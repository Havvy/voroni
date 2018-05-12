#![allow(unused_imports)]
use std::cell::{Ref, RefCell};
use std::cmp::Ordering;
use std::collections::binary_heap::BinaryHeap;
use std::collections::{BTreeMap, HashSet};
use std::f64::{NEG_INFINITY, INFINITY};
use std::hash::{Hash, Hasher};
use std::iter::FromIterator;
use std::ops::Index;
use std::rc::{Rc, Weak};

use ordered_float::{NotNaN};
use half_edge::{ConnectivityKernel};

use ::geometry::{LineSegment, Point, Triangle, TriangleOrientation, FocusDirectrixParabola, QuadraticParabola};
use self::Event::{Site, Circle};

type Number = NotNaN<f64>;
type Sweepline = Number;


/// An arc on the beachline.
///
/// A parabola with the specified focus and a directix of the scanline.
/// Since the scanline changes as the algorithm, this represents an entire
/// class of parabolas, though only one parabola at any point in time.
#[derive(Debug, Clone)]
struct Arc {
    /// The focus of the parabolas where the directrix is the sweepline.
    focus: Point,

    /// Associated circle event to the arc.
    circle_event: Option<Event>,

    /// Pointer up the tree.
    parent: Option<Weak<RefCell<ArcTreeComponent>>>
}

impl Arc {
    fn parabola_at_sweep_location(&self, sweep_location: &Sweepline) -> FocusDirectrixParabola {
        FocusDirectrixParabola::new_unwrap(self.focus, *sweep_location)
    }

    fn siblings(&self) -> (Option<&Arc>, Option<&Arc>) {
    // VParabola * lp = VParabola::GetLeftParent (b);
    // VParabola * rp = VParabola::GetRightParent(b);

    // VParabola * a  = VParabola::GetLeftChild (lp);
    // VParabola * c  = VParabola::GetRightChild(rp);
    unimplemented!()
    }

    fn compute_circle_event(&mut self) -> Option<Event> {
        unimplemented!()
    }
}

impl PartialEq for Arc {
    fn eq(&self, other: &Self) -> bool {
        self as *const Self as usize == other as *const Self as usize
    }
}

impl Eq for Arc {}

#[derive(Debug)]
struct Breakpoint {
    left: ArcTree,
    right: ArcTree,
    parent: Option<Weak<RefCell<ArcTreeComponent>>>
}

impl Breakpoint {
    fn new(left: ArcTree, right: ArcTree, parent: Option<Weak<RefCell<ArcTreeComponent>>>) -> Breakpoint {
        let breakpoint = Breakpoint {
            left: left,
            right: right,
            parent
        };

        // FIXME(Havvy)

        breakpoint
    }

    fn with_left_arc<F, Output>(&self, f: F) -> Output where F: Fn(&Arc) -> Output {
        fn right_descent<F, Output>(arc_tree: &ArcTree, f: F) -> Output where F: Fn(&Arc) -> Output {
            match *(arc_tree.0.borrow()) {
                ArcTreeComponent::Leaf(ref arc) => f(arc),
                ArcTreeComponent::Branch(ref breakpoint) => right_descent(&breakpoint.right, f),
            }
        }

         match *self.left.0.borrow() {
            ArcTreeComponent::Leaf(ref arc) => f(arc),
            ArcTreeComponent::Branch(ref breakpoint) => right_descent(&breakpoint.left, f),
        }
    }

    fn with_right_arc<F, Output>(&self, f: F) -> Output where F: Fn(&Arc) -> Output {
        fn left_descent<F, Output>(arc_tree: &ArcTree, f: F) -> Output where F: Fn(&Arc) -> Output {
            match *(arc_tree.0.borrow()) {
                ArcTreeComponent::Leaf(ref arc) => f(arc),
                ArcTreeComponent::Branch(ref breakpoint) => left_descent(&breakpoint.left, f),
            }
        }

         match *self.left.0.borrow() {
            ArcTreeComponent::Leaf(ref arc) => f(arc),
            ArcTreeComponent::Branch(ref breakpoint) => left_descent(&breakpoint.right, f),
        }
    }

    fn with_arcs<F, Output>(&self, f: F) -> Output where F: Fn(&Arc, &Arc) -> Output {
        self.with_left_arc(|left| {
            self.with_right_arc(|right| {
                f(left, right)
            })
        })
    }

    /// Calculate the 'x' coordinate where the breakpoint is at the
    /// specified sweep location.
    fn get_x(&self, sweep_location: &Sweepline) -> Number {
        // The breakpoint is the intersection of the two parabolas defined by
        // the foci of the arcs and the given sweep location.
        //
        // In general, we can find the intersection of two functions f(x) and 
        // g(x) by creating a function h(x) = f(x) - g(x) and then solving for
        // h(x) = 0.
        //
        // We do this by converting our focus-directrix description of the
        // parabolas into quadratic form (f(x) = ax^2 + bx + c), subtracting
        // them, and then using the quadratic equation to find zeros.
        //
        // There are two zeros, so we pick the one that corresponds to the
        // parabola segments we are looking at.
        // FIXME(Havvy): How does the `if` actually work? 
        self.with_arcs(|left, right| {
            let left_parabola: QuadraticParabola = left.parabola_at_sweep_location(sweep_location).into();
            let right_parabola: QuadraticParabola = right.parabola_at_sweep_location(sweep_location).into();

            let intersection_at_zero_parabola = &left_parabola - &right_parabola;

            let (q1, q2) = intersection_at_zero_parabola.quadratic_equation();

            if left.focus.y < right.focus.y {
                q1.max(q2)
            } else {
                q1.min(q2)
            }
        })
    }
}

impl PartialEq for Breakpoint {
    fn eq(&self, other: &Self) -> bool {
        self as *const Self as usize == other as *const Self as usize
    }
}

impl Eq for Breakpoint {}

/// Tree where the branches are arcs and the leaves are breakpoints.
#[derive(Debug)]
struct ArcTree(Rc<RefCell<ArcTreeComponent>>);

impl ArcTree {
    fn new_root(site: Point) -> ArcTree {
        ArcTree::new_from_component(ArcTreeComponent::Leaf(Arc { focus: site, circle_event: None, parent: None}))
    }

    fn new_triple((a1, a2, a3): (Arc, Arc, Arc)) -> ArcTree {
        let a1 = ArcTree::new_from_component(ArcTreeComponent::Leaf(a1));
        let a2 = ArcTree::new_from_component(ArcTreeComponent::Leaf(a2));
        let a3 = ArcTree::new_from_component(ArcTreeComponent::Leaf(a3));

        let mut right_branch = ArcTree::new_from_component(ArcTreeComponent::Branch(Breakpoint::new(a2, a3, None)));
        right_branch.set_childrens_parent();

        let mut left_branch = ArcTree::new_from_component(ArcTreeComponent::Branch(Breakpoint::new(a1, right_branch, None)));
        left_branch.set_childrens_parent();

        left_branch
    }

    fn new_from_component(tree_component: ArcTreeComponent) -> ArcTree {
        ArcTree(Rc::new(RefCell::new(tree_component)))
    }

    fn with_arc_above_site_mut<F, Output>
    (&mut self, site: Point, sweep_location: &Sweepline, mut f: F) -> Output
    where F: FnMut(&mut ArcTree) -> Output
    {
        match &mut *self.0.borrow_mut() {
            // What I want to do, but the borrow checker won't let me.
            // Instead, f(self) is at the bottom.
            // leaf @ &mut ArcTreeComponent::Leaf(..) => f(self),
            &mut ArcTreeComponent::Leaf(..) => {}
            &mut ArcTreeComponent::Branch(ref mut breakpoint) => {
                return if breakpoint.get_x(sweep_location) > site.x {
                    breakpoint.left.with_arc_above_site_mut(site, sweep_location, f)
                } else {
                    breakpoint.right.with_arc_above_site_mut(site, sweep_location, f)
                }
            },
        };

        f(self)
    }

    /// Set this ArcTree's children ArcTree's parents to this ArcTree. 
    fn set_childrens_parent(&mut self) {
        match *self.0.borrow_mut() {
            ArcTreeComponent::Branch(ref mut breakpoint) => {
                breakpoint.left.set_parent(Rc::downgrade(&self.0));
                breakpoint.right.set_parent(Rc::downgrade(&self.0));
            },
            ArcTreeComponent::Leaf(..) => { /* no children */ }
        };
    }

    fn set_parent_opt(&mut self, parent: Option<Weak<RefCell<ArcTreeComponent>>>) {
        match *self.0.borrow_mut() {
            ArcTreeComponent::Leaf(ref mut arc) => {
                arc.parent = parent;
            },
            ArcTreeComponent::Branch(ref mut breakpoint) => {
                breakpoint.parent = parent;
            }
        };
    }

    /// Set the parent of this ArcTree's contents to the specified parent.
    fn set_parent(&mut self, parent: Weak<RefCell<ArcTreeComponent>>) {
        self.set_parent_opt(Some(parent));
    }
}

/// The implemented tree for the beachline.
#[derive(Debug, PartialEq, Eq)]
enum ArcTreeComponent {
    Branch(Breakpoint),
    Leaf(Arc)
}

impl ArcTreeComponent {
    fn as_branch(&self) -> &Breakpoint {
        if let ArcTreeComponent::Branch(ref breakpoint) = *self {
            breakpoint
        } else {
            panic!("Method `as_branch` called on non-branch arc tree component.");
        }
    }

    fn parent(&self) -> Option<Weak<RefCell<ArcTreeComponent>>> {
        match *self {
            ArcTreeComponent::Branch(ref breakpoint) => breakpoint.parent.clone(),
            ArcTreeComponent::Leaf(ref arc) => arc.parent.clone(),
        }
    }

    // Get the fucking fucker fuck that fucking holds the fucking left path towards the fucking fuck parabola.
    fn with_left_parent<F>(&self, mut f: F) where F: FnMut(Option<&Breakpoint>) {
        match self.parent().and_then(|weak_parent| weak_parent.upgrade()) {
            None => { f(None); }
            Some(ref parent_rc) => {
                let parent_atc: Ref<ArcTreeComponent> = parent_rc.borrow();
                let parent_bp: Ref<Breakpoint> = Ref::map(Ref::clone(&parent_atc), |atc| atc.as_branch());
                let left: Ref<ArcTreeComponent> = parent_bp.left.0.borrow();

                if &*left == self {
                    parent_atc.with_left_parent(f);
                } else {
                    f(Some(self.as_branch()));
                }
            }
        };
    }

    fn with_right_parent<F>(&self, mut f: F) where F: FnMut(Option<&Breakpoint>) {
        match self.parent().and_then(|weak_parent| weak_parent.upgrade()) {
            None => { f(None); }
            Some(ref parent_rc) => {
                let parent_atc: Ref<ArcTreeComponent> = parent_rc.borrow();
                let parent_bp: Ref<Breakpoint> = Ref::map(Ref::clone(&parent_atc), |atc| atc.as_branch());
                let right: Ref<ArcTreeComponent> = parent_bp.right.0.borrow();

                if &*right == self {
                    parent_atc.with_left_parent(f);
                } else {
                    f(Some(self.as_branch()));
                }
            }
        };
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
struct ArcRef;

/// The beachline knows what all the arcs and breakpoints are.
///
/// It wraps an `ArcTree`, knowing the Voronoi specific operations, while
/// the `ArcTree` only knows about tree operations.
#[derive(Debug)]
struct Beachline(ArcTree);

impl Beachline {
    /// Construct a new Beachline. This is only done at the beginning of
    /// Fortune's Algorithm while handling the first point.
    fn new(site: Point) -> Beachline {
        Beachline(ArcTree::new_root(site))
    }

    /// Replace the arc currently located above `Point` as the specified sweep
    /// location with three arcs.
    ///
    /// The left and right arcs have the same focus as the one being removed,
    /// while the middle one has the specified site as its focus.
    ///
    /// Returns the invalid circle event associated with the old arc, if it
    /// existed.
    fn split_at_site(&mut self, site: Point, sweep_location: &Sweepline) -> Option<Event> {
        (self.0).with_arc_above_site_mut(site, sweep_location, |leaf| {
            let arc_triple = {
                let arc_ref = leaf.0.borrow();
                let arc: &Arc = match *arc_ref {
                    ArcTreeComponent::Leaf(ref arc) => arc,
                    _ => panic!("This should never happen.")
                };

                ((*arc).clone(), Arc { focus: site, circle_event: None, parent: None }, (*arc).clone())
            };

            let removed = ::std::mem::replace(leaf, ArcTree::new_triple(arc_triple));

            let (parent, invalid_circle_event) = match Rc::try_unwrap(removed.0).expect("There's only ever one stong Rc when splitting.").into_inner() {
                ArcTreeComponent::Leaf(Arc { parent, circle_event, ..}) => (parent, circle_event),
                _ => panic!("This should never happen.")
            };

            leaf.set_parent_opt(parent);

            invalid_circle_event
        })
    }
}

/// An event to be processed in Fortune's Algorithm.
///
/// Events are compared where the lowest `y` value is the greatest.
#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
enum Event {
    /// Handled when the sweepline is passing this point.
    Site(Point),

    /// Handled when the associated arc collapsed into a single point.
    /// May be removed from event queue from addition of another arc.
    Circle {
        /// The site of the circle's event.
        site: Point,

        /// The arc that has collapsed into a point when the site is handled.
        arc: ArcRef,

        /// The point the arc collasped at.
        vertex: Point
    }
}

impl Event {
    fn handle(self, event_queue: &mut EventQueue, beachline: &mut Beachline, sweep_location: &Sweepline, _voroni: &mut ConnectivityKernel) {
        match self {
            Site(site) => {
                let maybe_circle_event = beachline.split_at_site(site, sweep_location);

                if let Some(circle_event) = maybe_circle_event {
                    event_queue.remove(circle_event);
                }

                // let arc_left_event = arc_left.get_circle_event();
                // if let Some(ref event) = arc_left_event {
                //     event_queue.insert(event.clone())
                // };

                // // As a new 
                // beachline.insert(arc_center, None);

                // let arc_right_event = arc_right.get_circle_event();
                // if let Some(ref event) = arc_right_event {
                //     event_queue.insert(event.clone())
                // }
            },
            Circle { arc, site: _site, vertex: _vertex } => {
                // let (arc_left, arc_right) = beachline.remove_adjacent_arcs(&arc, event_queue);
                // beachline.remove_only(&arc);

                // breakpoints.remove(arc.left.expect("Left arc exists in circle event."));
                // breakpoints.remove(arc.right.expect("Right arc exists in circle event."));

                // VoronoiEdge e = new VoronoiEdge(ce.arc.left.site_left, ce.arc.right.site_right);
                // edgeList.add(e);

                // // Here we're trying to figure out if the Voronoi vertex we've found is the left
                // // or right point of the new edge.
                // // If the edges being traces out by these two beachline take a right turn then we know
                // // that the vertex is going to be above the current point
                // boolean turnsLeft = Point.ccw(arcLeft.right.edgeBegin, ce.p, arcRight.left.edgeBegin) == 1;
                // // So if it turns left, we know the next vertex will be below this vertex
                // // so if it's below and the slow is negative then this vertex is the left point
                // boolean isLeftPoint = (turnsLeft) ? (e.m < 0) : (e.m > 0);
                // if (isLeftPoint) {
                //     e.p1 = ce.vert;
                // }
                // else {
                //     e.p2 = ce.vert;
                // }
                // BreakPoint newBP = new BreakPoint(ce.arc.left.s1, ce.arc.right.s2, e, !isLeftPoint, this);
                // breakPoints.add(newBP);

                // arcRight.left = newBP;
                // arcLeft.right = newBP;

                // checkForCircleEvent(arcLeft);
                // checkForCircleEvent(arcRight); 
                unimplemented!()
            }
        };
    }

    fn y(&self) -> Number {
        match *self {
            Site(Point { y, .. }) => y,
            Circle { site: Point { y, .. }, .. } => y
        }
    }
}

/// Event Queue that allows removing events.
struct EventQueue {
    queue: BinaryHeap<Event>,
    removed_events: HashSet<Event>
}

impl EventQueue {
    fn new<I>(sites: I) -> EventQueue where I: Iterator<Item=Event> {
        EventQueue {
            queue: BinaryHeap::from_iter(sites),
            removed_events: HashSet::new()
        }
    }

    fn insert(&mut self, event: Event) {
        self.queue.push(event);
    }

    fn pop(&mut self) -> Option<Event> {
        loop {
            match self.queue.pop() {
                Some(event) => {
                    if !self.removed_events.remove(&event) {
                        break Some(event)
                    }
                },

                None => break None
            }
        }
    }

    fn remove(&mut self, event: Event) {
        self.removed_events.insert(event);
    }
}

pub fn fortune(sites: &[Point]) -> ConnectivityKernel {
    let mut voroni = ConnectivityKernel::new();
    let mut sweep_location: Sweepline;
    let mut event_queue = EventQueue::new(sites.into_iter().map(|&site| Site(site)));
    let mut beachline;

    // Handle the first point in its own way.
    match event_queue.pop() {
        Some(Site(site)) => {
            // sweep_location = site.y;
            beachline = Beachline::new(site);
        },
        None => {
            // No points were given.
            return voroni;
        },
        _ => panic!("First event must be a Site event if it exists.")
    };

    while let Some(event) = event_queue.pop() {
        sweep_location = event.y();
        event.handle(&mut event_queue, &mut beachline, &sweep_location, &mut voroni);
    }

    voroni
}