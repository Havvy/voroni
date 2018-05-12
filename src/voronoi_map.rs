use std::cell::Cell;
use std::cmp::Ordering;
use std::collections::binary_heap::BinaryHeap;
use std::collections::{BTreeMap, HashSet};
use std::f64::{NEG_INFINITY, INFINITY};
use std::hash::{Hash, Hasher};
use std::iter::FromIterator;
use std::ops::Index;
use std::rc::Rc;

use ordered_float::{NotNaN};
use half_edge::{ConnectivityKernel};

use ::geometry::{LineSegment, Point, Triangle, TriangleOrientation};
use self::Event::{Site, Circle};

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum EdgeDirection{Left, Right}

/// FIXME(Havvy): Transient struct? Should use the ConnectivityKernel instead?
#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Clone, Copy, Hash)]
struct VoroniEdge;

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
struct Sweepline(Cell<NotNaN<f64>>);

impl Sweepline {
    fn new(y: NotNaN<f64>) -> Sweepline {
        Sweepline(Cell::new(y))
    }

    fn set(&self, y: NotNaN<f64>) {
        self.0.set(y);
    }

    fn get(&self) -> NotNaN<f64> {
        self.0.get()
    }
}
/// The breakpoint of two coinciding arcs in the beachline
///
/// The struct only contains the information for calculating the breakpoint,
/// and the actual breakpoint is on the `point` method.
// TODO(Havvy): Derive PartialEq and Eq myself and ignore the cache.
#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Clone)]
struct ArcBreakpoint<'sl> {
    /// The site of the arc to the left on the breakpoint
    site_left: Point,
    /// The site of the arc to the right on the breakpoint
    site_right: Point,
    /// ???
    e: VoroniEdge,
    /// ???
    direction: EdgeDirection,
    /// ???
    edge_begin: Point,

    /// Where the sweep location is at all time. (This gets mutated despite being a non-mutable borrow.)
    sweep_location: &'sl Sweepline,

    /// Last location of the sweep location last time the breakpoint was queried.
    sweep_location_cache: Cell<NotNaN<f64>>,
    /// The last point returned last time the breakpoint was queried.
    point_cache: Cell<Point>
}

impl<'sl> ArcBreakpoint<'sl> {
    fn new(site_left: Point, site_right: Point, e: VoroniEdge, direction: EdgeDirection, sweep_location: &Sweepline) -> ArcBreakpoint {
        let zero = NotNaN::new(0f64).expect("0 is not NaN");
        let mut breakpoint = ArcBreakpoint {
            site_left, site_right, e, direction, sweep_location,

            edge_begin: Point { x: zero, y: zero },
            sweep_location_cache: Cell::new(zero),
            point_cache: Cell::new(Point { x: zero, y: zero })
        };
        breakpoint.edge_begin = breakpoint.point();
        breakpoint
    }

    /// Calculate the location of the breakpoint at the current sweep location.
    fn point(&self) -> Point {
        let sweep_location = self.sweep_location.get();
        // TODO(Havvy): Do with interior mutability?
        if sweep_location == self.sweep_location_cache.get() {
            return self.point_cache.get();
        }

        self.sweep_location_cache.set(sweep_location);

        // TODO(Havvy): Learn math and refactor.
        // TODO(Havvy): What happens if sweep location and y values are all the same when?
        let two = NotNaN::new(2f64).expect("Two is not NaN.");
        self.point_cache.set(if self.site_left.y == self.site_right.y {
            // Vertical line case.
            let x = (self.site_left.x + self.site_right.x) / two;

            let y_num = NotNaN::new((x - self.site_left.x).powi(2) + self.site_left.y.powi(2) - sweep_location.powi(2)).expect("Squaring cannot produce NaN.");
            let y_den = two * (self.site_left.y - *sweep_location);
            let y = y_num / y_den;

            Point { x, y }
        } else {
            // This method works by intersecting the line of the edge with the parabola of one of the higher site point.
            let site = if self.site_left > self.site_right { self.site_left } else { self.site_right };

            // TODO(Havvy): Math into geometry module
            let _d = two * (site.y - *sweep_location);

            //     double px = (s1.y > s2.y) ? s1.x : s2.x;
            //     double py = (s1.y > s2.y) ? s1.y : s2.y;
            //     double m = e.m;
            //     double b = e.b;

            //     double d = 2*(py - l);

            //     // Straight up quadratic formula
            //     double A = 1;
            //     double B = -2*px - d*m;
            //     double C = sq(px) + sq(py) - sq(l) - d*b;
            //     int sign = (s1.y > s2.y) ? -1 : 1;
            //     double det = sq(B) - 4 * A * C;
            //     // When rounding leads to a very very small negative determinant, fix it
            //     if (det <= 0) {
            //         x = -B / (2 * A);
            //     }
            //     else {
            //         x = (-B + sign * Math.sqrt(det)) / (2 * A);
            //     }
            //     y = m*x + b;

            unimplemented!()
        });

        // Recursive call that grabs the updated cached version.
        self.point()
    }
}

impl<'sl> Hash for ArcBreakpoint<'sl> {
    fn hash<H>(&self, hasher: &mut H) where H: Hasher  {
        // FIXME(Havvy): Some of these fields change?
        self.site_left.hash(hasher);
        self.site_right.hash(hasher);
        self.e.hash(hasher);
        self.direction.hash(hasher);
        self.edge_begin.hash(hasher);
    }
}

struct Breakpoints<'sl>(HashSet<Rc<ArcBreakpoint<'sl>>>);

impl<'sl> Breakpoints<'sl> {
    fn new() -> Breakpoints<'sl> {
        Breakpoints(HashSet::new())
    }

    fn insert(&mut self, breakpoint: ArcBreakpoint<'sl>) -> Rc<ArcBreakpoint<'sl>> {
        let rc = Rc::new(breakpoint);

        self.0.insert(rc.clone());

        rc
    }

    fn remove(&mut self, breakpoint: Rc<ArcBreakpoint<'sl>>) {
        self.0.remove(&breakpoint);
    }
}

/// An arc on the beachline.
///
/// The arc consists of the site that created it and the breakpoints with the arcs next to it.
///
/// The first arc cannot have breakpoints as no other arcs to break against.
/// Afterwards, only the leftmost arc's left breakpoint and rightmost arc's right breakpoint
/// will be None. 
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct Arc<'sl> {
    /// The location of the site event that created the arc.
    site: Point,
    left: Option<Rc<ArcBreakpoint<'sl>>>,
    right: Option<Rc<ArcBreakpoint<'sl>>>
}

impl<'sl> Arc<'sl> {
    fn first(site: Point) -> Arc<'sl> {
        Arc { site, left: None, right: None }
    }

    fn new(left: Option<Rc<ArcBreakpoint<'sl>>>, right: Option<Rc<ArcBreakpoint<'sl>>>) -> Arc<'sl> {
        let site = match (&left, &right) {
            (&Some(ref bp), _) => bp.site_right,
            (_, &Some(ref bp)) => bp.site_left,
            _ => { panic!("At least one breakpoint must exist."); }
        };

        Arc { left, right, site }
    }

    /// Whether or not the point is located under the arc.
    fn is_point_under(&self, point: &Point) -> bool {
        point.x >= self.left_breakpoint().x && point.x <= self.right_breakpoint().x
    }

    fn left_breakpoint(&self) -> Point {
        match self.left {
            Some(ref left) => left.point(),
            None => Point::new_unwrap(NEG_INFINITY, INFINITY)
        }
    }

    fn right_breakpoint(&self) -> Point {
        match self.left {
            Some(ref right) => right.point(),
            None => Point::new_unwrap(INFINITY, INFINITY)
        }
    }

    /// Get the center of the circle event's circle, if it exists.
    ///
    /// The outer arcs in the beachline do not have circle events.
    /// 
    /// Furthermore, the location of the arc's site w.r.t the adjacent sites
    /// can mean there is no circle event.
    fn circle_event_center(&self) -> Option<Point> {
        match *self {
            Arc { left: None, .. } => None,
            Arc { right: None, .. } => None,
            Arc { left: Some(ref left), right: Some(ref right), site } => {
                if Triangle(left.site_left, site, right.site_right).orientation() != TriangleOrientation::Counterclockwise {
                    None
                } else {
                    // FIXME(Havvy): return (this.left.getEdge().intersection(this.right.getEdge()));
                    unimplemented!()
                }
            }
        }
    }

    fn get_circle_event(&self) -> Option<Event<'sl>> {
        self.circle_event_center().map(|center| {
            let radius = (LineSegment { from: self.site, to: center }).length();
            Circle {
                site: Point { x: center.x, y: center.y - radius },
                arc: (*self).clone(),
                vertex: center
            }
        })
    }
}

impl<'sl> PartialOrd for Arc<'sl> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

// FIXME(Havvy): These are _false_ impls. Arcs cannot be ordered alone.
impl<'sl> Ord for Arc<'sl> {
    fn cmp(&self, other: &Self) -> Ordering {
        let my_left = self.left_breakpoint();
        let my_right = self.right_breakpoint();
        let other_left = other.left_breakpoint();
        let other_right = other.right_breakpoint();

        if my_left.x == other_left.x && my_right.x == other_right.x {
            Ordering::Equal
        } else if my_left.x >= other_right.x {
            Ordering::Greater
        } else if my_right.x <= other_left.x {
            Ordering::Less
        } else {
            let my_midpoint = (LineSegment { from: my_left, to: my_right }).midpoint();
            let other_midpoint = (LineSegment { from: other_left, to: other_right }).midpoint();

            my_midpoint.cmp(&other_midpoint)
        }
    }
}

/// The beachline is a map from Arcs to their associated Circle Events, if they have one.
#[derive(Debug)]
struct Beachline<'sl>(BTreeMap<Arc<'sl>, Option<Event<'sl>>>);

impl<'sl> Beachline<'sl> {
    fn new() -> Beachline<'sl> {
        Beachline(BTreeMap::new())
    }

    fn arcs<'bl>(&'bl self) -> ::std::collections::btree_map::Keys<'bl, Arc<'sl>, Option<Event<'sl>>> {
        self.0.keys()
    }

    fn remove_adjacent_arcs(&mut self, arc: &Arc<'sl>, event_queue: &mut EventQueue<'sl>) -> (Arc<'sl>, Arc<'sl>) {
        let arc_left = self
        .arcs()
        .filter(|beacharc| beacharc < &arc)
        .last()
        .expect("Only called with a circle event arc, so adjacent arcs must exist.")
        .clone();

        let arc_right = self
        .arcs()
        .filter(|beacharc| beacharc > &arc)
        .next()
        .expect("Only called with a circle event arc, so adjacent arcs must exist.")
        .clone();

        self.remove(&arc_left, event_queue);
        self.remove(&arc_right, event_queue);
        
        (arc_left, arc_right)
    }

    /// Put a new arc into the beachline with its associated circle event.
    /// When there aren't enough points for a circle event (less than 3?),
    /// then put None for the circle event.
    fn insert(&mut self, arc: Arc<'sl>, event: Option<Event<'sl>>) {
        self.0.insert(arc, event);
    }

    fn remove(&mut self, arc: &Arc<'sl>, event_queue: &mut EventQueue<'sl>) {
        let maybe_event = self.0.remove(arc).expect("Removed arc is in the beachline");

        if let Some(circle_event) = maybe_event {
            event_queue.remove(circle_event);
        }
    }

    fn remove_only(&mut self, arc: &Arc<'sl>) {
        self.0.remove(arc);
    }
}

impl<'a, 'sl> Index<&'a Arc<'sl>> for Beachline<'sl> {
    type Output = Option<Event<'sl>>;

    fn index(&self, arc: &'a Arc<'sl>) -> &Option<Event<'sl>> {
        &self.0[arc]
    }
}

/// An event to be processed in Fortune's Algorithm.
///
/// Events are compared where the lowest `y` value is the greatest.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
enum Event<'sl> {
    /// The sweepline is passing this point.
    Site(Point),

    /// An arc collapsed into a single point.
    Circle {
        /// The site of the circle's event.
        site: Point,

        /// The arc that has collapsed into a point when the site is handled.
        arc: Arc<'sl>,

        /// The point the arc collasped at.
        vertex: Point
    }
}

impl<'sl> Event<'sl> {
    fn handle(self, event_queue: &mut EventQueue<'sl>, beachline: &mut Beachline<'sl>, sweep_location: &'sl Sweepline, breakpoints: &mut Breakpoints<'sl>, _voroni: &mut ConnectivityKernel) {
        match self {
            Site(site) => {
                let arc_above_site = beachline
                .arcs()
                .find(|arc| arc.is_point_under(&site))
                .expect("There is always an arc above every non-initial site.")
                .clone();

                // Deal with the degenerate case where the first two points are at the same y value
                // ???(Havvy): Is this actually needed? The if check is wrong because size would never be 0.
                // ???(Havvy): If it is needed, should this be done outside of handle, perhaps in an e.g. handle_second?
                // if (arcs.size() == 0 && arcAbove.site.y == cur.p.y) {
                //     VoronoiEdge newEdge = new VoronoiEdge(arcAbove.site, cur.p);
                //     newEdge.p1 = new Point((cur.p.x + arcAbove.site.x)/2, Double.POSITIVE_INFINITY);
                //     BreakPoint newBreak = new BreakPoint(arcAbove.site, cur.p, newEdge, false, this);
                //     breakPoints.add(newBreak);
                //     this.edgeList.add(newEdge);
                //     Arc arcLeft = new Arc(null, newBreak, this);
                //     Arc arcRight = new Arc(newBreak, null, this);
                //     arcs.remove(arcAbove);
                //     arcs.put(arcLeft, null);
                //     arcs.put(arcRight, null);
                //     return;
                // }

                // Remove the circle event associated with this arc if there is one
                // Also remove the arc from the beachline, to be replaced by three other arcs below.
                beachline.remove(&arc_above_site, event_queue);

                // Create and insert 3 arcs to replace the removed arc ealier.
                // The center arc has 0 width at the current sweep line, but will grow out
                // towards the left and right arc as the sweep line moves. The left and right
                // arcs go from the respective left and right breakpoints of the removed arc to
                // the center arc.

                let left_arc_left_breakpoint = arc_above_site.left.clone();
                let right_arc_right_breakpoint = arc_above_site.right;
                let edge = VoroniEdge;
                // VoronoiEdge newEdge = new VoronoiEdge(arcAbove.site, cur.p);
                // this.edgeList.add(newEdge);
                let left_arc_right_breakpoint = ArcBreakpoint::new(arc_above_site.site, site, edge, EdgeDirection::Left, sweep_location);
                let left_arc_right_breakpoint = breakpoints.insert(left_arc_right_breakpoint);
                let left_arc_right_breakpoint = Some(left_arc_right_breakpoint);

                let right_arc_left_breakpoint = ArcBreakpoint::new(arc_above_site.site, site, edge, EdgeDirection::Right, sweep_location);
                let right_arc_left_breakpoint = breakpoints.insert(right_arc_left_breakpoint);
                let right_arc_left_breakpoint = Some(right_arc_left_breakpoint);

                let arc_left = Arc::new(left_arc_left_breakpoint, left_arc_right_breakpoint.clone());
                let arc_center = Arc::new(left_arc_right_breakpoint, right_arc_left_breakpoint.clone());
                let arc_right = Arc::new(right_arc_left_breakpoint, right_arc_right_breakpoint);

                let arc_left_event = arc_left.get_circle_event();
                if let Some(ref event) = arc_left_event {
                    event_queue.insert(event.clone())
                }
                beachline.insert(arc_left, arc_left_event);

                // As a new 
                beachline.insert(arc_center, None);

                let arc_right_event = arc_right.get_circle_event();
                if let Some(ref event) = arc_right_event {
                    event_queue.insert(event.clone())
                }
                beachline.insert(arc_right, arc_right_event);
            },
            Circle { arc, site: _site, vertex: _vertex } => {
                let (arc_left, arc_right) = beachline.remove_adjacent_arcs(&arc, event_queue);
                beachline.remove_only(&arc);

                breakpoints.remove(arc.left.expect("Left arc exists in circle event."));
                breakpoints.remove(arc.right.expect("Right arc exists in circle event."));

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

    fn y(&self) -> NotNaN<f64> {
        match *self {
            Site(Point{y, ..}) => y,
            Circle { site: Point { y, .. }, .. } => y
        }
    }
}

impl<'sl> PartialOrd for Event<'sl> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<'sl> Ord for Event<'sl> {
    fn cmp(&self, other: &Self) -> Ordering {
        // Note the ! operand.
        self.y().cmp(&other.y()).reverse()
    }
}

impl<'sl> From<Point> for Event<'sl> {
    fn from(p: Point) -> Event<'sl> {
        Site(p)
    }
}

struct EventQueue<'sl> {
    queue: BinaryHeap<Event<'sl>>,
    removed_events: HashSet<Event<'sl>>
}

impl<'sl> EventQueue<'sl> {
    fn new<I>(sites: I) -> EventQueue<'sl> where I: Iterator<Item=Event<'sl>> {
        EventQueue {
            queue: BinaryHeap::from_iter(sites),
            removed_events: HashSet::new()
        }
    }

    fn insert(&mut self, event: Event<'sl>) {
        self.queue.push(event);
    }

    fn pop(&mut self) -> Option<Event<'sl>> {
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

    fn remove(&mut self, event: Event<'sl>) {
        self.removed_events.insert(event);
    }
}


pub fn fortune(sites: &[Point]) -> ConnectivityKernel {
    let sweep_location: &Sweepline = &Sweepline::new(NotNaN::new(0f64).expect("Zero is not NaN"));
    let mut voroni = ConnectivityKernel::new();

    {
        let mut event_queue = EventQueue::new(sites.into_iter().map(|&site| Site(site)));
        let mut breakpoints = Breakpoints::new();
        let mut beachline = Beachline::new();

        // Deal with the first point specially.
        match event_queue.pop() {
            Some(Site(site)) => {
                sweep_location.set(site.y);
                beachline.insert(Arc::first(site), None);
            },
            None => {
                // No points were given.
                return voroni;
            },
            _ => panic!("First event must be a Site event if it exists.")
        };

        while let Some(event) = event_queue.pop() {
            sweep_location.set(event.y());
            event.handle(&mut event_queue, &mut beachline, &sweep_location, &mut breakpoints, &mut voroni);
        }
    }

    voroni
}