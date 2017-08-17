use ::std::collections::binary_heap::BinaryHeap;
use ::std::collections::{BTreeMap, HashSet};

use ordered_float::{FloatIsNaN, NotNaN};
use half_edge::{};

use ::geometry::{Point};

type EventQueue = BinaryHeap<Event>;

/// An event to be processed in Fortune's Algorithm.
/// Events are compared where the lowest `y` value is the greatest.
enum Event {
    Site(Point),
    Circle {
        arc: Arc, 
        point: NotNaN<f32>,
        vertex: Point
}

impl Event {
    fn handle(self, event_queue: &mut EventQueue, &mut beachline: Beachline) {
        match self {
            Site(site) => {
                let arc_above_site = beachline.keys().find(|&&arc| arc.is_point_under(site));

                // Deal with the degenerate case where the first two points are at the same y value
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
                // CircleEvent falseCE = arcEntryAbove.getValue();
                // if (falseCE != null) {
                //     events.remove(falseCE);
                if let Some(circle_event) = beachline[arc_above_site] {
                    event_queue.
                }
                

            },
            Circle { arc, point } => remove_parabola(arc, point)
        };
    }

    fn y(self) -> NotNaN<f32> {
        match self {
            Site(Point{y, ...}) => y,
            Circle { point: y, ... } => y
        }
    }
}

impl PartialOrd for Event {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Event {
    fn cmp(&self, other: &Self) -> Ordering {
        // Note the ! operand.
        !(self.y().partial_cmp(other.y()))
    }
}

impl From<Point> for Event {
    fn from(p: Point) -> Event {
        Site(p)
    }
}

struct Parabola {
    x_left: NotNaN<f32>,
    x_right: NotNaN<f64>
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum EdgeDirection{Left, Right}

struct Breakpoint {
    s1: Point,
    s2: Point,
    e: VoroniEdge,
    direction: EdgeDirection,
    edgeBig: Point,
    sweep_location_cache: NotNaN<f32>,
    point_cache: Point
}

/// The beachline is a map from Arcs to their associated Circle Events, if they have one.
#[derive(Debug)]
struct Beachline(BTreeMap<Arc, Option<Event>>);

impl Beachline {
    fn new() -> Beachline {
        Beachline(BTreeMap::new())
    }

    // TODO(Havvy): Is this useless by the way I refactored the code?
    fn is_arcless(&self) -> bool {
        self.0.is_empty()
    }

    fn insert(&mut self, arc: Arc, event: Option<Event>) {
        self.0.insert(arc, event);
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct Arc {
    site: Point
    left: Option<Breakpoint>,
    right: Option<Breakpoint>
}

impl Arc {
    fn first(site: Point) -> Arc {
        Arc { site, left: None, right: None }
    }

    /// Whether or not the point is located under the arc.
    fn is_point_under(&self, point: &Point) -> bool {
        if point.x >= self.left.expect("Left breakpoint exists").x && point.x <= self.right.expect("Right breakpoint exists").x {
            return true;
        }

        //if (myLeft.x == yourLeft.x && myRight.x == yourRight.x) return 0;
        //         if (myLeft.x >= yourRight.x) return 1;
        //         if (myRight.x <= yourLeft.x) return -1;

        // return Point.midpoint(myLeft, myRight).compareTo(Point.midpoint(yourLeft, yourRight));

        false
    }
}

//    function AddParabola ( point u )
//    {
//       par = arc under point u;
//       if (par has its circle event, when it is removed from the beachline)
//          remove this event from the queue
//       new arcs a, b, c;
//       b.site = u;
//       a.site = c.site = par.site; // site of arc is a focus of arc
//       xl, xr  = left and right edges, which comes from point on par under u
//       xl is a normal to  (a.site, b.site);
//       xr is a normal to (b.site, c.site);
//       replace par by the sequence a, xl, b, xr, c
//       CheckCircleEvent(a);
//       CheckCircleEvent(c);
//    }
fn add_parabola(site: Point, event_queue: &mut BinaryHeap) {
    
}

//    function RemoveParabola ( Parabola p )
//    {
//       l = an arc left to p;
//       r = an arc on the right from p;
//       if (l or r have their Circle events) remove these events from the queue
//       s = the circumcenter between l.site, p.site and r.site
//       x = new edge, starts at s, normal to (l.site, r.site)
//       finish two neighbour edges xl, xr at point s
//       replace a sequence xl, p, xr by new edge x
//       CheckCircleEvent(l);
//       CheckCircleEvent(r);
//    }
fn remove_parabola(parabola: Parabola, event_queue: &mut BinaryHeap) {

}

//    function CheckCircleEvent(Parabola p)
//    {
//       l = arc on the left to p;
//       r = arc on the right to p;
//       xl, xr = edges by the p
//       when there is no l  OR no r  OR  l.site=r.site  RETURN
//       s = middle point, where xl and xr cross each other
//       when there is no s (edges go like\ /) RETURN
//       r = distance between s an p.site (radius of curcumcircle)
//       if s.y + r is still under the sweepline  RETURN
//       e = new circle event
//       e.parabola = p;
//       e.y = s.y + r;
//       add e into queue 
//    }
fn check_circle_event(parabola: Parabola, event_queue: &mut BinaryHeap) {
    let arc_left = arc_left_of_parabola(parabola);
    let arc_right = arc_right_of_parabola(parabola);

    match (arc_left, arc_right) {
        (None, _) => { return; },
        (_, None) => { return; },
        (Arc{site: l_site, ..}, Arc{site: r_site, ..}) if l_site == r_site => { return; }
    }
}

fn fortune(sites: &[Point]) -> Vec<Edge> {
    // TODO: Consider BTreeSet?
    let mut event_queue = BinaryHeap::from_iter(sites.into_iter().map(|&site| Site(site)));
    let mut breakpoints = HashSet::new<Breakpoint>();
    let mut beachline = Beachline::new();
    let mut sweep_location: NotNaN<f32>;

    match event_queue.pop() {
        Some(Site(site)) => {
            beachline.insert(Arc::first(site), None);
        }
        None => return; // No points were given.
        _ => panic!("First event must be a Site event if it exists.")
    };

    while let Some(event) = event_queue.pop() {
        event.handle(&mut event_queue, &mut beachline);
    };
}