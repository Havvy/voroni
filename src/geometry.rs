use std::f64::{INFINITY, EPSILON};
use std::hash::{Hash, Hasher};
use std::ops::{Mul};
use std::{fmt};

use ::ordered_float::{NotNaN, FloatIsNaN};
use ::itertools;

/// Basic Point type for usage in the Voroni lib.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Point {
    pub x: NotNaN<f64>, 
    pub y: NotNaN<f64>
}

impl Point {
    /// Create a new point. Returns Err(FloatIsNaN) if either paramter is NaN.
    pub fn new(x: f64, y: f64) -> Result<Point, FloatIsNaN> {
        let x = NotNaN::new(x);
        let y = NotNaN::new(y);

        match (x, y) {
            (Ok(x), Ok(y)) => Ok(Point{x, y}),
            _ => x.and(y).map(|_| -> Point { unreachable!() })
        }
    }

    /// Like new, but panics instead of return an Err.
    pub fn new_unwrap(x: f64, y: f64) -> Point {
        Point::new(x, y).expect("Points cannot have NaN values.")
    }
}

impl<'a> Mul<&'a NotNaN<f64>> for Point {
    type Output = Point;
    fn mul(self, scalar: &NotNaN<f64>) -> Point {
        Point {
            x: self.x * *scalar,
            y: self.y * *scalar
        }
    }
}

impl fmt::Display for Point {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "({}, {})", self.x.into_inner(), self.y.into_inner())
    }
}

impl fmt::Debug for Point {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Point {{ x: {}, y: {} }}", self.x.into_inner(), self.y.into_inner())
    }
}

#[derive(Clone, Copy, Debug, Eq)]
pub struct LineSegment {
    pub from: Point, 
    pub to: Point
}

impl LineSegment {
    pub fn slope(&self) -> NotNaN<f64> {
        match self.dx() / self.dy() {
            x if x.into_inner() == -INFINITY => NotNaN::new(INFINITY).expect("Infinity is not NaN."),
            x => x
        }
    }

    fn perpendicular_bisector(&self) -> Line {
        let midpoint = self.midpoint();
        let to = Point { x: midpoint.x + self.dy(), y: midpoint.y - self.dx() };

        Line(LineSegment{from: midpoint, to})
    }

    pub fn length(&self) -> NotNaN<f64> {
        NotNaN::new((self.dx().powi(2) + self.dy().powi(2)).sqrt()).expect("Math works")
    }

    pub fn midpoint(&self) -> Point {
        Point {
            x: (self.from.x + self.to.x) / 2f64, 
            y: (self.from.y + self.to.y) / 2f64
        }
    }

    pub fn dx(&self) -> NotNaN<f64> {
        self.to.x - self.from.x
    }

    pub fn dy(&self) -> NotNaN<f64> {
        self.to.y - self.from.y
    }
}

impl PartialEq for LineSegment {
    fn eq(&self, other: &LineSegment) -> bool {
        (self.from == other.from && self.to == other.to)
        || (self.from == other.to && self.to == other.from)
    }
}

impl<'a> Mul<&'a NotNaN<f64>> for LineSegment {
    type Output = LineSegment;
    fn mul(self, scalar: &NotNaN<f64>) -> LineSegment {
        LineSegment {
            from: self.from * scalar,
            to: self.to * scalar
        }
    }
}


impl fmt::Display for LineSegment {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "\\{{({}, {}) -> ({}, {})}}", self.from.x, self.from.y, self.to.x, self.to.y)
    }
}

impl Hash for LineSegment {
    fn hash<H>(&self, hasher: &mut H) where H: Hasher  {
        if self.from < self.to {
            self.from.hash(hasher);
            self.to.hash(hasher);
        } else {
            self.to.hash(hasher);
            self.from.hash(hasher);
        }
    }
}

#[derive(Clone, Copy, Debug, Eq)]
struct Line(LineSegment);

impl Line {
    /// https://en.wikipedia.org/w/index.php?title=Line%E2%80%93line_intersection&oldid=764545912#Intersection_of_two_lines
    fn intersection(&self, other: &Line) -> Option<Point> {
        let selfs = self.0; // self as Segment
        let other = other.0;

        let denominator_add = (selfs.from.x - selfs.to.x) * (other.from.y - other.to.y);
        let denominator_sub = (selfs.from.y - selfs.to.y) * (other.from.x - other.to.x);
        let denominator = denominator_add - denominator_sub;

        if denominator.abs() < (2f64 * EPSILON) {
            return None;
        }

        let numerator_x_add = ((selfs.from.x * selfs.to.y) - (selfs.from.y * selfs.to.x)) * (other.from.x - other.to.x);
        let numerator_x_sub = ((other.from.x * other.to.y) - (other.from.y * other.to.x)) * (selfs.from.x - selfs.to.x);

        let numerator_y_add = ((selfs.from.x * selfs.to.y) - (selfs.from.y * selfs.to.x)) * (other.from.y - other.to.y);
        let numerator_y_sub = ((other.from.x * other.to.y) - (other.from.y * other.to.x)) * (selfs.from.y - selfs.to.y);

        Some(Point {
            x: (numerator_x_add - numerator_x_sub) / denominator,
            y: (numerator_y_add - numerator_y_sub) / denominator
        })
    }

    fn slope(&self) -> NotNaN<f64> {
        self.0.slope()
    }
}

impl PartialEq for Line {
    fn eq(&self, other: &Line) -> bool {
        self.slope() == other.slope() && self.intersection(other) == None
    }
}

impl fmt::Display for Line {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?})", self)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
struct Circle {
    center: Point,
    radius: NotNaN<f64>
}

impl Circle {
    fn contains_point(&self, point: Point) -> bool {
        (LineSegment {from: self.center, to: point}).length() < self.radius
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Triangle(pub Point, pub Point, pub Point);

impl Triangle {
    pub fn super_triangle_of(point: &[Point]) -> Triangle {
        unimplemented!();
    }

    pub fn edges(&self) -> TriangleEdges {
        TriangleEdges {
            triangle: self, 
            ix: 0
        }
    }

    pub fn vertexes(&self) -> TriangleVertexes {
        TriangleVertexes {
            triangle: self,
            ix: 0
        }
    }

    pub fn circumscribes_point(&self, point: Point) -> bool {
        self.circumcircle().contains_point(point)
    }

    fn circumcircle(&self) -> Circle {
        let mut edges = self.edges();
        let bisector_1 = edges.next().expect("Triangle has 3 edges.").perpendicular_bisector();
        let bisector_2 = edges.next().expect("Triangle has 3 edges.").perpendicular_bisector();

        // let center = bisector_1.intersection(&bisector_2).expect("There's no degenerate triangles.");
        let center = match bisector_1.intersection(&bisector_2) {
            Some(center) => center,
            None => {
                println!("degenerate triangle found:\n  {:?}", self);
                println!("bisectors:\n  {:?}\n  {:?}", bisector_1, bisector_2);
                panic!("degenerate triangle found.");
            }
        };

        let radius = (LineSegment{from: center, to: self.0}).length();

        Circle{center, radius}
    }

    pub fn shares_vertex(&self, other: &Triangle) -> bool {
        let Triangle(v0, v1, v2) = *self;
        let Triangle(o0, o1, o2) = *other;

        v0 == o0 || v0 == o1 || v0 == o2 ||
        v1 == o0 || v1 == o1 || v1 == o2 ||
        v2 == o0 || v2 == o1 || v2 == o2
    }
}

impl<'a> Mul<&'a NotNaN<f64>> for Triangle {
    type Output = Triangle;
    fn mul(self, scalar: &NotNaN<f64>) -> Triangle {
        Triangle(self.0 * scalar, self.1 * scalar, self.2 * scalar)
    }
}

impl fmt::Display for Triangle {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "△{{({}, {}), ({}, {}), ({}, {})}}", self.0.x, self.0.y, self.1.x, self.1.y, self.2.x, self.2.y)
    }
}

#[derive(Clone, Copy, Debug)]
pub struct TriangleEdges<'t>{
    triangle: &'t Triangle, 
    ix: u8
}

impl<'t> Iterator for TriangleEdges<'t> {
    type Item = LineSegment;

    fn next(&mut self) -> Option<<Self as Iterator>::Item> {
        if self.ix < 4 { self.ix += 1; }

        match self.ix {
            1 => Some(LineSegment{from: self.triangle.0, to: self.triangle.1}),
            2 => Some(LineSegment{from: self.triangle.1, to: self.triangle.2}),
            3 => Some(LineSegment{from: self.triangle.2, to: self.triangle.0}),
            _ => None
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct TriangleVertexes<'t>{
    triangle: &'t Triangle, 
    ix: u8
}

impl<'t> Iterator for TriangleVertexes<'t> {
    type Item = Point;

    fn next(&mut self) -> Option<<Self as Iterator>::Item> {
        if self.ix < 4 { self.ix += 1; }

        match self.ix {
            1 => Some(self.triangle.0),
            2 => Some(self.triangle.1),
            3 => Some(self.triangle.2),
            _ => None
        }
    }
}

#[derive(Clone, Debug, PartialEq, Hash)]
pub struct TriangleMesh(Vec<Triangle>);

impl TriangleMesh {
    pub fn add_triangle(&mut self, t: Triangle) {
        self.0.push(t);
    }

    pub fn remove_triangle(&mut self, t: Triangle) {
        let maybe_index = self.0.iter().position(|triangle_in_mesh| *triangle_in_mesh == t);
        
        if let Some(index) = maybe_index {
            self.0.remove(index);
        }
    }

    pub fn remove_triangles_with_vertexes_of_triangle(&mut self, pattern: Triangle) {
        self.0.retain(|triangle| !pattern.shares_vertex(triangle));
    }

    pub fn iter<'tm>(&'tm self) -> ::std::slice::Iter<'tm, Triangle> {
        self.0.iter()
    }
}

impl<'a> Mul<&'a NotNaN<f64>> for TriangleMesh {
    type Output = TriangleMesh;
    fn mul(self, scalar: &NotNaN<f64>) -> TriangleMesh {
        TriangleMesh(self.0.into_iter().map(|t| t * scalar).collect())
    }
}

impl Default for TriangleMesh {
    fn default() -> TriangleMesh {
        TriangleMesh(vec![])
    }
}

impl IntoIterator for TriangleMesh {
    type Item = Triangle;
    type IntoIter = ::std::vec::IntoIter<Triangle>;

    fn into_iter(self) -> ::std::vec::IntoIter<Triangle> {
        self.0.into_iter()
    }
}

impl fmt::Display for TriangleMesh {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "TriangleMesh{{\n    {}\n}}", itertools::join(self.0.iter(), ",\n    "))
    }
}

#[cfg(test)]
mod test {
    use super::*;
    
    fn generate_super_triangle() -> Triangle {
        Triangle(
            Point{x: NotNaN::new(-1f64).unwrap(), y: NotNaN::new(-1f64).unwrap()},
            Point{x: NotNaN::new(-1f64).unwrap(), y: NotNaN::new(3f64).unwrap()},
            Point{x: NotNaN::new(3f64).unwrap(), y: NotNaN::new(-1f64).unwrap()}
        )
    }

    #[test]
    fn perpendicular_bisector() {
        assert!(
            LineSegment{from: Point::new_unwrap(0.0, 0.0), to: Point::new_unwrap(0.0, 2.0)}.perpendicular_bisector()
            ==
            Line(LineSegment{from: Point::new_unwrap(0.0, 1.0), to: Point::new_unwrap(1.0, 1.0)})
        );

        assert!(
            LineSegment{from: Point::new_unwrap(0.0, 0.0), to: Point::new_unwrap(2.0, 0.0)}.perpendicular_bisector()
            ==
            Line(LineSegment{from: Point::new_unwrap(1.0, 0.0), to: Point::new_unwrap(1.0, 1.0)})
        );

        assert!(
            LineSegment{from: Point::new_unwrap(-1.0, -1.0), to: Point::new_unwrap(1.0, 1.0)}.perpendicular_bisector()
            ==
            Line(LineSegment{from: Point::new_unwrap(0.0, 0.0), to: Point::new_unwrap(-1.0, 1.0)})
        );

        assert!(
            LineSegment{from: Point::new_unwrap(0.0, 0.0), to: Point::new_unwrap(6.0, 4.0)}.perpendicular_bisector()
            ==
            Line(LineSegment{from: Point::new_unwrap(3.0, 2.0), to: Point::new_unwrap(7.0, -4.0)})
        );

        assert!(
            LineSegment{from: Point::new_unwrap(0.0, 0.0), to: Point::new_unwrap(6.0, 4.0)}.perpendicular_bisector()
            ==
            Line(LineSegment{from: Point::new_unwrap(3.0, 2.0), to: Point::new_unwrap(0.0, 6.5)})
        );
    }

    #[test]
    fn circumcircle_1() {
        assert!(generate_super_triangle().circumscribes_point(Point::new_unwrap(0.0, 0.0)));
    }

    #[test]
    fn circumcircle_2() {
        // △{(-1, -1), (-1, 4), (0, 0)}C(1, 0)=true
        let triangle = Triangle(
            Point::new_unwrap(0.0, 0.0),
            Point::new_unwrap(-2.0, -2.0),
            Point::new_unwrap(-2.0, 8.0)
        );

        let triangle_turned = Triangle(triangle.1, triangle.2, triangle.0);

        assert!(!triangle.circumscribes_point(Point::new_unwrap(2.0, 0.0)));
        assert!(!triangle_turned.circumscribes_point(Point::new_unwrap(2.0, 0.0)));
    }

    #[test]
    fn line_segment_equality() {
        let origin = Point::new_unwrap(0.0, 0.0);
        let x_hat = Point::new_unwrap(1.0, 0.0);
        // let y_hat = Point::new_unwrap(0.0, 1.0);

        assert!(LineSegment{from: origin, to: x_hat} == LineSegment{from: x_hat, to: origin});
        assert!(LineSegment{from: origin, to: x_hat} == LineSegment{from: origin, to: x_hat});
    }

    #[test]
    fn intersection_1() {
        let z = NotNaN::new(0f64).unwrap();
        let o = NotNaN::new(1f64).unwrap();
        let t = NotNaN::new(2f64).unwrap();

        let ls1 = LineSegment{from: Point{x: z, y: t}, to: Point{x: t, y: z}};
        let ls2 = LineSegment{from: Point{x: z, y: z}, to: Point{x: o, y: o}};

        let intersection = Point{x: o, y: o};

        assert!(Line(ls1).intersection(&Line(ls2)) == Some(intersection));
    }

    #[test]
    fn shares_vertex() {
        let nums: Vec<NotNaN<f64>> = vec![0f64, 1f64, 2f64, 3f64, 4f64].into_iter().map(NotNaN::new).map(|o| o.unwrap()).collect();
        let t1 = Triangle(Point{x: nums[0], y: nums[0]}, Point{x: nums[1], y: nums[0]}, Point{x: nums[0], y: nums[1]});
        let t2 = Triangle(Point{x: nums[2], y: nums[2]}, Point{x: nums[1], y: nums[0]}, Point{x: nums[3], y: nums[0]});

        assert!(t1.shares_vertex(&t2));
    }
}