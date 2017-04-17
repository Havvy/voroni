// function BowyerWatson (pointList)
//    // pointList is a set of coordinates defining the points to be triangulated
//    triangulation := empty triangle mesh data structure
//    add super-triangle to triangulation // must be large enough to completely contain all the points in pointList
//    for each point in pointList do // add all the points one at a time to the triangulation
//       badTriangles := empty set
//       for each triangle in triangulation do // first find all the triangles that are no longer valid due to the insertion
//          if point is inside circumcircle of triangle
//             add triangle to badTriangles
//       polygon := empty set
//       for each triangle in badTriangles do // find the boundary of the polygonal hole
//          for each edge in triangle do
//             if edge is not shared by any other triangles in badTriangles
//                add edge to polygon
//       for each triangle in badTriangles do // remove them from the data structure
//          remove triangle from triangulation
//       for each edge in polygon do // re-triangulate the polygonal hole
//          newTri := form a triangle from edge to point
//          add newTri to triangulation
//    for each triangle in triangulation // done inserting points, now clean up
//       if triangle contains a vertex from original super-triangle
//          remove triangle from triangulation
//    return triangulation

#![allow(unused_variables)]
#![allow(dead_code)]

extern crate itertools;
extern crate ordered_float;

use std::collections::{HashSet};
use std::hash::{Hash};
use std::{cmp, fmt};

pub use ordered_float::{FloatIsNaN, NotNaN};

/// Basic Point type for usage in the Voroni lib.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Point {
    pub x: NotNaN<f64>, 
    pub y: NotNaN<f64>
}

impl Point {
    /// Create a new point. Returns Err(FloatIsNaN) if either paramter is NaN.
    fn new(x: f64, y: f64) -> Result<Point, FloatIsNaN> {
        let x = NotNaN::new(x);
        let y = NotNaN::new(y);

        match (x, y) {
            (Ok(x), Ok(y)) => Ok(Point{x, y}),
            _ => x.and(y).map(|_| -> Point { unreachable!() })
        }
    }

    /// Like new, but panics instead of return an Err.
    fn new_unwrap(x: f64, y: f64) -> Point {
        Point::new(x, y).expect("Points cannot have NaN values.")
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
struct Line(LineSegment);

impl Line {
    /// https://en.wikipedia.org/w/index.php?title=Line%E2%80%93line_intersection&oldid=764545912#Intersection_of_two_lines
    fn intersection(&self, other: &Line) -> Option<Point> {
        let selfs = self.0; // self as Segment
        let other = other.0;

        let denominator_add = (selfs.from.x - selfs.to.x) * (other.from.y - other.to.y);
        let denominator_sub = (selfs.from.y - selfs.to.y) * (other.from.x - other.to.x);
        let denominator = denominator_add - denominator_sub;

        if denominator.abs() < (2f64 * std::f64::EPSILON) {
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

impl std::fmt::Display for Line {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?})", self)
    }
}

#[derive(Clone, Copy, Debug, Eq)]
struct LineSegment {
    from: Point, 
    to: Point
}

impl LineSegment {
    fn slope(&self) -> NotNaN<f64> {
        match self.dx() / self.dy() {
            x if x.into_inner() == -std::f64::INFINITY => NotNaN::new(std::f64::INFINITY).expect("Infinity is not NaN."),
            x => x
        }
    }

    fn perpendicular_bisector(&self) -> Line {
        let midpoint = self.midpoint();
        let to = Point { x: midpoint.x + self.dy(), y: midpoint.y - self.dx() };

        Line(LineSegment{from: midpoint, to})
    }

    fn length(&self) -> NotNaN<f64> {
        NotNaN::new((self.dx().powi(2) + self.dy().powi(2)).sqrt()).expect("Math works")
    }

    fn midpoint(&self) -> Point {
        Point {
            x: (self.from.x + self.to.x) / 2f64, 
            y: (self.from.y + self.to.y) / 2f64
        }
    }

    fn dx(&self) -> NotNaN<f64> {
        self.to.x - self.from.x
    }

    fn dy(&self) -> NotNaN<f64> {
        self.to.y - self.from.y
    }
}

impl PartialEq for LineSegment {
    fn eq(&self, other: &LineSegment) -> bool {
        (self.from == other.from && self.to == other.to)
        || (self.from == other.to && self.to == other.from)
    }
}


impl fmt::Display for LineSegment {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "\\{{({}, {}) -> ({}, {})}}", self.from.x, self.from.y, self.to.x, self.to.y)
    }
}

impl Hash for LineSegment {
    fn hash<H>(&self, hasher: &mut H) where H: std::hash::Hasher  {
        if self.from < self.to {
            self.from.hash(hasher);
            self.to.hash(hasher);
        } else {
            self.to.hash(hasher);
            self.from.hash(hasher);
        }
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
    fn super_triangle_of(point: &[Point]) -> Triangle {
        unimplemented!();
    }

    fn edges(&self) -> TriangleEdges {
        TriangleEdges {
            triangle: self, 
            ix: 0
        }
    }

    fn vertexes(&self) -> TriangleVertexes {
        TriangleVertexes {
            triangle: self,
            ix: 0
        }
    }

    fn circumscribes_point(&self, point: Point) -> bool {
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

    fn shares_vertex(&self, other: &Triangle) -> bool {
        let Triangle(v0, v1, v2) = *self;
        let Triangle(o0, o1, o2) = *other;

        v0 == o0 || v0 == o1 || v0 == o2 ||
        v1 == o0 || v1 == o1 || v1 == o2 ||
        v2 == o0 || v2 == o1 || v2 == o2
    }
}

impl std::fmt::Display for Triangle {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "△{{({}, {}), ({}, {}), ({}, {})}}", self.0.x, self.0.y, self.1.x, self.1.y, self.2.x, self.2.y)
    }
}

#[derive(Clone, Copy, Debug)]
struct TriangleEdges<'t>{
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
struct TriangleVertexes<'t>{
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
    fn add_triangle(&mut self, t: Triangle) {
        self.0.push(t);
    }

    fn remove_triangle(&mut self, t: Triangle) {
        let maybe_index = self.0.iter().position(|triangle_in_mesh| *triangle_in_mesh == t);
        
        if let Some(index) = maybe_index {
            self.0.remove(index);
        }
    }

    fn remove_triangles_with_vertexes_of_triangle(&mut self, pattern: Triangle) {
        self.0.retain(|triangle| !pattern.shares_vertex(triangle));
    }

    fn iter<'tm>(&'tm self) -> std::slice::Iter<'tm, Triangle> {
        self.0.iter()
    }
}

impl Default for TriangleMesh {
    fn default() -> TriangleMesh {
        TriangleMesh(vec![])
    }
}

impl IntoIterator for TriangleMesh {
    type Item = Triangle;
    type IntoIter = std::vec::IntoIter<Triangle>;

    fn into_iter(self) -> std::vec::IntoIter<Triangle> {
        self.0.into_iter()
    }
}

impl std::fmt::Display for TriangleMesh {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "TriangleMesh{{\n    {}\n}}", itertools::join(self.0.iter(), ",\n    "))
    }
}

/// Returns a vec with only non-duplicated elements.
fn unique_elements_only<T>(vec: Vec<T>) -> Vec<T>
where
    T: Eq + Hash + Clone,
{
    let mut duplicated = HashSet::new();
    let mut scanned = HashSet::new();

    for element in &vec {
        if scanned.contains(element) {
            duplicated.insert(element.clone());
        }

        scanned.insert(element.clone());
    }

    vec.into_iter().filter(|element| { !duplicated.contains(element) }).collect()
}

/// Translate a list of cordinates into the Delaunay triangulation via the Boyer-Watson algorithm.
///
/// https://en.wikipedia.org/wiki/Delaunay_triangulation
/// https://en.wikipedia.org/wiki/Bowyer%E2%80%93Watson_algorithm
pub fn bower_watson (point_list: &[Point]) -> TriangleMesh {
    let super_triangle = Triangle::super_triangle_of(point_list);

    bower_watson_with_super_triangle(point_list, super_triangle)
}

pub fn bower_watson_with_super_triangle(point_list: &[Point], super_triangle: Triangle) -> TriangleMesh {
    // Start with an empty Triangle Mesh.
    let mut triangulation = TriangleMesh::default();

    // Add a Triangle that contains all points in the point list.
    triangulation.add_triangle(super_triangle);

    println!("Added super triangle.\n{}\n", triangulation);

    for point in point_list {
        println!("Adding {}", point);

        let point = point.clone();

        println!("Testing triangles:");

        // Triangles that are no longer valid with the new point inserted.
        let bad_triangles: Vec<Triangle> = triangulation.iter().filter(|triangle| {
            let res = triangle.circumscribes_point(point);
            println!("    {} -> {}.", triangle, res);
            res
        }).map(|t| *t).collect();

        println!("Removing triangles:\n    {}", itertools::join(bad_triangles.iter(), "\n    "));

        for triangle in &bad_triangles {
            triangulation.remove_triangle(*triangle);
        }

        // Create a polygon of edges from the triangles where they don't touch.
        let polygon = unique_elements_only(bad_triangles.iter().flat_map(Triangle::edges).collect());

        println!("Polygon:\n    {}", itertools::join(&polygon, ",\n    "));
        println!("Adding triangles:");

        // Create new triangles going from the edge to the new point.
        for edge in polygon {
            let to_add = Triangle(edge.from, edge.to, point);
            println!("    {}", to_add);
            triangulation.add_triangle(to_add);
        }

        println!("");
    }

    println!("Removing triangles touching super triangle.\n");
    // Remove triangles touching the super-triangle
    triangulation.remove_triangles_with_vertexes_of_triangle(super_triangle);

    triangulation
}

#[cfg(test)]
mod tests {
    use ::{Point, LineSegment, Line, Triangle, TriangleMesh, unique_elements_only, bower_watson_with_super_triangle};
    use ordered_float::NotNaN;

    fn generate_super_triangle() -> Triangle {
        Triangle(
            Point{x: NotNaN::new(-1f64).unwrap(), y: NotNaN::new(-1f64).unwrap()},
            Point{x: NotNaN::new(-1f64).unwrap(), y: NotNaN::new(3f64).unwrap()},
            Point{x: NotNaN::new(3f64).unwrap(), y: NotNaN::new(-1f64).unwrap()}
        )
    }

    #[test]
    fn unique_empty() {
        let empty: Vec<u32> = vec![];
        assert!(unique_elements_only(empty) == vec![]);
    }

    #[test]
    fn unique_single_element() {
        let single_element = vec![0u32];
        assert!(unique_elements_only(single_element) == vec![0u32]);
    }

    #[test]
    fn unique_double_element_distinct() {
        let double_element_distinct = vec![0u32, 1u32];
        assert!(unique_elements_only(double_element_distinct) == vec![0u32, 1u32]);
    }

    #[test]
    fn unique_double_element_same() {
        let double_element_same = vec![0u32, 0u32];
        assert!(unique_elements_only(double_element_same) == vec![]);
    }

    #[test]
    fn unique_triple_element_distinct() {
        let triple_element_distinct = vec![0u32, 1u32, 2u32];
        assert!(unique_elements_only(triple_element_distinct) == vec![0, 1, 2]);
    }

    #[test]
    fn unique_triple_element_same() {
        let triple_element_same = vec![0u32, 0, 0];
        assert!(unique_elements_only(triple_element_same) == vec![]);
    }

    #[test]
    fn unique_triple_element_two_same() {
        let triple_element_two_same = vec![0u32, 0, 1];
        assert!(unique_elements_only(triple_element_two_same) == vec![1]);
    }

    #[test]
    fn unique_edges_of_touching_triangles() {
        // Triangles
        // △{(-1, 4), (4, -1), (0, 0)}
        // △{(4, -1), (-1, -1), (0, 0)}

        // Edges
        // \{(-1, 4) -> (4, -1)}
        // \{(4, -1) -> (0, 0)}
        // \{(0, 0) -> (-1, 4)}
        // \{(4, -1) -> (-1, -1)}
        // \{(-1, -1) -> (0, 0)}
        // \{(0, 0) -> (4, -1)}

        let points = [
            Point::new_unwrap(0.0, 0.0),
            Point::new_unwrap(-1.0, -1.0),
            Point::new_unwrap(-1.0, 4.0),
            Point::new_unwrap(4.0, -1.0)
        ];

        let triangles = vec![
            Triangle(points[2], points[3], points[0]),
            Triangle(points[3], points[1], points[0])
        ];

        let edges = triangles.iter().flat_map(Triangle::edges).collect::<Vec<_>>();

        println!("{}", ::itertools::join(&edges, "\n"));

        let uniques = unique_elements_only(edges);

        println!("\n\n{}", ::itertools::join(&uniques, "\n"));

        assert!(true && false);
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

    #[test]
    fn bower_watson_of_no_points_is_empty_triangle_mesh() {
        assert!(bower_watson_with_super_triangle(&[], generate_super_triangle()) == TriangleMesh::default())
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
}
