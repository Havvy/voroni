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

mod geometry;

use std::collections::{HashSet};
use std::hash::{Hash};

pub use ordered_float::{FloatIsNaN, NotNaN};

pub use ::geometry::{Point, Triangle, TriangleMesh};
// use geometry::{Triangle};

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
    use ::geometry::{Point, Triangle, TriangleMesh};
    use ::{unique_elements_only, bower_watson_with_super_triangle};
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

    // #[test]
    // fn unique_edges_of_touching_triangles() {
    //     // Triangles
    //     // △{(-1, 4), (4, -1), (0, 0)}
    //     // △{(4, -1), (-1, -1), (0, 0)}

    //     // Edges
    //     // \{(-1, 4) -> (4, -1)}
    //     // \{(4, -1) -> (0, 0)}
    //     // \{(0, 0) -> (-1, 4)}
    //     // \{(4, -1) -> (-1, -1)}
    //     // \{(-1, -1) -> (0, 0)}
    //     // \{(0, 0) -> (4, -1)}

    //     let points = [
    //         Point::new_unwrap(0.0, 0.0),
    //         Point::new_unwrap(-1.0, -1.0),
    //         Point::new_unwrap(-1.0, 4.0),
    //         Point::new_unwrap(4.0, -1.0)
    //     ];

    //     let triangles = vec![
    //         Triangle(points[2], points[3], points[0]),
    //         Triangle(points[3], points[1], points[0])
    //     ];

    //     let edges = triangles.iter().flat_map(Triangle::edges).collect::<Vec<_>>();

    //     println!("{}", ::itertools::join(&edges, "\n"));

    //     let uniques = unique_elements_only(edges);

    //     println!("\n\n{}", ::itertools::join(&uniques, "\n"));

    //     assert!(true && false);
    // }

    #[test]
    fn bower_watson_of_no_points_is_empty_triangle_mesh() {
        assert!(bower_watson_with_super_triangle(&[], generate_super_triangle()) == TriangleMesh::default())
    }
}
