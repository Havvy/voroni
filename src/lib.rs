extern crate itertools;
extern crate ordered_float;
extern crate half_edge;

pub mod geometry;
pub mod delaunay;
pub mod voroni;

use std::collections::{HashSet};
use std::hash::{Hash};

pub use ordered_float::{FloatIsNaN, NotNaN};

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

#[cfg(test)]
mod tests {
    use ::geometry::{Point, Triangle, TriangleMesh};
    use ::unique_elements_only;
    use :: delaunay::bower_watson_with_super_triangle;
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
    fn bower_watson_of_no_points_is_empty_triangle_mesh() {
        assert!(bower_watson_with_super_triangle(&[], generate_super_triangle()) == TriangleMesh::default())
    }
}