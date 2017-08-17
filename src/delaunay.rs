use ::itertools;

use ::geometry::{Point, Triangle, TriangleMesh};

/// Translate a list of cordinates into the Delaunay triangulation via the Boyer-Watson algorithm.
///
/// https://en.wikipedia.org/wiki/Delaunay_triangulation
/// https://en.wikipedia.org/wiki/Bowyer%E2%80%93Watson_algorithm
pub fn bower_watson (point_list: &[Point]) -> TriangleMesh {
    let super_triangle = Triangle::super_triangle_of(point_list);

    bower_watson_with_super_triangle_with_logger(point_list, super_triangle, DiscardLogger)
}

pub fn bower_watson_with_logger<L>(point_list: &[Point], logger: L) -> TriangleMesh where L : BowerWatsonLogger {
    let super_triangle = Triangle::super_triangle_of(point_list);

    bower_watson_with_super_triangle_with_logger(point_list, super_triangle, logger)
}

pub fn bower_watson_with_super_triangle(point_list: &[Point], super_triangle: Triangle) -> TriangleMesh {
    bower_watson_with_super_triangle_with_logger(point_list, super_triangle, DiscardLogger)
}

pub fn bower_watson_with_super_triangle_with_logger<L>(point_list: &[Point], super_triangle: Triangle, mut logger: L) -> TriangleMesh where L: BowerWatsonLogger {
    // Start with an empty Triangle Mesh.
    let mut triangulation = TriangleMesh::default();

    // Add a Triangle that contains all points in the point list.
    triangulation.add_triangle(super_triangle);

    logger.log(BowerWatsonLogMessage::AddSuperTriangle(&triangulation));

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
        let polygon = ::unique_elements_only(bad_triangles.iter().flat_map(Triangle::edges).collect());

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

#[derive(Debug)]
pub enum BowerWatsonLogMessage<'a> {
    AddSuperTriangle(&'a TriangleMesh),
    AddPoint(&'a Point),
    ValidTriangleTestBegin,
    ValidTriangleTest(&'a Triangle, bool)
}

impl<'a> ::std::fmt::Display for BowerWatsonLogMessage<'a> {
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> Result<(), ::std::fmt::Error> {
        match *self {
            BowerWatsonLogMessage::AddSuperTriangle(ref triangulation) => {
                write!(f, "Added super triangle.\n{}\n", triangulation)
            }
            _ => { write!(f, "") }
        }
    }
}

/// Logger trait for the BowerWatson algorithm.
pub trait BowerWatsonLogger {
    /// The function that gets a BowerWatsonLogMessage and does the actual logging.
    fn log(&mut self, message: BowerWatsonLogMessage);
}

#[derive(Debug)]
pub struct DiscardLogger;

impl BowerWatsonLogger for DiscardLogger {
    fn log(&mut self, _message: BowerWatsonLogMessage) {
        ()
    }
}

#[derive(Debug)]
pub struct PrintlnLogger;

impl BowerWatsonLogger for PrintlnLogger {
    fn log(&mut self, message: BowerWatsonLogMessage) {
        println!("{}", message);
    }
}