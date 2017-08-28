extern crate voroni;

use voroni::NotNaN;
use voroni::geometry::{Point, Triangle};
use voroni::delaunay::bower_watson_with_super_triangle_with_logger as bower_watson;
use voroni::delaunay::PrintlnLogger;

fn generate_corners() -> [Point; 4] {
    [
        Point::new_unwrap(0.0, 0.0),
        Point::new_unwrap(1.0, 0.0),
        Point::new_unwrap(1.0, 1.0),
        Point::new_unwrap(0.0, 1.0)
    ]
}

fn generate_super_triangle() -> Triangle {
    Triangle(
        Point::new_unwrap(-1.0, -1.0),
        Point::new_unwrap(-1.0, 4.0),
        Point::new_unwrap(4.0, -1.0)
    )
}

fn generate_grid(points: u8) -> Vec<Point> {
    let mut grid = Vec::with_capacity((points * points + 4) as usize);

    for x in 0..points {
        for y in 0..points {
            grid.push(Point::new_unwrap((x + 1) as f64 / (points + 1) as f64, (y + 1) as f64 / (points + 1) as f64));
        }
    }

    return grid;
}

fn main () {
    let points = {
        let mut points = generate_grid(3);
        for point in generate_corners().into_iter() {
            points.push(*point);
        }
        points
    };
    let out = bower_watson(&points, generate_super_triangle(), PrintlnLogger);
    println!("{}", out * &NotNaN::new(4f64).unwrap());
    println!("End of program.");
}