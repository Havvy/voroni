extern crate voroni;
extern crate ordered_float;

use voroni::{Point, Triangle, bower_watson_with_super_triangle as bower_watson, NotNaN};

/// Create a Point. Panics if given NaN.
fn point(x: f64, y: f64) -> Point {
    Point{x: NotNaN::new(x).unwrap(), y: NotNaN::new(y).unwrap()}
}

fn generate_corners() -> [Point; 4] {
    [
        point(0.0, 0.0),
        point(1.0, 0.0),
        point(1.0, 1.0),
        point(0.0, 1.0)
    ]
}

fn generate_super_triangle() -> Triangle {
    Triangle(
        point(-1.0, -1.0),
        point(-1.0, 4.0),
        point(4.0, -1.0)
    )
}

fn main () {
    let out = bower_watson(&generate_corners(), generate_super_triangle());
    println!("{}", out);
    println!("End of program.");
}