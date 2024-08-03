#![allow(unused)]

mod calculus;
mod linear_algebra;
mod types;


use calculus::integration::{
    integrate_left, integrate_middle, integrate_right, integrate_trapezoid
};


fn main() {
    let a = 0.;
    let b = 1.;
    let n = 1000;

    let result = integrate_left( f, a, b, n );
    println!( "ESQ: {}", result );

    let result = integrate_right( f, a, b, n );
    println!( "DIR: {}", result );

    let result = integrate_middle( f, a, b, n );
    println!( "CEN: {}", result );

    let result = integrate_trapezoid(f, a, b, n);
    println!( "TRA: {}", result );

    println!( "PI : {}", std::f64::consts::PI ); 
}


fn f(x: f64) -> f64 {
    4. / (x*x + 1.)
}
