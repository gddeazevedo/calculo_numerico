#![allow(dead_code, unused)]

mod numeric_methods;
mod functions;

#[macro_use]
use numeric_methods::function_roots::{
    bisec_method,
    // secant_method,
    sqrt,
    linear_iteration_method,
    newton_method
};


fn f(x: f64) -> f64
{
    f64::exp(x) + f64::sin(x*x) - 10.0
}

fn f1(x: f64) -> f64
{
    4.0 * f64::cos(x) - f64::exp(x)
}

fn main() {
    let x_bisec = bisec_method(f, 0.0, 10.0, 1e-12).unwrap_or(0.0);
    println!("Bissection: x = {} f(x) = {}", x_bisec, f(x_bisec));

    // let x_secant = secant_method(f, 0.0, 10.0, 1e-11).unwrap_or(0.0);
    // println!("Secant: x = {} f(x) = {}", x_secant, f(x_secant));

    println!("sqrt(144) = {}", sqrt(144.0).unwrap());

    println!("Root of f(x) = x*x - x - 2 is {}",
        linear_iteration_method(|x| f64::sqrt(2.0 + x), 2.5));

    let x_newton = newton_method(f1, 1.0);
    println!("Newton: x = {} f(x) = {}", x_newton, f1(x_newton));
}
