#![allow(dead_code, unused)]

use linear_systems::{exact_methods::cholesky_solver, helpers::{cholesky_method, print_matrix, transpose}, iterative_methods::{gauss_seidel_solver, jacobi_richardson_solver}};
use numeric_methods::function_roots::{linear_iteration_method, regula_falsi};

use crate::numeric_methods::function_roots::newton_method;

mod numeric_methods;
mod linear_systems;
mod types;


fn main() {
    let f = |x: f64| f64::exp(1.0833 * x) - 5.*x*x;
    let x = newton_method(f, 3.5);
    println!("x = {x}");
    println!("f({x}) = {}", f(x));
}
