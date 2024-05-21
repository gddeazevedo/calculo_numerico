#![allow(dead_code, unused)]


use linear_systems::{eigen::power_method, exact_methods::cholesky_solver, helpers::{cholesky_method, infinity_vecnorm, print_matrix, transpose}, iterative_methods::{gauss_seidel_solver, jacobi_richardson_solver}};
use numeric_methods::function_roots::{linear_iteration_method, regula_falsi};

use linear_systems::eigen::leverrier;

use crate::linear_systems::helpers::scalar_mul;

mod numeric_methods;
mod linear_systems;
mod types;


fn main() {
    let a = vec![
        vec![3., 0., 1.0],
        vec![2., 2., 2.],
        vec![4., 2., 5.]
    ];

    let l = power_method(&a);
    println!("Largest eigenvalue: {}", l);
}
