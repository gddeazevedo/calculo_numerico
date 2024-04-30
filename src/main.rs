#![allow(dead_code, unused)]

use eigen::eigen_values::{leverrier, tr};
use linear_systems::{exact_methods::cholesky_solver, helpers::{cholesky_method, print_matrix, transpose}, iterative_methods::{gauss_seidel_solver, jacobi_richardson_solver}};
use numeric_methods::function_roots::{linear_iteration_method, regula_falsi};

use crate::{linear_systems::helpers::matmat, numeric_methods::function_roots::newton_method};

mod numeric_methods;
mod linear_systems;
mod eigen;
mod types;


fn main() {
    let a = vec![
        vec![1., 1., -1.],
        vec![0., 0., 1.],
        vec![-1., 1., 0.]
    ];
    let p = leverrier(&a);
    println!("{:?}", p);
}
