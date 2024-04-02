#![allow(dead_code, unused)]

mod numeric_methods;
mod exact_methods;
mod types;

use exact_methods::linear_systems_solvers::{
    cholesky_solver, gaussian_elimination, gaussian_solver, lu_decomp
};
use types::Matrix;
use exact_methods::helpers::{ transpose, print_matrix, matmat };
use exact_methods::linear_systems_solvers::{lu_solver, solve_inf, solve_sup, lu_solver_solution_refinement};

use crate::exact_methods::helpers::vecnorm;



fn main() {
    let mut a = vec![
        vec![16.0, 5.0],
        vec![3.0,  2.5]
    ];

    let b = vec![21.0, 5.5];
    let x = lu_solver_solution_refinement(&a, &b);

    println!("{:?}", x);
}
