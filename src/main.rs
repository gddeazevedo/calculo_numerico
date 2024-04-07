#![allow(dead_code, unused)]

mod numeric_methods;
mod exact_methods;
mod types;

use exact_methods::linear_systems_solvers::{
    cholesky_solver, gaussian_elimination, gaussian_solver, inverse, lu_decomp
};
use types::Matrix;
use exact_methods::helpers::{ transpose, print_matrix, matmat };
use exact_methods::linear_systems_solvers::{lu_solver, solve_inf, solve_sup, lu_solver_solution_refinement};

use crate::exact_methods::helpers::vecnorm;



fn main() {
    let a = vec![
        vec![1.0, 3.0, 9.0],
        vec![1.0, 1.0, 3.0],
        vec![1.0, 1.0, 1.0],
    ];

    let inverse = inverse(&a);

    print_matrix(&inverse);
    print_matrix(&matmat(&a, &inverse));
}
