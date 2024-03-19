#![allow(dead_code, unused)]

mod numeric_methods;
mod exact_methods;
mod types;

use exact_methods::linear_systems_solvers::{
    cholesky_solver, matvec,
};

use exact_methods::helpers::{ transpose, print_matrix };

use crate::exact_methods::linear_systems_solvers::{gaussian_solver, lu_solver};


fn main() {
    let a = vec![
        vec![4.0, 2.0, -4.0],
        vec![2.0, 10.0, 4.0],
        vec![-4.0, 4.0, 9.0]
    ];

    let b = vec![0.0, 6.0, 5.0];
    let x = cholesky_solver(&a, &b);
    println!("{:?}", x);
}
