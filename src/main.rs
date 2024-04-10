#![allow(dead_code, unused)]

use linear_systems::iterative_methods::jacobi_richardson_solver;

use crate::linear_systems::iterative_methods::gauss_seidel_solver;

mod numeric_methods;
mod linear_systems;
mod types;





fn main() {
    let a = vec![
        vec![5.0, 1.0, 1.0],
        vec![3.0, 4.0, 1.0],
        vec![3.0, 3.0, 6.0],
    ];

    let b = vec![5.0, 6.0, 0.0];

    let x = jacobi_richardson_solver(&a, &b);
    println!("{:?}", x);
}
