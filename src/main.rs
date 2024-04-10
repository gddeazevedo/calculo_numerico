#![allow(dead_code, unused)]

use linear_systems::{exact_methods::gaussian_compact, iterative_methods::jacobi_richardson_solver};

use crate::linear_systems::{helpers::print_matrix, iterative_methods::gauss_seidel_solver};

mod numeric_methods;
mod linear_systems;
mod types;


fn main() {
    let mut a = vec![
        vec![5.0, 1.0, 1.0],
        vec![3.0, 4.0, 1.0],
        vec![3.0, 3.0, 6.0],
    ];

    let mut b = vec![5.0, 4.0, 6.0];

    gaussian_compact(&mut a, &mut b);

    println!("{:?}\n", b);
    print_matrix(&a);
}
