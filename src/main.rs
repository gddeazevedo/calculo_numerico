#![allow(dead_code, unused)]

mod numeric_methods;
mod linear_systems;
mod types;

use linear_systems::exact_methods::{inverse};
use linear_systems::helpers::{matmat, print_matrix};



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
