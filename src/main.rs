#![allow(unused)]

use crate::linear_algebra::helpers::print_matrix;

mod numeric_methods;
mod linear_algebra;
mod types;


fn main() {
    let a = vec![
        vec![2., 0., 1.],
        vec![0., 1., 0.],
        vec![1., 0., 1.]
    ];

    let eigenvalues = linear_algebra::eigen::rutishauser(&a);

    println!("Eigenvalues:");
    for eigenvalue in &eigenvalues {
        println!("{}", eigenvalue);
    }

    print_matrix(&a);
}
