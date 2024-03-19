#![allow(dead_code, unused)]

mod numeric_methods;
mod helpers;
mod exact_methods;
mod types;

use exact_methods::linear_systems_solvers::{
    gaussian_solver, lu_decomp, lu_solver, matmat, matvec, solve_inf, solve_sup
};

use crate::helpers::get_identity_matrix;


fn main() {

    // let a = vec![
    //     vec![1.0, -1.0, 2.0],
    //     vec![0.0,  3.0, 1.0],
    //     vec![0.0,  0.0, 1.0]
    // ];

    // let b = vec![1.0, 3.0, 6.0];

    // let x = solve_sup(&a, &b);

    // let b1 = matvec(&a, &x);

    // println!("x = {:?}", x);
    // println!("b = {:?}", b);
    // println!("A * x = {:?}", b1);

    // println!("{:?}", matmat(&a, &a));

    let a = vec![
        vec![5.0, 2.0, 1.0],
        vec![3.0, 1.0, 4.0],
        vec![1.0, 1.0, 3.0],
    ];

    let b = vec![0.0, -7.0, -5.0];
    let x = gaussian_solver(&a, &b, true);

    println!("x = {:?} (Gaussian solver partial pivot)", x);

    let x = lu_solver(&a, &b);
    println!("x = {:?} (LU solver)", x);

    let a = vec![
        vec![6.0, 2.0, -1.0],
        vec![2.0, 4.0,  1.0],
        vec![3.0, 2.0,  8.0],
    ];


    let b = vec![7.0, 7.0, 13.0];

    // let x: Vec<i32> = gaussian_solver(&a, &b, false)
    //     .into_iter()
    //     .map(|xi| f64::round(xi) as i32)
    //     .collect();

    let x = gaussian_solver(&a, &b, false);
    println!("x = {:?} (Gaussian solver no partial pivot)", x);

    let x = gaussian_solver(&a, &b, true);
    println!("x = {:?} (Gaussian solver partial pivot)", x);
}
