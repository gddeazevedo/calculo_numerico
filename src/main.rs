#![allow(dead_code, unused)]

mod numeric_methods;
mod functions;
mod exact_methods;

use exact_methods::linear_systems_solvers::{solve_inf, solve_sup, matvec, matmat};


fn main() {

    let a = vec![
        vec![1.0, -1.0, 2.0],
        vec![0.0,  3.0, 1.0],
        vec![0.0,  0.0, 1.0]
    ];

    let b = vec![1.0, 3.0, 6.0];

    let x = solve_sup(&a, &b);

    let b1 = matvec(&a, &x);

    println!("x = {:?}", x);
    println!("b = {:?}", b);
    println!("A * x = {:?}", b1);

    println!("{:?}", matmat(&a, &a));
}
