use std::vec;

use super::helpers::{divdot, infinity_vecnorm, lu_decomp, matmat, matvec, max_abs_value_in_inferior_triangle, scalar_mul, subvec, tr};
use crate::types::Matrix;


/**
 * The Leverrier method determines the coefficients of the characteristic polynomial of a matrix A
 */
pub fn leverrier(a: &Matrix<f64>) -> Vec<f64> {
    let n = a.len();
    let mut p = vec![0.0; n]; // stores the coefficients of the characteristic polynomial
    let mut s = vec![0.0; n];
    let mut a_cp = a.clone();

    for i in 0..n {
        s[i] = tr(&a_cp);
        a_cp = matmat(&a_cp, a);
    }

    for k in 1..=n {
        let mut sum = 0.0;
        
        for i in 1..k {
            sum += p[i-1] * s[k-i-1];
        }

        p[k-1] = (s[k-1] - sum) / (k as f64);
    }

    p
}


/**
 * The power method is an iterative method to find the largest eigenvalue of a matrix
 */
pub fn power_method(a: &Matrix<f64>) -> f64 {
    let n = a.len();
    let epsilon = 1e-100;

    let mut y = vec![1.; n]; // y[0]
    let mut z = matvec( a, &y ); // z[1]

    let mut l1 = z[0] / y[0];
    let mut laux = l1;

    loop {
        y = scalar_mul( &z, 1. / infinity_vecnorm( &z ) ); // y[k]
        z = matvec( a, &y ); // z[k + 1]

        l1 = z[0] / y[0];

        if f64::abs( l1 - laux ) / f64::abs( laux ) < epsilon {
            break;
        }

        laux = l1
    }

    l1
}


/**
 * Returns all the eigenvalues of a matrix using the LR algorithm
 */
pub fn rutishauser(a: &Matrix<f64>) -> Vec<f64> {
    let n = a.len();
    let mut a_ = a.clone();

    loop {
        let (l, r) = lu_decomp(&a_);
        a_ = matmat(&r, &l);

        if max_abs_value_in_inferior_triangle(&a_) < 1e-6 {
            break;
        }
    }

    let mut i: i64 = -1;

    a_.into_iter().map(|row| {
        i += 1;
        row[i as usize]
    }).collect()
}
