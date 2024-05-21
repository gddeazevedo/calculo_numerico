use std::vec;

use super::helpers::{divdot, infinity_vecnorm, matmat, matvec, scalar_mul, subvec, tr};
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
    let mut y = vec![1.0; n];
    let mut la = vec![0.0; n];
    let mut lp = vec![0.0; n];
    let epsilon = 1e-6;

    loop {
        let z = matvec(a, &y);
        lp = divdot(&z, &y);
        y = scalar_mul(&z, 1.0 / infinity_vecnorm(&z));

        if infinity_vecnorm(&subvec(&lp, &la)) / infinity_vecnorm(&lp) < epsilon {
            break;
        }
        la = lp;
    }
    
    lp[0]
}
