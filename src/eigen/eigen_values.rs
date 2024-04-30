use crate::linear_systems::helpers::matmat;
use crate::types::Matrix;


pub fn tr(a: &Matrix<f64>) -> f64 {
    let mut sum = 0.0;

    for i in 0..a.len() {
        sum += a[i][i];
    }

    sum
}

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
