use crate::types::Matrix;
use super::helpers::{lrb_star_decomp, matsum, matvec, subvec, vecnorm};


pub fn jacobi_richardson_solver(a: &Matrix<f64>, b: &Vec<f64>) -> Vec<f64> {
    let n = a.len();
    let epsilon = 1e-17;

    let (l_star, r_star, b_star) = lrb_star_decomp(a, b);
    let lr = matsum(&l_star, &r_star);

    let mut x0 = vec![0.0; n];
    let mut x1 = vec![0.0; n];

    loop {
        x1 = subvec(&b_star, &matvec(&lr, &x0));

        if vecnorm(&subvec(&x1, &x0)) / vecnorm(&x1) < epsilon {
            break;
        }

        x0 = x1.clone();
    }

    x1
}


pub fn gauss_seidel_solver(a: &Matrix<f64>, b: &Vec<f64>) -> Vec<f64> {
    let n = a.len();
    let epsilon = 1e-17;

    let (l_star, r_star, b_star) = lrb_star_decomp(a, b);

    let mut x0 = vec![0.0; n];
    let mut x1 = vec![0.0; n];

    loop {
        for i in 0..n {
            let mut r = 0.0;

            // lower triangular matrix i > j
            for j in 0..i {
                r -= l_star[i][j] * x1[j];
            }

            // upper triangular matrix i < j
            for j in (i + 1)..n {
                r -= r_star[i][j] * x0[j];
            }

            r += b_star[i];
            x1[i] = r;
        }

        if vecnorm(&subvec(&x1, &x0)) / vecnorm(&x1) < epsilon {
            break;
        }

        x0 = x1.clone();
    }

    x1
}
