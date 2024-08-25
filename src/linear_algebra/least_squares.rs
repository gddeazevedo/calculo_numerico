use crate::types::Matrix;
use super::{
    exact_methods::lu_solver,
    helpers::{matmat, matvec, mean, transpose}
};

/**
 * Finds the a and b that fits a line as close as possible to the points (x[i], y[i])
 */
pub fn linear_regression(x: &Vec<f64>, y: &Vec<f64>) -> (f64, f64) {
    if ( x.len() != y.len() ) {
        panic!("The vectors x and y must have the same length");
    }

    let n = x.len();
    let mean_x = mean(x);
    let mean_y = mean(y);

    let mut numerator   = 0.;
    let mut denominator = 0.;

    for i in 0..n {
        numerator   += x[i] * (y[i] - mean_y);
        denominator += x[i] * (x[i] - mean_x);
    }

    let a = numerator / denominator;
    let b = mean_y - a * mean_x;

    (a, b)
}

/**
 * Finds the vector beta that better fits the line to the points (x[1][i], x[2][i], ..., x[k][i], y[i]); i = 1, 2, ..., n
 */
pub fn linear_multiple_regression(x: &Matrix<f64>, y: &Vec<f64>, n: u8) -> Vec<f64> {
    let x_t   = transpose(x);
    let x_t_x = matmat( &x_t, x );
    let x_t_y = matvec( &x_t, y );

    lu_solver(&x_t_x, &x_t_y)
}
