use crate::types::Matrix;
use super::helpers::{
    choose_best_pivot,
    transpose,
    matvec,
    addvec,
    subvec,
    print_matrix,
    vecnorm
};


/**
 * 
 * a must be a square matrix where a[i][j] == 0 for i < j
 * b is the vector of independent terms
 * this functions returns x, which is the variables vector, solution vector
 * starts with the first item and goes until the diagonal item
 */
pub fn solve_inf(a: &Matrix<f64>, b: &Vec<f64>) -> Vec<f64>
{
    let n = b.len();
    let mut x = vec![0.0; n];

    x[0] = b[0] / a[0][0];

    for i in 1..n {
        let mut sum = 0.0;

        for j in 0..i {
            sum += a[i][j] * x[j];
        }

        x[i] = (b[i] - sum) / a[i][i];
    }

    x
}


/**
 *
 * a must be a square matrix where a[i][j] == 0 for i > j
 * b is the vector of independent terms
 * this functions returns x, which is the variables vector, solution vector
 * starts with the diagonal item and goes until the last item
 */
pub fn solve_sup(a: &Matrix<f64>, b: &Vec<f64>) -> Vec<f64>
{
    let n = b.len();
    let mut x = vec![0.0; n];

    x[n - 1] = b[n - 1] / a[n - 1][n - 1];

    for i in (0..(n - 1)).rev() {
        let mut sum = 0.0;

        for j in (i + 1)..n {
            sum += a[i][j] * x[j];
        }

        x[i] = (b[i] - sum) / a[i][i];
    }

    x
}


/**
 * Decompose a square matrix A (n x n) into two square triangular matrices
 * Returns:
 *  - L (lower triangular matrix)
 *  - U (upper triangular matrix)
 */
pub fn lu_decomp(a: &Matrix<f64>) -> (Matrix<f64>, Matrix<f64>)
{
    let n = a.len();

    let mut l = vec![
        vec![0.0; n]; n
    ];

    let mut u = vec![
        vec![0.0; n]; n
    ];

    for i in 0..n {
        l[i][i] = 1.0;

        for j in 0..n {
            let mut sum = 0.0;
            
            if i <= j {
                for k in 0..i {
                    sum += l[i][k] * u[k][j];
                }

                u[i][j] = a[i][j] - sum;
            } else {
                for k in 0..j {
                    sum += l[i][k] * u[k][j];
                }

                l[i][j] = (a[i][j] - sum) / u[j][j];
            }
        }
    }

    (l, u)
}


/**
 * Solve a square linear system of order n using de decomposition LU
 * A = LU
 * LUx = b
 * Ly = b
 * Ux = y
 */
pub fn lu_solver(a: &Matrix<f64>, b: &Vec<f64>) -> Vec<f64>
{
    let (l, u) = lu_decomp(a);
    let y = solve_inf(&l, b);
    solve_sup(&u, &y) // x
}


/**
 * Gaussian elimination method to create an upper triangular matrix
 * Returns a superior triangular matrix and a vector that are equivalent to the
 * original ones passed as arguments
 * This function does not change the original vector and matrix
 */
pub fn gaussian_elimination(a: &Matrix<f64>, b: &Vec<f64>) -> (Matrix<f64>, Vec<f64>)
{
    let mut a_ = a.clone();
    let mut b_ = b.clone();
    let n = b.len();

    for k in 0..(n - 1) {
        for i in (k + 1)..n {
            let p = a_[i][k] / a_[k][k];
            b_[i] = b_[i] - b_[k] * p;

            for j in k..n {
                a_[i][j] = a_[i][j] - a_[k][j] * p;
            }
        }
    }

    (a_, b_)
}


/**
 * Gaussian elimination method to create an upper triangular matrix, using partial pivot
 */
pub fn partial_pivot_gaussian_elimination(a: &Matrix<f64>, b: &Vec<f64>) -> (Matrix<f64>, Vec<f64>)
{
    let mut a_ = a.clone();
    let mut b_ = b.clone();
    let n = b.len();

    for k in 0..(n - 1) {
        choose_best_pivot(&mut a_, &mut b_, k);

        for i in (k + 1)..n {
            let p = a_[i][k] / a_[k][k];
            b_[i] = b_[i] - b_[k] * p;

            for j in k..n {
                a_[i][j] = a_[i][j] - a_[k][j] * p;
            }
        }

    }

    (a_, b_)
}


/**
 * Solves linear systems of order n using the gaussian elimination
 * partial_pivot controls whether the algorithm is going to use partial pivoting or not
 */
pub fn gaussian_solver(a: &Matrix<f64>, b: &Vec<f64>, partial_pivot: bool) -> Vec<f64>
{   
    let (a_, b_) = if partial_pivot {
        gaussian_elimination(a, b)
    } else {
        partial_pivot_gaussian_elimination(a, b)
    };

    solve_sup(&a_, &b_)
}


/**
 * 
 */
// pub fn gaussian_compact(a: &mut Matrix<f64>) -> Vec<f64>
// {
//     let n = a.len();
//     let mut m = vec![vec![0.0; n]; n];

//     m[0] = a[0].clone();

//     for i in 0..n {
//         for j in 0..n {
//             if i <= j {

//             } else {

//             }
//         }
//     }

//     m
// }


/**
 * Returns two matrices, G (lower triangular) and G transpose (upper triangular)
 * The cholesky method only works for matrices that are simetric and positive definite
 */
pub fn cholesky_method(a: &Matrix<f64>) -> (Matrix<f64>, Matrix<f64>)
{
    let n = a.len();

    let mut g = vec![
        vec![0.0; n]; n
    ];

    for i in 0..n {
        for j in 0..n {
            let mut s = 0.0;

            if i == j {
                for k in 0..i {
                    s += g[i][k] * g[i][k];
                }

                g[i][i] = f64::sqrt(a[i][i] - s);
            } else if i > j {
                for k in 0..j {
                    s += g[i][k] * g[j][k];
                }

                g[i][j] = (a[i][j] - s) / g[j][j];
            }
        }
    }

    let gt = transpose(&g);

    (g, gt)
}


pub fn cholesky_solver(a: &Matrix<f64>, b: &Vec<f64>) -> Vec<f64>
{
    let (g, gt) = cholesky_method(a);
    let y = solve_inf(&g, b);
    solve_sup(&gt, &y)
}

pub fn lu_solver_solution_refinement(a: &Matrix<f64>, b: &Vec<f64>) -> Vec<f64>
{
    let n = a.len();
    let mut r = b.clone();
    let mut x = vec![0.0; n];
    let epsilon = 1e-10;

    loop {
        let y = lu_solver(a, &r);
        x = addvec(&x, &y);
        r = subvec(b, &matvec(a, &x));

        if vecnorm(&r) < epsilon {
            break;
        }
    }

    x
}
