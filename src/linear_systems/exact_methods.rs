use crate::types::Matrix;
use super::helpers::{
    addvec,
    choose_best_pivot,
    infinity_norm, matvec,
    print_matrix,
    subvec,
    transpose,
    vecnorm,
    cholesky_method,
    lu_decomp,
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



fn lu_solver_refine(l: &Matrix<f64>, u: &Matrix<f64>, b: &Vec<f64>) -> Vec<f64> {
    let y = solve_inf(l, b);
    solve_sup(u, &y)
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
pub fn gaussian_compact(a: &mut Matrix<f64>, b: &mut Vec<f64>)
{
    let n = a.len();

    for k in 0..(n-1) {
        for i in (k + 1)..n {
            let p = a[i][k] / a[k][k];
            b[i] = b[i] - b[k] * p;

            for j in k..n {
                if i > j {
                    a[i][j] = p;
                } else {
                    a[i][j] = a[i][j] - a[k][j] * p;
                }
            }
        }
    }
}


pub fn cholesky_solver(a: &Matrix<f64>, b: &Vec<f64>) -> Vec<f64>
{
    let (g, gt) = cholesky_method(a);
    let y = solve_inf(&g, b);
    solve_sup(&gt, &y)
}


pub fn lu_solver_solution_refinement(a: &Matrix<f64>, b: &Vec<f64>) -> Vec<f64>
{
    let epsilon = 1e-10;
    let n = a.len();

    let (l, u) = lu_decomp(a);

    let mut x = lu_solver_refine(&l, &u, b);
    let mut r = subvec(b, &matvec(a, &x)); // r = b - Ax

    while vecnorm(&r) >= epsilon {
        let y = lu_solver_refine(&l, &u, b); // Ay = r
        x = addvec(&x, &y);
        r = subvec(b, &matvec(a, &x));
    }

    x
}


pub fn inverse(a: &Matrix<f64>) -> Matrix<f64> {
    let n = a.len();

    let mut inverse_t: Matrix<f64> = vec![]; // tranpose of the inverse matrix, its rows are the columns of the inverse

    let (l, u) = lu_decomp(a);

    for i in 0..n {
        let mut e: Vec<f64> = vec![0.0; n];

        for j in 0..n {
            e[j] = if i == j { 1.0 } else { 0.0 }; // column  j of identity matrix
        }

        let col = lu_solver_refine(&l, &u, &e);

        inverse_t.push(col);
    }

    transpose(&inverse_t)
}


pub fn cond(a: &Matrix<f64>) -> f64 {
    let a_inverse = inverse(a);
    infinity_norm(a) * infinity_norm(&a_inverse)
}
