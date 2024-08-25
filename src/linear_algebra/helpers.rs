use crate::types::Matrix;

/**
 * Returns the vector that results from the product between
 * a matrix n x n and a vector of dimension n
 */
pub fn matvec(a: &Matrix<f64>, v: &Vec<f64>) -> Vec<f64>
{
    let mut x: Vec<f64> = vec![0.0; v.len()];
    let n = v.len();

    for i in 0..n {
        let mut sum = 0.0;

        for j in 0..n {
            sum += a[i][j] * v[j];
        }

        x[i] = sum; 
    }

    x
}


/**
 * Return the product of two matrices
 */
pub fn matmat(a: &Matrix<f64>, b: &Matrix<f64>) -> Matrix<f64>
{
    let n = a.len();
    let mut c = vec![vec![0.0; n]; n];

    for i in 0..n {
        
        for j in 0..n {
            let mut sum = 0.0;
            
            for k in 0..n {
                sum += a[i][k] * b[k][j];
            }

            c[i][j] = sum;
        }
    }

    c
}


pub fn matsum(a: &Matrix<f64>, b: &Matrix<f64>) -> Matrix<f64> {
    let n = a.len();
    let mut c = vec![vec![0.0; n]; n];

    for i in 0..n {
        for j in 0..n {
            c[i][j] = a[i][j] + b[i][j];
        }
    }

    c
}


/**
 * Returns the sum of two vector
 */
pub fn addvec(u: &Vec<f64>, v: &Vec<f64>) -> Vec<f64>
{
    let n = u.len();
    let mut x = vec![0.0; n];

    for i in 0..n {
        x[i] = u[i] + v[i];
    }

    x
}


/**
 * Returns the subtraction of two vectors
 */
pub fn subvec(u: &Vec<f64>, v: &Vec<f64>) -> Vec<f64>
{
    let n = u.len();
    let mut x = vec![0.0; n];

    for i in 0..n {
        x[i] = u[i] - v[i];
    }

    x
}


/**
 * Returns the Euclidean norm of a vector
 */
pub fn vecnorm(v: &Vec<f64>) -> f64 {
    let norm: f64 = v.into_iter().map(|component| component * component).sum();
    f64::sqrt(norm)
}


/**
 * Returns the division of two vectors, component by component
 */
pub fn divdot(u: &Vec<f64>, v: &Vec<f64>) -> Vec<f64> {
    let n = u.len();
    let mut x = vec![0.0; n];

    for i in 0..n {
        x[i] = u[i] / v[i];
    }

    x
}


/**
 * Returns the infinity norm of a vector
 */
pub fn infinity_vecnorm(v: &Vec<f64>) -> f64 {
    let abs = |el: &f64| f64::abs(*el);
    v.into_iter().map(abs).fold(f64::NEG_INFINITY, f64::max)
}


/**
 * Returns the scalar multiplication of a vector
 */
pub fn scalar_mul(v: &Vec<f64>, scalar: f64) -> Vec<f64> {
    v.iter().map(|el| el * scalar).collect()
}


/**
 * Returns the infinity norm of a n x n matrix
 * ||A|| = max, 1 <= i<= n ( sum(a\[i\]\[j\]), 1 <= j <= n )
 */
pub fn infinity_norm(a: &Matrix<f64>) -> f64 {
    let n = a.len();
    let abs = |el: &f64| f64::abs(*el);

    let mut max: f64 = (&a[0]).into_iter().map(abs).sum();

    for i in 1..n {
        let sum = (&a[i]).into_iter().map(abs).sum();

        if sum > max {
            max = sum;
        }
    }

    max
}


/**
 * Returns the trace of a matrix
 */
pub fn tr(a: &Matrix<f64>) -> f64 {
    let mut sum = 0.0;

    for i in 0..a.len() {
        sum += a[i][i];
    }

    sum
}


/**
 * Returns the maximum absolute value in the inferior triangle of a matrix
 */
pub fn max_abs_value_in_inferior_triangle(a: &Matrix<f64>) -> f64 {
    let n = a.len();
    let mut max = f64::NEG_INFINITY;

    for i in 1..n {
        for j in 0..i {
            max = f64::max(max, f64::abs(a[i][j]));
        }
    }

    max
}


/**
 * Returns a n x n identity matrix
 */
pub fn get_identity_matrix(n: usize) -> Matrix<f64>
{
    let mut j = 0;

    let i = vec![
        vec![0.0; n]; n
    ]
        .into_iter()
        .map(|mut line| {
            line[j] = 1.0;
            j += 1;
            line
        })
        .collect();

    i
}


/**
 * Chooses the element with greatest absolute value as pivot of a column
 */
pub fn choose_best_pivot(a: &mut Matrix<f64>, b: &mut Vec<f64>, k: usize) {
    let n = b.len();
    let mut max_element = f64::abs(a[k][k]);
    let mut max_index = k;

    for line in (k + 1)..n {
        if f64::abs(a[line][k]) > max_element {
            max_element = f64::abs(a[line][k]);
            max_index   = line;
        }
    }

    if k != max_index {
        // swaps the line k with the line max_index which contains the max pivot

        let tmp = a[k].clone();
        a[k] = a[max_index].clone();
        a[max_index] = tmp;

        let tmp = b[k];
        b[k] = b[max_index];
        b[max_index] = tmp;
    }
}


/**
 * Returns the transpose of matrix a
 */
pub fn transpose(a: &Matrix<f64>) -> Matrix<f64>
{
    let lines = a.len();
    let cols = a[0].len();

    let mut at = vec![
        vec![0.0; lines]; cols
    ];

    for i in 0..cols {
        for j in 0..lines {
            at[i][j] = a[j][i];
        }
    }

    at
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


/**
 * Decompose a square matrix A in two matrices L* (upper) and R* (lower) and the vector b*
 * where A* = L* + I + R*
 * A* is the matrix A with each row divided by the correspondent main diagonal element
 */
pub fn lrb_star_decomp(a: &Matrix<f64>, b: &Vec<f64>) -> (Matrix<f64>, Matrix<f64>, Vec<f64>)
{
    let n = a.len();
    let mut l_star: Matrix<f64> = vec![vec![0.0; n]; n];
    let mut r_star: Matrix<f64> = vec![vec![0.0; n]; n];
    let mut b_star = vec![0.0; n];

    for i in 0..n {
        for j in 0..n {
            l_star[i][j] = if i > j { a[i][j] / a[i][i] } else { 0.0 };
            r_star[i][j] = if i < j { a[i][j] / a[i][i] } else { 0.0 };
        }

        b_star[i] = b[i] / a[i][i];
    }

    (l_star, r_star, b_star)
}


/**
 * Decompose a square matrix A in two matrices L* (upper) and R* (lower)
 * where A* = L* + I + R*
 * A* is the matrix A with each row divided by the correspondent main diagonal element
 */
pub fn lr_star_decomp(a: &Matrix<f64>) -> (Matrix<f64>, Matrix<f64>)
{
    let n = a.len();
    let mut l_star: Matrix<f64> = vec![vec![0.0; n]; n];
    let mut r_star: Matrix<f64> = vec![vec![0.0; n]; n];

    for i in 0..n {
        for j in 0..n {
            l_star[i][j] = if i > j { a[i][j] / a[i][i] } else { 0.0 };
            r_star[i][j] = if i < j { a[i][j] / a[i][i] } else { 0.0 };
        }
    }

    (l_star, r_star)
}


/**
 * Prints a matrix
 */
pub fn print_matrix(a: &Matrix<f64>) {
    for line in a {
        println!("{:?}", line);
    }
}

/**
 * Returns the mean of a collection of items
 */
pub fn mean(x: &Vec<f64>) -> f64 {
    x.iter().sum::<f64>() / x.len() as f64
}

