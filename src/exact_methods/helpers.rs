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


pub fn vecnorm(v: &Vec<f64>) -> f64 {
    let mut norm = 0.0;

    for i in v {
        norm += i * i;
    }

    f64::sqrt(norm)
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


pub fn print_matrix(a: &Matrix<f64>) {
    for line in a {
        println!("{:?}", line);
    }
}
