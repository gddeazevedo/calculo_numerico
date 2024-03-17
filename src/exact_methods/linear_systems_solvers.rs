/**
 * 
 * a must be a square matrix where a[i][j] == 0 for i < j
 * b is the vector of independent terms
 * this functions returns x, which is the variables vector, solution vector
 * starts with the first item and goes until the diagonal item
 */
pub fn solve_inf(a: &Vec<Vec<f64>>, b: &Vec<f64>) -> Vec<f64>
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
pub fn solve_sup(a: &Vec<Vec<f64>>, b: &Vec<f64>) -> Vec<f64>
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
 * Returns the vector that results from the product between
 * a matrix n x n and a vector of dimension n
 */
pub fn matvec(a: &Vec<Vec<f64>>, v: &Vec<f64>) -> Vec<f64>
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
pub fn matmat(a: &Vec<Vec<f64>>, b: &Vec<Vec<f64>>) -> Vec<Vec<f64>>
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
 * Decompose a square matrix A (n x n) into two square triangular matrices
 * L (lower triangular matrix)
 * U (upper triangular matrix)
 */
pub fn lu_decomp(a: &Vec<Vec<f64>>) -> (Vec<Vec<f64>>, Vec<Vec<f64>>)
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
 * 
 */
pub fn lu_solver(a: &Vec<Vec<f64>>, b: &Vec<f64>) -> Vec<f64>
{
    let (l, u) = lu_decomp(a);
    let y = solve_inf(&l, b);
    solve_sup(&u, &y)
}
