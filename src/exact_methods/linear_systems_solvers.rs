/**
 * 
 * a must be a square matrix where a[i][j] == 0 for i < j
 * b is the vector of independent terms
 * this functions returns x, which is the variables vector, solution vector
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
