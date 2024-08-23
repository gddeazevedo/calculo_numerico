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

    let mut numerator = 0.;
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
 * Returns the mean of a collection of items
 */
fn mean(x: &Vec<f64>) -> f64 {
    x.iter().sum::<f64>() / x.len() as f64
}
