/**
 * Starts the numeric integration from the left side of the rectangle.
 * x[i] = a + i * delta_x
 */
pub fn integrate_left(f: fn(f64) -> f64, a: f64, b: f64, n: u32) -> f64
{
    let delta_x = ( b - a ) / n as f64;
    let mut x = a;
    let mut sum = 0.;

    for i in 0..n {
        sum += f( x );
        x += delta_x;
    }

    sum * delta_x
}


/**
 * Starts the numeric integration from the right side of the rectangle.
 * x[i] = a + (i + 1) * delta_x
 */
pub fn integrate_right(f: fn(f64) -> f64, a: f64, b: f64, n: u32) -> f64
{
    let delta_x = ( b - a ) / n as f64;
    let mut x = a + delta_x;
    let mut sum = 0.;

    for _ in 0..n {
        sum += f( x );
        x += delta_x;
    }

    sum * delta_x
}

/**
 * Starts the numeric integration from the middle of the rectangle.
 * x[i] = a + (i + 0.5) * delta_x
 */
pub fn integrate_middle(f: fn(f64) -> f64, a: f64, b: f64, n: u32) -> f64
{
    let delta_x = ( b - a ) / n as f64;
    let mut x = a + delta_x / 2.;
    let mut sum = 0.;

    for i in 0..n {
        sum += f( x );
        x += delta_x;
    }

    sum * delta_x
}

/**
 * Integrates the function using the trapezoid method.
 * x[i] = a + i * delta_x
 * Area = [f(x) + f(x + delta_x)] * delta_x / 2
 */
pub fn integrate_trapezoid(f: fn(f64) -> f64, a: f64, b: f64, n: u32) -> f64
{
    let delta_x = ( b - a ) / n as f64;
    let mut x = a;
    let mut sum = 0.;

    for _ in 0..n {
        sum += f( x ) + f( x + delta_x );
        x += delta_x;
    }

    sum * delta_x / 2.
}
