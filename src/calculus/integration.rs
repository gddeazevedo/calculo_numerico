/**
 * Starts the numeric integration from a rectangle that starts from the left side of the point.
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
 * Starts the numeric integration from a rectangle that starts from the left side of the point.
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
 * Starts the numeric integration from a rectangle that starts from the middle of the point.
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


/**
 * Integrates the function f using the one third Simpon method
 * (h / 3) * ( f( x0 ) + f( x2N ) +  4*sum( f( x2i-1 ), 1, N ) + 2 * sum( f( x2i ), 1, N - 1 ) )
 * h = ( x2n - x0 ) / 2N 
 */
pub fn integreate_simpson(f: fn(f64) -> f64, x0: f64, x2n: f64, n: u32) -> f64 {
    let h = ( x2n - x0 ) / ( 2. * n as f64 );
    let mut sum = f( x0 ) + f( x2n );

    for i in 1..=n {
        let x_odd  = x0 + ( 2. * i as f64 - 1. ) * h;
        sum += 4. * f( x_odd );
    }

    for i in 1..n {
        let x_even = x0 + 2. * i as f64 * h;
        sum += 2. * f( x_even );
    }

    sum * h / 3.
}
