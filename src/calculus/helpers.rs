use crate::types::Matrix;


pub fn derivative(f: fn(f64) -> f64, x: f64) -> f64
{
    let delta_x = 1e-10;
    (f(x + delta_x) - f(x)) / delta_x
}


pub fn precision_error(a: f64, b: f64) -> f64
{
    f64::abs(b - a) / f64::max(1.0, f64::abs(b))
}


pub fn absolute_error(a: f64, b: f64) -> f64
{
    f64::abs(b - a)
}


pub fn relative_error(a: f64, b: f64) -> f64
{
    f64::abs(b - a) / f64::abs(b)
}
