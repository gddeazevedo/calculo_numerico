use crate::types::Matrix;


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
