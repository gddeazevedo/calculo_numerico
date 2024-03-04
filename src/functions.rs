pub fn f(x: f64) -> f64
{
    f64::exp(x) + f64::sin(x*x) - 10.0
}

pub fn derivative(f: fn(f64) -> f64, x: f64) -> f64
{
    let delta_x = 1e-10;
    (f(x + delta_x) - f(x)) / delta_x
}
