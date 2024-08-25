pub fn derivative(f: fn(f64) -> f64, x: f64) -> f64
{
    let delta_x = 1e-10;
    (f(x + delta_x) - f(x)) / delta_x
}

pub fn df_central(f: fn(f64) -> f64, x: f64) -> f64 {
    let h = 1e-5;
    (f(x + h) - f(x - h)) / (2. * h)
}

pub fn df_successor(f: fn(f64) -> f64, x: f64) -> f64 {
    let h = 1e-10;
    (f(x) - f(x - h)) / h
}

pub fn df_progressive(f: fn(f64) -> f64, x: f64) -> f64 {
    let h = 1e-10;
    (f(x + h) - f(x)) / h
}

pub fn df_taylor(f: fn(f64) -> f64, x: f64) -> f64 {
    let h = 1e-10;
    ((-3.) * f(x) + 4. * f(x + h) - f(x + 2. * h)) / (2. * h)
}
