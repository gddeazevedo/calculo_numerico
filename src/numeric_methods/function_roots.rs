use crate::functions::{derivative, precision_error, relative_error};


pub fn bisec_method(f: fn(f64) -> f64, mut a: f64, mut b: f64, epsilon: f64) -> Option<f64>
{
    if f(a) * f(b) > 0.0 {
        println!("f(a) must be negative and f(b) must be positive!");
        return None;        
    }

    let mut f_root: f64;

    loop {
        f_root = (b + a) / 2.0;

        if f(f_root) < 0.0 {
            a = f_root;
        } else {
            b = f_root;
        }

        if precision_error(a, b) < epsilon {
            break;
        }
    }

    Some(f_root)
}


pub fn secant_method()
{

}


/**
 * Calculate the square root of a floating point number with 64 bits
 * using the secant numeric method
 */
pub fn sqrt(n: f64) -> Option<f64>
{
    if n < 0.0 {
        println!("There is no square root for negative number!");
        return None;
    }

    let mut a = 0.0;
    let mut b = n;
    let epsilon = 1e-14;
    let mut root: f64;

    let f_sqrt: fn(f64, f64) -> f64 = |x, n| x*x - n;

    loop {
        let alpha = (f_sqrt(b, n) - f_sqrt(a, n)) / (b - a);
        let beta  = f_sqrt(a, n) - alpha * a;
        root = -beta / alpha;

        b = a;
        a = root;

        if precision_error(a, b) < epsilon {
            break;
        }
    }

    Some(root)
}


pub fn regula_falsi() -> f64
{
    0.0
}


pub fn linear_iteration_method(psi: fn(f64) -> f64, x0: f64) -> f64
{
    let mut xk = x0;
    let mut xk1 = 0.0;
    let epsilon = 1e-5;

    loop {
        xk1 = psi(xk);

        if precision_error(xk, xk1) < epsilon {
            break;
        }

        xk = xk1;
    }

    xk1
}


pub fn newton_method(f: fn(f64) -> f64, x0: f64) -> f64
{
    let mut xk = x0;
    let mut xk1 = 0.0;
    let epsilon = 1e-10;

    loop {
        xk1 = xk -  f(xk) / derivative(f, xk);

        if precision_error(xk, xk1) < epsilon {
            break;    
        }

        xk = xk1;
    }

    xk1
}
