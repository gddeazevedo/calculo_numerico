use crate::numeric_methods::helpers::{derivative, precision_error, relative_error};


/**
 * Finds the root of a one varible scalar function using the bisection numeric method
 * I = [a, b] -> f(a) * f(b) < 0, means that f has a root f(x) = 0
 */
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


/**
 * Finds the root of a one varible scalar function using the secant numeric method
 * 
 */
pub fn secant_method(f: fn(f64) -> f64, mut a: f64, mut b: f64) -> f64
{
    let e = 1e-14;
    let mut c: f64; // stores the function root

    loop {
        let alpha = (f(b) - f(a)) / (b - a);
        let beta = f(a) - alpha * a;
        c = -beta / alpha;
        b = a;
        a = c;

        if precision_error(a, b) < e {
            break
        }
    }

    c
}


/**
 * Finds the root of a one varible scalar function using the secant numeric method
 * 
 */
pub fn secant(f: fn(f64) -> f64, x0: f64, x1: f64) -> f64
{
    let mut xk = x0;
    let mut xk1 = x1;
    let mut xk2 = 0.0;
    let epsilon = 1e-14;

    loop {
        xk2 = (xk * f(xk1) - xk1 * f(xk)) / (f(xk1) - f(xk));

        if precision_error(xk1, xk2) < epsilon {
            break
        }

        xk = xk1;
        xk1 = xk2;
    }

    xk2
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


/**
 * Finds the root of a one varible scalar function using the Regula-Falsi numeric method
 * 
 */
pub fn regula_falsi(f: fn(f64) -> f64, x0: f64, x1: f64) -> f64
{
    if f(x0) * f(x1) > 0.0 {
        println!("f(a) must be negative and f(b) must be positive!");
        return 0.0;      
    }

    let epsilon = 1e-14;
    let mut xk  = x0;
    let mut xk1 = x1;
    let mut xk2 = 0.0;

    loop {
        xk2 = (xk * f(xk1) - xk1 * f(xk)) / (f(xk1) - f(xk)); // secant method formula

        if precision_error(xk, xk2) < epsilon || precision_error(xk1, xk2) < epsilon {
            break;
        }

        if f(xk2) * f(xk1) < 0.0 {
            xk = xk2;
        } else {
            xk1 = xk2;
        }
    }

    xk2
}


/**
 * Finds the root of a one varible scalar function using the linear iteration numeric method
 * 
 */
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


/**
 * Finds the root of a one varible scalar function using the Newton numeric method
 * 
 */
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
