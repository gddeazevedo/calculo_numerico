use crate::functions::f;


pub fn bisec_method(mut a: f64, mut b: f64, epsilon: f64) -> Option<f64>
{
    if f(a) > 0.0 && f(b) < 0.0 {
        let tmp = a;
        a = b;
        b = tmp;
    }

    if f(a) > 0.0 {
        println!("Error, invalid value for a");
        return None;
    }

    if f(b) < 0.0 {
        println!("Error, invalid value for a");
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

        if (b - a) / 2.0 < epsilon {
            break;
        }
    }

    Some(f_root)
}


pub fn secant_method(mut a: f64, mut b: f64, epsilon: f64) -> Option<f64>
{
    let mut f_root: f64;

    loop {
        let alpha = (f(b) - f(a)) / (b - a); // linear coefficient
        let beta  = f(a) - alpha * a; // x = a -> f(a) = alpha * a + beta -> beta = f(a) - alpha * a
        f_root = -beta / alpha; // alpha*x+beta = 0 -> x = -beta / alpha 

        b = a;
        a = f_root;

        if f64::abs(b - a) / f64::max(1 as f64, f64::abs(b)) < epsilon {
            break;
        }
    }

    Some(f_root)
}


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

    let f_sqrt: fn(f64, f64) -> f64 = |x:f64, n:f64| x*x - n;

    loop {
        let alpha = (f_sqrt(b, n) - f_sqrt(a, n)) / (b - a);
        let beta  = f_sqrt(a, n) - alpha * a;
        root = -beta / alpha;

        b = a;
        a = root;

        if f64::abs(b - a) / f64::max(1.0, f64::abs(b)) < epsilon {
            break;
        }
    }

    Some(root)
}
