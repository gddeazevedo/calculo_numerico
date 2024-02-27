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

    let mut c: f64;

    loop {
        c = (b + a) / 2.0;

        if f(c) < 0.0 {
            a = c;
        } else {
            b = c;
        }

        if (b - a) / 2.0 < epsilon {
            break;
        }
    }

    Some(c)
}


pub fn secant_method(mut a: f64, mut b: f64, epsilon: f64) -> Option<f64>
{
    let mut c: f64;

    loop {
        let alpha = (f(b) - f(a)) / (b - a); // linear coefficient
        let beta  = f(a) - alpha * a;
        c = -beta / alpha;

        b = a;
        a = c;

        if f64::abs(b - a) / f64::max(1 as f64, f64::abs(b)) < epsilon {
            break;
        }
    }

    Some(c)
}