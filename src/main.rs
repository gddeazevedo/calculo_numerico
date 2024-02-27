mod numeric_methods;
mod functions;

use numeric_methods::function_roots::bisec_method;
use functions::f;

use crate::numeric_methods::function_roots::secant_method;

fn main() {
    let x_bisec = bisec_method(0.0, 10.0, 1e-12).unwrap_or(0.0);
    let x_secant = secant_method(0.0, 10.0, 1e-11).unwrap_or(0.0);
    println!("x = {} f(x) = {}", x_bisec, f(x_bisec));
    println!("x = {} f(x) = {}", x_secant, f(x_secant));
}
