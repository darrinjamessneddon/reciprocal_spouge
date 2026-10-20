# SPECIAL FUNCTIONS USAGE

## The Spouge Approximation

The Spouge approximation takes two variables as input, a complex number z, and a parameter "a".

The latter must be an integer between 2 and 150, inclusive.

To improve the speed of computation a vector of Spouge Coefficients is pre-computed using the parameter "a".

The build script sets "a" to 100 by default, with Cargo features for 25, 50, 75, 125, and 150. The functions take "a" explicitly.

## Usage Examples

```rust
use f256::f256 as Float256;
use reciprocal_spouge::{
    digamma_c256, loggamma, spouge_c256, trigamma_c256, Complex256, ComplexOps,
};

const A: usize = 25;

fn assert_close(actual: Complex256, expected: Complex256, tolerance: f64) {
    let difference = actual.sub(expected).magnitude();
    assert!(
        difference <= Float256::from(tolerance),
        "expected {expected}, got {actual}, diff={difference}"
    );
}

fn main() {
    let z = Complex256::from_f64(1.0, 0.0);

    let gamma = spouge_c256(z, A).expect("valid Spouge parameter");
    assert_close(gamma, Complex256::from_f64(1.0, 0.0), 1e-12);

    let log_gamma_val = loggamma(z, A).expect("valid Spouge parameter");
    assert_close(log_gamma_val.exp(), gamma, 1e-12);

    let digamma_val = digamma_c256(z, A).expect("valid Spouge parameter");
    assert_close(
        digamma_val,
        Complex256::from_f64(-std::f64::consts::EULER_GAMMA, 0.0),
        1e-12,
    );

    let trigamma_val = trigamma_c256(z, A).expect("valid Spouge parameter");
    assert_close(
        trigamma_val,
        Complex256::from_f64(std::f64::consts::PI.powi(2) / 6.0, 0.0),
        1e-8,
    );
}
```


