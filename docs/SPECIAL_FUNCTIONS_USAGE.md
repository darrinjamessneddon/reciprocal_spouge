# SPECIAL FUNCTIONS USAGE

## The Spouge Approximation

The Spouge approximation takes two variables as input, a complex number z, and a parameter "a".

The latter is (for our purposes) an integer between 2 and 151.

To improve the speed of computation a vector of Spouge Coefficients is pre-computed using the parameter "a".

The default setting for "a" is 100, but the value can be re-set to 25, 50, 75, 125 or 150 as desired. These values have been stored in a build script.

## Usage Examples

let a = 150;

let z = Complex256::new(Float256::from(5.0), Float256::from(0.0));

let four_factorial = factorial(4);

// Add tolerance value
let tol = 1e-90;

let gamma = spouge(z, a).expect("valid Spouge parameter");

assert!(gamma =)

let log_gamma_val = loggamma(z, a).expect("valid Spouge parameter");

assert!(log_gamma_val = 

let digamma_val = digamma(z, a).expect("valid Spouge parameter");

assert!(digamma_val = 

let trigamma_val = trigamma(z, a).expect("valid Spouge parameter");

assert!(trigamma_val = 



