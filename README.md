**reciprocal_spouge**

Welcome and thank you for taking an interest in this project. 

This project aims create a library for f256 numbers and also a library of special functions including the reciprocal gamma, gamma, and log-gamma functions in the complex plane. These two sub-libraries have been bundled together as a single library crate.

This project was started because there are limits to the precision that can be achieved using Complex64 numbers based on f64 real and imaginary values. It is hoped that the construction of Complex256 numbers will achieve greater precision.

A major early goal was to use  re-arranged form of the Spouge approximation to compute the reciprocal gamma function, because it is a function that has no poles or singularities, and is analytic everywhere. Yet it was necessary to have a gamma function to compare its output with for testing purposes. So for this reason, the goals of the project were expanded to facilitate the calculation of other special functions.

This project is still very much a work in progress, and the error-handling module has a lot of work to be done on it.

It might be the case that the scope of the project may need to be broadened if the spouge, or reciprocal spouge functions do not give enough precision. If this is the case other ways to compute the gamma and reciprocal gamma function may need to be included to give users of this library more options.

The library is built on top of the f256 crate, and deals with special functions including the reciprocal gamma, gamma and log-gamma functions in the complex plane.

The Spouge approximation takes two inputs: the value for some complex number z, and the value of a positive integer parameter 'a'. The value of this parameter should be a >= 2.

**Features

Complex256 struct in containing a number of public functions:

add, sub, mul, div, abs, arg, conj, magnitude, powc, powi, powf, exp, ln , log10, recip, sqrt, sin, cos, tan.

**Functions for computing special functions:

rspouge(z, a) : takes a complex number, z, and a parameter 'a' with an integer value greater than or equal to 2 and returns `Result<String, MathError>`.

rspouge_c256(z, a): does the same thing but returns `Result<Complex256, MathError>` so that coefficient and range errors can be handled explicitly.

spouge(z, a): takes a complex number, z, and a parameter 'a' with an integer value greater than or equal to 2 and returns `Result<String, MathError>`.

spouge_c256(z, a): does the same thing but returns `Result<Complex256, MathError>` so that coefficient and range errors can be handled explicitly.

ln_gamma(z, a): returns `Result<Complex256, MathError>` for the principal-branch natural logarithm of the gamma function away from poles.

// Coefficient-dependent Spouge APIs now return Result values and report invalid parameters with MathError.

** Usage examples

* if linking to an executable file:

    cargo new special_functions --bin

    cargo build

in the Cargo.toml file add

    reciprocal_spouge = { git = https://github.com/darrinjamessneddon/reciprocal_spouge }

    f256 = "0.11.2"

in the main.rs file add the following:
  
    use reciprocal_spouge::{Complex256, ComplexOps, MathError, rspouge, rspouge_c256, spouge, spouge_c256};

    use f256::f256 as Float256;

in fn main() add:

    let z = Complex256::new(Float256::from(5.0), Float256::from(0.0));

    let a = 80_usize; // Can use a lower value if desired.

    let gamma = spouge(z, a).expect("valid Spouge parameter");

    let reciprocal_gamma = rspouge(z, a).expect("valid Spouge parameter");

    println!("gamma value for z: {}", gamma); // Returns gamma value as a string.

    println!("reciprocal gamma value for z: {}", reciprocal_gamma);// Returns rgamma value as a string.

    let gamma_256 = spouge_c256(z, a).expect("valid Spouge parameter");

    let rgamma_256 = rspouge_c256(z, a).expect("valid Spouge parameter");

    let gamma_64 = gamma_256.to_complex64();

    let rgamma_64 = rgamma_256.to_complex64();

    println!("gamma value as Complex64: {}", gamma_64);

    println!("reciprocal gamma value as Complex64: {}", rgamma_64);

* if linking to an executable file to use functions pertaining directly to complex numbers:

    cargo new my_app --bin

    cargo build
 
* add the dependencies to the Cargo.toml file in the same way as shown above.

in fn main() add:
  
    let z1 = Complex256::new(Float256::from(1.0), Float256::from(2.0));

    println!("z1: {}", z1);

    let z2 = Complex256::new(Float256::from(2.0), Float256::from(3.0));

    println!("z2: {}", z2);

    let z3 = z1.add(z2);

    println!("z1 plus z2: {}", z3);
  

**Installation Steps

git clone https://github.com/darrinjamessneddon/reciprocal_spouge

Cargo build

Create a binary application or test suite.

For contributor setup and workflow see [CONTRIBUTING.md](.github/CONTRIBUTING.md)

This project is licensed with the MIT license.
