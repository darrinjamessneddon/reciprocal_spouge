**reciprocal_spouge**

Welcome and thank you for taking an interest in this project. The original (unpublished) repository for this project was deleted, and the project re-built from scratch. I humbly apologise for any inconvenience this may have caused you.

This project aims create a library for f256 numbers and also a library of special functions including the reciprocal gamma, gamma, and log-gamma functions in the complex plane. These two sub-libraries have been bundled together as a single library crate.

This project was started because there are limits to the precision that can be achieved using Complex64 numbers based on f64 real and imaginary values. It is hoped that the construction of Complex256 numbers will achieve greater precision.

A major early goal was to use  re-arranged form of the Spouge approximation to compute the reciprocal gamma function, because it is a function that has no poles or singularities, and is analytic everywhere. Yet it was necessary to have a gamma function to compare its output with for testing purposes. So for this reason, the goals of the project were expanded to facilitate the calculation of other special functions.

This project is still very much a work in progress. At this point in time (the lngamma (log_gamma) submodule is yet to be created), and the error-handling module has a lot of work to be done on it.

This is a library created using Rust programming language because it is a language renowned for memory safety and concurrency.

The library uses the Spouge approximation because it is numerically stable.

The library is built on top of the f256 crate, and deals with special functions including the reciprocal gamma, gamma and log-gamma functions in the complex plane.

The Spouge approximation takes two inputs: the value for some complex number z, and the value of a positive integer parameter 'a'. The value of this parameter should be a > 2.

**Features

Complex256 struct in containing a number of public functions:

* {add, sub, mul, div} // standard arithmetical operations.
* {abs, arg, conj, magnitude}
* powc // z raised to complex power.
* powi // z raised to integer power.
* powf // z raised to non-whole number power.
* exp // e raised to the power of z.
* ln // natural logarithm of z
* log10 // log10 logarithm of z.
* recip // reciprocal of z.
* sqrt // square root of z.
* {sin, cos, tan} // sine, cosine and tangent of z.

**Functions for computing special functions:

rspouge(z, a) : takes a complex number, z, and a parameter 'a'with an integer value greater than 2 and returns the reciprocal gamma value as a string.

rspouge_c256(z, a): does the same thing but returns the reciprocal gamma value in f256 value form, so that it can be plugged directly into any other mathematical function you wish to create.

spouge(z, a): takes a complex number, z, and a parameter 'a' with an integer value greater than 2 and returns the gamma
value as a string.

spouge_c256(z, a): does the same thing but returns the gamma value function in f256 form, so that it can be plugged directly into any other mathematical function you wish to create.

// A log-gamma function is yet to be created.

// Seamless error-handling has yet to be added to these functions to prevent panics under certain conditions.

** Usage examples

let z1 = Complex256::new(5.0, 1.0);

let z2 = Complex256:new(3.0, 4.0);

let z3 = z1.add(&z2);

let z3_str = z3.to_string();

println!("z3 as string: {}", z3_str);

let z = Complex256::new(5.0, 0.0);

let a = 10;

let gamma = spouge(z, a);

let reciprocal_gamma = rspouge(z, a);

println!("gamma value for z: {}", gamma);

println!("reciprocal gamma value for z: {}", reciprocal_gamma);

**Installation Steps

git clone https://github.com/darrinjamessneddon/reciprocal_spouge

Cargo build

Create a binary application or test suite.

For contributor setup and workflow see [CONTRIBUTING.md](.github/CONTRIBUTING.md)

This project is licensed with the MIT license.
