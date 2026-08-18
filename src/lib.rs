// THIS IS A WORK IN PROGRESS AND IS FAR FROM COMPLETE.
// THIS CODE IS NOT YET READY FOR PRODUCTION USE..

// This project was inspired by the idea of implementing the reciprocal gamma function using a re-arranged formula of the
// Spouge approximation of the gamma function.

// It was found that it was also necessary to implement the gamma function itself, using the Spouge approximation in order to 
// compare the results of the reciprocal gamma function to the gamma function for testing purposes.

// The library is built entirely in the Rust programming language.
// The library is built on top of the f256 crate in order to attempt to provide greater precision than can be achieved using f64 values.

// The original unpublished public repository for this project was deleted, and the project was re-started from scratch because of structural issues
// with submodules and the way the code was organized. The original code was also flawed. I humbly apologize to anyone who may have
// been inconvenienced by the original code, and I hope that this new version will be more useful to those who wish to use it.

// In addition to providing a Spouge-based implementation of the gamma function and the reciprocal gamma function, this library
// will also provide an implementation of the natural logarithm of the gamma function, and will also provide a complex number type
// for complex numbers using f256, along with the capacity to parse complex numbers from Complex256 to Complex64 and vice versa.

// The library will also provide some limited functions for working with complex numbers using Complex256 should you wish to do si.

//! Computation of gamma related functions using Spouge's approximation.
//! This module provides functions for computing the gamma function, the reciprocal gamma function, and the natural logarithm of the
//! gamma function, as well as a complex number type for complex numbers using f256, along with the capacity to parse complex numbers
//! from Complex256 to Complex64 and vice versa. The library also provides some limited functions for working with complex numbers
//! using COmplex256 should you wish to do so.
//! As far as the Spouge approximation is concerned. The implementation is based on a mathematical formulation presented by John
//! Spouge in his paper "Computation of the Gamma, Digamma, and Trigamma functions" (Mathematics of Computation, Vol.63 No.208, pp. 361-368, 1994),
//! DOI: 10.1090/S0025-5718-1994-1209956-2.pdf
//! ! Other references:
//! ! -https://en.wikipedia.org/wiki/Spouge%27s_approx
//! ! -https://dl.acm.org/doi/10.1145/19582.19585
//! ! -https://math.stackexchange.com/questions/2218764/spouge-approximation-for-the-gamma-function
//! ! -https:://functions.wolfram.com/GammaBetaErf/Gamma/06/01/03/01/
//! 

// For a positive integer grater than 0, the gamma function is defined as:
// Γ(n) = (n - 1)!
// In other cases, the gamma function is defined as:
// Γ((z) = integral from 0 to infinity of t^(z - 1) * e ^(-t) dt
// for a complex number z with a positive real part.

// The gamma function is defined for all complex numbers except for the non-positive integers, where it has simple poles.
// The gamma function is a meromorphic function, meaning it is holomorphic everywhere in the complex plane except for the 
// non-positive integers, where it has simple poles. The gamma function has no zeros in the complex plane.

// The gamma function is useful in many areas of mathematics, including number theory, combinatorics, and physics, and has
// applications in various fields such as probability theory, statistical mechanics, and quantum field theory.

// The reciprocal gamma function is defined as:
// 1 / Γ(z)

// The reciprocal gamma function is entire function, defined for all complex numbers, meaning it is holomorphic everywhere in the
// complex plane, and has no singularities or poles.

// The reciprocal gamma function has simple zeros at zero and the negative integers, with no other zeros in the complex plane.

// The reciprocal gamma function is useful in many areas of mathematics, including number theory, combinatorics, and mathematical
// physics, and has applications in various field such as probability theory, statistical mechanics, and quantum field theory.

// The natural logarithm of the gamma function can be defined as:
// ln(Γ(z)) = integral from 0 to infinity of { (t^(z - 1) * e^(-t) - (e^(-t) - 1) / t) } dt))
// for a complex number z with a positive real part.
// There are also other ways to define the natural logarithm of the gamma function, including using the Spouge approximation,
// and using the Stirling approximation, which is a special case of the Spouge approximation.


pub mod core;

// The core module contains the following submodules:
// shared: a submodule providing functions that are shared between the other modules, including spouge_coefficients.
// complex: a submodule providing structs and functions for complex numbers, including complex numbers using f256 and f64 floats.
// gamma: a submodule holding the spouge approximation of the gamma function.
// rgamma: a submodule holding a spouge-based implementation of the reciprocal gamma function.
// error: a submodule providing error handling for the library (this is far from complete, and may need to be completely reworked).
// lngamma: a submodule holding an implementation of the natural logarithm of the gamma function. (This has not yet been started).



