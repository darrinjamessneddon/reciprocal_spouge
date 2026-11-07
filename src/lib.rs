#![allow(clippy::module_inception)]
// THIS IS A WORK IN PROGRESS AND IS FAR FROM COMPLETE.
// THIS CODE IS NOT YET READY FOR PRODUCTION USE..

// This project was inspired by the idea of implementing the reciprocal gamma function using a re-arranged formula of the
// Spouge approximation of the gamma function, and so the repository was given the name reciprocal_spouge.

// The rationale for this was that the reciprocal gamma function has the following key properties:
// 1. It is an entire function (analytic across the entire complex plane).
// 2. It has zeros instead of singularites and these are located at negative integers.
// 3. It has no essential singularity (no poles), and as such is a completely smooth and well-behaved function.
//! ! References:
//! Wikipedia: en.wikipedia.org
//! The Ohio State University: A Historical Profile of the Gamma Function: https://people.math.osu.edu/gautam.42/S20/DavisGammaFunction.pdf
//! Google AI

// It was envisaged that computation of the reciprocal gamma function would be a good starting point to compute other functions in the gamma function family.

// The original unpublished public repository for reciprocal_spouge was deleted, and the project was re-started from scratch because of structure h submodules and the way the code was organized. The original code was also flawed. I humbly apologize to anyone who may have
// been inconvenienced by the original code, and I hope that this new version will be more useful to those who wish to use it.

// This library is built entirely in the Rust programming language.
// When the library was re-built the scope of the project was expanded to include the following functions:
// 1. The gamma function.
// 2. The reciprocal gamma function.
// 3. The log-gamma function.
// 4. The digamma function.
// 5. The trigamma function.

// In order to facilitate the implementation of these functions listed above, a custom Complex256 struct was created and the f256 crate was used as a dependency.

// This Complex256 struct has two traits: ComplexOps and ComplexComparisons.
// ComplexOps deals with mathematical operations.
// ComplexComparisons deals with a small number of boolean comparison operations. This contains few functions, but can be extended if users of this library
// wish to make greater use of the Complex256 struct as a feature in itself.

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
//! !references:
//! ! en.wikipedia.org
//! ! Google AI.

// The natural logarithm of the gamma function can be defined as:
// ln(Γ(z)) = integral from 0 to infinity of { (t^(z - 1) * e^(-t) - (e^(-t) - 1) / t) } dt))
// for a complex number z with a positive real part.
// There are also other ways to define the natural logarithm of the gamma function, including using the Spouge approximation,
// and using the Stirling approximation, which is a special case of the Spouge approximation.
//! !references:
//! ! en.wikipedia.org
//! Google AI.

// The digamma function can be defined as:
// ψ(z) = d/dz lnΓ(z)

// Main uses and applications of the digamma function include:
// (i) exploiting its relationship with harmonic numbers, especially to evaluate and sum complex or alternating series.
// (ii) using it in maximum likelihood estimation.
// (iii) using it in the regularization of divergent integrals (particularly in mathematical physics and analysis where
// the function is used in analytical continuation and the regularization of certain divergent series and integrals.

// The digamma function is a meromorphic function on the complex plane with well-defined poles and zeros.
// It has poles at all non-positive integers.
// At each pole where z = -n (for n >= 0) the residue is -1.
// The function is analytic (holomorphic) everywhere on the complex plane. It has no branch points, branch cuts, or essential singularities
// in the finite complex plane.
// The zeros of the digamma function correspond to the saddle points of hte complex gamma function.
// All zeros of the digamma function are real and simple.
// There is a single zero on the positive real axis located near z = 1.46163.
// There is an infinite sequence of negative real zeros.
// There are no zeros off the real axis in the complex plane.
//! ! references:
//! Scribd: Zeros of Digamma and Barnes G Functions | PDF|
//! Digamma Function- an overview| ScienceDirect Topics (ScienceDirect.com)
//! The zeros of the digamma function- MathOverflow

// The trigamma function can be defined as:
// ψ'(z) = d^2/dz^2 lnΓ(z)

// The trigamma funciton is a meromorphic function on the complex plane.
// It is holomorphic (analytic) everywhere on the complex plane except at its poles.
// It has poles at all the non-positive integers.
// Near each non-positive integer, the trigamma function has a pole of order 2 (a double pole).
// Since it is a quotient/derivative derived from the gamma function (which itself is analytic everywhere on the complex
// plane except for isolated poles), the trigamma function fits the precise definition of a meromorphic function.
//! ! references:
//! GitHub: diofant/diofant/functions/special/gamma_functions.py at master
//! Google AI.

pub mod core;

pub use core::{
    digamma, digamma_c256, factorial, ln_gamma, loggamma, rspouge, rspouge_c256, spouge,
    spouge_c256, spouge_coefficients, trigamma, trigamma_c256, Complex256, ComplexComparisons,
    ComplexOps, MathError, ParameterOutOfRangeError, ParseError, ParseResult, Result,
};

// The core module contains the following submodules:
// shared: a submodule providing functions that are shared between the other modules, including spouge_coefficients.
// complex: a submodule providing structs and functions for complex numbers, including complex numbers using f256 and f64 floats.
// gamma: a submodule holding the spouge approximation of the gamma function.
// rgamma: a submodule holding a spouge-based implementation of the reciprocal gamma function.
// error: a submodule providing error handling for the library (this is far from complete, and may need to be completely reworked).
// lngamma: a submodule holding an implementation of the natural logarithm of the gamma function.
// digamma: a submodule holding an implementation of the digamma function.
// trigamma: a submodule holding an implementation of the trigamma function.
