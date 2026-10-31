# RECIPROCAL SPOUGE

## About
This rust library is to enable users to compute special functions with high precision and accuracy using the Spouge approximation.

Functions implemented include reciprocal gamma, gamma, log-gamma, digamma and trigamma.

The library also includes a custom struct Complex256 containing functions pertaining to complex numbers with f256 real and imaginary parts.

## Prerequisites & requirements

- Rust and Cargo installed
- A supported stable Rust toolchain
- Internet access to download crate dependencies
- This library depends on `num`, num-traits`, `num-complex`, `num-bigint`, `f256`, `astro_float`, and `rayon`
- For testing, `proptest` is used
- The project is still under active development and is not yet production-ready

## Installation

Installation guidelines are set out in docs/INSTALLATION.md

## Special Functions Usage Examples

Usage examples for the special functions mentioned above are set out in docs/SPECIAL_FUNCTIONS_USAGE.md

## Complex256 Usage Examples

The Complex256 struct contains two traits:

* ComplexOps deals with mathematical operations.
* ComplexComparisons deals with boolean operations.

Usage examples for Complex256::ComplexOps are set out in docs/Complex256/ComplexOps_USAGE.md

Usage examples for Complex256::ComplexComparisons are set in docs/Complex256/ComplexComparisons_USAGE.md

## Contributing Guidelines

Contributing guidelines are set out in the CONTRIBUTING.md

## Getting in touch:

[Join our Zulip Chat] https://reciprocal-spouge.zulipchat.com/join/542cjcee57n7a46r6tkbgud7/

## License
This project is licensed with the MIT License



