# reciprocal_spouge

`reciprocal_spouge` is a Rust library for experimenting with Spouge-based
approximations of the gamma function and the reciprocal gamma function in the
complex plane.

## Project status

This project is still a work in progress and is not ready for production use.
The current code is intended for learning, experimentation, and incremental
improvement.

## Current goals

The project currently aims to:

1. explore complex-number support built on the `f256` crate;
2. approximate the reciprocal gamma function directly with a modified Spouge
   formulation;
3. approximate the natural logarithm of the reciprocal gamma function; and
4. approximate the gamma function itself for values whose modulus is at most
   `10000`.

## Public API

The library currently exposes:

- `spouge_reciprocal::F256Complex` for complex values backed by `f256`;
- `r_spouge` for the reciprocal gamma approximation.

## Example

```rust
use f256::f256;
use reciprocal_spouge::spouge_reciprocal::r_spouge;
use reciprocal_spouge::spouge_reciprocal::F256Complex;

let gamma_input = F256Complex::new(f256::from(2.0), f256::from(0.0));
let reciprocal_gamma = r_spouge(gamma_input, f256::from(5.0));

println!("{:?}", reciprocal_gamma);
```

## Development

Run the existing test suite with:

```bash
cargo test
```

## References

- J. L. Spouge, *A New Approximation for the Gamma Function*, Mathematics of
  Computation 63(208), 159-176 (1994)
- <https://www.ams.org/journals/mcom/1994-63-208/S0025-5718-1994-1209956-2/S0025-5718-1994-1209956-2.pdf>
