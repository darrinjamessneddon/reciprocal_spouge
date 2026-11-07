# Installation

## Requirements
- Rust and Cargo
- MSRV = 1.86.0

## Add the crate
Once published:
```sh
cargo add reciprocal_spouge
```

## Or add it to Cargo.toml
```
[dependencies]
reciprocal_spouge = "0.1"
```

## If you are computing one or more special functions:
* Select a value for the spouge parameter "a": 80, 90 or 100.
* a = 100 is the default value if no selection is made.
### A larger value of the Spouge parameter:
[1] drastically reduces error
When the real part of z is positive (Re(z) > 0), the relative error of the approximation is bounded by:
Error <= a^-0.5 * (2* pi)^-(a + 0.5)
  
[2] increases computation cost
The summation is Spouge's formula runs from k = 1 to a - 1
* As the value of "a" becomes larger more terms must be calculated in the sum. This means you need to precompute
and store a larger vector of coefficients.

[3] modifies the convergence domain
* As Spouge's approximation is valid over the entire complex plane (except where z is a negative integer)
and the approximation is most accurate when the real part of z is greater than negative a, increasing the Spouge parameter,
effectively shifts the "safe" region of high accuracy further left along the negative real axis.
A larger sum 
[3] modifies the convergence domain
