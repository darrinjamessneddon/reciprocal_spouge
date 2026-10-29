# This document demonstrates Complex256 boolean operations using the ComplexComparisons trait

## Approx_Eq
// Compares whether two values of z are equal within a given absolute tolerance.

```rust
use reciprocal_spouge::{Complex256, ComplexComparisons, Float256};

fn main() {
    // Compare two z values where z2 is close to z1:
    let z1 = Complex256::new(Float256::from(1.0), Float256::from(2.0));
    let z2 = Complex256::new(Float256::from(1.0 + 1e-15), Float256::from(2.0));
    let tol = Float256::from(1e-10); // tolerance

    let result = z1.approx_eq(z2, tol);
    assert!(result);

    // Compare two z values where z3 == z1
    let z3 = Complex256(Float256::from(1.0), Float256::from(2.0));
    let tol = Float256::from(1e-10); // tolerance
    let result2 = z1.approx_eq(z3, tol);
    assert!(result2);
}
```

## Is_Zero
// Compares whether a value of z is zero to within a given absolute tolerance.

```rust
use reciprocal_spouge::{Complex256, ComplexComparisons, Float256};

fn main() {
    // Compare a z value to zero within a given tolerance:
    let z1 = Complex256::new(Float256::from(1e-60), Float256::from(0.0));
    let tol = Float256::from(1e-10);

    let result = z1.is_zero(tol);
    assert!(result);

    // Compare a z value to zero
    let z2 = Complex256::new(Float256::from(0.0), Float256::from(0.0));
    let tol = Float256::from(1e-10);

    let result2 = z2.is_zero(tol);
    assert!(result2);
}
```

## Is_i

```rust
use reciprocal_spouge::{Complex256, ComplexComparisons, Float256};

fn main() {
    let z1 = Complex256::new(Float256::from(1e-60), Float256::from(1.0));
    let tol = Float256::from(1e-10);

    let result = z1.is_i(tol);
    assert!(result);

    let z2 = Complex256::new(Float256::from(0.0), Float256::from(1.0));
    let tol = Float256::from(1e-10);

    let result2 = z2.is_i(tol);
    assert!(result2);
}
```

## Is_real

```rust
use reciprocal_spogue::{Complex256, ComplexComparisons, Float256};

fn main() {
    let z1 = Complex256::new(Float256::from(35.0), Float256::from(1e-60));
    let tol = Float256::from(1e-10);

    let result = z1.is_real(tol);
    assert!(result);

    let z2 = Complex256::new(Float256::from(35.0), Float256::from(0.0));
    let tol = Float256::from(1e-10);

    let result2 = z2.is_real(tol);
    assert!(result2);
}
```

## is_imaginary

```rust
use reciprocal_spouge::{Complex256, ComplexComparisons, Flaot256};

fn main() {
    let z1 = Complex256::new(Float256::from(1e-60), Float256::from(35.0));
    let tol = Float256::from(1e-10);

    let result = z1.is_imaginary(tol);
    assert!(result);

    let z2 = Complex256::new(FLoat256::from(0.0), Float256::from(35.0));
    let tol = Float256::from(1e-10);

    let result2 = z2.is_imaginary(tol);
    assert!(result2);
}
```
