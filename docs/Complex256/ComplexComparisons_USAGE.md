# This document demonstrates common Complex256 functions using the ComplexComparisons trait

## Approx_Eq

```rust
use reciprocal_spouge::{Complex256, ComplexComparisons, Float256};

fn main() {
    let z1 = Complex256::new(Float256::from(1.0), Float256::from(2.0));
    let z2 = Complex256::new(Float256::from(1.0 + 1e-15), Float256::from(2.0));
    let tol = Float256::from(1e-10);

    let result = z1.approx_eq(z2, tol);
    assert!(result);
}
```
