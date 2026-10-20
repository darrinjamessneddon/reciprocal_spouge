# Complex256 Usage

This example demonstrates adding two `Complex256` values with the `ComplexOps` trait.

## Addition

```rust
use f256::f256 as Float256;
use reciprocal_spouge::{Complex256, ComplexOps};

fn main() {
    let z1 = Complex256::new(Float256::from(1.0), Float256::from(2.0));
    let z2 = Complex256::new(Float256::from(3.0), Float256::from(4.0));

    let z3 = z1.add(z2);

    let expected = Complex256::new(Float256::from(4.0), Float256::from(6.0));
    assert!(z3 == expected);
}
```
