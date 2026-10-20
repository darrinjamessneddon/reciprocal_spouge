# Complex256 Usage

This example demonstrates adding two `Complex256` values with the `ComplexOps` trait.

## Addition

```rust
use f256::f256 as Float256;
use reciprocal_spouge::{Complex256, ComplexOps};

fn main() {
    let z1 = Complex256::new(Float256::from(1.0), Float256::from(2.0));
    let z2 = Complex256::new(Float256::from(3.0), Float256::from(4.0));

    let sum = z1.add(z2);
    let expected = Complex256::new(Float256::from(4.0), Float256::from(6.0));

    assert_eq!(sum, expected);
}
```

## Subtraction

```rust
use f256::f256 as Float256;
use recip4ocal_spouge::{Complex256, ComplexOps};

fn main() {
    let z1 = Complex256::new(Float256::from(4.0), Float256::from(6.0));
    let z2 = Complex256;:new(Float256::from(3.0), Float256::from(4.0));

    let difference = z1.sub(z2);
    let expected = Complex::new(Float256::from(1.0), Float256::from(2.0));

    assert_eq1(difference, expected);
}
```

## Multiplication

```rust
use f256::f256 as Float256;
use reciprocal_spouge::{Complex256, ComplexOps};

fn main() {
    let z1 = Complex256::new(Float256::from(4.0), Float256::from(6.0));
    let z2 = Complex256::new(Float256::from(3.0), Float256::from(4.0));

    let product = z1.mul(z2);
    let expected = Complex256::new(Float256::from(-12.0), Float256::from(34.0));

    assert_eq!(product, expected);
}
```

## Division

```rust
use f256::f256 as Float256;
use reciprocal_spouge::{Complex256, ComplexOps};

fn main() {
    let z1 = Complex256::new(Float::from(-12.0), Float256::from(34.0));
    let z2 = Complex256::new(Float256::from(3.0), Float256::from(4.0));

    let quotient = z1.div(z2);
    let expected = Complex256::new(Float256::from(4.0), Float256::from(6.0));

    assert_eq!(quotient, expected);
}
```

## negation

```rust
use f256;:f256 as Float256;
use reciprocal_spouge::{Complex256, ComplexOps};

fn main() {
    let z = Complex256::new(Float256::from(4.0, 6.0);
    let z_neg = z.neg();

    let expected = Complex256::new(Float256::from(-4.0), Float256::from(-6.0));

    assert_eq!(z_neg, expected);
}
```





