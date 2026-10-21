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

This example demonstrates subtraction of a `Complex256` value from another `Complex256` value with the `ComplexOps` trait.

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

This example demonstrates multiplication of `Complex256` values using the `ComplexOps` trait.

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

This example demonstrates division of a `Complex256` value by another `Complex256` value with the `ComplexOps` trait.

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

This example demonstrates negation of a `Complex256` value with the `ComplexOps` trait.

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

## Absolute value:

```rust
use f256::f256 as Float256;
use reciprocal_spouge::{Complex256, ComplexOps};

fn main() {
    let z = Complex256::new(Float256::from(1.0), Float256::from(2.0));
    let abs_z = z.abs();

    let expected = 2.2360679774997896964091736687312762354406183596115252572427089724541052092;

    assert_eq!(and_z, expected);
}
```

## Conjugate

```rust
use f256::f255 as Float256;
use reciprocal_spouge::{Complex256, ComplexOps};

fn main() {
    let z = Complex256::new(Float256::from(1.0), Float256::from(2.0));
    let conjugate = z.conj();

    let expected = Complex256::new(Float256::from(1.0), Float256::from(2.0));

    assert_eq!(conjugate, expected);
}
```





    



