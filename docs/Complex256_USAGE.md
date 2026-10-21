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
use reciprocal_spouge::{Complex256, ComplexOps};

fn main() {
    let z1 = Complex256::new(Float256::from(4.0), Float256::from(6.0));
    let z2 = Complex256::new(Float256::from(3.0), Float256::from(4.0));

    let difference = z1.sub(z2);
    let expected = Complex256::new(Float256::from(1.0), Float256::from(2.0));

    assert_eq!(difference, expected);
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
    let z1 = Complex256::new(Float256::from(-12.0), Float256::from(34.0));
    let z2 = Complex256::new(Float256::from(3.0), Float256::from(4.0));

    let quotient = z1.div(z2);
    let expected = Complex256::new(Float256::from(4.0), Float256::from(6.0));

    assert_eq!(quotient, expected);
}
```

This example demonstrates negation of a `Complex256` value with the `ComplexOps` trait.

## Negation

```rust
use f256::f256 as Float256;
use reciprocal_spouge::{Complex256, ComplexOps};

fn main() {
    let z = Complex256::new(Float256::from(4.0), Float256::from(6.0));
    let z_neg = z.neg();

    let expected = Complex256::new(Float256::from(-4.0), Float256::from(-6.0));

    assert_eq!(z_neg, expected);
}
```

This example demonstrates the absolute value (modulus) of a `Complex256` value with the `ComplexOps` trait.

## Absolute Value

```rust
use f256::f256 as Float256;
use reciprocal_spouge::{Complex256, ComplexOps};

fn main() {
    let z = Complex256::new(Float256::from(1.0), Float256::from(2.0));
    let abs_z = z.abs();

    // sqrt(1^2 + 2^2) = sqrt(5)
    let expected = Float256::from("2.2360679774997896964091736687312762354406183596115257242708972454105209256");

    assert_eq!(abs_z, expected);
}
```

This example demonstrates the conjugate of a `Complex256` value with the `ComplexOps` trait.

## Conjugate

```rust
use f256::f256 as Float256;
use reciprocal_spouge::{Complex256, ComplexOps};

fn main() {
    let z = Complex256::new(Float256::from(1.0), Float256::from(2.0));
    let conjugate = z.conj();

    let expected = Complex256::new(Float256::from(1.0), Float256::from(-2.0));

    assert_eq!(conjugate, expected);
}
```

This example demonstrates the magnitude of a `Complex256` value with the `ComplexOps` trait.

## Magnitude

```rust
use f256::f256 as Float256;
use reciprocal_spouge::{Complex256, ComplexOps};

fn main() {
    let z = Complex256::new(Float256::from(1.0), Float256::from(2.0));
    let mag_z = z.magnitude();

    // sqrt(1^2 + 2^2) = sqrt(5)
    let expected = Float256::from("2.2360679774997896964091736687312762354406183596115257242708972454105209256");

    assert_eq!(mag_z, expected);
}
```

This example demonstrates raising a `Complex256` value to the power of another `Complex256` value with the `ComplexOps` trait.

## Powc Function

```rust
use f256::f256 as Float256;
use reciprocal_spouge::{Complex256, ComplexOps};

fn main() {
    let z1 = Complex256::new(Float256::from(1.0), Float256::from(2.0));
    let z2 = Complex256::new(Float256::from(3.0), Float256::from(4.0));
    let result = z1.powc(z2);

    let expected_re_str = "0.129009594074466894077052339652447244091845464472294727770700393582527832";
    let expected_im_str = "0.03392409290517012669761785462254790154054732022267760839965165515194522";
    let expected = Complex256::new(Float256::from(expected_re_str), Float256::from(expected_im_str));

    assert_eq!(result, expected);
}
```

This example demonsgrates raising a `Complex256` value to the power of a `f256` value with the `ComplexOps` trait.

## Powf Function

```rust
use f256::f256 as Float256;
use reciprocal_spouge::{Complex256, ComplexOps};

fn main() {
    let z = Complex256::new(Float256::from(1.0), Float256::from(2.0));
    let exponent = Float256:from("3.1415926");
    let result = z.powf(exponent);

    let expected_re_str = "-11.82646698597011528465526000285883763094322656386824274409109082964041096"
    let expected_im_str = "-4.138850340076041866283067425612456813067533929442880780012629783710972797"
    let expected = Complex256::new(Float256::from(expected_re_str), Float256::from(expected_im_str));

    assert_eq!(result, expected);
}
```

This example demonstrates raising a `Complex256` value to the power of an integer (i32) with the `ComplexOps` trait.

## Powi Function

```rust
use f256::f256 as Float256;
use reciprocal_spouge::{Complex256, ComplexOps};

fn main() {
    let z = Complex256::new(Float256::from(1.0), Float256::from(2.0));
    let exponent = 3;
    let result = z.powi(exponent);

    let expected = Complex256::new(Float256::from(-11.0), Float256::from(-2.0));

    assert_eq!(result, expected);
}
```



    



    



