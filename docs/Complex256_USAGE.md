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

## Powc

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

## Powf

```rust
use f256::f256 as Float256;
use reciprocal_spouge::{Complex256, ComplexOps};

fn main() {
    let z = Complex256::new(Float256::from(1.0), Float256::from(2.0));
    let exponent = Float256:from("3.1415926");
    let result = z.powf(exponent);

    let expected_re_str = "-11.82646698597011528465526000285883763094322656386824274409109082964041096";
    let expected_im_str = "-4.138850340076041866283067425612456813067533929442880780012629783710972797";
    let expected = Complex256::new(Float256::from(expected_re_str), Float256::from(expected_im_str));

    assert_eq!(result, expected);
}
```

This example demonstrates raising a `Complex256` value to the power of an integer (i32) with the `ComplexOps` trait.

## Powi

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

This example demonstrates taking the exponent of a `Complex256` value with the `ComplexOps` trait.

## Exponential

```rust
use f256::f256 as Float256;
use reciprocal_spouge::{Complex256, ComplexOps};

fn main() {
    let z = Complex256::new(Float256::from(1.0), Float256::from(2.0));
    let result = z.exp();

    let expected_re_str = "-1.13120438375681363843125525551079471062886799582652575021772191041650193";
    let expected_im_str = "2.47172667200481892761693089355166453273619036924100818420075883527783967";
    let expected = Complex256::new(Float256::from(expected_re_str), Float256::from(expected_im_str));

    assert_eq!(result, expected);
}
```

This example demonstrates taking the natural logarithm of a `Complex256` value with the `ComplexOps` trait.

## Natural Logarithm

```rust
use f256::f256 as Float256;
use reciprocal_spouge::{Complex256, ComplexOps};

fn main() {
    let z = Complex256::new(Float256::from(1.0), Float256::from(2.0));
    let result = z.ln();

    let expected_re_str = "0.8047189562170501873003796666130938197628006771342588609563239457370849";
    let expected_im_str = "1.107148171779409050301706546017853704007004764540143264667653920743371034";
    let expected = Complex256::new(Float256::from(expected_re_str), Float256::from(Complex_im_str));

    assert_eq!(result, expected);
}
```

This example demonstrates taking the logarithm base10 of a `Complex256` value with the `ComplexOps` trait.

## Logarithm Base10

```rust
use f256::f256 as Float256;
use reciprocal_spouge::{Complex256, ComplexOps};

fn main() {
    let z = Complex256::new(Float256::from(1.0), Float256::from(2.0));
    let result = z.log10();

    let expected_re_str = "0.3494850021680093694455211138700073647695012827633090226305997180333825356";
    let expected_in_str = "0.480828578784234057373949959119329899707069251685754086023650078002816663";
    let expected = Complex256::new(Float256::from(expected_re_str), Float256::from(expected_im_str));

    assert_eq!(result, expected);
}
```

This example demonstrates taking the square root of a `Complex256` value with the `ComplexOps` traits.
It returns one of the square roots

## Square Root

```rust
use f256::F256 as Float256;
use reciprocal_spouge::{Complex256, ComplexOps};

fn main() {
    let z = Complex256::new(Float256:;from(1.0), Float256::from(2.0));
    let result = z.sqrt();

    let expected_re_str = "1.27201964951406896425242246173749149171560804184009624861664038253929758;
    let expected_im_str = "0.78615137775742328606955858584295892952312205783772323766490197010118205;
    let expected = Complex256::new(Float256::from(expected_re_str), Float256::from(expected_im_str));

    assert_eq!(result, expected);
}
```

This example demonstrates taking the reciprocal of a `Complex256` value with the `ComplexOps` trait.

## Reciprocal

```rust
use f256::f256 as Float256;
use reciprocal_spouge::{Complex256, ComplexOps};

fn main() {
    let z = Complex256::new(Float256::from(1.0), Float256::from(2.0));
    let result = z.recip();

    let expected = Complex256::new(Float256::from(0.2), Float256::from(-0.4));

    assert_eq!(result, expected);
}
```

This example demonstrates taking the sine of a `Complex256` value with the `ComplexOps` trait.

## Sine

```rust
use f256;:f256 as Float256;
use reciprocal_spouge::{Complex256, ComplexOps};

fn main() {
    let z = Complex256::new(Float256::from(1.0), Float256::from(2.0));
    let result = z.sin();
    let expected_re_str = "3.1657785132161681467407346171919055383791107678914689322893274643048881";
    let expected_im_str = "1.95960104142160589707035204998935827843632016018455965880190149916309425";
    let expected = Complex256::new(Float256::from(expected_re_str), Float256::from(expected_im_str));

    assert_eq!(result, expected);
}
```

This example demonstrates taking the sine of a `Complex256` value with the `ComplexOps` trait.

## Cosine

```rust
use f256::f256 as Float256;
use reciprocal_spouge::{Complex256, ComplexOps};

fn main() {
    let z = Complex256::new(Float256;:from(1.0), Float256::from(2.0));
    let result = z.cos();
    let expected_re_str = "2.03272300701966552943634344849951426373199040663875238194522463491499985";
    let expected_in_str = "3.0518977991518000575121156868951054528884376177333196446633774753865543";
    let expected = Complex256::new(Float256::from(expected_re_str), Float256::from(expected_im_str));

    assert_eq!(result, expected);
}
```

This example demonstrates taking the tangent of a `Complex256` value with the `ComplexOps` trait.

## Tangent

```rust
use f256::f256 as Float256,
use reciprocal_spouge::{Complex256, ComplexOps};

fn main() {
    let z = Complex256::new(Float256::from(1.0), Float256::from(2.0));
    let result = z.sin();
    let expected_re_str = "0.0338128260798966902843705597252873016401876693307558743667301120018375007;
    let expected_im_str = "1.014793616146633568117054175417967614163217471123533242464958508230529975;
    let expected = Complex256::new(Float256::from(expected_re_str), Float256::from(expected_im_str));

    assert_eq!(result, expected);
}
```

This example demonstrates taking the cotangent of a `Complex256` value with the `ComplexOps` trait.

## Cotangent

```rust
use f256::f256 as Float256;
use reciprocal_spouge::{Complex256, ComplexOps};

fn main() {
    let z = Complex256;:new(Float256::from(1.0), Float256::from(2.0));
    let result = z.cot();
    let expected_re_str = "0.032797755533752594062764565765830629341580233549321155658845723715218373;
    let expected_im_str = "-0.984329226458191029471888181689464448193903065068330455889400082143752147;
    let expected = Complex256::new(Float256::from(expected_re_str), Float256::from(expected_im_str));

    assert_eq!(result, expected);
}
```

These examples demonstrate conversion of `Complex256` values to `Complex64` values with the `ComplexOps` trait.

```rust
use f256::f256 as Float256;
use reciprocal_spouge::{Complex256, ComplexOps};

fn main() {
    let z = Complex256;:new(Float256::from(1.0), Float256::from(2.0));
    let z_complex64 = z.to_complex64();
    let z_c64_expected = Complex64::new(from(1.0), from(2.0));

    assert_eq!(z_complex64, z_c64_expected);

    let z_cot_c64 = z.cot().to_complex64();
    let z_cot_expected = Complex {re: 0.032797755552596, im; -0.984329226458191};
    assert_eq!(z_cot_c64, z_cot_expected);
}
```



    



    



    



