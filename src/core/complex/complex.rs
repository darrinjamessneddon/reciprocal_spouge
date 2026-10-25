pub mod complex {

    use f256::f256 as Float256;
    use num_complex::Complex64;

    use crate::{MathError, ParseError};

    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct Complex256 {
        pub re: Float256,
        pub im: Float256,
    }
    impl Complex256 {
        pub fn new(re: Float256, im: Float256) -> Self {
            Complex256 { re, im }
        }
        pub fn try_to_f64(self) -> Result<(f64, f64), ParseError> {
            let convert = |x: Float256| {
                x.to_string()
                    .parse::<f64>()
                    .ok()
                    .filter(|v| v.is_finite())
                    .ok_or(ParseError::ConversionFailure)
            };
            Ok((convert(self.re)?, convert(self.im)?))
        }
        /// Lossy conversion; returns `NaN` components if the conversion fails.
        pub fn to_f64(self) -> (f64, f64) {
            self.try_to_f64().unwrap_or((f64::NAN, f64::NAN))
        }
        pub fn from_f64(re: f64, im: f64) -> Self {
            Complex256 {
                re: Float256::from(re),
                im: Float256::from(im),
            }
        }
        pub fn try_to_complex64(self) -> Result<Complex64, ParseError> {
            let (re, im) = self.try_to_f64()?;
            Ok(Complex64::new(re, im))
        }
        /// Lossy conversion; returns `NaN` components if the conversion fails.
        pub fn to_complex64(self) -> Complex64 {
            let (re, im) = self.to_f64();
            Complex64::new(re, im)
        }
        pub fn from_complex64(c: Complex64) -> Self {
            Complex256 {
                re: Float256::from(c.re),
                im: Float256::from(c.im),
            }
        }
        #[allow(clippy::inherent_to_string_shadow_display)]
        pub fn from_string(s: &str) -> Option<Self> {
            s.parse().ok()
        }

        pub fn checked_add(self, other: Self) -> Option<Self> {
            let re = self.re + other.re;
            let im = self.im + other.im;

            if !re.is_finite() || !im.is_finite() {
                return None;
            }

            Some(Complex256 { re, im })
        }

        pub fn checked_sub(self, other: Self) -> Option<Self> {
            let re = self.re - other.re;
            let im = self.im - other.im;

            if !re.is_finite() || !im.is_finite() {
                return None;
            }

            Some(Complex256 { re, im })
        }

        fn zero() -> Self {
            Complex256::new(Float256::from(0.0), Float256::from(0.0))
        }

        fn is_exact_zero(self) -> bool {
            self.re == Float256::from(0.0) && self.im == Float256::from(0.0)
        }

        /// Validates a component pair: infinities map to `Overflow`, `NaN` to `Pole`.
        fn finite_or_err(re: Float256, im: Float256) -> Result<Self, MathError> {
            if re.is_nan() || im.is_nan() {
                Err(MathError::Pole)
            } else if !re.is_finite() || !im.is_finite() {
                Err(MathError::Overflow)
            } else {
                Ok(Complex256 { re, im })
            }
        }

        /// Like `finite_or_err`, but a zero result from non-zero inputs is `Underflow`.
        fn nonzero_result(self, re: Float256, im: Float256) -> Result<Self, MathError> {
            let out = Self::finite_or_err(re, im)?;
            if out.is_exact_zero() && !self.is_exact_zero() {
                return Err(MathError::Underflow);
            }
            Ok(out)
        }

        pub fn try_add(self, other: Self) -> Result<Self, MathError> {
            Self::finite_or_err(self.re + other.re, self.im + other.im)
        }

        pub fn try_sub(self, other: Self) -> Result<Self, MathError> {
            Self::finite_or_err(self.re - other.re, self.im - other.im)
        }

        pub fn try_mul(self, other: Self) -> Result<Self, MathError> {
            let terms = [
                self.re * other.re,
                self.im * other.im,
                self.re * other.im,
                self.im * other.re,
            ];
            if terms.iter().any(|t| t.is_infinite()) {
                return Err(MathError::Overflow);
            }
            let re = terms[0] - terms[1];
            let im = terms[2] + terms[3];
            if other.is_exact_zero() || self.is_exact_zero() {
                return Self::finite_or_err(re, im);
            }
            Self::finite_or_err(re, im)?;
            if re == Float256::from(0.0) && im == Float256::from(0.0) {
                return Err(MathError::Underflow);
            }
            Ok(Complex256 { re, im })
        }

        pub fn try_div(self, other: Self) -> Result<Self, MathError> {
            if other.is_exact_zero() {
                return Err(MathError::DivisionByZero);
            }
            let denom = other.re * other.re + other.im * other.im;
            if !denom.is_finite() {
                return Err(MathError::Overflow);
            }
            if denom == Float256::from(0.0) {
                return Err(MathError::Underflow);
            }
            self.nonzero_result(
                (self.re * other.re + self.im * other.im) / denom,
                (self.im * other.re - self.re * other.im) / denom,
            )
        }

        pub fn try_recip(self) -> Result<Self, MathError> {
            Complex256::new(Float256::from(1.0), Float256::from(0.0)).try_div(self)
        }

        pub fn try_exp(self) -> Result<Self, MathError> {
            let exp_re = self.re.exp();
            if exp_re == Float256::from(0.0) {
                return Err(MathError::Underflow);
            }
            if exp_re.is_infinite() {
                return Err(MathError::Overflow);
            }
            Self::finite_or_err(exp_re * self.im.cos(), exp_re * self.im.sin())
        }

        pub fn try_ln(self) -> Result<Self, MathError> {
            if self.is_exact_zero() {
                return Err(MathError::Pole);
            }
            let mag = self.magnitude();
            if !mag.is_finite() {
                return Err(MathError::Overflow);
            }
            if mag == Float256::from(0.0) {
                return Err(MathError::Underflow);
            }
            Self::finite_or_err(mag.ln(), self.arg())
        }

        pub fn try_powc(self, exp: Self) -> Result<Self, MathError> {
            if self.is_exact_zero() {
                return if exp.is_exact_zero() {
                    Err(MathError::Pole)
                } else if exp.re > Float256::from(0.0) {
                    Ok(Self::zero())
                } else {
                    Err(MathError::DivisionByZero)
                };
            }
            self.try_ln()?.try_mul(exp)?.try_exp()
        }

        pub fn try_powf(self, exp: Float256) -> Result<Self, MathError> {
            self.try_powc(Complex256::new(exp, Float256::from(0.0)))
        }

        pub fn try_powi(self, exp: i32) -> Result<Self, MathError> {
            let one = Complex256::new(Float256::from(1.0), Float256::from(0.0));

            if exp == 0 {
                return Ok(one);
            }

            let mut base = if exp < 0 { self.try_recip()? } else { self };
            let mut n: u32 = exp.unsigned_abs();
            let mut result = one;

            while n > 0 {
                if (n & 1) == 1 {
                    result = result.try_mul(base)?;
                }
                n >>= 1;
                if n > 0 {
                    base = base.try_mul(base)?;
                }
            }

            Ok(result)
        }

        pub fn try_sqrt(self) -> Result<Self, MathError> {
            let mag = self.magnitude();
            if !mag.is_finite() {
                return Err(MathError::Overflow);
            }
            let re_sqrt = ((self.re + mag) / Float256::from(2.0)).sqrt();
            let im_sqrt = ((mag - self.re) / Float256::from(2.0)).sqrt();
            let im = if self.im >= Float256::from(0.0) {
                im_sqrt
            } else {
                -im_sqrt
            };
            self.nonzero_result(re_sqrt, im)
        }

        fn try_exp_i(self, sign: f64) -> Result<Self, MathError> {
            self.try_mul(Complex256::from_f64(0.0, sign))?.try_exp()
        }

        pub fn try_sin(self) -> Result<Self, MathError> {
            let e_i = self.try_exp_i(1.0)?;
            let e_neg_i = self.try_exp_i(-1.0)?;
            e_i.try_sub(e_neg_i)?
                .try_div(Complex256::from_f64(0.0, 2.0))
        }

        pub fn try_cos(self) -> Result<Self, MathError> {
            let e_i = self.try_exp_i(1.0)?;
            let e_neg_i = self.try_exp_i(-1.0)?;
            e_i.try_add(e_neg_i)?
                .try_div(Complex256::from_f64(2.0, 0.0))
        }

        pub fn try_tan(self) -> Result<Self, MathError> {
            self.try_sin()?.try_div(self.try_cos()?)
        }

        pub fn try_cot(self) -> Result<Self, MathError> {
            self.try_cos()?.try_div(self.try_sin()?)
        }
    }

    pub trait ComplexOps {
        fn add(self, other: Self) -> Self;
        fn sub(self, other: Self) -> Self;
        fn mul(self, other: Self) -> Self;
        fn div(self, other: Self) -> Self;
        fn neg(self) -> Self;
        fn abs(self) -> Float256;
        fn arg(self) -> Float256;
        fn conj(self) -> Self;
        fn magnitude(self) -> Float256;
        fn powc(self, exp: Self) -> Self;
        fn powf(self, exp: Float256) -> Self;
        fn powi(self, exp: i32) -> Self;
        fn exp(self) -> Self;
        fn ln(self) -> Self;
        fn log10(self) -> Self;
        fn recip(self) -> Self;
        fn sqrt(self) -> Self;
        fn sin(self) -> Self;
        fn cos(self) -> Self;
        fn tan(self) -> Self;
        fn cot(self) -> Self;
    }

    #[allow(clippy::wrong_self_convention)]
    pub trait ComplexComparisons {
        fn approx_eq(self, other: Self, tol: Float256) -> bool;
        fn is_zero(self, tol: Float256) -> bool;
        fn is_one(self, tol: Float256) -> bool;
        fn is_i(self, tol: Float256) -> bool;
        fn is_real(self, tol: Float256) -> bool;
        fn is_imaginary(self, tol: Float256) -> bool;
    }

    fn has_valid_tolerance(tol: Float256) -> bool {
        tol.is_finite() && tol >= Float256::from(0.0)
    }

    impl ComplexOps for Complex256 {
        fn add(self, other: Self) -> Self {
            Complex256 {
                re: self.re + other.re,
                im: self.im + other.im,
            }
        }

        fn sub(self, other: Self) -> Self {
            Complex256 {
                re: self.re - other.re,
                im: self.im - other.im,
            }
        }

        fn mul(self, other: Self) -> Self {
            Complex256 {
                re: self.re * other.re - self.im * other.im,
                im: self.re * other.im + self.im * other.re,
            }
        }

        fn div(self, other: Self) -> Self {
            let denom = other.re * other.re + other.im * other.im;
            Complex256 {
                re: (self.re * other.re + self.im * other.im) / denom,
                im: (self.im * other.re - self.re * other.im) / denom,
            }
        }

        fn neg(self) -> Self {
            let re = -self.re;
            let im = -self.im;
            Complex256 { re, im }
        }

        fn abs(self) -> Float256 {
            (self.re * self.re + self.im * self.im).sqrt()
        }

        fn arg(self) -> Float256 {
            // f256's atan2 returns the wrong angle when |im| == |re| (non-zero).
            if self.re.abs() == self.im.abs() && self.re != Float256::from(0.0) {
                let quarter = f256::consts::PI / Float256::from(4.0);
                let angle = if self.re > Float256::from(0.0) {
                    quarter
                } else {
                    quarter * Float256::from(3.0)
                };
                return if self.im < Float256::from(0.0) {
                    -angle
                } else {
                    angle
                };
            }
            self.im.atan2(&self.re)
        }

        fn conj(self) -> Self {
            Complex256 {
                re: self.re,
                im: -self.im,
            }
        }

        fn magnitude(self) -> Float256 {
            (self.re * self.re + self.im * self.im).sqrt()
        }

        fn powc(self, exp: Self) -> Self {
            let ln_self = self.ln();
            let prod = ln_self.mul(exp);
            prod.exp()
        }

        fn powf(self, exp: Float256) -> Self {
            let exp_complex = Complex256 {
                re: exp,
                im: Float256::from(0.0),
            };
            self.powc(exp_complex)
        }

        fn powi(self, exp: i32) -> Self {
            self.try_powi(exp)
                .unwrap_or_else(|e| panic!("Complex256::powi failed: {}", e.message()))
        }

        fn exp(self) -> Self {
            let exp_re = self.re.exp();
            Complex256 {
                re: exp_re * self.im.cos(),
                im: exp_re * self.im.sin(),
            }
        }

        fn ln(self) -> Self {
            Complex256 {
                re: self.magnitude().ln(),
                im: self.arg(),
            }
        }

        fn log10(self) -> Self {
            self.ln().div(Complex256 {
                re: Float256::from(10.0_f64.ln()),
                im: Float256::from(0.0),
            })
        }

        fn recip(self) -> Self {
            let denom = self.re * self.re + self.im * self.im;
            Complex256 {
                re: self.re / denom,
                im: -self.im / denom,
            }
        }

        fn sqrt(self) -> Self {
            let mag = self.magnitude();
            let re_sqrt = ((self.re + mag) / Float256::from(2.0)).sqrt();
            let im_sqrt = ((mag - self.re) / Float256::from(2.0)).sqrt();
            Complex256 {
                re: re_sqrt,
                im: if self.im >= Float256::from(0.0) {
                    im_sqrt
                } else {
                    -im_sqrt
                },
            }
        }

        fn sin(self) -> Self {
            let e_i = self
                .mul(Complex256 {
                    re: Float256::from(0.0),
                    im: Float256::from(1.0),
                })
                .exp();
            let e_neg_i = self
                .mul(Complex256 {
                    re: Float256::from(0.0),
                    im: Float256::from(-1.0),
                })
                .exp();
            e_i.sub(e_neg_i).div(Complex256 {
                re: Float256::from(0.0),
                im: Float256::from(2.0),
            })
        }

        fn cos(self) -> Self {
            let e_i = self
                .mul(Complex256 {
                    re: Float256::from(0.0),
                    im: Float256::from(1.0),
                })
                .exp();
            let e_neg_i = self
                .mul(Complex256 {
                    re: Float256::from(0.0),
                    im: Float256::from(-1.0),
                })
                .exp();
            e_i.add(e_neg_i).div(Complex256 {
                re: Float256::from(2.0),
                im: Float256::from(0.0),
            })
        }

        fn tan(self) -> Self {
            self.sin().div(self.cos())
        }
        fn cot(self) -> Self {
            self.cos().div(self.sin())
        }
    }

    impl ComplexComparisons for Complex256 {
        fn approx_eq(self, other: Self, tol: Float256) -> bool {
            if !has_valid_tolerance(tol) {
                return false;
            }

            let re_diff = (self.re - other.re).abs();
            let im_diff = (self.im - other.im).abs();

            re_diff <= tol && im_diff <= tol
        }

        fn is_zero(self, tol: Float256) -> bool {
            self.approx_eq(
                Complex256 {
                    re: Float256::from(0.0),
                    im: Float256::from(0.0),
                },
                tol,
            )
        }

        fn is_one(self, tol: Float256) -> bool {
            self.approx_eq(
                Complex256 {
                    re: Float256::from(1.0),
                    im: Float256::from(0.0),
                },
                tol,
            )
        }

        fn is_i(self, tol: Float256) -> bool {
            self.approx_eq(
                Complex256 {
                    re: Float256::from(0.0),
                    im: Float256::from(1.0),
                },
                tol,
            )
        }

        fn is_real(self, tol: Float256) -> bool {
            has_valid_tolerance(tol) && self.im.abs() <= tol
        }

        fn is_imaginary(self, tol: Float256) -> bool {
            has_valid_tolerance(tol) && self.re.abs() <= tol && self.im.abs() > tol
        }
    }

    impl std::fmt::Display for Complex256 {
        fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            if self.im >= Float256::from(0.0) {
                write!(f, "{} + {}i", self.re, self.im)
            } else {
                write!(f, "{} - {}i", self.re, -self.im)
            }
        }
    }

    impl std::str::FromStr for Complex256 {
        type Err = ParseError;
        fn from_str(s: &str) -> Result<Self, Self::Err> {
            let s = s.trim();
            let s = s
                .strip_prefix('(')
                .and_then(|value| value.strip_suffix(')'))
                .unwrap_or(s);
            let body = s.strip_suffix('i').ok_or(ParseError::InvalidFormat)?.trim();

            let (re_str, im_str, imag_sign) = if let Some((re, im)) = body.rsplit_once(" + ") {
                (re, im, Float256::from(1.0))
            } else if let Some((re, im)) = body.rsplit_once(" - ") {
                (re, im, Float256::from(-1.0))
            } else {
                return Err(ParseError::InvalidFormat);
            };

            let re = re_str
                .trim()
                .parse::<Float256>()
                .map_err(|_| ParseError::InvalidFormat)?;
            let im = im_str
                .trim()
                .parse::<Float256>()
                .map_err(|_| ParseError::InvalidFormat)?
                * imag_sign;
            Ok(Complex256 { re, im })
        }
    }

    #[cfg(test)]
    mod tests {
        use super::{Complex256, ComplexComparisons, ComplexOps};
        use crate::MathError;
        use f256::f256 as Float256;
        use std::str::FromStr;

        fn c256(re: f64, im: f64) -> Complex256 {
            Complex256::from_f64(re, im)
        }

        #[test]
        fn approx_eq_accepts_componentwise_differences_within_tolerance() {
            let left = c256(1.0, -2.0);
            let right = c256(1.000_5, -2.000_5);

            assert!(left.approx_eq(right, Float256::from(0.001)));
        }

        #[test]
        fn approx_eq_rejects_componentwise_differences_outside_tolerance() {
            let left = c256(1.0, -2.0);
            let right = c256(1.002, -1.999);

            assert!(!left.approx_eq(right, Float256::from(0.001)));
        }

        #[test]
        fn approx_eq_rejects_invalid_tolerance() {
            let left = c256(1.0, 2.0);
            let right = c256(1.0, 2.0);

            assert!(!left.approx_eq(right, Float256::from(-0.001)));
            assert!(!left.approx_eq(right, Float256::INFINITY));
            assert!(!left.approx_eq(right, Float256::NAN));
        }

        #[test]
        fn comparison_helpers_match_expected_special_values() {
            let tol = Float256::from(0.001);

            assert!(c256(0.000_5, -0.000_5).is_zero(tol));
            assert!(c256(1.000_5, -0.000_5).is_one(tol));
            assert!(c256(0.000_5, 1.000_5).is_i(tol));
            assert!(c256(2.0, 0.000_5).is_real(tol));
            assert!(c256(0.000_5, 2.0).is_imaginary(tol));
        }

        #[test]
        fn comparison_helpers_reject_non_matching_values() {
            let tol = Float256::from(0.001);

            assert!(!c256(0.0, 0.002).is_zero(tol));
            assert!(!c256(1.0, 0.002).is_one(tol));
            assert!(!c256(0.002, 1.0).is_i(tol));
            assert!(!c256(2.0, 0.002).is_real(tol));
            assert!(!c256(0.000_5, 0.000_5).is_imaginary(tol));
        }

        #[test]
        fn try_div_and_recip_report_division_by_zero() {
            let zero = c256(0.0, 0.0);
            assert_eq!(c256(1.0, 1.0).try_div(zero), Err(MathError::DivisionByZero));
            assert_eq!(zero.try_recip(), Err(MathError::DivisionByZero));
            assert_eq!(zero.try_powi(-1), Err(MathError::DivisionByZero));
            assert!(c256(1.0, 0.0).try_tan().is_ok());
        }

        #[test]
        fn try_ops_report_overflow() {
            let big = Complex256::new(Float256::MAX, Float256::MAX);
            assert_eq!(big.try_add(big), Err(MathError::Overflow));
            assert_eq!(big.try_mul(big), Err(MathError::Overflow));
            assert_eq!(c256(1e6, 0.0).try_exp(), Err(MathError::Overflow));
            assert_eq!(c256(2.0, 0.0).try_powi(i32::MAX), Err(MathError::Overflow));
            assert_eq!(big.try_sqrt(), Err(MathError::Overflow));
        }

        #[test]
        fn try_ops_report_underflow() {
            let tiny = Complex256::new(Float256::MIN_POSITIVE, Float256::from(0.0));
            assert_eq!(tiny.try_mul(tiny), Err(MathError::Underflow));
            assert_eq!(c256(-1e6, 0.0).try_exp(), Err(MathError::Underflow));
        }

        #[test]
        fn try_ops_report_non_finite_inputs_and_poles() {
            let nan = Complex256::new(Float256::NAN, Float256::from(0.0));
            assert_eq!(nan.try_add(c256(1.0, 0.0)), Err(MathError::Pole));
            assert_eq!(c256(0.0, 0.0).try_ln(), Err(MathError::Pole));
            assert_eq!(
                c256(0.0, 0.0).try_powc(c256(0.0, 0.0)),
                Err(MathError::Pole)
            );
        }

        #[test]
        fn try_ops_match_unchecked_results_on_valid_input() {
            let a = c256(1.5, -2.0);
            let b = c256(0.5, 3.0);
            assert_eq!(a.try_add(b), Ok(a.add(b)));
            assert_eq!(a.try_mul(b), Ok(a.mul(b)));
            assert_eq!(a.try_div(b), Ok(a.div(b)));
            assert_eq!(a.try_powi(3), Ok(a.powi(3)));
            assert_eq!(a.try_sqrt(), Ok(a.sqrt()));
            assert_eq!(a.try_sin(), Ok(a.sin()));
        }

        #[test]
        #[should_panic(expected = "Complex256::powi failed")]
        fn powi_still_panics_on_overflow() {
            let _ = c256(2.0, 0.0).powi(i32::MAX);
        }

        #[test]
        fn parsing_round_trips_display_output_for_positive_and_negative_imaginary_parts() {
            let positive = c256(3.5, 4.25);
            let negative = c256(-3.5, -4.25);

            assert_eq!(
                Complex256::from_str(&positive.to_string()).ok(),
                Some(positive)
            );
            assert_eq!(
                Complex256::from_str(&negative.to_string()).ok(),
                Some(negative)
            );
            assert_eq!(
                Complex256::from_string(&format!("({})", negative)),
                Some(negative)
            );
        }
    }
}
