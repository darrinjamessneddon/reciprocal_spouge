pub mod complex {

    use f256::f256 as Float256;
    use num_complex::Complex64;

    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct Complex256 {
        pub re: Float256,
        pub im: Float256,
    }
    impl Complex256 {
        pub fn new(re: Float256, im: Float256) -> Self {
            Complex256 { re, im }
        }
        pub fn to_f64(self) -> (f64, f64) {
            let re_f256 = self.re;
            let im_f256 = self.im;
            let re_f256_str = re_f256.to_string();
            let im_f256_str = im_f256.to_string();
            let re_f64 = re_f256_str.parse::<f64>().unwrap_or(0.0);
            let im_f64 = im_f256_str.parse::<f64>().unwrap_or(0.0);
            (re_f64, im_f64)
        }
        pub fn from_f64(re: f64, im: f64) -> Self {
            Complex256 {
                re: Float256::from(re),
                im: Float256::from(im),
            }
        }
        pub fn to_complex64(self) -> Complex64 {
            let re_f64: f64 = self.re.to_string().parse().unwrap_or(0.0);
            let im_f64: f64 = self.im.to_string().parse().unwrap_or(0.0);
            Complex64::new(re_f64, im_f64)
        }
        pub fn from_complex64(c: Complex64) -> Self {
            Complex256 {
                re: Float256::from(c.re),
                im: Float256::from(c.im),
            }
        }
        #[allow(clippy::inherent_to_string_shadow_display)]
        pub fn from_string(s: &str) -> Option<Self> {
            let s = s.trim();
            if !s.ends_with(')') {
                return None;
            }
            let s = &s[..s.len() - 1];
            let parts: Vec<&str> = s.split('+').collect();
            if parts.len() != 2 {
                return None;
            }
            let re = parts[0].trim().parse::<Float256>().ok()?;
            let im = parts[1]
                .trim()
                .trim_end_matches('i')
                .parse::<Float256>()
                .ok()?;
            Some(Complex256 { re, im })
        }

        fn checked_mul(self, other: Self) -> Option<Self> {
            let re = (self.re * other.re) - (self.im * other.im);
            let im = (self.re * other.im) + (self.im * other.re);

            if !re.is_finite() || !im.is_finite() {
                return None;
            }

            Some(Complex256 { re, im })
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
        fn log(self) -> Self;
        fn log10(self) -> Self;
        fn recip(self) -> Self;
        fn sqrt(self) -> Self;
        fn sin(self) -> Self;
        fn cos(self) -> Self;
        fn tan(self) -> Self;
    }

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
            let one = Complex256 {
                re: Float256::from(1.0),
                im: Float256::from(0.0),
            };

            if exp == 0 {
                return one;
            }

            let mut base = if exp < 0 { self.recip() } else { self };
            let mut n: u32 = exp.unsigned_abs();
            let mut result = one;

            while n > 0 {
                if (n & 1) == 1 {
                    result = result.checked_mul(base).unwrap_or_else(|| {
                        panic!("Complex256::powi overflow/invalid value while multiplying result by base")
                    });
                }

                n >>= 1;

                if n > 0 {
                    base = base.checked_mul(base).unwrap_or_else(|| {
                        panic!("Complex256::powi overflow/invalid value while squaring base")
                    });
                }
            }

            result
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

        fn log(self) -> Self {
            self.ln().div(Complex256 {
                re: Float256::from(10.0_f64.ln()),
                im: Float256::from(0.0),
            })
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
        type Err = String;
        fn from_str(s: &str) -> Result<Self, Self::Err> {
            let s = s.trim();
            if !s.ends_with('i') {
                return Err("Invalid complex number format".to_string());
            }
            let s = &s[..s.len() - 1];
            let parts: Vec<&str> = s.split('+').collect();
            if parts.len() != 2 {
                return Err("Invalid complex number format".to_string());
            }
            let re = parts[0]
                .trim()
                .parse::<Float256>()
                .map_err(|_| "Invalid real part".to_string())?;
            let im = parts[1]
                .trim()
                .parse::<Float256>()
                .map_err(|_| "Invalid imaginary part".to_string())?;
            Ok(Complex256 { re, im })
        }
    }

    #[cfg(test)]
    mod tests {
        use super::{Complex256, ComplexComparisons};
        use f256::f256 as Float256;

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
    }
}
