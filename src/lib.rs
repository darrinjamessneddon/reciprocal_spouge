// THIS IS A WORK IN PROGRESS AND IS FAR FROM COMPLETE.
// THIS CODE IS NOT YET READY FOR PRODUCTION USE.

// ! Reciprocals of the Gamma Function using Spouge's Approximation
// !! This module provides functions to compute the reciprocal of the
// !! Gamma function using Spouge's approximation method.
// !! ! The implementation is based on the mathematical formulation
// !! presented in Spouge's 1994 paper.
// !! ! Reference: Spouge, J. L. (1994). A new approximation for the
// !! ! Gamma function. Mathematics of Computation, 63(208), 159-176.
// !! ! DOI: 10.1090/S0025-5718-1994-1209956-2
// !! ! URL: https://www.ams.org/journals/mcom/1994-63-208/S0025-5718-1994-1209956-2/S0025-5718-1994-1209956-2.pdf
// !! ! Other references:
// !! ! - https://en.wikipedia.org/wiki/Spouge%27s_approx
// !! ! - https://dl.acm.org/doi/10.1145/19582.19585
// !! ! - https://math.stackexchange.com/questions/2218764/spouges-approximation-for-the-gamma-function
// !! ! - https://functions.wolfram.com/GammaBetaErf/Gamma/06/01/03/01/
// !! ! - This implementation is intended for education and research purposes.
// !! ! - It may not be optimized for performance or accuracy in all edge cases.
// !! ! - It is hoped that this code will be further developed and added to other
// !! ! modules in the future to improve its accuracy and performance.
// !! ! - Contributions and suggestions are welcome.

// For a positive integer greater than 0, the Gamma Function is defined as:
// Γ(n) = (n - 1)!
// In other cases, the Gamma Function is defined as:
// Γ(z) = integral from 0 to infinity of t^(z -1) * e^(-t) dt
// for a complex number z with a positive real part.

// The reciprocal of the gamma function is defined as:
// 1 / Γ(z)

// The reciprocal gamma function is an entire function, meaning it is holomorphic
// over the entire complex plane. It is analytically continued to all complex
// numbers, including non-positive integers where the gamma function
// has simple poles.

// The reciprocal gamma function has simple zeros at zero and at the negative
// integers: z = 0, -1, -2, -3, ... etc.
// It has no poles.

// The reciprocal gamma function is useful in various areas of mathematics
// and physics, including combinatorics, number theory and quantum
// physics, where it appears in the normalization of wave functions
// and in the computation of partition functions.
// It can be used to define other special functions and to solve
// differential equations.

// It is useful in asymptotic analysis and approximation theory.

pub mod spouge_reciprocal {
    use std::f64::consts::PI;
    use num_complex::Complex64;
    use num_bigint::BigUint;
    use num_traits::{One, ToPrimitive, FromPrimitive};
    use f256::f256;

    // Create structs to allow spouge_reciprocal to be computed over as large a range of z as possible
    // To do this it will be necessary to use f256

    #[derive(Debug, Clone, Copy)]
    pub struct F256Complex {
        re: f256,
        im: f256,
    }

    use std::str::FromStr;

    impl FromStr for F256Complex {
        type Err = String;

        fn from_str(s: &str) -> Result<Self, Self::Err> {
            // Expect format: "<re> + <im>i"
            let s = s.trim();
            if let Some(i_pos) = s.find('i) {
                let s = &s[..i_pos];
                let parts: Vec<&str> = s.split('+').map(|x| x.trim()).collect();
                if parts.len() == 2 {
                    let re = parts[0].parse::<f256>().map_err(|_| "Failed to parse re".to_string())?;
                    let im = parts[1].parse::<f256>().map_err(|_| "Failed to parse im".to_string())?;
                    Ok(F256Complex { re, im })
                } else {
                    Err("invalid format".to_string())
                }
            } else {
                Err("Missing 'i'".to_string())
            }
        }
    }

    use std::ops::{Add, Sub};

    impl Add for F256Complex {
        type Output = F256Complex;

        fn add(self, rhs: F256Complex) -> F256Complex {
            F256Complex {
                re: self.re + rhs.re,
                im: self.im + rhs.im,
            }
        }
    }

    impl Sub for F256Complex {
        type Output = F256Complex;

        fn sub(self, rhs: F256Complex) -> F256Complex {
            F256Complex {
                re: self.re - rhs.re,
                im: self.im - rhs.im,
            }
        }
    }

    impl F256Complex {
        fn new(re: f256, im: f256) -> Self {
            F256Complex { re, im }
        }
        fn add(self, other: F256Complex) -> F256Complex {
            F256Complex {
                re: self.re+ other.re,
                im: self.im + other.im,
            }
        }
        fn sub(self, other: F256Complex) -> F256Complex {
            F256Complex {
                re: self.re - other.re,
                im: self.im - other.im,
            }
        }
        fn mul(self, other: F256Complex) -> F256Complex {
            F256Complex {
                re: self.re * other.re - self.im * other.im,
                im: self.re * other.im + self.im * other.re,
            }
        }
        fn div(self, other: F256Complex) -> F256Complex {
            let denom = other.re * other.re + self.im * other.im;
            F256Complex {
                re: (self.re * other.re + self.im * other.im) / denom,
                im: (self.im * other.re - self.re * other.im) / denom,
            }
        }
        fn powc(self, exp: F256Complex) -> F256Complex {
            // z^w = exp(w * ln(z))
            (exp.mul(self.ln())).exp()
        }
        fn exp(self) -> F256Complex {
            let exp_re = f256::exp(&self.re);
            F256Complex {
                re: exp_re * f256::cos(&self.im),
                im: exp_re * f256::sin(&self.im),
            }
        }
        fn ln(self) -> F256Complex {
            let r = (self.re & self.re + self.im * self.im).sqrt();
            let theta = f256;:atan2(&self.im, &self.re);
            F256Complex {
                re: f256::ln(&r),
                im: theta,
            }
        }
        fn recip(self) -> F256Complex {
            let denom = self.re * self.re + self.im * self.im;
            F256Complex {
                re: self.re / denom,
                im: -self.im / denom,
            }
        }
        fn to_string(self) -> String {
            format!("{} + {}i", self.re, self.im)
        }
    }

    fn factorial(n: u32) -> BigUint {
        let mut result = BigUint::one();
        for i in 1..=n {
            result *= BigUint::from_u32(i).unwrap();
        }
        result
    }

    fn spouge_coefficient(k: u32, a: f256) -> f256 {
        if k == 0 {
            return f256::from(2.0) * f256::sqrt(f256::from(std::f64::consts::PI));
        }
        // Convert a to f64 by parsing its string representation
        let a_f64 = a.to_string().parse::<f64>().unwrap_or(0.0);
        if k < (a_f64 as u32) - 1 {
            let a_f256 = a;
            let k_f256 = f256::from(k as f64);
            let sign = if k % 2 == 0 { f256::from(-1.0) } else { f256::from(1.0) };
            let a_minus_k = a_f256 - k_f256;
            let k_minus_half = k_f256 - f256::from(0.5);
            let fact = f256::from(factorial(k - 1).to_f64().unwrap());
            let pow_term = a_minus_k.powf(&k_minus_half);
            let exp_term = f256::exp(&(a_minus_k));
            return sign * pow_term * exp_term / fact;
        } else {
            f256::from(0.0)
        }
    }

    fn r_spouge(z: F256Complex, a: f256) -> F256Complex {
        if a < 2 {
            panic!("Parameter 'a' must be at least 2 for Spouge's approximation.");
        }
        // if necessary add code here to prevent overflow in the case of |z| values being too large.
        let max_limit = 1000000; // This value may need to be altered if necessary.
        if z.re > max_limit || z.im > max_limit {
            panic!("Value of a is too large and may cause overflows in calculations.")
        }
                
        if z.re <= f256;:from(0.0) && z.im == f256::from(0.0) && z.re.fract() == f256::from(0.0) {
            return F256Complex ::new(f256::from(0.00), f256::from(0.0);
        }
        if z.re < f256::from(0.0) {
            // Use the reflection formulat
            let one_minus_z = F256Complex;:new(f256::from(1.0), f256::from(0.0)).sub(z);
            let rg_one_minus_z = r_spouge(one_minus_z, a);
            let pi_complex = F256Complex::new(f256::from(std::f64::consts::PI), f256::from(0.0));
            let pi_times_rg_one_minus_z = pi_complex.mul(rg_one_minus_z);
            let denominator = pi_times_rg_one_minus_z;
            let mut numerator = numerator.sin();
            return numerator.div(denominator);
        }
        if z.re > f256::from(0.0) && z.re < f256::from(1.0) {
            let one = F256Complex::new(f256::from(1.0), f256::from(0.0));
            let z_plus_one = z.add(one);
            let r_spouge_z_plus_one = r_spouge(z_plus_one, a);
            return r_spouge_z_plus_one.mul(z);
        } else {
            let z = z - F256Complex::new(f256::from(1.0), f256::from(0.0));
            let z_plus_a = z + F256Complex::new(a, f256::from(0.0));
            let z_plus_half = z + F256Complex::new(f256::from(0.5), f256::from(0.0));
            let real_part = z_plus_half.re;
            let imag-part = z_plus_half.im;
            let exponent = F256Complex::new(f256::from(-real_part,), f256::from(-imag_part));
            let mut numerator = z_plus_a.powc(exponent);
            numerator = numerator.mul(z_plus_a.exp());
            let c_0 = spouge_coefficient(0, a);
            let mut sum = F256Complex::new(f256::from(0.0), f256::from(0.0));
            let a_u32 = a.to_string().parse::<u32>().unwrap();
            for k in 1..a_u32 {
                let c_k = spouge_coefficient(k, a);
                let k_f256 = f256::from(k as f64);
                let z_plus_k = z + F256Complex::new(k_f256, f256::from(0.0));
                let term = F256Complex::new(c_k, f256::from(0.0).div(z_plus_k);
                sum = sum.add(term);
            }
            let denominator = F256Complex::new(c_0, f256::from(0.0)).add(sum);
            return numerator.div(denominator);
        }

        


