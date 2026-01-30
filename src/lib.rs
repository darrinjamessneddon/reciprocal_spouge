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

use num_bigint::BigUint;
use num_traits::{One, FromPrimitive, ToPrimitive};

pub mod spouge_reciprocal {
    use f256::f256;

    // Create structs to allow spouge_reciprocal to be computed over as large a range of z as possible
    // To do this it will be necessary to use f256

    #[derive(Debug, Clone, Copy)]
    pub struct F256Complex {
        pub re: f256,
        pub im: f256,
    }

    // Move trait implementations and methods inside the module so F256Complex is in scope

    use std::str::FromStr;
    use std::ops::{Add, Sub, Neg};

    impl FromStr for F256Complex {
        type Err = String;

        fn from_str(s: &str) -> Result<Self, Self::Err> {
            // Expect format: "<re> + <im>i"
            let s = s.trim();
            if let Some(i_pos) = s.find('i') {
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

    impl Neg for F256Complex {
        type Output = F256Complex;

        fn neg(self) -> F256Complex {
            F256Complex {
                re: -self.re,
                im: -self.im,
            }
        }
    }

    impl F256Complex {
        pub fn new(re: f256, im: f256) -> Self {
            F256Complex { re, im }
        }
        pub fn add(self, other: F256Complex) -> F256Complex {
            F256Complex {
                re: self.re + other.re,
                im: self.im + other.im,
            }
        }
        pub fn sub(self, other: F256Complex) -> F256Complex {
            F256Complex {
                re: self.re - other.re,
                im: self.im - other.im,
            }
        }
        pub fn mul(self, other: F256Complex) -> F256Complex {
            F256Complex {
                re: self.re * other.re - self.im * other.im,
                im: self.re * other.im + self.im * other.re,
            }
        }
        pub fn div(self, other: F256Complex) -> F256Complex {
            let denom = other.re * other.re + self.im * other.im;
            F256Complex {
                re: (self.re * other.re + self.im * other.im) / denom,
                im: (self.im * other.re - self.re * other.im) / denom,
            }
        }
        pub fn powc(self, exp: F256Complex) -> F256Complex {
            // z^w = exp(w * ln(z))
            (exp.mul(self.ln())).exp()
        }
        pub fn exp(self) -> F256Complex {
            let exp_re = f256::exp(&self.re);
            F256Complex {
                re: exp_re * f256::cos(&self.im),
                im: exp_re * f256::sin(&self.im),
            }
        }
        pub fn ln(self) -> F256Complex {
            let r = (self.re * self.re + self.im * self.im).sqrt();
            let theta = f256::atan2(&self.im, &self.re);
            F256Complex {
                re: f256::ln(&r),
                im: theta,
            }
        }
        #[allow(dead_code)]
        fn recip(self) -> F256Complex {
            let denom = self.re * self.re + self.im * self.im;
            F256Complex {
                re: self.re / denom,
                im: -self.im / denom,
            }
        }
        pub fn sin(self) -> F256Complex {
            // sin(a + bi) = sin(a)cosh(b) + i*cos(a)sinh(b)
            let sin_re = f256::sin(&self.re);
            
            let cos_re = f256::cos(&self.re);
            
            // cosh(x) = (exp(x) + exp(-x)) / 2
            fn cosh(x: &f256) -> f256 {
                (f256::exp(x) + f256::exp(&(-*x))) / f256::from(2.0)
            }
            // sinh(x) = (exp(x) - exp(-x) / 2
            fn sinh(x: &f256) -> f256 {
                (f256::exp(x) - f256::exp(&(-*x))) / f256::from(2.0)
            }

            let sinh_im = sinh(&self.im);
            let cosh_im = cosh(&self.im);
            F256Complex {
                re: sin_re * cosh_im,
                im: cos_re * sinh_im,
            }
        }
        #[allow(dead_code)]
        fn to_string(self) -> String {
            format!("{} + {}i", self.re, self.im)
        }
    }

    } // end of mod spouge_reciprocal

    use crate::spouge_reciprocal::F256Complex;

    #[allow(dead_code)]
    fn factorial(n: u32) -> BigUint {
        let mut result = BigUint::one();
        for i in 1..=n {
            result *= BigUint::from_u32(i).unwrap();
        }
        result
    }

    // Helper function: sqrt for f256 using Newton's method
    #[allow(dead_code)]
    fn f256_sqrt(x: f256::f256) -> f256::f256 {
        if x <= f256::f256::from(0.0) {
            return f256::f256::from(0.0);
        }
        let mut guess = x / f256::f256::from(2.0);
        for _ in 0..20 {
            guess = (guess + x / guess) / f256::f256::from(2.0);
        }
        guess
    }

    #[allow(dead_code)]
    fn spouge_coefficient(k: u32, a: f256::f256) -> f256::f256 {
        if k == 0 {
            return f256_sqrt(f256::f256::from(2.0)) * f256_sqrt(f256::f256::from(std::f64::consts::PI));
        }
        // Convert a to f64 by parsing its string representation
        let a_f64 = a.to_string().parse::<f64>().unwrap_or(0.0);
        if k < (a_f64 as u32) - 1 {
            let a_f256 = a;
            let k_f256 = f256::f256::from(k as f64);
            let sign = if k % 2 == 0 { f256::f256::from(-1.0) } else { f256::f256::from(1.0) };
            let a_minus_k = a_f256 - k_f256;
            let k_minus_half = k_f256 - f256::f256::from(0.5);
            let fact = f256::f256::from(factorial(k - 1).to_f64().unwrap());
            let pow_term = a_minus_k.powf(&k_minus_half);
            let exp_term = a_minus_k.exp();
            return sign * pow_term * exp_term / fact;
        } else {
            f256::f256::from(0.0)
        }
    }

    pub fn r_spouge(z: F256Complex, a: f256::f256) -> F256Complex {
        if a < f256::f256::from(2.0) {
            panic!("Parameter 'a' must be at least 2 for Spouge's approximation.");
        }
        // if necessary add code here to prevent overflow in the case of |z| values being too large.
        let max_limit = f256::f256::from(10000.0); // This value may need to be altered if necessary.
        if z.re > max_limit || z.im > max_limit {
            panic!("Value of a is too large and may cause overflows in calculations.")
        }
                
        if z.re <= f256::f256::from(0.0) && z.im == f256::f256::from(0.0) && z.re.fract() == f256::f256::from(0.0) {
            return F256Complex::new(f256::f256::from(0.0), f256::f256::from(0.0));
        }
        if z.re < f256::f256::from(0.0) {
            // Use the reflection formula:
            // gamma(1 - z) = pi / (sin(pi * z) * gamma(z))
            // first calculate numerator: pi / sin(pi * z)
            let pi_complex = F256Complex::new(f256::f256::from(std::f64::consts::PI), f256::f256::from(0.0));
            let pi_times_z = pi_complex.mul(z);
            let numerator = pi_complex.div(pi_times_z.sin());
            // next calculate denominator: r_spouge(1 - z, a)

            let one_minus_z = F256Complex::new(f256::f256::from(1.0), f256::f256::from(0.0)).sub(z);
            let rg_one_minus_z = r_spouge(one_minus_z, a);
            let pi_complex = F256Complex::new(f256::f256::from(std::f64::consts::PI), f256::f256::from(0.0));
            let pi_times_rg_one_minus_z = pi_complex.mul(rg_one_minus_z);
            let denominator = pi_times_rg_one_minus_z;
            let numerator = numerator.sin();
            return numerator.div(denominator);
        }
        if z.re > f256::f256::from(0.0) && z.re < f256::f256::from(1.0) {
            let one = F256Complex::new(f256::f256::from(1.0), f256::f256::from(0.0));
            let z_plus_one = z.add(one);
            let r_spouge_z_plus_one = r_spouge(z_plus_one, a);
            return r_spouge_z_plus_one.mul(z);
        } else {
            let z = z - F256Complex::new(f256::f256::from(1.0), f256::f256::from(0.0));
            let z_plus_a = z + F256Complex::new(a, f256::f256::from(0.0));
            let z_plus_half = z + F256Complex::new(f256::f256::from(0.5), f256::f256::from(0.0));
            let real_part = z_plus_half.re;
            let imag_part = z_plus_half.im;
            let exponent = F256Complex::new(f256::f256::from(-real_part), f256::f256::from(-imag_part));
            let mut numerator = z_plus_a.powc(exponent);
            numerator = numerator.mul(z_plus_a.exp());
            let _c_0 = spouge_coefficient(0, a);
            let mut sum = F256Complex::new(f256::f256::from(0.0), f256::f256::from(0.0));
            let a_u32 = a.to_string().parse::<u32>().unwrap();
            for k in 1..a_u32 {
                let c_k = spouge_coefficient(k, a);
                let k_f256 = f256::f256::from(k as f64);
                let z_plus_k = z + F256Complex::new(k_f256, f256::f256::from(0.0));
                let term = F256Complex::new(c_k, f256::f256::from(0.0)).div(z_plus_k);
                sum = sum.add(term);
            }
            let denominator = F256Complex::new(_c_0, f256::f256::from(0.0)).add(sum);
                return numerator.div(denominator);
            }
        }
    
    
        #[allow(dead_code)]
        fn spouge(z: F256Complex, a: f256::f256) -> F256Complex {
            if a < f256::f256::from(2.0) {
                panic!("Parameter 'a' must be at least 2 for Spouge's approximation");
            }
            // code block to prevent overflow in the case of |z| being too large
            let max_limit = f256::f256::from(10000.0); // This value may need to be altered if necessary
            if z.re > max_limit || z.im > max_limit {
                panic!("Absolute value of z is too large and may cause overflows in calculation");
            }
            // Handle the pole at z = 0
            if z.re == f256::f256::from(0.0) && z.im == f256::f256::from(0.0) {
                return F256Complex::new(f256::f256::INFINITY, f256::f256::from(0.0))
            }
            // Handle the pole at negative integers
            if z.re < f256::f256::from(0.0) && z.im == f256::f256::from(0.0) {
                if z.re.fract() == f256::f256::from(0.0) {
                    return F256Complex::new(f256::f256::INFINITY, f256::f256::from(0.0))
                }
            }
            // Handle the case for negative non-integers
            if z.re < f256::f256::from(0.0) && z.re.fract() != f256::f256::from(0.0) {
                // Use the reflection formula:
                // gamma(1 - z)* gamma(z) = pi/ sin(pi * z)
                // first calculate gamma(1 - z):
                let one_minus_z = F256Complex::new(f256::f256::from(1.0), f256::f256::from(0.0)).sub(z);
                let denominator = spouge(one_minus_z, a);
                // next calculate pi/(sin(pi * z))
                let pi_complex = F256Complex::new(f256::f256::from(std::f64::consts::PI), f256::f256::from(0.0));
                let pi_times_z = pi_complex.mul(z);
                let numerator = pi_complex.div(pi_times_z);
                return numerator.div(denominator)
            }
            // Handle the case for z values with real z between 0 and 1:
            if z.re > f256::f256::from(0.0) && z.re < f256::f256::from(1.0) {
                // Use the relationsip: gamma(z + 1) = z * gamma(z);
                let one = F256Complex::new(f256::f256::from(1.0), f256::f256::from(0.0));
                let z_plus_one = z.add(one);
                let numerator = spouge(z_plus_one, a);
                let denominator = z;
                return numerator.div(denominator);           
            } else {
                // Handle the case for other z values by implementing Spouge's approximation
                let z = z - F256Complex::new(f256::f256::from(1.0), f256::f256::from(0.0));
                let z_plus_a = z + F256Complex::new(a, f256::f256::from(0.0));
                let z_plus_half = z + F256Complex::new(f256::f256::from(0.5), f256::f256::from(0.0));
                let real_part = z_plus_half.re;
                let imag_part = z_plus_half.im;
                let exponent = F256Complex::new(real_part, imag_part);
                let mut prefactor = z_plus_a.powc(exponent);
                prefactor = prefactor.mul((-z_plus_a).exp());
                let c_0 = spouge_coefficient(0, a);
                let mut sum = F256Complex::new(f256::f256::from(0.0), f256::f256::from(0.0));
                let a_u32 = a.to_string().parse::<u32>().unwrap();
                for k in 1..a_u32 {
                    let c_k = spouge_coefficient(k, a);
                    let k_f256 = f256::f256::from(k as f64);
                    let z_plus_k = z + F256Complex::new(k_f256, f256::f256::from(0.0));
                    let term = F256Complex::new(c_k, f256::f256::from(0.0)).div(z_plus_k);
                    sum = sum.add(term);
                }
                let sum = F256Complex::new(c_0, f256::f256::from(0.0)).add(sum);
                let factor = prefactor.mul(sum);
                return factor
            }
        }
    
        #[allow(dead_code)]
        fn ln_spouge(z: F256Complex, a: f256::f256) -> F256Complex {
            let spouge_gamma = spouge(z, a);
            let ln_spouge_gamma = spouge_gamma.ln();
            ln_spouge_gamma
        }




            

        


