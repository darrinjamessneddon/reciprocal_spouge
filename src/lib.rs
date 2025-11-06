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

    #[derive(Debug, Clone, Copy)]

    pub struct RSpouge {
        pub a: usize,// Spouge's parameter 'a'
        pub z: Complex64, // Input value 'z'
    }
    impl RSpouge {
        // Constructor for RSpouge
        pub fn new(z: Complex64, a: usize) -> Self {
            if a < 2 {
                panic!("Parameter 'a' must be at least 2 for Spouge's approximation.");
            }
            RSpouge { z, a }
        }
        pub fn factorial(n: usize) -> BigUint {
            let mut result = BigUint::one();
            for i in 2..=n {
                result *= BigUint::from_usize(i).unwrap();
            }
            result
        }
        pub fn compute(&self) -> Complex64 {
            // An adaptation of Spouge's approximation to compute
            // the reciprocal of the Gamma function.
            // Handle the zeros at negative integers and zero
            if self.z.re <= 0.0 && self.z.im == 0.0 && self.z.re.fract() == 0.0 {
                return Complex64::new(0.0, 0.0);
            }
            // Use the reflection formula for negative real parts
            if self.z.re < 0.0 {
                // Using the reflection formula:
                // Γ(z) * Γ(1 - z) =  π / sin(πz)
                // first compute Γ(1 - z)
                let one_minus_z = Complex64::new(1.0, 0.0) - self.z;
                let spouge_rg_one_minus_z = RSpouge::new(one_minus_z, self.a);
                // next compute sin(pi * z)
                let pi_z = Complex64::new(PI * self.z.re, PI * self.z.im);
                let sin_pi_z = pi_z.sin();
                let sin_pi_z_over_pi = sin_pi_z / Complex64::new(PI, 0.0);
                // finally compute the reciprocal gamma using the reflection formula
                return sin_pi_z_over_pi * spouge_rg_one_minus_z.compute();
            }
            // If 0 < Re(z) < 1, use the property Γ(z) = Γ(z + 1) / z
            if self.z.re > 0.0 && self.z.re < 1.0 {
                let spouge_rg_shifted = RSpouge::new(self.z + Complex64::new(1.0, 0.0), self.a);
                return spouge_rg_shifted.compute() * self.z;
            }
            // Spouge's approximation for Re(z) >= 1
            let a = self.a as f64;
            let z = self.z;
            let z_plus_a = z + Complex64::new(a, 0.0);
            let mut numerator = z_plus_a.powc(-(z + Complex64::new(0.5, 0.0)));
            numerator *= z_plus_a.exp();
            let c_0 = Complex64::new((2.0 * PI).sqrt(), 0.0);
            let mut sum = Complex64::new(0.0, 0.0);
            for k in 1..self.a {
                let k_f64 = k as f64;
                let c_k = ((-1.0f64).powf(k_f64 - 1.0)
                / RSpouge::factorial(k - 1).to_f64().unwrap())
                * (a - k_f64).powf(k_f64 - 0.5)
                * (a - k_f64).exp();
            let term = c_k / (z + Complex64::new(k_f64, 0.0));
            sum += term;
            }
            let denominator = c_0 + sum;
            numerator / denominator
        }
    }
}

