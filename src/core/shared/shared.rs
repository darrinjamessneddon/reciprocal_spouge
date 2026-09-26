pub mod shared {
    use core::str::FromStr;

    use ::f256::consts::TAU;
    use f256::f256 as Float256;
    use num_bigint::BigUint;

    use crate::core::MathError;

    pub fn factorial(n: u64) -> BigUint {
        (1..=n).fold(BigUint::from(1_u8), |acc, value| acc * BigUint::from(value))
    }

    /// Compute the Spouge coefficients for a given parameter 'a'. The coefficients are used in Spouge's approximation of the gamma function.
    /// allow for error handling by returning a Result type, which can either be Ok with the coefficients
    /// or Err with a MathError describing the error.
    pub fn spouge_coefficients(a: u64) -> Result<Vec<Float256>, MathError> {
        if a < 2 {
            return Err(MathError::ParameterOutOfRange);
        }

        let sqrt_two_pi = TAU.sqrt();
        let a_f256 = Float256::from(a);
        let mut coefficients = Vec::with_capacity(a as usize);
        coefficients.push(sqrt_two_pi);

        let mut factorial_k_minus_1 = BigUint::from(1_u8);
        for k in 1..a {
            let k_f256 = Float256::from(k);
            let sign = if k % 2 == 0 {
                Float256::from(-1.0)
            } else {
                Float256::from(1.0)
            };
            let fact_f256 = Float256::from_str(&factorial_k_minus_1.to_string())
                .map_err(|_| MathError::Overflow)?;
            if !fact_f256.is_finite() {
                return Err(MathError::Overflow);
            }
            let a_minus_k = a_f256 - k_f256;
            let k_minus_half = k_f256 - Float256::from(0.5);
            let pow_term = a_minus_k.powf(&k_minus_half);
            let exp_term = a_minus_k.exp();
            coefficients.push(sign * pow_term * exp_term / fact_f256);
            factorial_k_minus_1 *= BigUint::from(k);
        }
        Ok(coefficients)
    }
    /// Add structs and functions to facilitate the use of Kahan summation by feeding f256 terms into a Kahan accumulator using Neumaier's algorithm
    #[derive(Debug, Clone, Copy, Default)]
    pub struct KahanF256 {
        pub kahan_sum: Float256,
        pub compensation: Float256,
    }

    impl KahanF256 {
        /// Creates a compensated accumulator with the provided running sum and compensation.
        pub fn new(kahan_sum: Float256, compensation: Float256) -> Self {
            Self {
                kahan_sum,
                compensation,
            }
        }

        /// Code to feed a new f256 term into the accumulator using Neumaier's algorithm.
        pub fn add(&mut self, term: Float256) -> Self {
            let sum = self.kahan_sum + term;
            if self.kahan_sum.abs() >= term.abs() {
                self.compensation += (self.kahan_sum - sum) + term;
            } else {
                self.compensation += (term - sum) + self.kahan_sum;
            }
            self.kahan_sum = sum;
            *self
        }

        pub fn total(&self) -> Float256 {
            self.kahan_sum + self.compensation
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        fn assert_f256_close(actual: Float256, expected: Float256, tolerance: Float256) {
            let difference = (actual - expected).abs();
            assert!(
                difference <= tolerance,
                "expected {actual} to be within {tolerance} of {expected}, diff={difference}"
            );
        }

        #[test]
        fn factorial_handles_zero_and_small_values() {
            assert_eq!(factorial(0), BigUint::from(1_u8));
            assert_eq!(factorial(1), BigUint::from(1_u8));
            assert_eq!(factorial(5), BigUint::from(120_u16));
        }

        #[test]
        fn spouge_coefficients_follow_standard_formula() {
            let a = 12_u64;
            let coefficients_result = spouge_coefficients(a);
            let coefficients = match coefficients_result {
                Ok(value) => value,
                Err(error) => panic!("unexpected coefficient error: {error:?}"),
            };
            assert_eq!(coefficients.len(), a as usize);
            assert_f256_close(coefficients[0], TAU.sqrt(), Float256::from(1e-30));

            let mut factorial_k_minus_1 = BigUint::from(1_u8);
            for k in 1..a {
                let k_f256 = Float256::from(k);
                let a_minus_k = Float256::from(a) - k_f256;
                let k_minus_half = k_f256 - Float256::from(0.5);
                let fact_f256_result = Float256::from_str(&factorial_k_minus_1.to_string());
                let fact_f256 = match fact_f256_result {
                    Ok(value) => value,
                    Err(parse_error) => panic!("factorial parse failed: {parse_error:?}"),
                };
                let numerator = coefficients[k as usize] * fact_f256;
                let denominator = a_minus_k.powf(&k_minus_half) * a_minus_k.exp();
                let signed_ratio = numerator / denominator;
                let expected_sign = if k % 2 == 0 {
                    Float256::from(-1.0)
                } else {
                    Float256::from(1.0)
                };
                assert_f256_close(signed_ratio, expected_sign, Float256::from(1e-24));
                factorial_k_minus_1 *= BigUint::from(k);
            }
        }

        #[test]
        fn neumaier_accumulator_recovers_lost_low_order_term() {
            let large = Float256::from(1e100);
            let naive = (large + Float256::from(1.0)) + (-large);
            let mut accumulator = KahanF256::default();
            accumulator.add(large);
            accumulator.add(Float256::from(1.0));
            accumulator.add(-large);

            assert_f256_close(naive, Float256::from(0.0), Float256::from(0.0));
            assert_f256_close(
                accumulator.total(),
                Float256::from(1.0),
                Float256::from(0.0),
            );
        }
    }
}
