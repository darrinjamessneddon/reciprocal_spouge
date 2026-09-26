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

    use astro_float::{BigFloat, Consts, RoundingMode, ctx::Context};
    use rayon::prelude::*;

    pub const TARGET_PRECISION: usize = 476; // TO match a true 512-bit precision
    pub const WORKING_PRECISION: usize = TARGET_PRECISION + 128; // Total 604 bits
    pub const GUARD_DIGITS: usize = 128; // Number of guard digits added to the target precision

    pub fn create_context() -> Context {
        let consts = match Consts::new() {
            Ok(value) => value,
            Err(error) => panic!("failed to initialize astro-float constants for context: {error}"),
        };
        Context::new(
            WORKING_PRECISION,
            RoundingMode::ToEven,
            consts,
            -(GUARD_DIGITS as i32),
            GUARD_DIGITS as i32,
        )
    }

    pub fn create_consts(ctx: &Context) -> Consts {
        let _ = ctx;
        match Consts::new() {
            Ok(value) => value,
            Err(error) => panic!("failed to initialize astro-float constants: {error}"),
        }
    }

    pub fn create_context_and_consts() -> (Context, Consts) {
        let ctx = create_context();
        let consts = create_consts(&ctx);
        (ctx, consts)
    }

    pub fn precomputed_bigfloat_spouge_coefficients(a: usize) -> Vec<BigFloat> {
        let ctx = create_context();
        let mut consts = create_consts(&ctx);
        // Compute tau as 2pi using the mathematical constants object
        let two = BigFloat::from(2u32);
        let pi = consts.pi(WORKING_PRECISION, RoundingMode::ToEven);
        let tau_lib = two.mul(&pi, WORKING_PRECISION, RoundingMode::ToEven);
        // Compute the square root as exp(log(x) / 2), since this version of astro-float
        // does not provide a general `pow` method on `Context`.
        let tau_log = tau_lib.ln(WORKING_PRECISION, RoundingMode::ToEven, &mut consts);
        let tau_log_half = tau_log.mul(
            &BigFloat::from(0.5),
            WORKING_PRECISION,
            RoundingMode::ToEven,
        );
        let tau_sqrt = tau_log_half.exp(WORKING_PRECISION, RoundingMode::ToEven, &mut consts);
        let coefficients: Vec<BigFloat> = (0..a)
            .into_par_iter()
            .map(|k| {
                let mut local_consts = match Consts::new() {
                    Ok(value) => value,
                    Err(error) => {
                        panic!("failed to initialize local astro-float constants: {error}")
                    }
                };
                if k == 0 {
                    tau_sqrt.clone()
                } else {
                    let sign = if k % 2 == 0 {
                        -BigFloat::from(1u32)
                    } else {
                        BigFloat::from(1u32)
                    };
                    let k_big = BigFloat::from(k as u32);
                    let a_big = BigFloat::from(a as u32);
                    let a_minus_k = a_big.sub(&k_big, WORKING_PRECISION, RoundingMode::ToEven);
                    let k_big_minus_half = k_big.sub(
                        &BigFloat::from(0.5),
                        WORKING_PRECISION,
                        RoundingMode::ToEven,
                    );
                    let a_minus_k_ln =
                        a_minus_k.ln(WORKING_PRECISION, RoundingMode::ToEven, &mut local_consts);
                    let pow_exponent = a_minus_k_ln.mul(
                        &k_big_minus_half,
                        WORKING_PRECISION,
                        RoundingMode::ToEven,
                    );
                    let pow_term = pow_exponent.exp(
                        WORKING_PRECISION,
                        RoundingMode::ToEven,
                        &mut local_consts,
                    );
                    let exp_term =
                        a_minus_k.exp(WORKING_PRECISION, RoundingMode::ToEven, &mut local_consts);
                    let mut fact_big = BigFloat::from(1u32);
                    for i in 2..=(k - 1) {
                        fact_big = fact_big.mul(
                            &BigFloat::from(i as u32),
                            WORKING_PRECISION,
                            RoundingMode::ToEven,
                        );
                    }
                    let product = pow_term.mul(&exp_term, WORKING_PRECISION, RoundingMode::ToEven);
                    let numerator = sign.mul(&product, WORKING_PRECISION, RoundingMode::ToEven);
                    numerator.div(&fact_big, WORKING_PRECISION, RoundingMode::ToEven)
                }
            })
            .collect();
        coefficients
    }

    pub fn precomputed_coefficients(a: usize) -> Result<Vec<Float256>, MathError> {
        let coefficients = precomputed_bigfloat_spouge_coefficients(a);
        coefficients
            .iter()
            .map(|c| Float256::from_str(&c.to_string()).map_err(|_| MathError::Overflow))
            .collect()
    }
    // For testing processes we should be able to compute these precomputed coefficients as a string
    pub fn precomputed_coefficients_str(a: usize) -> Result<Vec<String>, MathError> {
        let result = precomputed_coefficients(a)?;
        Ok(result.iter().map(|c| c.to_string()).collect())
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
