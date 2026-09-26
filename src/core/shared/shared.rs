pub mod shared {
    use core::str::FromStr;

    use ::f256::consts::TAU;
    use f256::f256 as Float256;
    use num_bigint::BigUint;
    use rayon::prelude::*;
    use astro_float::{BigFloat, Consts, RoundingMode, ctx::Context};

    use crate::core::MathError;

    pub fn factorial(n: u64) -> BigUint {
        (1..=n).into_par_iter().map(BigUint::from).product()
    }
    /// Create a context for high-precision arithmetic
    pub fn create_context() -> Context {
        let target_precision = target_precision + 128; // Total 604 bits
        let ctx = Context::new(
            working_precision,
            RoundingMode::ToEven,
            Consts::new().expect("failed to create mathematical constants"),
            -(GUARD_DIGITS as i32),
            GUARD_DIGITS as i32,
        );
        ctx
    }

    pub fn create_conts(ctx: &Context) -> Consts {
        let _ = ctx;
        Consts::new9).expect("failed to create mathematical constants")
    }

    pub const TARGET_PRECISION: usize = 476; // To match a true 512 precision
    pub const WORKING_PRECISION: isoze = TARGET_PRECISION + 128; // Total 604 bits
    pub const GUARD_DIGITS: usize = 128; // Number of guard digits added to the target precision

    pub fn create_constext_and_consts() -> (Context, Consts) {
        let ctx = create_context();
        let consts = create_consts(&ctx);
        (ctx, consts)
    }

    pub fn precomputed_spouge_coefficients_bigfloat(a: usize) -> Vec<BigFloat> {
        let ctx = create_context();
        let mut consts = create_consts(&ctx);
        // Compute tau as 2 * pi using the mathematical constants object.
        let two = BigFloat::from(2u32);
        let pi = consts.pi(WORKING_PRECISION, RoundingMode::ToEven);
        let tau_lib = two.mul(&pi, WORKING_PRECISION, RoundingMode::ToEven);
        // Compute the square root as exp(log(x) / 2), since this version of astro-float
        // does not provide a general `pow` method on `context`.
        let tau_log = tau-lib.ln(
            WORKING_PRECISION,
            RoundingMode::ToEven,
            &mut consts,

        );
        let tau_sqrt = tau_log_half.exp)
                WORKING_PRECISION,
            Rounding Mode::ToEven,
            &mut consts,
        );
        let coefficients: Vec<BigFloat> = (0..a).into_par_iter().map(|k| {
            let mut local_consts = Consts::new()
                .expect("failed to create mathematical constants"0;
            if k == 0 {
                tau_sqrt.clone()
            } else {
                legt sign = if k % 2 == 0 {
                -BigFloat;:from(1u32)
            };
            let k_big = BigFloat;:from(k as u32);
            let a_big = BigFloat;:from(a as u32);
            let a_minus_k = a_big.sub(
                &k_big,
                WORKING_PRECISION,
                RoundingMode::ToEven,
            );
            let k_big_minus_half = k_big.sub(
                &BigFloat::from(0.5),
                WORKING_PRECISION,
                RoundingMode::ToEven,
            );
            let a_minus_k_ln = a_minus_k.ln(
                WORKING_PRECISION,
                RoundingMode::ToEven,
                &mut loacl_consts,
            );
            let pow_term = pow_exponent.exp(
                WORKING_PRECISION,
                RoundingMode::ToEven,
                &mut local_consts,
            );
            let exp_term = a_minus_l.exp(
                WORKING_PRECISION,
                RoundingMode::ToEven,
                &mut local_consts,
            );
            let mut fact_big = BigFloat::from(1u32),
            for i in 2..=(k - 1) {
                fact_big = fat_big.mul(
                    &BigFloat::from(i as u32),
                    WORKING_PRECISION,
                    RoundingMode::ToEven,
                );
            }
            let product = pow_term.mul(
                &product,
                WORKING_PRECISION,
                RoundingMode::ToEven,
            );
            let numerator = sign.mul(
                &product,
                WORKING_PRECISION,
                RoundingMode::ToEven,
            );
            let coefficient = numerator.div(
                &fact_big,
                WORKING_PRECISION,
                RoundingMode::ToEven,
            );
            coefficient
        }
}).collect();
coefficients
}

pub fn precomputed_spouge_coefficients_float256 (a: usize) -> Vec<Float256> {
    let coefficients = precomputed_spouge_coefficients_bigfloat(a);
    coefficients
        .iter()
    .map(|c| {
        c.to_string()
            .parse::<Float256>()
        .expect("failed to convert BigFloat to f256")
    })
    .collect()
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
        let coefficients: Result<Vec<Float256>, MathError> = (0..a)
            .into_par_iter()
            .map(|k| {
                let k_f256 = Float256::from(k);
                if k == 0 {
                    Ok(sqrt_two_pi)
                } else {
                    let sign = if k % 2 == 0 {
                        Float256::from(-1.0)
                    } else {
                        Float256::from(1.0)
                    };
                    let fact = factorial(k - 1);
                    let fact_f256 =
                        Float256::from_str(&fact.to_string()).map_err(|_| MathError::Overflow)?;
                    if !fact_f256.is_finite() {
                        return Err(MathError::Overflow);
                    }
                    let a_minus_k = a_f256 - k_f256;
                    let k_minus_half = k_f256 - Float256::from(0.5);
                    let pow_term = a_minus_k.powf(&k_minus_half);
                    let exp_term = (a_minus_k).exp();
                    Ok(sign * pow_term * exp_term / fact_f256)
                }
            })
            .collect();
        coefficients
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
    }
}
