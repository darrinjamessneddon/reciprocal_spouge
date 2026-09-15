pub mod shared {
    use core::str::FromStr;

    use ::f256::consts::TAU;
    use f256::f256 as Float256;
    use num_bigint::BigUint;
    use rayon::prelude::*;

    use crate::core::MathError;

    pub fn factorial(n: u64) -> BigUint {
        (1..=n).into_par_iter().map(BigUint::from).product()
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
}

/// Add structs and functions to facillitate the use of Kahan summation by feeding f256 terms into a Kahan accumulator using Neumaier's algorithm
#[derive(Debug, Clone, Copy, Default)]
pub struct KahanF256 {
    pub kahan_sum: Float256,
    pub compensation: Float256,
    // add implementation code here
}
impl KahanF256 {
    /// Creates a new compensated accumulator initialized to zero.
pub fn new(kahan_sum: Float256, compensation: Float256) -> Self {
    KahanF256{kahan_sum, compensation}
    }
}
/// Code to feed a new f256 term into the accumulator using Neumaier's algorithm.
pub fn add(&mut self, term: Float256) -> Self {
    let mut kahan_sum = Float256::from(0.0);
    let mut compensation = Float256::from(0.0);
    let mut terms = vec::new();
    for i in 0..terms.len() {
        let  mut t = kahan_sum + terms[i];
        if kahan_sum.abs() >= terms[i] {
            compensation += (kahan_sum - t) + terms[i]
        } else {
            compensation += (terms[i] - t) + kahan_sum
        }
        kahan_sum = t
    }
    return kahan_sum + compensation
}
