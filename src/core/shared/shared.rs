pub mod shared {
    use rayon::prelude::*;
    use num_bigint::BigUint;
    use num_traits::ToPrimitive;
    use f256::f256 as Float256;
    use ::f256::consts::TAU;

    
    pub fn factorial(n: u64) -> BigUint {
        (1..=n).into_par_iter().map(BigUint::from).product()
    }

/// Compute the Spouge coefficients for a given parameter 'a'. The coefficients are used in Spouge's approximation of the gamma function.
/// allow for error handling by returning a Result type, which can either be Ok with the coefficients
/// or Err with a String describing the error.
    pub fn spouge_coefficients(a: u64) -> Result<Vec<Float256>, String> {
        if a == 0 {
            return Err("Spouge coefficients require a > 0".to_string());
        }
        if a < 2 || a > 21 {
            return Err(format!("Parameter 'a' is out of range. Must be between 2 and 21, got {}", a));
        }

        let sqrt_two_pi = TAU.sqrt();
        let a_f256 = Float256::from(a as u64);
        let coefficients: Vec<Float256> = (0..a).into_par_iter().map(|k| {
            let k_f256 = Float256::from(k as u64);
            if k == 0 {
                sqrt_two_pi
            } else {
                let sign = if k % 2 == 0 { Float256::from(-1.0) } else { Float256::from(1.0) };
                let fact = factorial(k - 1);
                let fact_f256 = Float256::from(fact.to_u64().unwrap_or(0));
                let a_minus_k = a_f256 - k_f256;
                let k_minus_half = k_f256 - Float256::from(0.5);
                let pow_term = a_minus_k.powf(&k_minus_half);
                let exp_term = (a_minus_k).exp();
                sign * pow_term * exp_term / fact_f256
            }
        }).collect();
        Ok(coefficients)
    }
}