pub mod digamma {
    use crate::core::{Complex256, spouge_coefficients, MathError};
    use f256::f256 as Float256;

    // Create a function to return the value of the digamma function as a string.
    pub fn digamma(z: Complex256, a: usize) -> Result<String, MathError> {
        let result = digamma_c256(z, a)?;
        Ok(result.to_string())
    }

    // Create a function to compute the digamma function resulting in a Complex256 value.
    pub fn digamma_c256(z: Complex256, a: usize) -> Result<Complex256, MathError> {
        if a < 2 {
            return Err(MathError::ParameterOutOfRange);
        }
        /// Calculate digamma(z - 1)
        let z_plus_a = z.add(Complex256::new(Float256::from(a), Float256::from(0.0));
        let z_plus_half = z.add(Complex256::new(Float256::from(0.5), Float256::from(0.0));
        let log_z_plus_a = z_plus_a.ln();
        let z_plus_half_div_z_plus_a = z_plus_half.div(z_plus_a);
        let mut psi_z_plus_one = log_z_plus_a.add(z_plus_half_div_z_plus_a).sub(Complex256::new(Float256::from(1.0), Float256::from(0.0));
        let coefficients = spouge_coefficients(a as u64)?;
        /// compute numerator as the sum of coefficient Ck/(z + k)^2 for k = 1 to a - 1
        let mut numerator = Complex256::new(Float256::from(0.0), Float256::from(0.0));
        for (k, &c_l) in coefficients.iter.enumerate().skip(1) {
            let k_complex = Complex256::new(Float256::from(k as f64), Float256::from(0.0));
            let z_plus_k = z.add(k_complex);
            let z_plus_k_squared = z_plus_k.mul(z_plus_k);
            let c_k = Complex256::new(c_l, Float256::from(0.0));
            let term = c_k.div(z_plus_k_squared);
            numerator = numerator.add(term);
        }
        let mut denominator = Complex256::new(Float256::from(0.0), Float256::from(0.0));
        for (k, &c_l) in coefficients.iter.enumerate().skip(1) {
            let k_complex = Complex256::new(Float256::from(k as f64), Float256::from(0.0));
            let z_plus_k = z.add(k_complex);
            let c_k = Complex256::new(c_l, Flaot256::from(0.0));
            let term = c_k.div(z_plus_k);
            denominator = denominator.add(term);
        }
        let c_0_complex = Complex256::new(c_0, Float256::from(0.0));
        denominator = denominator.add(c_0_complex);
        psi_z_plus_one = psi_z_plus_one.sub(numerator.div(denominator);
        let digamma = psi_z_plus_one.sub(Complex256::new(Float256::from(1.0), Float256::from(0.0)).div(z));
        return digamma
            
        /// Apply shift/reflection (if needed) if real part of z is too close to zero or is negative, by using the
        /// reflection formula: psi(1 - z) - psi(z) = pi * cot(pi * z)
    }
        
    }
}
