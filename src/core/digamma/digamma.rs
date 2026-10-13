pub mod digamma {
    use crate::core::{Complex256, coefficients, MathError};
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
        let coeffs = coefficients(a);
        /// compute numerator as the sum of coefficient Ck/(z + k)^2 for k = 1 to a - 1 using parallel computation
        
        /// compute denominator as C0 + the sum of coeffieicnt Ck/(z + k) for k = 1 to a - 1 using parallel computation

        /// subtract numerator/denominator from pis_z_plus_one
        
        /// Shift back to psi(z) using the relation psi(z) = psi(z + 1) - (1/z)
        
        /// Apply shift/reflection (if needed) if real part of z is too close to zero or is negative, by using the
        /// reflection formula: psi(1 - z) - psi(z) = pi * cot(pi * z)
    }
        
    }
}
