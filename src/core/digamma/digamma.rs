pub mod digamma {
    use crate::core::{Complex256, MathError};
    use f256::f256 as Float256;

    // Create a function to return the value of the digamma function as a string.
    pub fn digamma(z: Complex256, a: usize) -> Result<String, MathError> {
        let result = digamma_c256(z, a)?;
        Ok(result.to_string())
    }

    // Create a function to compute the digamma function resulting in a Complex256 value.
    pub fn digamma_c256(_z: Complex256, a: usize) -> Result<Complex256, MathError> {
        if a < 2 {
            return Err(MathError::ParameterOutOfRange);
        }
        // for now return placeholder
        Ok(Complex256::new(Float256::from(0.0), Float256::from(0.0)))
    }
}
