pub mod lngamma {
    use crate::core::gamma::spouge_c256;
    use crate::core::{Complex256, ComplexOps, MathError};
    use f256::f256 as Float256;

    /// Compute the principal-branch natural logarithm of the gamma function.
    ///
    /// This delegates to the current Spouge gamma implementation and then applies
    /// the principal complex logarithm to the finite result. Poles at zero and the
    /// negative integers return `MathError::Pole`.
    pub fn ln_gamma(z: Complex256, a: usize) -> Result<Complex256, MathError> {
        if a < 2 {
            return Err(MathError::ParameterOutOfRange);
        }

        if z.im == Float256::from(0.0)
            && z.re <= Float256::from(0.0)
            && z.re.fract() == Float256::from(0.0)
        {
            return Err(MathError::Pole);
        }

        let gamma = spouge_c256(z, a)?;
        Ok(gamma.ln())
    }
}
