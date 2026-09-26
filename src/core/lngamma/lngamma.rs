pub mod lngamma {
    use crate::core::{Complex256, ComplexOps, MathError, spouge_coefficients};
    use f256::consts::PI;
    use f256::f256 as Float256;

    /// Compute the principal-branch natural logarithm of the gamma function.
    ///
    /// This uses the current Spouge approximation directly in logarithmic form and
    /// applies the reflection formula on the principal branch for inputs with a
    /// negative real part. Poles at zero and the negative integers return
    /// `MathError::Pole`.
    pub fn ln_gamma(z: Complex256, a: usize) -> Result<Complex256, MathError> {
        if a < 2 {
            return Err(MathError::ParameterOutOfRange);
        }

        let max_limit = Float256::from(10000.0);
        if z.re > max_limit || z.im > max_limit {
            return Err(MathError::Overflow);
        }
        if z.re < -max_limit || z.im < -max_limit {
            return Err(MathError::Underflow);
        }

        if z.im == Float256::from(0.0)
            && z.re <= Float256::from(0.0)
            && z.re.fract() == Float256::from(0.0)
        {
            return Err(MathError::Pole);
        }

        if z.re < Float256::from(0.0) {
            let pi_complex = Complex256::new(PI, Float256::from(0.0));
            let log_pi = Complex256::new(PI.ln(), Float256::from(0.0));
            let sin_pi_z = pi_complex.mul(z).sin();
            let one_minus_z = Complex256::new(Float256::from(1.0) - z.re, -z.im);
            return Ok(log_pi.sub(sin_pi_z.ln()).sub(ln_gamma(one_minus_z, a)?));
        }

        let coefficients = spouge_coefficients(a as u64)?;
        let mut sum = Complex256::new(coefficients[0], Float256::from(0.0));
        for (k, &coefficient) in coefficients.iter().enumerate().skip(1) {
            let k_complex = Complex256::new(Float256::from(k as f64), Float256::from(0.0));
            let z_plus_k = z.add(k_complex);
            let c_k = Complex256::new(coefficient, Float256::from(0.0));
            sum = sum.add(c_k.div(z_plus_k));
        }

        let z_plus_a = z.add(Complex256::new(
            Float256::from(a as f64),
            Float256::from(0.0),
        ));
        let z_plus_half = z.add(Complex256::new(Float256::from(0.5), Float256::from(0.0)));
        Ok(z_plus_half
            .mul(z_plus_a.ln())
            .sub(z_plus_a)
            .add(sum.ln())
            .sub(z.ln()))
    }
}
