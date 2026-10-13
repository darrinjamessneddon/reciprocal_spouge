pub mod digamma {
    use crate::core::{spouge_coefficients, Complex256, ComplexOps, MathError};
    use ::f256::consts::PI;
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
        if z.im == Float256::from(0.0)
            && z.re <= Float256::from(0.0)
            && z.re.fract() == Float256::from(0.0)
        {
            return Err(MathError::Pole);
        }
        if z.re < Float256::from(0.0) {
            let pi = Complex256::new(PI, Float256::from(0.0));
            let one_minus_z = Complex256::new(Float256::from(1.0) - z.re, -z.im);
            let reflected = digamma_c256(one_minus_z, a)?;
            return Ok(reflected.sub(pi.mul(z).cot().mul(pi)));
        }

        let z_plus_a = z.add(Complex256::new(
            Float256::from(a as f64),
            Float256::from(0.0),
        ));
        let z_plus_half = z.add(Complex256::new(Float256::from(0.5), Float256::from(0.0)));
        let log_z_plus_a = z_plus_a.ln();
        let z_plus_half_div_z_plus_a = z_plus_half.div(z_plus_a);
        let mut psi_z_plus_one = log_z_plus_a
            .add(z_plus_half_div_z_plus_a)
            .sub(Complex256::new(Float256::from(1.0), Float256::from(0.0)));
        let coefficients = spouge_coefficients(a as u64)?;
        let c_0 = coefficients[0];
        let mut numerator = Complex256::new(Float256::from(0.0), Float256::from(0.0));
        let mut denominator = Complex256::new(c_0, Float256::from(0.0));
        for (k, &c_k) in coefficients.iter().enumerate().skip(1) {
            let k_complex = Complex256::new(Float256::from(k as f64), Float256::from(0.0));
            let z_plus_k = z.add(k_complex);
            let z_plus_k_squared = z_plus_k.mul(z_plus_k);
            let c_k = Complex256::new(c_k, Float256::from(0.0));
            numerator = numerator.add(c_k.div(z_plus_k_squared));
            denominator = denominator.add(c_k.div(z_plus_k));
        }
        psi_z_plus_one = psi_z_plus_one.sub(numerator.div(denominator));
        Ok(psi_z_plus_one.sub(Complex256::new(Float256::from(1.0), Float256::from(0.0)).div(z)))
    }
}
