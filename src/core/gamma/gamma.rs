pub mod gamma {
    use crate::MathError;
    use crate::core::complex::complex::complex::{Complex256, ComplexOps};
    use crate::core::shared::shared::shared::spouge_coefficients;
    use ::f256::consts::PI;
    use f256::f256 as Float256;

    // Create a function to return the value of the gamma function as a string, using Spouge's
    // approximation for complex numbers.
    pub fn spouge(z: Complex256, a: usize) -> Result<String, MathError> {
        let result = spouge_c256(z, a)?;
        Ok(result.to_string)
    }

    // Create a function to perform implementation of Spouge's approximation for the gamma function
    // for complex numbers, resulting in a Complex256 value.

    pub fn spouge_c256(z: Complex256, a: usize) -> Result<Complex256, MathError> {
        if a < 2 {
            return Err(MathError::ParameterOutOfRange);
        }
        // The spouge approximation calculates takes a z value and a parameter a, and returns the value of
        // gamma(z + 1). To compensate for this, we will subtract 1 from the input value z before passing it to the approximation.
        let max_limit = Float256::from(10000.0);
        if z.re > max_limit || z.im > max_limit {
            return Err(MathError::Overflow);
        }
        if z.re < -max_limit || z.im < -max_limit {
            return Err(MathError::Underflow);
        }
        // Handle the singularity at z=0
        if z.re == Float256::from(0.0) && z.im == Float256::from(0.0) {
            return Ok(Complex256::new(Float256::INFINITY, Float256::from(0.0)));
        }
        // Handle the case for negative integers (negative real part and zero imaginary part)
        // Negative integers are poles of the gamma function, so we return infinity for those cases
        if z.re < Float256::from(0.0)
            && z.im == Float256::from(0.0)
            && z.re.fract() == Float256::from(0.0)
        {
            return Ok(Complex256::new(Float256::INFINITY, Float256::from(0.0)));
        }
        // Handle the case for negative non-integer values (negative real part and non-zero imaginary part)
        if z.re < Float256::from(0.0) && z.im != Float256::from(0.0) {
            // Use the reflection formula for the gamma function
            // gamma(z) = pi / (sin(pi * z) * gamma(1 - z))
            let one_minus_z = Complex256::new(Float256::from(1.0) - z.re, -z.im);
            let pi = PI;
            let pi_complex = Complex256::new(pi, Float256::from(0.0));
            let sin_pi_z = (pi_complex.mul(z)).sin();
            let gamma_one_minus_z = spouge_c256(one_minus_z, a)?;
            return Ok(pi_complex.div(sin_pi_z.mul(gamma_one_minus_z)));
        }
        // Handle the case for z values with negative non-integer real part and zero imaginary part
        if z.re < Float256::from(0.0) && z.im == Float256::from(0.0) {
            // Use the reflection formula for the gamma function
            // gamma(z) = pi / (sin(pi * z) * gamma(1 - z)
            let one_minus_z = Complex256::new(Float256::from(1.0) - z.re, Float256::from(0.0));
            let pi = PI;
            let pi_complex = Complex256::new(pi, Float256::from(0.0));
            let sin_pi_z = (pi_complex.mul(z)).sin();
            let gamma_one_minus_z = spouge_c256(one_minus_z, a)?;
            return Ok(pi_complex.div(sin_pi_z.mul(gamma_one_minus_z)));
        }
        // Handle the case for z values with real z between 0 and 1 (0 , Re(z) < 1)
        if z.re < Float256::from(0.0) && z.re.fract() != Float256::from(0.0) {
            // use the reflection formula for the gamma function
            // gamma(z) = pi / (sin(pi * z) * gamma(1 - z))
            let one_minus_z = Complex256::new(Float256::from(1.0) - z.re, -z.im);
            let pi = PI;
            let pi_complex = Complex256::new(pi, Float256::from(0.0));
            let sin_pi_z = (pi_complex.mul(z)).sin();
            let gamma_one_minus_z = spouge_c256(one_minus_z, a)?;
            Ok(pi_complex.div(sin_pi_z.mul(gamma_one_minus_z)))
        } else {
            // Handle the case for other z values by using Spouge's approximation directly
            let mut sum = Complex256::new(Float256::from(0.0), Float256::from(0.0));
            let z_plus_a = z.add(Complex256::new(
                Float256::from(a as f64),
                Float256::from(0.0),
            ));
            let z_plus_half = z.add(Complex256::new(Float256::from(0.5), Float256::from(0.0)));
            let pow_term = z_plus_a.powc(z_plus_half);
            let exp_term = z_plus_a
                .mul(Complex256::new(Float256::from(-1.0), Float256::from(0.0)))
                .exp();
            let coefficients = spouge_coefficients(a as u64)?;
            let c_0 = coefficients[0];
            // Compute the sum of c_l / (z + k) for k = 1 to a-1
            for (k, &c_l) in coefficients.iter().enumerate().skip(1) {
                let k_complex = Complex256::new(Float256::from(k as f64), Float256::from(0.0));
                let z_plus_k = z.add(k_complex);
                let c_k = Complex256::new(c_l, Float256::from(0.0));
                let term = c_k.div(z_plus_k);
                sum = sum.add(term);
            }
            let c_0_complex = Complex256::new(c_0, Float256::from(0.0));
            let result = pow_term.mul(exp_term).mul(c_0_complex.add(sum));
            Ok(result.div(z))
        }
    }
}
