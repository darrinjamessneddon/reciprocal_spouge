pub mod gamma {
    use crate::core::{spouge_coefficients, Complex256, ComplexOps, MathError};
    use ::f256::consts::PI;
    use f256::f256 as Float256;
    use num_complex::Complex;

    pub type C256 = Complex<f256::f256>;

    pub fn generic_gamma<T>(z: Complex<T>) -> Complex<T>
    where
        T: num_traits::Float + num_traits::NumCast,
    {
        let cast = |value| T::from(value).unwrap_or_else(T::nan);
        if z.re < cast(0.5) {
            let pi = cast(std::f64::consts::PI);
            let one = Complex::new(T::one(), T::zero());
            let pi_complex = Complex::new(pi, T::zero());
            return pi_complex / ((pi_complex * z).sin() * generic_gamma(one - z));
        }

        let coefficients = [
            676.5203681218851,
            -1259.1392167224028,
            771.3234287776531,
            -176.6150291621406,
            12.507343278686905,
            -0.13857109526572012,
            9.984369578019572e-6,
            1.5056327351493116e-7,
        ];
        let z_minus_one = z - Complex::new(T::one(), T::zero());
        let mut sum = Complex::new(cast(0.999_999_999_999_809_9), T::zero());
        for (index, coefficient) in coefficients.iter().enumerate() {
            let denominator = z_minus_one + Complex::new(cast((index + 1) as f64), T::zero());
            sum = sum + Complex::new(cast(*coefficient), T::zero()) / denominator;
        }

        let t = z_minus_one + Complex::new(cast(7.5), T::zero());
        let sqrt_two_pi = cast((2.0 * std::f64::consts::PI).sqrt());
        Complex::new(sqrt_two_pi, T::zero())
            * t.powc(z_minus_one + Complex::new(cast(0.5), T::zero()))
            * (-t).exp()
            * sum
    }

    // Create a function to return the value of the gamma function as a string, using Spouge's
    // approximation for complex numbers.
    pub fn spouge(z: Complex256, a: usize) -> Result<String, MathError> {
        let result = spouge_c256(z, a)?;
        Ok(result.to_string())
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
        // Handle poles at zero and the negative integers.
        if z.im == Float256::from(0.0)
            && z.re <= Float256::from(0.0)
            && z.re.fract() == Float256::from(0.0)
        {
            return Err(MathError::Pole);
        }
        if z.re < Float256::from(0.0) {
            // Use the reflection formula for the gamma function
            // gamma(z) = pi / (sin(pi * z) * gamma(1 - z))
            let one_minus_z = Complex256::new(Float256::from(1.0) - z.re, -z.im);
            let pi = PI;
            let pi_complex = Complex256::new(pi, Float256::from(0.0));
            let sin_pi_z = (pi_complex.mul(z)).sin();
            let gamma_one_minus_z = spouge_c256(one_minus_z, a)?;
            return Ok(pi_complex.div(sin_pi_z.mul(gamma_one_minus_z)));
        }
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

    #[cfg(test)]
    mod tests {
        use super::generic_gamma;
        use num_complex::Complex;

        #[test]
        fn generic_gamma_matches_real_reference_values() {
            let factorial_value = generic_gamma(Complex::new(5.0_f64, 0.0));
            assert!((factorial_value.re - 24.0).abs() < 1e-10);
            assert!(factorial_value.im.abs() < 1e-10);

            let half_value = generic_gamma(Complex::new(0.5_f64, 0.0));
            assert!((half_value.re - std::f64::consts::PI.sqrt()).abs() < 1e-10);
            assert!(half_value.im.abs() < 1e-10);
        }
    }
}
