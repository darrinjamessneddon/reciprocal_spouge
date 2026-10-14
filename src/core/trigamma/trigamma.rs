pub mod trigamma {
    use crate::core::{spouge_coefficients, Complex256, ComplexOps, MathError};
    use f256::consts::PI;
    use f256::f256 as Float256;

    pub fn trigamma(z: Complex256, a: usize) -> Result<String, MathError> {
        let result = trigamma_c256(z, a)?;
        Ok(result.to_string())
    }

    pub fn trigamma_c256(z: Complex256, a: usize) -> Result<Complex256, MathError> {
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
            let sin_pi_z = pi.mul(z).sin();
            let reflected = trigamma_c256(one_minus_z, a)?;
            return Ok(pi.mul(pi).div(sin_pi_z.mul(sin_pi_z)).sub(reflected));
        }

        let one = Complex256::new(Float256::from(1.0), Float256::from(0.0));
        let z_plus_a = z.add(Complex256::new(
            Float256::from(a as f64),
            Float256::from(0.0),
        ));
        let z_plus_a_squared = z_plus_a.mul(z_plus_a);
        let base = one.div(z_plus_a).add(
            Complex256::new(
                Float256::from(a as f64) - Float256::from(0.5),
                Float256::from(0.0),
            )
            .div(z_plus_a_squared),
        );
        let coefficients = spouge_coefficients(a as u64)?;
        let mut numerator = Complex256::new(Float256::from(0.0), Float256::from(0.0));
        let mut numerator_derivative = Complex256::new(Float256::from(0.0), Float256::from(0.0));
        let mut denominator = Complex256::new(coefficients[0], Float256::from(0.0));

        for (k, &coefficient) in coefficients.iter().enumerate().skip(1) {
            let k_complex = Complex256::new(Float256::from(k as f64), Float256::from(0.0));
            let z_plus_k = z.add(k_complex);
            let z_plus_k_squared = z_plus_k.mul(z_plus_k);
            let z_plus_k_cubed = z_plus_k_squared.mul(z_plus_k);
            let coefficient = Complex256::new(coefficient, Float256::from(0.0));
            numerator = numerator.add(coefficient.div(z_plus_k_squared));
            numerator_derivative = numerator_derivative.sub(
                Complex256::new(Float256::from(2.0), Float256::from(0.0))
                    .mul(coefficient)
                    .div(z_plus_k_cubed),
            );
            denominator = denominator.add(coefficient.div(z_plus_k));
        }

        let correction_derivative = numerator_derivative
            .div(denominator)
            .add(numerator.mul(numerator).div(denominator.mul(denominator)));
        Ok(base.sub(correction_derivative).add(one.div(z.mul(z))))
    }

    #[cfg(test)]
    mod tests {
        use super::trigamma_c256;
        use crate::core::{Complex256, MathError};

        #[test]
        fn trigamma_matches_known_values_and_reports_poles() {
            let at_one = trigamma_c256(Complex256::from_f64(1.0, 0.0), 25)
                .unwrap()
                .to_f64()
                .0;
            assert!((at_one - std::f64::consts::PI.powi(2) / 6.0).abs() < 1e-8);

            assert_eq!(
                trigamma_c256(Complex256::from_f64(0.0, 0.0), 25),
                Err(MathError::Pole)
            );
            assert_eq!(
                trigamma_c256(Complex256::from_f64(1.0, 0.0), 1),
                Err(MathError::ParameterOutOfRange)
            );
        }
    }
}
