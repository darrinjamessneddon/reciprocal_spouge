pub mod lngamma {
    use crate::core::{spouge_coefficients, Complex256, ComplexOps, MathError};
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

    /// Spouge series sum `c_0 + sum_k c_k / (z + k)`.
    fn spouge_sum(z: Complex256, coefficients: &[Float256]) -> Complex256 {
        let mut sum = Complex256::new(coefficients[0], Float256::from(0.0));
        for (k, &coefficient) in coefficients.iter().enumerate().skip(1) {
            let k_complex = Complex256::new(Float256::from(k as f64), Float256::from(0.0));
            let c_k = Complex256::new(coefficient, Float256::from(0.0));
            sum = sum.add(c_k.div(z.add(k_complex)));
        }
        sum
    }

    /// Continuous log-gamma for `Re z >= 0`. The only term of the Spouge form
    /// whose principal branch can jump is `ln(sum)`, so its argument is unwrapped
    /// along the vertical path from `Re z` to `z`.
    fn loggamma_right(z: Complex256, a: usize) -> Result<Complex256, MathError> {
        let coefficients = spouge_coefficients(a as u64)?;
        let zero = Float256::from(0.0);
        let pi = PI;
        let two_pi = PI + PI;

        let steps = (z.to_f64().1.abs().ceil() as usize)
            .saturating_mul(2)
            .max(1);
        let mut prev_arg = spouge_sum(Complex256::new(z.re, zero), &coefficients)
            .ln()
            .im;
        let mut unwrapped = prev_arg;
        let mut ln_sum = Complex256::new(zero, zero);
        for step in 1..=steps {
            let t = Float256::from(step as f64) / Float256::from(steps as f64);
            let point = if step == steps {
                z
            } else {
                Complex256::new(z.re, z.im * t)
            };
            ln_sum = spouge_sum(point, &coefficients).ln();
            let mut delta = ln_sum.im - prev_arg;
            if delta > pi {
                delta -= two_pi;
            } else if delta < -pi {
                delta += two_pi;
            }
            unwrapped += delta;
            prev_arg = ln_sum.im;
        }
        ln_sum.im = unwrapped;

        let z_plus_a = z.add(Complex256::new(Float256::from(a as f64), zero));
        let z_plus_half = z.add(Complex256::new(Float256::from(0.5), zero));
        Ok(z_plus_half
            .mul(z_plus_a.ln())
            .sub(z_plus_a)
            .add(ln_sum)
            .sub(z.ln()))
    }

    /// The function lngamma(z) computes the principal Log-Gamma function. Therefore it is necessary to add the loggamma function:
    /// The relationship between the two functions is given by:
    /// loggamma(z) = ln(gamma(z)) + 2 * pi * i * k(z), where k(z) is an integer that corrects for the winding number
    ///
    /// The result is the analytic (continuous) branch of log-gamma on the plane
    /// cut along the negative real axis, so `exp(loggamma(z)) == gamma(z)` and
    /// the imaginary part includes the winding-number term. Poles at zero and
    /// the negative integers return `MathError::Pole`.
    pub fn loggamma(z: Complex256, a: usize) -> Result<Complex256, MathError> {
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

        let zero = Float256::from(0.0);
        if z.im == zero && z.re <= zero && z.re.fract() == zero {
            return Err(MathError::Pole);
        }

        if z.re >= zero {
            return loggamma_right(z, a);
        }

        // Recurrence: lnGamma(z) = lnGamma(z + n) - sum_{k<n} Log(z + k), where
        // each principal Log is analytic off the negative real axis.
        let n = (-z.to_f64().0).floor() as usize + 1;
        let shifted = z.add(Complex256::new(Float256::from(n as f64), zero));
        let mut result = loggamma_right(shifted, a)?;
        for k in 0..n {
            let term = z.add(Complex256::new(Float256::from(k as f64), zero));
            result = result.sub(term.ln());
        }
        Ok(result)
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        const A: usize = 12;

        // Reference values of ln(Gamma(x)) (high precision, from known closed forms).
        const REFERENCES: [(f64, f64); 9] = [
            (0.5, 0.572_364_942_924_700_1),
            (1.0, 0.0),
            (1.5, -0.120_782_237_635_245_22),
            (2.0, 0.0),
            (3.0, std::f64::consts::LN_2),
            (4.0, 1.791_759_469_228_055),
            (5.0, 3.178_053_830_347_945_8),
            (10.0, 12.801_827_480_081_469),
            (20.0, 39.339_884_187_199_495),
        ];

        fn assert_close(input: f64, expected: f64, actual: f64) {
            let error = (actual - expected).abs();
            let tolerance = 1e-9 + 1e-9 * expected.abs();
            assert!(
                error <= tolerance,
                "loggamma({input}): expected {expected}, actual {actual}, error {error}"
            );
        }

        fn real_loggamma(x: f64) -> Complex256 {
            loggamma(Complex256::from_f64(x, 0.0), A).unwrap()
        }

        #[test]
        fn loggamma_matches_reference_values() {
            for (x, expected) in REFERENCES {
                let (re, im) = real_loggamma(x).to_f64();
                assert_close(x, expected, re);
                assert!(im.abs() <= 1e-9, "loggamma({x}): imaginary part {im}");
            }
        }

        #[test]
        fn loggamma_small_positive_input() {
            // ln(Gamma(x)) ~ -ln(x) - gamma_E * x for small x
            let x = 1e-3;
            let expected = 6.907_178_885_383_853; // ln(Gamma(0.001))
            assert_close(x, expected, real_loggamma(x).to_f64().0);
        }

        #[test]
        fn loggamma_satisfies_recurrence() {
            for x in [0.25, 0.75, 2.5, 7.3] {
                let lhs = real_loggamma(x + 1.0).to_f64().0;
                let rhs = real_loggamma(x).to_f64().0 + x.ln();
                assert_close(x, rhs, lhs);
            }
        }

        #[test]
        fn loggamma_rejects_poles_and_invalid_parameter() {
            assert_eq!(
                loggamma(Complex256::from_f64(0.0, 0.0), A),
                Err(MathError::Pole)
            );
            assert_eq!(
                loggamma(Complex256::from_f64(-2.0, 0.0), A),
                Err(MathError::Pole)
            );
            assert_eq!(
                loggamma(Complex256::from_f64(1.0, 0.0), 1),
                Err(MathError::ParameterOutOfRange)
            );
        }
    }
}
