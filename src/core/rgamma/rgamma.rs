pub mod rgamma {
    use crate::core::complex::complex::complex::{Complex256, ComplexOps};
    use crate::core::shared::shared::shared::spouge_coefficients;
    use f256::f256 as Float256;
    use f256::consts::PI;

    // Create a function to approximate the reciprocal of the gamma function using a re-arranged
    // version of the Spouge approximation for complex numbers, resulting a Complex256 value
    // that is then converted to a string and returned.
    pub fn rspouge(z: Complex256, a: usize) -> String {
        let result = rspouge_c256(z, a as i32);
        let result_re_str = result.re.to_string();
        let result_im_str = result.im.to_string();
        return format!("{} + {}", result_re_str, result_im_str);
    }


    // Create a function to approximate the gamma function using a re-arranged version of the Spouge approximation
    pub fn rspouge_c256(z: Complex256, a: i32) -> Complex256 {
        if a < 2 {
            panic!("Parameter 'a' must be greater than or equal to 2");
        }
        // The re-arranged approximation to calculate reciprocal gamma: 1/gamma(z + 1) = (z + a)^(-z-0.5) * exp(z + a) * (1/(c_0 + sum_(k=1 to a - 1) c_k/(z + k)))
        let max_limit = Float256::from(10000.0); // Set an arbitrary limit for the input value to avoid overflow or underflow issues
        if z.re > max_limit || z.im > max_limit {
            panic!("Input value is too large for Spouge's approximation");
        }
        if z.re < -max_limit || z.im < -max_limit {
            panic!("Input value is too small for Spouge's approximation");
        }
        // Handle the zero value case for z=0, which is a singularity for the gamma function
        if z.re == Float256::from(0.0) && z.im == Float256::from(0.0) {
            return Complex256::new(Float256::from(0.0), Float256::from(0.0));
        }
        // Handle the case for negative integers where the imaginary part is zero.
        if z.re < Float256::from(0.0) && z.im == Float256::from(0.0) {
            if z.re.fract() == Float256::from(0.0) {
                return Complex256::new(Float256::from(0.0), Float256::from(0.0));
            }
        }
        // Handle the case for negative non-integer values where the imaginary part is non-zero
        // use the reflection formula: gamma(z) = pi / (sin(pi * z) * gamma(1 - z))
        if z.re < Float256::from(0.0) && z.im != Float256::from(0.0) {
            let pi = PI;
            let pi_complex = Complex256::new(pi, Float256::from(0.0));
            let sin_pi_z = (pi_complex.mul(z)).sin();
            let rgamma_one_minus_z = rspouge_c256(Complex256::new(Float256::from(1.0) - z.re, -z.im), a);
            let numerator = sin_pi_z;
            let denominator = pi_complex.mul(rgamma_one_minus_z);
            return numerator.div(denominator);

        } else {
            // For positive values of z, we can use the re-arranged Spouge approximation directly
            let mut sum = Complex256::new(Float256::from(0.0), Float256::from(0.0));
            let coeffs = spouge_coefficients((a as u64).try_into().unwrap()).unwrap();
            let c_0 = Complex256::new(coeffs[0].into(), Float256::from(0.0));
            for k in 1..a as u64 {
                let c_k = Complex256::new(coeffs[k as usize].into(), Float256::from(0.0));
                let z_plus_k = Complex256::new(z.re + Float256::from(k as f64), z.im);
                let term = c_k.div(z_plus_k);
                sum = sum.add(term);
            }
            sum = sum.add(c_0);
            let z_plus_a = Complex256::new(z.re + Float256::from(a as f64), z.im);
            let z_plus_half = Complex256::new(z.re + Float256::from(0.5), z.im);
            let neg_z_plus_half = Complex256::new(-z_plus_half.re, -z_plus_half.im);
            let pow_term = z_plus_a.powc(neg_z_plus_half);
            let exp_term = (z_plus_a).exp();
            let numerator = pow_term.mul(exp_term);
            let denominator = sum;
            let mut rgamma = numerator.div(denominator);
            rgamma = rgamma.mul(z);
            return rgamma
        }
    }
}