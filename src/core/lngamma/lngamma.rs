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

    #[cfg(test)]
    mod stability_tests {
        use super::*;
        use proptest::prelude::*;
        use std::str::FromStr;

        const A_REF: usize = 40;

        fn f(s: &str) -> Float256 {
            Float256::from_str(s).unwrap()
        }

        fn c(re: f64, im: f64) -> Complex256 {
            Complex256::from_f64(re, im)
        }

        fn to64(x: Float256) -> f64 {
            x.to_string().parse::<f64>().unwrap_or(f64::NAN)
        }

        fn two_pi() -> Float256 {
            PI + PI
        }

        /// Relative error `|actual - expected| / max(|expected|, 1)`, as f64.
        fn rel_err(actual: Float256, expected: Float256) -> f64 {
            let denom = expected.abs().max(Float256::from(1.0));
            to64((actual - expected).abs() / denom)
        }

        /// Number of correct decimal digits implied by a relative error.
        fn digits(err: f64) -> f64 {
            if err <= 0.0 {
                f64::INFINITY
            } else {
                -err.log10()
            }
        }

        /// Distance between two angles modulo 2*pi.
        fn angle_diff(actual: Float256, expected: Float256) -> Float256 {
            let tau = two_pi();
            let mut d = (actual - expected) % tau;
            if d > PI {
                d -= tau;
            } else if d < -PI {
                d += tau;
            }
            d.abs()
        }

        fn assert_angle_close(actual: Float256, expected: Float256, tol: f64) {
            let d = to64(angle_diff(actual, expected));
            assert!(
                d <= tol,
                "angles differ by {d} (mod 2pi): {actual} vs {expected}"
            );
        }

        fn assert_re_im(z: Complex256, re: &str, im: &str, tol_re: f64, tol_im: f64, label: &str) {
            let er = rel_err(z.re, f(re));
            let ei = rel_err(z.im, f(im));
            eprintln!(
                "{label}: re digits {:.1}, im digits {:.1}",
                digits(er),
                digits(ei)
            );
            assert!(er <= tol_re, "{label}: re error {er}, got {}", z.re);
            assert!(ei <= tol_im, "{label}: im error {ei}, got {}", z.im);
        }

        #[test]
        fn reference_values_real_axis() {
            let refs = [
                ("1e-30", "69.0775527898213705205397436405"),
                ("1e-10", "23.0258509298827352736979859311"),
                ("0.5", "0.572364942924700087071713675677"),
                ("1", "0"),
                ("2", "0"),
                ("10", "12.8018274800814696112077178746"),
                ("100", "359.13420536957539877604401046"),
                ("1000", "5905.22042320918121182607691236"),
                ("9999.5", "82095.1123637576392281574846617"),
            ];
            for (x, expected) in refs {
                let z = Complex256::new(f(x), Float256::from(0.0));
                let r = loggamma(z, A_REF).unwrap();
                assert_re_im(r, expected, "0", 1e-25, 1e-25, x);
                let r = ln_gamma(z, A_REF).unwrap();
                assert_re_im(r, expected, "0", 1e-25, 1e-25, x);
            }
        }

        #[test]
        fn reference_values_complex() {
            let refs = [
                (
                    c(0.5, 10.0),
                    "-14.7890247347442934505328871802",
                    "13.0300200349110898508075452634",
                ),
                (
                    c(1.0, 1000.0),
                    "-1566.42351062220087796351437472",
                    "5908.5405938121983892516852624",
                ),
                (
                    c(100.0, 100.0),
                    "315.078044599493313233406035654",
                    "473.32107821888029677925879947",
                ),
            ];
            for (z, re, im) in refs {
                let r = loggamma(z, A_REF).unwrap();
                assert_re_im(r, re, im, 1e-22, 1e-22, &format!("{z}"));
            }
            let z = Complex256::new(f("1e-5"), f("1e-5"));
            let r = loggamma(z, A_REF).unwrap();
            assert_re_im(
                r,
                "11.1663461025336075514131807713",
                "-0.785403935389604719630712957404",
                1e-22,
                1e-22,
                "1e-5+1e-5i",
            );
        }

        #[test]
        fn reference_values_negative_real_axis() {
            // loggamma is continuous from above the cut, so Im = -pi * (number of
            // sign changes) as given by mpmath.
            let refs = [
                (
                    -0.5,
                    "1.26551212348464539648894579713",
                    "-3.14159265358979323846264338328",
                ),
                (
                    -1.5,
                    "0.86004701537648101451093268167",
                    "-6.28318530717958647692528676656",
                ),
                (
                    -10.5,
                    "-15.1472705907178411461011763955",
                    "-34.5575191894877256230890772161",
                ),
            ];
            for (x, re, im) in refs {
                let r = loggamma(c(x, 0.0), A_REF).unwrap();
                assert_re_im(r, re, im, 1e-22, 1e-22, &format!("{x}"));
            }
        }

        #[test]
        fn near_poles_behave_like_minus_ln_eps() {
            let eps = f("1e-20");
            let expected = 20.0 * std::f64::consts::LN_10;
            for n in [1.0, 2.0, 5.0] {
                for sign in [1.0, -1.0] {
                    let x = Float256::from(-n) + eps * Float256::from(sign);
                    let z = Complex256::new(x, Float256::from(0.0));
                    let r = loggamma(z, A_REF).unwrap();
                    let fact: f64 = (1..=(n as u64)).map(|k| k as f64).product();
                    let expected = expected - fact.ln();
                    let re = to64(r.re);
                    assert!(re.is_finite());
                    assert!((re - expected).abs() < 1e-6, "n={n} sign={sign}: {re}");
                }
            }
        }

        #[test]
        fn convergence_over_a() {
            let z = c(3.7, 0.0);
            let exact = loggamma(z, 100).unwrap();
            let mut errors = Vec::new();
            for a in [10_usize, 20, 40, 60] {
                let r = loggamma(z, a).unwrap();
                let e = rel_err(r.re, exact.re);
                eprintln!("a={a}: {:.1} digits", digits(e));
                errors.push(e);
            }
            for w in errors.windows(2) {
                assert!(w[1] <= w[0], "error did not shrink: {errors:?}");
            }
            // Spouge bound: a^(-1/2) (2 pi)^(-(a + 1/2)), times a generous constant.
            for (a, e) in [10_usize, 20, 40].iter().zip(&errors) {
                let a = *a as f64;
                let bound = 1e3 * a.powf(-0.5) * (2.0 * std::f64::consts::PI).powf(-(a + 0.5));
                assert!(*e <= bound, "a={a}: error {e} exceeds bound {bound}");
            }
            // f256 precision floor (~70 digits)
            assert!(errors[3] < 1e-45, "a=60 error {}", errors[3]);
        }

        #[test]
        fn known_values() {
            let a = 30;
            let tol = 1e-25;
            assert!(rel_err(loggamma(c(1.0, 0.0), a).unwrap().re, Float256::from(0.0)) < tol);
            assert!(rel_err(loggamma(c(2.0, 0.0), a).unwrap().re, Float256::from(0.0)) < tol);
            let half = loggamma(c(0.5, 0.0), a).unwrap();
            assert!(rel_err(half.re, PI.ln() / Float256::from(2.0)) < tol);
            let three = loggamma(c(3.0, 0.0), a).unwrap();
            assert!(rel_err(three.re, Float256::from(2.0).ln()) < tol);
            let mut fact = Float256::from(1.0);
            for n in 1..=20u32 {
                let r = loggamma(c(n as f64, 0.0), a).unwrap();
                assert!(rel_err(r.re, fact.ln()) < tol, "n={n}");
                fact *= Float256::from(n as f64);
            }
        }

        #[test]
        fn recurrence_holds() {
            for z in [c(0.3, 0.7), c(2.5, -3.0), c(10.0, 10.0), c(-2.3, 1.5)] {
                let lhs = loggamma(z.add(c(1.0, 0.0)), A_REF)
                    .unwrap()
                    .sub(loggamma(z, A_REF).unwrap())
                    .sub(z.ln());
                assert!(to64(lhs.re.abs()) < 1e-22, "{z}: {}", lhs.re);
                assert_angle_close(lhs.im, Float256::from(0.0), 1e-22);
            }
        }

        #[test]
        fn reflection_holds() {
            let pi_c = Complex256::new(PI, Float256::from(0.0));
            for z in [c(0.3, 0.7), c(0.5, 5.0), c(0.8, -2.0), c(0.25, 0.0)] {
                let lhs = ln_gamma(z, A_REF)
                    .unwrap()
                    .add(ln_gamma(c(1.0, 0.0).sub(z), A_REF).unwrap());
                let rhs = Complex256::new(PI.ln(), Float256::from(0.0)).sub(pi_c.mul(z).sin().ln());
                assert!(to64((lhs.re - rhs.re).abs()) < 1e-22, "{z}");
                assert_angle_close(lhs.im, rhs.im, 1e-22);
            }
        }

        #[test]
        fn conjugate_symmetry() {
            for z in [c(0.3, 0.7), c(5.0, 20.0), c(-2.3, 1.5), c(100.0, 3.0)] {
                let a = loggamma(z.conj(), A_REF).unwrap();
                let b = loggamma(z, A_REF).unwrap().conj();
                assert!(to64((a.re - b.re).abs()) < 1e-25, "{z}");
                assert!(to64((a.im - b.im).abs()) < 1e-25, "{z}");
            }
        }

        #[test]
        fn imaginary_axis_symmetry() {
            // Re lnGamma(x + iy) is even in y for the same x.
            for z in [c(0.5, 10.0), c(2.0, 7.0), c(1e-5, 3.0)] {
                let up = loggamma(z, A_REF).unwrap();
                let down = loggamma(z.conj(), A_REF).unwrap();
                assert!(to64((up.re - down.re).abs()) < 1e-25, "{z}");
            }
        }

        #[test]
        fn small_positive_reals_follow_minus_ln_x() {
            // lnGamma(x) = -ln x - gamma_E x + O(x^2)
            let gamma_e = 0.577_215_664_901_532_9;
            for x in [1e-30_f64, 1e-20, 1e-10, 1e-6] {
                let r = to64(loggamma(c(x, 0.0), A_REF).unwrap().re);
                let expected = -x.ln() - gamma_e * x;
                assert!((r - expected).abs() < 1e-9 * expected.abs(), "x={x}");
            }
        }

        #[test]
        fn large_reals_match_stirling() {
            let ln_two_pi = two_pi().ln();
            for x in [100.0_f64, 1e3, 5e3, 9999.0] {
                let xf = Float256::from(x);
                let x2 = xf * xf;
                let stirling = (xf - Float256::from(0.5)) * xf.ln() - xf
                    + ln_two_pi / Float256::from(2.0)
                    + Float256::from(1.0) / (Float256::from(12.0) * xf)
                    - Float256::from(1.0) / (Float256::from(360.0) * x2 * xf)
                    + Float256::from(1.0) / (Float256::from(1260.0) * x2 * x2 * xf)
                    - Float256::from(1.0) / (Float256::from(1680.0) * x2 * x2 * x2 * xf);
                let r = loggamma(c(x, 0.0), A_REF).unwrap();
                let e = rel_err(r.re, stirling);
                eprintln!("x={x}: {:.1} digits", digits(e));
                assert!(e < 1e-18, "x={x}: {e}");
            }
        }

        #[test]
        fn large_imaginary_part_does_not_overflow() {
            for z in [c(-0.5, 10000.0), c(-1.5, 5000.0), c(0.3, 10000.0)] {
                let r = ln_gamma(z, 20).unwrap();
                assert!(r.re.is_finite() && r.im.is_finite(), "{z}");
            }
        }

        #[test]
        fn continuity_along_vertical_path() {
            let h = 1e-3;
            let mut prev = loggamma(c(0.5, 0.0), 30).unwrap();
            let mut y = 0.0;
            for _ in 0..200 {
                y += h * 10.0;
                let cur = loggamma(c(0.5, y), 30).unwrap();
                let d = cur.sub(prev);
                assert!(to64(d.re.abs()) < 0.5 && to64(d.im.abs()) < 0.5, "y={y}");
                prev = cur;
            }
            // Across Re z = 0 into the left half-plane, off the cut.
            let mut prev = loggamma(c(0.05, 2.0), 30).unwrap();
            for i in 1..=40 {
                let x = 0.05 - 0.005 * i as f64;
                let cur = loggamma(c(x, 2.0), 30).unwrap();
                let d = cur.sub(prev);
                assert!(to64(d.re.abs()) < 0.5 && to64(d.im.abs()) < 0.5, "x={x}");
                prev = cur;
            }
        }

        #[test]
        fn error_handling() {
            for n in [0.0, -1.0, -2.0, -10.0] {
                assert_eq!(loggamma(c(n, 0.0), A_REF), Err(MathError::Pole));
                assert_eq!(ln_gamma(c(n, 0.0), A_REF), Err(MathError::Pole));
            }
            for a in [0, 1] {
                assert_eq!(
                    loggamma(c(1.0, 0.0), a),
                    Err(MathError::ParameterOutOfRange)
                );
                assert_eq!(
                    ln_gamma(c(1.0, 0.0), a),
                    Err(MathError::ParameterOutOfRange)
                );
            }
            assert_eq!(ln_gamma(c(10001.0, 0.0), 12), Err(MathError::Overflow));
            assert_eq!(ln_gamma(c(0.0, 10001.0), 12), Err(MathError::Overflow));
            assert_eq!(ln_gamma(c(-10001.0, 0.0), 12), Err(MathError::Underflow));
            assert_eq!(loggamma(c(10001.0, 0.0), 12), Err(MathError::Overflow));
            assert_eq!(loggamma(c(-10001.0, 0.0), 12), Err(MathError::Underflow));
        }

        proptest! {
            #![proptest_config(ProptestConfig::with_cases(24))]

            #[test]
            fn prop_recurrence(x in 0.1_f64..20.0, y in -10.0_f64..10.0) {
                let z = c(x, y);
                let lhs = loggamma(z.add(c(1.0, 0.0)), 30).unwrap()
                    .sub(loggamma(z, 30).unwrap())
                    .sub(z.ln());
                prop_assert!(to64(lhs.re.abs()) < 1e-15);
                prop_assert!(to64(angle_diff(lhs.im, Float256::from(0.0))) < 1e-15);
            }

            #[test]
            fn prop_reflection(x in 0.05_f64..0.95, y in -5.0_f64..5.0) {
                let z = c(x, y);
                let pi_c = Complex256::new(PI, Float256::from(0.0));
                let lhs = ln_gamma(z, 30).unwrap()
                    .add(ln_gamma(c(1.0, 0.0).sub(z), 30).unwrap());
                let rhs = Complex256::new(PI.ln(), Float256::from(0.0))
                    .sub(pi_c.mul(z).sin().ln());
                prop_assert!(to64((lhs.re - rhs.re).abs()) < 1e-15);
                prop_assert!(to64(angle_diff(lhs.im, rhs.im)) < 1e-15);
            }
        }
    }
}
