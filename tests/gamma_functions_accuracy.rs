use f256::f256 as Float256;
use reciprocal_spouge::{
    ln_gamma, rspouge, rspouge_c256, spouge, spouge_c256, Complex256, ComplexOps, MathError,
};
use std::str::FromStr;

const A: usize = 30;
const DEFAULT_TOLERANCE: f64 = 1e-12;

fn c(re: f64, im: f64) -> Complex256 {
    Complex256::from_f64(re, im)
}

fn assert_close(actual: Complex256, expected: Complex256, rel_tol: f64, label: &str) {
    let diff = actual.sub(expected).magnitude();
    let scale = Float256::from(1.0).max(expected.magnitude());
    let allowed = scale * Float256::from(rel_tol);
    assert!(
        diff <= allowed,
        "{label}: expected {expected}, got {actual}, diff={diff}, allowed={allowed}"
    );
}

fn real(x: f64) -> Complex256 {
    c(x, 0.0)
}

#[test]
fn spouge_matches_known_real_values() {
    let sqrt_pi = std::f64::consts::PI.sqrt();
    let cases = [
        (1.0, 1.0),
        (2.0, 1.0),
        (5.0, 24.0),
        (10.0, 362880.0),
        (0.5, sqrt_pi),
        (1.5, sqrt_pi / 2.0),
    ];
    for (x, expected) in cases {
        let actual = spouge_c256(real(x), A).unwrap();
        assert_close(actual, real(expected), 1e-12, &format!("gamma({x})"));
    }
}

#[test]
fn spouge_handles_negative_non_integers() {
    let sqrt_pi = std::f64::consts::PI.sqrt();
    let actual = spouge_c256(real(-0.5), A).unwrap();
    assert_close(actual, real(-2.0 * sqrt_pi), 1e-12, "gamma(-0.5)");
    let actual = spouge_c256(real(-1.5), A).unwrap();
    assert_close(actual, real(4.0 * sqrt_pi / 3.0), 1e-12, "gamma(-1.5)");
}

#[test]
fn spouge_matches_known_complex_value() {
    let actual = spouge_c256(c(1.0, 1.0), A).unwrap();
    let expected = c(0.498_015_668_118_356, -0.154_949_828_301_810_7);
    assert_close(actual, expected, 1e-12, "gamma(1+i)");
}

#[test]
fn spouge_satisfies_recurrence_relation() {
    for z in [c(0.5, 0.0), c(2.3, 0.7), c(3.1, -1.2), c(-0.5, 0.25)] {
        let lhs = spouge_c256(z.add(real(1.0)), A).unwrap();
        let rhs = z.mul(spouge_c256(z, A).unwrap());
        assert_close(lhs, rhs, DEFAULT_TOLERANCE, "gamma(z+1) = z*gamma(z)");
    }
}

#[test]
fn spouge_string_matches_complex256_output() {
    for z in [real(5.0), c(2.5, 0.5)] {
        let text = spouge(z, A).unwrap();
        let value = spouge_c256(z, A).unwrap();
        assert_eq!(text, value.to_string());
        let parsed = Complex256::from_str(&text).unwrap();
        assert_close(parsed, value, DEFAULT_TOLERANCE, "parsed gamma string");
    }
}

#[test]
fn spouge_reports_errors() {
    assert_eq!(spouge_c256(real(0.0), A), Err(MathError::Pole));
    assert_eq!(spouge_c256(real(-3.0), A), Err(MathError::Pole));
    assert_eq!(
        spouge_c256(real(1.0), 1),
        Err(MathError::ParameterOutOfRange)
    );
    assert_eq!(spouge(real(1.0), 0), Err(MathError::ParameterOutOfRange));
    assert_eq!(spouge_c256(real(20000.0), A), Err(MathError::Overflow));
}

#[test]
fn spouge_accuracy_improves_with_larger_a() {
    let z = real(5.0);
    let error = |a: usize| spouge_c256(z, a).unwrap().sub(real(24.0)).magnitude();
    let small = error(4);
    let large = error(30);
    assert!(
        large < small,
        "expected error with a=30 ({large}) to be below error with a=4 ({small})"
    );
    assert!(large <= Float256::from(1e-12));
}

#[test]
fn rspouge_matches_known_real_values() {
    let sqrt_pi = std::f64::consts::PI.sqrt();
    let cases = [
        (1.0, 1.0),
        (5.0, 1.0 / 24.0),
        (10.0, 1.0 / 362880.0),
        (0.5, 1.0 / sqrt_pi),
    ];
    for (x, expected) in cases {
        let actual = rspouge_c256(real(x), A).unwrap();
        assert_close(actual, real(expected), 1e-12, &format!("1/gamma({x})"));
    }
}

#[test]
fn rspouge_is_zero_at_zero_and_negative_integers() {
    for x in [0.0, -1.0, -2.0, -7.0] {
        let actual = rspouge_c256(real(x), A).unwrap();
        assert_eq!(actual.re, Float256::from(0.0), "1/gamma({x}) real part");
        assert_eq!(
            actual.im,
            Float256::from(0.0),
            "1/gamma({x}) imaginary part"
        );
    }
}

#[test]
fn rspouge_matches_known_complex_value() {
    let actual = rspouge_c256(c(1.0, 1.0), A).unwrap();
    let expected = spouge_c256(c(1.0, 1.0), A).unwrap().recip();
    assert_close(actual, expected, DEFAULT_TOLERANCE, "1/gamma(1+i)");
    let reference = c(0.498_015_668_118_356, -0.154_949_828_301_810_7).recip();
    assert_close(actual, reference, 1e-12, "1/gamma(1+i) reference");
}

#[test]
fn rspouge_is_inverse_of_spouge() {
    for z in [
        real(0.5),
        real(4.0),
        c(2.5, 0.25),
        c(-0.5, 0.0),
        c(1.2, -3.0),
    ] {
        let product = spouge_c256(z, A).unwrap().mul(rspouge_c256(z, A).unwrap());
        assert_close(product, real(1.0), DEFAULT_TOLERANCE, "gamma(z)*rgamma(z)");
    }
}

#[test]
fn rspouge_string_matches_complex256_output() {
    for z in [real(5.0), c(2.5, 0.5)] {
        let text = rspouge(z, A).unwrap();
        let value = rspouge_c256(z, A).unwrap();
        assert_eq!(text, value.to_string());
        let parsed = Complex256::from_str(&text).unwrap();
        assert_close(parsed, value, DEFAULT_TOLERANCE, "parsed rgamma string");
    }
}

#[test]
fn rspouge_reports_errors() {
    assert_eq!(
        rspouge_c256(real(1.0), 1),
        Err(MathError::ParameterOutOfRange)
    );
    assert_eq!(rspouge(real(1.0), 0), Err(MathError::ParameterOutOfRange));
    assert_eq!(rspouge_c256(real(20000.0), A), Err(MathError::Overflow));
    assert_eq!(rspouge_c256(real(-20000.0), A), Err(MathError::Underflow));
}

#[test]
fn ln_gamma_matches_known_values() {
    let sqrt_pi = std::f64::consts::PI.sqrt();
    let cases = [
        (1.0, 0.0),
        (2.0, 0.0),
        (5.0, 24.0f64.ln()),
        (10.0, 362880.0f64.ln()),
        (0.5, sqrt_pi.ln()),
    ];
    for (x, expected) in cases {
        let actual = ln_gamma(real(x), A).unwrap();
        assert_close(actual, real(expected), 1e-12, &format!("ln gamma({x})"));
    }
}

#[test]
fn ln_gamma_exponential_recovers_gamma() {
    for z in [real(3.5), c(2.5, 0.25), c(1.5, 1.0), c(-0.5, 0.25)] {
        let recovered = ln_gamma(z, A).unwrap().exp();
        let direct = spouge_c256(z, A).unwrap();
        assert_close(recovered, direct, DEFAULT_TOLERANCE, "exp(ln gamma(z))");
    }
}

#[test]
fn ln_gamma_satisfies_recurrence() {
    for z in [real(0.5), real(3.2), c(2.3, 0.7), c(1.5, -1.0)] {
        let lhs = ln_gamma(z.add(real(1.0)), A).unwrap();
        let rhs = z.ln().add(ln_gamma(z, A).unwrap());
        assert_close(lhs, rhs, DEFAULT_TOLERANCE, "ln gamma(z+1)");
    }
}

#[test]
fn ln_gamma_handles_negative_real_part_on_principal_branch() {
    let sqrt_pi = std::f64::consts::PI.sqrt();
    let actual = ln_gamma(real(-0.5), A).unwrap();
    // gamma(-0.5) = -2*sqrt(pi) < 0, so the principal log has imaginary part +/- pi.
    assert_close(
        Complex256::new(actual.re, Float256::from(0.0)),
        real((2.0 * sqrt_pi).ln()),
        1e-12,
        "Re ln gamma(-0.5)",
    );
    let pi = Float256::from(std::f64::consts::PI);
    assert!(
        (actual.im.abs() - pi).abs() <= Float256::from(1e-12),
        "expected |Im ln gamma(-0.5)| = pi, got {}",
        actual.im
    );
    let recovered = actual.exp();
    assert_close(
        recovered,
        real(-2.0 * sqrt_pi),
        1e-12,
        "exp(ln gamma(-0.5))",
    );
}

#[test]
fn ln_gamma_reports_errors() {
    assert_eq!(ln_gamma(real(0.0), A), Err(MathError::Pole));
    assert_eq!(ln_gamma(real(-4.0), A), Err(MathError::Pole));
    assert_eq!(ln_gamma(real(1.0), 1), Err(MathError::ParameterOutOfRange));
    assert_eq!(ln_gamma(real(20000.0), A), Err(MathError::Overflow));
    assert_eq!(ln_gamma(real(-20000.0), A), Err(MathError::Underflow));
}
