use f256::f256 as Float256;
use proptest::prelude::*;
use reciprocal_spouge::core::{
    Complex256, Complex256 as CoreComplex256, ComplexOps, ComplexOps as CoreComplexOps, MathError,
    ln_gamma as core_ln_gamma, rspouge as core_rspouge, rspouge_c256 as core_rspouge_c256,
    spouge as core_spouge, spouge_c256 as core_spouge_c256, spouge_coefficients,
};
use reciprocal_spouge::{
    Complex256 as RootComplex256, ComplexOps as RootComplexOps, ln_gamma as root_ln_gamma,
    rspouge as root_rspouge, rspouge_c256 as root_rspouge_c256, spouge as root_spouge,
    spouge_c256 as root_spouge_c256,
};
use std::str::FromStr;

/// Integration-test starter template for the current public API.
///
/// Notes for future test authors:
/// - Prefer approximate comparisons over exact equality when working with `f256` results.
/// - Keep the Spouge parameter `a` fixed inside a test so that failures are easier to reproduce.
/// - Start with identities that are mathematically stable, then add wider-input coverage.
/// - Use `prop_assume!` in property tests to skip poles and other unsupported regions.
const DEFAULT_SPOUGE_A: usize = 12;

fn c256(re: f64, im: f64) -> Complex256 {
    Complex256::from_f64(re, im)
}

fn assert_f256_close(actual: Float256, expected: Float256, tolerance: Float256) {
    let difference = (actual - expected).abs();
    assert!(
        difference <= tolerance,
        "expected {actual} to be within {tolerance} of {expected}, diff={difference}"
    );
}

fn assert_complex_close(actual: Complex256, expected: Complex256, tolerance: Float256) {
    assert_f256_close(actual.re, expected.re, tolerance);
    assert_f256_close(actual.im, expected.im, tolerance);
}

#[test]
fn complex_helper_example_uses_tolerances_for_high_precision_values() {
    // This is a small real test so the template is exercised by CI.
    // It also demonstrates the preferred pattern for comparing f256-based values.
    let value = c256(3.0, 4.0);
    assert_f256_close(
        value.magnitude(),
        Float256::from(5.0),
        Float256::from(1e-12),
    );
    assert_complex_close(value.conj(), c256(3.0, -4.0), Float256::from(0.0));
}

#[test]
fn complex_checked_add_returns_some_for_finite_inputs() {
    let left = c256(1.5, -2.0);
    let right = c256(2.25, 3.0);

    let sum = left
        .checked_add(right)
        .expect("checked_add should succeed for finite components");

    assert_complex_close(sum, c256(3.75, 1.0), Float256::from(1e-30));
}

#[test]
fn complex_checked_sub_returns_some_for_finite_inputs() {
    let left = c256(5.25, -1.5);
    let right = c256(2.0, 4.0);

    let difference = left
        .checked_sub(right)
        .expect("checked_sub should succeed for finite components");

    assert_complex_close(difference, c256(3.25, -5.5), Float256::from(1e-30));
}

#[test]
fn complex_checked_add_returns_none_for_non_finite_result() {
    let left = Complex256::new(Float256::INFINITY, Float256::from(1.0));
    let right = c256(1.0, 2.0);

    assert!(left.checked_add(right).is_none());
}

#[test]
fn complex_checked_sub_returns_none_for_non_finite_result() {
    let left = c256(1.0, 2.0);
    let right = Complex256::new(Float256::from(0.5), Float256::INFINITY);

    assert!(left.checked_sub(right).is_none());
}

#[test]
fn flattened_reexports_preserve_access_to_existing_core_api() {
    let root_value = RootComplex256::from_f64(3.0, 4.0);
    let core_value = CoreComplex256::from_f64(3.0, 4.0);

    assert_f256_close(
        RootComplexOps::magnitude(root_value),
        Float256::from(5.0),
        Float256::from(1e-12),
    );
    assert_complex_close(
        RootComplexOps::conj(RootComplexOps::conj(root_value)),
        Complex256::from_f64(3.0, 4.0),
        Float256::from(0.0),
    );
    assert_complex_close(
        CoreComplexOps::conj(CoreComplexOps::conj(core_value)),
        Complex256::from_f64(3.0, 4.0),
        Float256::from(0.0),
    );

    let input = c256(2.5, 0.25);

    assert_complex_close(
        root_spouge_c256(input, DEFAULT_SPOUGE_A).unwrap(),
        core_spouge_c256(input, DEFAULT_SPOUGE_A).unwrap(),
        Float256::from(0.0),
    );
    assert_complex_close(
        root_rspouge_c256(input, DEFAULT_SPOUGE_A).unwrap(),
        core_rspouge_c256(input, DEFAULT_SPOUGE_A).unwrap(),
        Float256::from(0.0),
    );
    assert_eq!(
        root_spouge(input, DEFAULT_SPOUGE_A).unwrap(),
        core_spouge(input, DEFAULT_SPOUGE_A).unwrap(),
    );
    assert_eq!(
        root_rspouge(input, DEFAULT_SPOUGE_A).unwrap(),
        core_rspouge(input, DEFAULT_SPOUGE_A).unwrap(),
    );
    assert_complex_close(
        root_ln_gamma(input, DEFAULT_SPOUGE_A).unwrap(),
        core_ln_gamma(input, DEFAULT_SPOUGE_A).unwrap(),
        Float256::from(0.0),
    );
}

#[test]
fn spouge_api_returns_parameter_out_of_range_for_small_a() {
    let input = c256(2.5, 0.25);

    assert_eq!(spouge_coefficients(1), Err(MathError::ParameterOutOfRange));
    assert_eq!(
        root_spouge_c256(input, 1),
        Err(MathError::ParameterOutOfRange)
    );
    assert_eq!(
        root_rspouge_c256(input, 1),
        Err(MathError::ParameterOutOfRange)
    );
    assert_eq!(root_spouge(input, 1), Err(MathError::ParameterOutOfRange));
    assert_eq!(root_rspouge(input, 1), Err(MathError::ParameterOutOfRange));
}

#[test]
fn spouge_api_still_supports_successful_coefficient_paths() {
    let gamma_of_one = root_spouge_c256(c256(1.0, 0.0), DEFAULT_SPOUGE_A).unwrap();
    let reciprocal_gamma_of_one = root_rspouge_c256(c256(1.0, 0.0), DEFAULT_SPOUGE_A).unwrap();
    let gamma_string = root_spouge(c256(1.0, 0.0), DEFAULT_SPOUGE_A).unwrap();
    let reciprocal_gamma_string = root_rspouge(c256(1.0, 0.0), DEFAULT_SPOUGE_A).unwrap();

    assert_complex_close(gamma_of_one, c256(1.0, 0.0), Float256::from(1e-15));
    assert_complex_close(
        reciprocal_gamma_of_one,
        c256(1.0, 0.0),
        Float256::from(1e-15),
    );
    assert_eq!(gamma_string, gamma_of_one.to_string());
    assert_eq!(reciprocal_gamma_string, reciprocal_gamma_of_one.to_string());
}

#[test]
fn coefficient_path_source_is_free_of_unwrap_and_expect() {
    for source in [
        include_str!("../src/core/shared/shared.rs"),
        include_str!("../src/core/gamma/gamma.rs"),
        include_str!("../src/core/rgamma/rgamma.rs"),
    ] {
        assert!(!source.contains("unwrap("));
        assert!(!source.contains("expect("));
    }
}

#[test]
fn gamma_matches_known_reference_values_and_rejects_poles() {
    let gamma_of_one = root_spouge_c256(c256(1.0, 0.0), DEFAULT_SPOUGE_A).unwrap();
    let gamma_of_five = root_spouge_c256(c256(5.0, 0.0), DEFAULT_SPOUGE_A).unwrap();

    assert_complex_close(gamma_of_one, c256(1.0, 0.0), Float256::from(1e-20));
    assert_complex_close(gamma_of_five, c256(24.0, 0.0), Float256::from(1e-12));
    assert_eq!(
        root_spouge_c256(c256(0.0, 0.0), DEFAULT_SPOUGE_A),
        Err(MathError::Pole)
    );
    assert_eq!(
        root_spouge_c256(c256(-2.0, 0.0), DEFAULT_SPOUGE_A),
        Err(MathError::Pole)
    );
}

#[test]
fn reciprocal_gamma_matches_inverse_relationship_and_zero_contract() {
    for pole in [c256(0.0, 0.0), c256(-1.0, 0.0), c256(-2.0, 0.0)] {
        assert_complex_close(
            root_rspouge_c256(pole, DEFAULT_SPOUGE_A).unwrap(),
            c256(0.0, 0.0),
            Float256::from(0.0),
        );
    }

    for input in [c256(2.5, 0.25), c256(-0.5, 0.0)] {
        let gamma = root_spouge_c256(input, DEFAULT_SPOUGE_A).unwrap();
        let reciprocal_gamma = root_rspouge_c256(input, DEFAULT_SPOUGE_A).unwrap();
        let product = gamma.mul(reciprocal_gamma);

        assert_complex_close(product, c256(1.0, 0.0), Float256::from(1e-12));
    }
}

#[test]
fn log_gamma_tracks_principal_logarithm_contract() {
    let log_gamma_of_one = root_ln_gamma(c256(1.0, 0.0), DEFAULT_SPOUGE_A).unwrap();
    let input = c256(2.5, 0.25);
    let log_gamma = root_ln_gamma(input, DEFAULT_SPOUGE_A).unwrap();
    let reconstructed_gamma = log_gamma.exp();

    assert_complex_close(log_gamma_of_one, c256(0.0, 0.0), Float256::from(1e-20));
    assert_complex_close(
        reconstructed_gamma,
        root_spouge_c256(input, DEFAULT_SPOUGE_A).unwrap(),
        Float256::from(1e-12),
    );
    assert_eq!(
        root_ln_gamma(c256(0.0, 0.0), DEFAULT_SPOUGE_A),
        Err(MathError::Pole)
    );
}

#[test]
fn complex_operations_cover_public_api_examples() {
    let z = c256(1.25, -0.75);
    let reciprocal = z.recip();
    let product = z.mul(reciprocal);
    let parsed_from_display = Complex256::from_str(&z.to_string()).unwrap();
    let parsed_from_parenthesized = Complex256::from_string(&format!("({})", z)).unwrap();

    assert_complex_close(product, c256(1.0, 0.0), Float256::from(1e-18));
    assert_complex_close(parsed_from_display, z, Float256::from(0.0));
    assert_complex_close(parsed_from_parenthesized, z, Float256::from(0.0));
    assert_complex_close(z.add(z.conj()), c256(2.5, 0.0), Float256::from(1e-18));
}

proptest! {
    #[test]
    #[ignore = "Template placeholder: enable after choosing property-test domains and tolerances"]
    fn gamma_and_reciprocal_gamma_stay_consistent_under_proptest(
        re in -4.5f64..4.5,
        im in -4.5f64..4.5,
    ) {
        // TODO: Narrow or widen these domains as the implementation matures.
        // TODO: Skip poles and numerically fragile regions with `prop_assume!`.
        // TODO: Record any counterexamples here as deterministic regression tests.
        prop_assume!(im.abs() > 1e-6 || re.fract().abs() > 1e-6 || re >= 0.0);

        let z = c256(re, im);
        let gamma = root_spouge_c256(z, DEFAULT_SPOUGE_A).unwrap();
        let reciprocal_gamma = root_rspouge_c256(z, DEFAULT_SPOUGE_A).unwrap();
        let product = gamma.mul(reciprocal_gamma);

        assert_f256_close(product.re, Float256::from(1.0), Float256::from(1e-12));
        assert_f256_close(product.im, Float256::from(0.0), Float256::from(1e-12));
    }
}
