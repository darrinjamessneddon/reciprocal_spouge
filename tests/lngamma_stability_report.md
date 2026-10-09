# `ln_gamma` and `loggamma` Numerical Stability Test Report

The deterministic stability tests live in the `stability_tests` module in
`src/core/lngamma/lngamma.rs`. Run them with:

```sh
cargo test --lib stability_tests
```

## Coverage and rationale

- Positive real inputs from `1e-30` through `9999.5` cover the high-curvature
  region near the pole at zero, ordinary values, and large results close to the
  supported input limit. Both functions are required to return finite values
  and agree on this shared branch.
- Negative non-integers `-0.5`, `-1.5`, and `-10.5` exercise reflection and
  branch handling. Their log values must remain finite, and exponentiating each
  result must recover the same gamma value even when the log branches differ.
- The pair `1e-30` and `2e-30` checks sensitivity near zero: doubling the
  positive input should change log-gamma by approximately `-ln(2)`, for both
  functions.
- Existing tests in the same module additionally cover high-precision real and
  complex reference values, values on both sides of negative-integer poles,
  large real and imaginary inputs, recurrence, reflection, conjugate symmetry,
  and overflow/underflow handling.

## Observed outcomes

The focused `stability_tests` suite passes. Outputs were finite for all
non-pole cases; both functions agreed across the tested positive range, their
exponentials agreed for the tested negative non-integers, and the near-zero
perturbation followed the expected `-ln(2)` change within the test tolerance.
