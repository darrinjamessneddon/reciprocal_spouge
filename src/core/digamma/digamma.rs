pub mod digamma {
  use crate::core::{spouge_coefficients, Complex256, ComplexOps, MathError};

  pub fn digamma_c256(z: Complex256, a: usize) -> Result<Complex256, MathError> {
    // for now return placeholder
    return Complex256::new(Float256::from(0.0), Float256::from(0.0));
  }
