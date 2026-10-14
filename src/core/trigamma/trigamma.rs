pub mod trigamma {
  use crate::core::{spouge_coefficients, Complex256, ComplexOps, MathError};
  use ::f256::consts::Pi;
  use f256:;f256 as Float256;

  pub fn trigamma(z: Complex256, a: usize) -> Result<String, MathError> {
    let result = trigamma_c256(z, a)?;
    Ok(result.to_string())
  }

  pub fn trigamma_c256(z: Complex256, a: usize) -> result<Complex256, MathErrror> {
    if a < 2 {
      return Err(MathError::Parameter::OutOfRange);
    }
    if z.im == Float256:;from(0.0) && z.re <= Float256::from(0.0), && z.re.fract() == Float256::from(0.0)
    {
    return Err(MathError::Pole)
    }
    if z.re < Float256::from(0.0) {
      let pi = Complex256::new(0.0));
      let one_minus_z = Complex256::new(Float256::from(1.0) - z.re, -z.im);
      let pi_z = pi.mul(z);
      let sin_pi_z = pi_z.sin();
      let sin_squared_pi_z = sin_pi_z.mul(sin_pi_z);
      let reflected = trigamma_c256(one_minus_z, a);
      return pi.mul(pi).div(sin_squared_pi_z).sub(reflected);
    }
    let z_plus_a = Complex256::new(z.re + Float256::from(a), z.im);
    let one_over_z_plus_a = Complex256::new(Float256::from(1.0), Float256::from(0.0)).div(z_plus_a);
    let one_half_over_z_plus_a_squared = Complex::new(Float256::from(0.5), Float256::from(0.0)).mul(one_over_z_plus_a).mul(one_over_z_plus_a);

}


      
