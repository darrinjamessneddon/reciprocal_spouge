pub mod complex; // a submodule providing structs and functions for complex numbers.
pub mod error; // a submodule providing error handling for the library.
pub mod gamma; // a submodule holding the spouge approximation of the gamma function.
pub mod lngamma;
pub mod rgamma; // a submodule holding a spouge-based implementation of the reciprocal gamma function.
pub mod shared; // a submodule providing functions that are shared between the other modules. // a submodule holding an implementation of the natural logarithm of the gamma function.

pub use complex::{Complex256, ComplexComparisons, ComplexOps};
pub use error::{ComplexError, MathError, ParameterOutOfRangeError, ParseError, ParseResult, Result};
pub use gamma::{spouge, spouge_c256};
pub use lngamma::ln_gamma;
pub use rgamma::{rspouge, rspouge_c256};
pub use shared::{factorial, spouge_coefficients};
