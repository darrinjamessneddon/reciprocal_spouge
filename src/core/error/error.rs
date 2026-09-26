pub mod error {

    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub enum MathError {
        ParameterOutOfRange,
        Pole,
        DivisionByZero,
        Overflow,
        Underflow,
    }

    impl MathError {
        pub fn message(&self) -> &'static str {
            match self {
                Self::ParameterOutOfRange => "Parameter 'a' is out of range",
                Self::Pole => "Function is undefined at this pole",
                Self::DivisionByZero => "Attempted division by zero",
                Self::Overflow => "Arithmetic overflow",
                Self::Underflow => "Arithmetic underflow",
            }
        }
    }

    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub struct ParameterOutOfRangeError {
        pub parameter: &'static str,
        pub min: u64,
    }
    impl ParameterOutOfRangeError {
        pub fn new(parameter: &'static str, min: u64) -> Self {
            Self { parameter, min }
        }

        pub fn error(a: u64) -> Option<Self> {
            let parameter = "a";
            let min = 2;

            if a < min {
                Some(Self::new(parameter, min))
            } else {
                None
            }
        }
    }

    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub enum ParseError {
        InvalidFormat,
        UnexpectedCharacter,
        InvalidType,
        UnknownError,
    }
    impl ParseError {
        pub fn message(&self) -> &'static str {
            match self {
                Self::InvalidFormat => "Invalid format",
                Self::UnexpectedCharacter => "Unexpected character",
                Self::InvalidType => "Invalid type",
                Self::UnknownError => "Unknown error",
            }
        }
    }
    pub type Result<T> = std::result::Result<T, MathError>;
    pub type ParseResult<T> = std::result::Result<T, ParseError>;

    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub enum ComplexError {
        DivisionByZero,
        Overflow,
        Underflow,
        ParseError(ParseError),
        OutOfRange,
        UnknownError,
    }
    impl ComplexError {
        pub fn message(&self) -> &'static str {
            match self {
                Self::DivisionByZero => "Attempted division by zero",
                Self::Overflow => "Arithmetic overflow",
                Self::Underflow => "Arithmetic underflow",
                Self::ParseError(err) => err.message(),
                Self::OutOfRange => "Complex number is out of range",
                Self::UnknownError => "Unknown error",
            }
        }
    }

    #[cfg(test)]
    mod tests {
        use super::ParameterOutOfRangeError;

        #[test]
        fn parameter_out_of_range_error_only_rejects_values_below_minimum() {
            assert_eq!(
                ParameterOutOfRangeError::error(1),
                Some(ParameterOutOfRangeError::new("a", 2))
            );
            assert_eq!(ParameterOutOfRangeError::error(2), None);
            assert_eq!(ParameterOutOfRangeError::error(200), None);
            assert_eq!(ParameterOutOfRangeError::error(u64::MAX), None);
        }
    }
}
