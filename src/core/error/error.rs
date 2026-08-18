pub mod error {


    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub enum MathError {
        ParameterOutOfRange,
        DivisionByZero,
        Overflow,
        Underflow,
    }

    impl MathError {
        pub fn message(&self) -> &'static str {
            match self {
                Self::ParameterOutOfRange => "Parameter 'a' is out of range",
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
        pub max: u64,
    }
    impl ParameterOutOfRangeError {
        pub fn new(parameter: &'static str, min: u64, max: u64) -> Self {
            Self {
                parameter,
                min,
                max,
            }
        }

        pub fn error(a: u64) -> Option<Self> {
            let parameter = "a";
            let min = 2;
            let max = 21;

            if a < min || a > max {
                Some(Self::new(parameter, min, max))
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
}