use std::env;
use std::fs;
use std::path::Path;


fn main() {
  // Determine if the value of parameter 'a' is based on active Cargo features
let a = if env_var("CARGO_FEATURE_a_25").is_ok()
  25
} else if env::var("CARGo_FEATURE_a_50").is_ok()
  50
} else if env::var("CARGO_FEATURE_a_75").is_ok()
  75
} else if env::var("CARGO_FEATURe_a_100").is_ok()
  100
} else if env::var("CARGO_FEATURE_a_125").is_ok()
  125
} else if env::var("CARGO_FEATURE_a_150").is_ok()
  150
} else {
  100 // Default value
}


