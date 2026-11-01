use std::env;

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-env-changed=SPOUGE_A");

    // Determine the value of parameter 'a' from the active Cargo features
    let a = if env::var("CARGO_FEATURE_A_50").is_ok() {
        50
    } else if env::var("CARGO_FEATURE_A_75").is_ok() {
        75
    } else {
        100
    };

    println!("cargo:rustc-env=SPOUGE_A={a}");
}
