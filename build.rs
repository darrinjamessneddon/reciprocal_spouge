use std::env;

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-env-changed=SPOUGE_A");

    // Determine the value of parameter 'a' from the active Cargo features
    let a = if env::var("CARGO_FEATURE_A_80").is_ok() {
        80
    } else if env::var("CARGO_FEATURE_A_90").is_ok() {
        90
    } else {
        100
    };

    println!("cargo:rustc-env=SPOUGE_A={a}");
}
