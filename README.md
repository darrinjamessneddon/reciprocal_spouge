A Rust library designed to directly calculate the reciprocal gamma function using a modified version of Spouge's
approximation for the gamma function.
The inputs for the Spouge approximation are a complex number z, and the Spouge parameter 'a'.
The Spouge parameter is an integer greater than 2. (Generally, in practice it is better that it be at least 7
or more. This parameter can be selected by the user.

Examples

let z = Complex64::new(5.0, 0.0);
let spouge_rg = RSpouge::new(z, 10);
let result = spouge_rg.compute();
println!("Approximate reciprocal gamma for {}: {}", z, result);

let gamma_result = Complex64::new(1.0, 0.0) / result;
println!("Approximate gamma function value for {}: {}", z, gamma_result);


let z = Complex64::new(2.5, 0.0);
let spouge_rg = RSpouge::new(z, 21);
let result = spouge_rg.compute();
println!("Reciprocal gamma for {}: {}", z, result);

let gamma_result = Complex64::new(1.0, 0.0) / result;
println!("Approximate gamma function value for {}: {}", z, gamma_result);

// or to calculate the reciprocal gamma function using an automatically generated Spouge parameter.

let z = Complex64::new(3.5, 1.5);
let reciprocal = RSpouge::with_auto_a(z);
let result = reciprocal.compute()
println!("Reciprocal gamma for {}: {}", z, result);

let gamma_result = Complex64::new(1.0, 0.0) / result;
println!("Gamma function result for {}: {}", z, gamma_result);


The library also contains a function to automatically generate the Spouge parameter.
