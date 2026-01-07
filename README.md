WELCOME!!...AND THANK YOU FOR VISITING THIS PROJECT

THIS IS MY FIRST PROJECT AND I NOT ONLY WELCOME COLLABORATORS (SHOULD YOU WISH TO GET INVOLVED),
I ALSO LOOK FORWARD TO ANY FEEDBACK ON HOW TO IMPROVE THINGS.

THIS PROJECT IS VERY MUCH A WORK IN PROGRESS AND IS NOT YET COMPLETE

A Rust library designed to directly calculate the reciprocal gamma function using a modified version of Spouge's
approximation for the gamma function.
The inputs for the Spouge approximation are a complex number z, and the Spouge parameter 'a'.
The Spouge parameter is an integer greater than 2. 
The library also has a function to calculate the reciprocal gamma function using an automatically generated Spouge
parameter.

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

Issues to be addressed include:

(1)  Documentation
(2)  Improving the accuracy of reciprocal gamma calculation for complex numbers where z is large.
(3)  Extending the range of complex numbers for which the reciprocal gamma value can be calculated.
(4)  Only the real part of a complex number is used in the function to automatically calculate 'a'
therefore this function may need to be tweaked if the imaginary part of the complex number has a substantial
impact on whether the optimal value for 'a' is being chosen. Otherwise it is best to keep the function as simple
as possible.
(5) Optimization of the computation of reciprocal gamma (if this can be done with relatively simple functions or algorithms.
(6)  Testing


