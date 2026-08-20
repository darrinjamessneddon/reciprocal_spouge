Welcome and thank you for taking an interest in this project.

This project is still very much a work in progress. At this point in time (the lngamma (log_gamma) submodule is yet to be created), and the error-handling module has a lot of work to be done on it.

The original (unpublished) repository for this project was deleted, and the project re-built from scratch.
I humbly apologise if this has caused you any inconvenience.

This is a library created using Rust programming language. Rust is a language renowned for memory safety and concurrency.
The library uses the Spouge approximation because it is numerically stable.
The library is built on top of the f256 crate, and deals with special functions including the reciprocal gamma, gamma and log-gamma functions in the complex plane.

The Spouge approximation takes two inputs: the value for z, and the value of a parameter 'a'. The value of this parameter should be a > 2.

The goals of the project are to allow the computation of the special functions mentioned above, giving a greater degree of precision than could be achieved by simply using complex numbers based on f64.

To achieve this, a Complex256 struct was created to work with complex numbers using f256 values.
