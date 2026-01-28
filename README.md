WELCOME!!...AND THANK YOU FOR VISITING THIS PROJECT

THIS IS MY FIRST PROJECT AND I NOT ONLY WELCOME COLLABORATORS (SHOULD YOU WISH TO GET INVOLVED),
I ALSO LOOK FORWARD TO ANY FEEDBACK ON HOW TO IMPROVE THINGS.

THIS PROJECT IS VERY MUCH A WORK IN PROGRESS AND IS NOT YET COMPLETE.

This project is a Rust library.

The goals of this project (at this preliminary stage) are as follows:
(1)  To use the f256 crate to create structs to enable the computation of special functions in the complex plane.
(2)  To use a modified version of Spouge's approximation to directly calculate values for the reciprocal gamma function for z values
with a modulus of <= 10000. To do this requires a simple rearrangement of Spouge's approximation formula.
(3)  To create a function to approximation the natural logarithm of the reciprocal gamma function.
(4)  To use Spouge's approximation to calculate values for the gamma function for z values with a modulus of <= 10000.
(5)  To have an accuracy of better than 1e-15. (The greater the accuracy that can be achieved without combining the Spouge approximation
with other methods at this stage, the better).
(6)  This library is intended not only to be used as a stand alone library, but should also be able to be used as a dependency for other libaries
seeking to optimize calculations of the functions mentioned above.

The inputs for the Spouge approximation are a complex number z, and the Spouge parameter 'a'.
The Spouge parameter is an integer greater than 2.

It was initially intended to attain a very high degree of accuracy, yet I have found it necessary to break those initial plans into smaller,
more easily achieved goals first.







