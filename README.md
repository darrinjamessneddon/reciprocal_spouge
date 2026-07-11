WELCOME!!...AND THANK YOU FOR VISITING THIS PROJECT

THIS IS MY FIRST PROJECT AND I NOT ONLY WELCOME COLLABORATORS (SHOULD YOU WISH TO GET INVOLVED),
I ALSO LOOK FORWARD TO ANY FEEDBACK ON HOW TO IMPROVE THINGS.

THIS PROJECT IS VERY MUCH A WORK IN PROGRESS AND IS NOT YET COMPLETE.

This project is a Rust library.

The goals of this project (at this preliminary stage) are as follows:
(1)  To use the f256 crate to create structs to enable the computation of special functions in the complex plane.
(2)  To facillitate the calculation of gamma and reciprocal gamma using f256 complex numbers.
(3)  To create a function to approximation the natural logarithm of the reciprocal gamma function.
(4)  To use Spouge's approximation to calculate values for the gamma function for z values with a modulus of <= 10000. (This is an arbitrarily chosen value).
(5)  To use the adapted Antoulas-Anderson (AAA) algorithm to further extend the range of complex numbers for which the values of these special functions can be computed.
(6)  This library is intended not only to be used as a stand alone library, but should also be able to be used as a dependency for other libaries
seeking to optimize calculations of the functions mentioned above.

The inputs for the Spouge approximation are a complex number z, and the Spouge parameter 'a'.
The Spouge parameter is an integer greater than 2.

The AAA algorithm has yet to be written (and may need to be in a separate module, or crate). Help with this aspect of the project would be particularly welcomed.







