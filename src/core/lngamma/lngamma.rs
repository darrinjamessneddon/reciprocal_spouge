pub mod lngamma {
    // The ln_gamma function for complex z is not continuous. It is important to consider zeros, singularities and discontinuities
    // of the function before creating it.

    // ln_gamma(z + 1) = (z + 0.5)* ln(z + a) - (z + a) + ln {co + (sum from k = 1 to a - 1 (ck/(z - 1 + k))

    // If the real part of z is small Re(z) < 0.5, we need to use the reflection formula:
    // ln(z) = ln(pi) - lnsin(pi * z) - ln_gamma(1 - z) before applying the approximation.

    pub fn ln_gamma() {
        // add code here
    }
}
