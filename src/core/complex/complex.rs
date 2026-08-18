pub mod complex {

    use f256::f256 as Float256;
    use num_complex::Complex64;

    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct Complex256 {
        pub re: Float256,
        pub im: Float256,
    }
    impl Complex256 {
        pub fn new(re: Float256, im: Float256) -> Self {
            Complex256 { re, im }
        }
        pub fn to_f64(self) -> (f64, f64) {
            let re_f256 = self.re;
            let im_f256 = self.im;
            let re_f256_str = re_f256.to_string();
            let im_f256_str = im_f256.to_string();
            let re_f64 = re_f256_str.parse::<f64>().unwrap_or(0.0);
            let im_f64 = im_f256_str.parse::<f64>().unwrap_or(0.0);
            (re_f64, im_f64)
        }
        pub fn from_f64(re: f64, im: f64) -> Self {
            Complex256 {
                re: Float256::from(re),
                im: Float256::from(im),
            }
        }
        pub fn to_complex64(self) -> Complex64 {
            let re_f64: f64 = self.re.to_string().parse().unwrap_or(0.0);
            let im_f64: f64 = self.im.to_string().parse().unwrap_or(0.0);
            Complex64::new(re_f64, im_f64)
        }
        pub fn from_complex64(c: Complex64) -> Self {
            Complex256 {
                re: Float256::from(c.re),
                im: Float256::from(c.im),
            }
        }
        pub fn to_string(self) -> String {
            let re_str = self.re.to_string();
            let im_str = self.im.to_string();
            format!("({} + {}i)", re_str, im_str)
        }
        pub fn from_string(s: &str) -> Option<Self> {
            let s = s.trim();
            if !s.ends_with(')') {
                return None;
            }
            let s = &s[..s.len() - 1];
            let parts: Vec<&str> = s.split('+').collect();
            if parts.len() != 2 {
                return None;
            }
            let re = parts[0].trim().parse::<Float256>().ok()?;
            let im = parts[1].trim().trim_end_matches('i').parse::<Float256>().ok()?;
            Some(Complex256 { re, im })
        }
    }

    pub trait ComplexOps {
        fn add(self, other: Self) -> Self;
        fn sub(self, other: Self) -> Self;
        fn mul(self, other: Self) -> Self;
        fn div(self, other: Self) -> Self;
        fn abs(self) -> Float256;
        fn arg(self) -> Float256;
        fn conj(self) -> Self;
        fn magnitude(self) -> Float256;
        fn powc(self, exp: Self) -> Self;
        fn powf(self, exp: Float256) -> Self;
        fn powi(self, exp: i32) -> Self;
        fn exp(self) -> Self;
        fn ln(self) -> Self;
        fn log(self) -> Self;
        fn log10(self) -> Self;
        fn recip(self) -> Self;
        fn sqrt(self) -> Self;
        fn sin(self) -> Self;
        fn cos(self) -> Self;
        fn tan(self) -> Self;

    }
    impl  ComplexOps for Complex256{
        fn add(self, other: Self) -> Self {
            Complex256 {
                re: self.re + other.re,
                im: self.im + other.im,
            }
        }
    
        fn sub(self, other: Self) -> Self {
            Complex256 {
                re: self.re - other.re,
                im: self.im - other.im,
            }
        }
    
        fn mul(self, other: Self) -> Self {
            Complex256 {
                re: self.re * other.re - self.im * other.im,
                im: self.re * other.im + self.im * other.re,
            }
        }
    
        fn div(self, other: Self) -> Self {
            let denom = other.re * other.re + other.im * other.im;
            Complex256 {
                re: (self.re * other.re + self.im * other.im) / denom,
                im: (self.im * other.re - self.re * other.im) / denom,
            }
        }
    
        fn abs(self) -> Float256 {
            (self.re * self.re + self.im * self.im).sqrt()
        }
    
        fn arg(self) -> Float256 {
            self.im.atan2(&self.re)
        }
    
        fn conj(self) -> Self {
            Complex256 {
                re: self.re,
                im: -self.im,
            }
        }
    
        fn magnitude(self) -> Float256 {
            (self.re * self.re + self.im * self.im).sqrt()
        }
    
        fn powc(self, exp: Self) -> Self {
            let ln_self = self.ln();
            let prod = ln_self.mul(exp);
            prod.exp()
        }
    
        fn powf(self, exp: Float256) -> Self {
            let exp_complex = Complex256 { re: exp, im: Float256::from(0.0) };
            self.powc(exp_complex)
        }
    
        fn powi(self, exp: i32) -> Self {
            let mut result = Complex256 { re: Float256::from(1.0), im: Float256::from(0.0) };
            let mut base = self;
            let mut n = exp;
            if n < 0 {
                base = self.div(Complex256 { re: Float256::from(1.0), im: Float256::from(0.0) });
                n = -n;
            }
            while n > 0 {
                if n % 2 == 1 {
                    result = result.mul(base);
                }
                base = base.mul(base);
                n /= 2;
            }
            result
        }
    
        fn exp(self) -> Self {
            let exp_re = self.re.exp();
            Complex256 {
                re: exp_re * self.im.cos(),
                im: exp_re * self.im.sin(),
            }
        }
    
        fn ln(self) -> Self {
            Complex256 {
                re: self.magnitude().ln(),
                im: self.arg(),
            }
        }
    
        fn log(self) -> Self {
            self.ln().div(Complex256 { re: Float256::from(10.0_f64.ln()), im: Float256::from(0.0) })
        }
    
        fn log10(self) -> Self {
            self.ln().div(Complex256 { re: Float256::from(10.0_f64.ln()), im: Float256::from(0.0) })
        }

        fn recip(self) -> Self {
            let denom = self.re * self.re + self.im * self.im;
            Complex256 {
                re: self.re / denom,
                im: -self.im / denom,
            }
        }
    
        fn sqrt(self) -> Self {
            let mag = self.magnitude();
            let re_sqrt = ((self.re + mag) / Float256::from(2.0)).sqrt();
            let im_sqrt = ((mag - self.re) / Float256::from(2.0)).sqrt();
            Complex256 {
                re: re_sqrt,
                im: if self.im >= Float256::from(0.0) { im_sqrt } else { -im_sqrt },
            }
        }
    
        fn sin(self) -> Self {
            let e_i = self.mul(Complex256 { re: Float256::from(0.0), im: Float256::from(1.0) }).exp();
            let e_neg_i = self.mul(Complex256 { re: Float256::from(0.0), im: Float256::from(-1.0) }).exp();
            e_i.sub(e_neg_i).div(Complex256 { re: Float256::from(0.0), im: Float256::from(2.0) })
        }
    
        fn cos(self) -> Self {
            let e_i = self.mul(Complex256 { re: Float256::from(0.0), im: Float256::from(1.0) }).exp();
            let e_neg_i = self.mul(Complex256 { re: Float256::from(0.0), im: Float256::from(-1.0) }).exp();
            e_i.add(e_neg_i).div(Complex256 { re: Float256::from(2.0), im: Float256::from(0.0) })
        }
    
        fn tan(self) -> Self {
            self.sin().div(self.cos())
        }
    }

    impl std::fmt::Display for Complex256 {
        fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            write!(f, "{} + {}i", self.re, self.im)
        }
    }

    impl std::str::FromStr for Complex256 {
        type Err = String;
        fn from_str(s: &str) -> Result<Self, Self::Err> {
            let s = s.trim();
            if !s.ends_with('i') {
                return Err("Invalid complex number format".to_string());
            }
            let s = &s[..s.len() - 1];
            let parts: Vec<&str> = s.split('+').collect();
            if parts.len() != 2 {
                return Err("Invalid complex number format".to_string());
            }
            let re = parts[0].trim().parse::<Float256>().map_err(|_| "Invalid real part".to_string())?;
            let im = parts[1].trim().parse::<Float256>().map_err(|_| "Invalid imaginary part".to_string())?;
            Ok(Complex256 { re, im })
        }
    }
}