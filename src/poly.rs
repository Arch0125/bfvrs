use crate::ntt::{intt, ntt};
use num_bigint::BigInt;
use num_traits::{One, ToPrimitive, Zero};
use rand::distributions::{Distribution, Uniform};
use rand::Rng;
use std::ops::Add;
use std::ops::Rem;

#[derive(Clone, Debug)]
pub struct Poly {
    pub n: usize,
    pub q: BigInt,
    pub np: [Vec<BigInt>; 4],
    pub f: Vec<BigInt>,
    pub in_ntt: bool,
}

impl Rem<BigInt> for Poly {
    type Output = Poly;

    fn rem(self, modulus: BigInt) -> Self::Output {
        let result_f = self.f.iter().map(|x| x % &modulus).collect();

        Poly {
            n: self.n,
            q: self.q.clone(),
            np: self.np.clone(),
            f: result_f,
            in_ntt: self.in_ntt,
        }
    }
}

impl Add for Poly {
    type Output = Result<Poly, String>;

    fn add(self, other: Poly) -> Self::Output {
        if self.n != other.n || self.q != other.q {
            return Err(
                "Polynomial Addition: Inputs must have the same degree and modulus.".to_string(),
            );
        }

        let result_f = self
            .f
            .iter()
            .zip(&other.f)
            .map(|(x, y)| (x + y) % &self.q)
            .collect();

        Ok(Poly {
            n: self.n,
            q: self.q.clone(),
            np: self.np.clone(),
            f: result_f,
            in_ntt: self.in_ntt,
        })
    }
}

impl Poly {
    pub fn new(n: usize, q: BigInt, np: [Vec<BigInt>; 4]) -> Self {
        Poly {
            n,
            q: q.clone(),
            np,
            f: vec![BigInt::zero(); n],
            in_ntt: false,
        }
    }

    pub fn randomize(&mut self, bound: &BigInt, domain: bool, dist_type: i32, mu: f64, sigma: f64) {
        if bound.is_zero() {
            self.f = vec![BigInt::zero(); self.n];
            self.in_ntt = domain;
            return;
        }

        let mut rng = rand::thread_rng();
        match dist_type {
            0 => {
                let half_bound = (bound / BigInt::from(2)).to_i64().unwrap_or(1);
                let range = Uniform::new(-half_bound, half_bound);
                self.f = (0..self.n)
                    .map(|_| BigInt::from(rng.sample(range)) % bound)
                    .collect();
            }
            1 => {
                if sigma == 0.0 {
                    self.f = vec![BigInt::zero(); self.n];
                } else {
                    self.f = (0..self.n)
                        .map(|_| {
                            let sample = rng.sample(rand_distr::Normal::new(mu, sigma).unwrap());

                            BigInt::from(sample.round() as i64) % bound
                        })
                        .collect();
                }
            }
            _ => panic!("Invalid distribution type"),
        }
        self.in_ntt = domain;
    }

    pub fn add(&self, other: &Poly) -> Result<Poly, String> {
        if self.in_ntt != other.in_ntt {
            return Err("Polynomial Addition: Inputs must be in the same domain.".to_string());
        }
        if self.q != other.q {
            return Err("Polynomial Addition: Inputs must have the same modulus.".to_string());
        }

        let result_f = self
            .f
            .iter()
            .zip(&other.f)
            .map(|(x, y)| (x + y) % &self.q)
            .collect();

        Ok(Poly {
            n: self.n,
            q: self.q.clone(),
            np: self.np.clone(),
            f: result_f,
            in_ntt: self.in_ntt,
        })
    }

    pub fn sub(&self, other: &Poly) -> Result<Poly, String> {
        if self.in_ntt != other.in_ntt {
            return Err("Polynomial Subtraction: Inputs must be in the same domain.".to_string());
        }
        if self.q != other.q {
            return Err("Polynomial Subtraction: Inputs must have the same modulus.".to_string());
        }

        let result_f = self
            .f
            .iter()
            .zip(&other.f)
            .map(|(x, y)| (x - y) % &self.q)
            .collect();

        Ok(Poly {
            n: self.n,
            q: self.q.clone(),
            np: self.np.clone(),
            f: result_f,
            in_ntt: self.in_ntt,
        })
    }

    pub fn mul(&self, other: &Poly) -> Result<Poly, String> {
        if self.in_ntt != other.in_ntt {
            return Err(
                "Polynomial Multiplication: Inputs must be in the same domain.".to_string(),
            );
        }
        if self.q != other.q {
            return Err(
                "Polynomial Multiplication: Inputs must have the same modulus.".to_string(),
            );
        }

        let result_f = if self.in_ntt {
            self.f
                .iter()
                .zip(&other.f)
                .map(|(x, y)| (x * y) % &self.q)
                .collect()
        } else {
            let w_table = &self.np[0];
            let wv_table = &self.np[1];
            let psi_table = &self.np[2];
            let psiv_table = &self.np[3];

            let s_p = self
                .f
                .iter()
                .enumerate()
                .map(|(pwr, x)| (x * &psi_table[pwr]) % &self.q)
                .collect();

            let b_p = other
                .f
                .iter()
                .enumerate()
                .map(|(pwr, x)| (x * &psi_table[pwr]) % &self.q)
                .collect();

            let s_n = ntt(s_p, w_table.clone(), &self.q);
            let b_n = ntt(b_p, w_table.clone(), &self.q);

            let sb_n = s_n
                .iter()
                .zip(&b_n)
                .map(|(x, y)| (x * y) % &self.q)
                .collect();

            let sb_p = intt(sb_n, wv_table.clone(), &self.q);

            sb_p.iter()
                .enumerate()
                .map(|(pwr, x)| (x * &psiv_table[pwr]) % &self.q)
                .collect()
        };

        Ok(Poly {
            n: self.n,
            q: self.q.clone(),
            np: self.np.clone(),
            f: result_f,
            in_ntt: self.in_ntt,
        })
    }

    pub fn coeff_mul(&self, other: &Poly) -> Result<Poly, String> {
        if self.n != other.n || self.q != other.q {
            return Err("Polynomials must have the same degree and modulus for coefficient-wise multiplication.".to_string());
        }

        let result_f = self
            .f
            .iter()
            .zip(&other.f)
            .map(|(x, y)| (x * y) % &self.q)
            .collect();

        Ok(Poly {
            n: self.n,
            q: self.q.clone(),
            np: self.np.clone(),
            f: result_f,
            in_ntt: false,
        })
    }

    pub fn modulo(&self, modulus: &BigInt) -> Poly {
        let result_f = self
            .f
            .iter()
            .map(|x| {
                let mut result = x % modulus;
                if result < BigInt::zero() {
                    result += modulus;
                }
                result
            })
            .collect();

        Poly {
            n: self.n,
            q: self.q.clone(),
            np: self.np.clone(),
            f: result_f,
            in_ntt: false,
        }
    }

    pub fn round(&self) -> Poly {
        let result_f = self
            .f
            .iter()
            .map(|x| {
                let x_f64 = x.to_f64().unwrap();
                BigInt::from(x_f64.round() as i64)
            })
            .collect();

        Poly {
            n: self.n,
            q: self.q.clone(),
            np: self.np.clone(),
            f: result_f,
            in_ntt: false,
        }
    }

    pub fn eq(&self, other: &Poly) -> bool {
        if self.n != other.n || self.q != other.q {
            return false;
        }
        self.f.iter().zip(&other.f).all(|(x, y)| x == y)
    }

    pub fn neg(&self) -> Poly {
        let result_f = self.f.iter().map(|x| (-x) % &self.q).collect();

        Poly {
            n: self.n,
            q: self.q.clone(),
            np: self.np.clone(),
            f: result_f,
            in_ntt: self.in_ntt,
        }
    }

    pub fn to_ntt(&self) -> Poly {
        if self.in_ntt {
            return self.clone();
        }

        let result_f = ntt(self.f.clone(), self.np[0].clone(), &self.q);

        Poly {
            n: self.n,
            q: self.q.clone(),
            np: self.np.clone(),
            f: result_f,
            in_ntt: true,
        }
    }

    pub fn to_pol(&self) -> Poly {
        if !self.in_ntt {
            return self.clone();
        }

        let result_f = intt(self.f.clone(), self.np[1].clone(), &self.q);

        Poly {
            n: self.n,
            q: self.q.clone(),
            np: self.np.clone(),
            f: result_f,
            in_ntt: false,
        }
    }
}

impl std::fmt::Display for Poly {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let tmp = std::cmp::min(self.n, 8);
        let mut pstr = format!("{}", self.f[0]);

        for i in 1..tmp {
            pstr = format!("{} + {}*x^{}", pstr, self.f[i], i);
        }

        if self.n > 8 {
            pstr = format!("{} + ...", pstr);
        }

        write!(f, "{}", pstr)
    }
}
