use num_bigint::BigInt;
use num_bigint::RandBigInt;
use num_traits::FromPrimitive;
use num_traits::ToPrimitive;
use num_traits::Signed;

use crate::ntt::intt;
use crate::ntt::ntt;

pub struct Poly {
    pub n: u32,
    pub q: BigInt,
    pub np: Vec<Vec<BigInt>>,
    pub f: Vec<BigInt>,
    pub in_ntt: bool,
}

impl Poly {
    pub fn new(n: u32, q: BigInt, np: Vec<Vec<BigInt>>) -> Poly {
        Poly {
            n: n,
            q: q,
            np: np,
            f: vec![BigInt::from(0); n as usize],
            in_ntt: false,
        }
    }

    pub fn randomize(&mut self, b: u32, domain: bool, type_p: u32, mu: u32, sigma: u32) {
        if type_p == 0 {
            self.f = (0..self.n).map(|_| {
                (rand::thread_rng().gen_bigint_range(&BigInt::from(-(b as i32) / 2), &BigInt::from(b as i32 / 2)) % self.q.clone()).abs()
            }).collect();
            self.in_ntt = domain;
        } else {
            self.f = (0..self.n).map(|_| {
                (rand::thread_rng().gen_bigint_range(&BigInt::from(mu) , &BigInt::from(sigma)) % self.q.clone()).abs()
            }).collect();
            self.in_ntt = domain;
        }
    }

    pub fn add(&self, b: &Poly) -> Poly {
        if self.in_ntt != b.in_ntt {
            panic!("Polynomial Addiditon: Inputs must be in the same domain.");
        } else if self.q != b.q {
            panic!("Polynomial Addiditon: Inputs must have the same modulus");
        } else {
            let mut c = Poly::new(self.n, self.q.clone(), self.np.clone());
            c.f = self.f.iter().zip(b.f.iter()).map(|(x, y)| (x.clone() + y.clone()) % self.q.clone()).collect();
            c.in_ntt = self.in_ntt;
            c
        }
    }

    pub fn sub(&self, b: &Poly) -> Poly {
        if self.in_ntt != b.in_ntt {
            panic!("Polynomial Subtraction: Inputs must be in the same domain.");
        } else if self.q != b.q {
            panic!("Polynomial Subtraction: Inputs must have the same modulus");
        } else {
            let mut c = Poly::new(self.n, self.q.clone(), self.np.clone());
            c.f = self.f.iter().zip(b.f.iter()).map(|(x, y)| (x.clone() - y.clone()) % self
                .q.clone()).collect();
            c.in_ntt = self.in_ntt;
            c
        }
    }

    pub fn mul(&self, b: &Poly) -> Poly {
        if self.in_ntt != b.in_ntt {
            panic!("Polynomial Multiplication: Inputs must be in the same domain.");
        } else if self.q != b.q {
            panic!("Polynomial Multiplication: Inputs must have the same modulus");
        } else {
            let mut c = Poly::new(self.n, self.q.clone(), self.np.clone());
            if self.in_ntt == true && b.in_ntt == true {
                c.f = self.f.iter().zip(b.f.iter()).map(|(x, y)| (x.clone() * y.clone()) % self.q.clone()).collect();
                c.in_ntt = true;
            } else {
                let w_table = self.np[0].clone();
                let wv_table = self.np[1].clone();
                let psi_table = self.np[2].clone();
                let psiv_table = self.np[3].clone();

                let s_p = self.f.iter().enumerate().map(|(pwr, x)| (x.clone() * psi_table[pwr].clone()) % self.q.clone()).collect();
                let b_p = b.f.iter().enumerate().map(|(pwr, x)| (x.clone() * psi_table[pwr].clone()) % self.q.clone()).collect();
                let s_n = ntt(&s_p, &w_table, &self.q);
                let b_n = ntt(&b_p, &w_table, &self.q);
                let sb_n = s_n.iter().zip(b_n.iter()).map(|(x, y)| (x.clone() * y.clone()) % self.q.clone()).collect();
                let sb_p = intt(&sb_n, &wv_table, &self.q);
                c.f = sb_p.iter().enumerate().map(|(pwr, x)| (x.clone() * psiv_table[pwr].clone()) % self.q.clone()).collect();
                c.in_ntt = false;
            }
            c
        }
    }

    pub fn modulo(&self, base: u32) -> Poly {
        let mut b = Poly::new(self.n, self.q.clone(), self.np.clone());
        b.f = self.f.iter().map(|x| x % base).collect();
        b.in_ntt = self.in_ntt;
        b
    }

    pub fn round(&self) -> Poly {
        let mut b = Poly::new(self.n, self.q.clone(), self.np.clone());
        b.f = self.f.iter().map(|x| x.clone()).collect();
        b.in_ntt = self.in_ntt;
        b
    }

    pub fn equal(&self, b: &Poly) -> bool {
        if self.n != b.n {
            return false;
        } else if self.q != b.q {
            return false;
        } else {
            for (i, j) in self.f.iter().zip(b.f.iter()) {
                if i != j {
                    return false;
                }
            }
            return true;
        }
    }

    pub fn neg(&self) -> Poly {
        let mut b = Poly::new(self.n, self.q.clone(), self.np.clone());
        b.f = self.f.iter().map(|x| (-x.clone() % self.q.clone())).collect();
        b.in_ntt = self.in_ntt;
        b
    }


    pub fn to_ntt(&self) -> Poly {
        let mut b = Poly::new(self.n, self.q.clone(), self.np.clone());
        if self.in_ntt == false {
            b.f = ntt(&self.f, &self.np[0], &self.q);
            b.in_ntt = true;
        } else {
            b.f = self.f.iter().map(|x| x.clone()).collect();
            b.in_ntt = true;
        }
        b
    }

    pub fn to_pol(&self) -> Poly {
        let mut b = Poly::new(self.n, self.q.clone(), self.np.clone());
        if self.in_ntt == false {
            b.f = self.f.iter().map(|x| x.clone()).collect();
            b.in_ntt = false;
        } else {
            b.f = intt(&self.f, &self.np[1], &self.q);
            b.in_ntt = false;
        }
        b
    }
}

