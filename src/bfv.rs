use crate::helper::{canonical_mod, div_round};
use crate::ntt::{intt, ntt};
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, ToPrimitive, Zero};
use poly::Poly;
use rand::distributions::{Distribution, Uniform};
use rand::Rng;
use std::ops::Add;
use std::os::unix::raw::pthread_t;

use crate::poly;

pub struct BFV {
    n: usize,
    q: BigInt,
    t: BigInt,
    T: BigInt,
    l: usize,
    p: BigInt,
    mu: f64,
    sigma: f64,
    qnp: [Vec<BigInt>; 4],
    sk: Poly,
    pk: Vec<Poly>,
    rlk1: Vec<Vec<Poly>>,
    rlk2: Vec<Poly>,
}

impl BFV {
    pub fn new(n: usize, q: BigInt, t: BigInt, mu: f64, sigma: f64, qnp: [Vec<BigInt>; 4]) -> Self {
        BFV {
            n,
            q: q.clone(),
            t: t.clone(),
            T: BigInt::zero(),
            l: 0,
            p: BigInt::zero(),
            mu,
            sigma,
            qnp: qnp.clone(),
            sk: Poly::new(n, q.clone(), qnp.clone()),
            pk: Vec::new(),
            rlk1: Vec::new(),
            rlk2: Vec::new(),
        }
    }

    pub fn int_decode(&self, m: &Poly) -> BigInt {
        let mut mr = BigInt::zero();

        let k = std::cmp::min(8, m.f.len());
        let thr = if self.t == BigInt::from(2) {
            BigInt::from(2)
        } else {
            (&self.t + BigInt::one()) / BigInt::from(2)
        };

        for (i, c) in m.f.iter().take(k).enumerate() {
            let c_ = if c >= &thr { -(&self.t - c) } else { c.clone() };
            let two_pow_i = BigInt::from(2).pow(i as u32);
            mr += c_ * two_pow_i;
        }
        mr
    }

    pub fn eval_key_gen_v2(&mut self, p: BigInt) {
        self.p = p.clone();

        let modulus = &self.p * &self.q;
        let mut a = Poly::new(self.n, modulus.clone(), self.qnp.clone());
        let mut e = Poly::new(self.n, modulus.clone(), self.qnp.clone());

        a.randomize(
            &BigInt::from(modulus.to_i32().unwrap_or(i32::MAX)),
            false,
            0,
            0.0,
            0.0,
        );
        e.randomize(&BigInt::from(0), false, 1, self.mu, self.sigma);

        let a_sk_coeffs: Vec<BigInt> =
            a.f.iter()
                .zip(&self.sk.f)
                .map(|(a_coeff, sk_coeff)| (a_coeff * sk_coeff) % &modulus)
                .collect();

        let a_sk = Poly {
            n: self.n,
            q: modulus.clone(),
            np: self.qnp.clone(),
            f: a_sk_coeffs,
            in_ntt: false,
        };

        let neg_a_sk_e = (a_sk.neg() + e).unwrap();

        let sk2_coeffs: Vec<BigInt> = self
            .sk
            .f
            .iter()
            .zip(&self.sk.f)
            .map(|(sk1_coeff, sk2_coeff)| (sk1_coeff * sk2_coeff) % &modulus)
            .collect();

        let sk2 = Poly {
            n: self.n,
            q: modulus.clone(),
            np: self.qnp.clone(),
            f: sk2_coeffs,
            in_ntt: false,
        };

        let p_poly = Poly {
            n: self.n,
            q: modulus.clone(),
            np: self.qnp.clone(),
            f: vec![self.p.clone(); self.n],
            in_ntt: false,
        };

        let p_sk2 = sk2.mul(&p_poly).unwrap();

        let c0 = (neg_a_sk_e + p_sk2).unwrap();

        self.rlk2 = vec![c0, a];
    }

    pub fn int_encode(&self, m: i32) -> Poly {
        let mut mr = Poly::new(self.n, self.t.clone(), self.qnp.clone());
        if m > 0 {
            let mut mt = m;
            for i in 0..self.n {
                mr.f[i] = BigInt::from(mt % 2);
                mt /= 2;
            }
        } else if m < 0 {
            let mut mt = -m;
            for i in 0..self.n {
                mr.f[i] = (&self.t - BigInt::from(mt % 2)) % &self.t;
                mt /= 2;
            }
        }
        mr
    }

    pub fn decrypt_v2(&self, ct: &[Poly]) -> Poly {
        let sk2 = self.sk.mul(&self.sk).unwrap();
        let term1 = ct[0].clone();
        let term2 = ct[1].mul(&self.sk).unwrap();
        let term3 = ct[2].mul(&sk2).unwrap();
        let tmp = (term1 + term2).unwrap();
        let m = (tmp + term3).unwrap();
        let m_scaled = Poly {
            n: self.n,
            q: self.t.clone(),
            np: self.qnp.clone(),
            f: m.f
                .iter()
                .map(|x| ((&self.t * x) / &self.q) % &self.t)
                .collect(),
            in_ntt: false,
        };
        m_scaled.round()
    }

    pub fn relinearization_v1(&self, ct: &[Poly]) -> Vec<Poly> {
        let c0 = ct[0].clone();
        let c1 = ct[1].clone();
        let c2 = ct[2].clone();

        let mut c0r = c0.clone();
        let mut c1r = c1.clone();

        for i in 0..=self.l {
            let c2i = c2.clone() % BigInt::from(self.T.pow(i as u32));
            c0r = (c0r + self.rlk1[i][0].mul(&c2i).unwrap()).unwrap();
            c1r = (c1r + self.rlk1[i][1].mul(&c2i).unwrap()).unwrap();
        }

        vec![c0r, c1r]
    }

    pub fn secret_key_gen(&mut self) {
        let mut s = Poly::new(self.n, self.q.clone(), self.qnp.clone());
        s.randomize(&BigInt::from(2), false, 0, 0.0, 0.0);
        self.sk = s;
    }

    pub fn public_key_gen(&mut self) {
        let mut a = Poly::new(self.n, self.q.clone(), self.qnp.clone());

        a.randomize(&self.q, false, 0, 0.0, 0.0);

        let mut e = Poly::new(self.n, self.q.clone(), self.qnp.clone());

        if self.sigma > 0.0 {
            e.randomize(&self.q, false, 1, self.mu, self.sigma);
        } else {
            e.randomize(&BigInt::from(1), false, 1, self.mu, self.sigma);
        }

        let pk0 = (a.mul(&self.sk).unwrap().neg() + e).unwrap();
        let pk1 = a;

        self.pk = vec![pk0, pk1];
    }

    pub fn div_round(numer: &BigInt, denom: &BigInt) -> BigInt {
        let (d, r) = numer.div_rem(denom);

        if &r * 2 >= *denom {
            d + BigInt::one()
        } else {
            d
        }
    }

    pub fn eval_key_gen_v1(&mut self, T: BigInt) {
        self.T = T.clone();
        self.l = (self.q.to_f64().unwrap().log(T.to_f64().unwrap()).floor()) as usize;

        let sk2 = self.sk.mul(&self.sk).unwrap();

        for i in 0..=self.l {
            let mut ai = Poly::new(self.n, self.q.clone(), self.qnp.clone());
            let mut ei = Poly::new(self.n, self.q.clone(), self.qnp.clone());

            ai.randomize(&BigInt::from(self.q.to_i32().unwrap()), false, 0, 0.0, 0.0);
            ei.randomize(&BigInt::from(2), false, 1, self.mu, self.sigma);

            let mut ts2 = Poly::new(self.n, self.q.clone(), self.qnp.clone());
            ts2.f = sk2
                .f
                .iter()
                .map(|x| (T.pow(i as u32) * x) % &self.q)
                .collect();

            let rlki0 = (ts2.sub(&ai.mul(&self.sk).unwrap()).unwrap() + ei).unwrap();
            let rlki1 = ai;

            self.rlk1.push(vec![rlki0, rlki1]);
        }
    }

    pub fn encrypt(&self, m: &Poly) -> Vec<Poly> {
        let delta = &self.q / &self.t;

        let mut u = Poly::new(self.n, self.q.clone(), self.qnp.clone());
        u.randomize(&BigInt::from(2), false, 0, 0.0, 0.0);

        let mut e1 = Poly::new(self.n, self.q.clone(), self.qnp.clone());
        let mut e2 = Poly::new(self.n, self.q.clone(), self.qnp.clone());
        e1.randomize(&BigInt::from(2), false, 1, self.mu, self.sigma);
        e2.randomize(&BigInt::from(2), false, 1, self.mu, self.sigma);

        let md = Poly {
            n: self.n,
            q: self.q.clone(),
            np: self.qnp.clone(),
            f: m.f.iter().map(|x| (delta.clone() * x) % &self.q).collect(),
            in_ntt: false,
        };

        let c0 = ((self.pk[0].mul(&u).unwrap() + e1).unwrap() + md).unwrap();
        let c1 = (self.pk[1].mul(&u).unwrap() + e2).unwrap();

        vec![c0, c1]
    }

    pub fn decrypt(&self, ct: &[Poly]) -> Poly {
        
        let m = (ct[1].mul(&self.sk).unwrap() + ct[0].clone()).unwrap();
    
        
        let scaled_coeffs = m.f.iter()
            .map(|x| div_round(&(&self.t * x), &self.q))
            .collect();
        
        let m_scaled = Poly {
            n: self.n,
            q: self.t.clone(),  
            np: self.qnp.clone(),
            f: scaled_coeffs,
            in_ntt: false,
        };
        
        m_scaled.modulo(&self.t)
    }
    
    
    

    pub fn homomorphic_addition(&self, ct0: &[Poly], ct1: &[Poly]) -> Vec<Poly> {
        let ct0_b = (ct0[0].clone() + ct1[0].clone()).unwrap();
        let ct1_b = (ct0[1].clone() + ct1[1].clone()).unwrap();
        vec![ct0_b, ct1_b]
    }

    pub fn homomorphic_multiplication(&self, ct0: &[Poly], ct1: &[Poly]) -> Vec<Poly> {
        let r0 = ct0[0].mul(&ct1[0]).unwrap();
        let r1 = ct0[0].mul(&ct1[1]).unwrap();
        let r2 = ct0[1].mul(&ct1[0]).unwrap();
        let r3 = ct0[1].mul(&ct1[1]).unwrap();

        vec![r0, r1, r2, r3]
    }
}

impl std::fmt::Display for BFV {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "--- Parameters:")?;
        writeln!(f, "n    : {}", self.n)?;
        writeln!(f, "q    : {}", self.q)?;
        writeln!(f, "t    : {}", self.t)?;
        writeln!(f, "T    : {}", self.T)?;
        writeln!(f, "l    : {}", self.l)?;
        writeln!(f, "p    : {}", self.p)?;
        writeln!(f, "mu   : {}", self.mu)?;
        writeln!(f, "sigma: {}", self.sigma)
    }
}
