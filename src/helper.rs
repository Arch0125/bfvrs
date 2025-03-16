use num_bigint::{BigInt, BigUint};
use num_rational::BigRational;
use num_traits::{Euclid, One, ToPrimitive, Zero};
use rand::Rng;
use std::{cmp::Ordering, ops::{Div, Mul, Rem, Sub}};
use num_integer::Integer;

use crate::generate_prime::is_prime;

pub fn modinv(a: BigInt, m: BigInt) -> Result<BigInt, String> {
    let (g, x, _) = extended_gcd(a.clone(), m.clone());

    if !g.is_one() {
        return Err("Modular inverse does not exist".to_string());
    } else {
        Ok(x.rem_euclid(&m))
    }
}

pub fn canonical_mod(x: &BigInt, modulus: &BigInt) -> BigInt {
    let mut r = x % modulus;
    if r < BigInt::zero() {
        r += modulus;
    }
    r
}

pub fn div_round(numer: &BigInt, denom: &BigInt) -> BigInt {
    let (d, r) = numer.div_rem(denom);
    if &r * 2 >= *denom {
        d + BigInt::one()
    } else {
        d
    }
}


fn extended_gcd(a: BigInt, b: BigInt) -> (BigInt, BigInt, BigInt) {
    if a.is_zero() {
        (b, BigInt::zero(), BigInt::one())
    } else {
        let (g, x, y) = extended_gcd(b.clone() % a.clone(), a.clone());
        (g, y - (b / a.clone()) * x.clone(), x)
    }
}

pub fn gcd(n1: BigInt, n2: BigInt) -> BigInt {
    let mut a = n1;
    let mut b = n2;
    while b != BigInt::zero() {
        let temp = b.clone();
        b = a % b;
        a = temp;
    }
    a
}

pub fn rational_round(r: BigRational) -> BigInt {
    // Compute quotient and remainder of numer/denom.
    let (quotient, remainder) = r.numer().div_rem(r.denom());
    // Multiply remainder by 2.
    let double_rem: BigInt = &remainder * 2;
    match double_rem.cmp(r.denom()) {
        Ordering::Greater => quotient + BigInt::one(),
        Ordering::Less => quotient,
        Ordering::Equal => {
            // Exactly half: round to even.
            if quotient.is_even() {
                quotient
            } else {
                quotient + BigInt::one()
            }
        }
    }
}

pub fn int_reverse(a: usize, n: usize) -> usize {
    let b = format!("{:0width$b}", a, width = n);
    usize::from_str_radix(&b.chars().rev().collect::<String>(), 2).unwrap()
}

pub fn index_reverse(a: Vec<BigInt>, r: usize) -> Vec<BigInt> {
    let n = a.len();
    let mut b = vec![BigInt::zero(); n];
    for i in 0..n {
        let rev_idx = int_reverse(i, r);
        b[rev_idx] = a[i].clone();
    }
    b
}

pub fn ref_pol_mul(a: Vec<BigInt>, b: Vec<BigInt>, modulus: BigInt) -> Vec<BigInt> {
    let n = a.len();
    let mut c = vec![BigInt::zero(); 2 * n];
    let mut d = vec![BigInt::zero(); n];

    for (index_a, elem_a) in a.iter().enumerate() {
        for (index_b, elem_b) in b.iter().enumerate() {
            c[index_a + index_b] =
                (c[index_a + index_b].clone() + elem_a.clone() * elem_b.clone()) % modulus.clone();
        }
    }

    for i in 0..n {
        d[i] = (c[i].clone() - c[i + n].clone()) % modulus.clone();
    }
    d
}

pub fn is_root_of_unity(w: &BigInt, m: &BigInt, q: &BigInt) -> bool {
    if w == &BigInt::zero() {
        false
    } else {
        let exponent = m.clone() / BigInt::from(2);
        let result = w.modpow(&exponent, q);
        result == q.clone() - BigInt::one()
    }
}

pub fn get_proper_prime(n: usize, log_q: usize) -> Result<BigUint, String> {
    let factor = BigUint::from(2 * n);
    let mut value = (BigUint::one() << log_q) - factor.clone() + BigUint::one();
    let lbound = BigUint::one() << (log_q - 1);

    while value > lbound {
        if is_prime(&value, 11) {
            return Ok(value);
        } else {
            value = value - factor.clone();
        }
    }
    Err("Failed to find a proper prime.".to_string())
}

pub fn find_primitive_root(m: &BigInt, q: &BigInt) -> Result<BigInt, String> {
    let g = (q.clone() - BigInt::one()) / m.clone();

    if (q.clone() - BigInt::one()) != g.clone() * m.clone() {
        return Err("Invalid parameters for primitive root.".to_string());
    }

    let mut attempt_ctr = 0;
    let attempt_max = 100;

    while attempt_ctr < attempt_max {
        let a = BigInt::from(rand::thread_rng().gen_range(2..q.to_u64().unwrap()));
        let b = a.modpow(&g, q);

        if is_root_of_unity(&b, m, q) {
            return Ok(b);
        } else {
            attempt_ctr += 1;
        }
    }
    Err("Failed to find a primitive root.".to_string())
}

pub fn param_gen(
    n: usize,
    log_q: usize,
) -> Result<(BigUint, BigInt, BigInt, BigInt, BigInt), String> {
    let mut pfound = false;
    let mut q = BigUint::zero();
    let mut psi = BigInt::zero();

    while !pfound {
        q = get_proper_prime(n, log_q)?;
        match find_primitive_root(&BigInt::from(2 * n), &BigInt::from(q.clone())) {
            Ok(root) => {
                psi = root;
                pfound = true;
            }
            Err(_) => continue,
        }
    }

    let psi_inv = modinv(psi.clone(), BigInt::from(q.clone()))?;
    let w = psi.modpow(&BigInt::from(2), &BigInt::from(q.clone()));
    let w_inv = modinv(w.clone(), BigInt::from(q.clone()))?;

    Ok((q, psi, psi_inv, w, w_inv))
}
