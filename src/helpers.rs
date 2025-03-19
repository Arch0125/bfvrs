use num_bigint::RandBigInt;
use num_traits::Num;
use rand::Rng;
use num_bigint::BigInt;
use num_traits::FromPrimitive;
use num_traits::ToPrimitive;
use crate::generate_prime::is_prime;

pub fn egcd(a: BigInt, b: BigInt) -> (BigInt, BigInt, BigInt) {
    if a == BigInt::from(0) {
        return (b, BigInt::from(0), BigInt::from(1));
    } else {
        let (g, y, x) = egcd(b.clone() % a.clone(), a.clone());
        return (g, x - (b.clone() / a.clone()) * y.clone(), y);
    }
}

pub fn modinv(a: BigInt, m: BigInt) -> BigInt {
    let (g, x, y) = egcd(a.clone(), m.clone());
    if g != BigInt::from(1) {
        panic!("Modular inverse does not exist");
    } else {
        return x % m;
    }
}

pub fn gcd(a: BigInt, b: BigInt) -> BigInt {
    let mut a = a;
    let mut b = b;
    while b != BigInt::from(0) {
        let temp = a.clone();
        a = b.clone();
        b = temp % b;
    }
    return a;
}

pub fn int_reverse(a: BigInt, n: u32) -> BigInt {
    let b = format!("{:0width$b}", a, width = n as usize);
    let c = b.chars().rev().collect::<String>();
    return BigInt::from_str_radix(&c, 2).unwrap();
}

pub fn ref_pol_mul(a: &Vec<BigInt>, b: &Vec<BigInt>, m: &BigInt) -> Vec<BigInt> {
    let mut c = vec![BigInt::from(0); 2 * a.len()];
    let mut d = vec![BigInt::from(0); a.len()];
    for (index_a, elem_a) in a.iter().enumerate() {
        for (index_b, elem_b) in b.iter().enumerate() {
            c[index_a + index_b] = (c[index_a + index_b].clone() + elem_a.clone() * elem_b.clone()) % m;
        }
    }
    for i in 0..a.len() {
        d[i] = (c[i].clone() - c[i + a.len()].clone()) % m;
    }
    return d;
}

pub fn ref_pol_mul_v2(a: &Vec<BigInt>, b: &Vec<BigInt>) -> Vec<BigInt> {
    let mut c = vec![BigInt::from(0); 2 * a.len()];
    let mut d = vec![BigInt::from(0); a.len()];
    for (index_a, elem_a) in a.iter().enumerate() {
        for (index_b, elem_b) in b.iter().enumerate() {
            c[index_a + index_b] = c[index_a + index_b].clone() + elem_a.clone() * elem_b.clone();
        }
    }
    for i in 0..a.len() {
        d[i] = c[i].clone() - c[i + a.len()].clone();
    }
    return d;
}

pub fn is_root_of_unity(w: BigInt, m: BigInt, q: BigInt) -> bool {
    if w == BigInt::from(0) {
        return false;
    } else if w.clone().modpow(&(m.clone() / BigInt::from(2)), &q) == q - BigInt::from(1) {
        return true;
    } else {
        return false;
    }
}

pub fn get_proper_prime(n: BigInt, logq: u32) -> BigInt {
    let factor = BigInt::from(2) * n.clone();
    let value = (BigInt::from(1) << logq) - factor.clone() + BigInt::from(1);
    let lbound = BigInt::from(1) << (logq - 1);
    let mut value = value;
    while value > lbound {
        if is_prime(value.clone(), BigInt::from(5)) {
            return value;
        } else {
            value = value - factor.clone();
        }
    }
    panic!("Failed to find a proper prime.");
}

pub fn find_primitive_root(m: BigInt, q: BigInt) -> (bool, BigInt) {
    let g = (q.clone() - BigInt::from(1)) / m.clone();
    if (q.clone() - BigInt::from(1)) != g.clone() * m.clone() {
        return (false, BigInt::from(0));
    }
    let mut attempt_ctr = 0;
    let attempt_max = 100;
    while attempt_ctr < attempt_max {
        let a = rand::thread_rng().gen_bigint_range(&BigInt::from(2), &(&q.clone() - BigInt::from(1)));
        let b = a.clone().modpow(&g, &q);
        if is_root_of_unity(b.clone(), m.clone(), q.clone()) {
            return (true, b);
        } else {
            attempt_ctr = attempt_ctr + 1;
        }
    }
    return (true, BigInt::from(0));
}

pub fn index_reverse(a: &Vec<BigInt>, r: u32) -> Vec<BigInt> {
    let n = a.len();
    let mut b = vec![BigInt::from(0); n];
    for i in 0..n {
        let rev_idx = int_reverse(BigInt::from(i), r);
        b[rev_idx.to_usize().unwrap()] = a[i].clone();
    }
    return b;
}