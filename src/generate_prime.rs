use num_bigint::{BigUint, RandBigInt};
use num_traits::{One, Zero};
use rand::Rng;
use std::ops::BitAnd;

pub fn miller_rabin(p: &BigUint, s: u32) -> bool {
    let mut r = p - BigUint::one();
    let mut u = 0;

    while (&r).bitand(&BigUint::one()) == BigUint::zero() {
        u += 1;
        r = r >> 1;
    }

    for _ in 0..s {
        let a = rand::thread_rng().gen_range(BigUint::from(2u32)..(p - BigUint::one()));
        let mut z = mod_pow(&a, &r, p);

        if z != BigUint::one() && z != (p - BigUint::one()) {
            for _ in 0..u - 1 {
                if z != p - BigUint::one() {
                    z = mod_pow(&z, &BigUint::from(2u32), p);
                    if z == BigUint::one() {
                        return false;
                    }
                } else {
                    break;
                }
            }
            if z != p - BigUint::one() {
                return false;
            }
        }
    }
    true
}

pub fn mod_pow(base: &BigUint, exp: &BigUint, modulus: &BigUint) -> BigUint {
    if modulus == &BigUint::one() {
        return BigUint::zero();
    }
    let mut result = BigUint::one();
    let mut base = base % modulus;
    let mut exp = exp.clone();

    while exp > BigUint::zero() {
        if &exp & BigUint::one() == BigUint::one() {
            result = (result * &base) % modulus;
        }
        exp >>= 1;
        base = (&base * &base) % modulus;
    }
    result
}

pub fn is_prime(n: &BigUint, s: u32) -> bool {
    let low_primes: [u32; 167] = [
        3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89,
        97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181,
        191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281,
        283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397,
        401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503,
        509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619,
        631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743,
        751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863,
        877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997,
    ];

    if n >= &BigUint::from(3u32) {
        if n.bitand(BigUint::one()) != BigUint::zero() {
            for &p in &low_primes {
                let p_big = BigUint::from(p);
                if n == &p_big {
                    return true;
                }
                if n % p_big == BigUint::zero() {
                    return false;
                }
            }
            return miller_rabin(n, s);
        }
    }
    false
}

pub fn generate_large_prime(k: usize, s: u32) -> Result<BigUint, String> {
    let mut r = (100.0 * ((k as f64).log2() + 1.0)) as i32;
    while r > 0 {
        let mut rng = rand::thread_rng();
        let n = rng.gen_biguint(k.try_into().unwrap());
        r -= 1;
        if is_prime(&n, s) {
            return Ok(n);
        }
    }
    Err(format!("Failure after {} tries.", r))
}

pub fn main() {
    match generate_large_prime(1024, 11) {
        Ok(prime) => println!("Generated prime: {}", prime),
        Err(e) => println!("{}", e),
    }
}
