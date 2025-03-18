use num_bigint::{BigInt, RandBigInt};
use rand::Rng;
use num_traits::ToPrimitive;

pub fn miller_rabin(p: BigInt, s: BigInt) -> bool {
    let mut r = p.clone() - BigInt::from(1);
    let mut u = BigInt::from(0);
    while r.to_u128().unwrap() & 1 == 0 {
        u += 1;
        r = r / 2;
    }
    
    let p_value = p.clone();
    let minus_p = p_value - BigInt::from(1);
    for _ in 0..s.to_u128().unwrap() {
        let a = rand::thread_rng().gen_bigint_range(&BigInt::from(2), &BigInt::from(minus_p.clone()));
        let mut z = a.pow(r.to_u32().unwrap() ) % p.clone();
        if z != BigInt::from(1) && z != minus_p {
            for _ in 0..u.to_u128().unwrap() - 1 {
                if z != minus_p {
                    z = z.pow(2) % p.clone();
                    if z == BigInt::from(1) {
                        return false;
                    }
                } else {
                    break;
                }
            }
            if z !=  minus_p {
                return false;
            }
        }
    }
    true
}

pub fn is_prime(n: BigInt, s: BigInt) -> bool {
    let low_primes = vec![3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97
                    ,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179
                    ,181,191,193,197,199,211,223,227,229,233,239,241,251,257,263,269
                    ,271,277,281,283,293,307,311,313,317,331,337,347,349,353,359,367
                    ,373,379,383,389,397,401,409,419,421,431,433,439,443,449,457,461
                    ,463,467,479,487,491,499,503,509,521,523,541,547,557,563,569,571
                    ,577,587,593,599,601,607,613,617,619,631,641,643,647,653,659,661
                    ,673,677,683,691,701,709,719,727,733,739,743,751,757,761,769,773
                    ,787,797,809,811,821,823,827,829,839,853,857,859,863,877,881,883
                    ,887,907,911,919,929,937,941,947,953,967,971,977,983,991,997];
    if n >= BigInt::from(3) {
        if n.clone() & BigInt::from(1) != BigInt::from(0) {
            for p in low_primes {
                if n == BigInt::from(p) {
                    return true;
                }
                if n.clone() % BigInt::from(p) == BigInt::from(0) {
                    return false;
                }
            }
            return miller_rabin(n, s);
        }
    }
    false
}

pub fn generate_large_prime(k: BigInt, s: BigInt) -> BigInt {
    let mut r = 100 * (k.to_f64().unwrap()).log2().ceil() as u32;
    while r > 0 {
        let n = rand::thread_rng().gen_bigint_range(&BigInt::from(2).pow(k.to_u32().unwrap() - 1), &BigInt::from(2).pow(k.to_u32().unwrap()));
        r -= 1;
        if is_prime(n.clone(), s.clone()) {
            return n;
        }
    }
    panic!("Failure after {} tries.", r);
}