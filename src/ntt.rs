use crate::helper::{index_reverse, modinv};
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::ToPrimitive;
use num_traits::{One, Zero};

pub fn ntt(a: Vec<BigInt>, w_table: Vec<BigInt>, q: &BigInt) -> Vec<BigInt> {
    let n = a.len();
    let mut b = a.clone();

    let v = (n as f64).log2() as usize;

    for i in 0..v {
        for j in 0..(1 << i) {
            for k in 0..(1 << (v - i - 1)) {
                let s = j * (1 << (v - i)) + k;
                let t = s + (1 << (v - i - 1));

                let w = &w_table[(1 << i) * k];

                let as_temp = b[s].clone();
                let at_temp = b[t].clone();

                b[s] = (&as_temp + &at_temp) % q;
                b[t] = ((&as_temp - &at_temp) * w) % q;
            }
        }
    }

    index_reverse(b, v)
}

pub fn intt(a: Vec<BigInt>, w_table: Vec<BigInt>, q: &BigInt) -> Vec<BigInt> {
    let n = a.len();
    let mut b = a.clone();

    let v = (n as f64).log2() as usize;

    for i in 0..v {
        for j in 0..(1 << i) {
            for k in 0..(1 << (v - i - 1)) {
                let s = j * (1 << (v - i)) + k;
                let t = s + (1 << (v - i - 1));

                let w = &w_table[(1 << i) * k];

                let as_temp = b[s].clone();
                let at_temp = b[t].clone();

                b[s] = (&as_temp + &at_temp) % q;
                b[t] = ((&as_temp - &at_temp) * w) % q;
            }
        }
    }

    let n_bigint = BigInt::from(n);

    let n_inv = modinv(n_bigint, q.clone()).expect("Failed to compute modular inverse");

    for i in 0..n {
        b[i] = (&b[i] * &n_inv) % q;
    }

    index_reverse(b, v)
}

pub fn precompute_ntt_tables(
    n: usize,
    q: &BigInt,
    psi: &BigInt,
    psi_inv: &BigInt,
    w: &BigInt,
    w_inv: &BigInt,
) -> [Vec<BigInt>; 4] {
    let mut w_table = vec![BigInt::one(); n];
    let mut wv_table = vec![BigInt::one(); n];
    let mut psi_table = vec![BigInt::one(); n];
    let mut psiv_table = vec![BigInt::one(); n];

    for i in 1..n {
        w_table[i] = (&w_table[i - 1] * w) % q;
        wv_table[i] = (&wv_table[i - 1] * w_inv) % q;
        psi_table[i] = (&psi_table[i - 1] * psi) % q;
        psiv_table[i] = (&psiv_table[i - 1] * psi_inv) % q;
    }

    [w_table, wv_table, psi_table, psiv_table]
}
