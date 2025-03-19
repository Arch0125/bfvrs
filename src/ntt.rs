use num_bigint::BigInt;

use crate::helpers::index_reverse;

pub fn ntt(a: &Vec<BigInt>, w_table: &Vec<BigInt>, q: &BigInt) -> Vec<BigInt> {
    let n = a.len();
    let mut b = vec![BigInt::from(0); n];
    let v = (n as f64).log2() as u32;

    for i in 0..v {
        for j in 0..(1 << i) {
            for k in 0..(1 << (v - i - 1)) {
                let s = j * (1 << (v - i)) + k;
                let t = s + (1 << (v - i - 1));

                let w = w_table[(1 << i) * k].clone();

                let as_temp = b[s].clone();
                let at_temp = b[t].clone();

                b[s] = (as_temp.clone() + at_temp.clone()) % q;
                b[t] = ((as_temp - at_temp) * w) % q;
            }
        }
    }

    index_reverse(&b, v)
}   

pub fn intt(a: &Vec<BigInt>, w_table: &Vec<BigInt>, q: &BigInt) -> Vec<BigInt> {
    let n = a.len();
    let mut b = vec![BigInt::from(0); n];
    let v = (n as f64).log2() as u32;

    for i in 0..v {
        for j in 0..(1 << i) {
            for k in 0..(1 << (v - i - 1)) {
                let s = j * (1 << (v - i)) + k;
                let t = s + (1 << (v - i - 1));

                let w = w_table[(1 << i) * k].clone();

                let as_temp = b[s].clone();
                let at_temp = b[t].clone();

                b[s] = (as_temp.clone() + at_temp.clone()) % q;
                b[t] = ((as_temp - at_temp) * w) % q;
            }
        }
    }

    b = index_reverse(&b, v);

    let n_inv = crate::helpers::modinv(BigInt::from(n), q.clone());
    for i in 0..n {
        b[i] = (b[i].clone() * n_inv.clone()) % q.clone();
    }

    b
}