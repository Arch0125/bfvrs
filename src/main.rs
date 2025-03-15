use bfv::bfv::BFV;
use bfv::helper::{modinv, param_gen};
use bfv::ntt::precompute_ntt_tables;
use num_bigint::{BigInt, ToBigInt};
use num_traits::{One, ToPrimitive, Zero};
use rand::Rng;

fn main() {
    let pd = 1;

    let (t, n, q, psi, psiv, w, wv) = if pd == 0 {
        let t = BigInt::from(16);
        let n = 1024;
        let q = BigInt::from(132120577);
        let psi = BigInt::from(73993);

        let psiv = modinv(psi.clone(), q.clone()).unwrap();
        let w = psi.modpow(&BigInt::from(2), &q);
        let wv = modinv(w.clone(), q.clone()).unwrap();

        (t, n, q, psi, psiv, w, wv)
    } else {
        let t = BigInt::from(16);
        let n = 1024;
        let logq = 30;
        let (q_u, psi, psiv, w, wv) = param_gen(n, logq).unwrap();
        let q = q_u.to_bigint().unwrap();

        (t, n, q, psi, psiv, w, wv)
    };

    let mu = 0.0;

    let sigma_values = vec![0.1, 1.0, 10.0, 10000.0];

    let t_relin = BigInt::from(256);
    let p = q.pow(3) + BigInt::one();

    let qnp = precompute_ntt_tables(n, &q, &psi, &psiv, &w, &wv);

    println!("--- Starting BFV Demo (Multiple Cycles) ---");
    println!("Parameters:");
    println!("* n: {}", n);
    println!("* q: {}", q);
    println!("* t: {}", t);
    println!("* psi: {}", psi);
    println!("* psiv: {}", psiv);
    println!("* w: {}", w);
    println!("* wv: {}", wv);
    println!("* mu: {}", mu);
    println!("* T (relin): {}", t_relin);
    println!("* p: {}", p);
    println!();

    let cycles = 100;

    for sigma in sigma_values {
        let mut successes = 0;
        for _ in 0..cycles {
            let mut evaluator = BFV::new(n, q.clone(), t.clone(), mu, sigma, qnp.clone());
            evaluator.secret_key_gen();
            evaluator.public_key_gen();
            evaluator.eval_key_gen_v1(t_relin.clone());

            let n1 = 6;
            let n2 = 5;
            let expected_sum = n1 + n2;

            let m1 = evaluator.int_encode(n1);
            let m2 = evaluator.int_encode(n2);

            let ct1 = evaluator.encrypt(&m1);
            let ct2 = evaluator.encrypt(&m2);

            let ct_add = evaluator.homomorphic_addition(&ct1, &ct2);
            let mt_add = evaluator.decrypt(&ct_add);
            let decoded = evaluator.int_decode(&mt_add);

            if decoded == BigInt::from(expected_sum) {
                successes += 1;
            }
        }
        let accuracy = (successes as f64 / cycles as f64) * 100.0;
        println!(
            "Sigma: {:.2} -> Accuracy: {}/{} = {:.2}%",
            sigma, successes, cycles, accuracy
        );
    }
}
