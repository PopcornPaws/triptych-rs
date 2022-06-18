use k256::Scalar;
use rand_core::OsRng;
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
pub fn generate_random_scalars(n: u32, m: u32) -> JsValue {
    // NOTE n = 2, i.e. N = 2^m
    // N - number of elements in the ring (should be padded to 2^m)
    let a_vec = (0..m)
        .map(|_| Scalar::generate_biased(OsRng))
        .collect::<Vec<Scalar>>();
    JsValue::from("Hello")
}
