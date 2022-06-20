#![feature(int_log)]
mod ring;
mod signature;

use k256::elliptic_curve::Field;
use k256::{AffinePoint, Scalar};
use rand_core::OsRng;
use wasm_bindgen::prelude::*;

type FrontendRing = Vec<String>;

#[wasm_bindgen]
pub fn sign(index: u32, ring: JsValue) -> Result<JsValue, JsValue> {
    let pedersen_h = (AffinePoint::GENERATOR * Scalar::random(OsRng)).to_affine();

    let frontend_ring: FrontendRing = ring.into_serde().map_err(|e| e.to_string())?;
    let mut ring = parse_ring(frontend_ring)?;

    Ok(JsValue::from("Hello"))
}

fn parse_ring(frontend_ring: FrontendRing) -> Result<ring::Ring, String> {
    // TODO
    Ok(vec![AffinePoint::GENERATOR; 3])
}
