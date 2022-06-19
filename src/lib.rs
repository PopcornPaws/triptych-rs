use k256::{AffinePoint, Scalar};
use rand_core::OsRng;
use wasm_bindgen::prelude::*;

type FrontendRing = Vec<String>;
type Ring = Vec<AffinePoint>;

// NOTE all these for loops should be lumped together to generate

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct VecElem {
    i_0: Scalar,
    i_1: Scalar,
}

// NOTE n = 2, i.e. N = 2^m
// N - number of elements in the ring (should be padded to 2^m)
#[wasm_bindgen]
pub fn sign(index: u32, ring: JsValue) -> Result<JsValue, JsValue> {
    let frontend_ring: FrontendRing = ring.into_serde().map_err(|e| e.to_string())?;
    let mut ring = parse_ring(frontend_ring)?;
    let m = pad_ring_to_2m(&mut ring)?;

    let a_vec = (0..m)
        .map(|_| {
            let i_0 = Scalar::generate_biased(OsRng);
            VecElem { i_0, i_1: -i_0 }
        })
        .collect::<Vec<VecElem>>();

    let sigma_vec = deltas(index, m);

    let c_vec = a_vec
        .iter()
        .zip(sigma_vec)
        .map(|(a, s)| {
            let i_0 = a.i_0 * (Scalar::ONE - (s.i_0 + s.i_0));
            let i_1 = a.i_1 * (Scalar::ONE - (s.i_1 + s.i_1));
            VecElem { i_0, i_1 }
        })
        .collect::<Vec<VecElem>>();

    let d_vec = a_vec
        .iter()
        .map(|a| {
            let i_0 = -(a.i_0 * a.i_0);
            let i_1 = -(a.i_1 * a.i_1);
            VecElem { i_0, i_1 }
        })
        .collect::<Vec<VecElem>>();

    // TODO generate commitments A, B, C, D

    Ok(JsValue::from("Hello"))
}

// NOTE max exponent of 2 is 32, i.e. max ring len = 2^32
fn deltas(index: u32, m: u8) -> Vec<VecElem> {
    assert!(m <= 32);
    assert!(1 << (m as u32) > index);

    let sigma_vec = (0..m)
        .map(|j| {
            if index & 1 << j == 0 {
                VecElem {
                    i_0: Scalar::ONE,
                    i_1: Scalar::ZERO,
                }
            } else {
                VecElem {
                    i_0: Scalar::ZERO,
                    i_1: Scalar::ONE,
                }
            }
        })
        .collect::<Vec<VecElem>>();

    sigma_vec
}

fn pad_ring_to_2m(ring: &mut Vec<AffinePoint>) -> Result<u8, String> {
    // TODO
    Ok(11)
}

fn parse_ring(frontend_ring: FrontendRing) -> Result<Ring, String> {
    // TODO
    Ok(vec![AffinePoint::GENERATOR; 3])
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    #[rustfmt::skip]
    fn kronecker_delta() {
        let index = 1234_u32; //0b0000010011010010
        let m = 11; // N = 2^11
        let d = deltas(index, m);

        assert_eq!(d.len(), 11);
        assert_eq!(d[0], VecElem { i_0: Scalar::ONE, i_1: Scalar::ZERO }); // 0
        assert_eq!(d[1], VecElem { i_0: Scalar::ZERO, i_1: Scalar::ONE }); // 1
        assert_eq!(d[2], VecElem { i_0: Scalar::ONE, i_1: Scalar::ZERO }); // 0
        assert_eq!(d[3], VecElem { i_0: Scalar::ONE, i_1: Scalar::ZERO }); // 0
        assert_eq!(d[4], VecElem { i_0: Scalar::ZERO, i_1: Scalar::ONE }); // 1
        assert_eq!(d[5], VecElem { i_0: Scalar::ONE, i_1: Scalar::ZERO }); // 0
        assert_eq!(d[6], VecElem { i_0: Scalar::ZERO, i_1: Scalar::ONE }); // 1
        assert_eq!(d[7], VecElem { i_0: Scalar::ZERO, i_1: Scalar::ONE }); // 1
        assert_eq!(d[8], VecElem { i_0: Scalar::ONE, i_1: Scalar::ZERO }); // 0
        assert_eq!(d[9], VecElem { i_0: Scalar::ONE, i_1: Scalar::ZERO }); // 0
        assert_eq!(d[10], VecElem { i_0: Scalar::ZERO, i_1: Scalar::ONE }); // 1
    }
}
