use crate::ring::{pad_ring_to_2n, Ring};
use k256::elliptic_curve::Field;
use k256::{AffinePoint, ProjectivePoint, Scalar};
use rand_core::OsRng;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct VecElem {
    i_0: Scalar,
    i_1: Scalar,
}

pub struct Signature {
    tag: AffinePoint,
}

// NOTE n = 2, i.e. N = 2^m
// N - number of elements in the ring (should be padded to 2^m)
impl Signature {
    pub fn new(index: u32, mut ring: Ring, pedersen_h: AffinePoint) -> Result<Self, String> {
        let m = pad_ring_to_2n(&mut ring)?;

        let a_vec = (0..m)
            .map(|_| {
                let i_0 = Scalar::random(OsRng);
                VecElem { i_0, i_1: -i_0 }
            })
            .collect::<Vec<VecElem>>();

        // NOTE max exponent of 2 is 32, i.e. max ring len = 2^32
        assert!(m <= 32);
        assert!(1 << (m as u32) > index);
        let b_vec = deltas(index, m); // sigma vec
                                      //
        let c_vec = a_vec
            .iter()
            .zip(b_vec.iter())
            .map(|(a, b)| {
                let i_0 = a.i_0 * (Scalar::ONE - (b.i_0 + b.i_0));
                let i_1 = a.i_1 * (Scalar::ONE - (b.i_1 + b.i_1));
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
        let mut generators = Vec::<[AffinePoint; 2]>::with_capacity(m);
        let mut a_com_gen = ProjectivePoint::IDENTITY;
        let mut b_com_gen = ProjectivePoint::IDENTITY;
        let mut c_com_gen = ProjectivePoint::IDENTITY;
        let mut d_com_gen = ProjectivePoint::IDENTITY;
        for i in 0..(m as usize) {
            let generator_0 = (AffinePoint::GENERATOR * Scalar::random(OsRng)).to_affine();
            let generator_1 = (AffinePoint::GENERATOR * Scalar::random(OsRng)).to_affine();
            a_com_gen += generator_0 * a_vec[i].i_0;
            a_com_gen += generator_1 * a_vec[i].i_1;
            b_com_gen += generator_0 * b_vec[i].i_0;
            b_com_gen += generator_1 * b_vec[i].i_1;
            c_com_gen += generator_0 * c_vec[i].i_0;
            c_com_gen += generator_1 * c_vec[i].i_1;
            d_com_gen += generator_0 * d_vec[i].i_0;
            d_com_gen += generator_1 * d_vec[i].i_1;
            generators.push([generator_0, generator_1]);
        }

        let a_com = (pedersen_h * Scalar::random(OsRng) + a_com_gen).to_affine();
        let b_com = (pedersen_h * Scalar::random(OsRng) + b_com_gen).to_affine();
        let c_com = (pedersen_h * Scalar::random(OsRng) + c_com_gen).to_affine();
        let d_com = (pedersen_h * Scalar::random(OsRng) + d_com_gen).to_affine();

        for k in 0..ring.len() {}

        todo!();
    }
}

fn deltas(num: u32, n: usize) -> Vec<VecElem> {
    (0..n)
        .map(|j| {
            if num & 1 << j == 0 {
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
        .collect::<Vec<VecElem>>()
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
