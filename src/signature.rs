use crate::ring::*;
use k256::elliptic_curve::group::GroupEncoding;
use k256::elliptic_curve::{Field, PrimeField};
use k256::{AffinePoint, ProjectivePoint, Scalar};
use rand_core::OsRng;
use sha3::{Digest, Keccak256};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct VecElem {
    i_0: Scalar,
    i_1: Scalar,
}

pub struct Signature {
    a_commitment: AffinePoint,
    b_commitment: AffinePoint,
    c_commitment: AffinePoint,
    d_commitment: AffinePoint,
    x_points: Vec<AffinePoint>,
    y_points: Vec<AffinePoint>,
    f_scalars: Vec<VecElem>,
    z_a_point: AffinePoint,
    z_c_point: AffinePoint,
    z_scalar: Scalar,
}

// NOTE n = 2, i.e. N = 2^m
// N - number of elements in the ring (should be padded to 2^m)
impl Signature {
    pub fn new(
        index: usize,
        mut ring: Ring,
        privkey: Scalar,
        h_point: AffinePoint,
        j_point: AffinePoint,
    ) -> Result<Self, String> {
        let mut hasher = Keccak256::new();
        let m = pad_ring_to_2n(&mut ring)?;
        let a_vec = (0..m)
            .map(|_| {
                let i_0 = Scalar::random(OsRng);
                VecElem { i_0, i_1: -i_0 }
            })
            .collect::<Vec<VecElem>>();

        // NOTE max exponent of 2 is 32, i.e. max ring len = 2^32
        assert!(m <= 32);
        assert!(1 << m > index);
        let b_vec = deltas(index as u32, m); // sigma vec
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

        let mut omegas = Vec::<Scalar>::with_capacity(m);
        let mut rho_vec = Vec::<Scalar>::with_capacity(m);
        let mut generators = Vec::<[AffinePoint; 2]>::with_capacity(m);
        let mut a_com_gen = ProjectivePoint::IDENTITY;
        let mut b_com_gen = ProjectivePoint::IDENTITY;
        let mut c_com_gen = ProjectivePoint::IDENTITY;
        let mut d_com_gen = ProjectivePoint::IDENTITY;
        for i in 0..m {
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
            omegas.push(Scalar::random(OsRng)); // x points for polynomial interpolation
            rho_vec.push(Scalar::random(OsRng));
        }

        let a_com = (h_point * Scalar::random(OsRng) + a_com_gen).to_affine();
        let b_com = (h_point * Scalar::random(OsRng) + b_com_gen).to_affine();
        let c_com = (h_point * Scalar::random(OsRng) + c_com_gen).to_affine();
        let d_com = (h_point * Scalar::random(OsRng) + d_com_gen).to_affine();

        hasher.update(a_com.to_bytes());
        hasher.update(b_com.to_bytes());
        hasher.update(c_com.to_bytes());
        hasher.update(d_com.to_bytes());

        let mut coeff_vecs = vec![Vec::<Scalar>::with_capacity(m); ring.len()];

        for k in 0..ring.len() {
            let mut highest_order = omegas[m - 1];
            let mut evals = vec![Scalar::ONE; m];
            for j in 0..m {
                highest_order *= omegas[m - 1];
                if k & (1 << j) == 0 {
                    evals[j] *= b_vec[j].i_0 * omegas[j] + a_vec[j].i_0;
                } else {
                    evals[j] *= b_vec[j].i_1 * omegas[j] + a_vec[j].i_1;
                }
            }
            if k == index {
                for eval in evals.iter_mut() {
                    *eval -= highest_order;
                }
            }

            // TODO unwrap
            let coeffs = interpolate(&omegas, &evals).unwrap();
            coeff_vecs.push(coeffs);
        }

        let mut x_points = Vec::<AffinePoint>::with_capacity(m);
        let mut y_points = Vec::<AffinePoint>::with_capacity(m);

        for (j, rho) in rho_vec.iter().enumerate() {
            let mut sum = ProjectivePoint::IDENTITY;
            for k in 0..ring.len() {
                sum += ring[k] * coeff_vecs[k][j];
            }
            let x = (sum + ProjectivePoint::GENERATOR * rho).to_affine();
            let y = (j_point * rho).to_affine();
            hasher.update(x.to_bytes());
            hasher.update(y.to_bytes());

            x_points.push(x);
            y_points.push(y);
        }

        // TODO unwrap
        let xi = Scalar::from_repr(hasher.finalize()).unwrap();

        let mut f_vec = Vec::<VecElem>::with_capacity(m);
        for j in 0..m {
            f_vec.push(VecElem {
                i_0: b_vec[j].i_0 * xi + a_vec[j].i_0,
                i_1: b_vec[j].i_1 * xi + a_vec[j].i_1,
            });
        }

        let z_a = b_com * xi + a_com;
        let z_c = c_com * xi + d_com;

        let mut xi_pow = Scalar::ONE;
        let z_sum = rho_vec.iter().fold(Scalar::ZERO, |acc, &r| {
            let next = r * xi_pow;
            xi_pow *= xi;
            acc + next
        });

        // note x_pow here should be xi^m due to stuff in fold
        let z_scalar = privkey * xi_pow - z_sum;

        Ok(Self {
            a_commitment: a_com,
            b_commitment: b_com,
            c_commitment: c_com,
            d_commitment: d_com,
            x_points,
            y_points,
            f_scalars: f_vec,
            z_a_point: z_a.to_affine(),
            z_c_point: z_c.to_affine(),
            z_scalar,
        })
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
