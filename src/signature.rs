use crate::ring::*;
use k256::elliptic_curve::group::GroupEncoding;
use k256::elliptic_curve::ops::Reduce;
use k256::elliptic_curve::{Field, PrimeField};
use k256::{AffinePoint, ProjectivePoint, Scalar, U256};
use rand_core::OsRng;
use sha3::{Digest, Keccak256};

// NOTE this is a public "private key" that determines the U point which tags the
// signatures
const U_SCALAR_U256: U256 =
    U256::from_le_hex("7c81a9587b8da43a9519bd50d96191fd8f2c4f66b8f1550e366e3c7f9ed18897");

pub struct Parameters {
    generators: Vec<VecElem<AffinePoint>>,
    h_point: AffinePoint,
    u_point: AffinePoint,
}

impl Parameters {
    pub fn new(n_generators: usize) -> Self {
        let mut generators = Vec::<VecElem<AffinePoint>>::with_capacity(n_generators);
        for _ in 0..n_generators {
            let pt_0 = AffinePoint::GENERATOR * Scalar::random(OsRng);
            let pt_1 = AffinePoint::GENERATOR * Scalar::random(OsRng);
            generators.push(VecElem {
                i_0: pt_0.to_affine(),
                i_1: pt_1.to_affine(),
            });
        }

        Self {
            generators,
            h_point: (AffinePoint::GENERATOR * Scalar::random(OsRng)).to_affine(),
            u_point: (AffinePoint::GENERATOR * Scalar::from_uint_reduced(U_SCALAR_U256))
                .to_affine(),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct VecElem<T: Clone + Copy + std::fmt::Debug + PartialEq + Eq> {
    i_0: T,
    i_1: T,
}

pub struct Signature {
    a_commitment: AffinePoint,
    b_commitment: AffinePoint,
    c_commitment: AffinePoint,
    d_commitment: AffinePoint,
    x_points: Vec<AffinePoint>,
    y_points: Vec<AffinePoint>,
    f_scalars: Vec<VecElem<Scalar>>,
    z_a_scalar: Scalar,
    z_c_scalar: Scalar,
    z_scalar: Scalar,
    tag: AffinePoint,
}

// NOTE n = 2, i.e. N = 2^m
// N - number of elements in the ring (should be padded to 2^m)
impl Signature {
    pub fn new(
        index: usize,
        ring: &Ring,
        ring_hash: &[u8],
        message_hash: &[u8],
        privkey: Scalar,
        parameters: &Parameters,
    ) -> Result<Self, String> {
        let mut hasher = Keccak256::new();
        hasher.update(message_hash);
        hasher.update(ring_hash);

        // TODO unwrap
        let j_point = parameters.u_point * privkey.invert().unwrap();

        let mut ring = ring.to_owned();
        let m = pad_ring_to_2n(&mut ring)?;
        let a_vec = (0..m)
            .map(|_| {
                let i_0 = Scalar::random(OsRng);
                VecElem { i_0, i_1: -i_0 }
            })
            .collect::<Vec<VecElem<Scalar>>>();

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
            .collect::<Vec<VecElem<Scalar>>>();

        let d_vec = a_vec
            .iter()
            .map(|a| {
                let i_0 = -(a.i_0 * a.i_0);
                let i_1 = -(a.i_1 * a.i_1);
                VecElem { i_0, i_1 }
            })
            .collect::<Vec<VecElem<Scalar>>>();

        let mut omegas = Vec::<Scalar>::with_capacity(m);
        let mut rho_vec = Vec::<Scalar>::with_capacity(m);
        let mut a_com_gen = ProjectivePoint::IDENTITY;
        let mut b_com_gen = ProjectivePoint::IDENTITY;
        let mut c_com_gen = ProjectivePoint::IDENTITY;
        let mut d_com_gen = ProjectivePoint::IDENTITY;
        for (i, gen) in parameters.generators.iter().take(m).enumerate() {
            a_com_gen += gen.i_0 * a_vec[i].i_0;
            a_com_gen += gen.i_1 * a_vec[i].i_1;
            b_com_gen += gen.i_0 * b_vec[i].i_0;
            b_com_gen += gen.i_1 * b_vec[i].i_1;
            c_com_gen += gen.i_0 * c_vec[i].i_0;
            c_com_gen += gen.i_1 * c_vec[i].i_1;
            d_com_gen += gen.i_0 * d_vec[i].i_0;
            d_com_gen += gen.i_1 * d_vec[i].i_1;
            omegas.push(Scalar::random(OsRng)); // x points for polynomial interpolation
            rho_vec.push(Scalar::random(OsRng));
        }

        let r_a = Scalar::random(OsRng);
        let r_b = Scalar::random(OsRng);
        let r_c = Scalar::random(OsRng);
        let r_d = Scalar::random(OsRng);

        let a_com = (parameters.h_point * r_a + a_com_gen).to_affine();
        let b_com = (parameters.h_point * r_b + b_com_gen).to_affine();
        let c_com = (parameters.h_point * r_c + c_com_gen).to_affine();
        let d_com = (parameters.h_point * r_d + d_com_gen).to_affine();

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
            coeff_vecs[k] = coeffs;
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

        let mut f_vec = Vec::<VecElem<Scalar>>::with_capacity(m);
        for j in 0..m {
            f_vec.push(VecElem {
                i_0: b_vec[j].i_0 * xi + a_vec[j].i_0,
                i_1: b_vec[j].i_1 * xi + a_vec[j].i_1,
            });
        }

        let z_a = r_b * xi + r_a;
        let z_c = r_c * xi + r_d;

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
            z_a_scalar: z_a,
            z_c_scalar: z_c,
            z_scalar,
            tag: j_point.to_affine(),
        })
    }

    pub fn verify(
        &self,
        ring: &Ring,
        ring_hash: &[u8],
        message_hash: &[u8],
        parameters: &Parameters,
    ) -> Result<(), String> {
        let mut hasher = Keccak256::new();
        hasher.update(message_hash);
        hasher.update(ring_hash);

        let mut ring = ring.to_owned();
        let m = pad_ring_to_2n(&mut ring)?;

        hasher.update(self.a_commitment.to_bytes());
        hasher.update(self.b_commitment.to_bytes());
        hasher.update(self.c_commitment.to_bytes());
        hasher.update(self.d_commitment.to_bytes());

        for (x, y) in self.x_points.iter().zip(self.y_points.iter()) {
            hasher.update(x.to_bytes());
            hasher.update(y.to_bytes());
        }

        let xi = Scalar::from_repr(hasher.finalize()).unwrap();

        let (f_0_0, f_1_0) = self
            .f_scalars
            .iter()
            .fold((Scalar::ZERO, Scalar::ZERO), |(acc_0, acc_1), elem| {
                (acc_0 + elem.i_0, acc_1 + elem.i_1)
            });

        let f_0_0 = xi - f_0_0;
        let f_1_0 = xi - f_1_0;

        // check commitments
        let mut gen_0_f = parameters.generators[0].i_0 * f_0_0;
        let mut gen_1_f = parameters.generators[0].i_1 * f_1_0;
        let mut gen_0_f_xi = parameters.generators[0].i_0 * (f_0_0 * (xi - f_0_0));
        let mut gen_1_f_xi = parameters.generators[0].i_1 * (f_1_0 * (xi - f_1_0));

        for (g, f) in parameters
            .generators
            .iter()
            .take(m)
            .skip(1)
            .zip(self.f_scalars.iter().skip(1))
        {
            gen_0_f += g.i_0 * f.i_0;
            gen_1_f += g.i_1 * f.i_1;
            gen_0_f_xi += g.i_0 * (f.i_0 * (xi - f.i_0));
            gen_1_f_xi += g.i_1 * (f.i_1 * (xi - f.i_1));
        }

        let com_f = (gen_0_f + gen_1_f + parameters.h_point * self.z_a_scalar).to_affine();
        let com_f_xi = (gen_0_f_xi + gen_1_f_xi + parameters.h_point * self.z_c_scalar).to_affine();

        if (ProjectivePoint::from(self.a_commitment) + self.b_commitment * xi).to_affine() != com_f
        {
            return Err("ab commitment mismatch".to_owned());
        }

        if (ProjectivePoint::from(self.d_commitment) + self.c_commitment * xi).to_affine()
            != com_f_xi
        {
            return Err("cd commitment mismatch".to_owned());
        }

        let mut sum_pk_prod_f = ProjectivePoint::IDENTITY;
        let mut sum_prod_f = Scalar::ZERO;
        for (k, &pk) in ring.iter().enumerate() {
            let prod_f = self
                .f_scalars
                .iter()
                .enumerate()
                .fold(Scalar::ZERO, |acc, (j, f)| {
                    if k & 1 << j == 0 {
                        acc * f.i_0
                    } else {
                        acc * f.i_1
                    }
                });

            sum_pk_prod_f += pk * prod_f;
            sum_prod_f += prod_f;
        }

        let mut xi_pow = Scalar::ONE; // xi^0
        let (x_sum, y_sum) = self.x_points.iter().zip(self.y_points.iter()).fold(
            (ProjectivePoint::IDENTITY, ProjectivePoint::IDENTITY),
            |(acc_0, acc_1), (&x, &y)| {
                let next_x = acc_0 + x * xi_pow;
                let next_y = acc_1 + y * xi_pow;
                xi_pow *= xi;
                (next_x, next_y)
            },
        );

        let first_zero = sum_pk_prod_f - x_sum - AffinePoint::GENERATOR * self.z_scalar;
        let second_zero = parameters.u_point * sum_prod_f - y_sum - self.tag * self.z_scalar;

        if first_zero.to_affine() != AffinePoint::IDENTITY {
            return Err("first constraint is nonzero".to_owned());
        }

        if second_zero.to_affine() != AffinePoint::IDENTITY {
            return Err("second constraint is nonzero".to_owned());
        }

        Ok(())
    }
}

fn deltas(num: u32, n: usize) -> Vec<VecElem<Scalar>> {
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
        .collect::<Vec<VecElem<Scalar>>>()
}

#[cfg(test)]
mod test {
    use super::*;

    struct Keypair {
        public: AffinePoint,
        private: Scalar,
    }

    impl Keypair {
        fn random() -> Self {
            let private = Scalar::random(OsRng);
            Self {
                public: (AffinePoint::GENERATOR * private).to_affine(),
                private,
            }
        }
    }

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

    #[test]
    fn valid_signature() {
        let parameters = Parameters::new(10);
        let keypair = Keypair::random();
        let ring = vec![
            (AffinePoint::GENERATOR * Scalar::random(OsRng)).to_affine(),
            (AffinePoint::GENERATOR * Scalar::random(OsRng)).to_affine(),
            keypair.public,
            (AffinePoint::GENERATOR * Scalar::random(OsRng)).to_affine(),
            (AffinePoint::GENERATOR * Scalar::random(OsRng)).to_affine(),
        ];
        let index = 2_usize;

        let mut hasher = Keccak256::new();
        for r in ring.iter() {
            hasher.update(r.to_bytes());
        }
        let ring_hash = hasher.finalize();

        let msg = b"hello there!";
        let mut hasher = Keccak256::new();
        hasher.update(msg);
        let msg_hash = hasher.finalize();

        let signature = Signature::new(
            index,
            &ring,
            &ring_hash,
            &msg_hash,
            keypair.private,
            &parameters,
        )
        .unwrap();
        assert_eq!(signature.f_scalars.len(), 3);
        assert_eq!(signature.x_points.len(), 3);
        assert_eq!(signature.y_points.len(), 3);
        let result = signature.verify(&ring, &ring_hash, &msg_hash, &parameters);
        println!("{:?}", result);
        assert!(result.is_ok());
    }
}
