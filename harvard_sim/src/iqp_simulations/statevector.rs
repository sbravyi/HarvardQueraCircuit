use super::simulation_params::SimulationParams;
use bitvec::prelude::*;
use rand::Rng;

pub fn generate_random_statevector(params: &SimulationParams) -> BitVec {
    #[cfg(debug_assertions)]
    return bits![0,0,1,1,1,1,0,0,1,0,0,0].to_bitvec();
    let mut rng = rand::thread_rng();
    let mut bv = bitvec![usize, Lsb0; 0; params.n_qubits as usize];
    let mut bits = params.n_qubits;
    while bits > 0 {
        let rand: u128 = rng.gen_range(0..(2_u128.pow(params.n_qubits)));
        bv.store(rand);
        bits = bits.saturating_sub(128);
        bv.shift_left(bits as usize);
    }
    bv
}
