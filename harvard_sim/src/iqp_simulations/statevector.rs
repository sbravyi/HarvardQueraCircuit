use super::simulation_params::SimulationParams;
use bitvec::prelude::*;
use rand::Rng;

pub fn generate_random_statevector(params: &SimulationParams) -> BitVec {
    let mut rng = rand::thread_rng();
    let rand: u32 = rng.gen_range(0..params.n_qubits);
    let mut bv = bitvec![usize, Lsb0; 0; params.n_qubits.ilog2() as usize ];
    bv.store(rand);
    bv
}
