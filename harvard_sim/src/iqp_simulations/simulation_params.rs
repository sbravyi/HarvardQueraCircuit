pub struct SimulationParams {
    pub boolean_cube_dimension: u32,
    pub nodes: u32,
    pub n_qubits: u32,
}

impl SimulationParams {
    pub fn new(boolean_cube_dimension: u32) -> Self {
        let nodes = 1 << boolean_cube_dimension;
        let n_qubits = 3 * nodes;
        Self {
            boolean_cube_dimension,
            nodes,
            n_qubits,
        }
    }
}
