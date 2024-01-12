use super::simulation::Simulation;
use crate::gray_code_flip_bit::GrayCodeFlipBit;
use anyhow::Result;

pub struct CPUSmallIntSimulation {
    boolean_cube_dimension: u32,
    nodes: u32,
    n: u32,
}

impl CPUSmallIntSimulation {
    pub fn new(boolean_cube_dimension: u32) -> Self {
        let nodes = 1 << boolean_cube_dimension;
        let n = 3 * nodes;
        Self {
            boolean_cube_dimension,
            nodes,
            n,
        }
    }
}

impl Simulation for CPUSmallIntSimulation {
    fn run(&mut self) -> Result<()> {
        let gc_flip_bit = GrayCodeFlipBit::new(self.nodes);
        for flip_bit in gc_flip_bit {}
        Ok(())
    }
}
