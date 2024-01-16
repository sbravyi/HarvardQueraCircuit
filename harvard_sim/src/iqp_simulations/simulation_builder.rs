use super::simulation::Simulation;
use super::small_int_simulation::CPUSmallIntSimulation;
use anyhow::{bail, Context, Result};
pub struct IQPSimulationBuilder {
    boolean_cube_dimension: u32,
}

impl IQPSimulationBuilder {
    pub fn new(boolean_cube_dimension: u32) -> Self {
        Self {
            boolean_cube_dimension,
        }
    }
}

impl IQPSimulationBuilder {
    pub fn run_appropriate_simulation_instance(self) -> Result<()> {
        if self.boolean_cube_dimension * 3 >= 128 {
            bail!(
                "no appropriate simulation implementation for IQP circuits ok size k = {}",
                self.boolean_cube_dimension
            );
        } else {
            let mut simulation = CPUSmallIntSimulation::new(self.boolean_cube_dimension);
            simulation.run().context("running small int simulation")?;
            Ok(())
        }
    }
}
