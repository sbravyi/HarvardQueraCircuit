use anyhow::Result;

pub trait Simulation {
    fn run(&mut self) -> Result<f64>;
}
