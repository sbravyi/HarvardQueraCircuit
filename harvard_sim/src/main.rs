use std::process;

use anyhow::{bail, Result};
use clap::Parser;
use env_logger::Env;
use harvard_sim::gray_code_flip_bit;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    #[arg(short, long, default_value_t = 2)]
    boolean_cube_dimension: u32,
}

trait Simulation {
    fn run(&mut self) -> Result<()>;
}

struct CPUSmallIntSimulation {
    boolean_cube_dimension: u32,
    nodes: u32,
    n: u32,
}

impl CPUSmallIntSimulation {
    fn new(boolean_cube_dimension: u32) -> Self {
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
        let gc_flip_bit = gray_code_flip_bit::GrayCodeFlipBit::new(self.nodes);
        for flip_bit in gc_flip_bit {}
        Ok(())
    }
}

struct IQPSimulationBuilder {
    boolean_cube_dimension: u32,
}

impl IQPSimulationBuilder {
    fn run_appropriate_simulation_instance(self) -> Result<()> {
        let mut appropriate_simulation = if self.boolean_cube_dimension * 3 >= 128 {
            bail!(
                "no appropriate simulation implementation for IQP circuits ok size k = {}",
                self.boolean_cube_dimension
            );
        } else {
            CPUSmallIntSimulation::new(self.boolean_cube_dimension)
        };
        appropriate_simulation.run()
    }
}

fn run_iqp_simulation(boolean_cube_dimension: u32) -> Result<()> {
    let simulation_factory = IQPSimulationBuilder {
        boolean_cube_dimension,
    };
    simulation_factory.run_appropriate_simulation_instance()
}

fn main() {
    let args = Args::parse();
    env_logger::Builder::from_env(Env::default().default_filter_or("debug")).init();
    log::debug!("Running with: {args:?}");
    if let Err(err) = run_iqp_simulation(args.boolean_cube_dimension) {
        eprintln!("ran into simulation error: {err:?}");
        process::exit(1)
    }
}
