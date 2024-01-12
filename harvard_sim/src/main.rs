use std::process;

use anyhow::Result;
use clap::Parser;
use env_logger::Env;
use harvard_sim::iqp_simulations::simulation_builder::IQPSimulationBuilder;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    #[arg(short, long, default_value_t = 2)]
    pub boolean_cube_dimension: u32,
}

fn run_iqp_simulation(boolean_cube_dimension: u32) -> Result<()> {
    let simulation_factory = IQPSimulationBuilder::new(boolean_cube_dimension);
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
