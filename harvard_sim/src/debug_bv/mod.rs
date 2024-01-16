use crate::{
    bit_matrix::{gauss_jordan::GaussJordan, matrix::BitMatrix},
    iqp_simulations::linear_systems::LinearSystems,
};
use bitvec::prelude::*;

pub fn debug_bitvec(bv: &BitVec) {
    println!("{}", bv);
}

pub fn debug_bitmatrix(u: &BitMatrix) {
    for bv in u.rows.iter() {
        println!("{}", bv);
    }
}

pub fn debug_gj(u: &GaussJordan) {
    for bv in u.rows.iter() {
        println!("{}", bv);
    }
}

pub fn debug_linear_system(ls: &LinearSystems) {
    // gamma: BitMatrix,
    // delta_b: BitVec,
    // delta_g: BitVec,
    // x_r: BitVec,
    // gj: GaussJordan,
    // solver: BackwardsSubstitution,
    // // caching allocations for
    // // often used intermediary vectors
    // sb_delta_b: BitVec,
    // sg_delta_g: BitVec,
    println!("DEBUGGING LINEAR SYSTEM!");
    println!("Gamma:");
    debug_bitmatrix(&ls.gamma);
    println!("delta_b:");
    debug_bitvec(&ls.delta_b);
    println!("delta_g:");
    debug_bitvec(&ls.delta_g);
    println!("x_r:");
    debug_bitvec(&ls.x_r);
    println!("gj:");
    debug_gj(&ls.solver.gj);
    println!("solver:");
    debug_bitvec(&ls.solver.solution);
    println!("sb_delta_b:");
    debug_bitvec(&ls.sb_delta_b);
    println!("sg_delta_g:");
    debug_bitvec(&ls.sg_delta_g);
}
