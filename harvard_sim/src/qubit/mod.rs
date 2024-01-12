use std::cmp::Ordering;

use strum_macros::EnumIter;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, EnumIter)]
pub enum Color {
    Red,
    Blue,
    Green,
}

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct Qubit {
    pub index: u32,
    pub color: Color,
}

impl Qubit {
    fn assign_color(index: u32) -> Color {
        let colors = [Color::Red, Color::Blue, Color::Green];
        colors[(index % 3) as usize]
    }

    pub fn new(index: u32) -> Self {
        Self {
            index,
            color: Self::assign_color(index),
        }
    }
}

impl Ord for Qubit {
    fn cmp(&self, other: &Self) -> Ordering {
        (self.index).cmp(&other.index)
    }
}

impl PartialOrd for Qubit {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
