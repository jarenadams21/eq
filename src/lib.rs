// src/lib.rs

use std::f64::consts::PI;
use serde::{Deserialize, Serialize};
use std::path::Path;
use std::fs::File;
use std::io::{BufReader, BufRead};

// Physical constants
pub const H_BAR: f64 = 1.054571e-34; // Planck's constant over 2π (J·s)
pub const E_CHARGE: f64 = 1.602176634e-19; // Electron charge (C)
pub const K_B: f64 = 1.380649e-23; // Boltzmann constant (J/K)

// Define the SimulationParams struct
#[derive(Deserialize, Serialize, Clone, Debug)]
pub struct SimulationParams {
    pub U: f64,         // Charging energy (eV)
    pub epsilon_0: f64, // Inter-shell spacing (eV)
    pub temperature: f64, // Temperature (K)
    pub eta: f64,       // Bias asymmetry
    pub gamma_l: f64,   // Tunneling rate at left lead (eV)
    pub gamma_r: f64,   // Tunneling rate at right lead (eV)
    pub delta_phi: f64, // Phase difference (radians)
    pub gamma_rel: f64, // Relaxation rate (eV)
    pub J: f64,         // Exchange interaction (eV)
    pub bias_voltage_range: (f64, f64), // Bias voltage range (V)
    pub gate_voltage_range: (f64, f64), // Gate voltage range (V)
    pub bias_steps: usize, // Number of bias voltage steps
    pub gate_steps: usize, // Number of gate voltage steps
    pub quantum_states_file: String, // Path to CSV file with valid quantum states
}

// Define a struct for QuantumState with Q-head and Q-tail
#[derive(Clone, Debug)]
pub struct QuantumState {
    pub head: f64, // Q-head encoding
    pub tail: f64, // Q-tail encoding
}

impl QuantumState {
    // Function to create a new QuantumState from given values
    pub fn new(head: f64, tail: f64) -> Self {
        QuantumState { head, tail }
    }

    // Function to simulate state transition
    pub fn transition(&self, delta: f64) -> Self {
        // Simulate state transition using a simple model
        QuantumState {
            head: self.head * delta.cos() - self.tail * delta.sin(),
            tail: self.head * delta.sin() + self.tail * delta.cos(),
        }
    }

    // Function to compute the amplitude (e.g., magnitude of the complex number)
    pub fn amplitude(&self) -> f64 {
        (self.head.powi(2) + self.tail.powi(2)).sqrt()
    }
}

// Function to load valid quantum states from CSV file
pub fn load_quantum_states(filename: &str) -> Result<Vec<QuantumState>, Box<dyn std::error::Error>> {
    let path = Path::new(filename);
    let file = File::open(&path)?;
    let reader = BufReader::new(file);

    let mut quantum_states = Vec::new();

    for line in reader.lines() {
        let line = line?;
        // Assume CSV format: head, tail
        let parts: Vec<&str> = line.trim().split(',').collect();
        if parts.len() != 2 {
            continue; // Skip invalid lines
        }
        let head = parts[0].parse::<f64>()?;
        let tail = parts[1].parse::<f64>()?;
        let state = QuantumState::new(head, tail);
        quantum_states.push(state);
    }
    Ok(quantum_states)
}

// Function to compute current according to the model, incorporating QuantumState
pub fn compute_current(
    V_b: f64,
    V_g: f64,
    params: &SimulationParams,
    quantum_states: &[QuantumState],
) -> f64 {
    // Unpack parameters
    let gamma_l = params.gamma_l; // Tunneling rate at left lead (eV)
    let gamma_r = params.gamma_r; // Tunneling rate at right lead (eV)
    let delta_phi = params.delta_phi; // Phase difference (rad)
    let gamma_rel = params.gamma_rel; // Relaxation rate (eV)
    let eta = params.eta;

    // Convert bias voltage to eV
    let V_b_eV = V_b;

    // Gate voltage effect
    let alpha_g = 1.0; // Lever arm factor (dimensionless)
    let E0 = -alpha_g * V_g; // Energy level position (eV)

    // Chemical potentials
    let mu_L = eta * V_b_eV;
    let mu_R = (eta - 1.0) * V_b_eV;

    // Determine if we are within the resonance window
    let resonance_width = 0.5e-3; // Resonance width in eV
    let in_resonance = (E0 - mu_L).abs() < resonance_width && (E0 - mu_R).abs() < resonance_width;

    // Simulate quantum state transitions
    let mut current = 0.0;
    for state in quantum_states {
        let new_state = state.transition(V_b_eV * V_g);
        let amplitude = new_state.amplitude();

        // Calculate current contribution from this state
        let numerator = gamma_l * gamma_r * delta_phi.sin().powi(2) * amplitude;
        let denominator = gamma_l + gamma_r + gamma_rel;
        let base_current = numerator / denominator;

        // Apply current suppression due to CPT
        let state_current = if in_resonance {
            base_current * 0.0 // Suppress current when CPT occurs
        } else {
            base_current
        };

        current += state_current;
    }

    // Normalize current by the number of states
    current /= quantum_states.len() as f64;

    // Return the current (in arbitrary units)
    current
}
