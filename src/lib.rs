// src/lib.rs
// NOTE: Deprecated, Needs to be Upgraded means it requires extensive proofreading and careful addition into the current system

use std::f64::consts::{PI, E};
use serde::{Deserialize, Serialize};
use std::path::Path;
use std::fs::File;
use std::io::{BufReader, BufRead};
use num_complex::Complex64;

// Physical constants
pub const H_BAR: f64 = 1.054571e-34;      // Planck's constant over 2π (J·s)
pub const E_CHARGE: f64 = 1.602176634e-19; // Electron charge (C)
pub const K_B: f64 = 1.380649e-23;        // Boltzmann constant (J/K)
pub const GOLDEN_RATIO: f64 = 1.618033988749895;

// Gas constant as per USSA1976
pub const R_STAR: f64 = 8.31432e3; // N·m·kmol⁻¹·K⁻¹ or J·K⁻¹·kmol⁻¹

// Adjusted gas constant (after 2019 SI redefinition)
pub const R_STANDARD: f64 = 8.31446261815324; // J·K⁻¹·mol⁻¹

// Mathematical constants and functions
pub const GAMMA_CONST: f64 = 0.5772156649; // Euler-Mascheroni constant

// Define the SimulationParams struct
#[derive(Deserialize, Serialize, Clone, Debug)]
pub struct SimulationParams {
    pub U: f64,         // Charging energy (eV)
    pub epsilon_0: f64, // Inter-shell spacing (eV)
    pub eta: f64,       // Bias asymmetry
    pub gamma_l: f64,   // Tunneling rate at left lead (eV)
    pub gamma_r: f64,   // Tunneling rate at right lead (eV)
    pub delta_phi: f64, // Phase difference (radians)
    pub gamma_rel: f64, // Relaxation rate (eV)
    pub pressure: f64,  // Pressure (Pa)
    pub volume: f64,    // Volume (m³)
    pub temperature: f64, // Temperature (K)
    pub quantum_states_file: String, // Path to CSV file with valid quantum states
}

// Define a struct for QuantumState with Q-head and Q-tail
#[derive(Clone, Debug)]
pub struct QuantumState {
    pub state: Complex64, // The quantum state as a complex number
    pub zeta_zero: f64,   // Non-trivial zero of the zeta function
    pub n_moles: f64,     // Amount of substance (kmol)
}

impl QuantumState {
    // Function to create a new QuantumState from given values
    pub fn new(real: f64, imag: f64, zeta_zero: f64, n_moles: f64) -> Self {
        QuantumState {
            state: Complex64::new(real, imag),
            zeta_zero,
            n_moles,
        }
    }

    // Function to simulate state transition incorporating ideal gas behavior and QED effects
    pub fn transition(&self, delta: f64, params: &SimulationParams) -> Self {
        // Ideal gas law: PV = nRT
        let pressure = params.pressure;
        let volume = params.volume;
        let temperature = params.temperature;

        // Use R* as per USSA1976 (in J·K⁻¹·kmol⁻¹)
        let n = self.n_moles; // Amount in kmol

        // Compute the discrepancy in pressure using different R values
        let pressure_standard = (n * R_STANDARD * temperature) / volume; // Using standard R
        let pressure_star = (n * R_STAR * temperature) / volume;         // Using R*
        let pressure_difference = pressure_standard - pressure_star;     // Discrepancy due to R value

        // QED effects: Incorporate fine-structure constant (α)
        let fine_structure_constant = 7.2973525693e-3; // Dimensionless

        // Adjust delta with pressure difference and fine-structure constant
        let delta_adjusted = delta + pressure_difference * fine_structure_constant;

        // Use the zeta zero to modulate the phase
        // This models the influence of the non-trivial zeros of the Riemann zeta function
        // on the quantum state. See:
        // - "The Riemann Zeta-Function: Theory and Applications" by A. Ivic
        // - "Quantum Chaos and the Riemann Zeta Function" for connections between quantum states and zeta zeros

        let energy_scale = 1.0e-3; // Scale factor to convert zeta zeros to eV
        let scaled_zeta_zero = self.zeta_zero * energy_scale;  
        let phase_shift = delta_adjusted * scaled_zeta_zero;

        // Update the state by rotating it in the complex plane
        let rotation = Complex64::from_polar(1.0, phase_shift);

        let new_state = self.state * rotation;

        QuantumState {
            state: new_state,
            zeta_zero: self.zeta_zero,
            n_moles: self.n_moles,
        }
    }

    // Function to compute the amplitude (magnitude of the complex number)
    pub fn amplitude(&self) -> f64 {
        self.state.norm()
    }
}

// Implement factorial function for u64
trait Factorial {
    fn factorial(self) -> Self;
}

impl Factorial for u64 {
    fn factorial(self) -> Self {
        (1..=self).product()
    }
}

// Gamma function approximation using Lanczos approximation
pub fn gamma(z: f64) -> f64 {
    if z < 0.5 {
        PI / ((PI * z).sin() * gamma(1.0 - z))
    } else {
        let z = z - 1.0;
        let x = 0.99999999999980993
            + 676.5203681218851 / (z + 1.0)
            - 1259.1392167224028 / (z + 2.0)
            + 771.32342877765313 / (z + 3.0)
            - 176.61502916214059 / (z + 4.0)
            + 12.507343278686905 / (z + 5.0)
            - 0.13857109526572012 / (z + 6.0)
            + 9.9843695780195716e-6 / (z + 7.0)
            + 1.5056327351493116e-7 / (z + 8.0);
        let t = z + 7.5;
        (2.0 * PI).sqrt() * t.powf(z + 0.5) * (-t).exp() * x
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
        let line = line.trim();

        // Skip comments and empty lines
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        // Assume CSV format: real, imag, zeta_zero, n_moles
        let parts: Vec<&str> = line.split(',').collect();
        if parts.len() != 4 {
            continue; // Skip invalid lines
        }

        let real = parts[0].parse::<f64>()?;
        let imag = parts[1].parse::<f64>()?;
        let zeta_zero = parts[2].parse::<f64>()?;
        let n_moles = parts[3].parse::<f64>()?;
        let state = QuantumState::new(real, imag, zeta_zero, n_moles);
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
    let gamma_l = params.gamma_l;
    let gamma_r = params.gamma_r;
    let delta_phi = params.delta_phi;
    let gamma_rel = params.gamma_rel;
    let eta = params.eta;

    // Convert bias voltage to eV
    let V_b_eV = V_b;

    // Gate voltage effect
    let alpha_g = 1.0;
    let E0 = -alpha_g * V_g;

    // Chemical potentials
    let mu_L = eta * V_b_eV;
    let mu_R = (eta - 1.0) * V_b_eV;

    // Determine if we are within the resonance window
    // Compute thermal energy in eV
    let k_b_t_eV = K_B * params.temperature / E_CHARGE; // Convert J to eV

    // Determine if we are within the resonance window
    let resonance_width = gamma_l + gamma_r + k_b_t_eV; // eV

    // Use logical OR for resonance condition
    let in_resonance = (E0 - mu_L).abs() < resonance_width || (E0 - mu_R).abs() < resonance_width;

    // Simulate quantum state transitions
    let mut current = 0.0;
    for state in quantum_states {
        // Incorporate ideal gas behavior into delta
        let delta = V_b_eV * V_g;

        let new_state = state.transition(delta, params);
        let amplitude = new_state.amplitude();

        // Calculate current contribution from this state
        let numerator = gamma_l * gamma_r * delta_phi.sin().powi(2) * amplitude;
        let denominator = gamma_l + gamma_r + gamma_rel;
        let base_current = numerator / denominator;

        // Apply current suppression due to CPT
        let state_current = if in_resonance {
            0.0 // Suppress current when CPT occurs
        } else {
            base_current
        };

        current += state_current;
    }

    // Normalize current by the number of states
    current /= quantum_states.len() as f64;

    current
}

// Riemann zeta function (ζ(s)) for Re(s) > 1
pub fn zeta(s: f64) -> f64 {
    let mut sum = 0.0;
    for n in 1..1000000 {
        sum += 1.0 / (n as f64).powf(s);
    }
    sum
}