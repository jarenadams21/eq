
// Probe Detuning Equation : See 'probeDetuning.png'

// Prob Phase Shift Iteration 1 Code
fn calculate_probe_phase_shift_1(num_excited_atoms: f64, probe_detuning: f64) -> f64 {
    // Susceptibility χ = N * μ² / (ε₀ ℏ (Δω + iΓ/2))
    let delta = probe_detuning;
    let gamma = GAMMA;

    let numerator = num_excited_atoms * DIP_TRANSITION.powi(2);
    let denominator = EPSILON_0 * H_BAR * Complex::new(delta, gamma / 2.0); // Denominator is complex

    let susceptibility = Complex::new(numerator, 0.0) / denominator;

    // Phase shift Δϕ = ω_p * Im[χ] * L / (2c)
    let path_length = 1e-3; // Length of the cloud (m)
    let probe_frequency = FREQUENCY_0 + probe_detuning;
    let delta_phase = probe_frequency * susceptibility.im * path_length / (2.0 * C);

    delta_phase
}

// Probe Phase Shift Iteration 2 Code
// Converts the probe detuning delta from Hz to rad/s
// Function to calculate probe phase shift
fn calculate_probe_phase_shift_2(num_excited_atoms: f64, probe_detuning: f64) -> f64 {
    // Convert detuning from Hz to rad/s
    let delta = probe_detuning * 2.0 * PI; // Δω in rad/s
    let gamma = GAMMA; // Γ in s⁻¹ (already in rad/s if GAMMA was defined correctly)

    let numerator = num_excited_atoms * DIP_TRANSITION.powi(2);
    let denominator = EPSILON_0 * H_BAR * Complex::new(delta, gamma / 2.0);

    let susceptibility = Complex::new(numerator, 0.0) / denominator;

    // Phase shift Δϕ = ω_p * Im[χ] * L / (2c)
    // Convert probe_frequency to rad/s
    let probe_frequency = (FREQUENCY_0 + probe_detuning) * 2.0 * PI;
    let path_length = 1e-3; // Length of the cloud (m)
    let delta_phase = probe_frequency * susceptibility.im * path_length / (2.0 * C);

    delta_phase
}



impl bbQuantumState {
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

use std::f64::consts::PI;
use std::f64::consts::SQRT_2; // Square root of 2
use std::f64::consts::E; // Euler's number

// Golden ratio
pub const GOLDEN_RATIO: f64 = 1.618033988749895;

impl bbbQuantumState {
    // Function to simulate state transition
    pub fn transition(&self, delta: f64) -> Self {
        // Use irrational numbers to introduce complex rotations
        let phi_rotation = (delta * GOLDEN_RATIO).sin_cos();
        let euler_rotation = (delta * E).sin_cos();

        QuantumState {
            head: self.head * phi_rotation.0 - self.tail * euler_rotation.1,
            tail: self.head * euler_rotation.0 + self.tail * phi_rotation.1,
        }
    }
}
