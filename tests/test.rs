// tests/test.rs

use std::f64::consts::PI;

// Import the functions and constants from your crate
use funhalf::{compute_current, SimulationParams, QuantumState};

#[test]
fn test_current_calculation() {
    let params = SimulationParams {
        U: 2.0e-3,             // Charging energy (eV)
        epsilon_0: 3.0e-3,     // Inter-shell spacing (eV)
        temperature: 2.7,      // Temperature (K)
        eta: 0.55,             // Bias asymmetry
        gamma_l: 0.2556e-6,    // Tunneling rate at left lead (eV)
        gamma_r: 0.0244e-6,    // Tunneling rate at right lead (eV)
        delta_phi: PI / (2.0_f64.sqrt()),   // Phase difference (radians)
        gamma_rel: 0.0,        // Relaxation rate (eV)
        J: 0.01e-3,            // Exchange interaction (eV)
        bias_voltage_range: (-5e-3, 5e-3), // Bias voltage range (V)
        gate_voltage_range: (-20e-3, 20e-3), // Gate voltage range (V)
        bias_steps: 200,       // Number of bias voltage steps
        gate_steps: 200,       // Number of gate voltage steps
        quantum_states_file: String::new(), // Not used in test
    };

    let V_b = 1e-3; // Bias voltage (V)
    let V_g = 0.0;  // Gate voltage (V)

    // Create a simple quantum state for testing
    // Interesting Note: Only (>= 1 in the Y causes current > 0)
    let quantum_states = vec![
        QuantumState::new(0.05, -0.10),
        QuantumState::new(0.05, -0.10),
    ];

    let current = compute_current(V_b, V_g, &params, &quantum_states);

    // Check that the current is non-zero
    assert!(current > 0.0);

    // Now test in resonance condition where CPT occurs
    let V_g_resonance = 0.0; // Assuming resonance at V_g = 0
    let V_b_resonance = 0.0; // Assuming resonance at V_b = 0
    let current_cpt = compute_current(V_b_resonance, V_g_resonance, &params, &quantum_states);

    // Check that the current is suppressed due to CPT
    // (!) Current test failures:
    //  1. left: 2.7910636024282933e-8
    //  2. right: 0.0
    assert_eq!(current_cpt, 0.0);
}
