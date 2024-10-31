// main.rs

use std::f64::consts::PI;
use plotters::prelude::*; // For plotting
use num_complex::Complex; // For complex numbers
use serde::{Deserialize, Serialize};
use log::{info, error};
use env_logger::Env;

// Define the SimulationParams struct
#[derive(Deserialize, Serialize, Clone, Debug)]
struct SimulationParams {
    pulse_duration: f64, // in seconds
    optical_depth: f64,
    probe_detuning: f64, // in Hz
}

// Physical constants
const H_BAR: f64 = 1.054571e-34; // Planck's constant over 2π (J·s)
const C: f64 = 299792458.0; // Speed of light in vacuum (m/s)
const EPSILON_0: f64 = 8.854187817e-12; // Vacuum permittivity (F/m)
const GAMMA: f64 = 1.0 / (26e-9); // Decay rate (s⁻¹)
const DIP_TRANSITION: f64 = 3.584e-29; // Dipole transition matrix element (C·m)
const WAVELENGTH_0: f64 = 780e-9; // Central wavelength (m)
const FREQUENCY_0: f64 = C / WAVELENGTH_0; // Central frequency (Hz)
const ATOMIC_LIFETIME: f64 = 26e-9; // Atomic lifetime (s)

fn main() {
    // Initialize the logger with default level "info"
    env_logger::Builder::from_env(Env::default().default_filter_or("info")).init();

    // Define static simulation parameters
    let params = SimulationParams {
        pulse_duration: 52e-9,    // Example
        optical_depth: 4.0,      // Example
        probe_detuning: -100e7,      // Example
    };

    // Log the start of the simulation
    info!("Starting simulation with parameters: {:?}", params);

    // Run the simulation
    match run_simulation(params) {
        Ok(_) => info!("Simulation completed successfully."),
        Err(e) => error!("Simulation failed: {:?}", e),
    }
}

// Function to run the simulation
fn run_simulation(params: SimulationParams) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    // Log start of simulation
    info!("Starting simulation with params: {:?}", params);

    // Unpack parameters
    let pulse_duration = params.pulse_duration; // Pulse duration (s)
    let optical_depth = params.optical_depth;   // Optical depth
    let probe_detuning = params.probe_detuning; // Probe detuning (Hz)

    // Simulation parameters based on Angulo's experiment
    let num_atoms = 1e6 as usize; // Number of atoms in the cloud
    let cloud_length = 1e-3;      // Length of the atomic cloud (m)
    let cloud_area = PI * (25e-6_f64).powi(2); // Cross-sectional area (m²)
    let _atom_density = num_atoms as f64 / (cloud_length * cloud_area); // Atom density (atoms/m³)

    // Signal parameters
    let signal_frequency = FREQUENCY_0; // Signal frequency (Hz)
    let signal_photon_energy = H_BAR * signal_frequency; // Signal photon energy (J)
    let mean_photon_number = 100.0;    // Mean photon number per pulse (from Angulo's experiment)

    // Time parameters
    let time_window = (-100e-9, 100e-9); // Time window for simulation (s)
    let time_steps = 1000;               // Number of time steps
    let delta_t = (time_window.1 - time_window.0) / (time_steps as f64); // Time step size (s)

    // Generate time vector
    let times: Vec<f64> = (0..time_steps)
        .map(|i| time_window.0 + i as f64 * delta_t)
        .collect();

    // Initialize vector for phase shifts
    let phase_shifts: Vec<f64> = times
        .iter()
        .map(|&time| {
            let signal_intensity = signal_pulse_intensity(
                time,
                pulse_duration,
                mean_photon_number,
                signal_photon_energy,
            );

            // Calculate atomic excitation probability
            let excitation_prob = calculate_excitation_probability(signal_intensity, delta_t);

            // Calculate number of excited atoms
            let num_excited_atoms = excitation_prob * num_atoms as f64;

            // Calculate probe phase shift due to excited atoms
            calculate_probe_phase_shift(num_excited_atoms, probe_detuning)
        })
        .collect();

    // Output results
    info!(
        "Pulse duration: {:.0} ns, Optical depth: {}, Probe detuning: {:.2} MHz",
        pulse_duration * 1e9,
        optical_depth,
        probe_detuning / 1e6,
    );

    // Plot the results
    let filename = format!(
        "phase_shift_{:.0}_ns_OD_{:.0}.png",
        pulse_duration * 1e9,
        optical_depth
    );
    let root_area = BitMapBackend::new(&filename, (1280, 720)).into_drawing_area();
    root_area.fill(&WHITE)?;
    let min_phase = phase_shifts.iter().cloned().fold(f64::INFINITY, f64::min);
    let max_phase = phase_shifts.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let mut chart = ChartBuilder::on(&root_area)
        .caption(
            format!(
                "Probe Phase Shift Over Time (Pulse: {:.0} ns, OD: {:.0})",
                pulse_duration * 1e9,
                optical_depth
            ),
            ("sans-serif", 30).into_font(),
        )
        .margin(20)
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d(
            time_window.0 * 1e9..time_window.1 * 1e9,
            min_phase..max_phase,
        )?;

    chart
        .configure_mesh()
        .x_desc("Time (ns)")
        .y_desc("Phase Shift (rad)")
        .draw()?;

    chart.draw_series(LineSeries::new(
        times.iter().map(|t| t * 1e9).zip(phase_shifts.iter().cloned()),
        &BLUE,
    ))?;

    // Save data to file
    let data_filename = format!(
        "data_{:.0}_ns_OD_{:.0}.csv",
        pulse_duration * 1e9,
        optical_depth
    );
    save_data(&data_filename, &times, &phase_shifts)?;

    info!("Simulation data saved to files: {}, {}", filename, data_filename);

    Ok(())
}

// Function for signal pulse intensity (Gaussian pulse)
fn signal_pulse_intensity(
    time: f64,
    pulse_duration: f64,
    mean_photon_number: f64,
    photon_energy: f64,
) -> f64 {
    let pulse_amplitude =
        (mean_photon_number * photon_energy) / (pulse_duration * (2.0 * PI).sqrt());
    pulse_amplitude * (-0.5 * (time / pulse_duration).powi(2)).exp()
}

// Function to calculate excitation probability
fn calculate_excitation_probability(signal_intensity: f64, delta_t: f64) -> f64 {
    // Rabi frequency Ω = μ * E / ℏ
    let electric_field = (2.0 * signal_intensity / (C * EPSILON_0)).sqrt();
    let rabi_frequency = DIP_TRANSITION * electric_field / H_BAR;

    // Excitation probability per atom P = (Ω * delta_t / 2)²
    let excitation_prob = (rabi_frequency * delta_t / 2.0).powi(2);
    excitation_prob.min(1.0) // Ensure probability does not exceed 1
}

// Function to calculate probe phase shift
fn calculate_probe_phase_shift(num_excited_atoms: f64, probe_detuning: f64) -> f64 {
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


// Function to save data to CSV
fn save_data(filename: &str, times: &[f64], phase_shifts: &[f64]) -> Result<(), std::io::Error> {
    use std::fs::File;
    use std::io::Write;

    let mut file = File::create(filename)?;
    writeln!(file, "Time(s),PhaseShift(rad)")?;
    for (&time, &phase_shift) in times.iter().zip(phase_shifts.iter()) {
        writeln!(file, "{},{}", time, phase_shift)?;
    }
    Ok(())
}
