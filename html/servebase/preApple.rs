use std::f64::consts::PI;
use rayon::prelude::*; // For parallel processing
use plotters::prelude::*; // For plotting
use num_complex::Complex; // For complex numbers

// Via research as of October 28, 2024 by Daniela Angulo
// Her research inspired this continued simulation for problem solving as you will
fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Physical constants
    const H_BAR: f64 = 1.054571817e-34; // Planck's constant over 2π (J·s)
    const C: f64 = 299792458.0; // Speed of light in vacuum (m/s)
    const EPSILON_0: f64 = 8.854187817e-12; // Vacuum permittivity (F/m)
    const GAMMA: f64 = 1.0 / (26e-9); // Decay rate (s⁻¹)
    const DIP_TRANSITION: f64 = 3.584e-29; // Dipole transition matrix element (C·m)
    const WAVELENGTH_0: f64 = 780e-9; // Central wavelength (m)
    const FREQUENCY_0: f64 = C / WAVELENGTH_0; // Central frequency (Hz)

    // Experimental parameters
    let pulse_durations = [10e-9, 18e-9, 27e-9, 36e-9]; // Pulse durations (s)
    let optical_depths = [2.0, 3.0, 4.0]; // Optical depths
    let num_atoms = 1e6 as usize; // Number of atoms in the cloud
    let cloud_length = 1e-3; // Length of the atomic cloud (m)
    let cloud_area = PI * (25e-6_f64).powi(2); // Cross-sectional area (m²)
    let atom_density = num_atoms as f64 / (cloud_length * cloud_area); // Atom density (atoms/m³)

    // Probe parameters
    let probe_detuning = -20e6; // Probe detuning (Hz)
    let probe_frequency = FREQUENCY_0 + probe_detuning; // Probe frequency (Hz)

    // Signal parameters
    let signal_frequency = FREQUENCY_0; // Signal frequency (Hz)
    let signal_photon_energy = H_BAR * signal_frequency; // Signal photon energy (J)
    let mean_photon_number = 100.0; // Mean photon number per pulse

    // Time parameters
    let time_window = (-100e-9, 100e-9); // Time window for simulation (s)
    let time_steps = 1000; // Number of time steps
    let delta_t = (time_window.1 - time_window.0) / (time_steps as f64); // Time step size (s)

    // Simulation loop over pulse durations and optical depths
    for &pulse_duration in &pulse_durations {
        for &optical_depth in &optical_depths {
            // Generate time vector
            let times: Vec<f64> = (0..time_steps)
                .map(|i| time_window.0 + i as f64 * delta_t)
                .collect();

            // Initialize vectors for phase shifts
            let phase_shifts: Vec<f64> = times
                .par_iter()
                .map(|&time| {
                    let signal_intensity = signal_pulse_intensity(
                        time,
                        pulse_duration,
                        mean_photon_number,
                        signal_photon_energy,
                    );

                    // Calculate atomic excitation probability
                    let excitation_prob =
                        calculate_excitation_probability(signal_intensity, delta_t);

                    // Calculate number of excited atoms
                    let num_excited_atoms = excitation_prob * num_atoms as f64;

                    // Calculate probe phase shift due to excited atoms
                    let delta_phase =
                        calculate_probe_phase_shift(num_excited_atoms, probe_detuning);

                    delta_phase
                })
                .collect();

            // Output results
            println!(
                "Pulse duration: {:.0} ns, Optical depth: {}",
                pulse_duration * 1e9,
                optical_depth
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
        }
    }

    Ok(())
}

// Function for signal pulse intensity
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
    const DIP_TRANSITION: f64 = 3.584e-29; // Dipole transition matrix element (C·m)
    const H_BAR: f64 = 1.054571817e-34; // Planck's constant over 2π (J·s)
    const C: f64 = 299792458.0; // Speed of light in vacuum (m/s)
    const EPSILON_0: f64 = 8.854187817e-12; // Vacuum permittivity (F/m)

    let electric_field = (2.0 * signal_intensity / (C * EPSILON_0)).sqrt();
    let rabi_frequency = DIP_TRANSITION * electric_field / H_BAR;

    // Excitation probability per atom P = (Ω * delta_t / 2)²
    let excitation_prob = (rabi_frequency * delta_t / 2.0).powi(2);
    excitation_prob.min(1.0) // Ensure probability does not exceed 1
}

// Function to calculate probe phase shift
fn calculate_probe_phase_shift(num_excited_atoms: f64, probe_detuning: f64) -> f64 {
    // Susceptibility χ = N * μ² / (ε₀ ℏ (Δω + iΓ/2))
    const DIP_TRANSITION: f64 = 3.584e-29; // Dipole transition matrix element (C·m)
    const EPSILON_0: f64 = 8.854187817e-12; // Vacuum permittivity (F/m)
    const H_BAR: f64 = 1.054571817e-34; // Planck's constant over 2π (J·s)
    const GAMMA: f64 = 1.0 / (26e-9); // Decay rate (s⁻¹)
    const FREQUENCY_0: f64 = C / 780e-9; // Central frequency (Hz)
    const C: f64 = 299792458.0; // Speed of light in vacuum (m/s)

    let delta = probe_detuning;
    let gamma = GAMMA;

    let numerator = num_excited_atoms * DIP_TRANSITION.powi(2);
    let denominator = EPSILON_0 * H_BAR
        * Complex::new(delta, gamma / 2.0); // Denominator is complex

    let susceptibility = Complex::new(numerator, 0.0) / denominator;

    // Phase shift Δϕ = ω_p * Im[χ] * L / (2c)
    let path_length = 1e-3; // Length of the cloud (m)
    let probe_frequency = FREQUENCY_0 + probe_detuning;
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