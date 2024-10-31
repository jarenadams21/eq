






// main.rs

use std::f64::consts::PI;
use plotters::prelude::*; // For plotting
use num_complex::Complex; // For complex numbers
use actix_web::{web, App, HttpResponse, HttpServer, Responder};
use serde::{Deserialize, Serialize};
use std::sync::Mutex;
use log::{info, error, warn};
use env_logger;

// Simulation parameters structure
#[derive(Deserialize, Serialize, Clone, Debug)]
struct SimulationParams {
    pulse_duration: f64,
    optical_depth: f64,
    probe_detuning: f64,
}

// Shared application state
struct AppState {
    simulation_params: Mutex<Option<SimulationParams>>,
}

#[actix_web::main]
async fn main() -> std::io::Result<()> {
    // Initialize the logger
    env_logger::init();

    // Initialize shared state
    let app_state = web::Data::new(AppState {
        simulation_params: Mutex::new(None),
    });

    info!("Starting server on 127.0.0.1:8080");

    HttpServer::new(move || {
        App::new()
            .app_data(app_state.clone())
            .route("/simulate", web::post().to(start_simulation))
            .route("/clear", web::post().to(clear_media))
            .route("/results", web::get().to(get_results))
    })
    .bind("0.0.0.0:8080")?
    .run()
    .await
}

// Handler to start simulation
async fn start_simulation(
    data: web::Data<AppState>,
    params: web::Json<SimulationParams>,
) -> impl Responder {
    // Log received parameters
    info!("Received simulation parameters: {:?}", params);

    // Validate simulation parameters
    if let Err(validation_error) = validate_params(&params) {
        error!("Validation error: {}", validation_error);
        return HttpResponse::BadRequest().body(format!("Invalid parameters: {}", validation_error));
    }

    // Update simulation parameters
    let mut sim_params = data.simulation_params.lock().unwrap();
    *sim_params = Some(params.into_inner());

    // Run simulation
    if let Some(sim_params) = data.simulation_params.lock().unwrap().clone() {
        // Run the simulation with the provided parameters
        if let Err(e) = run_simulation(sim_params) {
            error!("Simulation failed: {}", e);
            return HttpResponse::InternalServerError()
                .body(format!("Simulation failed: {}", e));
        }
    }

    info!("Simulation completed successfully");

    HttpResponse::Ok().body("Simulation completed")
}

// Function to validate simulation parameters
fn validate_params(params: &SimulationParams) -> Result<(), String> {
    if params.pulse_duration <= 0.0 {
        return Err("Pulse duration must be greater than zero".to_string());
    }
    if params.optical_depth <= 0.0 {
        return Err("Optical depth must be greater than zero".to_string());
    }
    if params.probe_detuning.is_nan() {
        return Err("Probe detuning must be a valid number".to_string());
    }
    Ok(())
}

// Handler to clear media directory
async fn clear_media() -> impl Responder {
    // Path to media directory
    let media_dir = "./";

    // Remove all PNG and CSV files in the directory
    if let Err(e) = remove_files_with_extension(media_dir, "png") {
        error!("Failed to remove PNG files: {}", e);
        return HttpResponse::InternalServerError()
            .body(format!("Failed to remove PNG files: {}", e));
    }
    if let Err(e) = remove_files_with_extension(media_dir, "csv") {
        error!("Failed to remove CSV files: {}", e);
        return HttpResponse::InternalServerError()
            .body(format!("Failed to remove CSV files: {}", e));
    }

    info!("Media directory cleared");

    HttpResponse::Ok().body("Media directory cleared")
}

// Handler to get results (list of files)
async fn get_results() -> impl Responder {
    use std::fs;

    let mut results = Vec::new();
    let entries = match fs::read_dir("./") {
        Ok(entries) => entries,
        Err(e) => {
            error!("Failed to read directory: {}", e);
            return HttpResponse::InternalServerError()
                .body(format!("Failed to read directory: {}", e));
        }
    };

    for entry in entries {
        let entry = match entry {
            Ok(entry) => entry,
            Err(e) => {
                warn!("Failed to read directory entry: {}", e);
                continue;
            }
        };
        let path = entry.path();
        if let Some(ext) = path.extension() {
            if ext == "png" || ext == "csv" {
                if let Some(name) = path.file_name() {
                    results.push(name.to_string_lossy().into_owned());
                }
            }
        }
    }

    info!("Retrieved results: {:?}", results);

    HttpResponse::Ok().json(results)
}

// Function to remove files with specific extension
fn remove_files_with_extension(dir: &str, extension: &str) -> std::io::Result<()> {
    use std::fs;
    let entries = fs::read_dir(dir)?;

    for entry in entries {
        let entry = entry?;
        let path = entry.path();
        if let Some(ext) = path.extension() {
            if ext == extension {
                fs::remove_file(path)?;
            }
        }
    }
    Ok(())
}

// Function to run the simulation
fn run_simulation(params: SimulationParams) -> Result<(), Box<dyn std::error::Error>> {
    // Physical constants
    const H_BAR: f64 = 1.054571e-34; // Planck's constant over 2π (J·s)
    const C: f64 = 299792458.0; // Speed of light in vacuum (m/s)
    const EPSILON_0: f64 = 8.854187817e-12; // Vacuum permittivity (F/m)
    const GAMMA: f64 = 1.0 / (26e-9); // Decay rate (s⁻¹)
    const DIP_TRANSITION: f64 = 3.584e-29; // Dipole transition matrix element (C·m)
    const WAVELENGTH_0: f64 = 780e-9; // Central wavelength (m)
    const FREQUENCY_0: f64 = C / WAVELENGTH_0; // Central frequency (Hz)

    // Experimental parameters
    let num_atoms = 1e6 as usize; // Number of atoms in the cloud
    let cloud_length = 1e-3; // Length of the atomic cloud (m)
    let cloud_area = PI * (25e-6_f64).powi(2); // Cross-sectional area (m²)
    let _atom_density = num_atoms as f64 / (cloud_length * cloud_area); // Atom density (atoms/m³)

    // Unpack parameters
    let pulse_duration = params.pulse_duration; // Pulse duration (s)
    let optical_depth = params.optical_depth; // Optical depth
    let probe_detuning = params.probe_detuning; // Probe detuning (Hz)

    // Signal parameters
    let signal_frequency = FREQUENCY_0; // Signal frequency (Hz)
    let signal_photon_energy = H_BAR * signal_frequency; // Signal photon energy (J)
    let mean_photon_number = 100.0; // Mean photon number per pulse

    // Time parameters
    let time_window = (-100e-9, 100e-9); // Time window for simulation (s)
    let time_steps = 1000; // Number of time steps
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
            let delta_phase = calculate_probe_phase_shift(num_excited_atoms, probe_detuning);

            delta_phase
        })
        .collect();

    // Output results
    info!(
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

    info!("Simulation data saved to files: {}, {}", filename, data_filename);

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
    const H_BAR: f64 = 1.054571e-34; // Planck's constant over 2π (J·s)
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
    const DIP: f64 = 3.584e-29; // Dipole transition matrix element (C·m)
    const EPSILON_0: f64 = 8.854187817e-12; // Vacuum permittivity (F/m)
    const H_BAR: f64 = 1.054571e-34; // Planck's constant over 2π (J·s)
    const GAMMA: f64 = 1.0 / (26e-9); // Decay rate (s⁻¹)
    const C: f64 = 299792458.0; // Speed of light in vacuum (m/s)
    const WAVELENGTH_0: f64 = 780e-9; // Central wavelength (m)
    const FREQUENCY_0: f64 = C / WAVELENGTH_0; // Central frequency (Hz)

    let delta = probe_detuning;
    let gamma = GAMMA;

    let numerator = num_excited_atoms * DIP.powi(2);
    let denominator = EPSILON_0 * H_BAR * Complex::new(delta, gamma / 2.0); // Denominator is complex

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
