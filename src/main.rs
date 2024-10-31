// main.rs

use std::f64::consts::PI;
use plotters::prelude::*; // For plotting
use log::{info, error};
use env_logger::Env;

// Import structs, constants, and functions from lib.rs
use funhalf::{
    SimulationParams, compute_current, load_quantum_states, QuantumState,
};

fn main() {
    // Initialize the logger with default level "info"
    env_logger::Builder::from_env(Env::default().default_filter_or("info")).init();

    // Define simulation parameters
    let params = SimulationParams {
        U: 2.0e-3,             // Charging energy (eV)
        epsilon_0: 3.0e-3,     // Inter-shell spacing (eV)
        temperature: 0.5,      // Temperature (K)
        eta: 0.55,             // Bias asymmetry
        gamma_l: 0.3556e-6,    // Tunneling rate at left lead (eV)
        gamma_r: 0.0444e-6,    // Tunneling rate at right lead (eV)
        delta_phi: PI / 2.0,   // Phase difference (radians)
        gamma_rel: 0.0,        // Relaxation rate (eV)
        J: 0.01e-3,            // Exchange interaction (eV)
        bias_voltage_range: (-5e-3, 5e-3), // Bias voltage range (V)
        gate_voltage_range: (-20e-3, 20e-3), // Gate voltage range (V)
        bias_steps: 200,       // Number of bias voltage steps
        gate_steps: 200,       // Number of gate voltage steps
        quantum_states_file: "valid_carbon_readings.csv".to_string(), // Path to CSV file
    };

    // Log the start of the simulation
    info!("Starting simulation with parameters: {:?}", params);

    // Load quantum states from CSV file
    let quantum_states = match load_quantum_states(&params.quantum_states_file) {
        Ok(states) => {
            info!("Loaded {} quantum states from {}", states.len(), params.quantum_states_file);
            states
        },
        Err(e) => {
            error!("Failed to load quantum states: {}", e);
            return;
        }
    };

    // Run the simulation
    match run_simulation(params, &quantum_states) {
        Ok(_) => info!("Simulation completed successfully."),
        Err(e) => error!("Simulation failed: {:?}", e),
    }
}

// Function to run the simulation
fn run_simulation(params: SimulationParams, quantum_states: &[QuantumState]) -> Result<(), Box<dyn std::error::Error>> {
    // Unpack parameters
    let bias_voltage_range = params.bias_voltage_range;
    let gate_voltage_range = params.gate_voltage_range;
    let bias_steps = params.bias_steps;
    let gate_steps = params.gate_steps;

    // Generate bias voltage and gate voltage vectors
    let bias_voltages: Vec<f64> = (0..bias_steps)
        .map(|i| bias_voltage_range.0 + i as f64 * (bias_voltage_range.1 - bias_voltage_range.0) / (bias_steps as f64))
        .collect();

    let gate_voltages: Vec<f64> = (0..gate_steps)
        .map(|i| gate_voltage_range.0 + i as f64 * (gate_voltage_range.1 - gate_voltage_range.0) / (gate_steps as f64))
        .collect();

    // Initialize a 2D vector for currents
    let mut currents: Vec<Vec<f64>> = Vec::with_capacity(gate_steps);

    for &V_g in &gate_voltages {
        let current_row: Vec<f64> = bias_voltages
            .iter()
            .map(|&V_b| {
                compute_current(V_b, V_g, &params, quantum_states)
            })
            .collect();
        currents.push(current_row);
    }

    // Output results
    info!("Simulation data calculated.");

    // Plot the results
    let filename = format!("current_map.png");
    let root_area = BitMapBackend::new(&filename, (800, 600)).into_drawing_area();
    root_area.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root_area)
        .caption("Current as a function of Bias Voltage and Gate Voltage", ("sans-serif", 20))
        .margin(20)
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d(
            bias_voltage_range.0..bias_voltage_range.1,
            gate_voltage_range.0..gate_voltage_range.1,
        )?;

    chart
        .configure_mesh()
        .x_desc("Bias Voltage (V)")
        .y_desc("Gate Voltage (V)")
        .draw()?;

    // Prepare data for plotting
    let data: Vec<(f64, f64, f64)> = gate_voltages
        .iter()
        .zip(currents.iter())
        .flat_map(|(&V_g, current_row)| {
            bias_voltages.iter().zip(current_row.iter()).map(move |(&V_b, &I)| (V_b, V_g, I))
        })
        .collect();

    chart.draw_series(
        data.iter().map(|&(x, y, z)| {
            Circle::new((x, y), 2, RGBColor(0, (z * 1e9) as u8, 0).filled())
        }),
    )?;

    // Save data to file
    let data_filename = "current_data.csv";
    save_data(data_filename, &bias_voltages, &gate_voltages, &currents)?;
    info!("Simulation data saved to files: {}, {}", filename, data_filename);

    Ok(())
}

// Function to save data to CSV
fn save_data(
    filename: &str,
    bias_voltages: &[f64],
    gate_voltages: &[f64],
    currents: &[Vec<f64>],
) -> Result<(), std::io::Error> {
    use std::fs::File;
    use std::io::Write;

    let mut file = File::create(filename)?;
    // Write header
    write!(file, "Gate Voltage (V)/Bias Voltage (V)")?;
    for &V_b in bias_voltages {
        write!(file, ",{}", V_b)?;
    }
    writeln!(file)?;

    // Write data
    for (i, &V_g) in gate_voltages.iter().enumerate() {
        write!(file, "{}", V_g)?;
        for &I in &currents[i] {
            write!(file, ",{}", I)?;
        }
        writeln!(file)?;
    }
    Ok(())
}
