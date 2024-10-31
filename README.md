# Quantum Optics

1. Reference: https://www.physicsforums.com/threads/whats-the-meaning-of-mean-photon-number.947078/
 - Typically, introductory texts to quantum optics do not consider free space cw fields, but light fields inside a perfectly reflecting cavity. There the meaning of the mean photon number is quite clear as no additional energy enters the system and no light leaves the cavity. The mean photon number is just given by the light field inside the cavity. As you incorporate losses to this cavity mode, you will usually create a balance of pump and decay. Some additional light enters the cavity and some leaves the cavity. A steady state will form and the mean photon number of the state inside the cavity is still well defined. The light field outside the cavity, however, is defined by the leakage of photons through the cavity walls. This is given by a leakage rate (some number of photons per second). For the light field inside the cavity, the notion of a mean photon number makes sense as you can measure the total number of photons inside the cavity and this does not depend on integration time. For the light field outside the cavity, such a notion does not make sense as you will always measure a photon flux: The number of photons that pass through a certain detector diameter in a certain amount of time. To get to a total detected photon number you will integrate the photon flux over the detection time. Based on that, you can calculate a mean photon flux or a mean photon number for a certain given integration time.
 - (1) You can calculate a kind of mean photon number, but it is not necessarily related to the concept of a mean photon number that one would use for the state of the field inside the cavity. The difference is most obvious for a single photon state. There, the maximum photon number you can detect is 1 and it will not increase with time. (2) For a single photon state of 1 W, you would need a REALLY high-energy single photon. However, the typical meaning of a single photon state in real devices and experiments is somewhat different. Usually you take a single photon emitter such as a single atom or a quantum dot and place it inside a cavity. You excite it using a pulsed laser and after some time, it will emit a single photon. The cavity has a finite lifetime (say, 1 ns), so the single photon will leave the cavity within this time range and you will have a single photon state outside the cavity. Then you simply excite the system again using the next pulse. A typical time between consecutive pulses often found in experiments is 13 ns. So 13 ns after the first photon another single photon will leave the cavity within a time range of 1 ns. This way, you get a consecutive train of single photons. So here, the power in the light beam is mainly governed by the pulse repetition rate. Accordingly, your question about "mean photon numbers" is meaningless. You can define a mean photon number per arbitrarily chosen integration time (which is a mean count rate) or (more common) for pulsed light fields you can define a mean photon number per pulse, but for free space light fields the notion of a mean photon number does not make much sense without reference to some fixed time range and volume. Mentz114 said: The mean photon number per unit time is given by the Poisson distribution. Only for coherent light fields or classical light fields with a coherence time shorter than the integration time (or of course horrible detector efficiency).
2. Daniela Angulo: https://arxiv.org/pdf/2409.03680
 - Negative time spent in an atom cloud


 ## Parameter Indexes
1. // Experimental parameters
    - 1.1) let pulse_durations = [10e-9, 18e-9, 27e-9, 36e-9]; // Pulse durations (s)
    - 1.2) let optical_depths = [2.0, 3.0, 4.0]; // Optical depths
    - 1.3) let num_atoms = 1e6 as usize; // Number of atoms in the cloud
    - 1.4) let cloud_length = 1e-3; // Length of the atomic cloud(m)
    - 1.5) let cloud_area = PI * (25e-6_f64).powi(2);
    - 1.6) // Cross-sectional area (m²) let atom_density = num_atoms as f64 / (cloud_length * cloud_area); // Atom density (atoms/m³)

2. // Probe parameters
    - 2.1) let probe_detuning = -20e6; // Probe detuning (Hz)
    - 2.2) let probe_frequency = FREQUENCY_0 + probe_detuning; //Probe frequency (Hz)


    