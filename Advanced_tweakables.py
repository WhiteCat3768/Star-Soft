class Advancedtweakbles():  # Класс, содержащий специфические параметры симуляции.

    Neutrino_spectrum_time = 1e-200  # Point in time in years from beginning of simulation at which you want to see a neutrino spectrum graph.
    Distance_neutrino_spectrum = 1.5e+8  # Distance from the center of the stellar core to the imagenary neutrino detector in kilometes (should be at least 8-10 bigger the the raduim of the core itself for accurate results).
    Max_up_concentration_change = 1e+0  # Maximum relative concentration increase in one step of the simulation.
    Max_time_step = 1e+6  # Single time step limit in years.
    Max_relative_time_step_change = 0.05  # Limits single time step to a fraction to time that had passed from the beginning of the simulation.
    Resonant_reactions_relative_temperature_update = 0.01  # A maximum relative change of temperature at which resonant reactions cross-sections have been updated last time, for them to be updated again.
    Speed_of_light_velocity_limit = 0.01  # A fraction of speed of light up to which velocities will be accounted for while calculating integrals to update resonant reactions cross-sections, should be kept below 0.05 - 0.1 to avoid relativistic effects, however must be increased while working with ultra-hot plasma.
    Resonance_width_energy_range = 500  # A parameter that determines the number of full resonance widths that make up half the width of the energy integration range when calculating cross-sections of resonant reactions.
    Annihilation = True  # Allows for electrons and positrons to annihilate releasing energy (its assumed to be instant), dont know why you may need to turn it off.
    Neutrino_cooling = True  # Allows neutrino to actually take energy from stellar core when they are emmited, slightly affects performance for a slight accuracy improvement.
