import numpy as np
import matplotlib.pyplot as plt
import os

# Create docs folder if not exists
os.makedirs("docs", exist_ok=True)

# =============================================================================
# 0. CONSTANTS AND MATERIAL DATA
# =============================================================================

# Fundamental constants
k0 = 8.9875517923e9         # Coulomb constant
e_charge = 1.60217663e-19   # Elementary charge
m_e = 9.10938356e-31        # Electron mass
c_light = 2.99792458e8      # Speed of light

# Unit conversions
MeV_to_Joule = 1e6 * e_charge
Joule_to_MeV = 1.0 / MeV_to_Joule

# Material parameters
rho_water = 1.0             # g/cm^3
n_water = 3.3428e29         # 1/m^3
I_water = 75.0              # eV

rho_h2 = 8.988e-5           # g/cm^3
n_h2 = 5.3745e25            # 1/m^3
I_h2 = 19.0                 # eV

# Particle masses
m_proton = 1.6726219e-27    # kg
m_alpha = 6.644657e-27      # kg


# =============================================================================
# 1. RELATIVISTIC KINEMATICS AND BETHE–BLOCH
# =============================================================================

def calculate_beta_gamma(energy_MeV: float, mass_kg: float):
    e_kin_J = energy_MeV * MeV_to_Joule
    e_rest_J = mass_kg * c_light**2
    e_total_J = e_kin_J + e_rest_J
    gamma = e_total_J / e_rest_J
    if gamma < 1.0: gamma = 1.0
    beta = np.sqrt(1.0 - 1.0 / gamma**2)
    return beta, gamma

def bethe_bloch_exact(energy_MeV, particle_Z, particle_mass_kg, n_m3, I_eV):
    beta, gamma = calculate_beta_gamma(energy_MeV, particle_mass_kg)
    if beta < 1e-4: return 0.0

    I_Joule = I_eV * e_charge
    prefactor = (4.0 * np.pi * k0**2 * particle_Z**2 * e_charge**4 * n_m3 / (m_e * c_light**2 * beta**2))
    argument = (2.0 * m_e * c_light**2 * beta**2) / (I_Joule * (1.0 - beta**2))

    if argument <= 1.0: return 0.0

    bracket = np.log(argument) - beta**2
    dE_dx_J_m = prefactor * bracket
    dE_dx_MeV_cm = (dE_dx_J_m * Joule_to_MeV) / 100.0
    return dE_dx_MeV_cm

def calculate_beta(energy_MeV, mass_kg):
    beta, _ = calculate_beta_gamma(energy_MeV, mass_kg)
    return beta

def bethe_formula(energy_MeV, particle_Z, particle_mass_kg, n_m3, I_eV):
    beta = calculate_beta(energy_MeV, particle_mass_kg)
    if beta < 1e-5: return 0.0
    I_Joule = I_eV * e_charge
    prefactor = (4.0 * np.pi * k0**2 * particle_Z**2 * e_charge**4 * n_m3 / (m_e * c_light**2 * beta**2))
    argument = (2.0 * m_e * c_light**2 * beta**2) / (I_Joule * (1.0 - beta**2))
    if argument <= 1.0: return 0.0
    bracket = np.log(argument) - beta**2
    return prefactor * bracket


# =============================================================================
# 2. HIGH-PRECISION PHASE-CHANGE SIMULATION
# =============================================================================

def simulate_high_precision(start_energy_MeV, mode, transition_cm):
    current_E = start_energy_MeV
    current_x = 0.0
    positions = [0.0]
    energies = [current_E]
    dedx_values = []

    if mode == 'vaporization':
        current_n, current_I = n_water, I_water
        target_n, target_I = n_h2, I_h2
    else:
        current_n, current_I = n_h2, I_h2
        target_n, target_I = n_water, I_water

    switched = False
    initial_sp = bethe_bloch_exact(current_E, 1, m_proton, current_n, current_I)
    dedx_values.append(initial_sp)

    while current_E > 0.01:
        if (not switched) and (current_x >= transition_cm):
            current_n, current_I = target_n, target_I
            switched = True
            new_sp = bethe_bloch_exact(current_E, 1, m_proton, current_n, current_I)
            positions.append(current_x)
            energies.append(current_E)
            dedx_values.append(new_sp)

        dedx = bethe_bloch_exact(current_E, 1, m_proton, current_n, current_I)
        if dedx <= 0.0: break

        fractional_loss = 0.001
        delta_E = current_E * fractional_loss
        if delta_E < 1e-5: delta_E = 1e-5
        delta_x = delta_E / dedx

        if (not switched) and (current_x + delta_x > transition_cm):
            delta_x = transition_cm - current_x
            if delta_x < 1e-12:
                delta_x = 0.0
                current_x = transition_cm
            else:
                delta_E = delta_x * dedx

        current_x += delta_x
        current_E -= delta_E
        positions.append(current_x)
        energies.append(current_E)
        dedx_values.append(dedx)

    return np.array(positions), np.array(energies), np.array(dedx_values)

# =============================================================================
# 3. RANGE AND STOPPING POWER UTILITIES
# =============================================================================

def get_stopping_power_arrays(energies_MeV, Z, mass_kg, n_m3, I_eV, rho_g_cm3):
    sp = []
    msp = []
    for E in energies_MeV:
        dE_dx_J_m = bethe_formula(E, Z, mass_kg, n_m3, I_eV)
        dE_dx_MeV_cm = (dE_dx_J_m * Joule_to_MeV) / 100.0
        sp.append(dE_dx_MeV_cm)
        msp.append(dE_dx_MeV_cm / rho_g_cm3)
    return np.array(sp), np.array(msp)

def simulate_trajectory(start_E_MeV, Z, mass_kg, n_m3, I_eV):
    distances = [0.0]
    energies = [start_E_MeV]
    current_E = start_E_MeV
    current_x = 0.0
    while current_E > 0.05:
        dE_dx_J_m = bethe_formula(current_E, Z, mass_kg, n_m3, I_eV)
        dE_dx_MeV_cm = (dE_dx_J_m * Joule_to_MeV) / 100.0
        if dE_dx_MeV_cm <= 0.0: break
        dE = current_E * 0.005
        dx = dE / dE_dx_MeV_cm
        current_x += dx
        current_E -= dE
        distances.append(current_x)
        energies.append(current_E)
    return np.array(distances), np.array(energies)

# =============================================================================
# 4. MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    print("Running Bethe-Bloch Energy Loss Simulation...")

    # --- SIMULATION 1: Phase Change ---
    dist_vap, en_vap, sp_vap = simulate_high_precision(200.0, "vaporization", 1e-30)
    dist_cond, en_cond, sp_cond = simulate_high_precision(200.0, "condensation", 500.0)

    # Plot 1: Phase Change Analysis (MAIN VISUAL)
    fig_pc, axs_pc = plt.subplots(2, 2, figsize=(14, 10))
    
    # Vaporization
    ax = axs_pc[0, 0]
    ax.plot(dist_vap, en_vap, 'b', linewidth=1.5)
    ax.axvline(x=0.0, color='r', linestyle='--', label='Phase change')
    ax.set_title("Vaporization: Energy vs Distance")
    ax.set_ylabel("Energy (MeV)")
    ax.legend()
    ax.grid(True, alpha=0.3)

    ax = axs_pc[0, 1]
    ax.plot(dist_vap, sp_vap, 'orange', linewidth=1.5)
    ax.axvline(x=0.0, color='r', linestyle='--')
    ax.set_title("Vaporization: dE/dx vs Distance")
    ax.set_ylabel("Stopping Power (MeV/cm)")
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3)

    # Condensation
    ax = axs_pc[1, 0]
    ax.plot(dist_cond, en_cond, 'g', linewidth=1.5)
    ax.axvline(x=500.0, color='r', linestyle='--', label='Phase change')
    ax.set_title("Condensation: Energy vs Distance")
    ax.set_ylabel("Energy (MeV)")
    ax.set_xlabel("Distance (cm)")
    ax.legend()
    ax.grid(True, alpha=0.3)

    ax = axs_pc[1, 1]
    ax.plot(dist_cond, sp_cond, 'purple', linewidth=1.5)
    ax.axvline(x=500.0, color='r', linestyle='--')
    ax.set_title("Condensation: dE/dx vs Distance")
    ax.set_ylabel("Stopping Power (MeV/cm)")
    ax.set_xlabel("Distance (cm)")
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("docs/fig1_phase_change.png")
    print("Saved: docs/fig1_phase_change.png")
    
    # (Optional: Code calculates stopping power and range curves but we only save Fig 1 for README)
