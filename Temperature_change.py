from World_constants import Worldconstants as Wc  # Импорт класса, содержащего мировые константы.

from Neutrino_spectrum import *

def Temperature_change_speed_count(Unresonant_reactions, Resonant_reactions, Decay_reactions, Processes, T_surface, R_star, R_core, Concentrations, Neutrino_unresonant_reactions_matrix, Neutrino_resonant_reactions_matrix, Neutrino_decay_reactions_matrix, Neutrino_processes_matrix, Distance_neutrino_spectrum):
    Total_energy_impact = 0  # Функция, считающая изменение температуры ядра звезды по модели одноатомного ИГ (хоть плазма и состоит из ядер, но это плазма, поэтому гамма не 5/3, а 4/3. Значит i=6, а не i=3.)
    Moles_per_unit = 0
    for t in range(len(Unresonant_reactions[0])):
        Total_energy_impact += Unresonant_reactions[16][t] * Unresonant_reactions[21][t]
    for t in range(len(Resonant_reactions[0])):
        Total_energy_impact += Resonant_reactions[25][t] * Resonant_reactions[21][t]
    for t in range(len(Decay_reactions[0])):
        Total_energy_impact += Decay_reactions[19][t] * Decay_reactions[16][t]
    for t in range(len(Processes[0])):
        Total_energy_impact += Processes[9][t] * Processes[7][t]  # Блок, отвечающий за рассчёт общего энергетического эффекта всех реакций в ядре звезды в SAДжоулях.
    Total_energy_impact *= (4 / 3) * Wc.Pi * R_core ** 3
    Total_energy_lose = Wc.St_Bol * 4 * Wc.Pi * R_star ** 2 * T_surface ** 4  # Рассчёт суммарных потерь энергии звездой при свечении по модели абсолютно чёрного тела по формуле Стефана-Больцмана.

    Neutrino_cooling = True
    if Neutrino_cooling:
        Neutrino_unresonant_reactions_matrix, Neutrino_resonant_reactions_matrix, Neutrino_decay_reactions_matrix, Neutrino_processes_matrix = Neutrino_spectrum_create(Neutrino_unresonant_reactions_matrix, Unresonant_reactions, Neutrino_resonant_reactions_matrix, Resonant_reactions, Neutrino_decay_reactions_matrix, Decay_reactions, Neutrino_processes_matrix, Processes, R_core, Distance_neutrino_spectrum)
        Total_energy_lose += Neutrino_total_energy_count(Neutrino_unresonant_reactions_matrix, Neutrino_resonant_reactions_matrix, Neutrino_decay_reactions_matrix, Neutrino_processes_matrix) * (4 / 3) * Wc.Pi * R_core ** 3

    Energy_change = Total_energy_impact - Total_energy_lose
    C_v = (Wc.i / 2) * Wc.R  # Так как давление в центре может меняться при жизни звезды, то возьму изохорическую теплоёмкость, так как объём при нахождении звезды на ГП не меняется.
    for t1 in range(len(Concentrations[1])):  # Надеюсь, что то, что частицы имеют разную массу, не повлияет на точность расчётов.
        Moles_per_unit += Concentrations[1][t1] / Wc.N_a
    Temperature_change = Energy_change / (C_v * Moles_per_unit * (4 / 3) * Wc.Pi * R_core ** 3)
    return Temperature_change


def Temperature_change_time_count(T, Temperature_change_speed):
    if Temperature_change_speed > 0:
        Temperature_change_time = float("+inf")
    elif T == 0 or Temperature_change_speed == 0:
            Temperature_change_time = float("+inf")
    else:
        Temperature_change_time = abs(T / Temperature_change_speed)
    return Temperature_change_time


def El_pos_annihilation_temperature_change(Annihilation_energy, Concentrations, R_core):
    Moles_per_unit = 0
    C_v = (Wc.i / 2) * Wc.R  # Так как давление в центре может меняться при жизни звезды, то возьму изохорическую теплоёмкость, так как объём при нахождении звезды на ГП не меняется.
    for t1 in range(len(Concentrations[1])):  # Надеюсь, что то, что частицы имеют разную массу, не повлияет на точность расчётов.
        Moles_per_unit += Concentrations[1][t1] / Wc.N_a
    Temperature_change = Annihilation_energy / (C_v * Moles_per_unit * (4 / 3) * Wc.Pi * R_core ** 3)
    return Temperature_change
