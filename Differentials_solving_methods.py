from datetime import datetime

from World_constants import Worldconstants as Wc  # Импорт класса, содержащего мировые константы.

from Unresonant_reactions import *
from Resonant_reactions import *
from Decay_reactions import *
from Processes import *
from Element_operations import *
from Neutrino_spectrum import *
from Electron_positron_operations import *
from Temperature_change import *


Max_up_concentration_change = 1e+0

def Curant_diff_solve(Unresonant_reactions, Resonant_reactions, Decay_reactions, Processes, Speed, Time, Smart, Concentrations, Mass_fractions, R_core, T, ro, Average_per_particle_weight, Electron_unresonant_reaction_effect, Positron_unresonant_reaction_effect, Electron_resonant_reaction_effect, Positron_resonant_reaction_effect, Electron_decay_reaction_effect, Positron_decay_reaction_effect, Electron_process_effect, Positron_process_effect, Neutrino_unresonant_reactions_matrix, Neutrino_resonant_reactions_matrix, Neutrino_decay_reactions_matrix, Neutrino_processes_matrix, Elements_burning_time, Burning_speed, Elements, Electrons_concentration, Positrons_concentration, T_surface, R_star, Distance_neutrino_spectrum):
    
    Unresonant_reactions = Unresonant_reactions_reactions_count(Unresonant_reactions, Concentrations, Mass_fractions, ro, T, Average_per_particle_weight)
    Resonant_reactions = Resonant_reactions_reactions_count(Resonant_reactions, T, Elements, Concentrations)
    Decay_reactions = Decay_reactions_count(Decay_reactions, Concentrations)
    Processes = Processes_count(Processes, Concentrations)  # Блок, отвечающий за подсчёт нерезонансных, резонансных реакций, реакций распада и процессов по данным из матрицы.
    
    Elements_burning_time[1] = []
    Burning_speed[1] = []
    for t2 in Concentrations[0]:
        Elements_burning_time[1].append(0)
        Burning_speed[1].append(0)  # Блок, отвечающий за обнуление матриц скосроти и времени выгорания элементов перед новым подсчётом.
            
    for t2 in range(len(Burning_speed[0])):
        for t3 in range(len(Unresonant_reactions[0])):
            for t4 in [Unresonant_reactions[0][t3], Unresonant_reactions[1][t3]]:
                if Burning_speed[0][t2] == t4:
                    Burning_speed[1][t2] -= float(Unresonant_reactions[16][t3])
            for t4 in [Unresonant_reactions[10][t3], Unresonant_reactions[11][t3], Unresonant_reactions[12][t3]]:
                if Burning_speed[0][t2] == t4:
                    Burning_speed[1][t2] += float(Unresonant_reactions[16][t3])
        for t3 in range(len(Resonant_reactions[0])):
            for t4 in [Resonant_reactions[0][t3], Resonant_reactions[1][t3]]:
                if Burning_speed[0][t2] == t4:
                    Burning_speed[1][t2] -= float(Resonant_reactions[25][t3])
            for t4 in [Resonant_reactions[9][t3], Resonant_reactions[10][t3]]:
                if Burning_speed[0][t2] == t4:
                    Burning_speed[1][t2] += float(Resonant_reactions[25][t3])
        for t3 in range(len(Decay_reactions[0])):
            for t4 in [Decay_reactions[0][t3]]:
                if Burning_speed[0][t2] == t4:
                    Burning_speed[1][t2] -= float(Decay_reactions[19][t3])
            for t4 in [Decay_reactions[4][t3], Decay_reactions[5][t3]]:
                if Burning_speed[0][t2] == t4:
                    Burning_speed[1][t2] += float(Decay_reactions[19][t3])
        for t3 in range(len(Processes[0])):
            for t4 in [Processes[0][t3], Processes[1][t3], Processes[2][t3]]:
                if Burning_speed[0][t2] == t4:
                    Burning_speed[1][t2] -= float(Processes[9][t3])
            for t4 in [Processes[3][t3], Processes[4][t3], Processes[5][t3], Processes[6][t3]]:
                if Burning_speed[0][t2] == t4:
                    Burning_speed[1][t2] += float(Processes[9][t3])  # Блок, отвечающий за пересчёт скосротей протекания реакций, в скорости выгорания отдельных элементов в зависимости от их участия и роли (продукт-реагент) в реакции.
                    
    Electrons_burning_speed, Positrons_burning_speed = Electron_positron_burning_speed(Electron_unresonant_reaction_effect, Positron_unresonant_reaction_effect, Unresonant_reactions, Electron_resonant_reaction_effect, Positron_resonant_reaction_effect, Resonant_reactions, Electron_decay_reaction_effect, Positron_decay_reaction_effect, Decay_reactions, Electron_process_effect, Positron_process_effect, Processes)
    # Блок, отвечающий за расчёт скосроти выгорания электронов и позитронов.
    
    for t2 in range(len(Elements_burning_time[0])):
        if Burning_speed[1][t2] == 0 or Concentrations[1][t2] == 0:
            Elements_burning_time[1][t2] = float("+inf")
        elif Burning_speed[1][t2] < 0:
            Elements_burning_time[1][t2] = abs(Concentrations[1][t2] / Burning_speed[1][t2])
        else:
            # Elements_burning_time[1][t2] = abs((min((Concentrations[1][t2] * Max_up_concentration_change), (Max_concentration_of_element(Concentrations[0][t2], Elements, ro) - Concentrations[1][t2]))) /Burning_speed[1][t2])  # Блок, отвечающий за пересчёт скорости выгорания элементов во время их выгорания при текущей скорости выгорания.
            Elements_burning_time[1][t2] = abs((Max_concentration_of_element(Concentrations[0][t2], Elements, ro) - Concentrations[1][t2]) / Burning_speed[1][t2])  # Блок, отвечающий за пересчёт скорости выгорания элементов во время их выгорания при текущей скорости выгорания.

#    print(Concentrations, Burning_speed, Elements_burning_time)
    Electrons_burning_time, Positrons_burning_time = El_pos_burning_time_count(Electrons_burning_speed, Positrons_burning_speed, Electrons_concentration, Positrons_concentration)  # Блок, отвечающий за пересчёт скорости выгорания электронов и позитронов во время их выгорания при текущей скорости выгорания.

    Neutrino_unresonant_reactions_matrix, Neutrino_resonant_reactions_matrix, Neutrino_decay_reactions_matrix, Neutrino_processes_matrix = Neutrino_spectrum_create(Neutrino_unresonant_reactions_matrix, Unresonant_reactions, Neutrino_resonant_reactions_matrix, Resonant_reactions, Neutrino_decay_reactions_matrix, Decay_reactions, Neutrino_processes_matrix, Processes, R_core, Distance_neutrino_spectrum)

    Temperature_change_speed = Temperature_change_speed_count(Unresonant_reactions, Resonant_reactions, Decay_reactions, Processes, T_surface, R_star, R_core, Concentrations, Neutrino_unresonant_reactions_matrix, Neutrino_resonant_reactions_matrix, Neutrino_decay_reactions_matrix, Neutrino_processes_matrix, Distance_neutrino_spectrum)
    
    Temperature_change_time = Temperature_change_time_count(T, Temperature_change_speed)
            
    dt = min(min(Elements_burning_time[1]), Temperature_change_time, Electrons_burning_time) * Speed # шаг по времени, с которым будем интегрировать
    if Time > 1e-9:
        dt = min(dt, Time * 0.05)
    else:
        dt = min(dt, 1e-10)
    
    return Burning_speed, Speed, Electrons_burning_speed, Positrons_burning_speed, dt, Elements_burning_time, Temperature_change_speed, Temperature_change_time


def Runge_Kutta_4(Unresonant_reactions, Resonant_reactions, Decay_reactions, Processes, Speed, Time, Concentrations, Mass_fractions, R_core, T, ro, Average_per_particle_weight, Electron_unresonant_reaction_effect, Positron_unresonant_reaction_effect, Electron_resonant_reaction_effect, Positron_resonant_reaction_effect, Electron_decay_reaction_effect, Positron_decay_reaction_effect, Electron_process_effect, Positron_process_effect, Neutrino_unresonant_reactions_matrix, Neutrino_resonant_reactions_matrix, Neutrino_decay_reactions_matrix, Neutrino_processes_matrix, Elements_burning_time, Burning_speed, Elements, Electrons_concentration, Positrons_concentration, T_surface, R_star, Temperature_change_conf, Distance_neutrino_spectrum):

    Unresonant_reactions = Unresonant_reactions_reactions_count(Unresonant_reactions, Concentrations, Mass_fractions,
                                                                ro, T, Average_per_particle_weight)
    Resonant_reactions = Resonant_reactions_reactions_count(Resonant_reactions, T, Elements, Concentrations)
    Decay_reactions = Decay_reactions_count(Decay_reactions, Concentrations)
    Processes = Processes_count(Processes,
                                Concentrations)  # Блок, отвечающий за подсчёт нерезонансных, резонансных реакций, реакций распада и процессов по данным из матрицы.


    Neutrino_unresonant_reactions_matrix, Neutrino_resonant_reactions_matrix, Neutrino_decay_reactions_matrix, Neutrino_processes_matrix = Neutrino_spectrum_create(
        Neutrino_unresonant_reactions_matrix, Unresonant_reactions, Neutrino_resonant_reactions_matrix,
        Resonant_reactions, Neutrino_decay_reactions_matrix, Decay_reactions, Neutrino_processes_matrix, Processes,
        R_core, Distance_neutrino_spectrum)

    Burning_speed_1 = Burning_speed[:]

    Electrons_burning_speed_1, Positrons_burning_speed_1 = Electron_positron_burning_speed(
        Electron_unresonant_reaction_effect, Positron_unresonant_reaction_effect, Unresonant_reactions,
        Electron_resonant_reaction_effect, Positron_resonant_reaction_effect, Resonant_reactions,
        Electron_decay_reaction_effect, Positron_decay_reaction_effect, Decay_reactions, Electron_process_effect,
        Positron_process_effect, Processes)
    # Блок, отвечающий за расчёт скосроти выгорания электронов и позитронов.

    Temperature_change_speed_1 = Temperature_change_speed_count(Unresonant_reactions, Resonant_reactions, Decay_reactions,
                                                              Processes, T_surface, R_star, R_core, Concentrations,
                                                              Neutrino_unresonant_reactions_matrix,
                                                              Neutrino_resonant_reactions_matrix,
                                                              Neutrino_decay_reactions_matrix,
                                                              Neutrino_processes_matrix, Distance_neutrino_spectrum)

    Burning_speed_1[1] = []
    for t2 in Concentrations[0]:
        Burning_speed_1[1].append(0)  # Блок, отвечающий за обнуление матриц скосроти и времени выгорания элементов перед новым подсчётом.

    for t2 in range(len(Burning_speed_1[0])):
        for t3 in range(len(Unresonant_reactions[0])):
            for t4 in [Unresonant_reactions[0][t3], Unresonant_reactions[1][t3]]:
                if Burning_speed_1[0][t2] == t4:
                    Burning_speed_1[1][t2] -= float(Unresonant_reactions[16][t3])
            for t4 in [Unresonant_reactions[10][t3], Unresonant_reactions[11][t3], Unresonant_reactions[12][t3]]:
                if Burning_speed_1[0][t2] == t4:
                    Burning_speed_1[1][t2] += float(Unresonant_reactions[16][t3])
        for t3 in range(len(Resonant_reactions[0])):
            for t4 in [Resonant_reactions[0][t3], Resonant_reactions[1][t3]]:
                if Burning_speed_1[0][t2] == t4:
                    Burning_speed_1[1][t2] -= float(Resonant_reactions[25][t3])
            for t4 in [Resonant_reactions[9][t3], Resonant_reactions[10][t3]]:
                if Burning_speed_1[0][t2] == t4:
                    Burning_speed_1[1][t2] += float(Resonant_reactions[25][t3])
        for t3 in range(len(Decay_reactions[0])):
            for t4 in [Decay_reactions[0][t3]]:
                if Burning_speed_1[0][t2] == t4:
                    Burning_speed_1[1][t2] -= float(Decay_reactions[19][t3])
            for t4 in [Decay_reactions[4][t3], Decay_reactions[5][t3]]:
                if Burning_speed_1[0][t2] == t4:
                    Burning_speed_1[1][t2] += float(Decay_reactions[19][t3])
        for t3 in range(len(Processes[0])):
            for t4 in [Processes[0][t3], Processes[1][t3], Processes[2][t3]]:
                if Burning_speed_1[0][t2] == t4:
                    Burning_speed_1[1][t2] -= float(Processes[9][t3])
            for t4 in [Processes[3][t3], Processes[4][t3], Processes[5][t3], Processes[6][t3]]:
                if Burning_speed_1[0][t2] == t4:
                    Burning_speed_1[1][t2] += float(Processes[9][t3])  # Блок, отвечающий за пересчёт скосротей протекания реакций, в скорости выгорания отдельных элементов в зависимости от их участия и роли (продукт-реагент) в реакции.

    for t2 in range(len(Elements_burning_time[0])):
        if Burning_speed_1[1][t2] == 0 or Concentrations[1][t2] == 0:
            Elements_burning_time[1][t2] = float("+inf")
        elif Burning_speed_1[1][t2] < 0:
            Elements_burning_time[1][t2] = abs(Concentrations[1][t2] / Burning_speed_1[1][t2])
        else:
            Elements_burning_time[1][t2] = abs(
                (min((Concentrations[1][t2] * Max_up_concentration_change), (Max_concentration_of_element(Concentrations[0][t2], Elements, ro) - Concentrations[1][t2]))) /
                Burning_speed_1[1][t2])  # Блок, отвечающий за пересчёт скорости выгорания элементов во время их выгорания при текущей скорости выгорания.
            # Elements_burning_time[1][t2] = abs((Max_concentration_of_element(Concentrations[0][t2], Elements, ro) - Concentrations[1][t2]) / Burning_speed[1][t2])  # Блок, отвечающий за пересчёт скорости выгорания элементов во время их выгорания при текущей скорости выгорания.


    #    print(Concentrations, Burning_speed, Elements_burning_time)
    Electrons_burning_time, Positrons_burning_time = El_pos_burning_time_count(Electrons_burning_speed_1,
                                                                               Positrons_burning_speed_1,
                                                                               Electrons_concentration,
                                                                               Positrons_concentration)  # Блок, отвечающий за пересчёт скорости выгорания электронов и позитронов во время их выгорания при текущей скорости выгорания.

    Temperature_change_time = Temperature_change_time_count(T, Temperature_change_speed_1)

    dt = min(min(Elements_burning_time[1]), Temperature_change_time, Electrons_burning_time) * Speed  # шаг по времени, с которым будем интегрировать
    # if Time > 1e-9:
    #     dt = min(dt, Time * 0.01)

    '''
    После определения шага дифференцирования можно применить сам метод Рунге-Кутты 4 порядка.  
    '''

    Concentrations_1 = Concentrations[:]
    Mass_fractions_1 = Mass_fractions[:]
    # Вычисляем скорости убывания элементов:
    for t2 in range(len(Concentrations_1[0])):
        Concentrations_1[1][t2] += Burning_speed_1[1][t2] * dt / 2

    for t2 in range(len(Mass_fractions_1[0])):  # пересчёт концентраций в массовые доли.
        Mass_fractions_1[1][t2] = Concentrations_1[1][t2] * Wc.M_nuc * Mass_of_element(Mass_fractions_1[0][t2], Elements) / ro

    Electrons_concentration_1 = Electrons_concentration + Electrons_burning_speed_1 * dt / 2
    Positrons_concentration_1 = Positrons_concentration + Positrons_burning_speed_1 * dt / 2  # Расчёт изменения концетрации электронов и позитронов.
    Average_per_particle_weight = Average_per_particle_weight_count(Electrons_concentration, Positrons_concentration,
                                                                    Elements,
                                                                    Concentrations_1)  # Рассчёт среднего молекулярного веса (вся масса в эдинице объёма, делённая на концентрацию всех частиц в объёме).

    if Temperature_change_conf:
        T_1 = T + Temperature_change_speed_1 * dt / 2
    else:
        T_1 = T  # Блок, отвечающий за расчёт изменения температуры (если оно включено).

    '''
    Вторая итерация метода.
    '''

    Unresonant_reactions = Unresonant_reactions_reactions_count(Unresonant_reactions, Concentrations_1, Mass_fractions_1,
                                                                ro, T_1, Average_per_particle_weight)
    Resonant_reactions = Resonant_reactions_reactions_count(Resonant_reactions, T_1, Elements, Concentrations_1)
    Decay_reactions = Decay_reactions_count(Decay_reactions, Concentrations_1)
    Processes = Processes_count(Processes,
                                Concentrations_1)  # Блок, отвечающий за подсчёт нерезонансных, резонансных реакций, реакций распада и процессов по данным из матрицы.

    Burning_speed_2 = Burning_speed_1[:]

    Electrons_burning_speed_2, Positrons_burning_speed_2 = Electron_positron_burning_speed(
        Electron_unresonant_reaction_effect, Positron_unresonant_reaction_effect, Unresonant_reactions,
        Electron_resonant_reaction_effect, Positron_resonant_reaction_effect, Resonant_reactions,
        Electron_decay_reaction_effect, Positron_decay_reaction_effect, Decay_reactions, Electron_process_effect,
        Positron_process_effect, Processes)
    # Блок, отвечающий за расчёт скосроти выгорания электронов и позитронов.

    Temperature_change_speed_2 = Temperature_change_speed_count(Unresonant_reactions, Resonant_reactions,
                                                                Decay_reactions,
                                                                Processes, T_surface, R_star, R_core, Concentrations_1,
                                                                Neutrino_unresonant_reactions_matrix,
                                                                Neutrino_resonant_reactions_matrix,
                                                                Neutrino_decay_reactions_matrix,
                                                                Neutrino_processes_matrix, Distance_neutrino_spectrum)

    Burning_speed_2[1] = []
    for t2 in Concentrations_1[0]:
        Burning_speed_2[1].append(0)  # Блок, отвечающий за обнуление матриц скосроти и времени выгорания элементов перед новым подсчётом.

    for t2 in range(len(Burning_speed_2[0])):
        for t3 in range(len(Unresonant_reactions[0])):
            for t4 in [Unresonant_reactions[0][t3], Unresonant_reactions[1][t3]]:
                if Burning_speed_2[0][t2] == t4:
                    Burning_speed_2[1][t2] -= float(Unresonant_reactions[16][t3])
            for t4 in [Unresonant_reactions[10][t3], Unresonant_reactions[11][t3], Unresonant_reactions[12][t3]]:
                if Burning_speed_2[0][t2] == t4:
                    Burning_speed_2[1][t2] += float(Unresonant_reactions[16][t3])
        for t3 in range(len(Resonant_reactions[0])):
            for t4 in [Resonant_reactions[0][t3], Resonant_reactions[1][t3]]:
                if Burning_speed_2[0][t2] == t4:
                    Burning_speed_2[1][t2] -= float(Resonant_reactions[25][t3])
            for t4 in [Resonant_reactions[9][t3], Resonant_reactions[10][t3]]:
                if Burning_speed_2[0][t2] == t4:
                    Burning_speed_2[1][t2] += float(Resonant_reactions[25][t3])
        for t3 in range(len(Decay_reactions[0])):
            for t4 in [Decay_reactions[0][t3]]:
                if Burning_speed_2[0][t2] == t4:
                    Burning_speed_2[1][t2] -= float(Decay_reactions[19][t3])
            for t4 in [Decay_reactions[4][t3], Decay_reactions[5][t3]]:
                if Burning_speed_2[0][t2] == t4:
                    Burning_speed_2[1][t2] += float(Decay_reactions[19][t3])
        for t3 in range(len(Processes[0])):
            for t4 in [Processes[0][t3], Processes[1][t3], Processes[2][t3]]:
                if Burning_speed_2[0][t2] == t4:
                    Burning_speed_2[1][t2] -= float(Processes[9][t3])
            for t4 in [Processes[3][t3], Processes[4][t3], Processes[5][t3], Processes[6][t3]]:
                if Burning_speed_2[0][t2] == t4:
                    Burning_speed_2[1][t2] += float(Processes[9][t3])  # Блок, отвечающий за пересчёт скосротей протекания реакций, в скорости выгорания отдельных элементов в зависимости от их участия и роли (продукт-реагент) в реакции.

    Concentrations_2 = Concentrations[:]
    Mass_fractions_2 = Mass_fractions[:]
    # Вычисляем скорости убывания элементов:
    for t2 in range(len(Concentrations_1[0])):
        Concentrations_2[1][t2] += Burning_speed_2[1][t2] * dt / 2

    for t2 in range(len(Mass_fractions_2[0])):  # пересчёт концентраций в массовые доли.
        Mass_fractions_2[1][t2] = Concentrations_2[1][t2] * Wc.M_nuc * Mass_of_element(Mass_fractions_2[0][t2], Elements) / ro

    Electrons_concentration_2 = Electrons_concentration + Electrons_burning_speed_2 * dt / 2
    Positrons_concentration_2 = Positrons_concentration + Positrons_burning_speed_2 * dt / 2  # Расчёт изменения концетрации электронов и позитронов.
    Average_per_particle_weight = Average_per_particle_weight_count(Electrons_concentration, Positrons_concentration,
                                                                    Elements,
                                                                    Concentrations_2)  # Рассчёт среднего молекулярного веса (вся масса в эдинице объёма, делённая на концентрацию всех частиц в объёме).

    if Temperature_change_conf:
        T_2 = T + Temperature_change_speed_2 * dt / 2
    else:
        T_2 = T  # Блок, отвечающий за расчёт изменения температуры (если оно включено).

    '''
    Третья итерация метода.
    '''

    Unresonant_reactions = Unresonant_reactions_reactions_count(Unresonant_reactions, Concentrations_2,
                                                                Mass_fractions_2,
                                                                ro, T_2, Average_per_particle_weight)
    Resonant_reactions = Resonant_reactions_reactions_count(Resonant_reactions, T_2, Elements, Concentrations_2)
    Decay_reactions = Decay_reactions_count(Decay_reactions, Concentrations_2)
    Processes = Processes_count(Processes,
                                Concentrations_2)  # Блок, отвечающий за подсчёт нерезонансных, резонансных реакций, реакций распада и процессов по данным из матрицы.

    Burning_speed_3 = Burning_speed_2[:]

    Electrons_burning_speed_3, Positrons_burning_speed_3 = Electron_positron_burning_speed(
        Electron_unresonant_reaction_effect, Positron_unresonant_reaction_effect, Unresonant_reactions,
        Electron_resonant_reaction_effect, Positron_resonant_reaction_effect, Resonant_reactions,
        Electron_decay_reaction_effect, Positron_decay_reaction_effect, Decay_reactions, Electron_process_effect,
        Positron_process_effect, Processes)
    # Блок, отвечающий за расчёт скосроти выгорания электронов и позитронов.

    Temperature_change_speed_3 = Temperature_change_speed_count(Unresonant_reactions, Resonant_reactions,
                                                                Decay_reactions,
                                                                Processes, T_surface, R_star, R_core, Concentrations_2,
                                                                Neutrino_unresonant_reactions_matrix,
                                                                Neutrino_resonant_reactions_matrix,
                                                                Neutrino_decay_reactions_matrix,
                                                                Neutrino_processes_matrix, Distance_neutrino_spectrum)

    Burning_speed_3[1] = []
    for t2 in Concentrations_2[0]:
        Burning_speed_3[1].append(0)  # Блок, отвечающий за обнуление матриц скосроти и времени выгорания элементов перед новым подсчётом.

    for t2 in range(len(Burning_speed_3[0])):
        for t3 in range(len(Unresonant_reactions[0])):
            for t4 in [Unresonant_reactions[0][t3], Unresonant_reactions[1][t3]]:
                if Burning_speed_3[0][t2] == t4:
                    Burning_speed_3[1][t2] -= float(Unresonant_reactions[16][t3])
            for t4 in [Unresonant_reactions[10][t3], Unresonant_reactions[11][t3], Unresonant_reactions[12][t3]]:
                if Burning_speed_3[0][t2] == t4:
                    Burning_speed_3[1][t2] += float(Unresonant_reactions[16][t3])
        for t3 in range(len(Resonant_reactions[0])):
            for t4 in [Resonant_reactions[0][t3], Resonant_reactions[1][t3]]:
                if Burning_speed_3[0][t2] == t4:
                    Burning_speed_3[1][t2] -= float(Resonant_reactions[25][t3])
            for t4 in [Resonant_reactions[9][t3], Resonant_reactions[10][t3]]:
                if Burning_speed_3[0][t2] == t4:
                    Burning_speed_3[1][t2] += float(Resonant_reactions[25][t3])
        for t3 in range(len(Decay_reactions[0])):
            for t4 in [Decay_reactions[0][t3]]:
                if Burning_speed_3[0][t2] == t4:
                    Burning_speed_3[1][t2] -= float(Decay_reactions[19][t3])
            for t4 in [Decay_reactions[4][t3], Decay_reactions[5][t3]]:
                if Burning_speed_3[0][t2] == t4:
                    Burning_speed_3[1][t2] += float(Decay_reactions[19][t3])
        for t3 in range(len(Processes[0])):
            for t4 in [Processes[0][t3], Processes[1][t3], Processes[2][t3]]:
                if Burning_speed_3[0][t2] == t4:
                    Burning_speed_3[1][t2] -= float(Processes[9][t3])
            for t4 in [Processes[3][t3], Processes[4][t3], Processes[5][t3], Processes[6][t3]]:
                if Burning_speed_3[0][t2] == t4:
                    Burning_speed_3[1][t2] += float(Processes[9][t3])  # Блок, отвечающий за пересчёт скосротей протекания реакций, в скорости выгорания отдельных элементов в зависимости от их участия и роли (продукт-реагент) в реакции.

    Concentrations_3 = Concentrations[:]
    Mass_fractions_3 = Mass_fractions[:]
    # Вычисляем скорости убывания элементов:
    for t2 in range(len(Concentrations_2[0])):
        Concentrations_3[1][t2] += Burning_speed_3[1][t2] * dt

    for t2 in range(len(Mass_fractions_3[0])):  # пересчёт концентраций в массовые доли.
        Mass_fractions_3[1][t2] = Concentrations_3[1][t2] * Wc.M_nuc * Mass_of_element(Mass_fractions_3[0][t2],
                                                                                       Elements) / ro

    Electrons_concentration_3 = Electrons_concentration + Electrons_burning_speed_3 * dt
    Positrons_concentration_3 = Positrons_concentration + Positrons_burning_speed_3 * dt  # Расчёт изменения концетрации электронов и позитронов.
    Average_per_particle_weight = Average_per_particle_weight_count(Electrons_concentration, Positrons_concentration,
                                                                    Elements,
                                                                    Concentrations_3)  # Рассчёт среднего молекулярного веса (вся масса в эдинице объёма, делённая на концентрацию всех частиц в объёме).

    if Temperature_change_conf:
        T_3 = T + Temperature_change_speed_3 * dt
    else:
        T_3 = T  # Блок, отвечающий за расчёт изменения температуры (если оно включено).

    '''
    Четвёртая итерация метода.
    '''

    Unresonant_reactions = Unresonant_reactions_reactions_count(Unresonant_reactions, Concentrations_3,
                                                                Mass_fractions_3,
                                                                ro, T_3, Average_per_particle_weight)
    Resonant_reactions = Resonant_reactions_reactions_count(Resonant_reactions, T_3, Elements, Concentrations_3)
    Decay_reactions = Decay_reactions_count(Decay_reactions, Concentrations_3)
    Processes = Processes_count(Processes,
                                Concentrations_3)  # Блок, отвечающий за подсчёт нерезонансных, резонансных реакций, реакций распада и процессов по данным из матрицы.

    Burning_speed_4 = Burning_speed_3[:]

    Electrons_burning_speed_4, Positrons_burning_speed_4 = Electron_positron_burning_speed(
        Electron_unresonant_reaction_effect, Positron_unresonant_reaction_effect, Unresonant_reactions,
        Electron_resonant_reaction_effect, Positron_resonant_reaction_effect, Resonant_reactions,
        Electron_decay_reaction_effect, Positron_decay_reaction_effect, Decay_reactions, Electron_process_effect,
        Positron_process_effect, Processes)
    # Блок, отвечающий за расчёт скосроти выгорания электронов и позитронов.

    Temperature_change_speed_4 = Temperature_change_speed_count(Unresonant_reactions, Resonant_reactions,
                                                                Decay_reactions,
                                                                Processes, T_surface, R_star, R_core, Concentrations_3,
                                                                Neutrino_unresonant_reactions_matrix,
                                                                Neutrino_resonant_reactions_matrix,
                                                                Neutrino_decay_reactions_matrix,
                                                                Neutrino_processes_matrix, Distance_neutrino_spectrum)

    Burning_speed_4[1] = []
    for t2 in Concentrations_3[0]:
        Burning_speed_4[1].append(0)  # Блок, отвечающий за обнуление матриц скосроти и времени выгорания элементов перед новым подсчётом.

    for t2 in range(len(Burning_speed_4[0])):
        for t3 in range(len(Unresonant_reactions[0])):
            for t4 in [Unresonant_reactions[0][t3], Unresonant_reactions[1][t3]]:
                if Burning_speed_4[0][t2] == t4:
                    Burning_speed_4[1][t2] -= float(Unresonant_reactions[16][t3])
            for t4 in [Unresonant_reactions[10][t3], Unresonant_reactions[11][t3], Unresonant_reactions[12][t3]]:
                if Burning_speed_4[0][t2] == t4:
                    Burning_speed_4[1][t2] += float(Unresonant_reactions[16][t3])
        for t3 in range(len(Resonant_reactions[0])):
            for t4 in [Resonant_reactions[0][t3], Resonant_reactions[1][t3]]:
                if Burning_speed_4[0][t2] == t4:
                    Burning_speed_4[1][t2] -= float(Resonant_reactions[25][t3])
            for t4 in [Resonant_reactions[9][t3], Resonant_reactions[10][t3]]:
                if Burning_speed_4[0][t2] == t4:
                    Burning_speed_4[1][t2] += float(Resonant_reactions[25][t3])
        for t3 in range(len(Decay_reactions[0])):
            for t4 in [Decay_reactions[0][t3]]:
                if Burning_speed_4[0][t2] == t4:
                    Burning_speed_4[1][t2] -= float(Decay_reactions[19][t3])
            for t4 in [Decay_reactions[4][t3], Decay_reactions[5][t3]]:
                if Burning_speed_4[0][t2] == t4:
                    Burning_speed_4[1][t2] += float(Decay_reactions[19][t3])
        for t3 in range(len(Processes[0])):
            for t4 in [Processes[0][t3], Processes[1][t3], Processes[2][t3]]:
                if Burning_speed_4[0][t2] == t4:
                    Burning_speed_4[1][t2] -= float(Processes[9][t3])
            for t4 in [Processes[3][t3], Processes[4][t3], Processes[5][t3], Processes[6][t3]]:
                if Burning_speed_4[0][t2] == t4:
                    Burning_speed_4[1][t2] += float(Processes[9][t3])  # Блок, отвечающий за пересчёт скосротей протекания реакций, в скорости выгорания отдельных элементов в зависимости от их участия и роли (продукт-реагент) в реакции.

    '''
    Вычисление финальных производных и значений.
    '''

    Elements_burning_time[1] = []
    Burning_speed[1] = []
    for t2 in Concentrations[0]:
        Elements_burning_time[1].append(0)
        Burning_speed[1].append(0)  # Блок, отвечающий за обнуление матриц скосроти и времени выгорания элементов перед новым подсчётом.

    Burning_speed[1] = []
    for t2 in range(len(Concentrations[0])):
        Burning_speed[1].append((Burning_speed_1[1][t2] + 2 * Burning_speed_2[1][t2] + 2 * Burning_speed_3[1][t2] + Burning_speed_4[1][t2]) / 6)

    Electrons_burning_speed = (Electrons_burning_speed_1 + 2 * Electrons_burning_speed_2 + 2 * Electrons_burning_speed_3 + Electrons_burning_speed_4) / 6
    Positrons_burning_speed = (Positrons_burning_speed_1 + 2 * Positrons_burning_speed_2 + 2 * Positrons_burning_speed_3 + Positrons_burning_speed_4) / 6

    Temperature_change_speed = (Temperature_change_speed_1 + 2 * Temperature_change_speed_2 + 2 * Temperature_change_speed_3 + Temperature_change_speed_4) / 6

    for t2 in range(len(Elements_burning_time[0])):
        if Burning_speed[1][t2] == 0 or Concentrations[1][t2] == 0:
            Elements_burning_time[1][t2] = float("+inf")
        elif Burning_speed[1][t2] < 0:
            Elements_burning_time[1][t2] = abs(Concentrations[1][t2] / Burning_speed[1][t2])
        else:
            Elements_burning_time[1][t2] = abs(
                (min((Concentrations[1][t2] * Max_up_concentration_change), (Max_concentration_of_element(Concentrations[0][t2], Elements, ro) - Concentrations[1][t2]))) /
                Burning_speed[1][t2])  # Блок, отвечающий за пересчёт скорости выгорания элементов во время их выгорания при текущей скорости выгорания.
            # Elements_burning_time[1][t2] = abs((Max_concentration_of_element(Concentrations[0][t2], Elements, ro) - Concentrations[1][t2]) / Burning_speed[1][t2])  # Блок, отвечающий за пересчёт скорости выгорания элементов во время их выгорания при текущей скорости выгорания.


    #    print(Concentrations, Burning_speed, Elements_burning_time)
    Electrons_burning_time, Positrons_burning_time = El_pos_burning_time_count(Electrons_burning_speed,
                                                                               Positrons_burning_speed,
                                                                               Electrons_concentration,
                                                                               Positrons_concentration)  # Блок, отвечающий за пересчёт скорости выгорания электронов и позитронов во время их выгорания при текущей скорости выгорания.

    Temperature_change_time = Temperature_change_time_count(T, Temperature_change_speed)

    # dt = min(min(Elements_burning_time[1]), Temperature_change_time, Electrons_burning_time) * Speed  # шаг по времени, с которым будем интегрировать

    return Burning_speed, Speed, Electrons_burning_speed, Positrons_burning_speed, dt, Elements_burning_time, Temperature_change_speed, Temperature_change_time


def Maxwells_energy_distribution_count(Pi, e, k, T):  # На данный момент функция нигде не используется напрямую.
    # По оси x - энергия, по оси y - часть частиц с такой энергией.
    E_min = 0
    E_max = 100 * k * T
    Point_of_calculation = 100
    E_step = (E_max - E_min) / (Point_of_calculation - 1)
    
    E_values = []
    f_E_values = []
    for t1 in range(Point_of_calculation):
        E = E_min + t1 * E_step
        E_values.append(E)
        f_E_values.append(2 * Pi * E ** 0.5 * e ** (-E / (k * T)) / (Pi * k * T) ** 1.5)
    
    # Нормирование графика методом трапеций (численное интегрирование).
    Norm_factor = 0
    for t1 in range(1, len(E_values) - 1):
        Norm_factor += (E_values[t1] - E_values[t1 - 1]) * (f_E_values[t1] + f_E_values[t1 + 1]) / 2
    Maxwells_energy_distribution = [f_E / Norm_factor for f_E in f_E_values] # Распределение именно по энергии одинаково для всех элементов, дальше на его основе для каждой реакции можно получить распределение по столкновениям на основе концентраций.
    
    return Maxwells_energy_distribution
    
    
