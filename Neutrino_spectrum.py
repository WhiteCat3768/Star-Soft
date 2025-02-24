import matplotlib.pyplot as plt  # Библиотека для построения графиков.
import numpy as np  # Использована для частичного упрощения вычислений.

from World_constants import Worldconstants as Wc  # Импорт класса, содержащего мировые константы.

def Is_neutrino(Components):
    for t1 in Components:
        if t1 == 'v':
            return True
    return False

def Neutrino_flux(Distance, R_core, R):
        V = (4 / 3) * np.pi * R_core ** 3
        S = 4 * np.pi * (Distance * 1e+5) ** 2
        return V * R / S


def Neutrino_energy(Energy, E_max):
    X = Energy / E_max
    return X ** 2 * (1 - 2 * X + X ** 2) * (1 - X ** 2) ** 0.5


def Neutrino_total_energy_count(Neutrino_unresonant_reactions_matrix, Neutrino_resonant_reactions_matrix, Neutrino_decay_reactions_matrix, Neutrino_processes_matrix):

    Total_neutrino_energy = 0

    for t1 in range(len(Neutrino_unresonant_reactions_matrix[0])):
        if Neutrino_unresonant_reactions_matrix[0][t1] == 1:
            for t2 in range(len(Neutrino_unresonant_reactions_matrix[6][t1][0])):
                Total_neutrino_energy += Neutrino_unresonant_reactions_matrix[6][t1][0][t2] * Neutrino_unresonant_reactions_matrix[6][t1][1][t2]
        elif Neutrino_unresonant_reactions_matrix[0][t1] == 2:
            Total_neutrino_energy += Neutrino_unresonant_reactions_matrix[2][t1] * Neutrino_unresonant_reactions_matrix[3][t1] * Neutrino_unresonant_reactions_matrix[7][t1]
            if Neutrino_unresonant_reactions_matrix[4][t1] != '0':  # Некоторые реакции нейтринной эмиссии типа 2 имеют всего один канал.
                Total_neutrino_energy += Neutrino_unresonant_reactions_matrix[4][t1] * Neutrino_unresonant_reactions_matrix[5][t1] * Neutrino_unresonant_reactions_matrix[8][t1]

    for t1 in range(len(Neutrino_resonant_reactions_matrix[0])):
        if Neutrino_resonant_reactions_matrix[0][t1] == 1:
            for t2 in range(len(Neutrino_resonant_reactions_matrix[6][t1][0])):
                Total_neutrino_energy += Neutrino_resonant_reactions_matrix[6][t1][0][t2] * Neutrino_resonant_reactions_matrix[6][t1][1][t2]
        elif Neutrino_resonant_reactions_matrix[0][t1] == 2:
            Total_neutrino_energy += Neutrino_resonant_reactions_matrix[2][t1] * Neutrino_resonant_reactions_matrix[3][t1] * Neutrino_resonant_reactions_matrix[7][t1]
            if Neutrino_resonant_reactions_matrix[4][t1] != '0':  # Некоторые реакции нейтринной эмиссии типа 2 имеют всего один канал.
                Total_neutrino_energy += Neutrino_resonant_reactions_matrix[4][t1] * Neutrino_resonant_reactions_matrix[5][t1] * Neutrino_resonant_reactions_matrix[8][t1]

    for t1 in range(len(Neutrino_decay_reactions_matrix[0])):
        if Neutrino_decay_reactions_matrix[0][t1] == 1:
            for t2 in range(len(Neutrino_decay_reactions_matrix[6][t1][0])):
                Total_neutrino_energy += Neutrino_decay_reactions_matrix[6][t1][0][t2] * Neutrino_decay_reactions_matrix[6][t1][1][t2]
        elif Neutrino_decay_reactions_matrix[0][t1] == 2:
            Total_neutrino_energy += Neutrino_decay_reactions_matrix[2][t1] * Neutrino_decay_reactions_matrix[3][t1] * Neutrino_decay_reactions_matrix[7][t1]
            if Neutrino_decay_reactions_matrix[4][t1] != '0':  # Некоторые реакции нейтринной эмиссии типа 2 имеют всего один канал.
                Total_neutrino_energy += Neutrino_decay_reactions_matrix[4][t1] * Neutrino_decay_reactions_matrix[5][t1] * Neutrino_decay_reactions_matrix[8][t1]

    for t1 in range(len(Neutrino_processes_matrix[0])):
        if Neutrino_processes_matrix[0][t1] == 1:
            for t2 in range(len(Neutrino_processes_matrix[6][t1][0])):
                Total_neutrino_energy += Neutrino_processes_matrix[6][t1][0][t2] * Neutrino_processes_matrix[6][t1][1][t2]
        elif Neutrino_processes_matrix[0][t1] == 2:
            Total_neutrino_energy += Neutrino_processes_matrix[2][t1] * Neutrino_processes_matrix[4][t1] * Neutrino_processes_matrix[7][t1]
            if Neutrino_processes_matrix[4][t1] != '0':  # Некоторые процессы нейтринной эмиссии типа 2 имеют всего один канал.
                Total_neutrino_energy += Neutrino_processes_matrix[4][t1] * Neutrino_processes_matrix[5][t1] * Neutrino_processes_matrix[8][t1]

    return Total_neutrino_energy


def Neutrino_unresonant_reactions_matrix_create():

    with open("Simulation_input.txt", 'r') as f1:
        Ammount_of_unresonant_neutrino_channels = 0
        for t1 in f1:
            if t1[:-1] == 'Реакции:':
                break
        for t1 in f1:
            if t1[:-1] == '+Resonant reactions:':  # Блок, отвечающий за выбор правильных строк из файла Simulation_input.txt.
                break
            elif t1[:17] == '+                ' and t1[17] != ' ':
                if Is_neutrino(t1.split('_')[0].split(' ')):
                    for t3 in t1.split('_')[3][1:-2].split('/'):
                        Ammount_of_unresonant_neutrino_channels = max(Ammount_of_unresonant_neutrino_channels, int(t3.split('(')[1][0:-1]))

    Neutrino_unresonant_reactions_matrix = [[], [], [], [], [], [], [], [], [], []]  # Тип реакции (1 - континуум (нейтрино с разной энергией), 2 - нейтрино с фиксированной энергией (энергиями)), максимальная энергия нейтрино (для типа 1), (энергия канала нейтрино, вероятность канала нейтрино) 1 - 2 (для типа 2), значения спектра нейтринной эмиссии объёмом вещеста (для типа 1), значения нейтринной эмиссии фиксированных энергий по каждому из каналов распада объёмом вещества 1 - 2 (для типа 2), (нерезонансные реакции, имеющте этот нейтринный канал).
    for t1 in range(Ammount_of_unresonant_neutrino_channels):
        for i in range(0, 9):
            Neutrino_unresonant_reactions_matrix[i].append('0')
        Neutrino_unresonant_reactions_matrix[9].append([])

    t0 = 0
    with open("Simulation_input.txt", 'r') as f1_:
        for t1 in f1_:
            if t1[:-1] == 'Реакции:':
                break
        for t1 in f1_:
            if t1[:-1] == '+Resonant reactions:':  # Блок, отвечающий за выбор правильных строк из файла Simulation_input.txt.
                break
            elif t1[:17] == '+                ' and t1[17] != ' ':

                for t3 in t1.split('_')[3][1:-2].split('/'):
                    if Is_neutrino(t1.split('_')[0].split(' ')):
                        t = int(t3.split('(')[1][0:-1]) - 1
                        if t3[1:6] == 'Const':
                            Neutrino_unresonant_reactions_matrix[0][t] = 2
                        elif t3[1:5] == 'Cont':
                            Neutrino_unresonant_reactions_matrix[0][t] = 1
        
                        if Neutrino_unresonant_reactions_matrix[0][t] == 1:
                            Neutrino_unresonant_reactions_matrix[1][t] = float(t3.split('(')[0][6:-1])
        
                        elif Neutrino_unresonant_reactions_matrix[0][t] == 2:
                            for t4 in range(len(t3[7:-4].split(';'))):
                                Neutrino_unresonant_reactions_matrix[2 + 2 * t4][t], Neutrino_unresonant_reactions_matrix[3 + 2 * t4][t] = float(t3.split('(')[0][7:-1].split(';')[t4].split('-')[0]), float(t3.split('(')[0][7:-1].split(';')[t4].split('-')[1])

                        Neutrino_unresonant_reactions_matrix[9][t].append(t0)

                t0 += 1

        return Neutrino_unresonant_reactions_matrix


def Neutrino_resonant_reactions_matrix_create():

    with open("Simulation_input.txt", 'r') as f1:
        Ammount_of_resonant_neutrino_channels = 0
        for t1 in f1:
            if t1[:-1] == '+Resonant reactions:':
                break
        for t1 in f1:
            if t1[
               :-1] == '+Decay reactions:':  # Блок, отвечающий за выбор правильных строк из файла Simulation_input.txt.
                break
            elif t1[:17] == '+                ' and t1[17] != ' ':
                if Is_neutrino(t1.split('_')[0].split(' ')):
                    for t3 in t1.split('_')[2][1:-2].split('/'):
                        Ammount_of_resonant_neutrino_channels = max(Ammount_of_resonant_neutrino_channels, int(t3.split('(')[1][0:-1]))

    Neutrino_resonant_reactions_matrix = [[], [], [], [], [], [], [], [], [], []]  # Тип реакции (1 - континуум (нейтрино с разной энергией), 2 - нейтрино с фиксированной энергией (энергиями)), максимальная энергия нейтрино (для типа 1), (энергия канала нейтрино, вероятность канала нейтрино) 1 - 2 (для типа 2), значения спектра нейтринной эмиссии объёмом вещеста (для типа 1), значения нейтринной эмиссии фиксированных энергий по каждому из каналов распада объёмом вещества 1 - 2 (для типа 2), (резонансные реакции, имеющте этот нейтринный канал).
    for t1 in range(Ammount_of_resonant_neutrino_channels):
        for i in range(0, 9):
            Neutrino_resonant_reactions_matrix[i].append('0')
        Neutrino_resonant_reactions_matrix[9].append([])

    t0 = 0
    with open("Simulation_input.txt", 'r') as f1_:
        for t1 in f1_:
            if t1[:-1] == '+Resonant reactions:':
                break
        for t1 in f1_:
            if t1[
               :-1] == '+Decay reactions:':  # Блок, отвечающий за выбор правильных строк из файла Simulation_input.txt.
                break
            elif t1[:17] == '+                ' and t1[17] != ' ':

                for t3 in t1.split('_')[2][1:-2].split('/'):
                    if Is_neutrino(t1.split('_')[0].split(' ')):
                        t = int(t3.split('(')[1][0:-1]) - 1
                        if t3[1:6] == 'Const':
                            Neutrino_resonant_reactions_matrix[0][t] = 2
                        elif t3[1:5] == 'Cont':
                            Neutrino_resonant_reactions_matrix[0][t] = 1

                        if Neutrino_resonant_reactions_matrix[0][t] == 1:
                            Neutrino_resonant_reactions_matrix[1][t] = float(t3.split('(')[0][6:-1])

                        elif Neutrino_resonant_reactions_matrix[0][t] == 2:
                            for t4 in range(len(t3[7:-4].split(';'))):
                                Neutrino_resonant_reactions_matrix[2 + 2 * t4][t], Neutrino_resonant_reactions_matrix[3 + 2 * t4][t] = float(t3.split('(')[0][7:-1].split(';')[t4].split('-')[0]), float(t3.split('(')[0][7:-1].split(';')[t4].split('-')[1])

                        Neutrino_resonant_reactions_matrix[9][t].append(t0)

                t0 += 1

        return Neutrino_resonant_reactions_matrix

        
def Neutrino_decay_reactions_matrix_create():

    with open("Simulation_input.txt", 'r') as f1:
        Ammount_of_decay_neutrino_channels = 0
        for t1 in f1:
            if t1[:-1] == '+Decay reactions:':
                break
        for t1 in f1:
            if t1[
               :-1] == '+Processes:':  # Блок, отвечающий за выбор правильных строк из файла Simulation_input.txt.
                break
            elif t1[:17] == '+                ' and t1[17] != ' ':
                if Is_neutrino(t1.split('_')[0].split(' ')):
                    for t3 in t1.split('_')[2][1:-2].split('/'):
                        Ammount_of_decay_neutrino_channels = max(Ammount_of_decay_neutrino_channels, int(t3.split('(')[1][0:-1]))

    Neutrino_decay_reactions_matrix = [[], [], [], [], [], [], [], [], [], []]  # Тип реакции (1 - континуум (нейтрино с разной энергией), 2 - нейтрино с фиксированной энергией (энергиями)), максимальная энергия нейтрино (для типа 1), (энергия канала нейтрино, вероятность канала нейтрино) 1 - 2 (для типа 2), значения спектра нейтринной эмиссии объёмом вещеста (для типа 1), значения нейтринной эмиссии фиксированных энергий по каждому из каналов распада объёмом вещества 1 - 2 (для типа 2), (реакции распада, имеющте этот нейтринный канал).
    for t1 in range(Ammount_of_decay_neutrino_channels):
        for i in range(0, 9):
            Neutrino_decay_reactions_matrix[i].append('0')
        Neutrino_decay_reactions_matrix[9].append([])

    t0 = 0
    with open("Simulation_input.txt", 'r') as f1_:
        for t1 in f1_:
            if t1[:-1] == '+Decay reactions:':
                break
        for t1 in f1_:
            if t1[
               :-1] == '+Processes:':  # Блок, отвечающий за выбор правильных строк из файла Simulation_input.txt.
                break
            elif t1[:17] == '+                ' and t1[17] != ' ':

                for t3 in t1.split('_')[2][1:-2].split('/'):
                    if Is_neutrino(t1.split('_')[0].split(' ')):
                        t = int(t3.split('(')[1][0:-1]) - 1
                        if t3[1:6] == 'Const':
                            Neutrino_decay_reactions_matrix[0][t] = 2
                        elif t3[1:5] == 'Cont':
                            Neutrino_decay_reactions_matrix[0][t] = 1

                        if Neutrino_decay_reactions_matrix[0][t] == 1:
                            Neutrino_decay_reactions_matrix[1][t] = float(t3.split('(')[0][6:-1])

                        elif Neutrino_decay_reactions_matrix[0][t] == 2:
                            for t4 in range(len(t3[7:-4].split(';'))):
                                Neutrino_decay_reactions_matrix[2 + 2 * t4][t], Neutrino_decay_reactions_matrix[3 + 2 * t4][t] = float(t3.split('(')[0][7:-1].split(';')[t4].split('-')[0]), float(t3.split('(')[0][7:-1].split(';')[t4].split('-')[1])

                        Neutrino_decay_reactions_matrix[9][t].append(t0)

                t0 += 1

        return Neutrino_decay_reactions_matrix

def Neutrino_processes_matrix_create():

    with open("Simulation_input.txt", 'r') as f1:
        Ammount_of_process_neutrino_channels = 0
        for t1 in f1:
            if t1[:-1] == '+Processes:':
                break
        for t1 in f1:
            if t1[:17] == '+                ' and t1[17] != ' ':
                if Is_neutrino(t1.split('_')[0].split(' ')):
                    for t3 in t1.split('_')[2][1:-2].split('/'):
                        Ammount_of_process_neutrino_channels = max(Ammount_of_process_neutrino_channels, int(t3.split('(')[1][0:-1]))

    Neutrino_processes_matrix = [[], [], [], [], [], [], [], [], [], []]  # Тип реакции (1 - континуум (нейтрино с разной энергией), 2 - нейтрино с фиксированной энергией (энергиями)), максимальная энергия нейтрино (для типа 1), (энергия канала нейтрино, вероятность канала нейтрино) 1 - 2 (для типа 2), значения спектра нейтринной эмиссии объёмом вещеста (для типа 1), значения нейтринной эмиссии фиксированных энергий по каждому из каналов распада объёмом вещества 1 - 2 (для типа 2), (реакции распада, имеющте этот нейтринный канал).
    for t1 in range(Ammount_of_process_neutrino_channels):
        for i in range(0, 9):
            Neutrino_processes_matrix[i].append('0')
        Neutrino_processes_matrix[9].append([])

    t0 = 0
    with open("Simulation_input.txt", 'r') as f1_:
        for t1 in f1_:
            if t1[:-1] == '+Processes:':
                break
        for t1 in f1_:
            if t1[:17] == '+                ' and t1[17] != ' ':

                for t3 in t1.split('_')[2][1:-2].split('/'):
                    if Is_neutrino(t1.split('_')[0].split(' ')):
                        t = int(t3.split('(')[1][0:-1]) - 1
                        if t3[1:6] == 'Const':
                            Neutrino_processes_matrix[0][t] = 2
                        elif t3[1:5] == 'Cont':
                            Neutrino_processes_matrix[0][t] = 1

                        if Neutrino_processes_matrix[0][t] == 1:
                            Neutrino_processes_matrix[1][t] = float(t3.split('(')[0][6:-1])

                        elif Neutrino_processes_matrix[0][t] == 2:
                            for t4 in range(len(t3[7:-4].split(';'))):
                                Neutrino_processes_matrix[2 + 2 * t4][t], Neutrino_processes_matrix[3 + 2 * t4][t] = float(t3.split('(')[0][7:-1].split(';')[t4].split('-')[0]), float(t3.split('(')[0][7:-1].split(';')[t4].split('-')[1])

                        Neutrino_processes_matrix[9][t].append(t0)

                t0 += 1

        return Neutrino_processes_matrix


def Neutrino_spectrum_create(Neutrino_unresonant_reactions_matrix, Unresonant_reactions, Neutrino_resonant_reactions_matrix, Resonant_reactions, Neutrino_decay_reactions_matrix, Decay_reactions, Neutrino_processes_matrix, Processes, R_core, Distance_neutrino_spectrum):

    for t1 in range(len(Neutrino_unresonant_reactions_matrix[0])):
        if Neutrino_unresonant_reactions_matrix[0][t1] == 1:
            Energy_range = []
            E_min = 0
            E_max = Neutrino_unresonant_reactions_matrix[1][t1]
            Points_of_calculation = 1000
            E_step = (E_max - E_min) / (Points_of_calculation - 1)
            for t2 in range(Points_of_calculation):
                Energy_range.append(E_min + t2 * E_step)
            
            Neutrino_unresonant_reactions_matrix[6][t1] = [Energy_range, [Neutrino_energy(Energy, E_max) * Neutrino_flux(Distance_neutrino_spectrum, R_core, sum([Unresonant_reactions[16][Reaction_number] for Reaction_number in Neutrino_unresonant_reactions_matrix[9][t1]])) / (365.2421897 * 24 * 60 * 60) for Energy in Energy_range]]
        
        elif Neutrino_unresonant_reactions_matrix[0][t1] == 2:
            if Neutrino_unresonant_reactions_matrix[2][t1] != '0':
                Neutrino_unresonant_reactions_matrix[7][t1] = Neutrino_unresonant_reactions_matrix[3][t1] * Neutrino_flux(Distance_neutrino_spectrum, R_core, sum([Unresonant_reactions[16][Reaction_number] for Reaction_number in Neutrino_unresonant_reactions_matrix[9][t1]])) / (365.2421897 * 24 * 60 * 60)
            if Neutrino_unresonant_reactions_matrix[4][t1] != '0':
                Neutrino_unresonant_reactions_matrix[8][t1] = Neutrino_unresonant_reactions_matrix[5][t1] * Neutrino_flux(Distance_neutrino_spectrum, R_core, sum([Unresonant_reactions[16][Reaction_number] for Reaction_number in Neutrino_unresonant_reactions_matrix[9][t1]])) / (365.2421897 * 24 * 60 * 60)

    for t1 in range(len(Neutrino_resonant_reactions_matrix[0])):
        if Neutrino_resonant_reactions_matrix[0][t1] == 1:
            Energy_range = []
            E_min = 0
            E_max = Neutrino_resonant_reactions_matrix[1][t1]
            Points_of_calculation = 1000
            E_step = (E_max - E_min) / (Points_of_calculation - 1)
            for t2 in range(Points_of_calculation):
                Energy_range.append(E_min + t2 * E_step)
            
            Neutrino_resonant_reactions_matrix[6][t1] = [Energy_range, [Neutrino_energy(Energy, E_max) * Neutrino_flux(Distance_neutrino_spectrum, R_core, sum([Resonant_reactions[16][Reaction_number] for Reaction_number in Neutrino_resonant_reactions_matrix[9][t1]])) / (365.2421897 * 24 * 60 * 60) for Energy in Energy_range]]
        
        elif Neutrino_resonant_reactions_matrix[0][t1] == 2:
            if Neutrino_resonant_reactions_matrix[2][t1] != '0':
                Neutrino_resonant_reactions_matrix[7][t1] = Neutrino_resonant_reactions_matrix[3][t1] * Neutrino_flux(Distance_neutrino_spectrum, R_core, sum([Resonant_reactions[16][Reaction_number] for Reaction_number in Neutrino_resonant_reactions_matrix[9][t1]])) / (365.2421897 * 24 * 60 * 60)
            if Neutrino_spectrum_create[4][t1] != '0':
                Neutrino_resonant_reactions_matrix[8][t1] = Neutrino_resonant_reactions_matrix[5][t1] * Neutrino_flux(Distance_neutrino_spectrum, R_core, sum([Resonant_reactions[16][Reaction_number] for Reaction_number in Neutrino_resonant_reactions_matrix[9][t1]])) / (365.2421897 * 24 * 60 * 60)
    
    for t1 in range(len(Neutrino_decay_reactions_matrix[0])):
        if Neutrino_decay_reactions_matrix[0][t1] == 1:
            Energy_range = []
            E_min = 0
            E_max = Neutrino_decay_reactions_matrix[1][t1]
            Points_of_calculation = 1000
            E_step = (E_max - E_min) / (Points_of_calculation - 1)
            for t2 in range(Points_of_calculation):
                Energy_range.append(E_min + t2 * E_step)
            
            Neutrino_decay_reactions_matrix[6][t1] = [Energy_range, [Neutrino_energy(Energy, E_max) * Neutrino_flux(Distance_neutrino_spectrum, R_core, sum([Decay_reactions[19][Reaction_number] for Reaction_number in Neutrino_decay_reactions_matrix[9][t1]])) / (365.2421897 * 24 * 60 * 60) for Energy in Energy_range]]
        
        elif Neutrino_decay_reactions_matrix[0][t1] == 2:
            if Neutrino_decay_reactions_matrix[2][t1] != '0':
                Neutrino_decay_reactions_matrix[7][t1] = Neutrino_decay_reactions_matrix[3][t1] * Neutrino_flux(Distance_neutrino_spectrum, R_core, sum([Decay_reactions[19][Reaction_number] for Reaction_number in Neutrino_decay_reactions_matrix[9][t1]])) / (365.2421897 * 24 * 60 * 60)
            if Neutrino_decay_reactions_matrix[4][t1] != '0':
                Neutrino_decay_reactions_matrix[8][t1] = Neutrino_decay_reactions_matrix[5][t1] * Neutrino_flux(Distance_neutrino_spectrum, R_core, sum([Decay_reactions[19][Reaction_number] for Reaction_number in Neutrino_decay_reactions_matrix[9][t1]])) / (365.2421897 * 24 * 60 * 60)
    
    for t1 in range(len(Neutrino_processes_matrix[0])):
        if Neutrino_processes_matrix[0][t1] == 1:
            Energy_range = []
            E_min = 0
            E_max = Neutrino_processes_matrix[1][t1]
            Points_of_calculation = 1000
            E_step = (E_max - E_min) / (Points_of_calculation - 1)
            for t2 in range(Points_of_calculation):
                Energy_range.append(E_min + t2 * E_step)
            
            Neutrino_processes_matrix[6][t1] = [Energy_range, [Neutrino_energy(Energy, E_max) * Neutrino_flux(Distance_neutrino_spectrum, R_core, sum([Processes[9][Reaction_number] for Reaction_number in Neutrino_processes_matrix[9][t1]])) / (365.2421897 * 24 * 60 * 60) for Energy in Energy_range]]
        
        elif Neutrino_processes_matrix[0][t1] == 2:
            if Neutrino_processes_matrix[2][t1] != '0':
                Neutrino_processes_matrix[7][t1] = Neutrino_processes_matrix[3][t1] * Neutrino_flux(Distance_neutrino_spectrum, R_core, sum([Processes[9][Reaction_number] for Reaction_number in Neutrino_processes_matrix[9][t1]])) / (365.2421897 * 24 * 60 * 60)
            if Neutrino_processes_matrix[4][t1] != '0':
                Neutrino_processes_matrix[8][t1] = Neutrino_processes_matrix[5][t1] * Neutrino_flux(Distance_neutrino_spectrum, R_core, sum([Processes[9][Reaction_number] for Reaction_number in Neutrino_processes_matrix[9][t1]])) / (365.2421897 * 24 * 60 * 60)


    return Neutrino_unresonant_reactions_matrix, Neutrino_resonant_reactions_matrix, Neutrino_decay_reactions_matrix, Neutrino_processes_matrix

def Neutrino_spectrum_visualize(Neutrino_unresonant_reactions_matrix, Neutrino_resonant_reactions_matrix, Neutrino_decay_reactions_matrix, Neutrino_processes_matrix,  R_core, Distance_neutrino_spectrum):

    fig, ax = plt.subplots()

    for t1 in range(len(Neutrino_unresonant_reactions_matrix[0])):
        if Neutrino_unresonant_reactions_matrix[0][t1] == 1:
            ax.plot(Neutrino_unresonant_reactions_matrix[6][t1][0], Neutrino_unresonant_reactions_matrix[6][t1][1])

        elif Neutrino_unresonant_reactions_matrix[0][t1] == 2:
            if Neutrino_unresonant_reactions_matrix[2][t1] != '0':
                ax.plot([Neutrino_unresonant_reactions_matrix[2][t1], Neutrino_unresonant_reactions_matrix[2][t1]], [0, Neutrino_unresonant_reactions_matrix[7][t1]])
            if Neutrino_unresonant_reactions_matrix[4][t1] != '0':
                ax.plot([Neutrino_unresonant_reactions_matrix[4][t1], Neutrino_unresonant_reactions_matrix[4][t1]], [0, Neutrino_unresonant_reactions_matrix[8][t1]])

    for t1 in range(len(Neutrino_resonant_reactions_matrix[0])):
        if Neutrino_resonant_reactions_matrix[0][t1] == 1:
            ax.plot(Neutrino_resonant_reactions_matrix[6][t1][0], Neutrino_resonant_reactions_matrix[6][t1][1])

        elif Neutrino_resonant_reactions_matrix[0][t1] == 2:
            if Neutrino_resonant_reactions_matrix[2][t1] != '0':
                ax.plot([Neutrino_resonant_reactions_matrix[2][t1], Neutrino_resonant_reactions_matrix[2][t1]], [0, Neutrino_resonant_reactions_matrix[7][t1]])
            if Neutrino_spectrum_create[4][t1] != '0':
                ax.plot([Neutrino_resonant_reactions_matrix[4][t1], Neutrino_resonant_reactions_matrix[4][t1]], [0, Neutrino_resonant_reactions_matrix[8][t1]])

    for t1 in range(len(Neutrino_decay_reactions_matrix[0])):
        if Neutrino_decay_reactions_matrix[0][t1] == 1:
            ax.plot(Neutrino_decay_reactions_matrix[6][t1][0], Neutrino_decay_reactions_matrix[6][t1][1])

        elif Neutrino_decay_reactions_matrix[0][t1] == 2:
            if Neutrino_decay_reactions_matrix[2][t1] != '0':
                ax.plot([Neutrino_decay_reactions_matrix[2][t1], Neutrino_decay_reactions_matrix[2][t1]], [0, Neutrino_decay_reactions_matrix[7][t1]])
            if Neutrino_decay_reactions_matrix[4][t1] != '0':
                ax.plot([Neutrino_decay_reactions_matrix[4][t1], Neutrino_decay_reactions_matrix[4][t1]], [0, Neutrino_decay_reactions_matrix[8][t1]])

    for t1 in range(len(Neutrino_processes_matrix[0])):
        if Neutrino_processes_matrix[0][t1] == 1:
            ax.plot(Neutrino_processes_matrix[6][t1][0], Neutrino_processes_matrix[6][t1][1])

        elif Neutrino_processes_matrix[0][t1] == 2:
            if Neutrino_processes_matrix[2][t1] != '0':
                ax.plot([Neutrino_processes_matrix[2][t1], Neutrino_processes_matrix[2][t1]], [0, Neutrino_decay_reactions_matrix[7][t1]])
            if Neutrino_processes_matrix[4][t1] != '0':
                ax.plot([Neutrino_processes_matrix[4][t1], Neutrino_processes_matrix[4][t1]], [0, Neutrino_decay_reactions_matrix[8][t1]])

    ax.set_title('Neutrino spectrum')
    ax.set_xlabel('Energy, MeV')
    ax.set_ylabel('Neutrino flux')
    ax.set_xscale('log')
    ax.set_yscale('log')
    #     ax.set_ylim(1e+1, 10e+12)
    #     ax.set_xlim(1e-1, 25)
    ax.legend()
    plt.show()
