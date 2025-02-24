from Element_operations import *

from World_constants import Worldconstants as Wc  # Импорт класса, содержащего мировые константы.

def Electron_positron_unresonant_matrix_create(R):  # Функция для парсинга, которая читает уравнения нерезонансных реакций и создаёт матрицу, хранящую эффект на количество электронов / позитронов от одной такой реакции.
    with open("Simulation_input.txt", 'r') as f1:
        Electron_unresonant_reaction_effect = [0] * R
        Positron_unresonant_reaction_effect = [0] * R
        t = 0
        for t1 in f1:
            if t1[:-1] == 'Реакции:':
                break
        for t1 in f1:
            if t1[:-1] == '+Resonant reactions:':  # Блок, отвечающий за выбор правильных строк из файла Simulation_input.txt.
                break
            elif t1[:17] == '+                ' and t1[17] != ' ':
                Reaction = t1.split('_')[0]
                Electrons_add, Positrons_add = Read_the_reaction_el_pos(Reaction)
                Electron_unresonant_reaction_effect[t] += Electrons_add
                Positron_unresonant_reaction_effect[t] += Positrons_add  # Блок, отвечающий за заполнение матрицы.
                t += 1
        return Electron_unresonant_reaction_effect, Positron_unresonant_reaction_effect


def Electron_positron_resonant_matrix_create(R):  # Функция для парсинга, которая читает уравнения резонансных реакций и создаёт матрицу, хранящую эффект на количество электронов / позитронов от одной такой реакции.
    with open("Simulation_input.txt", 'r') as f1:
        Electron_resonant_reaction_effect = [0] * R
        Positron_resonant_reaction_effect = [0] * R
        t = 0
        for t1 in f1:
            if t1[:-1] == '+Resonant reactions:':
                break
        for t1 in f1:
            if t1[:-1] == '+Decay reactions:':  # Блок, отвечающий за выбор правильных строк из файла Simulation_input.txt.
                break
            elif t1[:17] == '+                ' and t1[17] != ' ':
                Reaction = t1.split('_')[0]
                Electrons_add, Positrons_add = Read_the_reaction_el_pos(Reaction)
                Electron_resonant_reaction_effect[t] += Electrons_add
                Positron_resonant_reaction_effect[t] += Positrons_add  # Блок, отвечающий за заполнение матрицы.
                t += 1
        return Electron_resonant_reaction_effect, Positron_resonant_reaction_effect


def Electron_positron_decay_matrix_create(R):  # Функция для парсинга, которая читает уравнения реакций распада и создаёт матрицу, хранящую эффект на количество электронов / позитронов от одной такой реакции.
    with open("Simulation_input.txt", 'r') as f1:
        Electron_decay_reaction_effect = [0] * R
        Positron_decay_reaction_effect = [0] * R
        t = 0
        for t1 in f1:
            if t1[:-1] == '+Decay reactions:':
                break
        for t1 in f1:
            if t1[:-1] == '+Processes:':  # Блок, отвечающий за выбор правильных строк из файла Simulation_input.txt.
                break
            elif t1[:17] == '+                ' and t1[17] != ' ':
                Reaction = t1.split('_')[0]
                Electrons_add, Positrons_add = Read_the_reaction_el_pos(Reaction)
                Electron_decay_reaction_effect[t] += Electrons_add
                Positron_decay_reaction_effect[t] += Positrons_add  # Блок, отвечающий за заполнение матрицы.
                t += 1
        return Electron_decay_reaction_effect, Positron_decay_reaction_effect


def Electron_positron_processes_matrix_create(R):  # Функция для парсинга, которая читает уравнения процессов и создаёт матрицу, хранящую эффект на количество электронов / позитронов от одного такого процесса.
    with open("Simulation_input.txt", 'r') as f1:
        Electron_process_effect = [0] * R
        Positron_process_effect = [0] * R
        t = 0
        for t1 in f1:
            if t1[:-1] == '+Processes:':  # Блок, отвечающий за выбор правильных строк из файла Simulation_input.txt.
                break
        for t1 in f1:
            if t1[:17] == '+                ' and t1[17] != ' ':
                Reaction = t1.split('_')[0]
                Electrons_add, Positrons_add = Read_the_reaction_el_pos(Reaction)
                Electron_process_effect[t] += Electrons_add
                Positron_process_effect[t] += Positrons_add  # Блок, отвечающий за заполнение матрицы.
                t += 1
        return Electron_process_effect, Positron_process_effect


def Electron_positron_initial_concentration_count(Elements, Concentrations):  # Функция, отвечающая за поиск изначальной концентрации электронов, исходя из предположения общей электрической нейтральности звезды.
    Electrons_concentration = 0
    Positrons_concentration = 0
    for t1 in range(len(Elements[0])):
        Charge = Concentrations[1][t1] * Elements[2][t1]
        if Charge > 0:
            Electrons_concentration += Charge
        elif Charge < 0:
            Positrons_concentration += -Charge  # Блок, отвечающий за поиск концентраций на основе подсчёта суммарного заряда элементов.
    return Electrons_concentration, Positrons_concentration


def Electron_positron_burning_speed(Electron_unresonant_reaction_effect, Positron_unresonant_reaction_effect, Unresonant_reactions, Electron_resonant_reaction_effect, Positron_resonant_reaction_effect, Resonant_reactions, Electron_decay_reaction_effect, Positron_decay_reaction_effect, Decay_reactions, Electron_process_effect, Positron_process_effect, Processes):
    Electrons_burning_speed = 0  # Функция, отвечающая за подсчёт скорости изменения концентрации электронов и позитронов по матрицам электронно-позитронного эффекта для всех типов реакций.
    Positrons_burning_speed = 0
    for t in range(len(Unresonant_reactions[0])):
        Electrons_burning_speed += Electron_unresonant_reaction_effect[t] * float(Unresonant_reactions[16][t])
        Positrons_burning_speed += Positron_unresonant_reaction_effect[t] * float(Unresonant_reactions[16][t])
    for t in range(len(Resonant_reactions[0])):
        Electrons_burning_speed += Electron_resonant_reaction_effect[t] * float(Resonant_reactions[25][t])
        Positrons_burning_speed += Positron_resonant_reaction_effect[t] * float(Resonant_reactions[25][t])
    for t in range(len(Decay_reactions[0])):
        Electrons_burning_speed += Electron_decay_reaction_effect[t] * float(Decay_reactions[19][t])
        Positrons_burning_speed += Positron_decay_reaction_effect[t] * float(Decay_reactions[19][t])
    for t in range(len(Processes[0])):
        Electrons_burning_speed += Electron_process_effect[t] * float(Processes[9][t])
        Positrons_burning_speed += Positron_process_effect[t] * float(Processes[9][t])
    return Electrons_burning_speed, Positrons_burning_speed


def Read_the_reaction_el_pos(Reaction):  # Функция для парсинга, которая непосредственно просматривает компоненты реакции (реагенты и продукты) на наличие электронов или позитронов, считая электронно-позитронные эффекты конкретной реакции.
    Components = Reaction.split()[1:]
    Electrons_add = 0
    Positrons_add = 0
    t = 0
    Results = False
    while t <= len(Components) - 1:
        if Components[t] == '=':
            Results = True
            t += 1
        if Components[t] == 'e-':
            if Results:
                Electrons_add += 1
            else:
                Electrons_add -= 1
        if Components[t] == 'e+':
            if Results:
                Positrons_add += 1
            else:
                Positrons_add -= 1  # Для электрона / позитрона в продуктах - +1, в реагентах - -1.
        t += 1
    return Electrons_add, Positrons_add


def El_pos_burning_time_count(Electrons_burning_speed, Positrons_burning_speed, Electrons_concentration, Positrons_concentration):
    if Electrons_concentration == 0 or Electrons_burning_speed <= 0:
        Electrons_burning_time = float('+inf')
    else:
        Electrons_burning_time = abs(Electrons_concentration / Electrons_burning_speed)
    if Positrons_concentration == 0 or Positrons_burning_speed <= 0:
        Positrons_burning_time = float('+inf')
    else:
        Positrons_burning_time = abs(Positrons_concentration / Positrons_burning_speed)
    return(Electrons_burning_time, Positrons_burning_time)


def El_pos_annihilation(Electrons_concentration, Positrons_concentration):  # Анигилятся электронов и позитронов считается мгновенной.
    if Electrons_concentration >= Positrons_concentration:
        Annihilation_energy = 2 * Positrons_concentration * 1e+22 / 6.24e+13
        Positrons_concentration = 0
        Electrons_concentration -= Positrons_concentration
    else:
        Annihilation_energy = 2 * Electrons_concentration * 1e+22 / 6.24e+13
        Electrons_concentration = 0
        Positrons_concentration -= Electrons_concentration
    return Electrons_concentration, Positrons_concentration, Annihilation_energy
