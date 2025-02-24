# -*- coding: utf-8 -*-

from datetime import datetime  # Метод для вывода времени в Logs.txt.
import matplotlib.pyplot as plt  # Библиотека для построения графиков.
import matplotlib.cm as cm
import mplcursors  # Enables interactive tooltips
import numpy as np  # Использована для частичного упрощения вычислений.
# import scipy

from World_constants import Worldconstants as Wc  # Импорт класса, содержащего мировые константы.
from Advanced_tweakables import Advancedtweakbles as Ac  # Импорт класса, содержащего специфические параметры настройки симуляции.

from Element_operations import *
from Electron_positron_operations import *
from Temperature_change import *
from Neutrino_spectrum import *
from Differentials_solving_methods import *
from Displays import *
from Unresonant_reactions import *
from Resonant_reactions import *
from Decay_reactions import *
from Processes import *  # Импорт методов из других модулей программы.


Iterations = 1
Physical_characteristics = [0] * 12
with open("Simulation_input.txt", 'r') as f1:
    for t1 in range(1, 25):
        line = f1.readline()
        if t1 % 2 == 0:
            Physical_characteristics[int(t1 / 2) - 1] = line
R_core = int(Physical_characteristics[0]) * 1e+5
ro = int(Physical_characteristics[1])
T = int(Physical_characteristics[2])
T_surface = int(Physical_characteristics[3])
R_star = int(Physical_characteristics[4]) * 1e+5
Mass_fractions_input = Physical_characteristics[5][:-1].split('; ')
Time_limit = float(Physical_characteristics[6]) * 1e+6
Time_step_limit = float(Physical_characteristics[7])
if Physical_characteristics[8].split()[0] == 'Curant':
    Curant = True
    RK4 = False
    if Physical_characteristics[8].split()[1] == 'Force':
        Speed = float(Physical_characteristics[8].split()[2])
        Smart = False
    else:
        Smart = True
elif Physical_characteristics[8].split()[0] == 'Runge-Kutta_4':
    RK4 = True
    Curant = False
    Speed = float(Physical_characteristics[8].split()[1])
if Physical_characteristics[9][:-1] == 'True':
    Temperature_change_conf = True
else:
    Temperature_change_conf = False
Frequency = int(Physical_characteristics[10])
Test = 1
if Physical_characteristics[11][:-1] == 'File':
    File = True
else:
    File = False  # Это блок отвечает за считывание основных параметров из Simulation_input.txt.

with open("Simulation_input.txt", 'r') as f1, open("Simulation_output.txt", 'w') as f2, open("Elements_data.txt", 'r') as f3, open("Logs.txt", 'w') as f4:  # Начало блока, отвечающего за подготовку к запуску основного цикла вычиселний.
#     if T < 14000000:
#         print('Слишком низкая темпераутра в ядре, термоядерные реакции невозможны. Минимум - 14000000 K  .', file = f2)
#     else:

    Max_up_concentration_change = 1e+0

    Elements = Elements_list_create()  # Создание массива содержащего данные про элементы в формате: [Элемент, масса ядра, заряд ядра, спин ядра, чётность ядра].
    Amount_of_elements = len(Elements[0])
    Burning_speed = [[], []]
    Elements_burning_time = [[], []]

    # Система измерений: грамм, год, сантиметр, моль, Кельвин, элементарный заряд, Гаусс.
    # Остальные единицы измерений: SA[name] SA - Star Astro:
    # SAДжоуль = 1e-22 Джоуль;
    # SAПаскаль = 1e-16 Паскаль;
    # SAВатт = 1e-29 Ватт;
    #

    Distance_neutrino_spectrum = 1.5e+8  # В километрах.
    Neutrino_spectrum_shown = False

    Mo = 4 / 3 * Wc.Pi * R_core ** 3 * ro  # Вычисление массы ядра звезды.

    Mass_fractions = [[], []]  # Распаковка входного списка массовых долей элементов из Simulation_input.txt в формате [[Элементы](программа использует эту часть массива для навигации), [Массовые доли в соответствующем элементам порядке]].
    for t2 in Mass_fractions_input:
        Mass_fractions[0].append(t2.split(': ')[0])
        Mass_fractions[1].append(float(t2.split(': ')[1]))

    Concentrations = [[], []]
    Concentrations[0][:] = Mass_fractions[0][:]
    for t2 in range(len(Mass_fractions[0])):  # По массиву Mass_fractions создаётся массив с концентрациями с тем же форматом хранения данных.
        Concentrations[1].append(float(ro * float(Mass_fractions[1][t2]) / (Mass_of_element(Mass_fractions[0][t2], Elements) * Wc.M_nuc)))  # Комплекс стабильности элементов.

    for t2 in Concentrations[0]:
        Elements_burning_time[0].append(t2)
        Elements_burning_time[1].append(0)
        Burning_speed[0].append(t2)
        Burning_speed[1].append(0)  # В том же формате по массиву Concentrations для будущих рассчётов создаются массивы со скоростью выгорания элементов, а также

    print(f"[{datetime.now().time()}]>>>Creating unresonant reaction matrix...", file = f4)
    Unresonant_reactions = Unresonant_reactions_matrix_create(Elements)  # Создание матрицы, хранящей информацию про все нерезонансные реакции.
    print(f"[{datetime.now().time()}]>>>Succesfully created unresonant reaction matrix with {len(Unresonant_reactions[0])} reactions.", file = f4)
    Unresonant_reactions = Unresonant_reactions_matrix_update(Unresonant_reactions, T)  # Обновление матрицы, хранящей информацию про все нерезонансные реакции, под определённое значение температуры.
    Unresonant_reactions_matrix_print(Unresonant_reactions)
#     print(Unresonant_reactions)

    print(f"[{datetime.now().time()}]>>>Creating resonant reactions matrix...", file = f4)
    Resonant_reactions = Resonant_reactions_matrix_create(Elements, T)  # Создание матрицы, хранящей информацию про все резонансные реакции.
    print(f"[{datetime.now().time()}]>>>Successfully created resonant reaction matrix with {len(Resonant_reactions[0])} reactions.", file = f4)
    Resonant_reactions = Resonant_reactions_matrix_update(Resonant_reactions, Concentrations, T)  # Обновление матрицы, хранящей информацию про все резонансные реакции, под определённое значение температуры.
#     Resonant_reactions_matrix_print(Resonant_reactions)
    print(Resonant_reactions)

    print(f"[{datetime.now().time()}]>>>Creating decay reactions matrix...", file = f4)
    Decay_reactions = Decay_reactions_matrix_create(Elements, T)  # Создание матрицы, хранящей информацию про все реакции ядерного распада.
    print(f"[{datetime.now().time()}]>>>Successfully created decay reaction matrix with {len(Decay_reactions[0])} reactions.", file = f4)
    Decay_reactions = Decay_reactions_matrix_update(Decay_reactions, ro, T)  # Обновление матрицы, хранящей информацию про все реакции ядерного распада, под определённые значение температуры и плотности.
#     print(Decay_reactions)

    print(f"[{datetime.now().time()}]>>>Creating processes matrix...", file = f4)
    Processes = Processes_matrix_create(Elements)  # Создание матрицы, хранящей информацию про все реакции процессы.
    print(f"[{datetime.now().time()}]>>>Succesfully created processes matrix with {len(Processes[0])} processes.", file = f4)
    Processes = Processes_matrix_update(Processes, T)  # Обновление матрицы, хранящей информацию про все процессы, под определённое значение температуры

    Electrons_concentration = 0
    Positrons_concentration = 0
    Electrons_concentration, Positrons_concentration = Electron_positron_initial_concentration_count(Elements, Concentrations)  # На основе того, что суммарный заряд звезды равен нулю, рассчитывание изначальной концентрации электронов.
    Electron_unresonant_reaction_effect, Positron_unresonant_reaction_effect = Electron_positron_unresonant_matrix_create(len(Unresonant_reactions[0]))
    Electron_resonant_reaction_effect, Positron_resonant_reaction_effect = Electron_positron_resonant_matrix_create(len(Resonant_reactions[0]))
    Electron_decay_reaction_effect, Positron_decay_reaction_effect = Electron_positron_decay_matrix_create(len(Decay_reactions[0]))
    Electron_process_effect, Positron_process_effect = Electron_positron_processes_matrix_create(len(Processes[0]))  # По уравнениям соответствующего типа для электронов и позтронов создаётся матрица, хранящая эффект на их количество от одной такой реакции.
    Average_per_particle_weight = Average_per_particle_weight_count(Electrons_concentration, Positrons_concentration, Elements, Concentrations)

    Neutrino_unresonant_reactions_matrix = Neutrino_unresonant_reactions_matrix_create()
    Neutrino_resonant_reactions_matrix = Neutrino_resonant_reactions_matrix_create()
    Neutrino_decay_reactions_matrix = Neutrino_decay_reactions_matrix_create()
    Neutrino_processes_matrix = Neutrino_processes_matrix_create()

    Total_mass_fraction = 0
    for t1 in Mass_fractions[1]:
        Total_mass_fraction += t1
    Mass_fraction_fix = 1 / Total_mass_fraction  # Важно! Этот кусок кода - не костыль, а вынужденное нормирование концентраций. На каждой итерации из-за неточностей работы с числами типа float звезда теряет миллиардную долю от массы, но на масштабах моделирования это приводит к потере ею её заметной части, этот код это предотвращает.
    for t1 in range(len(Mass_fractions[0])):
        Mass_fractions[1][t1] *= Mass_fraction_fix

    if File is True:
#         print(Simulation_state_display(Mass_fractions) + Time_display(1, Time_limit), file = f2)
        print(Time_display(0, Time_limit) + Simulation_state_display_simpler(Mass_fractions) + ' ' + El_pos_state_display(Electrons_concentration, Positrons_concentration) + ' ' + Temperature_display(T), file = f2)
    else:
#         print(Simulation_state_display(Mass_fractions) + Time_display(1, Time_limit))
        print(Time_display(0, Time_limit) + Simulation_state_display_simpler(Mass_fractions) + ' ' + El_pos_state_display(Average_per_electron_atom_mass, Average_per_positron_atom_mass) + ' ' + Temperature_display(T))

    print(f"[{datetime.now().time()}]>>>Starting simulation...", file = f4)
#         print(f"[{datetime.now().time()}]>>>{Time_display(0, Time_limit)} {Simulation_state_display_simpler(Mass_fractions)} {El_pos_state_display(Electrons_concentration, Positrons_concentration)} {Temperature_display(T)}", file = f4)
    Time = 0
    while Time <= Time_limit:

        Concentrations[0][:] = Mass_fractions[0][:]
        Concentrations[1] = []
        for t2 in range(len(Mass_fractions[0])):  # По массиву Mass_fractions обновляем массив с концентрациями с тем же форматом хранения данных.
            Concentrations[1].append(float(ro * float(Mass_fractions[1][t2]) / (Mass_of_element(Mass_fractions[0][t2], Elements) * Wc.M_nuc)))

        Unresonant_reactions = Unresonant_reactions_matrix_update(Unresonant_reactions, T)
        if len(Resonant_reactions[0]) > 0 and (Resonant_reactions[22][0] / T < 0.99 or Resonant_reactions[22][0] / T > 1.01):  # Пока для оптимизацйии вставлен такой костыль, так как обновление матрицы резонансных реакций - очень дорогая операция.
            Resonant_reactions = Resonant_reactions_matrix_update(Resonant_reactions, Concentrations, T)
        Decay_reactions = Decay_reactions_matrix_update(Decay_reactions, ro, T)
        Processes = Processes_matrix_update(Processes, T)

        Test += 1
        if Iterations == 40:
            Test *= 2
            print(Test)

        if Curant:  # Ядро всей симуляции - метод дифференцирования, определящий шаг времени, а также изменения всех величин.
            Burning_speed, Speed, Electrons_burning_speed, Positrons_burning_speed, dt, Elements_burning_time, Temperature_change_speed, Temperature_change_time = Curant_diff_solve(Unresonant_reactions, Resonant_reactions, Decay_reactions, Processes, Speed, Time, Smart, Concentrations, Mass_fractions, R_core, T, ro, Average_per_particle_weight, Electron_unresonant_reaction_effect, Positron_unresonant_reaction_effect, Electron_resonant_reaction_effect, Positron_resonant_reaction_effect, Electron_decay_reaction_effect, Positron_decay_reaction_effect, Electron_process_effect, Positron_process_effect, Neutrino_unresonant_reactions_matrix, Neutrino_resonant_reactions_matrix, Neutrino_decay_reactions_matrix, Neutrino_processes_matrix, Elements_burning_time, Burning_speed, Elements, Electrons_concentration, Positrons_concentration, T_surface, R_star, Distance_neutrino_spectrum)

        if RK4:
            Burning_speed, Speed, Electrons_burning_speed, Positrons_burning_speed, dt, Elements_burning_time, Temperature_change_speed, Temperature_change_time = Runge_Kutta_4(Unresonant_reactions, Resonant_reactions, Decay_reactions, Processes, Speed, Time, Concentrations, Mass_fractions, R_core, T, ro, Average_per_particle_weight, Electron_unresonant_reaction_effect, Positron_unresonant_reaction_effect, Electron_resonant_reaction_effect, Positron_resonant_reaction_effect, Electron_decay_reaction_effect, Positron_decay_reaction_effect, Electron_process_effect, Positron_process_effect, Neutrino_unresonant_reactions_matrix, Neutrino_resonant_reactions_matrix, Neutrino_decay_reactions_matrix, Neutrino_processes_matrix, Elements_burning_time, Burning_speed, Elements, Electrons_concentration, Positrons_concentration, T_surface, R_star, Temperature_change_conf, Distance_neutrino_spectrum)

        Time += dt

        # Вычисляем скорости убывания элементов:
        for t2 in range(len(Concentrations[0])):
            Concentrations[1][t2] += Burning_speed[1][t2] * dt

        for t2 in range(len(Mass_fractions[0])):  # пересчёт концентраций в массовые доли.
            Mass_fractions[1][t2] = Concentrations[1][t2] * Wc.M_nuc * Mass_of_element(Mass_fractions[0][t2], Elements) / ro

        Total_mass_fraction = 0
        for t1 in Mass_fractions[1]:
            Total_mass_fraction += t1
        Mass_fraction_fix = 1 / Total_mass_fraction  # Важно! Этот кусок кода - не костыль, а вынужденное нормирование концентраций. На каждой итерации из-за неточностей работы с числами типа float звезда теряет миллионную долю от массы, но на масштабах моделирования это приводит к потере ею около четверти массы, этот код это предотвращает.
        for t1 in range(len(Mass_fractions[0])):
            Mass_fractions[1][t1] *= Mass_fraction_fix

        Electrons_concentration += Electrons_burning_speed * dt
        Positrons_concentration += Positrons_burning_speed * dt  # Расчёт изменения концетрации электронов и позитронов.
        Electrons_concentration, Positrons_concentration, Annihilation_energy = El_pos_annihilation(Electrons_concentration, Positrons_concentration)
        Average_per_particle_weight = Average_per_particle_weight_count(Electrons_concentration, Positrons_concentration, Elements, Concentrations)  # Рассчёт среднего молекулярного веса (вся масса в эдинице объёма, делённая на концентрацию всех частиц в объёме).

        if Temperature_change_conf:
            T += Temperature_change_speed * dt
            T += El_pos_annihilation_temperature_change(Annihilation_energy, Concentrations, R_core)
        else:
            T += 0  # Блок, отвечающий за расчёт изменения температуры (если оно включено).

        if Time > Ac.Neutrino_spectrum_time and not Neutrino_spectrum_shown:
            Neutrino_unresonant_reactions_matrix, Neutrino_resonant_reactions_matrix, Neutrino_decay_reactions_matrix, Neutrino_processes_matrix = Neutrino_spectrum_create(Neutrino_unresonant_reactions_matrix, Unresonant_reactions, Neutrino_resonant_reactions_matrix, Resonant_reactions, Neutrino_decay_reactions_matrix, Decay_reactions, Neutrino_processes_matrix, Processes, R_core, Distance_neutrino_spectrum)
            Neutrino_spectrum_visualize(Neutrino_unresonant_reactions_matrix, Neutrino_resonant_reactions_matrix, Neutrino_decay_reactions_matrix, Neutrino_processes_matrix,  R_core, Distance_neutrino_spectrum)
            Neutrino_spectrum_shown = True

        Iterations += 1

        if Iterations % Frequency == 0:  # Вывод информации в консоль и файлы Simulation_output.txt и Logs.txt.
            if File is True:
#                 print(Simulation_state_display(Mass_fractions) + Time_display(Time, Time_limit), file = f2)
                    print(Time_display(Time, Time_limit) + ' ' + Simulation_state_display_simpler(Mass_fractions) + ' ' + El_pos_state_display(Electrons_concentration, Positrons_concentration) + ' ' + Temperature_display(T), file = f2)
                    print(Time_display(Time, Time_limit))
                    print(min(Elements_burning_time[1]), Elements_burning_time[0][Elements_burning_time[1].index(min(Elements_burning_time[1]))])
                    print(Concentrations)
                    print(Burning_speed)
                    print(Elements_burning_time)
                    print(f"{T} K, {Temperature_change_time}")
#                 Debug_display(Unresonant_reactions, Concentrations, Electrons_concentration, Positrons_concentration)
#                 print(Test, file = f2)
            else:
#                 print(Simulation_state_display(Mass_fractions) + Time_display(Time, Time_limit))
                print(Time_display(Time, Time_limit) + ' ' + Simulation_state_display_simpler(Mass_fractions) + ' ' + El_pos_state_display(Electrons_concentration, Positrons_concentration) + ' ' + Temperature_display(T))
#                 Debug_display(Unresonant_reactions, Concentrations, Electrons_concentration, Positrons_concentration)
                print(dt)
            print(f"[{datetime.now().time()}]>>>Iteration: {Iterations}, Time passed: {Time_display(Time, Time_limit)}, Time step: {dt}y",file=f4)
            print(f"[{datetime.now().time()}]>>>General info: Mass fractions: {Simulation_state_display(Mass_fractions)} Electron/positron consentrations: {El_pos_state_display(Electrons_concentration, Positrons_concentration)} Temperature: {Temperature_display(T)}", file = f4)
            print(f"[{datetime.now().time()}]>>>Slowest element: {min(Elements_burning_time[1])}, {Elements_burning_time[0][Elements_burning_time[1].index(min(Elements_burning_time[1]))]}", file = f4)
            print(f"[{datetime.now().time()}]>>>Reactions speed: Unresonant: {Unresonant_reactions[16]} Resonant: {Resonant_reactions[25]} Decay: {Decay_reactions[19]} Processes: {Processes[9]}", file = f4)
            print(f"[{datetime.now().time()}]>>>Concentrations: {Concentrations_state_display(Concentrations)}", file = f4)
            print(f"[{datetime.now().time()}]>>>Elements_burning_speed: {Elements_burning_speed_state_display(Burning_speed)}", file = f4)
            print(f"[{datetime.now().time()}]>>>Elements_burning_time: {Elements_burning_time_state_display(Elements_burning_time)}", file = f4)
            if Electrons_burning_speed == 0 and Positrons_burning_speed == 0:
                print(f"[{datetime.now().time()}]>>>Simulation condition: {General_condition_dispay(Concentrations, Elements_burning_time, 1e+300, 1e+300)}", file = f4)
            elif Electrons_burning_speed == 0:
                print(f"[{datetime.now().time()}]>>>Simulation condition: {General_condition_dispay(Concentrations, Elements_burning_time, 1e+300, Positrons_concentration / Positrons_burning_speed)}", file = f4)
            elif Positrons_burning_speed == 0:
                print(f"[{datetime.now().time()}]>>>Simulation condition: {General_condition_dispay(Concentrations, Elements_burning_time, Electrons_concentration / Electrons_burning_speed, 1e+300)}", file = f4)
            else:
                print(f"[{datetime.now().time()}]>>>Simulation condition: {General_condition_dispay(Concentrations, Elements_burning_time, Electrons_concentration / Electrons_burning_speed, Positrons_concentration / Positrons_burning_speed)}", file = f4)

#                 print(Test)
#                 Test += 1
    print('Simulation finished!')


if File is True:
    print('Visualizing data...')

lines = []
with open('Simulation_output.txt', 'r') as f2:
    # Читаем все строки из файла
    for line in f2:
        lines.append(line)

# Инициализируем список для хранения данных
Time_values = []
Element_values = []  # List of lists for element data
Electron_values = []
Positron_values = []
Temperature_values = []

# Итерируемся по нечетным строкам и извлекаем значения
for t1 in range(0, len(lines)):
    values = lines[t1].strip().split(',')

    # Extracting the time value
    Time_values.append(float(values[0][:-1]))  # Извлекаем время без единиц измерения

    # Extracting element values dynamically
    Element_values_row = [float(x) for x in values[2:Amount_of_elements + 2]]
    Element_values.append(Element_values_row)  # Store all element values in a row

    # Extracting electron and positron values
    El_pos_values = [float(x) for x in values[Amount_of_elements + 2:Amount_of_elements + 4]]
    Electron_values.append(El_pos_values[0])
    Positron_values.append(El_pos_values[1])  # If needed

    # Extracting temperature values
    Temperature_values.append(float(values[Amount_of_elements + 4][:-2]) / T)

# Create the plot
plt.figure(figsize=(10, 6))
plt.xscale('log')
plt.yscale('log')

Element_values_transposed = list(zip(*Element_values))
Colors = cm.viridis(np.linspace(0, 1, Amount_of_elements))  # Use a colormap

# Store line objects for interactive legend
lines = []

# Iterate through the elements and plot each one dynamically
for i in range(len(Elements[0])):  # Iterate over the elements
    line, = plt.plot(Time_values, Element_values_transposed[i], color=Colors[i], label=Elements[0][i])
    lines.append(line)

# Plot temperature values
temp_line, = plt.plot(Time_values, Temperature_values, 'g*', label='T')
lines.append(temp_line)

# Add labels, title
plt.xlabel('Time')
plt.ylabel('Element/el_pos/temperature values')
plt.title('Star core evolution')

# Enable interactive hover tooltips
cursor = mplcursors.cursor(lines, hover=True)
cursor.connect("add", lambda sel: sel.annotation.set_text(sel.artist.get_label()))

# Show the plot
plt.show()
