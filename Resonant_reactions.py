import matplotlib.pyplot as plt  # Библиотека для построения графиков.
import numpy as np  # Использована для частичного упрощения вычислений.

from World_constants import Worldconstants as Wc  # Импорт класса, содержащего мировые константы.

from Element_operations import *
from Differentials_solving_methods import *


def Reaction_index_define(Resonant_reactions, Reaction_key):
    for t2 in range(len(Resonant_reactions[0])):
        if Reaction_key == Resonant_reactions[26][t2]:
            return t2, True
    return -1, False

def Integration_trapz(x, y):  # Техническая функция для численного интегрирования площади под графиками методом трапеций (используется для нормирвоания).
    Integral = 0
    for t1 in range(1, len(x)):
        Integral += (x[t1] - x[t1 - 1]) * (y[t1] + y[t1 - 1]) / 2
    return Integral


def Resonant_reactions_matrix_create(Elements, T):  # Функция для создания матрицы, хранящей информацию про резонансные реакции.
    
    with open("Simulation_input.txt", 'r') as f1, open("Elements_data.txt", 'r') as f2, open("Resonant_reactions_resonance_parameters.txt", 'r') as f3:
        
        print('Creating resonant reaction matrix...')

        Resonant_reactions = []
        for i in range(0, 27):
            Resonant_reactions.append([])  # Элемент 1, элемент 2, масса 1, масса 2, заряд 1, заряд 2, спин 1, спин 2, приведённая масса реакции, (элемент-результат, масса, заряд, спин) 1 - 2, [парцианая ширина входных каналов] 1 - n, [парцианая ширина выходных каналов] 1 - n, [полная ширина резонансов] 1 - n, [энергии резонансов] 1 - n, энергитический эффект реакции, температура для которой посчитано сечение образования составного ядра, [зависимости сечений от энергии ядер] 1 - n, [<сечение * скорость> образования составного ядра для каждого из резонансов] 1 - n, количество реакций в единицу времени в единицу времени, ключ реакции (запись в виде строки в Resonant_reactions_resonance_parameters.txt).
        Reaction_number = 0

        for t1 in f1:
            if t1[:-1] == '+Resonant reactions:':
                break
        for t1 in f1:
            if t1[:-1] == '+Decay reactions:':
                break
            elif t1[:17] == '+                ' and t1[17] != ' ':
                for i in range(0, 17):
                    Resonant_reactions[i].append('0')
                for i in range(17, 21):
                    Resonant_reactions[i].append([])
                for i in range(21, 23):
                    Resonant_reactions[i].append('0')
                for i in range(23, 25):
                    Resonant_reactions[i].append([])
                for i in range(25, 27):
                    Resonant_reactions[i].append('0')
                Resonant_reactions = Read_the_resonant_reaction(t1, Reaction_number, Elements, Resonant_reactions)
                Reaction_number += 1  # Блок, отвечающий за чтение непосредственно уравнений реакций из файла Simulation_input.txt с помощью соответствующей функции.
        
        for t1 in range(Reaction_number):
            Resonant_reactions[8][t1] = Reduced_mass(int(Resonant_reactions[2][t1]), int(Resonant_reactions[3][t1]))  # Блок, отвечающий за расчёт приведённой массы реакции.

        Reading_resonances = False
        if len(Resonant_reactions[0]) > 0:
            for t1 in f3:
                if t1[0] != '[':
                    Reaction_key = t1[:-2]
                    Reaction_index, Reading_resonances = Reaction_index_define(Resonant_reactions, Reaction_key)
                else:
                    if Reading_resonances and t1 != '':
                        E_res, Wd_in, Wd_out, Wd_res = [float(x) for x in t1[1:-2].split("; ")]
                        Resonant_reactions[17][Reaction_index].append(E_res * 1e+22 / 1.602e+13)
                        Resonant_reactions[18][Reaction_index].append(Wd_in * 1e+22 / 1.602e+13)
                        Resonant_reactions[19][Reaction_index].append(Wd_out * 1e+22 / 1.602e+13)
                        Resonant_reactions[20][Reaction_index].append(E_res * 1e+22 / 1.602e+13)
                        Resonant_reactions[23][Reaction_index].append('0')
                        Resonant_reactions[24][Reaction_index].append('0')

        for t1 in range(len(Resonant_reactions[0])):

            for t3 in range(len(Resonant_reactions[17][t1])):

                n = 200  # Параметр, определяющий количество полных ширин резонанса, которые составляют половину ширины диапозона интегрирования энергии при подсчёте сечений резонансных реакций.

                E_min = Resonant_reactions[20][t1][t3] - Resonant_reactions[19][t1][t3] * n
                E_max = Resonant_reactions[20][t1][t3] + Resonant_reactions[19][t1][t3] * n
                Points_of_calculation_energy = 10
                Collision_energy_range = []
                E_step = (E_max - E_min) / (Points_of_calculation_energy - 1)
                for t2 in range(Points_of_calculation_energy):
                    Collision_energy_range.append(E_min + t2 * E_step)  # Блок, отвечающий за диапазона рассматриваемых суммарных энергий системы из двух ядер для дальнейшего вычисления значения сечений реакции для значений энергии из диапазона.

                Collision_energy_range_specified = Collision_energy_range[:]  # Блок, отвечающий за выбор из диапазона энергий окресности энергии резонанса (все остальные значения энергии обнуляются) Это сделано для оптимизации, а также из-за границ применимости формулы Брейта-Вигнера.
                for t2 in range(len(Collision_energy_range_specified)):
                    if Collision_energy_range_specified[t2] < (Resonant_reactions[20][t1][t3] - n * Resonant_reactions[19][t1][t3]) or Collision_energy_range_specified[t2] > (Resonant_reactions[20][t1][t3] + n * Resonant_reactions[19][t1][t3]):
                        Collision_energy_range_specified[t2] = 0  # Выделяем область резонанса (50 полных ширин резонанса в каждую сторону от энергии резонанса.)
            
                Resonant_reactions[23][t1][t3] = [Collision_energy_range_specified, []]
                # Сечение по формуле Брейта-Вигнера.
                if Resonant_reactions[10][t1] == '0':
                    Spin_factor = (2 * Spin_of_element(Resonant_reactions[9][t1], Elements) + 1) / ((2 * Spin_of_element(Resonant_reactions[0][t1], Elements) + 1) * (2 * Spin_of_element(Resonant_reactions[1][t1], Elements) + 1))
                else:
                    Spin_factor = (2 * Spin_of_element(Compound_core(Resonant_reactions[9][t1], Resonant_reactions[10][t1], Elements), Elements) + 1) / ((2 * Spin_of_element(Resonant_reactions[0][t1], Elements) + 1) * (2 * Spin_of_element(Resonant_reactions[1][t1], Elements) + 1))
                for Energy in Resonant_reactions[23][t1][t3][0]:
                    if Energy == 0:
                        Resonant_reactions[23][t1][t3][1].append(0)
                    else:
                        Constant = np.pi * (Wc.h_ ** 2 / (2 * Resonant_reactions[8][t1] * Wc.M_nuc * Energy))
                        Resonance = Resonant_reactions[17][t1][t3] * Resonant_reactions[18][t1][t3] / ((Energy - Resonant_reactions[20][t1][t3]) ** 2 + (Resonant_reactions[19][t1][t3] / 2) ** 2)
                        Cross_section = Constant * Spin_factor * Resonance
                        Resonant_reactions[23][t1][t3][1].append(Cross_section)  # Блок, отвечающий за расчёт и запись в матрицу значения сечений резонансной реакции для соответствующих значений энергии после выбора из диапозона окресностей энергии резонанса.
        
        print(f"Succesfully created resonant reaction matrix with {Reaction_number} reactions.")

    return Resonant_reactions


def Resonant_reactions_matrix_update(Resonant_reactions, Concentrations, T):  # Функция для обновления матрицы резонансных реакций под текущее значение температуры ядра.

    Alpha = 0.01

    Amount_of_resonances = 0
    for t1 in range(len(Resonant_reactions[0])):
        Amount_of_resonances += len(Resonant_reactions[17][t1])

    Percentage = 0
    for t1 in range(len(Resonant_reactions[0])):

        for t5 in range(len(Resonant_reactions[17][t1])):

            Average_speed_cross_section = 0
            Relative_speed_distribution = []
            V_min = 0
            V_max = Wc.c * Alpha
            Points_of_calculation_speed = 2000
            Relative_speed_range = []
            V_step = (V_max - V_min) / (Points_of_calculation_speed - 1)  # Блок, отвечающий за создание диапазона рассматриваемых относительных скоростей ядер.

            for t2 in range(Points_of_calculation_speed):
                Relative_speed_range.append(V_min + t2 * V_step)
            for Relative_speed in Relative_speed_range:
                Relative_speed_distribution.append(4 * np.pi * (Resonant_reactions[8][t1] * Wc.M_nuc / (2 * np.pi * Wc.k * T)) ** (3 / 2) * Relative_speed ** 2 * np.exp(-Resonant_reactions[8][t1] * Wc.M_nuc * Relative_speed ** 2 / (2 * Wc.k * T)))
            Norm_factor = Integration_trapz(Relative_speed_range, Relative_speed_distribution)
            Relative_speed_distribution_normalized = []
            for t3 in Relative_speed_distribution:
                Relative_speed_distribution_normalized.append(t3 / Norm_factor)  # Блок, отвечающий за расчёт и нормализацию распределения плотности вероятности относительной скорости ядер по рассматриваемому диапазону.

            t3 = 0
            for Relative_speed in Relative_speed_range:

                n = 500  # Параметр, определяющий количество полных ширин резонанса, которые составляют половину ширины диапозона интегрирования энергии при подсчёте сечений резонансных реакций.

                Energy_distribution = []
                E_min = Resonant_reactions[20][t1][t5] - Resonant_reactions[19][t1][t5] * n
                E_max = Resonant_reactions[20][t1][t5] + Resonant_reactions[19][t1][t5] * n
                Points_of_calculation_energy = 10
                Energy_range = []
                Energy_1_range = []
                E_step = (E_max - E_min) / (Points_of_calculation_energy - 1)
                for t4 in range(Points_of_calculation_energy):
                    Energy_range.append(E_min + t4 * E_step)  # Блок, отвечающий за диапазона рассматриваемых суммарных энергий системы из двух ядер для конкретного значения из диапазона рассматривавемых относительных скоростей.

                for Energy in Energy_range:
                    Energy_relative = Resonant_reactions[8][t1] * Wc.M_nuc * Relative_speed ** 2 / 2
                    Energy_additional = Energy - Energy_relative
                    if Energy_additional <= 0:
                        Energy_distribution.append(0)
                    else:
                        Energy_distribution.append(Energy_additional ** 0.5 * np.exp(-Energy_additional / (Wc.k * T)))
                Norm_factor = Integration_trapz(Energy_range, Energy_distribution)
                Energy_distribution_normalized = []
                if Norm_factor != 0:
                    for t6 in Energy_distribution:
                        Energy_distribution_normalized.append(t6 / Norm_factor)  # Блок, отвечающий за расчёт и нормализацию распределения плотности вероятности суммарной энергии системы из двух ядер по рассматриваемому диапазону для конкретного значения относительной скорости.
                else:
                    Energy_distribution_normalized = Energy_distribution[:]

                Cross_section = 0
                for Energy in Energy_distribution:
                    Cross_section += Resonant_reactions[23][t1][t5][1][t4] * Energy  # Блок, отвечающий за подсчёт сечения согласно значению энергии и соответствующей ей плотности вероятности по распределению.

                Average_speed_cross_section += Relative_speed * Cross_section * Relative_speed_distribution[Relative_speed_range.index(Relative_speed)]
                t3 += 1
                Percentage += 1 / (Amount_of_resonances * Points_of_calculation_speed)
                print(f"{Percentage * 100}% updated...")  # Блок, отвечающий за вывод в консоль процента обновления матрицы резонансных реакций (так как это очень дорогая по времени операция).

            Resonant_reactions[24][t1][t5] = Average_speed_cross_section
        Resonant_reactions[22][t1] = T  # Блок, отвечающий за запись обновлённого под новое значение температуры ядра сечение реакции, а также значения температуры, для которого оно было посчитано в матрицу резонансных реакций.

    return Resonant_reactions
                
def Resonant_reactions_reactions_count(Resonant_reactions, T, Elements, Concentrations):  # Функция для подсчёта скорости протекания резонансых реакций по данным из матрицы резонансных реакций.
    
    for t1 in range(len(Resonant_reactions[0])):
        Resonant_reactions[25][t1] = (365.2421897 * 24 * 60 * 60) * sum(Resonant_reactions[24][t1]) * Concentration_of_element(Resonant_reactions[0][t1], Concentrations) * Concentration_of_element(Resonant_reactions[1][t1], Concentrations)
    return Resonant_reactions


def Resonant_reactions_matrix_print(Resonant_reactions):  # Функция для вывода матрицы резонансных реакций в консоль (пока без нормального интерфейса).
    print('Resonant_reactions matrix:')
    print(Resonant_reactions)
