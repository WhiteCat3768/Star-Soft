import numpy as np  # Использована для частичного упрощения вычислений.

from World_constants import Worldconstants as Wc  # Импорт класса, содержащего мировые константы.

from Element_operations import *
from Differentials_solving_methods import *


def Decay_reactions_matrix_create(Elements, T):  # Функция для создания матрицы, хранящей информацию про реакции распада.
    
    with open("Simulation_input.txt", 'r') as f1, open("Simulation_input.txt", 'r') as f1_, open("Elements_data.txt", 'r') as f2:
        
        print('Creating decay reactions matrix...')
        Decay_reactions = []
        for i in range(0, 20):
            Decay_reactions.append([])  # Элемент 1, масса 1, заряд 1, спин 1, (элемент-результат, масса, заряд, спин) 1 - 2, нормальный период полураспада, константа, тип реакции (1 - альфа-распад, 2 - бетта-распад (+/-)), энергитический эффект реакции, нормальная температура, нормальная плотность вещества из элемента с изначальным ядром, период полураспада, количество реакций в единицу времени в единицу времени.
        Reaction_number = 0
        H_life = [[], []]
        Const = [[], []]  # Созданеи матриц для факторов, определяющих скорости реакций распада.
        
        for t1 in f1:
            if t1[:3] == '===':
                break
        for t1 in f1:
            if t1[:6] == 'H_life':
                H_life[0].append(t1.split(' ')[0].split('_')[2])
                H_life[1].append((float(t1.split(' ')[2]) * 10 ** float(t1.split(' ')[6][1:-2])))
            elif t1[:5] == 'Const':
                Const[0].append(t1.split(' ')[0].split('_')[1])
                Const[1].append((float(t1.split(' ')[2]) * 10 ** float(t1.split(' ')[6][1:-2])))
            elif t1[:4] == '===':
                break  # Блок, отвечающий за чтение факторов, определяющих скорость реакции из файла Simulation_output.txt.
        
        for t1 in f1_:
            if t1[:-1] == '+Decay reactions:':
                break     
        for t1 in f1_:
            if t1[:-1] == '+Processes:':
                break   
            elif t1[:17] == '+                ' and t1[17] != ' ':
                for i in range(0, 20):
                    Decay_reactions[i].append('0')
                Decay_reactions = Read_the_decay_reaction(t1, Reaction_number, Elements, Decay_reactions)
                Reaction_number += 1  # Блок, отвечающий за чтение непосредственно уравнений реакций из файла Simulation_input.txt с помощью соответствующей функции.
        
        for t2 in range(len(H_life[0])):
            for t3 in range(len(Decay_reactions[0])):
                if H_life[0][t2] == Decay_reactions[0][t3]:
                    Decay_reactions[12][t3] = float(H_life[1][t2])
                    
        for t2 in range(len(Const[0])):
            for t3 in range(len(Decay_reactions[0])):
                if Const[0][t2] == Decay_reactions[0][t3]:
                    Decay_reactions[13][t3] = float(Const[1][t2])  # Блок, отвечающий за запись константы зависимости периода полураспада от внешних условий для реакций распада.
        
        for t1 in range(len(Decay_reactions[0])):
            Decay_reactions[16][t1] = Wc.T0  # Блок, отвечающий за запись нормальной температуры.
        
        for t1 in range(len(Decay_reactions[0])):
            Decay_reactions[17][t1] = Wc.p0 * Mass_of_element(Decay_reactions[0][t1], Elements) / (Wc.R * Wc.T0)  # Блок, отвечающий за запись нормальной плотности (при атмосферном давлении) для вещества, состоящего из атомов с ядрами, для которых происходит соответствующая реакция распада.
        
        print(f"Succesfully created decay reactions matrix with {Reaction_number} reactions.")
        
    return Decay_reactions


def Decay_reactions_matrix_update(Decay_reactions, ro, T):  # Функция для обновления матрицы реакций распада под текущие значения температуры и плотности ядра.
    
    for t1 in range(len(Decay_reactions[0])):
        if Decay_reactions[14][t1] == 1:
            Lambda_normal = np.log(2) / Decay_reactions[12][t1]
            Lambda_current = Lambda_normal * (1 + Decay_reactions[13][t1] * (T / Decay_reactions[16][t1]))
            Decay_reactions[18][t1] = np.log(2) / Lambda_current  # Блок, отвечающий за подсчёт изменённого периода полураспада для реакций альфа-распада от температуры ядра.
        else:
            Lambda_normal = np.log(2) / Decay_reactions[12][t1]
            Lambda_current = Lambda_normal * (1 + Decay_reactions[13][t1] * (T / Decay_reactions[16][t1]) * (ro / Decay_reactions[17][t1]) ** -(3 / 2))
            Decay_reactions[18][t1] = np.log(2) / Lambda_current  # Блок, отвечающий за подсчёт изменённого периода полураспада для реакций бетта-распада от температуры и плотности ядра.
            
    return Decay_reactions


def Decay_reactions_count(Decay_reactions, Concentrations):  # Функция для подсчёта скорости протекания реакций распада по данным из матрицы реакций распада.
    
    for t1 in range(len(Decay_reactions[0])):
        Decay_reactions[19][t1] = 365.2421897 * 24 * 60 * 60 * (Concentration_of_element(Decay_reactions[0][t1], Concentrations) * np.log(2) / Decay_reactions[18][t1])
    
    return Decay_reactions


        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        