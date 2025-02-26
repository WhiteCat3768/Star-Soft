from Element_operations import *

from World_constants import Worldconstants as Wc  # Импорт класса, содержащего мировые константы.


def Unresonant_reactions_matrix_create(
        Elements):  # Функция для создания матрицы, хранящей информацию про нерезонансные реакции.
    with open("Simulation_input.txt", 'r') as f1, open("Simulation_input.txt", 'r') as f1_, open("Elements_data.txt",
                                                                                                 'r') as f2:

        print('Creating unresonant reaction matrix...')

        Unresonant_reactions = []
        for i in range(0, 22):
            Unresonant_reactions.append(
                [])  # Элемент 1, элемент 2, масса 1, масса 2, заряд 1, заряд 2, приведённая масса реакции, S0, символ Кронекера, тау, (элемент-результат, масса) 1 - 3, количество реакций в единицу времени, зависимость от электронов (позитронов), 1 электронный (позитронный) коэффициент пропорциональности, 2 электронный (позитронный) коэффициент пропорциональности, номер определяющей реакции (если есть), энергетический эффект реакции.
        Reaction_number = 0
        S0pp = [[], [], []]
        Elpp = [[], [], [], []]  # Созданеи матриц для факторов, определяющих скорости нерезонансных реакций.

        for t1 in f1:
            if t1[:3] == '===':
                break
        for t1 in f1:
            if t1[:4] == 'S0pp':
                S0pp[0].append(t1.split(' ')[0].split('_')[1])
                S0pp[1].append(t1.split(' ')[0].split('_')[2])
                S0pp[2].append((float(t1.split(' ')[2]) * 10 ** float(t1.split(' ')[6][1:-2])) * (
                            365.2421897 * 24 * 60 * 60))  # S0 - для секунды, а мы сичтаем для года, поэтому его надо увеличить для системы SA.
            elif t1[:4] == 'Elpp':
                Elpp[0].append(t1.split('; ')[0].split(' ')[0].split('_')[1])
                Elpp[1].append(t1.split('; ')[0].split(' ')[0].split('_')[2])
                Elpp[2].append(
                    (float(t1.split('; ')[0].split(' ')[2]) * 10 ** float(t1.split('; ')[0].split(' ')[6][1:-1])) * (
                                365.2421897 * 24 * 60 * 60))  # El_factor1 - для секунды, а мы сичтаем для года, поэтому его надо увеличить для системы SA.
                Elpp[3].append(float(t1.split('; ')[1]))
            elif t1[:4] == '===':
                break  # Блок, отвечающий за чтение факторов, определяющих скорость реакции из файла Simulation_output.txt.

        for t1 in f1_:
            if t1[:-1] == '+Unresonant reactions:':
                break
        for t1 in f1_:
            if t1[:-1] == '+Resonant reactions:':
                break
            elif t1[:17] == '+                ' and t1[17] != ' ':
                for i in range(0, 22):
                    Unresonant_reactions[i].append('0')
                Unresonant_reactions = Read_the_unresonant_reaction(t1, Reaction_number, Elements, Unresonant_reactions)
                Reaction_number += 1  # Блок, отвечающий за чтение непосредственно уравнений реакций из файла Simulation_output.txt с помощью соответствующей функции.

        for t1 in range(Reaction_number):
            Unresonant_reactions[6][t1] = str(Reduced_mass(int(Unresonant_reactions[2][t1]), int(
                Unresonant_reactions[3][t1])))  # Блок, отвечающий за расчёт приведённой массы реакции.

        for t2 in range(len(S0pp[0])):
            for t3 in range(len(Unresonant_reactions[0])):
                if not (Unresonant_reactions[17][t3]) and S0pp[0][t2] == Unresonant_reactions[0][t3] and S0pp[1][t2] == \
                        Unresonant_reactions[1][t3]:
                    Unresonant_reactions[7][t3] = str(S0pp[2][
                                                          t2])  # Блок, отвечающий за запись астрофизических факторов S0 в данные о соответствующих нерезонансных реакциях, не зависящих от электронов.

        for t2 in range(len(Elpp[0])):
            for t3 in range(len(Unresonant_reactions[0])):
                if Unresonant_reactions[17][t3] and Elpp[0][t2] == Unresonant_reactions[0][t3] and Elpp[1][t2] == \
                        Unresonant_reactions[1][t3]:
                    Unresonant_reactions[18][t3] = str(Elpp[2][t2])
                    Unresonant_reactions[19][t3] = str(Elpp[3][
                                                           t2])  # Блок, отвечающий за запись факторов El1 и El2 в данные о соответствующих нерезонансных реакциях, зависящих от электронов.

        for t2 in range(len(Unresonant_reactions[0])):
            if Unresonant_reactions[0][t2] == Unresonant_reactions[1][t2]:
                Unresonant_reactions[8][t2] = '1'
            else:
                Unresonant_reactions[8][t2] = '0'  # Блок, отвечающий за определение символа Кронекера.

        print(f"Succesfully created unresonant reaction matrix with {Reaction_number} reactions.")

    return Unresonant_reactions


def Unresonant_reactions_matrix_update(Unresonant_reactions,
                                       T):  # Функция для обновления матрицы нерезонансных реакций под текущее значение температуры ядра.
    for t2 in range(len(Unresonant_reactions[0])):
        if not (Unresonant_reactions[17][t2]):
            Unresonant_reactions[9][t2] = str(42.49 * (
                        float(Unresonant_reactions[4][t2]) ** 2 * float(Unresonant_reactions[5][t2]) ** 2 * float(
                    Unresonant_reactions[6][t2]) / (T / 10 ** 6)) ** (1 / 3))
    return Unresonant_reactions


def Unresonant_reactions_reactions_count(Unresonant_reactions, Concentrations, Mass_fractions, ro, T,
                                         Average_per_particle_weight):  # Функция для подсчёта скорости протекания нерезонансых реакций по данным из матрицы нерезонансных реакций.
    T6 = T / 10 ** 6
    for t2 in range(len(Unresonant_reactions[0])):
        if not (Unresonant_reactions[17][t2]):
            Unresonant_reactions[16][t2] = 7.2 * 10 ** (-19) * (Concentration_of_element(Unresonant_reactions[0][t2],
                                                                                         Concentrations) * Concentration_of_element(
                Unresonant_reactions[1][t2], Concentrations) * float(Unresonant_reactions[7][t2]) * float(
                Unresonant_reactions[9][t2]) ** 2 * Wc.e ** (-float(Unresonant_reactions[9][t2]))) / (
                                                       (1 + float(Unresonant_reactions[8][t2])) * float(
                                                   Unresonant_reactions[4][t2]) * float(
                                                   Unresonant_reactions[5][t2]) * float(Unresonant_reactions[6][t2]))
            if Unresonant_reactions[20][t2] != -1:
                Unresonant_reactions[16][t2] = Unresonant_reactions[16][t2] * Unresonant_reactions[16][
                    Unresonant_reactions[20][t2] - 1] / (
                                                           365.2421897 * 24 * 60 * 60)  # Блок, отвечающий за подсчёт скорости протекания нерезонансных реакций, не зависящих от электронов.
        else:
            if Unresonant_reactions[1][t2] == '0':
                Unresonant_reactions[16][t2] = float(Unresonant_reactions[18][t2]) * Concentration_of_element(
                    Unresonant_reactions[0][t2], Concentrations) * ro * (
                                                           1 + Mass_fraction_of_element('H1', Mass_fractions)) * T6 ** (
                                                   -0.5)
            else:
                Unresonant_reactions[16][t2] = float(Unresonant_reactions[18][t2]) * ro * T6 ** (-0.5) * (
                            1 + float(Unresonant_reactions[19][t2]) * (T6 - 16)) / Average_per_particle_weight
            if Unresonant_reactions[20][t2] != -1:
                Unresonant_reactions[16][t2] = Unresonant_reactions[16][t2] * Unresonant_reactions[16][
                    Unresonant_reactions[20][t2] - 1] / (
                                                           365.2421897 * 24 * 60 * 60)  # Блок, отвечающий за подсчёт скорости протекания нерезонансных реакций, зависящих от электронов.
    return Unresonant_reactions


def Unresonant_reactions_matrix_print(Unresonant_reactions):
    print('Unresonant_reactions matrix:')
    print('================')
    for t1 in range(len(Unresonant_reactions[0])):
        print(f"Reaction {t1 + 1}:")
        print(
            f"Element 1: {Unresonant_reactions[0][t1]} Mass 1: {Unresonant_reactions[2][t1]} Charge 1: {Unresonant_reactions[4][t1]}")
        print(
            f"Element 2: {Unresonant_reactions[1][t1]} Mass 2: {Unresonant_reactions[3][t1]} Charge 2: {Unresonant_reactions[5][t1]}")
        print(f"Reduced mass: {Unresonant_reactions[6][t1]}, The Kroneker symbol: {Unresonant_reactions[8][t1]}")
        print(f"Result 1: {Unresonant_reactions[10][t1]}, Mass 1: {Unresonant_reactions[13][t1]}")
        print(f"Result 2: {Unresonant_reactions[11][t1]}, Mass 2: {Unresonant_reactions[14][t1]}")
        print(f"Result 3: {Unresonant_reactions[12][t1]}, Mass 3: {Unresonant_reactions[15][t1]}")
        #         print(f"Reactions ammount per sec per cm^3: {Unresonant_reactions[16][t1]}")
        if not (Unresonant_reactions[17][t1]):
            print(f"S0 factor: {Unresonant_reactions[7][t1]}, Initial tau: {Unresonant_reactions[9][t1]}")
        else:
            print(
                f"Electron factor 1: {Unresonant_reactions[18][t1]}, Electron factor 2: {Unresonant_reactions[19][t1]}")
        if Unresonant_reactions[20][t1] != '0':
            print(f"Defining reaction number: {Unresonant_reactions[20][t1]}")
        print('================')
