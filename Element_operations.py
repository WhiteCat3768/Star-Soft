from World_constants import Worldconstants as Wc  # Импорт класса, содержащего мировые константы.


def Elements_list_create():  # Создание массива содержащего данные про элементы в формате: [Элемент, масса ядра, заряд ядра, спин ядра, чётность ядра].
    with open("Elements_data.txt", 'r') as f3:
        print('Creating list of elements...')
        Ammount_of_elements = 0
        Elements = [[], [], [], [], []]
        for t1 in f3:
            Ammount_of_elements += 1
            Element = t1.split(' ')
            Elements[0].append(Element[0])
            Elements[1].append(float(Element[1]))
            Elements[2].append(int(Element[2]))
            Elements[3].append(float(Element[3][:-2]))
            if Element[3][-2:-1] == '+':
                Elements[4].append(1)
            if Element[3][-2:-1] == '-':
                Elements[4].append(-1)
        print(f"Successfully created list of elements with {Ammount_of_elements} elements.")
        print(Elements)
    return Elements


def Is_element_new(x, Elements):  # Низкоуровневая техническая функция, использующаяся при создании матрицы элементов.
    for t1 in Elements[0]:
        if x == t1:
            return False
    return True


def Mass_of_element(x, Elements):  # Техническая функция для поиска в массиве Elements массы ядра заданного элемента.
    for t1 in range(len(Elements[0])):
        if Elements[0][t1] == x:
            return Elements[1][t1]


def Charge_of_element(x, Elements):  # Техническая функция для поиска в массиве Elements заряда ядра заданного элемента.
    for t1 in range(len(Elements[0])):
        if Elements[0][t1] == x:
            return Elements[2][t1]


def Spin_of_element(x, Elements):  # Техническая функция для поиска в массиве Elements спина ядра заданного элемента.
    for t1 in range(len(Elements[0])):
        if Elements[0][t1] == x:
            return Elements[3][t1]


def Parity_of_element(x,
                      Elements):  # Техническая функция для поиска в массиве Elements чётности ядра заданного элемента.
    for t1 in range(len(Elements[0])):
        if Elements[0][t1] == x:
            return Elements[4][t1]


def Concentration_of_element(x,
                             Concentrations):  # Техническая функция для поиска в массиве Concentrations концентрации ядер заданного элемента.
    if x == '0':
        return 1
    else:
        for t in range(len(Concentrations[0])):
            if Concentrations[0][t] == x:
                return float(Concentrations[1][t])


def Mass_fraction_of_element(x,
                             Mass_fractions):  # Техническая функция для поиска в массиве Mass_fractions массовой доли ядер заданного элемента от вещества ядра.
    for t in range(len(Mass_fractions[0])):
        if Mass_fractions[0][t] == x:
            return float(Mass_fractions[1][t])


def Max_concentration_of_element(x, Elements,
                                 ro):  # Техническая функция для рассчёта максимально возможной концентрации элеметнов. Используется для опредления времени выгорания элемента, если его концентрация увеличивается.
    Max_concentration = float(ro * 1 / (Mass_of_element(x, Elements) * Wc.M_nuc))
    return Max_concentration


def Compound_core(x, y,
                  Elements):  # Техническая функция для резонансных реакций, определяющая условное ядро, возникшее до моментального альфа-распада в некоторых резонансных реакциях с выходом He4.
    Mass = round(Mass_of_element(x, Elements)) + int(Mass_of_element(y, Elements))
    Charge = Charge_of_element(x, Elements) + Charge_of_element(y, Elements)
    for t1 in range(len(Elements[0])):
        if Mass == round(Elements[1][t1]) and Charge == Elements[2][t1]:
            return Elements[0][t1]


def Reduced_mass(m1, m2):  # Техническая функция для резонансных реакций, рассчитывающая приведённую массу реакции.
    return m1 * m2 / (m1 + m2)


def Average_per_particle_weight_count(El_con, Pos_con, Elements,
                                      Concentrations):  # Функция для рассчёта среднего молекулярного веса (средняя масса любой частицы).
    Total_atom_mass_per_volume = 0
    for t in range(len(Elements[0])):
        Total_atom_mass_per_volume += Mass_of_element(Elements[0][t], Elements) * Concentration_of_element(
            Elements[0][t], Concentrations)
    Total_concentration_per_volume = 0
    for t in range(len(Elements[0])):
        Total_concentration_per_volume += Concentration_of_element(Elements[0][t], Concentrations)
    Total_concentration_per_volume += El_con
    #    Total_concentration_per_volume += Pos_con
    Average_per_particle_weight = Total_atom_mass_per_volume / Total_concentration_per_volume
    return Average_per_particle_weight


def Read_element(x,
                 Elements):  # Функция для парсинга, объединяющая результаты функций для нахождения массы, заряда и спина ядра элемента по заданному элементу.
    with open("Elements_data.txt", 'r') as f1:
        for t1 in range(len(Elements[0])):
            if x == Elements[0][t1]:
                Mass = Elements[1][t1]
                Charge = Elements[2][t1]
                Spin = Elements[3][t1]
                Parity = Elements[4][t1]
                break
    return x, Mass, Charge, Spin, Parity


def Read_the_unresonant_reaction(line, Reaction_number, Elements,
                                 Unresonant_reactions):  # Функция для парсинга, предназначенная для записи части данных про реакцию в матрицу нерезонансных реакций.
    with open("Elements_data.txt", 'r') as f1:
        Reaction, Defining_reaction_number, Temperature_effect, *_ = line.split('_')
        Unresonant_reactions[20][Reaction_number] = int(Defining_reaction_number[
                                                        1:-1])  # Номер реакции в матрице, от которой зависит скосроть этой реакции, если такой реакции нет, то это значение равно -1.
        Unresonant_reactions[21][Reaction_number] = float(
            Temperature_effect[1:-1]) * 1e+22 / 6.24e+13  # Энергетический эффект реакции в SAДжоулях.

        Components = Reaction.split()[1:]
        Component = 0
        Result = 0
        t = 0
        Results = False
        Unresonant_reactions[17][Reaction_number] = False
        while t <= len(Components) - 1:  # Цикл, отдельно считывающий каждый компонент реакции.
            if Components[t] == '=':
                Results = True
                t += 1
            if Components[t] != '+':
                if Components[t] == 'e+' or Components[
                    t] == 'e-':  # Определение того, зависит ли реакция от электронов / позитронов (есть ли они среди реагентов).
                    if not (Results):
                        Unresonant_reactions[17][Reaction_number] = True
                if Components[t] != 'e+' and Components[t] != 'e-':
                    if Components[t] != 'y' and Components[t] != 'v':
                        if Results is False:
                            Unresonant_reactions[Component][Reaction_number], Unresonant_reactions[Component + 2][
                                Reaction_number], Unresonant_reactions[Component + 4][
                                Reaction_number], *_ = Read_element(Components[t],
                                                                    Elements)  # Считываем нужную информацию для элемента-реагента реакции.
                            Component += 1
                        else:
                            Element_info = Read_element(Components[t],
                                                        Elements)  # Считываем нужную информацию для элемента-продукта реакции.
                            Unresonant_reactions[Result + 10][Reaction_number] = Element_info[0]
                            Unresonant_reactions[Result + 13][Reaction_number] = Element_info[1]
                            Result += 1
            t += 1
        return Unresonant_reactions


def Read_the_resonant_reaction(line, Reaction_number, Elements,
                               Resonant_reactions):  # Функция для парсинга, предназначенная для записи части данных про реакцию в матрицу резонансных реакций.
    with open("Elements_data.txt", 'r') as f1:
        Reaction, Temperature_effect, *_ = line.split('_')
        Resonant_reactions[21][Reaction_number] = float(
            Temperature_effect[1:-1]) * 1e+22 / 6.24e+13  # Энергетический эффект реакции в SAДжоулях.
        Resonant_reactions[26][Reaction_number] = Reaction[17:]

        Components = Reaction.split()[1:]
        Component = 0
        Result = 0
        t = 0
        Results = False
        while t <= len(Components) - 1:
            if Components[t] == '=':
                Results = True
                t += 1
            if Components[t] != '+':
                if Components[t] != 'e+' and Components[t] != 'e-':
                    if Components[t] != 'y' and Components[t] != 'v':
                        if Results is False:
                            Resonant_reactions[Component][Reaction_number], Resonant_reactions[Component + 2][
                                Reaction_number], Resonant_reactions[Component + 4][Reaction_number], \
                            Resonant_reactions[Component + 6][Reaction_number], *_ = Read_element(Components[t],
                                                                                                  Elements)  # Считываем нужную информацию для элемента-реагента реакции.
                            Component += 1
                        else:
                            Element_info = Read_element(Components[t],
                                                        Elements)  # Считываем нужную информацию для элемента-продукта реакции.
                            Resonant_reactions[Result + 9][Reaction_number] = Element_info[0]
                            Resonant_reactions[Result + 11][Reaction_number] = Element_info[1]
                            Resonant_reactions[Result + 13][Reaction_number] = Element_info[2]
                            Resonant_reactions[Result + 15][Reaction_number] = Element_info[3]
                            Result += 1
            t += 1

        return Resonant_reactions


def Read_the_decay_reaction(line, Reaction_number, Elements,
                            Decay_reactions):  # Функция для парсинга, предназначенная для записи части данных про реакцию в матрицу реакций распада.
    with open("Elements_data.txt", 'r') as f1:
        Reaction, Temperature_effect, *_ = line.split('_')
        Decay_reactions[15][Reaction_number] = float(
            Temperature_effect[1:-1]) * 1e+22 / 6.24e+13  # Энергетический эффект реакции в SAДжоулях.
        print(Reaction)

        Components = Reaction.split()[1:]
        Component = 0
        Result = 0
        Alpha_decay = False  # Если значение остаётся False, значит реакция бетта-распада.
        t = 0
        Results = False
        while t <= len(Components) - 1:
            if Components[t] == '=':
                Results = True
                t += 1
            if Components[t] != '+':
                if Components[t] == 'e+' or Components[t] == 'e-':
                    if Results:
                        Electrons_in_results = True
                if Components[t] != 'e+' and Components[t] != 'e-':
                    if Components[t] != 'y' and Components[t] != 'v':
                        if Results is False:
                            Decay_reactions[Component][Reaction_number], Decay_reactions[Component + 1][
                                Reaction_number], Decay_reactions[Component + 2][Reaction_number], \
                            Decay_reactions[Component + 3][Reaction_number], *_ = Read_element(Components[t],
                                                                                               Elements)  # Считываем нужную информацию для элемента-реагента реакции.
                            Component += 1
                        else:
                            Element_info = Read_element(Components[t],
                                                        Elements)  # Считываем нужную информацию для элемента-продукта реакции.
                            Decay_reactions[Result + 4][Reaction_number] = Element_info[0]
                            Decay_reactions[Result + 6][Reaction_number] = Element_info[1]
                            Decay_reactions[Result + 8][Reaction_number] = Element_info[2]
                            Decay_reactions[Result + 10][Reaction_number] = Element_info[3]
                            if Element_info[0] == 'He4':
                                Alpha_decay = True
                            Result += 1
            t += 1

        if Alpha_decay:  # Запись информации про тип распада (от этого зависит корректировка времени полураспада).
            Decay_reactions[14][Reaction_number] = 1
        else:
            Decay_reactions[14][Reaction_number] = 2

        return Decay_reactions


def Read_the_process(line, Reaction_number, Elements,
                     Processes):  # Функция для парсинга, предназначенная для записи части данных про процесс в матрицу процессов.
    with open("Elements_data.txt", 'r') as f1:
        Reaction, Temperature_effect, *_ = line.split('_')
        print(Reaction)
        Processes[7][Reaction_number] = float(
            Temperature_effect[1:-1]) * 1e+22 / 6.24e+13  # Энергетический эффект реакции в SAДжоулях.

        Components = Reaction.split()[1:]
        Component = 0
        Result = 0
        t = 0
        Results = False
        while t <= len(Components) - 1:
            if Components[t] == '=':
                Results = True
                t += 1
            if Components[t] != '+':
                if Components[t] != 'e+' and Components[t] != 'e-':
                    if Components[t] != 'y' and Components[t] != 'v':
                        if Results is False:
                            Processes[Component][Reaction_number], *_ = Read_element(Components[t],
                                                                                     Elements)  # Считываем нужную информацию для элемента-реагента реакции.
                            Component += 1
                        else:
                            Element_info = Read_element(Components[t],
                                                        Elements)  # Считываем нужную информацию для элемента-продукта реакции.
                            Processes[Result + 3][Reaction_number] = Element_info[0]
                            Result += 1
            t += 1

        return Processes
