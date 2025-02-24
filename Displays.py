def Simulation_state_display(Mass_fractions):  # Функция для отображения массовых долей элементов.
    Answer = ''
    for t1 in range(len(Mass_fractions[0])):
        Answer = Answer + f"{Mass_fractions[0][t1]}: {Mass_fractions[1][t1]} \t"
    return Answer


def Concentrations_state_display(Concentrations):  # Функция для отображения концентраций элементов.
    Answer = ''
    for t1 in range(len(Concentrations[0])):
        Answer = Answer + f"{Concentrations[0][t1]}: {Concentrations[1][t1]} \t"
    return Answer


def Elements_burning_speed_state_display(Elements_burning_speed):  # Функция для отображения скоростей выгорания элементов.
    Answer = ''
    for t1 in range(len(Elements_burning_speed[0])):
        Answer = Answer + f"{Elements_burning_speed[0][t1]}: {Elements_burning_speed[1][t1]} \t"
    return Answer


def Elements_burning_time_state_display(Elements_burning_time):  # Функция для отображения времени выгорания элементов.
    Answer = ''
    for t1 in range(len(Elements_burning_time[0])):
        Answer = Answer + f"{Elements_burning_time[0][t1]}: {Elements_burning_time[1][t1]} \t"
    return Answer



def Simulation_state_display_simpler(Mass_fractions):  # Функция для отображения массовых долей элементов в упрощённом виде.
    Answer = ''
    for t1 in range(len(Mass_fractions[0])):
#         if str(Mass_fractions[1][t1]).find('e') != -1:
#             Answer = Answer + str(Mass_fractions[1][t1]).split('e')[0][:-6] + 'e' + str(Mass_fractions[1][t1]).split('e')[1] + '\t'
#         else:
#             Answer = Answer + str(Mass_fractions[1][t1]) + '\t'
            Answer = Answer + str(Mass_fractions[1][t1]) + ','
    return Answer


def El_pos_state_display(Electrons_concentration, Positrons_concentration):  # Функция для отображения среднего морлекулярного веса на один электрон и позитрон.
    return f" {Electrons_concentration}, {Positrons_concentration},"


def Time_display(Time, Time_limit):  # Функция для отображения прошедшего времени симуляции и процентного отношения прошедшего времени симуляции к заданному.
    Percentage = (Time / Time_limit) * 100
    Years = Time
#     Time -= Years * 365.2425 * 24 * 60 * 60
#     Days = Time // (24 * 60 * 60)
#     Time -= Days * 24 * 60 * 60
#     Hours = Time // (60 * 60)
#     Time -= Hours * 60 * 60
#     Minutes = Time // 60
#     Time -= Minutes * 60
#     Seconds = Time

#     return (f"{Years}y, {Percentage}%")
    return (f"{Years}y, {Percentage}%,")

#    return (f"{Years}y {Days}d {Hours}h {Minutes} m {Seconds}s, {Percentage}%")


def Temperature_display(Temperature):  # Функция для отображения температуры.
    return(f"{Temperature} K")


def General_condition_dispay(Concentrations, Elements_burning_time, Electrons_burning_time, Positrons_burning_time):  # Функция для отображения состояния стабильности симуляции. Наиболее частые признаки нестабильности - уменьшение какой-либо величины из концентраций или времён выгорания (элементов, электронов или позитронов) до слишком низких значений.
    Concentrations_cleared = []
    for t1 in Concentrations[1]:
        if t1 != 0:
            Concentrations_cleared.append(abs(t1))
    Burning_time_cleared = []
    for t1 in Elements_burning_time[1]:
        if t1 != 0:
            Burning_time_cleared.append(abs(t1))
    if Electrons_burning_time != 0:
        Burning_time_cleared.append(abs(Electrons_burning_time))
    if Positrons_burning_time != 0:
        Burning_time_cleared.append(abs(Positrons_burning_time))
    if min(Concentrations_cleared) < 10e-10 or min(Burning_time_cleared) < 1e-15:
        return 'Critical!'
    elif min(Concentrations_cleared) < 10e0 or min(Burning_time_cleared) < 1e-10:
        return 'Suspicious...'
    else:
        return 'Stable'


def Debug_display(Unresonant_reactions, Concentrations, Electrons_concentration, Positrons_concentration):  # Неиспользуемая функция для отображения температуры дополнительных параметров.
    print(Concentrations, Electrons_concentration, Positrons_concentration)