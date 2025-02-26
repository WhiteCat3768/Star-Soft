import matplotlib.cm as cm
import matplotlib.pyplot as plt
import mplcursors  # Enables interactive tooltips
import numpy as np

from Element_operations import Elements_list_create

Elements = Elements_list_create()
Amount_of_elements = len(Elements[0])

T = 1.57e+7

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
