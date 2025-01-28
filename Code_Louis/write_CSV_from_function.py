import numpy as np
import csv
from scipy.interpolate import interp1d

# Définition de la fonction
mag_profile = lambda x, profile_factor: np.exp(-profile_factor * (x - 1)**2)

# Paramètres
profile_factor = 4.0  # Valeur à ajuster selon tes besoins
x_values = np.linspace(0, 2, 100)  # Générer 100 valeurs de x entre 0 et 2

# Calcul des valeurs de y en utilisant la fonction
y_values = mag_profile(x_values, profile_factor)

# Écriture dans le fichier CSV
with open('output.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['x', 'y'])  # En-têtes
    for x, y in zip(x_values, y_values):
        writer.writerow([x, y])

print("Les données ont été écrites dans 'output.csv'.")

import matplotlib.pyplot as plt

# Interpolation
f_interpolated = interp1d(x_values, y_values, kind='cubic')
x_interpolated = np.linspace(0, 2, 400)
y_interpolated = f_interpolated(x_interpolated)

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(x_values, y_values, 'o', label='Original Data')
plt.plot(x_interpolated, y_interpolated, '-', label='Interpolated Data')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Original Function and Interpolation')
plt.legend()
plt.grid(True)
plt.show()