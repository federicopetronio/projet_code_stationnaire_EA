import numpy as np
import subprocess
import matplotlib.pyplot as plt
import sys
from hall_papier import plot_densities, plot_electron_temperature, plot_ion_velocity, plot_electric_field_and_ionization_source, plot_magnetic_field, plot_fluxes

# Exécuter les deux scripts pour générer les fichiers de données
subprocess.run([sys.executable, "hall1.py"])
subprocess.run([sys.executable, "hall_papier.py"])

# Charger les données de chaque fichier
data1 = np.loadtxt("values.csv", delimiter=',')
data_anna = np.loadtxt("values_anna.csv", delimiter=',')
data_anna2 = np.loadtxt("values_anna_Te.csv", delimiter=',')

# Extraire les données X et Y

x, n_i, n_g, T_e, u_i, E_x, S_iz, BB, vng, vni = data1.T
z, z_bar, n_i_anna, n_g_anna, Te_end, u_i_anna, Ez, S_iz_anna, BB_anna, vng_anna, vni_anna = data_anna.T
z_bar_array, T_bar_array = data_anna2.T

####### Comparaison des données #######
#Tracer les deux courbes ensemble
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()

# Tracé des premières données
plot_densities(z_bar, n_i_anna, n_g_anna, 'n_i_anna', 'n_g_anna', 'blue', 'green', ax1=ax1, ax2=ax2)

# Tracé des deuxièmes données (sur les mêmes axes)
plot_densities(x, n_i, n_g, 'n_i', 'n_g', 'red', 'orange', ax1=ax1, ax2=ax2)

# Affichage final
plt.xlim([0, 1])  # Limites communes pour l'axe x
plt.tight_layout()  # Ajustement automatique des espaces
plt.show()

################
plt.figure()
#plot_electron_temperature(z_bar_array, T_bar_array, 'blue')
plot_electron_temperature(z_bar, Te_end, 'blue')
plot_electron_temperature(x, T_e, 'red')
plt.legend(["Anna", "v1"])

# ################# Siz et Ez
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()

# Tracé des premières données
plot_electric_field_and_ionization_source(z_bar, Ez, S_iz_anna, 'Ez', 'S_iz_anna', 'blue', 'green', ax1=ax1, ax2=ax2)

# Tracé des deuxièmes données (sur les mêmes axes)
plot_electric_field_and_ionization_source(x, E_x, S_iz, 'E_x', 'S_iz_anna', 'red', 'orange', ax1=ax1, ax2=ax2)

# Affichage final
plt.xlim([0, 1])  # Limites communes pour l'axe x
plt.tight_layout()  # Ajustement automatique des espaces
plt.show()


plt.figure()
plot_ion_velocity(z_bar, u_i_anna, color='blue')
plot_ion_velocity(x, u_i, color='red')
plt.legend(["Anna", "v1"])

# plt.figure()
# plot_magnetic_field(z_bar, BB_anna)
# plot_magnetic_field(x, BB)
# plt.legend(["Anna", "v1"])

plt.figure()
plot_fluxes(z_bar, vng_anna, vni_anna, 'blue','green')
plot_fluxes(x, vng, vni, 'red', 'orange')
plt.legend(["vng_Anna", "vng_anna", "vng_v1", "vni_v1"])

plt.show()