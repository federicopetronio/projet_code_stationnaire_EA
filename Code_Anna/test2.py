import numpy as np
import subprocess
import matplotlib.pyplot as plt
import sys

################ Plot functions ################

def plot_densities(z_bar, n_i, n_g, label_i, label_g, color_i, color_g, ax1=None, ax2=None):
    if ax1 is None or ax2 is None:  # Si les axes ne sont pas fournis, créez-en
        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()
    
    # Tracé des densités ioniques sur l'axe principal
    line1, = ax1.plot(z_bar, n_i, color=color_i, linewidth=1, label=label_i)
    ax1.set_xlabel(r'$\overline{z}$', fontsize=14)
    ax1.set_ylabel('$n_i$ (m$^{-3}$)', fontsize=14, color=color_i)
    ax1.tick_params(axis='y', labelcolor=color_i)
    
    # Tracé des densités neutres sur l'axe secondaire
    line2, = ax2.plot(z_bar, n_g, color=color_g, linewidth=1, label=label_g)
    ax2.set_ylabel('$n_g$ (m$^{-3}$)', fontsize=14, color=color_g)
    ax2.tick_params(axis='y', labelcolor=color_g)
    
    # Titre commun
    ax1.set_title('Ion and Neutral Densities', fontsize=16)
    
    ax1.legend(handles=[line1], loc='upper left', fontsize=12)  # Légende pour ax1
    ax2.legend(handles=[line2], loc='upper right', fontsize=12)  # Légende pour ax2
    
    return ax1, ax2  # Renvoie les axes pour des ajustements supplémentaires


#def plot_electron_temperature(z_bar_array, T_bar_array, color):
def plot_electron_temperature(z_bar, Temp, color):
    plt.plot(z_bar, Temp, linewidth=1, color = color)
    plt.xlim([0, 1])
    plt.xlabel(r'$\overline{z}$', fontsize=14)
    plt.ylabel('$T_e$ (V)', fontsize=14)
    plt.gca().tick_params(labelsize=14)
    plt.title('Electron temperature')

def plot_ion_velocity(z_bar, u_i, color):
    plt.plot(z_bar, u_i, linewidth=1, color = color)
    plt.xlim([0, 1])
    plt.xlabel(r'$\overline{z}$', fontsize=14)
    plt.ylabel('$u_i$ (m/s)', fontsize=14)
    plt.gca().tick_params(labelsize=14)
    plt.title('Ion velocity')

def plot_electric_field_and_ionization_source(z_bar, Ez, S_iz, label_ez, label_siz, color_ez, color_siz, ax1=None, ax2=None):
    if ax1 is None or ax2 is None:  
        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()
    
    # trace electric field on the main axis
    line1, = ax1.plot(z_bar, Ez, color=color_ez, linewidth=1, label=label_ez)
    ax1.set_xlabel(r'$\overline{z}$', fontsize=14)
    ax1.set_ylabel('$E_z$ (V/m)', fontsize=14, color=color_ez)
    ax1.tick_params(axis='y', labelcolor=color_ez)
    
    # trace ionization
    line2, = ax2.plot(z_bar, S_iz, color=color_siz, linewidth=1, label=label_siz)
    ax2.set_ylabel('$S_{\\rm iz}$ (m$^3$/s$^{-1}$)', fontsize=14, color=color_siz)
    ax2.tick_params(axis='y', labelcolor=color_siz)
    
    # common title
    ax1.set_title('Electric field and ionization source', fontsize=16)
    
    ax1.legend(handles=[line1], loc='upper left', fontsize=12)  
    ax2.legend(handles=[line2], loc='upper right', fontsize=12) 

    return ax1, ax2 

def plot_magnetic_field(z_bar, BB):
    plt.plot(z_bar, BB * 1e4, 'b', linewidth=1)
    plt.xlim([0, 1])
    plt.xlabel(r'$\overline{z}$', fontsize=14)
    plt.ylabel('$B$ (Gauss)', fontsize=14)
    plt.gca().tick_params(labelsize=14)
    plt.title('Magnetic field')  

def plot_fluxes(z_bar, vng, vni, color_g, color_i):
    plt.plot(z_bar, vng, label='Neutrals', color = color_g)
    plt.plot(z_bar, vni, label='Ions', color = color_i)
    plt.xlim([0, 1])
    plt.xlabel(r'$\overline{z}$', fontsize=14)
    plt.ylabel('$\\Gamma$ (m$^{-2}$s$^{-1}$)', fontsize=14)
    plt.legend(fontsize=14)
    plt.title('Mass fluxes')
    plt.gca().tick_params(labelsize=14)

# Charger les données de chaque fichier
data1 = np.loadtxt("values1.csv", delimiter=',')
data_anna = np.loadtxt("values_anna.csv", delimiter=',')
# data_anna2 = np.loadtxt("values_anna_Te.csv", delimiter=',')

# Extraire les données X et Y

x, n_i, n_g, T_e, u_i, E_x, S_iz, BB, vng, vni = data1.T
z_bar, n_i_anna, n_g_anna, Te_end, u_i_anna, Ez, S_iz_anna, BB_anna, vng_anna, vni_anna = data_anna.T
# z_bar_array, T_bar_array = data_anna2.T

####### Comparaison des données #######

plt.clf()

plt.subplot(3, 2, 1)
plt.plot(z_bar, Te_end, linewidth=1, color = "blue")
plt.plot(x, T_e, linewidth=1, color = "red")
plt.xlim([0, 1])
plt.xlabel(r'$\overline{z}$', fontsize=14)
plt.ylabel('$T_e$ (V)', fontsize=14)
plt.gca().tick_params(labelsize=14)
plt.title('Electron temperature', fontsize=14)

plt.subplot(3, 2, 2)
plt.plot(z_bar, u_i_anna, linewidth=1, color = "blue")
plt.plot(x, u_i, linewidth=1, color = "red")
plt.xlim([0, 1])
plt.xlabel(r'$\overline{z}$', fontsize=14)
plt.ylabel('$u_i$ (m/s)', fontsize=14)
plt.gca().tick_params(labelsize=14)
plt.title('Ion velocity', fontsize=14)


plt.subplot(3, 2, 3)
plt.plot(z_bar, BB_anna * 1e4, 'b', linewidth=1)
plt.plot(x, BB * 1e4, 'r', linewidth=1)
plt.xlim([0, 1])
plt.xlabel(r'$\overline{z}$', fontsize=14)
plt.ylabel('$B$ (Gauss)', fontsize=14)
plt.gca().tick_params(labelsize=14)
plt.title('Magnetic field', fontsize=14)

plt.subplot(3, 2, 4)
plt.plot(z_bar, vng_anna, label='Neutrals', color = "green")
plt.plot(z_bar, vni_anna, label='Ions', color = "red")
plt.plot(x, vng, label='Neutrals', color = "orange")
plt.plot(x, vni, label='Ions', color = "blue")
plt.xlim([0, 1])
plt.xlabel(r'$\overline{z}$', fontsize=14)
plt.ylabel('$\\Gamma$ (m$^{-2}$s$^{-1}$)', fontsize=14)
plt.legend(fontsize=14)
plt.title('Mass fluxes', fontsize=14)
plt.gca().tick_params(labelsize=14)

ax1 = plt.subplot(3, 2, 5)
ax2 = ax1.twinx()  # Create a twin y-axis

color_i = "green"
color_g = "red"
# Tracé des densités ioniques sur l'axe principal
line1, = ax1.plot(z_bar, n_i_anna, color=color_i, linewidth=1, label="label_i")
ax1.set_xlabel(r'$\overline{z}$', fontsize=14)
ax1.set_ylabel('$n_i$ (m$^{-3}$)', fontsize=14, color=color_i)
ax1.tick_params(axis='y', labelcolor=color_i)

# Tracé des densités neutres sur l'axe secondaire
line2, = ax2.plot(z_bar, n_g_anna, color=color_g, linewidth=1, label="label_g")
ax2.set_ylabel('$n_g$ (m$^{-3}$)', fontsize=14, color=color_g)
ax2.tick_params(axis='y', labelcolor=color_g)

# Titre commun
ax1.set_title('Ion and Neutral Densities', fontsize=14)

ax1.legend(handles=[line1], loc='upper left', fontsize=12)  # Légende pour ax1
ax2.legend(handles=[line2], loc='upper right', fontsize=12)  # Légende pour ax2

color_i = "blue"
color_g = "orange"
# Tracé des densités ioniques sur l'axe principal
line1, = ax1.plot(x, n_i, color=color_i, linewidth=1, label="label_i")
ax1.set_xlabel(r'$\overline{z}$', fontsize=14)
ax1.set_ylabel('$n_i$ (m$^{-3}$)', fontsize=14, color=color_i)
ax1.tick_params(axis='y', labelcolor=color_i)

# Tracé des densités neutres sur l'axe secondaire
line2, = ax2.plot(x, n_g, color=color_g, linewidth=1, label="label_g")
ax2.set_ylabel('$n_g$ (m$^{-3}$)', fontsize=14, color=color_g)
ax2.tick_params(axis='y', labelcolor=color_g)

        
ax1 = plt.subplot(3, 2, 6)
ax2 = ax1.twinx()

# trace electric field on the main axis
line1, = ax1.plot(z_bar, Ez, color="blue", linewidth=1, label="label_ez")
ax1.set_xlabel(r'$\overline{z}$', fontsize=14)
ax1.set_ylabel('$E_z$ (V/m)', fontsize=14, color="blue")
ax1.tick_params(axis='y', labelcolor="blue")

# trace ionization
line2, = ax2.plot(z_bar, S_iz_anna, color="red", linewidth=1, label="label_siz")
ax2.set_ylabel('$S_{\\rm iz}$ (m$^3$/s$^{-1}$)', fontsize=14, color="red")
ax2.tick_params(axis='y', labelcolor="red")

# common title
ax1.set_title('Electric field and ionization source', fontsize=14)

ax1.legend(handles=[line1], loc='upper left', fontsize=12)  
ax2.legend(handles=[line2], loc='upper right', fontsize=12)

# trace electric field on the main axis
line1, = ax1.plot(x, E_x, color="green", linewidth=1, label="label_ez")
ax1.set_xlabel(r'$\overline{z}$', fontsize=14)
ax1.set_ylabel('$E_z$ (V/m)', fontsize=14, color="blue")
ax1.tick_params(axis='y', labelcolor="blue")

# trace ionization
line2, = ax2.plot(x, S_iz, color="orange", linewidth=1, label="label_siz")
ax2.set_ylabel('$S_{\\rm iz}$ (m$^3$/s$^{-1}$)', fontsize=14, color="red")
ax2.tick_params(axis='y', labelcolor="red")



# Adjust layout to avoid overlapping
plt.tight_layout()

plt.show()





# plt.figure(figsize=(12, 10))

# #Tracer les deux courbes ensemble
# fig, ax1 = plt.subplots()
# ax2 = ax1.twinx()

# # Tracé des premières données
# plt.subplot(3, 2, 1)
# plot_densities(z_bar, n_i_anna, n_g_anna, 'n_i_anna', 'n_g_anna', 'blue', 'red', ax1=ax1, ax2=ax2)

# # Tracé des deuxièmes données (sur les mêmes axes)
# plot_densities(x, n_i, n_g, 'n_i', 'n_g', 'red', 'orange', ax1=ax1, ax2=ax2)

# # Affichage final
# plt.xlim([0, 1])  # Limites communes pour l'axe x
# plt.tight_layout()  # Ajustement automatique des espaces
# plt.show()

# ################
# # plt.figure()
# #plot_electron_temperature(z_bar_array, T_bar_array, 'blue')
# plt.subplot(3, 2, 2)
# plot_electron_temperature(z_bar, Te_end, 'blue')
# plot_electron_temperature(x, T_e, 'red')
# plt.legend(["Anna", "v1"])

# # ################# Siz et Ez
# fig, ax1 = plt.subplots()
# ax2 = ax1.twinx()

# plt.subplot(3, 2, 3)
# # Tracé des premières données
# plot_electric_field_and_ionization_source(z_bar, Ez, S_iz_anna, 'Ez', 'S_iz_anna', 'blue', 'green', ax1=ax1, ax2=ax2)

# # Tracé des deuxièmes données (sur les mêmes axes)
# plot_electric_field_and_ionization_source(x, E_x, S_iz, 'E_x', 'S_iz_anna', 'red', 'orange', ax1=ax1, ax2=ax2)

# # Affichage final
# plt.xlim([0, 1])  # Limites communes pour l'axe x
# plt.tight_layout()  # Ajustement automatique des espaces
# plt.show()

# plt.subplot(3, 2, 4)
# # plt.figure()
# plot_ion_velocity(z_bar, u_i_anna, color='blue')
# plot_ion_velocity(x, u_i, color='red')
# plt.legend(["Anna", "v1"])

# # plt.figure()
# # plot_magnetic_field(z_bar, BB_anna)
# # plot_magnetic_field(x, BB)
# # plt.legend(["Anna", "v1"])

# # plt.figure()
# plot_fluxes(z_bar, vng_anna, vni_anna, 'blue','green')
# plot_fluxes(x, vng, vni, 'red', 'orange')
# plt.legend(["vng_Anna", "vng_anna", "vng_v1", "vni_v1"])

# plt.show()