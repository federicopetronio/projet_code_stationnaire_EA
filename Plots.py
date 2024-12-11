import numpy as np
import subprocess
import matplotlib.pyplot as plt
import sys

################ Plot functions ################

def plot_densities(z_bar, n_i, n_g, Gammab, Gb, label_i, label_g, color_i, color_g,ax1=None, ax2=None):
    if ax1 is None or ax2 is None:  # Si les axes ne sont pas fournis, créez-en
        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()
    
    # Tracé des densités ioniques sur l'axe principal
    line1, = ax1.plot(z_bar, n_i(Gammab, Gb), color=color_i, linewidth=1, label=label_i)
    ax1.set_xlabel('$z_{bar}$', fontsize=14)
    ax1.set_ylabel('$n_i$ (m$^{-3}$)', fontsize=14, color=color_i)
    ax1.tick_params(axis='y', labelcolor=color_i)
    
    # Tracé des densités neutres sur l'axe secondaire
    line2, = ax2.plot(z_bar, n_g(Gammab), color=color_g, linewidth=1, label=label_g)
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
    plt.xlabel('$x/L_{\\rm ch}$', fontsize=14)
    plt.ylabel('$T_e$ (V)', fontsize=14)
    plt.gca().tick_params(labelsize=14)
    plt.title('Electron temperature')

def plot_ion_velocity(z_bar, u_i, Gammab, Gb, color):
    plt.plot(z_bar, u_i(Gammab, Gb), linewidth=1, color = color)
    plt.xlim([0, 1])
    plt.ylabel('$u_i$ (m/s)', fontsize=14)
    plt.gca().tick_params(labelsize=14)
    plt.title('Ion velocity')

def plot_electric_field_and_ionization_source(z_bar, Ez, S_iz, Gammab, Gb, label_ez, label_siz, color_ez, color_siz, ax1=None, ax2=None):
    if ax1 is None or ax2 is None:  
        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()
    
    # trace electric field on the main axis
    line1, = ax1.plot(z_bar, Ez(Gammab, Gb), color=color_ez, linewidth=1, label=label_ez)
    ax1.set_xlabel('$z_{bar}$', fontsize=14)
    ax1.set_ylabel('$E_x$ (V/m)', fontsize=14, color=color_ez)
    ax1.tick_params(axis='y', labelcolor=color_ez)
    
    # trace ionization
    line2, = ax2.plot(z_bar, S_iz(Gammab, Gb), color=color_siz, linewidth=1, label=label_siz)
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
    plt.xlabel('$z/L_{\\rm ch}$', fontsize=14)
    plt.ylabel('$B$ (Gauss)', fontsize=14)
    plt.gca().tick_params(labelsize=14)
    plt.title('Magnetic field')  

def plot_fluxes(z_bar, vng, vni, color_g, color_i):
    plt.plot(z_bar, vng, label='Neutrals', color = color_g)
    plt.plot(z_bar, vni, label='Ions', color = color_i)
    plt.xlim([0, 1])
    plt.xlabel('$z_bar$', fontsize=14)
    plt.ylabel('$\\Gamma$ (m$^{-2}$s$^{-1}$)', fontsize=14)
    plt.legend(fontsize=14)
    plt.title('Mass fluxes')
    plt.gca().tick_params(labelsize=14)
