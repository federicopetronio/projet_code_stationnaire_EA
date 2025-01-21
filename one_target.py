import physics
from physics import test_run
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

def run_one_target(B_max_values, Q_mgs_values, beta_mag, I_bar, thruster, propellant, var1=None, val1=None):
    
    n = len(B_max_values)
    m = len(Q_mgs_values)
    print("One target")

    results = [[0 for j in range(m)] for i in range(n)]
    
    for i in range(n):
        for j in range(m):
            # test_run est supposé être une fonction calculant [ISP, Thrust]
            results[i][j] = (test_run(B_max_values[i]/10000, beta_mag, Q_mgs_values[j], I_bar, thruster, propellant)) # /10000 pour convertir en tesla

    print(results)

    # Exemple de matrice indexée par i, j contenant des listes [A, B]
    data = np.array(results)

    # Séparer les valeurs de A (ISP) et B (Thrust)
    A = data[:, :, 0]  # Matrice de A
    B = data[:, :, 1]  # Matrice de B

    sigma = 1  # Écart type du filtre gaussien
    A_smoothed = gaussian_filter(A, sigma=sigma)
    B_smoothed = gaussian_filter(B, sigma=sigma)

    # Tracer la carte de couleur pour A
    plt.figure(figsize=(12, 5))

    # Carte de couleur pour A (ISP)
    plt.subplot(1, 2, 1)
    im_A = plt.imshow(A_smoothed, cmap='viridis', origin='lower', aspect='auto', 
               extent=[min(Q_mgs_values), max(Q_mgs_values), min(B_max_values)/10000, max(B_max_values)/10000])
    plt.colorbar(im_A, label='ISP (s)')
    plt.title('Specific impulse')
    plt.xlabel('Q_mgs (mg/s)')
    plt.ylabel('B_max (T)')

    # Ajout d'une courbe isobare pour A (ISP)

    if var1 == "ISP" and val1 is not None:
        isobar_value_A = val1
        Q_mgs_grid, B_max_grid = np.meshgrid(Q_mgs_values, B_max_values)
        CS_A = plt.contour(Q_mgs_grid, B_max_grid/10000, A_smoothed, levels=[isobar_value_A], colors='red')
        plt.clabel(CS_A, inline=True, fontsize=10, fmt=f'ISP = {isobar_value_A} s')

    # Carte de couleur pour B (Thrust)
    plt.subplot(1, 2, 2)
    im_B = plt.imshow(B_smoothed, cmap='plasma', origin='lower', aspect='auto', 
               extent=[min(Q_mgs_values), max(Q_mgs_values), min(B_max_values)/10000, max(B_max_values)/10000])
    plt.colorbar(im_B, label='Thrust (mN)')
    plt.title('Thrust')
    plt.xlabel('Q_mgs (mg/s)')
    plt.ylabel('B_max (T)')

    # Ajout d'une courbe isobare pour B (Thrust)
    if var1 == "Thrust" and val1 is not None:
        isobar_value_B = val1
        Q_mgs_grid, B_max_grid = np.meshgrid(Q_mgs_values, B_max_values)
        CS_B = plt.contour(Q_mgs_grid, B_max_grid/10000, B_smoothed, levels=[isobar_value_B], colors='blue')
        plt.clabel(CS_B, inline=True, fontsize=10, fmt=f'Thrust = {isobar_value_B} mN')

    plt.tight_layout()
    plt.show()
