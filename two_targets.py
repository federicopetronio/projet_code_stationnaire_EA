import test
from test import test_run
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter


def run_two_targets(B_max_values, Q_mgs_values, beta_mag, I_bar, thruster, propellant, var1, val1, var2, val2):
    
    var_index = {"ISP": 0, "Thrust": 1}

    n = len(B_max_values)
    m = len(Q_mgs_values)
    print("Run two targets")

    results = [[0 for j in range(m)] for i in range(n)]
    
    for i in range(n):
        for j in range(m):
            # test_run est supposé être une fonction calculant [ISP, Thrust]
            results[i][j] = (test_run(B_max_values[i]/10000, beta_mag, Q_mgs_values[j], I_bar, thruster, propellant)) # /10000 pour convertir en tesla

    print(results)

    M_A = np.array(results)[:, :, var_index[var1]]  # Matrice de A
    M_B = np.array(results)[:, :, var_index[var2]]  # Matrice de B

    Q_optimal, B_optimal = find_isobars_intersection(M_A, M_B, Q_mgs_values, B_max_values, val1, val2)
    print(f"Point d'intersection: Q = {Q_optimal}, B = {B_optimal}")

    # Trouver les indices du point d'intersection
    i_closest = np.argmin(np.abs(B_max_values - B_optimal))
    j_closest = np.argmin(np.abs(Q_mgs_values - Q_optimal))

    # Tracer la carte de couleur pour A
    plt.figure(figsize=(12, 5))

    # Carte de couleur pour A (ISP)
    plt.subplot(1, 2, 1)
    im_A = plt.imshow(M_A, cmap='viridis', origin='lower', aspect='auto', 
               extent=[min(Q_mgs_values), max(Q_mgs_values), min(B_max_values)/10000, max(B_max_values)/10000])
    plt.colorbar(im_A, label=var1)
    plt.title('-')
    plt.xlabel('Q_mgs (mg/s)')
    plt.ylabel('B_max (T)')

    # Ajout d'une courbe isobare pour A    
    isobar_value_A = val1
    Q_mgs_grid, B_max_grid = np.meshgrid(Q_mgs_values, B_max_values)
    CS_A = plt.contour(Q_mgs_grid, B_max_grid/10000, M_A, levels=[isobar_value_A], colors='red')
    plt.clabel(CS_A, inline=True, fontsize=10, fmt='Constraint 1')

    # Ajouter un point à la position i_closest, j_closest
    plt.plot(Q_mgs_values[j_closest], B_max_values[i_closest]/10000, 'go', markersize=10, label='Solution')
    plt.legend()

    # Carte de couleur pour B
    plt.subplot(1, 2, 2)
    im_B = plt.imshow(M_B, cmap='plasma', origin='lower', aspect='auto', 
               extent=[min(Q_mgs_values), max(Q_mgs_values), min(B_max_values)/10000, max(B_max_values)/10000])
    plt.colorbar(im_B, label=var2)
    plt.title('-')
    plt.xlabel('Q_mgs (mg/s)')
    plt.ylabel('B_max (T)')

    # Ajout d'une courbe isobare pour B
    isobar_value_B = val2
    Q_mgs_grid, B_max_grid = np.meshgrid(Q_mgs_values, B_max_values)
    CS_B = plt.contour(Q_mgs_grid, B_max_grid/10000, M_B, levels=[isobar_value_B], colors='blue')
    plt.clabel(CS_B, inline=True, fontsize=10, fmt='Constraint 2')

    # Ajouter un point à la position i_closest, j_closest
    plt.plot(Q_mgs_values[j_closest], B_max_values[i_closest]/10000, 'go', markersize=10, label='Solution')
    plt.legend()

    plt.tight_layout()
    plt.show()

    return B, Q



from scipy.interpolate import griddata

def find_isobars_intersection(M_A, M_B, Q_mgs_values, B_max_values, target_A, target_B):
    # Créer une grille pour l'interpolation
    Q_mgs_grid, B_max_grid = np.meshgrid(Q_mgs_values, B_max_values)

    # Trouver les valeurs des isobares pour A et B
    isobar_A = griddata((Q_mgs_grid.flatten(), B_max_grid.flatten()), M_A.flatten(), (Q_mgs_grid, B_max_grid), method='linear')
    isobar_B = griddata((Q_mgs_grid.flatten(), B_max_grid.flatten()), M_B.flatten(), (Q_mgs_grid, B_max_grid), method='linear')

    # Trouver les indices des isobares
    diff_A = np.abs(isobar_A - float(target_A))
    diff_B = np.abs(isobar_B - float(target_B))
    
    # Créer une matrice des différences
    total_diff = diff_A + diff_B

    # Trouver les indices du minimum dans la matrice des différences
    min_index = np.unravel_index(np.argmin(total_diff), total_diff.shape)

    # Obtenir les coordonnées des isobares
    Q_optimal = Q_mgs_grid[min_index]
    B_optimal = B_max_grid[min_index]

    return Q_optimal, B_optimal




 
