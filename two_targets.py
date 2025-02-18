from Calculations import test_run
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter


def run_two_targets(B_max_values, Q_mgs_values, magProfile, I_bar, thruster, propellant, var1, val1, var2, val2):
    
    var_index = {"ISP": 0, "Thrust": 1, "Thrust_power": 2, "Mass_utilization": 3, "Thrust_to_power_mN_kW": 4, "Total_efficiency": 5, "Elec_efficiency": 6}

    n = len(B_max_values)
    m = len(Q_mgs_values)
    print("Run two targets")

    results = [[0 for j in range(m)] for i in range(n)]
    
    magProfileSignature = '_'.join([str(magProfile(x)) for x in [0, 0.25, 0.5, 0.75, 1]])

    for i in range(n):
        for j in range(m):
            # Check if results are already in the file
            with open('results.txt', 'r') as file:
                lines = file.readlines()
                found = False
                for line in lines:
                    t, p, i_bar, q, b, beta, result = line.strip().split('| ')
                    if [t, p, float(i_bar), float(q), float(b), beta] == [thruster, propellant, I_bar, Q_mgs_values[j], B_max_values[i], magProfileSignature]:
                        results[i][j] = eval(result)
                        found = True
                        print("Déjà calculé")
                        break

            if not found:
                try:
                    results[i][j] = test_run(B_max_values[i]/10000, magProfile, Q_mgs_values[j], I_bar, thruster, propellant)
                except Exception as e:
                    results[i][j] = [0, 0, 0, 0, 0, 0, 0]
                    print("Erreur, ", e)
                print("Calcul")
                with open('results.txt', 'a') as file:
                    file.write(f"{thruster}| {propellant}| {I_bar}| {Q_mgs_values[j]}| {B_max_values[i]}| {magProfileSignature}| {results[i][j]}\n")

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

    return Q_optimal, B_optimal



import numpy as np
from scipy.interpolate import griddata
from scipy.optimize import fsolve

def find_isobars_intersection(M_A, M_B, Q_mgs_values, B_max_values, target_A, target_B):
    # Create a grid for interpolation
    Q_mgs_grid, B_max_grid = np.meshgrid(Q_mgs_values, B_max_values)

    # Define the function to find the intersection
    def intersection_func(x):
        Q, B = x
        A_value = griddata((Q_mgs_grid.flatten(), B_max_grid.flatten()), M_A.flatten(), (Q, B), method='linear')
        B_value = griddata((Q_mgs_grid.flatten(), B_max_grid.flatten()), M_B.flatten(), (Q, B), method='linear')
        
        # Handle potential None values
        if A_value is None or B_value is None:
            return [float('inf'), float('inf')]  # Return a large number if not found
        
        return [float(A_value) - float(target_A), float(B_value) - float(target_B)]

    # Initial guess: take a point near the grid center or a known approximate intersection
    initial_guess = [np.mean(Q_mgs_values), np.mean(B_max_values)]

    # Use fsolve to find the intersection point
    Q_optimal, B_optimal = fsolve(intersection_func, initial_guess)

    return Q_optimal, B_optimal

# Example usage (with your actual data for M_A, M_B, Q_mgs_values, and B_max_values)
# Q, B = find_isobars_intersection(M_A, M_B, Q_mgs_values, B_max_values, target_A, target_B)


 
