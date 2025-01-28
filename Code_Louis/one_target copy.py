import physics
from Calculations import test_run
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

def run_one_target(B_max_values, Q_mgs_values, magProfile, I_bar, thruster, propellant, var1=None, val1=None):
    n = len(B_max_values)
    m = len(Q_mgs_values)
    print("One target")

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
                results[i][j] = test_run(B_max_values[i]/10000, magProfile, Q_mgs_values[j], I_bar, thruster, propellant)
                print("Calcul")
                with open('results.txt', 'a') as file:
                    file.write(f"{thruster}| {propellant}| {I_bar}| {Q_mgs_values[j]}| {B_max_values[i]}| {magProfileSignature}| {results[i][j]}\n")


    # Convert results to numpy arrays
    data = np.array(results)
    A, B, C, D, E, F, G = [data[:, :, i] for i in range(7)]
    matrixes = [A, B, C, D, E, F, G]

    titles = ['ISP', 'Thrust', 'Thrust_power', 'Mass_utilization', 'Thrust_to_power_mN_kW', 'Total_efficiency', 'Elec_efficiency']
    units = ['s', 'mN', 'W', '', 'mN/kW', '', '']

    # Adjust for min/max constraints
    if var1 == "min":
        val1 = np.min(data[:, :, titles.index(var1)]) * 1.1
        print("Min: ", val1)
    elif var1 == "max":
        val1 = np.max(data[:, :, titles.index(var1)]) * 0.9
        print("Max: ", val1)

    # First window: Overview of all matrices
    fig1, axs1 = plt.subplots(2, 4, figsize=(12, 5))
    axs1 = axs1.flatten()

    for idx, matrix in enumerate(matrixes):
        ax = axs1[idx]
        im = ax.imshow(matrix, aspect='auto', origin='lower',
                       extent=[Q_mgs_values[0], Q_mgs_values[-1], B_max_values[0], B_max_values[-1]])
        ax.set_title(f"{titles[idx]} ({units[idx]})")
        ax.set_xlabel('Q_mgs_values')
        ax.set_ylabel('B_max_values')
        plt.colorbar(im, ax=ax)

        if var1 and titles[idx] == var1:
            CS = ax.contour(matrix, levels=[val1], colors='white',
                            extent=[Q_mgs_values[0], Q_mgs_values[-1], B_max_values[0], B_max_values[-1]])
            ax.clabel(CS, inline=1, fontsize=10, fmt='%1.1f')

    fig1.delaxes(axs1[-1])  # Remove empty subplot
    plt.tight_layout()

    # Second window: Detailed variation along the isobar
    if var1 and len(CS.allsegs[0]) > 1:
        fig2, axs2 = plt.subplots(4, 2, figsize=(12, 10))
        axs2 = axs2.flatten()
        isobar_coords = CS.allsegs[0][0]
        curvilinear_abscissa = np.linspace(0, 1, len(isobar_coords))

        for idx, matrix in enumerate(matrixes):
            ax = axs2[idx]
            values_along_isobar = [
                matrix[int((y - B_max_values[0]) / (B_max_values[-1] - B_max_values[0]) * (matrix.shape[0] - 1)),
                       int((x - Q_mgs_values[0]) / (Q_mgs_values[-1] - Q_mgs_values[0]) * (matrix.shape[1] - 1))]
                for x, y in isobar_coords
            ]
            ax.plot(curvilinear_abscissa, values_along_isobar, label=titles[idx])
            ax.set_title(f'Variation of {titles[idx]} along the isobar')
            ax.set_xlabel('Curvilinear abscissa')
            ax.set_ylabel('Values')
            ax.legend()

        for ax in axs2[len(matrixes):]:  # Hide unused subplots
            fig2.delaxes(ax)

        plt.tight_layout()

    # Handle click events
    def on_click(event):
        if event.inaxes:
            x, y = event.xdata, event.ydata
            print(f"Clicked at x={x}, y={y} in plot '{event.inaxes.title.get_text()}'")
            search(x, y)

    def search(x, y):
        if x < Q_mgs_values[0] or x > Q_mgs_values[-1] or y < B_max_values[0] or y > B_max_values[-1]:
            print("Click coordinates are out of bounds.")
            return
        print(f"Function launched with coordinates: Q={x} mg/s, B_max={y} Gauss.")
        plt.figure(figsize=(12, 10))
        test_run(y/10000, magProfile, x, I_bar, thruster, propellant, plotting=True)
        plt.show()

    fig1.canvas.mpl_connect("button_press_event", on_click)

    plt.show()
