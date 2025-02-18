import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
from scipy.ndimage import gaussian_filter
from Calculations import test_run

global axs1, fig1, line_active_B, line_active_Q

# Function for horizontal line (Fixer B)
def on_move_B(event):
    if event.inaxes and line_active_B:
        for ax in axs1:
            # Get the current y-axis limits
            ymin, ymax = ax.get_ylim()
            # Remove any previous horizontal lines
            for line in ax.get_lines():
                line.set_visible(False)
            # Draw a new horizontal line at the current cursor position
            ax.axhline(event.ydata, color='red', linestyle='--', linewidth=1)
            # Ensure the y-limits remain fixed
            ax.set_ylim(ymin, ymax)
        fig1.canvas.draw_idle()  # Update the figure canvas

# Function for vertical line (Fixer Q)
def on_move_Q(event):
    if event.inaxes and line_active_Q:
        for ax in axs1:
            # Get the current x-axis limits
            xmin, xmax = ax.get_xlim()
            # Remove any previous vertical lines
            for line in ax.get_lines():
                line.set_visible(False)
            # Draw a new vertical line at the current cursor position
            ax.axvline(event.xdata, color='blue', linestyle='--', linewidth=1)
            # Ensure the x-limits remain fixed
            ax.set_xlim(xmin, xmax)
        fig1.canvas.draw_idle()  # Update the figure canvas

def on_button_click_point(event):
    print("Point unique")
    global line_active_B
    global line_active_Q
    line_active_B = False
    line_active_Q = False

# Toggle button for "Fixer B"
def on_button_click_B(event):
    global line_active_B
    global line_active_Q
    if line_active_B: line_active_B = False
    if not line_active_B: line_active_B, line_active_Q = True, False
    
    if line_active_B:
        fig1.canvas.mpl_connect('motion_notify_event', on_move_B)  # Track cursor movement for horizontal line
    else:
        for ax in axs1:
            # Remove all horizontal lines when deactivating
            for line in ax.get_lines():
                line.set_visible(False)
        fig1.canvas.draw_idle()  # Update the figure canvas

# Toggle button for "Fixer Q"
def on_button_click_Q(event):
    global line_active_Q
    global line_active_B
    if line_active_Q: line_active_Q = False
    if not line_active_Q: line_active_Q, line_active_B = True, False

    if line_active_Q:
        fig1.canvas.mpl_connect('motion_notify_event', on_move_Q)  # Track cursor movement for vertical line
    else:
        for ax in axs1:
            # Remove all vertical lines when deactivating
            for line in ax.get_lines():
                line.set_visible(False)
        fig1.canvas.draw_idle()  # Update the figure canvas

def run_one_target(B_max_values, Q_mgs_values, magProfile, I_bar, thruster, propellant, voltage, var1=None, val1=None):
    global axs1, fig1, line_active_B, line_active_Q
    line_active_B = False
    line_active_Q = False

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
                try:
                    results[i][j] = test_run(B_max_values[i]/10000, magProfile, Q_mgs_values[j], I_bar, thruster, propellant)
                except Exception as e:
                    results[i][j] = [0, 0, 0, 0, 0, 0, 0, 0]
                    print("Erreur, ", e)
                print("Calcul")
                with open('results.txt', 'a') as file:
                    file.write(f"{thruster}| {propellant}| {I_bar}| {Q_mgs_values[j]}| {B_max_values[i]}| {magProfileSignature}| {results[i][j]}\n")


    # Convert results to numpy arrays
    print(results)
    data = np.array(results)
    A, B, C, D, E, F, G, H = [data[:, :, i] for i in range(8)]
    matrixes = [A, B, C, D, E, F, G, H]

    titles = ['ISP', 'Thrust', 'Thrust_power', 'Mass_utilization', 'Thrust_to_power_mN_kW', 'Total_efficiency', 'Elec_efficiency', 'Voltage']
    units = ['s', 'mN', 'W', '', 'mN/kW', '', '', 'V']

    # First window: Overview of all matrices
    fig1, axs1 = plt.subplots(2, 5, figsize=(15, 6), gridspec_kw={'width_ratios': [1, 1, 1, 1, 0.8]})

    # Désactivation des axes de la colonne réservée aux boutons
    axs1[0,4].axis('off')
    axs1[1,4].axis('off')

    axs1 = axs1.flatten()

    saut = 0 # Pour sauter la case vide

    for idx, matrix in enumerate(matrixes):
        if idx == 4 and saut == 0: saut = 1
        ax = axs1[idx+saut]
        im = ax.imshow(matrix, aspect='auto', origin='lower',
                       extent=[Q_mgs_values[0], Q_mgs_values[-1], B_max_values[0], B_max_values[-1]])
        ax.set_title(f"{titles[idx]} ({units[idx]})")
        ax.set_xlabel('Q_mgs_values')
        ax.set_ylabel('B_max_values')
        plt.colorbar(im, ax=ax)

        global isobar_coords

        if var1 and titles[idx] == var1:
            isobar_coords = []
            CS = ax.contour(matrix, levels=[val1], colors='white', extent=[Q_mgs_values[0], Q_mgs_values[-1], B_max_values[0], B_max_values[-1]])
            ax.clabel(CS, inline=1, fontsize=10, fmt='%1.1f')
            for path in CS.collections[0].get_paths():
                isobar_coords.extend(path.vertices)

        if voltage and titles[idx] == 'Voltage':
            isobar_coords = []
            CS = ax.contour(matrix, levels=[voltage], colors='yellow', extent=[Q_mgs_values[0], Q_mgs_values[-1], B_max_values[0], B_max_values[-1]])
            ax.clabel(CS, inline=1, fontsize=10, fmt='%1.1fV')
            
            for path in CS.collections[0].get_paths():
                isobar_coords.extend(path.vertices)

    # Handle click events for interactivity
    def on_click(event):
        if event.inaxes:
            x, y = event.xdata, event.ydata
            print(f"Clicked at x={x}, y={y} in plot '{event.inaxes.title.get_text()}'")
            if line_active_B:
                print("Cutting horizontal line at y = ", y)
                curve_coords = [(x, y) for x in Q_mgs_values]
                coupe(curve_coords, Q_mgs_values, matrixes, absissa_title='Mass flow rate (mg/s)', type='horizontal')

            elif line_active_Q:
                print("Cutting vertical line at x = ", x)
                curve_coords = [(x, y) for y in B_max_values]
                coupe(curve_coords, B_max_values, matrixes, absissa_title='Magnetic field (gauss)', type='vertical')

            else:
                search(x, y)

    def search(x, y):
        if x < Q_mgs_values[0] or x > Q_mgs_values[-1] or y < B_max_values[0] or y > B_max_values[-1]:
            print("Click coordinates are out of bounds.")
            return
        print(f"Function launched with coordinates: Q={x} mg/s, B_max={y} Gauss.")
        plt.figure(figsize=(12, 10))
        test_run(y/10000, magProfile, x, I_bar, thruster, propellant, plotting=True)
        plt.show()

    def coupe(curve_coords, result_abscissa, matrixes, absissa_title='Curvilinear abscissa', type='isobar'):
        """
            Plot the variation of the results along the curve defined by curve_coords.
            The absiciissa of the curve is given by result_abscissa.
            EX: Vertical cut: curve_coords: liste de points verticaux, result_abscissa: Q_mgs_values
        """
        fig3, axs3 = plt.subplots(4, 2, figsize=(12, 10))
        axs3 = axs3.flatten()

        for idx, matrix in enumerate(matrixes):
            ax = axs3[idx]
            values_along_isobar = [
                matrix[int((y - B_max_values[0]) / (B_max_values[-1] - B_max_values[0]) * (matrix.shape[0] - 1)),
                    int((x - Q_mgs_values[0]) / (Q_mgs_values[-1] - Q_mgs_values[0]) * (matrix.shape[1] - 1))] 
                for x, y in curve_coords
            ]
            ax.plot(result_abscissa, values_along_isobar, label=titles[idx])
            ax.set_title(f'Variation of {titles[idx]} along the {type}')
            ax.set_xlabel(absissa_title)
            ax.legend()

        plt.tight_layout(pad=2.0)
        plt.show()


    def on_button_click_iso(event):
        global isobar_coords
        curvilinear_abscissa = np.linspace(0, 1, len(isobar_coords))
        coupe(isobar_coords, curvilinear_abscissa, matrixes)

    def on_move(event):
        if event.inaxes:  # Vérifie que la souris est bien dans un subplot
            x, y = event.xdata, event.ydata
            
            # Calcul des indices correspondants dans les matrices
            x_idx = int((x - Q_mgs_values[0]) / (Q_mgs_values[-1] - Q_mgs_values[0]) * (matrixes[0].shape[1] - 1))
            y_idx = int((y - B_max_values[0]) / (B_max_values[-1] - B_max_values[0]) * (matrixes[0].shape[0] - 1))

            # Clamping des indices pour éviter des erreurs hors limites
            x_idx = max(0, min(x_idx, matrixes[0].shape[1] - 1))
            y_idx = max(0, min(y_idx, matrixes[0].shape[0] - 1))

            # Mise à jour de tous les subplots avec la valeur correspondante
            saut = 0
            for idx, matrix in enumerate(matrixes):
                if idx == 4 and saut == 0: saut = 1
                value = matrix[y_idx, x_idx]  # Récupère la valeur de la matrice
                axs1[idx+saut].set_title(f"{titles[idx]} ({units[idx]})\nValue: {value:.2f}")

            fig1.canvas.draw_idle()  # Met à jour la figure

    # Connexion de l'événement
    fig1.canvas.mpl_connect('motion_notify_event', on_move)


    # Ajustement de la mise en page pour laisser de la place à droite
    # Ici, les sous-graphiques occupent de 5% à 75% de la largeur, et l'espace entre 75% et 100% est réservé aux boutons.
    fig1.subplots_adjust(left=0.05, right=0.75, top=0.9, bottom=0.1)

    # Création des boutons dans la zone libre à droite (coordonnées en fraction de la figure)
    ax_button_B = fig1.add_axes([0.85, 0.7, 0.12, 0.07])
    button_B = Button(ax_button_B, 'Coupe à B fixe')
    button_B.on_clicked(on_button_click_B)

    ax_button_Q = fig1.add_axes([0.85, 0.6, 0.12, 0.07])
    button_Q = Button(ax_button_Q, 'Coupe à Q fixe')
    button_Q.on_clicked(on_button_click_Q)

    ax_button_point = fig1.add_axes([0.85, 0.5, 0.12, 0.07])
    button_point = Button(ax_button_point, 'Analyser un point')
    button_point.on_clicked(on_button_click_point)

    ax_button_iso = fig1.add_axes([0.85, 0.4, 0.12, 0.07])
    button_iso = Button(ax_button_iso, "Coupe sur l'isobare")
    button_iso.on_clicked(on_button_click_iso)

    #fig1.delaxes(axs1[-1])  # Remove the last empty subplot
    plt.tight_layout(pad=2.0)

    fig1.canvas.mpl_connect("button_press_event", on_click)
    plt.show()