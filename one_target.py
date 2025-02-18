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

def run_one_target(B_max_values, Q_mgs_values, magProfile, I_bar, thruster, propellant, var1=None, val1=None):
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

    # Create button for "Fixer B"
    ax_button_B = fig1.add_axes([0.85, 0.05, 0.1, 0.075])  # Position for button
    button_B = Button(ax_button_B, 'Fixer B')
    button_B.on_clicked(on_button_click_B)

    # Create button for "Fixer Q"
    ax_button_Q = fig1.add_axes([0.85, 0.15, 0.1, 0.075])  # Position for "Fixer Q" button
    button_Q = Button(ax_button_Q, 'Fixer Q')
    button_Q.on_clicked(on_button_click_Q)

    fig1.delaxes(axs1[-1])  # Remove the last empty subplot
    plt.tight_layout(pad=2.0)

    # Second window: Detailed variation along the isobar
    fig2, axs2 = plt.subplots(4, 2, figsize=(12, 10))
    axs2 = axs2.flatten()
    isobar_coords = np.array([[0.2, 0.4], [0.5, 0.6]])  # Example coordinates of the isobar
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

    plt.tight_layout(pad=2.0)

    # Handle click events for interactivity
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