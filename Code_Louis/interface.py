import tkinter as tk
from tkinter import ttk
from tkinter import filedialog, messagebox
import numpy as np
import pandas as pd
import time
import threading
import one_target
from one_target import run_one_target
from two_targets import run_two_targets
#from physics import test_run
from Calculations import test_run
from RangeSlider.RangeSlider import RangeSliderH

# Create the main window
root = tk.Tk()
root.title("COMHET - Operational optimizer")
root.geometry("920x725")  # Increased height for new section

# Set a nice font and style for labels and entries
style = ttk.Style()
style.configure("TLabel", font=("Arial", 12))
style.configure("TEntry", font=("Arial", 12))
style.configure("TButton", font=("Arial", 12, "bold"), foreground="white", background="blue")
style.configure("TCombobox", font=("Arial", 12))

file_path = ""

def calculate_for_one_target(thrust, propellant, magProfile, var1, val1, B_max_min, B_max_max, Q_mgs_min, Q_mgs_max, precision):
    print(precision)
    B_length = int(precision[0])
    Q_length = int(precision[1])
    B_max_values = [B_max_min + ((B_max_max-B_max_min)/B_length)*i for i in range(B_length)]
    Q_mgs_values = [Q_mgs_min + ((Q_mgs_max-Q_mgs_min)/Q_length)*i for i in range(Q_length)]

    run_one_target(B_max_values, Q_mgs_values, magProfile, 1.34, thrust, propellant, var1, val1)
    return "Consultez les graphiques"

def calculate_for_two_targets(thrust, propellant, magProfile, var1, val1, var2, val2, B_max_min, B_max_max, Q_mgs_min, Q_mgs_max, precision):
    B_length = int(precision[0])
    Q_length = int(precision[1])
    B_max_values = [B_max_min + ((B_max_max-B_max_min)/B_length)*i for i in range(B_length)]
    Q_mgs_values = [Q_mgs_min + ((Q_mgs_max-Q_mgs_min)/Q_length)*i for i in range(Q_length)]

    B, Q = run_two_targets(B_max_values, Q_mgs_values, magProfile, 1.34, thrust, propellant, var1, val1, var2, val2)
    return B, Q

# Function to mimic a long-running calculation
def perform_calculation(thrust, propellant, profile_factor, b_max_min, b_max_max, q_min, q_max, precision, var1=None, val1=None, var2=None, val2=None): 
    # Show progress bar in the output frame
    output_frame.grid(row=4, column=0, columnspan=2)  # Show the output panel using grid
    output_b_value.grid_forget()  # Hide B and Q labels initially
    output_q_value.grid_forget()
    progress_bar.grid(row=0, column=0, columnspan=2)  # Show the progress bar using grid
    
    print(file_path)
    if file_path == "":
        magProfile = lambda x: np.exp(-float(profile_factor) * (x - 1)**2)
        print("gaussienne")
    else:
        magProfile = interpolate_csv(file_path)
        print("custom")

    print(magProfile)
    # Simulate calculation time
    progress_bar.start(10)  # Start progress bar animation

    # Perform calculation
    if var2 and val2:  # If both target variables are filled
        Q, B = calculate_for_two_targets(thrust, propellant, magProfile, var1, val1, var2, val2, b_max_min, b_max_max, q_min, q_max, precision)
        output_b_value.config(text=f"B: {int(B)} Gauss")  # Show calculated B and Q
        output_q_value.config(text=f"Q: {int(Q)} mg/s")
        output_b_value.grid(row=1, column=0)
        output_q_value.grid(row=1, column=1)
        test_run(B, profile_factor, Q, 1.34, thrust, propellant, plotting=False)

    else:  # Only the first target variable is filled
        calculate_for_one_target(thrust, propellant, magProfile, var1, val1, b_max_min, b_max_max, q_min, q_max, precision)

    # After calculation is done, display the result
    progress_bar.stop()  # Stop progress bar animation
    progress_bar.grid_forget()  # Hide the progress bar

# Function to be executed when 'Run' is clicked
def run_command():
    # Get the inputs from the interface
    thrust = thrust_dropdown.get()
    propellant = propellant_dropdown.get()
    profile_factor = profile_entry.get()
    var1 = variable_1_dropdown.get()
    val1 = value_1_entry.get()
    var2 = variable_2_dropdown.get()
    val2 = value_2_entry.get()

    # Get the B and Q min/max values from the scales
    b_min = b_range_slider.getValues()[0]
    b_max = b_range_slider.getValues()[1]
    q_min = q_range_slider.getValues()[0]
    q_max = q_range_slider.getValues()[1]
   
    precision = [int(precision_B_entry.get()), int(precision_Q_entry.get())]

    # Run the calculation in a separate thread
    calculation_thread = threading.Thread(
        target=perform_calculation, 
        args=(thrust, propellant, profile_factor, b_min, b_max, q_min, q_max, precision, var1, val1, var2, val2)
    )
    calculation_thread.start()

# Title Frame to ensure the title is at the very top
title_frame = tk.Frame(root)
title_frame.grid(row=0, column=0, columnspan=2, sticky="n", pady=(10, 0))  # sticky="n" aligns to the top

# Add a title label inside the title frame
title_label = tk.Label(title_frame, text="COMHET - Operational optimizer", font=("Arial", 24, "bold"))
title_label.pack()

# Engineering Inputs section
engineering_frame = tk.LabelFrame(root, text="Engineering Inputs", padx=20, pady=20, font=("Arial", 14, "bold"))
engineering_frame.grid(row=1, column=0, padx=20, pady=20, sticky="n")

# Dropdown for thrust input
thrust_label = ttk.Label(engineering_frame, text="Thruster:")
thrust_label.grid(row=0, column=0, sticky="w", pady=5)
thrust_options = ["PPSX00", "PPS1350", "PPS5000", "PPSX000_Hi", "PPS20k"]
thrust_dropdown = ttk.Combobox(engineering_frame, values=thrust_options)
thrust_dropdown.grid(row=0, column=1, pady=5)

def modify_ranges(event):
    b_ranges = {"": (50, 500), "PPSX00": (60, 200), "PPS1350": (50, 500), "PPS5000": (50, 500), "PPSX000_Hi": (50, 500), "PPS20k": (50, 500)}
    q_ranges = {"": (1, 20), "PPSX00": (2, 10), "PPS1350": (1, 20), "PPS5000": (1, 20), "PPSX000_Hi": (1, 20), "PPS20k": (1, 20)}

    selected_thrust = thrust_dropdown.get()
    b_range = b_ranges.get(selected_thrust, (50, 500))
    q_range = q_ranges.get(selected_thrust, (1, 20))

    b_range_slider.forceValues([b_range[0], b_range[1]])
    q_range_slider.forceValues([q_range[0], q_range[1]])


thrust_dropdown.bind("<<ComboboxSelected>>", modify_ranges)

# Dropdown for propellant input
propellant_label = ttk.Label(engineering_frame, text="Propellant:")
propellant_label.grid(row=1, column=0, sticky="w", pady=5)
propellant_options = ["xenon", "krypton"]
propellant_dropdown = ttk.Combobox(engineering_frame, values=propellant_options)
propellant_dropdown.grid(row=1, column=1, pady=5)

def choose_profile():
    # Enable or disable the appropriate input fields
    if profile_var.get() == "Gaussian":
        profile_entry.config(state="normal")  # Enable manual entry
        upload_button.config(state="disabled")  # Disable upload button
    elif profile_var.get() == "Custom":
        profile_entry.config(state="disabled")  # Disable manual entry
        upload_button.config(state="normal")  # Enable upload button

def upload_csv():
    global file_path
    # Open a file dialog for the user to select a CSV file
    file_path = filedialog.askopenfilename(
        filetypes=[("CSV Files", "*.csv"), ("All Files", "*.*")]
    )
    if file_path:
        # You can handle the uploaded file here
        messagebox.showinfo("File Selected", f"Custom profile loaded: {file_path}")
    else:
        messagebox.showwarning("No File Selected", "Please select a valid CSV file.")

def interpolate_csv(file_path, tolerance=0.001):
    """
    Lit un fichier .csv contenant des valeurs x et y et renvoie une fonction polynomiale 
    qui interpole ces points avec une précision spécifiée (par défaut 0.1%).
    
    Args:
        file_path (str): Chemin vers le fichier CSV.
        tolerance (float): Tolérance relative (ex : 0.001 pour 0.1%).

    Returns:
        callable: Fonction polynomiale qui interpole les données.
    """
    # Lire le fichier CSV
    data = pd.read_csv(file_path)
    if "x" not in data.columns or "y" not in data.columns:
        raise ValueError("Le fichier doit contenir des colonnes 'x' et 'y'.")

    x = data['x'].values
    y = data['y'].values

    # Vérification que les données sont triées par x
    if not np.all(np.diff(x) > 0):
        raise ValueError("Les valeurs de x doivent être strictement croissantes.")

    # Interpolation avec ajustement du degré
    for degree in range(1, len(x)):  # Tester des degrés croissants
        coefficients = np.polyfit(x, y, degree)  # Ajuste un polynôme de degré 'degree'
        polynomial = np.poly1d(coefficients)    # Crée une fonction polynomiale
        interpolated_y = polynomial(x)          # Calcule les y interpolés
        relative_error = np.abs((y - interpolated_y) / y)  # Erreur relative

        if np.all(relative_error <= tolerance):  # Si toutes les erreurs respectent la tolérance
            # Interpolation
            import matplotlib.pyplot as plt

            # Plot the raw data points
            plt.scatter(x, y, color='red', label='Raw Data')

            # Plot the interpolated polynomial
            x_dense = np.linspace(min(x), max(x), 500)  # Create a dense range of x values
            y_dense = polynomial(x_dense)  # Evaluate the polynomial at these x values
            plt.plot(x_dense, y_dense, label=f'Polynomial Degree {degree}')

            plt.xlabel('x')
            plt.ylabel('y')
            plt.title('Interpolation of CSV Data')
            plt.legend()
            plt.grid(True)
            plt.show()
            return polynomial

    # Si aucun polynôme ne satisfait la tolérance
    raise RuntimeError("Impossible d'interpoler avec la précision demandée.")

# Create a frame for the engineering options
profile_frame = ttk.Frame(engineering_frame, padding="10")
profile_frame.grid(row=2, column=0)

# Label for Profile Factor
profile_label = ttk.Label(profile_frame, text="Profile Factor:")
profile_label.grid(row=0, column=0, sticky="w")

# Radio buttons to select profile type
profile_var = tk.StringVar(value="Gaussian")  # Default selection

gaussian_button = ttk.Radiobutton(
    profile_frame, text="Gaussian Profile", variable=profile_var, value="Gaussian", command=choose_profile
)
gaussian_button.grid(row=1, column=0, sticky="w", pady=5)

custom_button = ttk.Radiobutton(
    profile_frame, text="Custom Profile", variable=profile_var, value="Custom", command=choose_profile
)
custom_button.grid(row=2, column=0, sticky="w", pady=5)

# Entry for Gaussian Profile Factor
profile_entry = ttk.Entry(profile_frame, state="normal")  # Default enabled
profile_entry.grid(row=1, column=1, pady=5)

# Button to upload a CSV file for Custom Profile
upload_button = ttk.Button(profile_frame, text="Upload CSV", state="disabled", command=upload_csv)
upload_button.grid(row=2, column=1, pady=5)

# Target section
target_frame = tk.LabelFrame(root, text="Target", padx=20, pady=20, font=("Arial", 14, "bold"))
target_frame.grid(row=1, column=1, padx=20, pady=20, sticky="n")

# First variable and value
variable_1_label = ttk.Label(target_frame, text="1) Variable:")
variable_1_label.grid(row=0, column=0, sticky="w", pady=5)
variable_1_options = ["ISP", "Thrust", "Thrust_power", "Mass_utilization", "Thrust_to_power_mN_kW", "Total_efficiency", "Elec_efficency"]
variable_1_dropdown = ttk.Combobox(target_frame, values=variable_1_options)
variable_1_dropdown.grid(row=0, column=1, pady=5)

value_1_label = ttk.Label(target_frame, text="Value: ")
value_1_label.grid(row=1, column=0, sticky="w", pady=5)
value_1_entry = ttk.Entry(target_frame)
value_1_entry.grid(row=1, column=1, pady=5)

# Second variable and value
variable_2_label = ttk.Label(target_frame, text="2) Variable:")
variable_2_label.grid(row=2, column=0, sticky="w", pady=5)
variable_2_dropdown = ttk.Combobox(target_frame, values=variable_1_options)
variable_2_dropdown.grid(row=2, column=1, pady=5)

value_2_label = ttk.Label(target_frame, text="Value: ")
value_2_label.grid(row=3, column=0, sticky="w", pady=5)
value_2_entry = ttk.Entry(target_frame)
value_2_entry.grid(row=3, column=1, pady=5)

# Ranges for B and Q section
range_frame = tk.LabelFrame(root, text="Ranges for B_max and Q", padx=20, pady=20, font=("Arial", 14, "bold"))
range_frame.grid(row=3, column=0, padx=20, pady=20, columnspan=2, sticky="n")

b_ranges = {"":(50,500), "PPSX00": (60, 200), "PPS1350": (50, 500), "PPS5000": (50, 500), "PPSX000_Hi": (50, 500), "PPS20k": (50, 500)}
q_ranges = {"":(1,20), "PPSX00": (2, 10), "PPS1350": (1, 20), "PPS5000": (1, 20), "PPSX000_Hi": (1, 20), "PPS20k": (1,20)}

b_range = b_ranges[thrust_dropdown.get()]
q_range = q_ranges[thrust_dropdown.get()]

# B_min and B_max scales
b_range_label = ttk.Label(range_frame, text="B range (gauss):")
b_range_label.grid(row=0, column=0, pady=0)
b_range_min = tk.DoubleVar(value = b_range[0])  #left handle variable initialised to value 0.2
b_range_max = tk.DoubleVar(value = b_range[1])  #right handle variable initialised to value 0.85
b_range_slider = RangeSliderH(range_frame , [b_range_min, b_range_max] , min_val = b_range[0], max_val = b_range[1], padX = 17, bgColor="#f0f0f0")   #horizontal slider, [padX] value might be needed to be different depending on system, font and handle size. Usually [padX] = 12 serves,
b_range_slider.grid(row=1, column=0, pady=0)

def on_b_range_change(event):
    print("B range changed")
    b_min, b_max = b_range_slider.getValues()
    b_range = b_ranges.get(thrust_dropdown.get(), (50, 500))
    
    if b_min < b_range[0]:
        b_range_slider.forceValues([b_range[0], b_max])
    elif b_max > b_range[1]:
        b_range_slider.forceValues([b_min, b_range[1]])

b_range_min.trace_add('write', on_b_range_change)

# Q_min and Q_max scales
q_range_label = ttk.Label(range_frame, text="Q range (mg/s):")
q_range_label.grid(row=2, column=0, pady=0)
q_range_min = tk.DoubleVar(value = q_range[0])  #left handle variable initialised to value 0.2
q_range_right = tk.DoubleVar(value = q_range[1])  #right handle variable initialised to value 0.85
q_range_slider = RangeSliderH( range_frame , [q_range_min, q_range_right], min_val = q_range[0],max_val = q_range[1], padX = 17, bgColor="#f0f0f0")   #horizontal slider, [padX] value might be needed to be different depending on system, font and handle size. Usually [padX] = 12 serves,
q_range_slider.grid( row=3, column=0, pady=0)

# Precision for B and Q text entries
precision_B_label = ttk.Label(range_frame, text="Precision B:")
precision_B_label.grid(row=4, column=0, pady=5)
precision_B_entry = ttk.Entry(range_frame, width=5)
precision_B_entry.grid(row=4, column=1, pady=5)

precision_Q_label = ttk.Label(range_frame, text="Precision Q:")
precision_Q_label.grid(row=5, column=0, pady=5)
precision_Q_entry = ttk.Entry(range_frame, width=5)
precision_Q_entry.grid(row=5, column=1, pady=5)

# Output section with B, Q values and progress bar
output_frame = tk.Frame(root)  # Frame that will hold either progress bar or results
output_frame.grid(row=4, column=0, columnspan=2, pady=10)  # Initially hidden

# Progress bar
progress_var = tk.DoubleVar()
progress_bar = ttk.Progressbar(output_frame, variable=progress_var, maximum=100, mode="indeterminate")

# Output B and Q labels, initially hidden
output_b_value = ttk.Label(output_frame, text="", font=("Arial", 14))
output_q_value = ttk.Label(output_frame, text="", font=("Arial", 14))

# Run button to initiate calculation
style = ttk.Style()
style.configure("Custom.TButton", foreground="white", background="blue")
run_button = ttk.Button(root, text="Run", style="Custom.TButton", command=run_command)
run_button.grid(row=0, column=2, columnspan=2, pady=20)

# Start the GUI event loop
root.mainloop()