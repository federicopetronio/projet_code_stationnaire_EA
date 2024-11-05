import tkinter as tk
from tkinter import ttk
import time
import threading
import one_target
from one_target import run_one_target
from two_targets import run_two_targets
# Create the main window
root = tk.Tk()
root.title("Engineering Inputs Interface")
root.geometry("800x800")  # Increased height for new section

# Set a nice font and style for labels and entries
style = ttk.Style()
style.configure("TLabel", font=("Arial", 12))
style.configure("TEntry", font=("Arial", 12))
style.configure("TButton", font=("Arial", 12, "bold"), foreground="white", background="#4CAF50")
style.configure("TCombobox", font=("Arial", 12))

loading_texts = ["Calculating", "Calculating.", "Calculating..", "Calculating..."]
loading_index = 0
is_calculating = False

# Function to update the loading animation
def loading_animation():
    global loading_index
    if is_calculating:
        output_b_value.config(text=loading_texts[loading_index])
        output_q_value.config(text=loading_texts[loading_index])
        loading_index = (loading_index + 1) % len(loading_texts)
        root.after(500, loading_animation)

# Dummy functions to simulate different types of calculations
def calculate_for_one_target(thrust, propellant, profile_factor, var1, val1, B_max_min, B_max_max, Q_mgs_min, Q_mgs_max, precision):
    print(precision)
    B_length = int(precision[0])
    Q_length = int(precision[1])
    B_max_values = [B_max_min + ((B_max_max-B_max_min)/B_length)*i for i in range(B_length)]
    Q_mgs_values = [Q_mgs_min + ((Q_mgs_max-Q_mgs_min)/Q_length)*i for i in range(Q_length)]

    run_one_target(B_max_values, Q_mgs_values, profile_factor, 1.34, thrust, propellant, var1, val1)
    return "Consultez les graphiques"

def calculate_for_two_targets(thrust, propellant, profile_factor, var1, val1, var2, val2, B_max_min, B_max_max, Q_mgs_min, Q_mgs_max, precision):
    B_length = int(precision[0])
    Q_length = int(precision[1])
    B_max_values = [B_max_min + ((B_max_max-B_max_min)/B_length)*i for i in range(B_length)]
    Q_mgs_values = [Q_mgs_min + ((Q_mgs_max-Q_mgs_min)/Q_length)*i for i in range(Q_length)]

    B, Q = run_two_targets(B_max_values, Q_mgs_values, profile_factor, 1.34, thrust, propellant, var1, val1, var2, val2)
    return B, Q

# Function to mimic a long-running calculation
def perform_calculation(thrust, propellant, profile_factor, b_max_min, b_max_max, q_min, q_max, precision, var1=None, val1=None, var2=None, val2=None): 
    global is_calculating
    is_calculating = True
    loading_animation()  # Start the loading animation
    
    if var2 and val2:  # If both target variables are filled
        B, Q = calculate_for_two_targets(thrust, propellant, profile_factor, var1, val1, var2, val2, b_max_min, b_max_max, q_min, q_max, precision)
    else:  # Only the first target variable is filled
        B, Q = calculate_for_one_target(thrust, propellant, profile_factor, var1, val1, b_max_min, b_max_max, q_min, q_max, precision)

    # After the calculation is done, display the result
    output_b_value.config(text=B)
    output_q_value.config(text=Q)
    
    is_calculating = False  # Stop the loading animation

# Function to be executed when 'Run' is clicked
def run_command():
    # Get the inputs from the interface
    thrust = thrust_dropdown.get()
    propellant = propellant_dropdown.get()
    profile_factor = float(profile_entry.get())
    var1 = variable_1_dropdown.get()
    val1 = value_1_entry.get()
    var2 = variable_2_dropdown.get()
    val2 = value_2_entry.get()

    # Get the B and Q min/max values from the scales
    b_min = b_max_min_scale.get()
    b_max = b_max_max_scale.get()
    q_min = q_min_scale.get()
    q_max = q_max_scale.get()
   
    precision = [int(precision_B_entry.get()), int(precision_Q_entry.get())]
    print(precision)
    
    # Check if the second target is filled or not
    if var2 and val2:
        # Run the calculation for two targets in a separate thread
        calculation_thread = threading.Thread(target=perform_calculation, args=(thrust, propellant, profile_factor, b_min, b_max, q_min, q_max, precision,  var1, val1, var2, val2))
    else:
        # Run the calculation for one target in a separate thread
        calculation_thread = threading.Thread(target=perform_calculation, args=(thrust, propellant, profile_factor, b_min, b_max, q_min, q_max, precision, var1, val1, precision))
    
    calculation_thread.start()

# Engineering Inputs section
engineering_frame = tk.LabelFrame(root, text="Engineering Inputs", padx=20, pady=20, font=("Arial", 14, "bold"))
engineering_frame.grid(row=0, column=0, padx=20, pady=20, sticky="n")

# Dropdown for thrust input
thrust_label = ttk.Label(engineering_frame, text="Thruster:")
thrust_label.grid(row=0, column=0, sticky="w", pady=5)
thrust_options = ["PPSX00", "PPS1350", "PPS5000", "PPSX000_Hi", "PPS20k"]
thrust_dropdown = ttk.Combobox(engineering_frame, values=thrust_options)
thrust_dropdown.grid(row=0, column=1, pady=5)

# Dropdown for propellant input
propellant_label = ttk.Label(engineering_frame, text="Propellant:")
propellant_label.grid(row=1, column=0, sticky="w", pady=5)
propellant_options = ["xenon", "krypton"]
propellant_dropdown = ttk.Combobox(engineering_frame, values=propellant_options)
propellant_dropdown.grid(row=1, column=1, pady=5)

# Profile Factor text entry
profile_label = ttk.Label(engineering_frame, text="Profile Factor:")
profile_label.grid(row=2, column=0, sticky="w", pady=5)
profile_entry = ttk.Entry(engineering_frame)
profile_entry.grid(row=2, column=1, pady=5)

# Target section
target_frame = tk.LabelFrame(root, text="Target", padx=20, pady=20, font=("Arial", 14, "bold"))
target_frame.grid(row=0, column=1, padx=20, pady=20, sticky="n")

# First variable and value
variable_1_label = ttk.Label(target_frame, text="1) Variable:")
variable_1_label.grid(row=0, column=0, sticky="w", pady=5)
variable_1_options = ["ISP", "Thrust"]
variable_1_dropdown = ttk.Combobox(target_frame, values=variable_1_options)
variable_1_dropdown.grid(row=0, column=1, pady=5)

value_1_label = ttk.Label(target_frame, text="Value: (s)")
value_1_label.grid(row=1, column=0, sticky="w", pady=5)
value_1_entry = ttk.Entry(target_frame)
value_1_entry.grid(row=1, column=1, pady=5)

# Second variable and value
variable_2_label = ttk.Label(target_frame, text="2) Variable:")
variable_2_label.grid(row=2, column=0, sticky="w", pady=5)
variable_2_dropdown = ttk.Combobox(target_frame, values=variable_1_options)
variable_2_dropdown.grid(row=2, column=1, pady=5)

value_2_label = ttk.Label(target_frame, text="Value: (mN)")
value_2_label.grid(row=3, column=0, sticky="w", pady=5)
value_2_entry = ttk.Entry(target_frame)
value_2_entry.grid(row=3, column=1, pady=5)

# Ranges for B and Q section
range_frame = tk.LabelFrame(root, text="Ranges for B_max and Q", padx=20, pady=20, font=("Arial", 14, "bold"))
range_frame.grid(row=1, column=0, padx=20, pady=20, columnspan=2, sticky="n")

# Fonction pour mettre à jour les valeurs affichées
def update_b_min_max_value(*args):
    b_min_value_label.config(text=f"{b_max_min_scale.get():.1f} Gauss")
    b_max_value_label.config(text=f"{b_max_max_scale.get():.1f} Gauss")

def update_q_min_max_value(*args):
    q_min_value_label.config(text=f"{q_min_scale.get():.1f} mg/s")
    q_max_value_label.config(text=f"{q_max_scale.get():.1f} mg/s")


# B_min and B_max scales
b_max_min_label = ttk.Label(range_frame, text="B_min:")
b_max_min_label.grid(row=0, column=0, pady=5)
b_max_min_scale = ttk.Scale(range_frame, from_=100, to=300, orient="horizontal", command=update_b_min_max_value)
b_max_min_scale.grid(row=0, column=1, pady=5)
b_min_value_label = ttk.Label(range_frame, text="100.0 Gauss")  # Valeur initiale
b_min_value_label.grid(row=0, column=2, padx=5)

b_max_max_label = ttk.Label(range_frame, text="B_max:")
b_max_max_label.grid(row=1, column=0, pady=5)
b_max_max_scale = ttk.Scale(range_frame, from_=100, to=300, orient="horizontal", command=update_b_min_max_value)
b_max_max_scale.grid(row=1, column=1, pady=5)
b_max_value_label = ttk.Label(range_frame, text="100.0 Gauss")  # Valeur initiale
b_max_value_label.grid(row=1, column=2, padx=5)

# Q_min and Q_max scales
q_min_label = ttk.Label(range_frame, text="Q_min:")
q_min_label.grid(row=2, column=0, pady=5)
q_min_scale = ttk.Scale(range_frame, from_=1, to=50, orient="horizontal", command=update_q_min_max_value)
q_min_scale.grid(row=2, column=1, pady=5)
q_min_value_label = ttk.Label(range_frame, text="1.0 mg/s")  # Valeur initiale
q_min_value_label.grid(row=2, column=2, padx=5)

q_max_label = ttk.Label(range_frame, text="Q_max:")
q_max_label.grid(row=3, column=0, pady=5)
q_max_scale = ttk.Scale(range_frame, from_=1, to=50, orient="horizontal", command=update_q_min_max_value)
q_max_scale.grid(row=3, column=1, pady=5)
q_max_value_label = ttk.Label(range_frame, text="1.0 mg/s")  # Valeur initiale
q_max_value_label.grid(row=3, column=2, padx=5)

precision_B_label = ttk.Label(range_frame, text="Number of B steps:")
precision_B_label.grid(row=4, column=0, pady=5)
precision_B_entry = ttk.Entry(range_frame)
precision_B_entry.grid(row=4, column=1, pady=5)

precision_Q_label = ttk.Label(range_frame, text="Number of Q steps:")
precision_Q_label.grid(row=5, column=0, pady=5)
precision_Q_entry = ttk.Entry(range_frame)
precision_Q_entry.grid(row=5, column=1, pady=5)

# Outputs section
output_frame = tk.LabelFrame(root, text="Outputs", padx=20, pady=20, font=("Arial", 14, "bold"))
output_frame.grid(row=2, column=0, padx=20, pady=20, columnspan=3, sticky="n")

# Output labels
output_b_label = ttk.Label(output_frame, text="B:")
output_b_label.grid(row=0, column=0, sticky="w", pady=5)
output_b_value = ttk.Label(output_frame, text="---")
output_b_value.grid(row=0, column=1, pady=5)

output_q_label = ttk.Label(output_frame, text="Q:")
output_q_label.grid(row=1, column=0, sticky="w", pady=5)
output_q_value = ttk.Label(output_frame, text="---")
output_q_value.grid(row=1, column=1, pady=5)

# Run button with styling
run_button = ttk.Button(root, text="Run", command=run_command)
run_button.grid(row=3, column=1, pady=20)

# Add padding around all elements to give more space
for widget in root.winfo_children():
    widget.grid_configure(padx=10, pady=10)

# Start the GUI event loop
root.mainloop()
