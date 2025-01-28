import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from scipy.constants import elementary_charge, electron_mass, Boltzmann, proton_mass
import matplotlib.pyplot as plt
# import Plots as plt_hall
import csv
from Km_plot import Km_Kiz

# 1D model of a Hall thruster in xenon and krypton

#################### CONSTANTS ####################
#Fundamental constants
ee = elementary_charge
m_e = electron_mass
k_B = Boltzmann
T_g = 500 # noter la temperature du gaz

# Engineering inputs
B_max = 0.0186  # Bmax in tesla
C_mag = 4  # Profile factor
Q_mgs = 2.0  # Mass flow rate in mg/s

##################### THRUSTER #####################
# choose thruster and propellant
thruster = 'PPSX00'
propellant = 'xenon'

# Model input is a normalized current
I_bar = 1.34

if thruster == 'PPSX00':
    L_ch = 0.024  # channel length in meters
    d_ave = 0.056  # Thruster average diameter
    Delta_r = 0.015  # channel gap
elif thruster == 'PPS1350':
    L_ch = 0.029  # channel length in meters
    d_ave = 0.0845  # Thruster average diameter
    Delta_r = 0.0155  # channel gap
elif thruster == 'PPS5000':
    L_ch = 0.033  # channel length in meters
    d_ave = 0.1275  # Thruster average diameter
    Delta_r = 0.0225  # channel gap
elif thruster == 'PPSX000_Hi':
    L_ch = 0.032  # channel length in meters
    d_ave = 0.1365  # Thruster average diameter
    Delta_r = 0.0135  # channel gap
elif thruster == 'PPS20k':
    L_ch = 0.071  # channel length in meters
    d_ave = 0.2700  # Thruster average diameter
    Delta_r = 0.05  # channel gap

if propellant == 'xenon':
    MM = 131 * proton_mass
    v_g = np.sqrt(k_B * T_g / MM)
    E_iz = 12.127
    Km_star = 2.37e-13
    K_iz_star = 5.45e-14
    delta_anom = 3.5e-3  # anomalous transport
elif propellant == 'krypton':
    MM = 83.8 * proton_mass
    v_g = np.sqrt(k_B * T_g / MM)
    E_iz = 14.00  # approx ionization rate
    Km_star = 1.91e-13
    K_iz_star = 4.62e-14
    delta_anom = 5e-4  # anomalous transport
elif propellant == 'argon':
    MM = 39.95 * proton_mass
    v_g = np.sqrt(k_B * T_g / MM)
    E_iz = 15.76  # approx ionization rate
    Km_star = 2.07e-13
    K_iz_star = 3.57e-14
    delta_anom = 5e-4  ###### TO FIND 
    

# Geometry of the Thruster
R_out = (d_ave + Delta_r) / 2  # external radius in meters
R_in = (d_ave - Delta_r) / 2  # internal radius in meters
A_ch = np.pi * (R_out**2 - R_in**2)  # channel cross sectional area
omega_max = ee * B_max / m_e # maximum electron cyclotron frequency
Q_m = Q_mgs * 1e-6 # mass flow rate in kg/s
f_epsilon = 2 # see Eq. (34) and Table I
h_R = 0.5 # edge-to-center density ratio
sigma_scl = 1 - 8.3*np.sqrt(m_e / MM) # sigma critical, see Eq. (13)

# Similarity parameters
v_star = np.sqrt(ee * E_iz / MM) # just above (21)
alpha = K_iz_star * L_ch * Q_m / (MM * A_ch * v_g * v_star) #(23)
beta = m_e * v_g * omega_max**2 * A_ch * L_ch / (v_star * Km_star * Q_m) #(24)
gamma = delta_anom * MM * v_g * omega_max * A_ch / (Km_star * Q_m) #(25)
lmbd = 2 * h_R * L_ch / Delta_r #(26)
# f_m = 1 # see Table I

I_d = I_bar * ee * Q_m / MM # denormalized current (27)
Gamma_d = I_d / (ee * A_ch) # after (17) in the text
Gamma_m = Q_m / (MM * A_ch) # after (15) in the text

################# FUNCTIONS ##################

def F_Te(Temp, dist, Gb, Gammab):
    '''This function calculates the electron temperature Te_bar at a given distance z_bar along the channel'''
    '''temp = Te_bar , dist = z_bar, Gb = G_bar, Gammab = Gamma_bar'''
    
    f_ce = np.exp(-C_mag * (dist - 1)**2) #(45)

    # verif_interpol(Temp * E_iz, K_iz, K_iz_x, K_iz_y)
    f_iz = K_iz(Temp * E_iz)/K_iz_star #(30)
    sigma_see = 0.207 * (Temp * E_iz)**.549
    sigma = np.minimum(sigma_scl, sigma_see) #(14)
    M_bar = MM / m_e
    
    # We use the interpolation to find the value of f_m depending on T_e
    # verif_interpol(Temp * E_iz, Km, Km_x, Km_y)
    f_m = Km(Temp * E_iz) / Km_star
    
    #(32) the unknown is Te
    numerateur = beta * f_ce**2 * Gb**2 * (1 - Gammab)**2
    denominateur = Gammab**4 * (f_m * (1 - I_bar * Gammab) + gamma * f_ce)
    term1 = alpha * f_iz * f_epsilon * (1 - I_bar * Gammab)
    term2 = lmbd * Temp**(3/2)*(2/(1 - sigma) + np.log((1 - sigma)*np.sqrt(M_bar/(2*np.pi))))

    return numerateur/denominateur - term1 - term2

def find_sign_change(f, a, b, steps=500):
    '''This function finds the interval [a, b] where the function f changes sign'''
    x_values = [a + j * (b - a) / steps for j in range(steps + 1)]
    for k in range(len(x_values) - 1):
        if f(x_values[k]) * f(x_values[k + 1]) < 0:
            return x_values[k], x_values[k + 1]
    raise ValueError("Pas de changement de signe trouvé dans l'intervalle donné")

def secant_method(Funct, dist, Gb, Gammab, x0, x1):
    '''This function finds the root of the function Funct using the secant method'''
    '''Funct = F_Te, dist = z_bar, Gb = G_bar, Gammab= Gamma_bar, x0 = Te_0_bar, x1 = Te_1_bar'''
    
    F0 = Funct(x0, dist, Gb, Gammab)
    F1 = Funct(x1, dist, Gb, Gammab)
        
    while np.abs(F0) > 1e-3:
        x2 = x1 - F1 * (x1 - x0) / (F1 - F0)        
        x0, x1 = x1, x2
        F0, F1 = F1, Funct(x2, dist, Gb, Gammab)
    xf = x1
    return xf  # Return the final approximation after max iterations

#################### ODEs ######################

z_bar_array = []
Te_bar_array = []


def RHS_Anna_HET(dist, y):
    '''This function calculates the right-hand side of the system of ODEs'''
    '''dist = z_bar, y = [Gamma_bar, G_bar]'''
        
    # z corresponds to the distance along the channel, y is the vector of the unknownss
    dy = np.zeros(2)
    Gamma_bar, G_bar = y # we calculate Te thanks to an initialisation value of G and Gamma
    
    ######## Secant method to find the root of F_Te
    [Te_0_bar, Te_1_bar] = find_sign_change(lambda Te_bar: F_Te(Te_bar, dist, G_bar, Gamma_bar), 0.5/E_iz, 50/E_iz) # initial interval for the secant method
    Te_bar = secant_method(F_Te, dist, G_bar, Gamma_bar, Te_0_bar, Te_1_bar)
    z_bar_array.append(dist) # save the values of z_bar
    Te_bar_array.append(Te_bar) # save the values of Te_bar
    
    f_m = Km(Te_bar * E_iz) / Km_star
    f_iz = K_iz(Te_bar * E_iz)/K_iz_star
    f_ce = np.exp(-C_mag * (dist - 1)**2) #(45)
        
    dy[0] = (alpha * f_iz * Gamma_bar**2 * (1 - I_bar * Gamma_bar))/(G_bar) 
    #- (lmbd * Gamma_bar**2 * np.sqrt(Te_bar))/G_bar #(21) dGamma_bar/dz_bar
    dy[1] = (beta * f_ce**2 * (1 - Gamma_bar))/(f_m * (1 - I_bar * Gamma_bar) + gamma * f_ce) 
    #- (lmbd * Gamma_bar * np.sqrt(Te_bar)) #(22) dG_bar/dz_bar
    
    return dy

#################### INTERPOLATION ######################

## Interpolation : Km and K_iz ##

# prepares the interpolation to find Km in the solver => to reduce the complexity of the code
Km_Kiz(propellant)

Km_x = []
Km_y = []

with open('//Users/annamasset/Desktop/X/Cours/3A/EA/code_statinnaire/Fichiers_CSV/Km_Kiz.csv', 'r') as file:
    reader = csv.reader(file)
    next(reader)
    for row in reader:
        Km_x.append(float(row[0]))
        Km_y.append(float(row[1]))

Km = interp1d(Km_x, Km_y, fill_value='extrapolate')

K_iz_x = []
K_iz_y = []

with open('//Users/annamasset/Desktop/X/Cours/3A/EA/code_statinnaire/Fichiers_CSV/Km_Kiz.csv', 'r') as file:
    reader = csv.reader(file)
    next(reader)
    for row in reader:
        K_iz_x.append(float(row[0]))
        K_iz_y.append(float(row[2]))
                    
K_iz = interp1d(K_iz_x, K_iz_y, fill_value='extrapolate')

# def verif_interpol(x, function, x_list, y_list):
#     if x < x_list[0] or x > x_list[-1]:  # Si x est hors de la plage
#         print(f"Valeur {x} hors de la plage d'interpolation. Mise à jour...")
#         # Ajouter x dans la liste des abscisses
#         x_list.append(x)
#         x_list.sort()  # Trier après ajout
#         # Ajouter une valeur par défaut (par exemple, 0 ou extrapoler depuis la fonction existante)
#         if x < x_list[0]:
#             y_list.insert(0, function(x))  # Ajouter au début
#         elif x > x_list[-1]:
#             y_list.append(function(x))  # Ajouter à la fin
#         # Recréer la fonction d'interpolation
#         return interp1d(x_list, y_list, fill_value='extrapolate')
#     return function  # Si x est dans la plage, pas besoin de mise à jour

            

#################### INTEGRATION ######################

# initial conditions
Gamma_0_bar = 0.005 
G_0_bar = 0.1 * Gamma_0_bar * f_epsilon**(1/2)

# solution
print("Solving ODEs")
sol = solve_ivp(RHS_Anna_HET, [0, 1], [Gamma_0_bar, G_0_bar], method='LSODA', rtol=1e-7, atol=[1e-7, 1e-7]) #solution of the system of ODEs with initial conditions Gamma_0 and G_0
z_bar, Y = sol.t, sol.y # z is the vector of the distances along the channel, Y is the matrix of the unknowns
print("ODEs solved")
z_bar_array = np.array(z_bar_array)
Te_array = np.array(Te_bar_array) * E_iz

###################### RETURN ########################

N0 = len(z_bar) # number of points in the solution
Gamma_bar, G_bar = Y[0, :], Y[1, :] # Gamma and G are the vectors of the unknowns

# Recalculation of Te
Te_bar_end = np.zeros(N0)
print("Recalculating Te")
for i in range(N0):
    [Te_0_bar, Te_1_bar] = find_sign_change(lambda Te_bar: F_Te(Te_bar, z_bar[i], G_bar[i], Gamma_bar[i]), 0.5/E_iz, 50/ E_iz)
    Te_bar_end[i] = secant_method(F_Te, z_bar[i], G_bar[i], Gamma_bar[i], Te_0_bar, Te_1_bar)
print("Te recalculated")
Te_end = Te_bar_end * E_iz

f_ce = np.exp(-C_mag * (z_bar - 1)**2)
BB = B_max * f_ce # (29) magnetic field
omega_ce = ee * BB / m_e # electron cyclotron frequency

def u_i(Gammab, Gb):
    '''This function calculates the ion velocity'''    
    u_i = v_star * Gb / Gammab
    return u_i

def n_i(Gammab, Gb):
    '''This function calculates the ion density'''
    n_i = Gammab * Gamma_d / u_i(Gammab, Gb) # (15) ion density
    return n_i

def n_g(Gammab):
    '''This function calculates the neutral density'''
    n_g = (Gamma_m - Gammab * Gamma_d) / v_g # (15) neutral density
    return n_g

def Ez(Gammab, Gb):
    '''This function calculates the electric field'''
    nu_eff = Km_star * n_g(Gammab) + delta_anom * omega_ce # (6) effective collision frequency
    mu_eff = ee * nu_eff / (m_e * omega_ce**2) # (16) effective mobility
    Ez = Gamma_d * (1 - Gammab) / (n_i(Gammab, Gb) * mu_eff) # (17) electric field
    return Ez

def P_d(Gammab, Gb):
    '''This function calculates the normalized discharge voltage'''
    f_m = Km(Te_end) / Km_star
    Phi_d_bar = np.trapz((beta * f_ce**2 * Gb*(1 - Gammab))/(Gammab**2 * (f_m * (1 - I_bar * Gammab) + gamma * f_ce)), z_bar) # (36) normalized discharge voltage
    P_d = Phi_d_bar * E_iz * I_d # (35under) (41) power
    return P_d

def S_iz(Gammab, Gb):
    '''This function calculates the ionization source'''
    S_iz = n_i(Gammab, Gb) * n_g(Gammab) * K_iz(Te_end)
    return S_iz

# Engineering outputs

def mass_utilization(Gamma_bar):
    '''This function calculates the mass utilization
    Gamma_L = Gamma_bar at the end of the channel, I_bar = normalized current'''
    Gamma_L = Gamma_bar[N0 - 1]
    return Gamma_L * I_bar

def thrust_mN_calcul(Gamma_bar, G_bar):
    '''This function calculates the thrust in (mN)'''
    thrust = n_i(Gamma_bar, G_bar)[N0 - 1] * u_i(Gamma_bar, G_bar)[N0 - 1] * MM * A_ch * u_i(Gamma_bar, G_bar)[N0 - 1] * 1e3  # in mN
    return thrust

def thrust_power_calcul(Gamma_bar, G_bar):
    '''This function calculates the thrust power'''
    thrust_power = 0.5 * MM * n_i(Gamma_bar, G_bar)[N0 - 1] * u_i(Gamma_bar, G_bar)[N0 - 1]**3 * A_ch
    return thrust_power

def thrust_to_power_ratio(Gamma_bar, G_bar):
    '''This function calculates the thrust to power ratio mN_kW'''
    thrust_to_power = thrust_mN_calcul(Gamma_bar, G_bar) / P_d(Gamma_bar, G_bar)
    return thrust_to_power

def I_sp(Gamma_bar, G_bar):
    '''This function calculates the specific impulse'''
    I_sp = thrust_mN_calcul(Gamma_bar, G_bar) * 1e-3 / (Q_m * 9.81)  # in s
    return I_sp

def total_efficiency(Gamma_bar, G_bar):
    '''This function calculates the total efficiency'''
    total_efficiency = (thrust_mN_calcul(Gamma_bar, G_bar) * 1e-3)**2 / (2 * Q_m * P_d(Gamma_bar, G_bar))
    return total_efficiency

def power_efficiency(Gamma_bar, G_bar):
    '''This function calculates the power efficiency'''
    power_efficiency = thrust_power_calcul(Gamma_bar, G_bar) / P_d(Gamma_bar, G_bar)
    return power_efficiency

vng = n_g(Gamma_bar) * v_g
vni = n_i(Gamma_bar, G_bar) * u_i(Gamma_bar, G_bar)


#saves the values of the variables in a csv file
data = np.array([z_bar, n_i(Gamma_bar, G_bar), n_g(Gamma_bar), Te_end, u_i(Gamma_bar, G_bar), Ez(Gamma_bar, G_bar), S_iz(Gamma_bar, G_bar), BB, vng, vni])
np.savetxt("values_anna.csv", data.T, delimiter=',')


#################### PLOTS ######################
'''For now we let the plots in this file but they will incorporated in the interface'''

# plt_hall.plot_densities(z_bar, n_i, n_g, Gamma_bar, G_bar, "n_i", "n_g", 'red', 'blue')

# plt.figure()
# plt_hall.plot_electron_temperature(z_bar_array, Te_array, 'black')

# plt.figure()
# plt_hall.plot_ion_velocity(z_bar, u_i, Gamma_bar, G_bar,'red')

# plt_hall.plot_electric_field_and_ionization_source(z_bar, Ez, S_iz, Gamma_bar, G_bar, "Ez"," S_iz",'black', 'green')

# plt.figure()
# plt_hall.plot_magnetic_field(z_bar, BB)

# plt.figure()
# plt_hall.plot_fluxes(z_bar, vng, vni, 'blue', 'red')

# plt.show()