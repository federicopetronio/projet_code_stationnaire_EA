import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar

# 1D model of a Hall thruster in xenon and krypton
# Trevor Lafleur and Pascal Chabert, June 2024

#measures the time it takes to run the code
import time as ti

start_time = ti.time()

global sigma_scl, C_mag, M, f_epsilon, E_iz, m_e, alpha, lmbd, beta, gamma, I_bar

#################### CONSTANTS ####################
#Fundamental constants
e = 1.6e-19
m_e = 9e-31
k_B = 1.38e-23
T_g = 500

# Engineering inputs
B_max = 300e-4  # Bmax in tesla
C_mag = 4  # Profile factor
Q_mgs = 5  # Mass flow rate in mg/s

##################### THRUSTER #####################
# choose thruster and propellant
thruster = 'PPS1350'
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
    M = 131 * 1.67e-27
    v_g = np.sqrt(k_B * T_g / M)
    K_iz_star = 1.8 * 1e-13
    E_iz = 12.127
    kexc0 = 1.2921e-13
    E_exc = 11.6
    Km_star = 2.5e-13
    delta_anom = 3.3e-3  # anomalous transport
elif propellant == 'krypton':
    M = 83.8 * 1.67e-27
    v_g = np.sqrt(k_B * T_g / M)
    K_iz_star = 1.7 * 1e-13
    E_iz = 14  # approx ionization rate
    kexc0 = 0.6e-13
    E_exc = 12.75  # excitation rate
    Km_star = 1.8e-13
    delta_anom = 5e-4  # anomalous transport

# Geometry of the Thruster
R_out = (d_ave + Delta_r) / 2  # external radius in meters
R_in = (d_ave - Delta_r) / 2  # internal radius in meters
A_ch = np.pi * (R_out**2 - R_in**2)  # channel cross sectional area
omega_max = e * B_max / m_e # maximum electron cyclotron frequency
Q_m = Q_mgs * 1e-6 # mass flow rate in kg/s
f_epsilon = 2 # see Eq. (34) and Table I
h_R = 0.5 # edge-to-center density ratio
sigma_scl = 1 - 8.3*np.sqrt(m_e / M) # sigma critical, see Eq. (13)


# Similarity parameters
E_c = E_iz * f_epsilon # (34)
v_star = np.sqrt(e * E_iz / M) # just above (21)
alpha = K_iz_star * L_ch * Q_m / (M * A_ch * v_g * v_star) #(23)
beta = m_e * v_g * omega_max**2 * A_ch * L_ch / (v_star * Km_star * Q_m) #(24)
gamma = delta_anom * M * v_g * omega_max * A_ch / (Km_star * Q_m) #(25)
lmbd = 2 * h_R * L_ch / Delta_r #(26)


n_g0 = alpha * v_star / (K_iz_star * L_ch) # Nominal propellant gas density (55)
I_d = I_bar * e * Q_m / M # denormalized current (27)
Gamma_d = I_d / (e * A_ch) # after (17) in the text
Gamma_m = Q_m / (M * A_ch) # after (15) in the text

# integration starts here
Gamma_0_bar = 0.005 
#Gamma_0 = Gamma_0_bar * Gamma_d
#G_0_bar = Gamma_0_bar * np.sqrt((k_B * T_g)/(e * E_iz))  
G_0_bar = 0.1 * Gamma_0_bar
#G_0 = (v_star * Gamma_d) * G_0_bar

#Km = 2.6e-13
#f_m = Km/Km_star #(31)
f_m = 1 # Pour l'instant on prend f_m = 1, see Table I

def F_Te(Te_bar, z_bar, G_bar, Gamma_bar):
    '''This function calculates the electron temperature Te_bar at a given distance z_bar along the channel'''
    # global f_ce, sigma, M_bar
    
    K_iz = K_iz_star * ((3/2 * Te_bar)**0.25) * np.exp(-4/(3 * Te_bar)) # see Eq. (30) paper BM
    # z_bar = z/L_ch #(21) just above
    f_ce = np.exp(-C_mag * (z_bar - 1)**2) #(45)

    f_iz = K_iz/K_iz_star #(30)
    #sigma_see = 0.207 * Te**0.549 #(12)
    # sigma_see = 1/25 * Te_bar * E_iz
    sigma_see = 0.207 * (Te_bar * E_iz)**.549
    sigma = np.minimum(sigma_scl, sigma_see) #(14)
    M_bar = M / m_e
    #(32) the unknown is Te
    numerateur = beta * f_ce**2 * G_bar**2 * (1 - Gamma_bar)**2
    denominateur = Gamma_bar**4 * (f_m*(1 - I_bar * Gamma_bar) + gamma * f_ce)
    term1 = alpha * f_iz * f_epsilon * (1 - I_bar * Gamma_bar)
    term2 = lmbd * Te_bar**(3/2)*(2/(1 - sigma) + np.log((1 - sigma)*np.sqrt(M_bar/(2*np.pi))))
    
    return numerateur/denominateur - term1 - term2


def F_Te_terms(Te_bar, z_bar, G_bar, Gamma_bar):
    '''This function calculates the electron temperature Te_bar at a given distance z_bar along the channel'''
    # global f_ce, sigma, M_bar

    K_iz = K_iz_star * ((3 / 2 * Te_bar) ** 0.25) * np.exp(-4 / (3 * Te_bar))  # see Eq. (30) paper BM
    # z_bar = z/L_ch #(21) just above
    f_ce = np.exp(-C_mag * (z_bar - 1) ** 2)  # (45)

    f_iz = K_iz / K_iz_star  # (30)
    # sigma_see = 0.207 * Te**0.549 #(12)
    # sigma_see = 1/25 * Te_bar * E_iz
    sigma_see = 0.207 * (Te_bar * E_iz) ** .549
    sigma = np.minimum(sigma_scl, sigma_see)  # (14)
    M_bar = M / m_e
    # (32) the unknown is Te
    numerateur = beta * f_ce ** 2 * G_bar ** 2 * (1 - Gamma_bar) ** 2
    denominateur = Gamma_bar ** 4 * (f_m * (1 - I_bar * Gamma_bar) + gamma * f_ce)
    term1 = alpha * f_iz * f_epsilon * (1 - I_bar * Gamma_bar)
    term2 = lmbd * Te_bar ** (3 / 2) * (2 / (1 - sigma) + np.log((1 - sigma) * np.sqrt(M_bar / (2 * np.pi))))

    return numerateur / denominateur, term1, term2


def find_sign_change(f, a, b, steps=100):
    '''This function finds the interval [a, b] where the function f changes sign'''
    x_values = [a + j * (b - a) / steps for j in range(steps + 1)]
    for k in range(len(x_values) - 1):
        if f(x_values[k]) * f(x_values[k + 1]) < 0:
            return x_values[k], x_values[k + 1]
    raise ValueError("Pas de changement de signe trouvé dans l'intervalle donné")

def secant_method(F_Te, z_bar, G_bar, Gamma_bar, Te_0_bar, Te_1_bar):
    '''This function finds the root of the function F_Te using the secant method'''
    F_Te_0 = F_Te(Te_0_bar, z_bar, G_bar, Gamma_bar)
    F_Te_1 = F_Te(Te_1_bar, z_bar, G_bar, Gamma_bar)
        
    while np.abs(F_Te_0) > 1e-3:
        Te_2_bar = Te_1_bar - F_Te_1 * (Te_1_bar - Te_0_bar) / (F_Te_1 - F_Te_0)        
        Te_0_bar, Te_1_bar = Te_1_bar, Te_2_bar
        F_Te_0, F_Te_1 = F_Te_1, F_Te(Te_2_bar, z_bar, G_bar, Gamma_bar)
    Te_bar = Te_1_bar
        
    return Te_bar  # Return the final approximation after max iterations

z_bar_array = []
T_bar_array = []

def RHS_Anna_HET(z_bar, y):
    '''This function calculates the right-hand side of the system of ODEs'''
    #z corresponds to the distance along the channel, y is the vector of the unknowns

    dy = np.zeros(2)
    Gamma_bar, G_bar = y # we calculate Te thanks to an initialisation value of G and Gamma

    # z = z_bar * L_ch #(21) just above
    
    ######## Secant method to find the root of F_Te
    # Initial guess
    [Te_0_bar, Te_1_bar] = find_sign_change(lambda Te_bar: F_Te(Te_bar, z_bar, G_bar, Gamma_bar), 0.5/E_iz, 50/E_iz)
    Te_bar = secant_method(F_Te, z_bar, G_bar, Gamma_bar, Te_0_bar, Te_1_bar)
    z_bar_array.append(z_bar)
    T_bar_array.append(Te_bar)
    
    
    K_iz = K_iz_star * ((3/2 * Te_bar)**0.25) * np.exp(-4/(3 * Te_bar))    
    f_iz = K_iz/K_iz_star
    f_ce = np.exp(-C_mag * (z_bar - 1)**2) #(45)
        
    dy[0] = (alpha * f_iz * Gamma_bar**2 * (1 - I_bar * Gamma_bar))/(G_bar)
    # - (lmbd * Gamma_bar**2 * np.sqrt(Te_bar))/G_bar #(21) dGamma_bar/dz_bar
    dy[1] = (beta * f_ce**2 * (1 - Gamma_bar))/(f_m * (1 - I_bar * Gamma_bar) + gamma * f_ce)
             # - lmbd * Gamma_bar * np.sqrt(Te_bar)) #(22) dG_bar/dz_bar
    
    return dy

#z_eval = np.linspace(0, 1, 1000)
# sol = solve_ivp(RHS_Anna_HET, [0, 1], [Gamma_0_bar, G_0_bar], method='RK45', rtol=1e-7, atol=[1e-7, 1e-7]) #solution of the system of ODEs with initial conditions Gamma_0 and G_0
sol = solve_ivp(RHS_Anna_HET, [0, 1], [Gamma_0_bar, G_0_bar], method='LSODA', rtol=1e-7, atol=[1e-7, 1e-7]) #solution of the system of ODEs with initial conditions Gamma_0 and G_0
z_bar, Y = sol.t, sol.y #z is the vector of the distances along the channel, Y is the matrix of the unknowns

z_bar_array = np.array(z_bar_array)
T_bar_array = np.array(T_bar_array)*E_iz




#plot the solution
# plt.figure(1)
# plt.plot(z_bar, Y[0, :], 'b', linewidth=1)
# plt.plot(z_bar, Y[1, :], 'r', linewidth=1)
# plt.xlabel('$x/L_{\\rm ch}$', fontsize=14)
# plt.ylabel('$\\Gamma$, $G$', fontsize=14)
# plt.legend(['$\\Gamma$', '$G$'], fontsize=14)
# plt.gca().tick_params(labelsize=14)
# plt.show()

N0 = len(z_bar) #number of points in the solution
Gamma_bar, G_bar = Y[0, :], Y[1, :] #Gamma and G are the vectors of the unknowns
u_i = v_star * G_bar/Gamma_bar #ion velocity ######################################(also look at 19 under)
n_i = Gamma_bar * Gamma_d / u_i #(15) in text under, ion density
n_g = (Gamma_m - Gamma_bar * Gamma_d) / v_g #(15) neutral density


# z = z_bar/L_ch #(21) just above
print ("z_bar: "+ str(z_bar))

f_ce = np.exp(-C_mag * (z_bar - 1)**2) #radial magnetic field profile
print("f_ce: "+ str(f_ce))

intB2 = np.trapz(f_ce**2, z_bar) #integral of the square of the magnetic field
BB = B_max * f_ce # (29) magnetic field
omega_ce = e * BB / m_e #electron cyclotron frequency

alpha_bar = beta * f_ce**2 #normalized electron mobility
 
nu_eff = Km_star * n_g + delta_anom * omega_ce # (6) effective collision frequency
mu_eff = e * nu_eff / (m_e * omega_ce**2) # (16) effective mobility
Ez = Gamma_d * (1 - Gamma_bar) / (n_i * mu_eff) # (17) electric field

Phi_d_bar = np.trapz((beta * f_ce**2 * G_bar*(1-Gamma_bar))/(Gamma_bar**2*(f_m*(1-I_bar*Gamma_bar) + gamma*f_ce)), z_bar) # (36) normalized discharge voltage
P_d = Phi_d_bar * E_iz * I_d # (35under) (41) power

# Recalculation of Te
Te_bar_end = np.zeros(N0)
for i in range(N0):
    print(i, 'of', N0)
    [Te_0_bar, Te_1_bar] = find_sign_change(lambda Te_bar: F_Te(Te_bar, z_bar[i], G_bar[i], Gamma_bar[i]), 0.5/E_iz, 50/ E_iz)
    Te_bar_end[i] = secant_method(F_Te, z_bar[i], G_bar[i], Gamma_bar[i], Te_0_bar, Te_1_bar)
Te_end = Te_bar_end * E_iz

K_iz = K_iz_star * ((3 * Te_end / (2 * E_iz))**0.25) * np.exp(-4 * E_iz / (3 * Te_end))
S_iz = n_i * n_g * K_iz

# Engineering output
Gamma_L = Gamma_bar[N0 - 1]
thrust_mN = n_i[N0 - 1] * u_i[N0 - 1] * M * A_ch * u_i[N0 - 1] * 1e3  # in mN
thrust_power = 0.5 * M * n_i[N0 - 1] * u_i[N0 - 1]**3 * A_ch
I_sp = thrust_mN * 1e-3 / (Q_m * 9.81)  # in s
mass_utilization = Gamma_L * I_bar
thrust_to_power_mN_kW = 1e3 * thrust_mN / P_d
total_efficiency = (thrust_mN * 1e-3)**2 / (2 * Q_m * P_d)
elec_efficency = thrust_power / (I_d * Phi_d_bar * E_iz)


print("Thrust: "+ str(thrust_mN)+ " mN")
print("ISP: "+ str(I_sp)+ " s")
print("Thrust to power ratio: "+ str(thrust_to_power_mN_kW)+ " mN/kW")
print("Time: "+ str(ti.time() - start_time) + " s")
#Te_end = np.zeros(N0)
z = z_bar * L_ch
np.savetxt("values_anna.csv", np.column_stack((z, z_bar, n_i, n_g, Te_end, u_i, Ez, S_iz, BB, n_g * v_g, n_i * u_i)), delimiter=',')
np.savetxt("values_anna_Te.csv", np.column_stack((z_bar_array, T_bar_array)), delimiter=',')

# Plotting

#plots the function F_Te to find the root
# Te = np.linspace(0 * E_iz, 3 * E_iz, 2500)  
# x = z_bar[N0-1]
# F = np.zeros(2500)
# for i in range(2500):
#     F[i] = F_Te(Te[i], x, G_bar[N0-1], Gamma_bar[N0-1])

# plt.plot(Te/E_iz, F)
# plt.xlabel('Te_Bar')
# plt.ylabel('F_Te')
# plt.show()

# sigma_see = 1/25 * Te_end
# sigma_scl_array = np.array([1 - 8.3*np.sqrt(m_e / M) for i in range(N0)])
# sigma = np.minimum(sigma_scl, sigma_see) #(14) 
# plt.plot(Te_end, sigma_see, label='sigma_see')
# plt.plot(Te_end, sigma_scl_array, label='sigma_scl')
# plt.plot(Te_end,sigma, label='sigma')
# plt.xlabel('Te_end')
# plt.ylabel('sigma')
# plt.legend()
# plt.show()


def plot_densities(z_bar, n_i, n_g, label_i, label_g, color_i, color_g,ax1=None, ax2=None):
    if ax1 is None or ax2 is None:  # Si les axes ne sont pas fournis, créez-en
        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()
    
    # Tracé des densités ioniques sur l'axe principal
    line1, = ax1.plot(z_bar, n_i, color=color_i, linewidth=1, label=label_i)
    ax1.set_xlabel('$z_{bar}$', fontsize=14)
    ax1.set_ylabel('$n_i$ (m$^{-3}$)', fontsize=14, color=color_i)
    ax1.tick_params(axis='y', labelcolor=color_i)
    
    # Tracé des densités neutres sur l'axe secondaire
    line2, = ax2.plot(z_bar, n_g, color=color_g, linewidth=1, label=label_g)
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
    #plot T_bar_array to see the evolution of Te
    #plt.plot(z_bar_array, T_bar_array, marker='o', color='b', ls ='')
    plt.xlim([0, 1])
    # plt.ylim([0, 25])
    plt.xlabel('$x/L_{\\rm ch}$', fontsize=14)
    plt.ylabel('$T_e$ (V)', fontsize=14)
    plt.gca().tick_params(labelsize=14)
    plt.title('Electron temperature')

def plot_ion_velocity(z_bar, u_i, color):
    plt.plot(z_bar, u_i, linewidth=1, color = color)
    plt.xlim([0, 1])
    plt.ylabel('$u_i$ (m/s)', fontsize=14)
    plt.gca().tick_params(labelsize=14)
    plt.title('Ion velocity')


def plot_electric_field_and_ionization_source(z_bar, Ez, S_iz, label_ez, label_siz, color_ez, color_siz, ax1=None, ax2=None):
    if ax1 is None or ax2 is None:  # Si les axes ne sont pas fournis, créez-en
        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()
    
    # Tracé des densités ioniques sur l'axe principal
    line1, = ax1.plot(z_bar, Ez, color=color_ez, linewidth=1, label=label_ez)
    ax1.set_xlabel('$z_{bar}$', fontsize=14)
    ax1.set_ylabel('$E_x$ (V/m)', fontsize=14, color=color_ez)
    ax1.tick_params(axis='y', labelcolor=color_ez)
    
    # Tracé des densités neutres sur l'axe secondaire
    line2, = ax2.plot(z_bar, S_iz, color=color_siz, linewidth=1, label=label_siz)
    ax2.set_ylabel('$S_{\\rm iz}$ (m$^3$/s$^{-1}$)', fontsize=14, color=color_siz)
    ax2.tick_params(axis='y', labelcolor=color_siz)
    
    # Titre commun
    ax1.set_title('Electric field and ionization source', fontsize=16)
    
    ax1.legend(handles=[line1], loc='upper left', fontsize=12)  # Légende pour ax1
    ax2.legend(handles=[line2], loc='upper right', fontsize=12)  # Légende pour ax2
    
    return ax1, ax2  # Renvoie les axes pour des ajustements supplémentaires

def plot_magnetic_field(z_bar, BB):
    plt.plot(z_bar, BB * 1e4, 'b', linewidth=1)
    plt.xlim([0, 1])
    plt.xlabel('$z/L_{\\rm ch}$', fontsize=14)
    plt.ylabel('$B$ (Gauss)', fontsize=14)
    plt.gca().tick_params(labelsize=14)
    plt.title('Magnetic field')
    

vng = n_g * v_g
vni = n_i * u_i

def plot_fluxes(z_bar, vng, vni, color_g, color_i):
    plt.plot(z_bar, vng, label='Neutrals', color = color_g)
    plt.plot(z_bar, vni, label='Ions', color = color_i)
    plt.xlim([0, 1])
    plt.xlabel('$z_bar$', fontsize=14)
    plt.ylabel('$\\Gamma$ (m$^{-2}$s$^{-1}$)', fontsize=14)
    plt.legend(fontsize=14)
    plt.title('Mass fluxes')
    plt.gca().tick_params(labelsize=14)

# plt.figure()
# plot_densities(z_bar, n_i, n_g, "n_i", "n_g", 'red', 'blue')
# # plt.show()
#
plt.figure()
plot_electron_temperature(z_bar, Te_end, 'black')
plt.axhline(y=(sigma_scl / .207) ** (1 / .549), color='g', linestyle='--', label='$T_{scl}$')
plt.legend()
# # plt.show()
#
# plt.figure()
# plot_ion_velocity(z_bar, u_i, 'red')
# plt.show()
#
# plt.figure()
# plot_electric_field_and_ionization_source(z_bar, Ez, S_iz, "Ez"," S_iz",'black', 'green')
# # plt.show()
#
# plt.figure()
# plot_magnetic_field(z_bar, BB)
# # plt.show()
#
# plt.figure()
# plot_fluxes(z_bar, vng, vni, 'blue', 'red')
plt.show()
