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
omega_max = e * B_max / m_e 
Q_m = Q_mgs * 1e-6
f_epsilon = 2
h_R = 0.5 # edge-to-center density ratio
sigma_scl = 1 - 8.3*np.sqrt(m_e / M) #(13)


# Similarity parameters
E_c = E_iz * f_epsilon #(34)
v_star = np.sqrt(e * E_iz / M) #(21) just above
alpha = K_iz_star * L_ch * Q_m / (M * A_ch * v_g * v_star) #(23)
beta = m_e * v_g * omega_max**2 * A_ch * L_ch / (v_star * Km_star * Q_m) #(24)
gamma = delta_anom * M * v_g * omega_max * A_ch / (Km_star * Q_m) #(25)
lmbd = 2 * h_R * L_ch / Delta_r #(26)


n_g0 = alpha * v_star / (K_iz_star * L_ch) #(55)
I_d = I_bar * e * Q_m / M #(27)
Gamma_d = I_d / (e * A_ch) #(15) in the text
Gamma_m = Q_m / (M * A_ch) #(15) in the text

# integration starts here
Gamma_0_bar = 0.005 
#Gamma_0 = Gamma_0_bar * Gamma_d
G_0_bar = Gamma_0_bar * np.sqrt((k_B * T_g)/(e * E_iz))  
#G_0 = (v_star * Gamma_d) * G_0_bar

####### TESTTTTT
#Km = 2.6e-13
#f_m = Km/Km_star #(31)
f_m = 1 # Pour l'instant on prend f_m = 1

def F_Te(Te, z, G_bar, Gamma_bar):
    global f_ce, sigma, M_bar
    
    Te_bar = Te/E_iz #(28)
    K_iz = K_iz_star * ((3/2 * Te_bar)**0.25) * np.exp(-4/(3 * Te_bar))
    z_bar = z/L_ch #(21) just above
    f_ce = np.exp(-C_mag * (z_bar - 1)**2) #(45)

    f_iz = K_iz/K_iz_star #(30)
    sigma_see = 0.207 * Te**0.549 #(12)
    sigma = np.minimum(sigma_scl, sigma_see) #(14) 
    M_bar = M/m_e #(33)

    #Gamma_bar = Gamma/Gamma_d #(17) in the text below
    #G_bar = G/(v_star * Gamma_d) #(22) in the text just above

    #(32) the unknown is Te
    numerateur = beta * f_ce**2 * G_bar**2 * (1 - Gamma_bar)**2
    denominateur = Gamma_bar**4 * (f_m*(1 - I_bar * Gamma_bar) + gamma * f_ce)
    term1 = alpha * f_iz * f_epsilon * (1 - I_bar * Gamma_bar)
    term2 = lmbd * Te_bar**(3/2)*(2/(1 - sigma) + np.log((1 - sigma)*np.sqrt(M_bar/(2*np.pi))))
    
    return numerateur/denominateur - term1 - term2

#plots the function F_Te to find the root
# Te = np.linspace(0 * E_iz, 3 * E_iz, 2500)    
# z = 1
# F = np.zeros(2500)
# for i in range(2500):
#     F[i] = F_Te(Te[i], z, G_0_bar, Gamma_0_bar)
    
# plt.plot(Te/E_iz, F)
# plt.xlabel('Te_Bar')
# plt.ylabel('F_Te')
# plt.show()

Te_bar_global = None

def RHS_Anna_HET(z_bar, y): #z corresponds to the distance along the channel, y is the vector of the unknowns
    global Te_bar_global, K_iz, f_iz
    
    dy = np.zeros(2)
    Gamma_bar, G_bar = y #we calculate Te thanks to an initialisation value of G and Gamma

    z = z_bar * L_ch #(21) just above
    
    ######## Secant method to find the root of F_Te
    #Initial guess
    Te_0 = 0.5 * E_iz
    Te_1 = Te_0 + 0.1 * E_iz
    
    F_Te_0 = F_Te(Te_0, z, G_bar, Gamma_bar)
    F_Te_1 = F_Te(Te_1, z, G_bar, Gamma_bar)
        
    while np.abs(F_Te_0) > 1e-5:
        Te_2 = Te_1 - F_Te_1 * (Te_1 - Te_0) / (F_Te_1 - F_Te_0)        
        Te_0, Te_1 = Te_1, Te_2
        F_Te_0, F_Te_1 = F_Te_1, F_Te(Te_2, z, G_bar, Gamma_bar)
    
    # Te = np.linspace(0.5, 25, 2500)
    
    # i_min = np.argmin(np.abs(F_Te(Te, z, G_bar, Gamma_bar)))
    # Te_0 = Te[i_min]

    
    Te_0_bar = Te_0/E_iz #(28)
    Te_bar_global = Te_0_bar
    
    #Gamma_bar = Gamma/Gamma_d #(17) in the text below
    #G_bar = G/(v_star * Gamma_d) #(22) in the text just above
    
    K_iz = K_iz_star * ((3/2 * Te_0_bar)**0.25) * np.exp(-4/(3 * Te_0_bar))    
    f_iz = K_iz/K_iz_star
    f_ce = np.exp(-C_mag * (z_bar - 1)**2) #(45)
        
    dy[0] = (alpha * f_iz * Gamma_bar**2 * (1 - I_bar * Gamma_bar))/(G_bar) 
    #- (lmbd * Gamma_bar**2 * np.sqrt(Te_0_bar))/G_bar #(21) dGamma_bar/dz_bar
    dy[1] = (beta * f_ce**2 * (1 - Gamma_bar))/(f_m * (1 - I_bar * Gamma_bar) + gamma * f_ce)
    #- lmbd * Gamma_bar * np.sqrt(Te_0_bar) #(22) dG_bar/dz_bar
    
    return dy
t_eval = np.linspace(0, 1, 20)
sol = solve_ivp(RHS_Anna_HET, [0, 1], [Gamma_0_bar, G_0_bar], method='RK45', rtol=1e-5, atol=[1e-5, 1e-5], t_eval=t_eval) #solution of the system of ODEs with initial conditions Gamma_0 and G_0
z_bar, Y = sol.t, sol.y #z is the vector of the distances along the channel, Y is the matrix of the unknowns

print("z: "+ str(z_bar))
print("Y: "+ str(Y))


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
n_i = Gamma_bar * Gamma_d / u_i #ion density
n_g = (Gamma_m - Gamma_bar * Gamma_d) / v_g #(15) neutral density


z = z_bar/L_ch #(21) just above
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

#V_d = E_c * np.trapz(alpha_bar * G_bar * (1 - Gamma_bar) / (Gamma_bar**2 * (1 + gamma * f_ce - I_bar * Gamma_bar)), z) #voltage drop
Phi_d_bar = np.trapz((beta * f_ce**2 * G_bar*(1-Gamma_bar))/(Gamma_bar**2*(f_m*(1-I_bar*Gamma_bar) + gamma*f_ce)), z_bar) # (36) normalized discharge voltage
P_d = Phi_d_bar * E_iz * I_d # (35under) (41) power

#print toutes les valeurs
# print("Gamma_bar: "+ str(Gamma_bar))
# print("G_bar: "+ str(G_bar))
#print("u_i: "+ str(u_i))
#print("n_i: "+ str(n_i))
#print("n_g: "+ str(n_g))
#print("f_ce: "+ str(f_ce))
# print("intB2: "+ str(intB2))
# print("BB: "+ str(BB))
# print("omega_ce: "+ str(omega_ce))
# print("alpha_bar: "+ str(alpha_bar))
# print("nu_eff: "+ str(nu_eff))
# print("mu_eff: "+ str(mu_eff))
# print("Ez: "+ str(Ez))
# print("Phi_d_bar: "+ str(Phi_d_bar))
# print("P_d: "+ str(P_d))


# E_w = np.zeros(N0)
# for kk in range(N0):
#     E_w[kk] = Te_bar_global[kk] * (2 + (1 - sigma) * np.log((1 - sigma) * np.sqrt(M / (2 * np.pi * m_e)))) #(11)

K_iz = K_iz_star * ((3 * Te_bar_global / (2 * E_iz))**0.25) * np.exp(-4 * E_iz / (3 * Te_bar_global))
#OMEGA = L_ch * K_iz * Q_m / (M * A_ch * v_g * v_star)
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

# n_imax, i_nimax = np.max(n_i), np.argmax(n_i)
# zip_nmax = z[i_nimax]
# Siz_max, i_Sizmax = np.max(S_iz), np.argmax(S_iz)
# z_infl = z[i_Sizmax]
# E_max, i_Emax = np.max(Ez), np.argmax(Ez)
# z_Emax = z[i_Emax]
# Te_max, i_Temax = np.max(Te_bar_global), np.argmax(Te_bar_global)
# z_Temax = z[i_Temax]

print("Thrust: "+ str(thrust_mN)+ " mN")
print("ISP: "+ str(I_sp)+ " s")
print("Thrust to power ratio: "+ str(thrust_to_power_mN_kW)+ " mN/kW")
print("Time: "+ str(ti.time() - start_time) + " s")


# Plotting

plt.figure(1)
plt.plot(z_bar, n_i, 'k', linewidth=1)
plt.twinx()
plt.plot(z_bar, n_g, 'r', linewidth=1)
plt.legend(['$n_i$', '$n_g$'], fontsize=14)
plt.xlim([0, 1])
plt.xlabel('$z_bar$', fontsize=14)
plt.ylabel('$n_i$ m$^{-3}$', fontsize=14)
plt.ylabel('$n_g$ m$^{-3}$', fontsize=14)
plt.gca().tick_params(labelsize=14)
#
# #plt.savefig('densities.pdf')

# plt.figure(2)
# plt.plot(z, Te_bar_global, 'b', linewidth=1)
# plt.xlim([0, 1])
# plt.ylim([0, 25])
# plt.xlabel('$x/L_{\\rm ch}$', fontsize=14)
# plt.ylabel('$T_e$ (V)', fontsize=14)
# plt.gca().tick_params(labelsize=14)
# #plt.savefig('Te.pdf')

plt.figure(3)
plt.plot(z_bar, u_i, 'b', linewidth=1)
plt.xlim([0, 1])
plt.ylabel('$u_i$ (m/s)', fontsize=14)
plt.gca().tick_params(labelsize=14)
#plt.savefig('ion-velocity.pdf')

plt.figure(4)
plt.plot(z_bar, Ez, 'k', linewidth=1)
plt.twinx()
plt.plot(z_bar, S_iz, 'r', linewidth=1)
plt.xlim([0, 1])
plt.xlabel('$z_bar$', fontsize=14)
plt.ylabel('$E_x$ (V/m)', fontsize=14)
plt.ylabel('$S_{\\rm iz}$ (m$^3$/s$^{-1}$)', fontsize=14)
plt.gca().tick_params(labelsize=14)
#plt.savefig('ExSiz.pdf')

plt.figure(6)
plt.plot(z_bar, BB * 1e4, 'b', linewidth=1)
plt.xlim([0, 1])
plt.xlabel('$z/L_{\\rm ch}$', fontsize=14)
plt.ylabel('$B$ (Gauss)', fontsize=14)
plt.gca().tick_params(labelsize=14)
#plt.savefig('Bfield.pdf')

plt.figure(8)
plt.plot(z_bar, n_g * v_g, 'b', label='Neutrals')
plt.plot(z_bar, n_i * u_i, 'r', label='Ions')
plt.xlim([0, 1])
plt.xlabel('$z_bar$', fontsize=14)
plt.ylabel('$\\Gamma$ (m$^{-2}$s$^{-1}$)', fontsize=14)
plt.legend(fontsize=14)
plt.gca().tick_params(labelsize=14)
plt.show()
#plt.savefig('Fluxes.pdf')