import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar

# 1D model of a Hall thruster in xenon and krypton
# Trevor Lafleur and Pascal Chabert, June 2024
# Code from Mananas Assez

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

def test_run(B_max, C_mag, Q_mgs, I_bar, thruster, propellant, plotting=False):

    global sigma_scl, M, f_epsilon, E_iz, m_e, alpha, lmbd, beta, gamma

    print("Plotting: "+ str(plotting))
    
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
    alpha = K_iz_star * L_ch * Q_m / (M * A_ch * v_g * v_star)
    lmbd = 2 * h_R * L_ch / Delta_r #(26)
    gamma = delta_anom * M * v_g * omega_max * A_ch / (Km_star * Q_m) #(25)
    beta = m_e * v_g * omega_max**2 * A_ch * L_ch / (v_star * Km_star * Q_m) #(24)

    n_g0 = Q_mgs * 1e-6 / (M * v_g * A_ch)
    I_d = I_bar * e * Q_m / M #(27)
    Gamma_d = I_d / (e * A_ch) #(15) in the text
    Gamma_m = Q_m / (M * A_ch) #(15) in the text

    # integration starts here
    Gamma_0_bar = 0.01
    Gamma_0 = Gamma_0_bar * Gamma_d
    G_0_bar = 0.001  # initial velocity is then v_star/10
    G_0 = (v_star * Gamma_d) * G_0_bar

    Km = 2.6e-13

    def F_Te(Te, z_bar, G, Gamma):
        global f_ce, f_m, f_iz, sigma, M_bar
        
        Te_bar = Te/E_iz #(28)
        K_iz = K_iz_star * ((3/2 * Te_bar)**0.25) * np.exp(-4/(3 * Te_bar))
        f_ce = np.exp(-C_mag * (z_bar - 1)**2) #(45)
        f_m = Km/Km_star #(31)
        f_m = 1
        f_iz = K_iz/K_iz_star #(30)
        #sigma_see = 0.207 * Te**0.549 #(12)
        sigma = np.minimum(sigma_scl, Te_bar*E_iz/25) #(14) POURQUOI MET ÇA POUR SIGMA SEE ??
        M_bar = M/m_e #(33)

        Gamma_bar = Gamma/Gamma_d #(17) in the text below
        G_bar = G/(v_star * Gamma_d) #(22) in the text just above

        #(32) the unknown is Te
        numerateur = beta * f_ce**2 * G_bar**2 * (1 - Gamma_bar)**2
        denominateur = Gamma_bar**4 * (f_m*(1 - I_bar * Gamma_bar) + gamma * f_ce)
        term1 = alpha * f_iz * f_epsilon * (1 - I_bar * Gamma_bar)
        term2 = lmbd * Te_bar**(3/2)*(2/(1 - sigma) + np.log((1 - sigma)*np.sqrt(M_bar/(2*np.pi))))
        
        return numerateur/denominateur - term1 - term2

    #plots the function F_Te to find the root
    # Te = np.linspace(0.5, 25, 2500)    
    # z = 0
    # F = np.zeros(2500)
    # for i in range(2500):
    #     F[i] = F_Te(Te[i], z, G_0, Gamma_0)
    # plt.plot(Te/E_iz, F)
    # plt.xlabel('Te_Bar')
    # plt.ylabel('F_Te')
    # plt.show()

    Te_bar_global = None

    def RHS_Anna_HET(z_bar, y): #z corresponds to the distance along the channel, y is the vector of the unknowns
        global sigma_scl, C_mag, M, f_epsilon, E_iz, m_e, alpha, lmbd, beta, gamma, I_bar, Te_bar_global
        
        dy = np.zeros(2)
        Gamma, G = y #we calculate Te thanks to an initialisation value of G and Gamma

        # Secant method to find the root of F_Te
        # Initial guess
        Te_0 = 9
        Te_1 = 11
        
        F_Te_0 = F_Te(Te_0, z_bar, G, Gamma)
        F_Te_1 = F_Te(Te_1, z_bar, G, Gamma)
            
        while np.abs(F_Te_0) > 1e-5:
            Te_2 = Te_1 - F_Te_1 * (Te_1 - Te_0) / (F_Te_1 - F_Te_0)        
            Te_0, Te_1 = Te_1, Te_2
            F_Te_0, F_Te_1 = F_Te_1, F_Te(Te_2, z_bar, G, Gamma)
        
        Te_0_bar = Te_0/E_iz #(28)
        Te_bar_global = Te_0_bar
        
        Gamma_bar = Gamma/Gamma_d #(17) in the text below
        G_bar = G/(v_star * Gamma_d) #(22) in the text just above
        
        K_iz = K_iz_star * ((3/2 * Te_0_bar)**0.25) * np.exp(-4/(3 * Te_0_bar))    
        f_iz = K_iz/K_iz_star
        f_ce = np.exp(-C_mag * (z_bar - 1)**2) #(45)
            
        dy[0] = (alpha * f_iz * Gamma_bar**2 * (1 - I_bar * Gamma_bar))/(G_bar) 
        #- (lmbd * Gamma_bar**2 * np.sqrt(Te_0_bar))/G_bar #(21) dGamma_bar/dz_bar
        dy[1] = (beta * f_ce**2 * (1 - Gamma_bar))/(f_m * (1 - I_bar * Gamma_bar) + gamma * f_ce)
        #- lmbd * Gamma_bar * np.sqrt(Te_0_bar) #(22) dG_bar/dz_bar
        
        return dy
            
    # print('Trevor Te = ', RHS_Trevor_HET(0, [chi_0, g_0])[1])
    # print('Anna Te = ', RHS_Anna_HET(0, [chi_0, g_0])[1])

    # print(solve_ivp(RHS_Trevor_HET, [0, 1], [chi_0, g_0], method='RK45', rtol=1e-5, atol=[1e-5, 1e-5]))
    # print(solve_ivp(RHS_Anna_HET, [0, 1], [chi_0, g_0], method='RK45', rtol=1e-5, atol=[1e-5, 1e-5]))

    ##### potentiellement à reprendre !!
    sol = solve_ivp(RHS_Anna_HET, [0, 1/L_ch], [Gamma_0, G_0], method='RK45', rtol=1e-5, atol=[1e-5, 1e-5]) #solution of the system of ODEs with initial conditions Gamma_0 and G_0
    z_bar, Y = sol.t, sol.y #z is the vector of the distances along the channel, Y is the matrix of the unknowns
    z = z_bar/L_ch #(21) just above

    print('z:', z)
    print('z_bar:', z_bar)
    print('Y:', Y)

    N0 = len(z) #number of points in the solution
    Gamma_bar, G_bar = Y[0, :], Y[1, :] #Gamma and G are the vectors of the unknowns
    u_i = ((v_star * Gamma_d) * G_bar / (Gamma_bar * Gamma_d)) * v_star #ion velocity
    n_i = Gamma_d * Gamma_bar / u_i #ion density
    n_g = (Gamma_m - Gamma_d * Gamma_bar) / v_g #(15) neutral density
    f_ce = np.exp(-C_mag * (z/L_ch - 1)**2) #radial magnetic field profile
    beta_bar = gamma * f_ce #normalized electron mobility
    intB2 = np.trapz(f_ce**2, z) #integral of the square of the magnetic field
    BB = B_max * f_ce #magnetic field
    omega_ce = e * BB / m_e #electron cyclotron frequency
    alpha_bar = beta * f_ce**2 #normalized electron mobility
    nu_eff = Km_star * n_g + delta_anom * omega_ce #effective collision frequency
    mu_eff = e * nu_eff / (m_e * omega_ce**2) #effective mobility
    E_x = Gamma_d * (1 - Gamma_bar) / (n_i * mu_eff) #electric field
    V_d = E_c * np.trapz(alpha_bar * G_bar * (1 - Gamma_bar) / (Gamma_bar**2 * (1 + beta_bar - I_bar * Gamma_bar)), z) #voltage drop
    Power = V_d * I_d #power


    # E_w = np.zeros(N0)
    # for kk in range(N0):
    #     E_w[kk] = Te_bar_global[kk] * (2 + (1 - sigma) * np.log((1 - sigma) * np.sqrt(M / (2 * np.pi * m_e)))) #(11)

    # K_iz = K_iz_star * ((3 * Te_bar_global / (2 * E_iz))**0.25) * np.exp(-4 * E_iz / (3 * Te_bar_global))
    # OMEGA = L_ch * K_iz * Q_m / (M * A_ch * v_g * v_star)
    # S_iz = n_i * n_g * K_iz

    # Engineering output
    Gamma_L = Gamma_bar[N0 - 1]
    thrust_mN = n_i[N0 - 1] * u_i[N0 - 1] * M * A_ch * u_i[N0 - 1] * 1e3  # in mN
    thrust_power = 0.5 * M * n_i[N0 - 1] * u_i[N0 - 1]**3 * A_ch
    I_sp = thrust_mN * 1e-3 / (Q_m * 9.81)  # in s
    mass_utilization = Gamma_L * I_bar
    thrust_to_power_mN_kW = 1e3 * thrust_mN / Power
    total_efficiency = (thrust_mN * 1e-3)**2 / (2 * Q_m * Power)
    elec_efficency = thrust_power / (I_d * V_d)


    print("Thrust: "+ str(thrust_mN)+ " mN")
    print("ISP: "+ str(I_sp)+ " s")
    print("Thrust to power ratio: "+ str(thrust_to_power_mN_kW)+ " mN/kW")
    print("Time: "+ str(ti.time() - start_time) + " s")

    # Plotting
    if plotting:
        plt.figure(1)
        plt.plot(z, n_i, 'k', linewidth=1)
        plt.twinx()
        plt.plot(z, n_g, 'r', linewidth=1)
        plt.xlim([0, 1])
        plt.xlabel('$x$ (m)', fontsize=14)
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
        plt.plot(z, u_i, 'b', linewidth=1)
        plt.xlim([0, 1])
        plt.ylabel('$u_i$ (m/s)', fontsize=14)
        plt.gca().tick_params(labelsize=14)
        #plt.savefig('ion-velocity.pdf')

        plt.figure(4)
        plt.plot(z, E_x, 'k', linewidth=1)
        plt.twinx()
        plt.plot(z, S_iz, 'r', linewidth=1)
        plt.xlim([0, 1])
        plt.xlabel('$x$(m)', fontsize=14)
        plt.ylabel('$E_x$ (V/m)', fontsize=14)
        plt.ylabel('$S_{\\rm iz}$ (m$^3$/s$^{-1}$)', fontsize=14)
        plt.gca().tick_params(labelsize=14)
        #plt.savefig('ExSiz.pdf')

        plt.figure(6)
        plt.plot(z, BB * 1e4, 'b', linewidth=1)
        plt.xlim([0, 1])
        plt.xlabel('$x/L_{\\rm ch}$', fontsize=14)
        plt.ylabel('$B$ (Gauss)', fontsize=14)
        plt.gca().tick_params(labelsize=14)
        #plt.savefig('Bfield.pdf')

        plt.figure(8)
        plt.plot(z, n_g * v_g, 'b', label='Neutrals')
        plt.plot(z, n_i * u_i, 'r', label='Ions')
        plt.xlim([0, 1])
        plt.xlabel('$x$(m)', fontsize=14)
        plt.ylabel('$\\Gamma$ (m$^{-2}$s$^{-1}$)', fontsize=14)
        plt.legend(fontsize=14)
        plt.gca().tick_params(labelsize=14)
        plt.show()
        #plt.savefig('Fluxes.pdf')

    return ([I_sp, thrust_mN])

test_run(B_max, C_mag, Q_mgs, I_bar, thruster, propellant, plotting=True)
