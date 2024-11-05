import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar

# 1D model of a Hall thruster in xenon and krypton
# Trevor Lafleur and Pascal Chabert, June 2024

#measures the time it takes to run the code
import time as ti
start_time = ti.time()

global sigma_scl, beta_mag, M, gamma_bar, E_iz, m_e, A_bar, B_bar, alpha, beta, I_bar


e = 1.6e-19
m_e = 9e-31
k_B = 1.38e-23
T_g = 500

# choose thruster and propellant
thruster = 'PPS1350'
propellant = 'xenon'

# Engineering inputs
B_max = 300e-4  # Bmax in tesla
beta_mag = 4  # Profile factor
Q_mgs = 5  # Mass flow rate in mg/s

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
    kiz0 = 1.8 * 1e-13
    E_iz = 12.127
    kexc0 = 1.2921e-13
    E_exc = 11.6
    kel = 2.5e-13
    delta_anom = 3.3e-3  # anomalous transport
elif propellant == 'krypton':
    M = 83.8 * 1.67e-27
    v_g = np.sqrt(k_B * T_g / M)
    kiz0 = 1.7 * 1e-13
    E_iz = 14  # approx ionization rate
    kexc0 = 0.6e-13
    E_exc = 12.75  # excitation rate
    kel = 1.8e-13
    delta_anom = 5e-4  # anomalous transport

R_out = (d_ave + Delta_r) / 2  # external radius in meters
R_in = (d_ave - Delta_r) / 2  # internal radius in meters
A_ch = np.pi * (R_out**2 - R_in**2)  # channel cross sectional area
omega_cemax = e * B_max / m_e
Q_m = Q_mgs * 1e-6
gamma_bar = 2
h_R = 0.5
sigma_scl = 1 - np.sqrt(m_e / M)

# Similarity parameters
E_c = E_iz * gamma_bar
v_star = np.sqrt(e * E_c / M)
A_bar = kiz0 * L_ch * Q_m / (M * A_ch * v_g * v_star)
B_bar = 2 * h_R * L_ch / Delta_r
beta = delta_anom * M * v_g * omega_cemax * A_ch / (kel * Q_m)
alpha = m_e * v_g * omega_cemax**2 * A_ch * L_ch / (v_star * kel * Q_m)

n_g0 = Q_mgs * 1e-6 / (M * v_g * A_ch)
I_d = I_bar * e * Q_m / M
Gamma_d = I_d / (e * A_ch)
Gamma_m = Q_m / (M * A_ch)

# integration starts here
chi_0 = 5e-3
g_0 = 0.1 * chi_0  # initial velocity is then v_star/10

def RHS_Trevor_HET(x, y): #x corresponds to the distance along the channel, y is the vector of the unknowns
    global sigma_scl, beta_mag, M, gamma_bar, E_iz, m_e, A_bar, B_bar, alpha, beta, I_bar

    dy = np.zeros(2)

    chi, g = y

    f = np.exp(-beta_mag * (x - 1)**2) #(45) f effectively represents the normalized radial magnetic field profile with Bmax
    beta_bar = beta * f  
    alpha_bar = alpha * f**2 

    Te = np.linspace(0.5, 25, 2500)
    sigma = np.minimum(sigma_scl, (1/25) * Te)
    E_w = 2 * Te + Te * (1 - sigma) * np.log((1 - sigma) * np.sqrt(M / (2 * np.pi * m_e)))
    term1 = A_bar * (3 * Te / (2 * E_iz))**0.25 * np.exp(-4 * E_iz / (3 * Te)) + B_bar * np.sqrt(Te / (gamma_bar * E_iz)) * (E_w / (E_iz * gamma_bar)) / (1 - sigma) / (1 - I_bar * chi)
    term2 = (1 - chi)**2 * g**2 * alpha_bar / (chi**4 * (1 - I_bar * chi) * (1 + beta_bar - I_bar * chi))
    
    i_min = np.argmin(np.abs(term2 - term1))
    T_e = Te[i_min]
    
    OMEGA = A_bar * (3 * T_e / (2 * E_iz))**0.25 * np.exp(-4 * E_iz / (3 * T_e))

    dy[0] = OMEGA * chi**2 * (1 - I_bar * chi) / g
    dy[1] = alpha_bar * (1 - chi) / (1 + beta_bar - I_bar * chi)

    return dy


sol = solve_ivp(RHS_Trevor_HET, [0, 1], [chi_0, g_0], method='RK45', rtol=1e-5, atol=[1e-5, 1e-5])
x = sol.t
Y = sol.y


# plt.figure(1)
# plt.plot(x, Y[0, :], 'b', linewidth=1)
# plt.plot(x, Y[1, :], 'r', linewidth=1)
# plt.xlabel('$x/L_{\\rm ch}$', fontsize=14)
# plt.ylabel('$\\Gamma$, $G$', fontsize=14)
# plt.legend(['$\\Gamma$', '$G$'], fontsize=14)
# plt.gca().tick_params(labelsize=14)
# plt.show()

N0 = len(x)
chi = Y[0, :]
g = Y[1, :]
u_i = (g / chi) * v_star
n_i = chi * Gamma_d / u_i
n_g = (Gamma_m - chi * Gamma_d) / v_g

print("x: "+ str(x))
f = np.exp(-beta_mag * (x - 1)**2)
beta_bar = beta * f
intB2 = np.trapz(f**2, x)
BB = B_max * f
omega_ce = e * BB / m_e
alpha_bar = alpha * f**2
nu_eff = kel * n_g + delta_anom * omega_ce
mu_eff = e * nu_eff / (m_e * omega_ce**2)
E_x = Gamma_d * (1 - chi) / (n_i * mu_eff)
V_d = E_c * np.trapz(alpha_bar * g * (1 - chi) / (chi**2 * (1 + beta_bar - I_bar * chi)), x)
Power = V_d * I_d

#print toutes les valeurs
print("g: "+ str(g))
print("u_i: "+ str(u_i))
print("n_i: "+ str(n_i))
print("n_g: "+ str(n_g))
print("f: "+ str(f))
print("beta_bar: "+ str(beta_bar))
print("intB2: "+ str(intB2))
print("BB: "+ str(BB))
print("omega_ce: "+ str(omega_ce))
print("alpha_bar: "+ str(alpha_bar))
print("nu_eff: "+ str(nu_eff))
print("mu_eff: "+ str(mu_eff))
print("E_x: "+ str(E_x))
print("V_d: "+ str(V_d))
print("Power: "+ str(Power))


T_e = np.zeros(N0)
E_w = np.zeros(N0)
for kk in range(N0):
    T_e[kk] = 0.5
    for ii in range(2500):
        sigma = min(sigma_scl, (1 / 25) * T_e[kk]**1)
        E_w[kk] = 2 * T_e[kk] + T_e[kk] * (1 - sigma) * np.log((1 - sigma) * np.sqrt(M / (2 * np.pi * m_e)))
        term1 = A_bar * (3 * T_e[kk] / (2 * E_iz))**0.25 * np.exp(-4 * E_iz / (3 * T_e[kk])) + B_bar * np.sqrt(T_e[kk] / (gamma_bar * E_iz)) * (E_w[kk] / (E_iz * gamma_bar)) / (1 - sigma) / (1 - I_bar * chi[kk])
        term2 = (1 - chi[kk])**2 * g[kk]**2 * alpha_bar[kk] / (chi[kk]**4 * (1 - I_bar * chi[kk]) * (1 + beta_bar[kk] - I_bar * chi[kk]))
        if term1 < term2:
            T_e[kk] = 0.5 + 0.01 * ii

K_iz = kiz0 * ((3 * T_e / (2 * E_iz))**0.25) * np.exp(-4 * E_iz / (3 * T_e))
OMEGA = L_ch * K_iz * Q_m / (M * A_ch * v_g * v_star)
S_iz = n_i * n_g * K_iz

# Engineering output
chi_L = chi[N0 - 1]
thrust_mN = n_i[N0 - 1] * u_i[N0 - 1] * M * A_ch * u_i[N0 - 1] * 1e3  # in mN
thrust_power = 0.5 * M * n_i[N0 - 1] * u_i[N0 - 1]**3 * A_ch
I_sp = thrust_mN * 1e-3 / (Q_m * 9.81)  # in s
mass_utilization = chi_L * I_bar
thrust_to_power_mN_kW = 1e3 * thrust_mN / Power
total_efficiency = (thrust_mN * 1e-3)**2 / (2 * Q_m * Power)
elec_efficency = thrust_power / (I_d * V_d)

n_imax, i_nimax = np.max(n_i), np.argmax(n_i)
x_nmax = x[i_nimax]
Siz_max, i_Sizmax = np.max(S_iz), np.argmax(S_iz)
x_infl = x[i_Sizmax]
E_max, i_Emax = np.max(E_x), np.argmax(E_x)
x_Emax = x[i_Emax]
Te_max, i_Temax = np.max(T_e), np.argmax(T_e)
x_Temax = x[i_Temax]

print("Thrust: "+ str(thrust_mN)+ " mN")
print("ISP: "+ str(I_sp)+ " s")
print("Thrust to power ratio: "+ str(thrust_to_power_mN_kW)+ " mN/kW")
print("Time: "+ str(ti.time() - start_time) + " s")


# Plotting

plt.figure(1)
plt.plot(x, n_i, 'k', linewidth=1)
plt.twinx()
plt.plot(x, n_g, 'r', linewidth=1)
plt.xlim([0, 1])
plt.xlabel('$x$ (m)', fontsize=14)
plt.ylabel('$n_i$ m$^{-3}$', fontsize=14)
plt.ylabel('$n_g$ m$^{-3}$', fontsize=14)
plt.gca().tick_params(labelsize=14)
#
# #plt.savefig('densities.pdf')

plt.figure(2)
plt.plot(x, T_e, 'b', linewidth=1)
plt.xlim([0, 1])
plt.ylim([0, 25])
plt.xlabel('$x/L_{\\rm ch}$', fontsize=14)
plt.ylabel('$T_e$ (V)', fontsize=14)
plt.gca().tick_params(labelsize=14)
#plt.savefig('Te.pdf')

plt.figure(3)
plt.plot(x, u_i, 'b', linewidth=1)
plt.xlim([0, 1])
plt.ylabel('$u_i$ (m/s)', fontsize=14)
plt.gca().tick_params(labelsize=14)
#plt.savefig('ion-velocity.pdf')

plt.figure(4)
plt.plot(x, E_x, 'k', linewidth=1)
plt.twinx()
plt.plot(x, S_iz, 'r', linewidth=1)
plt.xlim([0, 1])
plt.xlabel('$x$(m)', fontsize=14)
plt.ylabel('$E_x$ (V/m)', fontsize=14)
plt.ylabel('$S_{\\rm iz}$ (m$^3$/s$^{-1}$)', fontsize=14)
plt.gca().tick_params(labelsize=14)
#plt.savefig('ExSiz.pdf')

plt.figure(6)
plt.plot(x, BB * 1e4, 'b', linewidth=1)
plt.xlim([0, 1])
plt.xlabel('$x/L_{\\rm ch}$', fontsize=14)
plt.ylabel('$B$ (Gauss)', fontsize=14)
plt.gca().tick_params(labelsize=14)
#plt.savefig('Bfield.pdf')

plt.figure(8)
plt.plot(x, n_g * v_g, 'b', label='Neutrals')
plt.plot(x, n_i * u_i, 'r', label='Ions')
plt.xlim([0, 1])
plt.xlabel('$x$(m)', fontsize=14)
plt.ylabel('$\\Gamma$ (m$^{-2}$s$^{-1}$)', fontsize=14)
plt.legend(fontsize=14)
plt.gca().tick_params(labelsize=14)
plt.show()
#plt.savefig('Fluxes.pdf')