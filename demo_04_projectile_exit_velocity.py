import numpy as np
import matplotlib.pyplot as plt

# http://persweb.wabash.edu/facstaff/madsenm/publications/AJP_80_24_rohrbach_air_cannon.pdf
# https://www.am1.us/wp-content/uploads/Documents/U24146_Potato_Cannon_and_Chronograph_Part_I.pdf
# https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&cad=rja&uact=8&ved=2ahUKEwjT87rUten1AhU3IjQIHRM0BlMQFnoECAcQAQ&url=http%3A%2F%2Fpersweb.wabash.edu%2Ffacstaff%2Fmadsenm%2Fpublications%2FAJP_80_24_rohrbach_air_cannon.pdf&usg=AOvVaw3cb64OCiWA0-ciEPyqPfQR
# http://persweb.wabash.edu/facstaff/madsenm/publications/AJP_80_24_rohrbach_air_cannon.pdf

# Physical constants
G = 9.81  # acceleration due to gravity (m/s^2)

r = 0.05  # projectile radius (m)
m = 0.2  # mass (kg)

Cd = 0.47  # drag coefficient of a sphere
A = np.pi * r ** 2  # area (m^2)
rho_air = 1.205  # Air density (kg.m-3) @ NTP
rho_hydrogen = 0.0899  # Air density (kg.m-3) @ STP


def run_demo():
    # Initial Conditions -----------------------------------------------------------------------------------------------
    m = 19.40 * 1e-3  # mass [kg]
    h = 4.8019 * 1e-2  # length of projectile [m]

    x0 = 3.35 * 1e-2  # initial projectile position before launch [m]
    xbarrel = 91.6 * 1e-2  # length of barrel [m]
    L = xbarrel - x0  # travel distance [m]

    D = 1.913 * 1e-2  # diameter of barrel [m]
    A = np.pi * (D / 2) ** 2  # cross-sectional area of the barrel [m^2]

    # cannon reservoir volume increases as the projectile moves down the barrel: V0 = A * x0
    V0 = 4.196 * 1e-3  # initial air reservoir volume [m^3]
    P0 = 200  # initial pressure of air reservoir [kPa]

    gamma = 7 / 5  # diatomic ideal gasses
    Patm = 101.325  # pressure at atmosphere [kPa]
    Ff = 0  # force due to friction [N] or [kg·m·s-2] (acetal copolymer)

    # VELOCITY AS A FUNCTION OF BARREL LENGTH --------------------------------------------------------------------------
    x = np.linspace(0, L, 1000)  # [m]
    yad_x = exit_velocity_adiabatic(x, P0, V0, A, m, gamma, Patm, Ff)  # [m/s]
    yiso_x = exit_velocity_isothermal(x, P0, V0, A, m, Patm, Ff)  # [m/s]

    # VELOCITY AS A FUNCTION OF RESERVOIR PRESSURE ---------------------------------------------------------------------
    Pmin_ad = get_root_adiabatic(L, V0, A, gamma, Patm, Ff)  # minimum pressure in adiabatic expansion
    Pmin_iso = get_root_isothermal(L, V0, A, Patm, Ff)  # minimum pressure in isothermal expansion

    xfixed = L - x0
    P0 = np.linspace(max(Pmin_ad, Pmin_iso), 800, 1000)  # x grid of pressures

    yad_p = exit_velocity_adiabatic(xfixed, P0, V0, A, m, gamma, Patm, Ff)  # [m/s]
    yiso_p = exit_velocity_isothermal(xfixed, P0, V0, A, m, Patm, Ff)  # [m/s]

    # PLOT RESULTS -----------------------------------------------------------------------------------------------------
    plot(x, [yad_x, yiso_x], P0, [yad_p, yiso_p])


def exit_velocity_adiabatic(x, P0, V0, A, m, g, Patm, Ff):
    # for adiabatic expansion, has no heat transfer and temperature is reduced
    P0, Patm = (P0 * 1e3), (Patm * 1e3)  # convert to base pascal units
    return np.sqrt((2 / m) * (((P0 * V0) / (g - 1)) * (1 - (V0 / (V0 + A * x)) ** (g - 1)) - (Patm * A * x) - Ff * x))


def exit_velocity_isothermal(x, P0, V0, A, m, Patm, Ff):
    # for isothermal expansion, temperature remains constant and pressure is reduced
    # The gas from a pressurized reservoir expands isothermally
    P0, Patm = (P0 * 1e3), (Patm * 1e3)  # convert to base pascal units
    return np.sqrt((2 / m) * (P0 * V0 * np.log(1 + ((A * x) / V0)) - (Patm * A * x) - Ff * x))


def exit_velocity_corrected(x, P0, V0, A, m, Patm, Ff):
    # the corrected model takes into account the flow rate of air through the valve.
    P0, Patm = (P0 * 1e3), (Patm * 1e3)  # convert to base pascal units
    return np.sqrt((2 / m) * (P0 * V0 * np.log(1 + ((A * x) / V0)) - (Patm * A * x) - Ff * x))


def get_root_adiabatic(x, V0, A, g, Patm, Ff):
    # minimum pressure in adiabatic expansion
    return ((1-g) * x * (-Patm * A - Ff)) / (V0 - (((V0 / (V0 + A * x)) ** g) * (A * x + V0)))


def get_root_isothermal(x, V0, A, Patm, Ff):
    # minimum pressure in isothermal expansion
    return (x * (A * Patm + Ff)) / (V0 * np.log((A * x) / V0 + 1))


def plot(x, y, x2, y2):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8))
    fig.suptitle('Compressed Air Cannon', fontsize=18)

    ax1.plot(x * 1e2, np.asarray(y).T)

    ax1.grid()
    ax1.set_title('Exit Velocity as a function of barrel travel length')
    ax1.set_xlabel('barrel length (travelled) (cm)')
    ax1.set_ylabel('v (m/s)')
    ax1.legend(['adiabatic', 'isothermal'])

    ax2.plot(x2, np.asarray(y2).T)

    ax2.grid()
    ax2.set_title('Exit Velocity as a function of the initial pressure')
    ax2.set_xlabel('initial pressure of air reservoir $P_0 (m)$')
    ax2.set_ylabel('v (m/s)')

    plt.tight_layout()
    plt.show()


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    run_demo()
