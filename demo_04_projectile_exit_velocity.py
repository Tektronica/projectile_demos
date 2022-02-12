import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# http://persweb.wabash.edu/facstaff/madsenm/publications/AJP_80_24_rohrbach_air_cannon.pdf
# https://arxiv.org/pdf/1106.2803.pdf
# https://aircannonplans.com/pdf/air-cannon-velocity.pdf

# https://www.mathworks.com/matlabcentral/answers/307969-ode45-doesn-t-run
# https://physics.stackexchange.com/q/518845

# https://www.am1.us/wp-content/uploads/Documents/U24146_Potato_Cannon_and_Chronograph_Part_I.pdf
# (section 2.2.6) - https://www.emerson.com/documents/automation/manual-flow-measurement-user-manual-ras-en-133856.pdf
# https://www.cpp.edu/~pbsiegel/phy409/compnotes/diffeq.pdf
# http://faculty.olin.edu/bstorey/Notes/DiffEq.pdf

# Physical Conditions --------------------------------------------------------------------------------------------------
G = 9.81  # acceleration due to gravity (m/s^2)
Cd = 0.47  # drag coefficient of a sphere
rho_air = 1.205  # Air density (kg.m-3) @ NTP

Gg = 1  # Specific gravity of air
T = 293  # temperature, [k]
K = 1.380e-23  # boltzmann constant, [J/K]
Z = 1  # compressibility Factor
Cv = 1.93  # flow coefficient
B = 3.11e19  # engineering Flow Constant converts quantity into units of molecules per time, sqrt(k)/(Pa*s)
gamma = 7 / 5  # diatomic ideal gasses


def run_demo():
    # Initial Conditions -----------------------------------------------------------------------------------------------
    m = 19.40 * 1e-3  # mass [kg]
    h = 4.8019 * 1e-2  # length of projectile [m]

    x0 = 3.35 * 1e-2  # initial projectile position before launch [m]
    xbarrel = 91.6 * 1e-2  # length of barrel [m]
    L = xbarrel - x0  # travel distance [m]
    d = xbarrel - L  # initial distance from the valve to the slug [m]
    D = 1.913 * 1e-2  # diameter of barrel [m]
    A = np.pi * (D / 2) ** 2  # cross-sectional area of the barrel [m^2]

    # cannon reservoir volume increases as the projectile moves down the barrel: V0 = A * x0
    V0 = 4.196 * 1e-3  # initial air reservoir volume [m^3]
    Pt0 = 200 * 1e3  # initial pressure of air reservoir [Pa]
    Patm = 101.325 * 1e3  # pressure at atmosphere [kPa]
    rmax = 0.8  # Flow Ratio (For chocked vs unchocked)

    Ff = 0  # force due to friction [N] or [kg·m·s-2] (acetal copolymer)

    # VELOCITY AS A FUNCTION OF BARREL LENGTH --------------------------------------------------------------------------
    x = np.linspace(0, L, 1000)  # [m]
    yad_x = exit_velocity_adiabatic(x, Pt0, V0, A, m, gamma, Patm, Ff)  # [m/s]
    yiso_x = exit_velocity_isothermal(x, Pt0, V0, A, m, Patm, Ff)  # [m/s]

    # VELOCITY AS A FUNCTION OF RESERVOIR PRESSURE ---------------------------------------------------------------------
    Pmin_ad = get_root_adiabatic(L, V0, A, gamma, Patm, Ff)  # minimum pressure in adiabatic expansion
    Pmin_iso = get_root_isothermal(L, V0, A, Patm, Ff)  # minimum pressure in isothermal expansion

    xfixed = L - x0
    Px = np.linspace(max(Pmin_ad, Pmin_iso), 800 * 1e3, 100)  # x grid of pressures

    yad_p = exit_velocity_adiabatic(xfixed, Px, V0, A, m, gamma, Patm, Ff)  # [m/s]
    yiso_p = exit_velocity_isothermal(xfixed, Px, V0, A, m, Patm, Ff)  # [m/s]

    x0 = 0.0
    v0 = 0.0
    xc, yc = exit_velocity_corrected_ode(x0, v0, V0, Pt0, Patm, A, d, m, rmax)  # x0, v0, V0, P0, Patm, A, d
    yc_p = exit_velocity_corrected(L, x0, v0, V0, Px, Patm, A, d, m, rmax)  # x0, v0, V0, P0, Patm, A, d

    # PLOT RESULTS -----------------------------------------------------------------------------------------------------
    plot([x, x, xc], [yad_x, yiso_x, yc], [Px, Px, Px], [yad_p, yiso_p, yc_p])  # [x], [y] , [x2], [y2]


def exit_velocity_adiabatic(x, P0, V0, A, m, g, Patm, Ff):
    # for adiabatic expansion, has no heat transfer and temperature is reduced
    return np.sqrt((2 / m) * (((P0 * V0) / (g - 1)) * (1 - (V0 / (V0 + A * x)) ** (g - 1)) - (Patm * A * x) - Ff * x))


def get_root_adiabatic(x, V0, A, g, Patm, Ff):
    # minimum pressure in adiabatic expansion
    return ((1 - g) * x * (-Patm * A - Ff)) / (V0 - (((V0 / (V0 + A * x)) ** g) * (A * x + V0)))


def exit_velocity_isothermal(x, P0, V0, A, m, Patm, Ff):
    # for isothermal expansion, temperature remains constant and pressure is reduced
    # The gas from a pressurized reservoir expands isothermally
    return np.sqrt((2 / m) * (P0 * V0 * np.log(1 + ((A * x) / V0)) - (Patm * A * x) - Ff * x))


def get_root_isothermal(x, V0, A, Patm, Ff):
    # minimum pressure in isothermal expansion
    return (x * (A * Patm + Ff)) / (V0 * np.log((A * x) / V0 + 1))


def exit_velocity_corrected_ode(x0, v0, V0, P0, Patm, A, d, m, rmax):
    # the corrected model takes into account the flow rate of air through the valve.
    Pt0 = P0
    Pb0 = Patm  # initial barrel pressure
    N0 = Pt0 * V0 / (K * T)  # initial number of molecules in the tank
    Nb0 = Pb0 * A * d / (K * T)  # initial number of molecules in the barrel

    # Interval of integration by ODE method up to tf -------------------------------------------------------------------
    t0 = 0
    tf = 0.04
    steps = int(1e4)

    # Initial conditions -----------------------------------------------------------------------------------------------
    u0 = x0, v0, N0, Nb0
    args = (V0, Patm, A, d, m, rmax)

    # STOP CONDITIONS --------------------------------------------------------------------------------------------------
    # Stop the integration when we hit the target.
    # end_barrel.terminal = True

    # COMPUTE ODE ------------------------------------------------------------------------------------------------------
    # scipy.integrate.solve_ivp(func, t_span, y0, args=() ...)
    # each event occurs at the zeros of a continuous function of time and state
    soln = solve_ivp(deriv, (t0, tf), u0, args=args, dense_output=True)

    # A fine grid of time points from 0 until impact time.
    t = np.linspace(0, tf, steps)

    # Retrieve the solution and append
    sol = soln.sol(t)

    position = sol[0]
    velocity = sol[1]

    return position, velocity


def exit_velocity_corrected(L, x0, v0, V0, Px, Patm, A, d, m, rmax):
    # the corrected model takes into account the flow rate of air through the valve.
    ycp = []

    # Interval of integration by ODE method up to tf -------------------------------------------------------------------
    t0 = 0.0
    tf = 0.04
    steps = 1e5  # step size must be small enough to avoid "NaN" "values due to overstepping
    dt = (tf - t0) / steps

    for i, Pt in enumerate(Px[:]):
        # Initial conditions -------------------------------------------------------------------------------------------
        x = x0
        v = v0
        t = t0

        Pb0 = Patm  # initial barrel pressure
        Nt = Pt * V0 / (K * T)  # initial number of molecules in the tank
        Nb = Pb0 * A * d / (K * T)  # initial number of molecules in the barrel

        # iterative loop -----------------------------------------------------------------------------------------------
        while x <= L:
            u = x, v, Nt, Nb

            # ODE loop -------------------------------------------------------------------------------------------------
            v, a, dNt, dNb = deriv(t, u, V0, Patm, A, d, m, rmax)

            # increment state by dt ------------------------------------------------------------------------------------
            x += v * dt
            v += a * dt

            Nt += dNt * dt
            Nb += dNb * dt

            t += t * dt
        ycp.append(v)

    return ycp


def get_Q(Pt, Pb, rmax):
    # regime selection -------------------------------------------------------------------------------------------------
    # the molecular flow rate Q through the valve is a function of the ratio
    r = (Pt - Pb) / Pt

    if r <= rmax:
        # non-choked regime
        Q = B * Pt * Cv * (1 - (r / (3 * rmax))) * np.sqrt(r / (Gg * T * Z))
    else:
        # choked regime
        Q = (2 / 3) * B * Pt * Cv * np.sqrt(rmax / (Gg * T * Z))

    return Q


def deriv(t, u, *args):
    # STATE VARIABLES --------------------------------------------------------------------------------------------------
    x, v, Nt, Nb = u

    # CONSTANTS --------------------------------------------------------------------------------------------------------
    V0, Patm, A, d, m, rmax = args

    # FLOW RATE --------------------------------------------------------------------------------------------------------
    Pt = (Nt * K * T) / V0  # Tank Pressure
    Pb = (Nb * K * T) / (A * (d + x))  # Barrel Pressure
    Q = get_Q(Pt, Pb, rmax)

    # SYSTEM OF DIFFERENTIAL EQUATIONS ---------------------------------------------------------------------------------
    a = A * (Pb - Patm) / m  # acceleration [m/s^2]
    # k = 0.5 * Cd * rho_air * A
    # a = a + (-k / m * v * Vx)
    dNt = -Q  # Tank Molecules number Differential
    dNb = Q  # Barrel Molecules Number Differential

    return v, a, dNt, dNb


# def end_barrel(t, u, *args):
#     # projectile has reached the end of barrel
#     return u[5]
#
#
# def low_pressure(t, u, *args):
#     # tank pressure is less than barrel pressure
#     return u[0]


def plot(x1, y1, x2, y2):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8))
    fig.suptitle('Compressed Air Cannon', fontsize=18)

    for i in range(len(x1)):
        ax1.plot(x1[i] * 1e2, y1[i])

    ax1.grid()
    ax1.set_title('Exit Velocity as a function of barrel travel length (200 kPa)')
    ax1.set_xlabel('barrel length (travelled) (cm)')
    ax1.set_ylabel('v (m/s)')
    ax1.legend(['Adiabatic Expansion', 'Isothermal Expansion', 'Valve flow rate corrected'])

    for i in range(len(x2)):
        ax2.plot(x2[i] * 1e-3, y2[i])

    ax2.grid()
    ax2.set_title('Exit Velocity as a function of the initial pressure (L = 88 cm)')
    ax2.set_xlabel('Initial pressure of air reservoir (kPa)')
    ax2.set_ylabel('v (m/s)')

    plt.tight_layout()
    plt.show()


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    run_demo()
