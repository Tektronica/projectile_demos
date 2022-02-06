import numpy as np
import matplotlib.pyplot as plt

# https://physics.stackexchange.com/a/127994
# https://en.wikipedia.org/wiki/Secant_method

# Physical constants
G = 9.81  # acceleration due to gravity (m/s^2)

r = 0.05  # projectile radius (m)
m = 0.2  # mass (kg)

Cd = 0.47  # drag coefficient of a sphere
A = np.pi * r ** 2  # area (m^2)
rho_air = 1.205  # Air density (kg.m-3) @ NTP
rho_hydrogen = 0.0899  # Air density (kg.m-3) @ STP


def run_demo():
    # Initial speed and launch angle (from the horizontal).
    launch_angle = 65  # deg
    launch_velocity = 50  # m/s
    launch_height = 1  # m

    # run_projectile_resistance_demo(launch_angle, launch_velocity, launch_height)
    root = secant_iterative_method(fun2, (-0.33, 1.25), 5)
    print(f"Root: {root}")  # Root: 24.738633748750722


def secant_iterative_method(fun, u0, iterations, atol=1e-05):
    # https://en.wikipedia.org/wiki/Secant_method

    x0, x1 = u0
    x2 = np.nan

    """Return the root calculated using the secant method."""
    for i in range(iterations):
        x2 = x1 - fun(x1) * ((x1 - x0) / float(fun(x1) - fun(x0)))
        x0, x1 = x1, x2

    return x2


def fun1(x):
    # (x0, x1) = (10, 30)
    # Root: 24.738633748750722
    return x ** 2 - 612


def fun2(x):
    # (x0, x1) = (-0.33, 1.25)
    # Root: 0.451583
    return np.sin(np.cos(np.exp(x)))


def run_projectile_resistance_demo(launch_angle, v0, h):
    phi0 = np.radians(launch_angle)

    # Initial conditions
    x0 = 0
    y0 = h
    Vx0 = v0 * np.cos(phi0)
    Vy0 = v0 * np.sin(phi0)

    # Interval of integration by iteration up to tf
    steps = 1000
    tf = (Vy0 + np.sqrt(Vy0 ** 2 + 2 * G * h)) / G  # time of flight
    dt = tf / steps

    # No drag ----------------------------------------------------------------------------------------------------------
    X_ND = list()
    Y_ND = list()

    for t in range(steps + 1):
        X_ND.append(x0 + Vx0 * dt * t)
        Y_ND.append(y0 + Vy0 * dt * t - 0.5 * G * (dt * t) ** 2)

    # With drag --------------------------------------------------------------------------------------------------------
    X_WD = list()
    Y_WD = list()

    for loop, rho in enumerate([rho_air, rho_hydrogen]):
        sol = iterative((x0, y0, Vx0, Vy0), dt, steps, rho)

        # Retrieve the solution and append
        X_WD.append(sol[0])
        Y_WD.append(sol[1])
        tof = sol[4]
        ttm = sol[5]

        print(f'Time to target = {tof} s')
        print(f'Time to highest point = {ttm} s')
        print(f'Range to target, xmax = {X_WD[-1]} m')
        print(f'Maximum height, zmax = {max(Y_WD[loop])} m')

    # Plot results
    x = [X_ND] + X_WD
    y = [Y_ND] + Y_WD

    plot(x, y)


def iterative(u0, dt, steps, rho):
    x0, y0, Vx0, Vy0 = u0

    x = list()
    y = list()
    Vx = list()
    Vy = list()

    x.append(x0)
    y.append(y0)
    Vx.append(Vx0)
    Vy.append(Vy0)

    stop = 0  # stop condition flag to end for loop
    tof = 0  # time of flight
    ttm = 0  # time to max
    last_smallest_Vy = 1e05

    k = 0.5 * Cd * rho * A  # convenience constant

    for t in range(1, steps + 1):
        if stop != 1:
            Vxy = np.hypot(Vx[t - 1], Vy[t - 1])

            # First calculate velocity
            Vx.append(Vx[t - 1] * (1.0 - k / m * Vxy * dt))
            Vy.append(Vy[t - 1] + (- G - k / m * Vy[t - 1] * Vxy) * dt)

            # Now calculate position
            x.append(x[t - 1] + Vx[t - 1] * dt)
            y.append(y[t - 1] + Vy[t - 1] * dt)

            # log event - reached highest point why Vy=0 (note: discrete dt step misses 0.0)
            if np.abs(Vy[t]) < last_smallest_Vy:
                last_smallest_Vy = Vy[t]
                ttm = t * dt

            # stop event - hit target
            if y[t] <= 0.0:
                tof = t * dt
                stop = 1

    return x, y, Vx, Vy, tof, ttm


def plot(x, y):
    fig, ax1 = plt.subplots()
    if isinstance(x[0], list):
        for i in range(len(x)):
            plt.plot(x[i], y[i], label='id %s' % i)
    else:
        plt.plot(x, y)

    ax1.set_title('Projectile Demo with air resistance')
    ax1.set_xlabel('distance (m)')
    ax1.set_ylabel('height (m)')
    ax1.legend(['$\\rho=0$', f'$\\rho={rho_air}$ (air)', f'$\\rho={rho_hydrogen}$ (hydrogen)'])

    plt.grid()
    plt.show()


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    run_demo()
