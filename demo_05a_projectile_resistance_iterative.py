import numpy as np
import matplotlib.pyplot as plt

# https://physics.stackexchange.com/a/134863

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

    run_projectile_resistance_demo(launch_angle, launch_velocity, launch_height)


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

    for loop, rho in enumerate([rho_hydrogen, rho_air]):
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

    x_list = list()
    y_list = list()
    Vx_list = list()
    Vy_list = list()

    x_list.append(x0)
    y_list.append(y0)
    Vx_list.append(Vx0)
    Vy_list.append(Vy0)

    # Interval of integration by ODE method up to tf -------------------------------------------------------------------
    stop = 0  # stop condition flag to end for loop
    tof = 0  # time of flight
    ttm = 0  # time to max
    last_smallest_Vy = 1e05

    # Initial conditions -------------------------------------------------------------------------------------------
    x, y, Vx, Vy = x0, y0, Vx0, Vy0

    for t in range(1, steps + 1):
        if stop != 1:
            u = x, y, Vx, Vy

            # log data -------------------------------------------------------------------------------------------------
            Vx, Vy, ax, ay = deriv(t, u, rho)

            # increment state by dt ------------------------------------------------------------------------------------
            # position
            x += Vx * dt
            y += Vy * dt

            # velocity
            Vx += ax * dt
            Vy += ay * dt

            t += t * dt

            # log data -------------------------------------------------------------------------------------------------
            x_list.append(x)
            y_list.append(y)

            Vx_list.append(Vx)
            Vy_list.append(Vy)

            # log event - reached highest point why Vy=0 (note: discrete dt step misses 0.0)
            if np.abs(Vy) < last_smallest_Vy:
                last_smallest_Vy = Vy
                ttm = t * dt

            # stop event - hit target
            if y <= 0.0:
                tof = t * dt
                stop = 1

    return x_list, y_list, Vx_list, Vy_list, tof, ttm


def deriv(t, u, rho):
    x, y, Vx, Vy = u
    k = 0.5 * Cd * rho * A  # convenience constant

    Vxy = np.hypot(Vx, Vy)
    ax = -k / m * Vxy * Vx  # acceleration in x direction
    ay = -k / m * Vxy * Vy - G  # acceleration in y direction
    return Vx, Vy, ax, ay


def plot(x, y):
    fig, ax1 = plt.subplots(figsize=(10, 5))
    if isinstance(x[0], list):
        for i in range(len(x)):
            plt.plot(x[i], y[i], label='id %s' % i)
    else:
        plt.plot(x, y)

    ax1.set_title('Projectile Demo with air resistance')
    ax1.set_xlabel('distance (m)')
    ax1.set_ylabel('height (m)')
    ax1.legend(['$\\rho=0.0$ (vacuum)', f'$\\rho={rho_hydrogen}$ (hydrogen)', f'$\\rho={rho_air}$ (air)'])

    plt.grid()
    plt.show()


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    run_demo()
