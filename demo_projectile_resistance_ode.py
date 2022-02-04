import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# https://scipython.com/book2/chapter-8-scipy/examples/a-projectile-with-air-resistance/
# https://en.wikipedia.org/wiki/Projectile_motion
# http://orca.phys.uvic.ca/~tatum/classmechs/class7.pdf
# https://ipython-books.github.io/123-simulating-an-ordinary-differential-equation-with-scipy/
# https://en.wikipedia.org/wiki/Ordinary_differential_equation
# https://github.com/rossant/awesome-math/#ordinary-differential-equations
# http://djm.cc/library/Differential_Equations_Phillips_edited.pdf
# https://gamedev.stackexchange.com/a/149612

G = 9.81  # acceleration due to gravity (m/s^2)

r = 0.05  # projectile radius (m)
m = 0.2  # mass (kg)

Cd = 0.47  # drag coefficient
A = np.pi * r ** 2  # area (m^2)
rho_air = 1.28  # Air density (kg.m-3)

k = 0.5 * Cd * rho_air * A  # convenience constant.


def run_demo():
    # Initial speed and launch angle (from the horizontal).
    launch_angle = 65  # deg
    launch_velocity = 50  # m/s
    launch_height = 1.  # m

    run_projectile_resistance_demo(launch_angle, launch_velocity, launch_height)


def run_projectile_resistance_demo(launch_angle, v0, h):
    phi0 = np.radians(launch_angle)

    # Initial conditions
    x0 = 0.
    y0 = h
    Vx0 = v0 * np.cos(phi0)
    Vy0 = v0 * np.sin(phi0)

    u0 = x0, y0, Vx0, Vy0

    # Interval of integration by ODE method up to tf
    t0 = 0
    tf = (Vy0 + np.sqrt(Vy0**2 + 2*G*h)) / G  # time of flight
    steps = 1000
    dt = tf / steps

    # No drag ----------------------------------------------------------------------------------------------------------
    X_ND = list()
    Y_ND = list()

    for t in range(steps + 1):
        X_ND.append(x0 + Vx0 * dt * t)
        Y_ND.append(y0 + Vy0 * dt * t - 0.5 * G * (dt * t)**2)

    # With drag --------------------------------------------------------------------------------------------------------
    # Stop the integration when we hit the target.
    hit_target.terminal = True
    # We must be moving downwards (don't stop before we begin moving upwards!)
    hit_target.direction = -1

    soln = solve_ivp(deriv, (t0, tf), u0, dense_output=True, events=(hit_target, max_height))

    # A fine grid of time points from 0 until impact time.
    t = np.linspace(0, soln.t_events[0][0], 100)

    # Retrieve the solution for the time grid
    sol = soln.sol(t)
    X_WD, Y_WD = sol[0], sol[1]

    # print(soln)
    print(f'Time to target = {soln.t_events[0][0]} s')
    print(f'Time to highest point = {soln.t_events[1][0]} s')
    print(f'Range to target, xmax = {X_WD[-1]} m')
    print(f'Maximum height, zmax = {max(Y_WD)} m')

    # plot the trajectory.
    x = [X_ND, X_WD]
    y = [Y_ND, Y_WD]

    plot(x, y)


def deriv(t, u):
    x, y, Vx, Vy = u

    Vxy = np.hypot(Vx, Vy)
    ax = -k / m * Vxy * Vx  # acceleration in x direction
    ay = -k / m * Vxy * Vy - G  # acceleration in y direction
    return Vx, Vy, ax, ay


def hit_target(t, u):
    # We've hit the target if the z-coordinate is 0.
    return u[1]


def max_height(t, u):
    # The maximum height is obtained when the y-velocity is zero.
    return u[3]


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

    plt.grid()
    plt.show()


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    run_demo()
