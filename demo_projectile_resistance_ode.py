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

G = 9.81  # acceleration due to gravity (m/s^2)

r = 0.05  # projectile radius (m)
m = 0.2  # mass (kg)

c = 0.47  # drag coefficient
A = np.pi * r ** 2  # area (m^2)
rho_air = 1.28  # Air density (kg.m-3)

k = 0.5 * c * rho_air * A  # convenience constant.


def run_demo():
    # Initial speed and launch angle (from the horizontal).
    launch_angle = 65  # deg
    launch_velocity = 50  # m/s
    launch_height = 1  # m

    run_projectile_resistance_demo(launch_angle, launch_velocity, launch_height)


def run_projectile_resistance_demo(launch_angle, v0, h):
    phi0 = np.radians(launch_angle)

    # Initial conditions: x0, v0_x, y0, v0_y.
    u0 = 0, v0 * np.cos(phi0), 0., v0 * np.sin(phi0)
    # Interval of integration - Integrate up to tf unless we hit the target sooner.
    t0, tf = 0, 50

    # Stop the integration when we hit the target.
    hit_target.terminal = True
    # We must be moving downwards (don't stop before we begin moving upwards!)
    hit_target.direction = -1

    soln = solve_ivp(deriv, (t0, tf), u0, dense_output=True, events=(hit_target, max_height))

    # A fine grid of time points from 0 until impact time.
    t = np.linspace(0, soln.t_events[0][0], 100)

    # Retrieve the solution for the time grid
    sol = soln.sol(t)
    x, y = sol[0], sol[2]

    # print(soln)
    print(f'Time to target = {soln.t_events[0][0]} s')
    print(f'Time to highest point = {soln.t_events[1][0]} s')
    print(f'Range to target, xmax = {x[-1]} m')
    print(f'Maximum height, zmax = {max(y)} m')

    # plot the trajectory.
    plot(x, y)


def deriv(t, u):
    x, Vx, y, Vy = u
    Vo = np.hypot(Vx, Vy)
    ax = -k / m * Vo * Vx  # acceleration in x direction
    ay = -k / m * Vo * Vy - G  # acceleration in y direction
    return Vx, ax, Vy, ay


def hit_target(t, u):
    # We've hit the target if the z-coordinate is 0.
    return u[2]


def max_height(t, u):
    # The maximum height is obtained when the y-velocity is zero.
    return u[3]


def plot(x, y):
    fig, ax1 = plt.subplots()
    plt.plot(x, y)

    ax1.set_title('Projectile Demo with air resistance')
    ax1.set_xlabel('distance (m)')
    ax1.set_ylabel('height (m)')

    plt.grid()
    plt.show()


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    run_demo()
