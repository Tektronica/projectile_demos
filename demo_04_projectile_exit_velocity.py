import numpy as np
import matplotlib.pyplot as plt

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
    # Initial speed and launch angle (from the horizontal).
    launch_angle = 65  # deg
    launch_velocity = 50  # m/s
    launch_height = 1  # m

    run_projectile_exit_velocity(launch_angle, launch_velocity, launch_height)


def run_projectile_exit_velocity(launch_angle, v0, h):
    phi0 = np.radians(launch_angle)

    # Initial conditions
    x0 = 0
    y0 = h
    Vx0 = v0 * np.cos(phi0)
    Vy0 = v0 * np.sin(phi0)

    x, y = [1, 2], np.array([[1, 2], [2, 3]])
    plot(x, y.T)


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
