import matplotlib.pyplot as plt
import numpy as np

G = 9.81  # m/s^2


def run_demo():
    # velocity initial is tangential to the
    Vo = 10  # m/s
    launch_angle = 45  # deg
    launch_height = 2  # m

    run_projectile_demo(launch_height, launch_angle, Vo)


def run_projectile_demo(h, theta, Vo):
    rad = np.radians(theta)
    Vx = Vo * np.cos(rad)

    # print the governing projectile equation
    print(f'Equation: y = {h} + {round(np.tan(rad), 5)}x - {round(G / (Vx ** 2 * 2.0), 5)}x^2')

    # compute final x position
    xfinal = find_final_position(h, theta, Vo)
    print(f'Distance = {round(xfinal, 3)} m')

    # compute max height of projectile
    ymax = find_height_max(h, theta, Vo)
    print(f'Peak Height = {round(ymax, 3)} m')

    # compute time elapsed
    tfinal = time_of_flight(h, theta, Vo)
    print(f'Time of flight = {round(tfinal, 3)} s')

    # plot
    x = np.linspace(0, xfinal, 100)
    y = h + np.tan(rad) * x - (G / (2.0 * Vx ** 2)) * (x ** 2)

    plot_items = {'xdata': x, 'ydata': y,
                  'xlim': (0, xfinal), 'ylim': (),
                  'xlabel': 'distance (m)', 'ylabel': 'height (m)', 'title': 'Projectile Demo'}
    plot(plot_items)


def find_position_at_time(h, theta, Vo, t):
    rad = np.radians(theta)
    Vy = Vo * np.sin(rad)
    return h + Vy*t - ((G / 2.0) * (t ** 2))


def find_final_position(h, theta, Vo):
    rad = np.radians(theta)
    Vx = Vo * np.cos(rad)
    return time_of_flight(h, theta, Vo) * Vx


def time_of_flight(h, theta, Vo):
    rad = np.radians(theta)
    Vy = Vo * np.sin(rad)
    return (Vy + np.sqrt(Vy ** 2 + 2 * G * h)) / G


def find_height_max(h, theta, Vo):
    rad = np.radians(theta)
    Vy = Vo * np.sin(rad)
    return h + (Vy ** 2) / (2 * G)


def plot(plot_items):
    fig, ax1 = plt.subplots()

    ax1.plot(plot_items['xdata'], plot_items['ydata'], linewidth=2.0)
    plt.grid()

    ax1.set(xlim=plot_items['xlim'])

    ax1.set_title(plot_items['title'])
    ax1.set_xlabel(plot_items['xlabel'])
    ax1.set_ylabel(plot_items['ylabel'])

    plt.show()


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    run_demo()
