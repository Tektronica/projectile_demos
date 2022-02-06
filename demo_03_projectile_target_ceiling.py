import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import numpy as np

G = 9.81  # m/s^2


def run_demo():
    launch_height = 3  # m
    target_distance = 10  # m
    target_height = 1  # m
    ceiling_height = 5  # m

    launch_angle = get_launch_angle(launch_height, target_distance, target_height, ceiling_height)
    target_angle = get_target_angle(launch_height, target_distance, ceiling_height, launch_angle)
    Vo = get_launch_velocity(launch_angle, target_angle, target_distance)
    print(f'Launch Angle: {launch_angle}')
    print(f'Target Angle: {target_angle}')
    print(f'Project Velocity: {Vo}')
    run_projectile_target_ceiling_demo(launch_height, target_height, target_distance, launch_angle, Vo)


# def get_target_angle(y0, xfinal, yfinal, ymax):
#     # y = k - a(x-h)^2
#     h = xfinal/2
#     k = ymax
#     a = ((ymax - y0) - (yfinal - y0))
#     return 43-np.rad2deg(np.arctan(2*a*(xfinal - h)))


def get_target_angle(h, xfinal, ymax, theta_launch):
    rad_launch = np.radians(theta_launch)
    return np.rad2deg(np.arctan(np.tan(rad_launch)-((xfinal)/(2*(ymax - h)))*np.tan(rad_launch)**2))


def get_launch_angle(h, xfinal, yfinal, ymax):
    return np.rad2deg(np.arctan((2 * ((ymax-h) + (np.sqrt((ymax-h)*(ymax-yfinal)))) / xfinal)))


def get_launch_velocity(theta_initial, theta_final, xfinal):
    rad_initial = np.radians(theta_initial)
    rad_final = np.radians(theta_final)
    return (1 / np.cos(rad_initial)) * np.sqrt((G * xfinal) / (np.tan(rad_initial) - np.tan(rad_final)))


def run_projectile_target_ceiling_demo(h, hfinal, xfinal, theta, Vo):
    rad = np.radians(theta)
    Vx = Vo * np.cos(rad)

    # print the governing projectile equation
    print(f'Equation: y = {h} + {round(np.tan(rad), 5)}x - {round(G / (Vx ** 2 * 2.0), 5)}x^2')

    # compute final x position
    # xfinal = find_final_position(h, theta, Vo)
    print(f'Distance = {round(xfinal, 3)} m')

    # compute max height of projectile
    ymax = find_height_max(h, theta, Vo)
    print(f'Peak Height = {round(ymax, 3)} m')

    # compute time elapsed
    tfinal = time_of_flight(h, theta, Vo)
    print(f'Time of Flight = {round(tfinal, 3)} s')

    # plot
    x = np.linspace(0, xfinal, 100)
    y = h + np.tan(rad) * x - (G / (2.0 * Vx ** 2)) * (x ** 2)

    plot_items = {'xdata': x, 'ydata': y,
                  'xlim': (-0.2, xfinal+0.2), 'ylim': (),
                  'box_left': (-1, 0), 'box_left_height': h, 'box_right': (xfinal-1, 0), 'box_right_height': hfinal,
                  'box_ceiling_anchor': (-0.2, ymax), 'box_ceiling_width': xfinal+.4,
                  'xlabel': 'distance (m)', 'ylabel': 'height (m)', 'title': 'Projectile Demo'}
    plot(plot_items)


def find_position_at_time(h, theta, Vo, t):
    rad = np.radians(theta)
    Vy = Vo * np.sin(rad)
    return h + Vy * t - ((G / 2.0) * (t ** 2))


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
    xdata = plot_items['xdata']
    ydata = plot_items['ydata']
    ax1.plot(xdata, ydata, linewidth=2.0)

    # rectangles
    left_anchor = plot_items['box_left']
    right_anchor = plot_items['box_right']
    box_left_height = plot_items['box_left_height']
    box_right_height = plot_items['box_right_height']
    ceiling_anchor = plot_items['box_ceiling_anchor']
    ceiling_width = plot_items['box_ceiling_width']

    boxes = [Rectangle(left_anchor, 2, box_left_height),
             Rectangle(right_anchor, 2, box_right_height),
             Rectangle(ceiling_anchor, ceiling_width, 3)]

    # Create patch collection with specified colour/alpha
    pc = PatchCollection(boxes, facecolor='lightgray', hatch='///', edgecolor='dimgray', linewidth=2)
    ax1.add_collection(pc)

    # limits
    ax1.set(xlim=plot_items['xlim'])

    ax1.set_title(plot_items['title'])
    ax1.set_xlabel(plot_items['xlabel'])
    ax1.set_ylabel(plot_items['ylabel'])

    plt.grid()
    plt.show()


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    run_demo()
