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
    iterations = 5
    u0 = (-0.5, 0)
    boundary = (-2, 2)

    root, x_guess, y_guess = secant_iterative_method(fun2, u0, iterations)
    print(f"Root: {root}")

    x_curve = np.linspace(boundary[0], boundary[1], 100)
    y_curve = fun2(x_curve)

    plot(x_curve, y_curve, x_guess, y_guess)


def secant_iterative_method(fun, u0, iterations, atol=1e-05):
    # https://en.wikipedia.org/wiki/Secant_method

    x0, x1 = u0
    x2 = np.nan

    x_guess = [None] * (iterations + 2)
    y_guess = [None] * (iterations + 2)

    x_guess[0], x_guess[1] = x0, x1
    y_guess[0], y_guess[1] = fun(x0), fun(x1)

    """Return the root calculated using the secant method."""
    for i in range(iterations):
        x2 = x1 - fun(x1) * ((x1 - x0) / float(fun(x1) - fun(x0)))
        x0, x1 = x1, x2

        x_guess[i + 2] = x2
        y_guess[i + 2] = fun(x2)
    print(x_guess)
    return x2, x_guess, y_guess


def fun1(x):
    # (x0, x1) = (10, 30)
    # Root: 24.738633748750722
    return x ** 2 - 612


def fun2(x):
    # (x0, x1) = (-0.33, 1.25)
    # Root: 0.451583
    return np.sin(np.cos(np.exp(x)))


def plot(x, y, x_guess, y_guess):
    fig, (ax1, ax2) = plt.subplots(2, 1)
    fig.suptitle('Secant Method', fontsize=18)

    ax1.plot(x, y)
    ax1.scatter(x_guess, y_guess, marker='^', c='tab:orange')
    ax1.hlines(0, x[0], x[-1], colors='r', linestyles='--')

    for i in range(len(x_guess)):
        ax1.annotate(f'x{i}', (x_guess[i], y_guess[i]))

    ax1.grid()
    ax1.set_title('Iterative method for converging on root')
    ax1.set_xlabel('distance (m)')
    ax1.set_ylabel('height (m)')
    ax1.legend(['curve', 'guess'])

    ax2.plot(x, y)
    ax2.scatter(x_guess, y_guess, marker='^', c='tab:orange')
    ax2.hlines(0, x[0], x[-1], colors='r', linestyles='--')
    ax2.set_xlim(min(x_guess)*(1+0.25), max(x_guess)*(1+0.25))

    for i in range(len(x_guess)):
        ax2.annotate(f'x{i}', (x_guess[i], y_guess[i]))

    ax2.grid()
    ax2.set_title('zoomed')
    ax2.set_xlabel('distance (m)')
    ax2.set_ylabel('height (m)')

    plt.tight_layout()
    plt.show()


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    run_demo()
