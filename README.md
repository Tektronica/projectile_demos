# Projectile Demonstrations:

# Table of Contents
1. [Introduction](#introduction)
2. [The Standard Projectile](#the-standard-projectile)
3. [Projectile given velocity and distance](#projectile-given-velocity-and-distance)
4. [Projectile given target angle and distance](#projectile-given-target-angle-and-distance)
5. [Projectile given max height](#projectile-given-max-height)
6. [Calculating Exit Velocity](#calculating-exit-velocity)
7. [Projectile with air resistance using iterative method](#projectile-with-air-resistance-using-iterative-method)
8. [Projectile with air resistance using ODE scipy method](#projectile-with-air-resistance-using-ODE-scipy-method)
9. [Projectile with air resistance to target](#projectile-with-air-resistance-to-target)

## Introduction

This repository is a collection of demo scripts that analyze various concepts concerning projectiles both with and
without air resistance. To a limited degree, each section discusses the basic mechanics and derivation of each new
concept introduced.

## The Standard Projectile

![trajectory](https://latex.codecogs.com/svg.image?\bg_white&space;y&space;=&space;y_0&space;&plus;&space;tan(\theta_{i})&space;x&space;-&space;\frac{g}{2V_{0}^2cos(\theta_{i})}x&space;^2)

## Projectile given velocity and distance

The standard projectile equation is used to find the launch angle given a fixed velocity and the distance to the target.
This example contends with a fixed muzzle velocity to adjust the launch angle to achieve the appropriate distance. Note,
however the maximum distance achieved by a fixed velocity launch is at 45 degrees.

<!-- \theta_{i}=tan^{-1} \left(\frac{V_{0}^{2} + \sqrt{V_{0}^{4} - g(g x_{final}^2 + 2(y_{final} - y_{0})V_{0}^{2})}}{g x_{final}}\right) -->
![find launch angle with fixed velocity](https://latex.codecogs.com/svg.image?\bg_white&space;\theta_{i}=tan^{-1}&space;\left(\frac{V_{0}^{2}&space;&plus;&space;\sqrt{V_{0}^{4}&space;-&space;g(g&space;x_{final}^2&space;&plus;&space;2(y_{final}&space;-&space;y_{0})V_{0}^{2})}}{g&space;x_{final}}\right))

## Projectile given target angle and distance

The standard projectile equation is used to find the launch angle and velocity given the target distance and angle. This
example maximizes the angle of attack (AoA) to the target to minimize bounce or to shoot through a hoop.

<!-- \theta_{i} = tan^{-1}\left(\frac{-x_{final}tan(\theta_{final}) - 2(y_{0} + y_{final})}{x_{final}}\right) -->
![find launch angle given target angle](https://latex.codecogs.com/svg.image?\bg_white&space;\theta_{i}&space;=&space;tan^{-1}\left(\frac{-x_{final}tan(\theta_{final})&space;-&space;2(y_{0}&space;&plus;&space;y_{final})}{x_{final}}\right))

<!-- V_{0}=\frac{1}{cos(\theta_{i})} \sqrt{\frac{g x_{final}}{tan(\theta_{i}) - tan(\theta_{final})}} -->
![find launch velocity given target angle](https://latex.codecogs.com/svg.image?\bg_white&space;V_{0}=\frac{1}{cos(\theta_{i})}&space;\sqrt{\frac{g&space;x_{final}}{tan(\theta_{i})&space;-&space;tan(\theta_{final})}})

## Projectile given max height

The standard projectile equation is used to find the launch angle and velocity given the target distance constrained to
a max ceiling. This example limits the height of the trajectory to avoid a ceiling.

<!-- \theta_{final}=tan^{-1}\left({tan(\theta_{i})-\frac{x_{final}}{2(y_{max}- y_{0})}tan^2(\theta_{i})}\right) -->
![find launch velocity given target angle](https://latex.codecogs.com/svg.image?\bg_white&space;\theta_{final}=tan^{-1}\left({tan(\theta_{i})-\frac{x_{final}}{2(y_{max}-&space;y_{0})}tan^2(\theta_{i})}\right))

<!-- \theta_{i}=tan^{-1}\left(\frac{2(y_{max}-h) + \sqrt{(y_{max}-y_{0})(y_{max}-y_{final}})}{x_{final}}\right) -->
![find launch velocity given target angle](https://latex.codecogs.com/svg.image?\bg_white&space;\theta_{i}=tan^{-1}\left(\frac{2(y_{max}-h)&space;&plus;&space;\sqrt{(y_{max}-y_{0})(y_{max}-y_{final}})}{x_{final}}\right))

<!-- V_{0}=\frac{1}{cos(\theta_{i})} \sqrt{\frac{g x_{final}}{tan(\theta_{i}) - tan(\theta_{final})}} -->
![find launch velocity given target angle](https://latex.codecogs.com/svg.image?\bg_white&space;V_{0}=\frac{1}{cos(\theta_{i})}&space;\sqrt{\frac{g&space;x_{final}}{tan(\theta_{i})&space;-&space;tan(\theta_{final})}})

## Calculating Exit Velocity

This demo introduces the mass of the projectile, which requires the exit velocity. This example provides an initial pressure of an air cannon. Currently, the model is over-simplified by not accounting for the flow rate of air through the valve connecting the reservoir and barrel. While the gas from a pressurized reservoir expands isothermally, adiabatic model produces a similar result since the temperature drop associated with the adiabatic expansion is so small.

<!-- v_{ab}={\sqrt{\frac{2}{m}\left(\frac{P_{0}V_{0}}{{\gamma} - 1}\left(1 - \left(\frac{V_{0}}{V_{0} + Ax}\right)^{{\gamma} - 1}\right) - xAP_{atm} - xF_{friction}\right)}}  -->
![exit velocity under adiabatic expansion](https://latex.codecogs.com/png.image?\dpi{110}%20\bg_white%20v_{ab}={\sqrt{\frac{2}{m}\left(\frac{P_{0}V_{0}}{{\gamma}%20-%201}\left(1%20-%20\left(\frac{V_{0}}{V_{0}%20+%20Ax}\right)^{{\gamma}%20-%201}\right)%20-%20xAP_{atm}%20-%20xF_{friction}\right)}})

<!-- v_{iso}={\sqrt{\frac{2}{m}\left(P_{0}V_{0}\ln{\left(1+\frac{xA}{V_{0}}\right)} - xAP_{atm} - xF_{friction}\right)}} -->
![exit velocity under isothermal expansion](https://latex.codecogs.com/png.image?\dpi{110}%20\bg_white%20v_{iso}={\sqrt{\frac{2}{m}\left(P_{0}V_{0}\ln{\left(1+\frac{xA}{V_{0}}\right)}%20-%20xAP_{atm}%20-%20xF_{friction}\right)}})

## Projectile with air resistance using iterative method

This demo introduces air resistance affecting the trajectory of the projectile. The path is plotted using an iterative loop to solve the ODE system of equations.

<!--
\frac{\mathrm{d}}{\mathrm{d}x}
\begin{pmatrix}
 x\\
 y\\
 V_{x}\\
 V_{y}
\end{pmatrix}
=
\begin{pmatrix}
\begin{array}{l}
 V_{x}\\
 V_{y}\\
 -{\mu}V_{x}\sqrt{V_{x}^2+V_{y}^2}\\
 -{\mu}V_{y}\sqrt{V_{x}^2+V_{y}^2}-g
\end{array}
\end{pmatrix}
=
\begin{pmatrix}
\begin{array}{l}
 V_{x}\\
 V_{y}\\
 -{\mu}V_{x}V_{0}\\
 -{\mu}V_{y}V_{0}-g
\end{array}
\end{pmatrix}
-->
![find launch velocity given target angle](https://latex.codecogs.com/svg.image?\bg_white&space;\frac{\mathrm{d}&space;}{\mathrm{d}&space;x}\begin{pmatrix}&space;x&space;\\\\&space;y&space;\\\\&space;V_{x}&space;\\\\&space;V_{y}&space;\end{pmatrix}=\begin{pmatrix}\begin{array}&space;{l}V_{x}&space;\\\\&space;V_{y}&space;\\\\-{\mu}V_{x}\sqrt{V_{x}^2&plus;V_{y}^2}&space;\\\\&space;-{\mu}V_{y}\sqrt{V_{x}^2&plus;V_{y}^2}-g&space;\end{array}\end{pmatrix}=\begin{pmatrix}\begin{array}&space;{l}V_{x}&space;\\\\&space;V_{y}&space;\\\\&space;-{\mu}V_{x}V_{0}&space;\\\\&space;-{\mu}V_{y}V_{0}-g&space;\end{array}\end{pmatrix})

<!--
\begin{matrix}{\mu}=\frac{k}{m}
 &\text{and}&k=\frac{1}{2}C_{d}{\rho}_{air}A\\ 
\end{matrix}
-->
![find launch velocity given target angle](https://latex.codecogs.com/svg.image?\bg_white&space;\begin{matrix}{\mu}=\frac{k}{m}&space;&\text{and}&k=\frac{1}{2}C_{d}{\rho}_{air}A\\&space;\end{matrix})


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


## Projectile with air resistance using ODE scipy method

This demo improves on the previous section by using the scipy built-in function to solve for the ODE. Events are used to trigger when the projectile has reached its target.

![trajectory](https://latex.codecogs.com/svg.image?\bg_white&space;y&space;=&space;y_0&space;&plus;&space;tan(\theta_{i})&space;x&space;-&space;\frac{g}{2V_{0}^2cos(\theta_{i})}x&space;^2)

    def run_projectile_resistance_demo(launch_angle, v0, h):
        phi0 = np.radians(launch_angle)
        rho = rho_air

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
    
        # No drag (ND) -------------------------------------------------------------------------------------------------
        X_ND = []
        Y_ND = []
    
        for t in range(steps + 1):
            X_ND.append(x0 + Vx0 * dt * t)
            Y_ND.append(y0 + Vy0 * dt * t - 0.5 * G * (dt * t)**2)
    
        # With drag (WD) -----------------------------------------------------------------------------------------------
        X_WD = []
        Y_WD = []
    
        # Stop the integration when we hit the target.
        hit_target.terminal = True
        # We must be moving downwards (don't stop before we begin moving upwards!)
        hit_target.direction = -1
    
        # scipy.integrate.solve_ivp(func, t_span, y0, args=() ...)
        soln = solve_ivp(deriv, (t0, tf), u0, args=(rho,), dense_output=True, events=(hit_target, max_height))

        # A fine grid of time points from 0 until impact time.
        t = np.linspace(0, soln.t_events[0][0], steps)

        # Retrieve the solution and append
        sol = soln.sol(t)
        X_WD = sol[0]
        Y_WD = sol[2]
    
        # plot the trajectories.
        x = [X_ND, X_WD]
        y = [Y_ND, Y_WD]
    
        plot(x, y)
    
    
    def deriv(t, u, rho):
        x, y, Vx, Vy = u
        k = 0.5 * Cd * rho * A  # convenience constant
    
        Vxy = np.hypot(Vx, Vy)
        ax = -k / m * Vxy * Vx  # acceleration in x direction
        ay = -k / m * Vxy * Vy - G  # acceleration in y direction
        return Vx, Vy, ax, ay
    
    
    def hit_target(t, u, *args):
        # We've hit the target if the z-coordinate is 0.
        return u[1]
    
    
    def max_height(t, u, *args):
        # The maximum height is obtained when the y-velocity is zero.
        return u[3]


## Projectile with air resistance to target

This demo brings all the concepts together in order to solve the launch angle of a trajectory for a projectile to reach its target while factoring air resistance. Since a closed form solution does not exist, the secant method is used to algorithmically resolve the roots of the ODE system of equations.

![trajectory](https://latex.codecogs.com/svg.image?\bg_white&space;y&space;=&space;y_0&space;&plus;&space;tan(\theta_{i})&space;x&space;-&space;\frac{g}{2V_{0}^2cos(\theta_{i})}x&space;^2)
