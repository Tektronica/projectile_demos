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

This demo introduces the mass of the projectile, which requires the exit velocity. This example provides an initial pressure of an air cannon. 

![trajectory](https://latex.codecogs.com/svg.image?\bg_white&space;y&space;=&space;y_0&space;&plus;&space;tan(\theta_{i})&space;x&space;-&space;\frac{g}{2V_{0}^2cos(\theta_{i})}x&space;^2)

## Projectile with air resistance using iterative method

This demo introduces air resistance affecting the trajectory of the projectile. The path is plotted using an iterative loop to solve the ODE system of equations.

![trajectory](https://latex.codecogs.com/svg.image?\bg_white&space;y&space;=&space;y_0&space;&plus;&space;tan(\theta_{i})&space;x&space;-&space;\frac{g}{2V_{0}^2cos(\theta_{i})}x&space;^2)

## Projectile with air resistance using ODE scipy method

This demo improves on the previous section by using the scipy built-in function to solve for the ODE. Events are used to trigger when the projectile has reached its target.

![trajectory](https://latex.codecogs.com/svg.image?\bg_white&space;y&space;=&space;y_0&space;&plus;&space;tan(\theta_{i})&space;x&space;-&space;\frac{g}{2V_{0}^2cos(\theta_{i})}x&space;^2)

## Projectile with air resistance to target

This demo brings all the concepts together in order to solve the launch angle of a trajectory for a projectile to reach its target while factoring air resistance. Since a closed form solution does not exist, the secant method is used to algorithmically resolve the roots of the ODE system of equations.

![trajectory](https://latex.codecogs.com/svg.image?\bg_white&space;y&space;=&space;y_0&space;&plus;&space;tan(\theta_{i})&space;x&space;-&space;\frac{g}{2V_{0}^2cos(\theta_{i})}x&space;^2)
