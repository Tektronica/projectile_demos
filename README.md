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

![trajectory](https://latex.codecogs.com/svg.image?\bg_white&space;y&space;=&space;y_0&space;&plus;&space;tan(\theta_{i})&space;x&space;-&space;\frac{g}{2V_{0}^2cos(\theta_{i})}x&space;^2)

## Projectile given target angle and distance

The standard projectile equation is used to find the launch angle and velocity given the target distance and angle. This
example maximizes the angle of attack (AoA) to the target to minimize bounce or to shoot through a hoop.

![trajectory](https://latex.codecogs.com/svg.image?\bg_white&space;y&space;=&space;y_0&space;&plus;&space;tan(\theta_{i})&space;x&space;-&space;\frac{g}{2V_{0}^2cos(\theta_{i})}x&space;^2)

## Projectile given max height

The standard projectile equation is used to find the launch angle and velocity given the target distance constrained to
a max ceiling. This example limits the height of the trajectory to avoid a ceiling.

![trajectory](https://latex.codecogs.com/svg.image?\bg_white&space;y&space;=&space;y_0&space;&plus;&space;tan(\theta_{i})&space;x&space;-&space;\frac{g}{2V_{0}^2cos(\theta_{i})}x&space;^2)

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
