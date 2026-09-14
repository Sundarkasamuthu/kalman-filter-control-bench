# Discrete-Time Luenberger State Observer for Double Integrator System

A MATLAB implementation of a 2nd-order discrete-time Luenberger state observer designed for state estimation and noise filtering using continuous-to-discrete system transformation and dual-method pole placement.

## Project Overview

This repository implements Part I & Part II of a state estimation suite for an electro-mechanical servo system modeled as a double-integrator plant ($A_2, b_2$). 

Key features include:
* **Exact Discretization:** Closed-form continuous-to-discrete transformation ($F_2, g_2$) verified analytically via inverse Laplace transforms $\mathcal{L}^{-1}\{(sI - A)^{-1}\}$ and numerically via `expm()`.
* **Observer Gain Design:** Feedback matrix $h_2$ derived using Ackermann's matrix polynomial expansion and verified against direct algebraic pole placement and MATLAB's `acker()` function.
* **Real-Time Data Estimation:** Reconstruction of angular velocity ($\dot{x}$) from noisy position measurements ($x$) obtained from physical hardware logfiles.

## Observer Dynamics & Mathematical Formulation

The continuous-time plant model:
$$\dot{x}(t) = \begin{bmatrix} 0 & 1 \\ 0 & 0 \end{bmatrix} x(t) + \begin{bmatrix} 0 \\ k_{ga} \end{bmatrix} u(t)$$

Discretized system over sampling time $T$:
$$x[k+1] = F_2 x[k] + g_2 u[k]$$
$$y[k] = c_2^T x[k], \quad c_2 = \begin{bmatrix} 1 \\ 0 \end{bmatrix}$$

Luenberger Observer state update equation:
$$\hat{x}[k+1] = F_2 \hat{x}[k] + g_2 u[k] + h_2 \left( y[k] - c_2^T \hat{x}[k] \right)$$

## Results & Verification

The observer was evaluated against experimental hardware log data (`xPosMeas`, `xPosEstLog`, `xDotEstLog`).

![Observer Performance vs Hardware Log Data](docs/observer_verification.png)

* **Position Tracking ($x$):** The observer output (`xPosSim2Obs`) tracks physical camera position measurements seamlessly with zero phase delay.
* **Velocity Estimation ($\dot{x}$):** The state observer (`xDotSim2Obs`) acts as an optimal estimator, successfully reconstructing velocity while smoothing out high-frequency differentiation noise visible in raw numerical differences (`xDotDifLog`).

## Repository Structure

* `src/main_luenberger_obs.m` - Main MATLAB simulation script containing discretization, gain verification, open-loop tests, and observer loop.
* `data/` - Hardware test drive logfiles.
* `docs/` - Performance plots and analytical verification notes.

## How to Run

1. Clone this repository:
   ```bash
   git clone [https://github.com/YOUR_USERNAME/kalman-filter-control-bench.git](https://github.com/YOUR_USERNAME/kalman-filter-control-bench.git)
