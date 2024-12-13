import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Constants
G = 6.67430e-11  # Gravitational constant
c = 299792458
m = 1.989e30  # solar mass
m1 = 1.4 * m
m2 = 1.4 * m
M = m1 + m2
Mc = (m1 * m2) ** (0.6) / M ** 0.2
µ = m1 * m2 / M
AU = 150 * 10 ** 9
a = AU / 10 ** 4.5  # in AU
e = 0  # eccentricity
T = 2 * np.pi * np.sqrt(a ** 3 / (G * M))
tm = 50  # merger time
fusion_distance = a / 27  # distance below which stars will merge
R = 500 * 3.086 * 10 ** (19)  # distance of observation 500Mpc
α = 4.1
β = 4
γ = 0.5

# Functions for the system
def distance(x1, y1, x2, y2):
    return np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)

def friction_coefficient(r, t):
    return β / ((t + γ * tm) ** (3 / 5) * np.exp(α * r / a))  # Avoid division by zero

def equation(r, t):
    x1, y1, vx1, vy1, x2, y2, vx2, vy2 = r

    # Gravitational forces
    R1 = np.sqrt(x1 ** 2 + y1 ** 2)
    R2 = np.sqrt(x2 ** 2 + y2 ** 2)

    dvx1dt = -G * m2 * x1 / (R1 ** 3 + 1e-10)  # Avoid division by zero
    dvy1dt = -G * m2 * y1 / (R1 ** 3 + 1e-10)

    dvx2dt = -G * m1 * x2 / (R2 ** 3 + 1e-10)
    dvy2dt = -G * m1 * y2 / (R2 ** 3 + 1e-10)

    # Friction forces
    v_rel = np.array([vx1 - vx2, vy1 - vy2])
    dist = distance(x1, y1, x2, y2)
    k = friction_coefficient(dist, tm)
    dvx1dt -= k * v_rel[0]
    dvy1dt -= k * v_rel[1]
    dvx2dt += k * v_rel[0]
    dvy2dt += k * v_rel[1]

    return [vx1, vy1, dvx1dt, dvy1dt, vx2, vy2, dvx2dt, dvy2dt]

def h(t_values, distances):
    f_GW = np.sqrt(G * M / (distances + 1e-10) ** 3) / (np.pi)
    h0 = 4 * G * Mc / (c ** 2 * R) * (G / c ** 3 * np.pi * f_GW * Mc) ** (2 / 3)
    dfdt = 96 / 5 * c ** 3 / G * f_GW / Mc * (G / c ** 3 * np.pi * f_GW * Mc) ** (8 / 3)
    ϕ = 2 * np.pi * (f_GW * t_values + 0.5 * dfdt * t_values ** 2)

    h_values = np.where(distances < fusion_distance, 0, h0 * np.cos(ϕ))
    return h_values

# Initial conditions
initial_conditions = [a * (1 - e) * m1 / M, 0, 0, np.sqrt(G * M / a), -a * (1 - e) * m2 / M, 0, 0, -np.sqrt(G * M / a)]
t_values = np.arange(0.0, tm / 2, 0.01)
solution = odeint(equation, initial_conditions, t_values)

# Compute distances and strain values
distances = np.array([distance(solution[i, 0], solution[i, 1], solution[i, 4], solution[i, 5]) for i in range(len(t_values))])
h_values = h(t_values, distances)

# Plot and animation
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 9))
fig.set_facecolor("black")

# Upper: strain vs time
ax1.set_xlim(0, tm / 2)
ax1.set_ylim(h_values.min(), h_values.max())
ax1.set_facecolor("black")
ax1.set_title("Gravitational Wave Strain vs Time", color="white")
ax1.set_xlabel("Time (s)", color="white")
ax1.set_ylabel("Strain", color="white")
ax1.tick_params(axis="x", colors="white")
ax1.tick_params(axis="y", colors="white")
ax1.spines["bottom"].set_color("white")
ax1.spines["left"].set_color("white")
line, = ax1.plot([], [], color="white")

# Lower: dynamic orbit plot
ax2.set_xlim(-a, a)
ax2.set_ylim(-a, a)
ax2.set_facecolor("black")
star1, = ax2.plot([], [], "o", markersize=20, color="lightblue")
star2, = ax2.plot([], [], "o", markersize=20, color="white")
star_merged, = ax2.plot([], [], "o", markersize=30, color="gray")

def update(frame):
    # Update strain plot
    line.set_data(t_values[:frame], h_values[:frame])

    # Update orbit plot
    dist = distance(solution[frame, 0], solution[frame, 1], solution[frame, 4], solution[frame, 5])
    if dist < fusion_distance:
        star1.set_data([], [])
        star2.set_data([], [])
        star_merged.set_data(0, 0)
    else:
        star1.set_data(solution[frame, 0], solution[frame, 1])
        star2.set_data(solution[frame, 4], solution[frame, 5])
        star_merged.set_data([], [])

    return line, star1, star2, star_merged

ax2.set_position([0.1, 0.35, 0.8, 0.63])  # Increase the size of the upper plot
ax1.set_position([0.1, 0.07, 0.8, 0.20])  # Decrease the size of the lower plot

ani = FuncAnimation(fig, update, frames=len(t_values), interval=1, blit=True)

plt.show()
