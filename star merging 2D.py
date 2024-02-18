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
a = AU / 10 ** 3.2  # in AU
e = 0  # eccentricity
T = 2 * np.pi * np.sqrt(a ** 3 / (G * M))
tm = 1000  # merger time
fusion_distance = a / 27  # distance below which stars will merge
R = 500 * 3.086 * 10 ** (19)  # distance of observation 500Mpc
α = 2
β = 1


def v(a):
    return 2 * np.pi * a / T


def relative_velocity(v1, v2):
    return np.sqrt((v1[0] - v2[0]) ** 2 + (v1[1] - v2[1]) ** 2)


def distance(x1, y1, x2, y2):
    return np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)


def friction_coefficient(r, t):
    return β / (t * np.exp(α * r / a))  # Avoid division by zero


def equation(r, t):
    x1, y1, vx1, vy1, x2, y2, vx2, vy2 = r

    R1 = np.sqrt(x1 ** 2 + y1 ** 2)
    R2 = np.sqrt(x2 ** 2 + y2 ** 2)

    # Gravitational forces
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


def h(t):
    f_GW = [np.sqrt(G * M / distances[i] ** 3) / (np.pi) for i, dist in enumerate(distances)]
    h0 = [
        4 * G * Mc / (c ** 2 * R) * (G / c ** 3 * np.pi * f_GW[i] * Mc) ** (2 / 3) for i, dist in enumerate(distances)
    ]
    dfdt = [
        96 / 5 * c ** 3 / G * f_GW[i] / Mc * (G / c ** 3 * np.pi * f_GW[i] * Mc) ** (8 / 3)
        for i, dist in enumerate(distances)
    ]
    ϕ = [2 * np.pi * (f_GW[i] * t[i] + 0.5 * dfdt[i] * t[i] ** 2) for i, dist in enumerate(distances)]

    h = [0 if dist < fusion_distance else h0[i] * np.cos(ϕ[i]) for i, dist in enumerate(distances)]
    return h


def update(frame):
    dist = distance(solution[frame, 0], solution[frame, 1], solution[frame, 4], solution[frame, 5])

    # Check distance and merge if below fusion_distance
    if dist < fusion_distance:
        star1.set_data(0, 0)
        star2.set_data(0, 0)
        star_merged.set_data(0, 0)
        star_merged.set_markersize(33)  # Adjust size
        star_merged.set_color("gray")  # Set color to gray

    else:
        star1.set_data(solution[frame, 0], solution[frame, 1])
        star2.set_data(solution[frame, 4], solution[frame, 5])
        star_merged.set_data([], [])  # Clear merged star data

    return star1, star2, star_merged


initial_conditions = [a * (1 - e) * m1 / M, 0, 0, v(a), -a * (1 - e) * m2 / M, 0, 0, -v(a)]

t_values = np.arange(0.01, α * tm / β, 1)  # One orbital period

solution = odeint(equation, initial_conditions, t_values)

# Calcul de la distance entre les deux étoiles pour chaque instant
distances = [distance(solution[i, 0], solution[i, 1], solution[i, 4], solution[i, 5]) for i in range(len(t_values))]

# Create the figure and axes
fig = plt.figure(figsize=(10, 12))
fig.set_facecolor("black")
gs = fig.add_gridspec(2, hspace=0.20)

# Upper plot for animation
ax1 = fig.add_subplot(gs[0])
ax1.set_xlim(-1.1 * a / 2, 1.1 * a / 2)
ax1.set_ylim(-1.1 * a / 2, 1.1 * a / 2)
ax1.set_facecolor("black")

star1, = ax1.plot([], [], "o", markersize=20, color="lightblue")
star2, = ax1.plot([], [], "o", markersize=20, color="white")
star_merged, = ax1.plot([], [], "o", markersize=33, color="gray")

# Lower plot for signal
ax2 = fig.add_subplot(gs[1])
ax2.set_xlim(0, α * tm / β)
ax2.set_ylim(min(h(t_values)), max(h(t_values)))
ax2.set_facecolor("black")

ax2.plot(t_values, h(t_values), color="white")
ax2.set_title("Gravitational wave signal (Newtonian approximation)", color="white")
ax2.set_xlabel("Time", color="white")
ax2.set_ylabel("Strain", color="white")

# Customize ticks and spines for the signal plot
ax2.tick_params(axis="x", colors="white")
ax2.tick_params(axis="y", colors="white")
ax2.spines["bottom"].set_color("white")
ax2.spines["left"].set_color("white")

# Adjust the size of the plots
ax1.set_position([0.1, 0.35, 0.8, 0.63])  # Increase the size of the upper plot
ax2.set_position([0.1, 0.07, 0.8, 0.20])  # Decrease the size of the lower plot

# Animation
ani = FuncAnimation(fig, update, frames=len(t_values), interval=1, blit=True)

plt.show()
