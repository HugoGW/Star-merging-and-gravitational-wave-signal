# Star-merging-and-gravitational-wave-signal

Here is a code that simulates in 2D (the 3D version is coming soon) the fusion of stars into a black hole with the associated gravitational wave signal (in Newtonian approximation).

I assume in my simulation that the stars are immersed in a "fluid" filling all space and slowing down their movement through friction. For this fluid, I choose the characteristic fusion time tmtm at which the stars will merge. I also choose the distance at which my stars merge. If this distance is reached, then the stars turn into a black hole (shown in gray).

Let's start with pen and paper to numerically simulate the Newtonian equations for our compact binary.

We use the 2nd Newton's law : $\displaystyle \sum_k \vec{F}_k = m \vec{a}$ which gives us the equation $\displaystyle m_1 \vec{a} = G \frac{m_1 m_2}{r^3} \vec{r} - k \vec{v}$

In Cartesian coordinates, the system of equations for star 1 is quite the same as the one for the solar system:

 - $\displaystyle \frac{d^2x}{dt^2} = -\frac{Gm_2}{r^2} \cos(\theta) - k(r) \frac{dx}{dt}$
 - $\displaystyle \frac{d^2y}{dt^2} = -\frac{Gm_2}{r^2} \sin(\theta) - k(r) \frac{dy}{dt}$

where $r = \sqrt{x^2 + y^2}$, $\cos(\theta) = \frac{x}{r}$ and $\sin(\theta) = \frac{y}{r}$. 
The friction force $\vec{f} = -k(r) \vec{v}$  introduces the term $k(r) \sim \frac{1}{t_m} e^{-r/a}$ which is the coefficient of friction.
We do not take into account the mass in the coefficient of friction. We assume that  $\displaystyle \frac{k(r)}{m_i} \rightarrow k(r)$,  otherwise the merger will take too long.

We numerically solve this differential equation with odeint from the scipy.integrate library according to our initial conditions $r_0 = \pm a(1-e) \frac{m_{1/2}}{m_1+m_2}$ and $v_0 = \pm \sqrt{\frac{GM}{a}}$. 

    initial_conditions = [a * (1 - e)*m1/M, 0, 0, v(a), -a * (1 - e)*m2/M, 0, 0, -v(a)]

It gives us 4 different arrays, 2 for the positions $(x,y)$ and 2 for the velocities $(v_x, v_y)$. We simulate and plot the position $r[i] = \sqrt{x^2[i] + y^2[i]}$ depending on the time. I took a light blue star and a white one. 
As the distance decreases because of the friction, the velocity increases until the two stars are close enough to merge. I arbitrarily chose a distance equal to $a/27$. When this distance is reached by our two stars, they vanish and a gray black hole appears for a few seconds.

    if dist < fusion_distance:
        star1.set_data(0, 0)
        star2.set_data(0, 0)
        star_merged.set_data(0, 0)
        star_merged.set_markersize(33)
        star_merged.set_color('gray')  # Set color to gray

    else:
        star1.set_data(solution[frame, 0], solution[frame, 1])
        star2.set_data(solution[frame, 4], solution[frame, 5])
        star_merged.set_data([], [])  # Clear merged star data


Then, for the gravitational signal, we use the same newtonian approximation that we used before. 
For a compact binary with masses $m_1$ and $m_2$ in a circular orbit with gravitational wave frequency $\displaystyle f_{GW} = 2f = \frac{1}{\pi}\sqrt{\frac{GM}{a^3}}$
then :

 - the Chirp Mass : $\displaystyle \mathcal{M}_c = \frac{(m_1m_2)^{3/5}}{(m_1 + m_2)^{1/5}}$
 - the Scaling Amplitude : $\displaystyle h_0 = 4 \frac{G}{c^2} \frac{\mathcal{M}_c}{R} \big(\frac{G\pi \mathcal{M}_c}{c^3} f \big) ^{2/3}$
 - Chirp $\displaystyle \dot{f} = \frac{96c^3 f}{5G\mathcal{M}_c} \big(\frac{G\pi \mathcal{M}_c}{c^3} f \big) ^{8/3}$
 - the gravitational phase : $\displaystyle \phi(t) = 2\pi (ft + \frac{1}{2} \dot{f}t^2)$


The signal is given by the equation $\displaystyle h(t) = h_0 \cos(\phi(t))$ where :

 - $\displaystyle h_0 = 4 \frac{G}{c^2} \frac{\mathcal{M}_c}{R} \big(\frac{G\pi \mathcal{M}_c}{c^3} f \big) ^{2/3}$
 - $\displaystyle \phi = 2\pi (ft + \frac{1}{2} \dot{f}t^2)$

The waveform is given by : 

$\displaystyle h(t) = h_0 \cos(\phi(t)) = h_0 \cos(2\pi ft + \pi  \dot{f}t^2)$

    def h(t):
    f_GW = [np.sqrt(G*M/distances[i]**3)/(np.pi) for i, dist in enumerate(distances)]
    h0 = [4*G*Mc/(c**2*R) * (G/c**3 * np.pi * f_GW[i]* Mc)**(2/3) for i, dist in enumerate(distances)]
    dfdt = [96/5 * c**3/G * f_GW[i]/Mc * (G/c**3 * np.pi * f_GW[i] * Mc)**(8/3) for i, dist in enumerate(distances)]
    ϕ = [2*np.pi*(f_GW[i]*t[i] + 0.5 * dfdt[i] * t[i]**2) for i, dist in enumerate(distances)]
    
    h = [0 if dist < fusion_distance else h0[i]*np.cos(ϕ[i]) for i, dist in enumerate(distances)]
    return h

Finally, we animate the movement of the two stars as they approach each other and plot the expected gravitational signal using the data collected from the simulation. We take the distance $r[i]$ between the two stars and implement it into the signal $h[i]$. Thus, we end up with a gravitational signal reflecting the dynamics of our star fusion.

<img width="486" alt="image" src="https://github.com/HugoGW/Star-merging-and-gravitational-wave-signal/assets/140922475/5b090e9f-bd35-4a13-833f-987991ee6bab">

https://github.com/HugoGW/Star-merging-and-gravitational-wave-signal/assets/140922475/e2458a05-a4ff-4916-b206-112dd749745d








