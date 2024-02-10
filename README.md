# Star-merging-and-gravitational-wave-signal

Here is a code that simulates in 2D (the 3D version is coming soon) the fusion of stars into a black hole with the associated gravitational wave signal (in Newtonian approximation).

I assume in my simulation that the stars are immersed in a "fluid" filling all space and slowing down their movement through friction. For this fluid, I choose the characteristic fusion time $tm$ at which the stars will merge. I also choose the distance at which my stars merge. If this distance is reached, then the stars turn into a black hole (shown in gray).

Let's start with pen and paper to numerically simulate the Newtonian equations for our compact binary.

We use the 2nd Newton's law : $\displaystyle \sum_k \vec{F}_k = m \vec{a}$ which gives us the equation $\displaystyle m \vec{a} = G \frac{Mm}{r^3} \vec{r} + k \vec{v}$

In cartesian coordinates, the system of equation is quite the same than the one from the solar system : 

 - $\displaystyle \frac{d^2x}{dt^2} = -\frac{GM}{r^2} \cos(\theta) - k(r) \frac{dx}{dt}$
 - $\displaystyle \frac{d^2y}{dt^2} = -\frac{GM}{r^2} \sin(\theta) - k(r) \frac{dy}{dt}$

where $r = \sqrt{x^2 + y^2}$, $\cos(\theta) = \frac{x}{r}$ and $\sin(\theta) = \frac{y}{r}$. 
The friction force $\vec{f} = -k(r) \vec{v}$ makes appear the term $k(r) ~ \frac{1}{t_m} e^{-r/a}$ which is the coefficient of friction.
