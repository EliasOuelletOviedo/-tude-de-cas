import numpy as np
from graphique_effort_interne import graph
import matplotlib.pyplot as plt
from scipy.integrate import cumtrapz

plt.rcParams['text.usetex'] = True


n_points_zone = 100
length = 0.95

# Calcul de T
end_zone1 = 0.35
end_zone2 = 0.75
end_zone3 = 0.95


T_zone1 = 0
T_zone2 = 150
T_zone3 = 0


T = [0, 0, T_zone2, T_zone2, 0, 0]
x =   [0, end_zone1, end_zone1, end_zone2,  end_zone2, end_zone3]



graph(x, T, r'Position $x$ [m]', r'Moment de torsion interne $T(x)$ [N$\cdot$m]')
plt.show()