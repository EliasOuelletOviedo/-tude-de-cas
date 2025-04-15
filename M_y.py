import numpy as np
from graphique_effort_interne import graph
import matplotlib.pyplot as plt
from scipy.integrate import cumtrapz


plt.rcParams['text.usetex'] = True

n_points_zone = 100
end_zone1 = 0.75
end_zone2 = 0.95
length = 0.95

v_y_zone1 = 632 
v_y_zone2 = -2368 

v_y = [0,v_y_zone1, v_y_zone1, v_y_zone2, v_y_zone2, v_y_zone2, 0]
x =   [0,0, end_zone1, end_zone1, end_zone2, end_zone2, length]

M_y = cumtrapz(v_y, x, initial=0)

x_upsampled = np.linspace(x[0], x[-1], 200)
M_y_upsampled = np.interp(x_upsampled, x, M_y)

graph(x_upsampled, M_y_upsampled, r'Position $x$ [m]', r'Moment de flexion $M_y(x)$ [N$\cdot$m]')
plt.show()