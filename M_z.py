import numpy as np
from graphique_effort_interne import graph
import matplotlib.pyplot as plt
from scipy.integrate import cumtrapz

plt.rcParams['text.usetex'] = True


n_points_zone = 100
length = 0.95

v_z_zone1 = 1263
v_z_zone2 = -737

end_zone1 = 0.35
end_zone2 = 0.95

v_z = [0,v_z_zone1, v_z_zone1, v_z_zone2, v_z_zone2, v_z_zone2, 0]
x =   [0,0, end_zone1, end_zone1, end_zone2, end_zone2, length]

M_z = cumtrapz(v_z, x, initial=0)

x_upsampled = np.linspace(x[0], x[-1], 200)
M_z_upsampled = np.interp(x_upsampled, x, M_z)

graph(x_upsampled, M_z_upsampled, r'Position $x$ [m]', r'Moment de flexion $M_z(x)$ [N$\cdot$m]')
plt.show()