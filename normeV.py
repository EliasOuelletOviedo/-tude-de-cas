import numpy as np
from graphique_effort_interne import graph
import matplotlib.pyplot as plt
from scipy.integrate import cumtrapz

plt.rcParams['text.usetex'] = True


n_points_zone = 100
length = 0.95

# Calcul de V_y
end_zone1 = 0.75
end_zone2 = 0.95

v_y_zone1 = 632 
v_y_zone2 = -2368 

v_y = [0,v_y_zone1, v_y_zone1, v_y_zone2, v_y_zone2, v_y_zone2, 0]
x =   [0,0, end_zone1, end_zone1, end_zone2, end_zone2, length]

x_upsampled = np.linspace(x[0], x[-1], 500)
v_y_upsampled = np.interp(x_upsampled, x, v_y)


# ---------------------------

# Calcul de V_z
v_z_zone1 = 1263
v_z_zone2 = -737

end_zone1 = 0.35
end_zone2 = 0.95

v_z = [0,v_z_zone1, v_z_zone1, v_z_zone2, v_z_zone2, v_z_zone2, 0]
x =   [0,0, end_zone1, end_zone1, end_zone2, end_zone2, length]


x_upsampled = np.linspace(x[0], x[-1], 500)
v_z_upsampled = np.interp(x_upsampled, x, v_z)

normeV = np.sqrt(v_z_upsampled**2 + v_y_upsampled**2)


graph(x_upsampled, normeV, r'Position $x$ [m]', r'Norme de $V(x)$ [N]')
plt.plot([0,0],[0,normeV[0]], color="orange")
plt.show()