import numpy as np
from graphique_effort_interne import graph
import matplotlib.pyplot as plt
from scipy.integrate import cumtrapz

plt.rcParams['text.usetex'] = True


n_points_zone = 100
length = 0.95

# Calcul de M_y
end_zone1 = 0.75
end_zone2 = 0.95

v_y_zone1 = 632 
v_y_zone2 = -2368 

v_y = [0,v_y_zone1, v_y_zone1, v_y_zone2, v_y_zone2, v_y_zone2, 0]
x =   [0,0, end_zone1, end_zone1, end_zone2, end_zone2, length]

M_y = cumtrapz(v_y, x, initial=0)

x_upsampled = np.linspace(x[0], x[-1], 200)
M_y_upsampled = np.interp(x_upsampled, x, M_y)


# ---------------------------

# Calcul de M_z
v_z_zone1 = 1263
v_z_zone2 = -737

end_zone1 = 0.35
end_zone2 = 0.95

v_z = [0,v_z_zone1, v_z_zone1, v_z_zone2, v_z_zone2, v_z_zone2, 0]
x =   [0,0, end_zone1, end_zone1, end_zone2, end_zone2, length]

M_z = cumtrapz(v_z, x, initial=0)

x_upsampled = np.linspace(x[0], x[-1], 200)
M_z_upsampled = np.interp(x_upsampled, x, M_z)

normeM = np.sqrt(M_z_upsampled**2 + M_y_upsampled**2)

graph(x_upsampled, normeM, r'Position $x$ [m]', r'Norme du moment $M(x)$ [N$\cdot$m]')
# graph(x, M_z ,r'Position $x$ [m]', r'Norme du moment $M(x)$ [N$\cdot$m]')
# graph(x, M_y ,r'Position $x$ [m]', r'Norme du moment $M(x)$ [N$\cdot$m]')
plt.show()