from graphique_effort_interne import graph
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['text.usetex'] = True

n_points_zone = 100
length = 0.95

v_z_zone1 = 1263
v_z_zone2 = -737

end_zone1 = 0.35
end_zone2 = 0.95

v_z = [0,v_z_zone1, v_z_zone1, v_z_zone2, v_z_zone2, v_z_zone2, 0]
x =   [0,0, end_zone1, end_zone1, end_zone2, end_zone2, length]

graph(x, v_z, r'Position $x$ [m]', r'Effort tranchant $V_z(x)$ [N]')
plt.show()