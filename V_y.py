from graphique_effort_interne import graph
import matplotlib.pyplot as plt
import numpy as np


plt.rcParams['text.usetex'] = True

n_points_zone = 100
end_zone1 = 0.75
end_zone2 = 0.95
length = 0.95

v_y_zone1 = 632 
v_y_zone2 = -2368 

v_y = [0,v_y_zone1, v_y_zone1, v_y_zone2, v_y_zone2, v_y_zone2, 0]
x =   [0,0, end_zone1, end_zone1, end_zone2, end_zone2, length]

print(max(np.abs(v_y)))
graph(x, v_y, r'Position $x$ [m]', r'Effort tranchant $V_y(x)$ [N]')
plt.show()