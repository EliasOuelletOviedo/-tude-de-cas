import numpy as np
import matplotlib.pyplot as plt

def idx(array, value):
    array = np.asarray(array)
    idx_val = (np.abs(array - value)).argmin()
    return idx_val

class Poutre:
    def __init__(self, N=10000):
        self.E = 205e9  # Young's modulus in Pa
        self.G = 80e9   # Shear modulus in Pa
        self.d = 0.3    # Diameter in m
        self.L = 0.95   # Length in m

        self.N = N
        self.x, self.dx = np.linspace(0, self.L, self.N, retstep=True)

        self.z_forces = {}
        self.y_forces = {}
        self.z_loads = {}
        self.y_loads = {}
        self.torsion = {}
        self.supports = {0 : 0, 1: 0.95}

    def add_force(self, F, l, angle=0):
        self.z_forces[l] = F*np.cos(np.deg2rad(angle))
        self.y_forces[l] = F*np.sin(np.deg2rad(angle))

    def add_load(self, F, l_i, l_f, angle=0):
        self.z_loads[(l_i, l_f)] = F*np.cos(np.deg2rad(angle))
        self.y_loads[(l_i, l_f)] = F*np.sin(np.deg2rad(angle))

    def add_torsion(self, T, l):
        self.torsion[l] = T

    def mod_support(self, i, l):
        self.supports[i] = l

    def moment_graph(self, V):
        M = np.cumsum(V) * self.dx

        M[0] = 0
        M[-1] = 0

        return M

    def shear_graph(self, supports_reaction, forces, loads):
        V = np.zeros(self.N)

        V[idx(self.x, self.supports[0]):] += supports_reaction[0]
        V[idx(self.x, self.supports[1]):] += supports_reaction[1]

        for pos, F in forces.items():
            V[idx(self.x, pos):] -= F

        for l, F in loads.items():
            for x_0 in self.x[idx(self.x, l[0]):idx(self.x, l[1])]:
                V[idx(self.x, x_0):] -= F / (l[1]-l[0]) * self.dx

        V[0] = 0
        V[-1] = 0

        M = self.moment_graph(V)

        return V, M
    
    def torsion_graph(self, torsion):
        T = np.zeros(self.N)

        for pos, t in torsion.items():
            T[idx(self.x, pos):] += t

        return T

    def solve(self, graph = False):
        pos = self.supports[1] - self.supports[0]

        # Forces en x
        moment_z = 0

        for l, F in self.z_forces.items():
                moment_z += F * (l - self.supports[0])

        for l, F in self.z_loads.items():
                moment_z += F * ((l[0]+l[1])/2 - self.supports[0])

        A_z = np.array([[1, 1], [0, pos]])
        b_z = np.array([sum(self.z_forces.values()) + sum(self.z_loads.values()), moment_z])

        self.supports_reaction_z = np.linalg.solve(A_z, b_z)

        # self.max_V_z, self.max_M_z, M_z = self.shear_graph(self.supports_reaction_z, self.z_forces, graph)
        self.V_z, self.M_y = self.shear_graph(self.supports_reaction_z, self.z_forces, self.z_loads)

        # Forces en y
        moment_y = 0

        for l, F in self.y_forces.items():
                moment_y += F * (l - self.supports[0])

        for l, F in self.y_loads.items():
                moment_y += F * ((l[0]+l[1])/2 - self.supports[0])

        A_y = np.array([[1, 1], [0, pos]])
        b_y = np.array([sum(self.y_forces.values()) + sum(self.y_loads.values()), moment_y])

        self.supports_reaction_y = np.linalg.solve(A_y, b_y)

        # self.max_V_y, self.max_M_y, M_y = self.shear_graph(self.supports_reaction_y, self.y_forces, graph)
        self.V_y, self.M_z = self.shear_graph(self.supports_reaction_y, self.y_forces, self.y_loads)

        self.V = np.sqrt(self.V_z**2 + self.V_y**2)
        self.theta_V = np.degrees(np.arctan2(self.V_z, self.V_y))
        # self.theta_V[self.theta_V<0] = 0
        # self.theta_V[self.theta_V<0] += 360
        self.M = np.sqrt(self.M_z**2 + self.M_y**2)
        self.theta_M = np.degrees(np.arctan2(self.M_z, self.M_y))
        self.theta_M[self.theta_M<0] = 0
        # self.theta_M[self.theta_M<0] += 360

        fig, axs = plt.subplots(ncols=2, nrows=4, figsize=(16, 8),layout="constrained")

        for ax in axs.flatten():
            ax.axhline(0, color='gray', linewidth=1)
            ax.axvline(0.35, color='gray', linewidth=1, linestyle='--')
            ax.axvline(0.75, color='gray', linewidth=1, linestyle='--')
        
        axs[0, 0].plot(self.x, self.V_z, label="V_z(x)", color='blue')
        axs[0, 0].tick_params(axis='x', labelbottom=False)
        axs[0, 0].set_ylabel(r"Effort tranchant en $z$ (N)")

        axs[1, 0].plot(self.x, self.V_y, label="M_z(x)", color='blue')
        axs[1, 0].tick_params(axis='x', labelbottom=False)
        axs[1, 0].set_ylabel(r"Effort tranchant en $y$ (N)")

        axs[2, 0].plot(self.x, self.V, label="V_y(x)", color='blue')
        axs[2, 0].tick_params(axis='x', labelbottom=False)
        axs[2, 0].set_ylabel(r"Effort tranchant (N)")

        axs[3, 0].plot(self.x, self.theta_V, label="theta_V(x)", color='blue')
        axs[3, 0].set_xlabel(r"Position $x$ (m)")
        axs[3, 0].set_ylabel(r"Angle de $V(x)$ (°)")

        axs[0, 1].plot(self.x, self.M_z, label="M_z(x)", color='blue')
        axs[0, 1].tick_params(axis='x', labelbottom=False)
        axs[0, 1].set_ylabel(r"Moment de flexion en $z$ (N$\cdot$m)")

        axs[1, 1].plot(self.x, self.M_y, label="M_y(x)", color='blue')
        axs[1, 1].tick_params(axis='x', labelbottom=False)
        axs[1, 1].set_ylabel(r"Moment de flexion en $y$ (N$\cdot$m)")

        axs[2, 1].plot(self.x, self.M, label="M(x)", color='blue')
        axs[2, 1].tick_params(axis='x', labelbottom=False)
        axs[2, 1].set_ylabel(r"Moment de flexion (N$\cdot$m)")

        axs[3, 1].plot(self.x, self.theta_M, label="theta_M(x)", color='blue')
        axs[3, 1].set_xlabel(r"Position $x$ (m)")
        axs[3, 1].set_ylabel(r"Angle de $M(x)$ (°)")

        plt.show()

        # Torsion
        self.max_T = self.torsion_graph(self.torsion)

        # return self.max_V_z, self.max_M_z, self.max_V_y, self.max_M_y, self.max_T
        return np.max(self.V), np.max(self.M), np.max(self.max_T)

P = 8500
rpm = 541
T = P / (2 * np.pi * rpm / 60)

h_poutre = 0.03
L_poutre = 0.95

e_gear = 0.032
h_B = 0.05
h_C = 0.075
F_B = T/h_B
F_C = T/h_C

rho = 7850

P_poutre = 7850*np.pi*h_poutre**2*L_poutre
P_B = 7850*np.pi*h_B**2*e_gear
P_C = 7850*np.pi*h_C**2*e_gear

# print(P_B, P_C, P_poutre)
# print(F_B, F_C)

test = Poutre()

test.L = 0.95

test.add_force(F_B, 0.75, 90)
test.add_force(F_C, 0.35)

# Si on considère la largeur des engrenages et le poids de la poutre
# test.add_load(P_poutre, 0, 0.95)
# test.add_load(P_B, 0.75-e_gear/2, 0.75+e_gear/2)
# test.add_load(P_C, 0.35-e_gear/2, 0.35+e_gear/2)
# test.add_load(F_B, 0.75-e_gear/2, 0.75+e_gear/2, 90)
# test.add_load(F_C, 0.35-e_gear/2, 0.35+e_gear/2)

test.add_torsion(h_B*F_B, 0.2)
test.add_torsion(-h_C*F_C, 0.6)

# test.solve(graph = True)
print(test.solve(graph = True))
