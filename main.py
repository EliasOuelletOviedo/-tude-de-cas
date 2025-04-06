import numpy as np
import matplotlib.pyplot as plt

def idx(array, value):
    array = np.asarray(array)
    idx_val = (np.abs(array - value)).argmin()
    return idx_val

class Poutre:
    def __init__(self):
        self.E = 205e9  # Young's modulus in Pa
        self.G = 80e9   # Shear modulus in Pa
        self.d = 0.3    # Diameter in m
        self.L = 0.95   # Length in m

        self.x_forces = {}
        self.y_forces = {}
        self.torsion = {}
        self.supports = {0 : 0, 1: 0.95}

    def add_force(self, F, l, angle=0):
        self.x_forces[l] = F*np.cos(np.deg2rad(angle))
        self.y_forces[l] = F*np.sin(np.deg2rad(angle))

    def add_torsion(self, T, l):
        self.torsion[l] = T

    def mod_support(self, i, l):
        self.supports[i] = l

    def moment_graph(self, V, x, graph):
        dx = x[1] - x[0]
        M = np.cumsum(V) * dx

        if graph:
            plt.figure(figsize=(10, 5))
            plt.plot(x, M, label="Moment de torsion M(x)", color='blue')
            plt.show()

        return np.max(np.abs(M))

    def shear_graph(self, supports_reaction, forces, graph):
        x = np.linspace(0, self.L, 1000)
        V = np.zeros(1000)

        V[idx(x, self.supports[0]):] += supports_reaction[0]
        V[idx(x, self.supports[1]):] += supports_reaction[1]

        for pos, F in forces.items():
            V[idx(x, pos):] -= F

        V[0] = 0
        V[-1] = 0

        if graph:
            plt.figure(figsize=(10, 5))
            plt.plot(x, V, label="Effort tranchant V(x)", color='blue')
            plt.show()

        max_V = np.max(np.abs(V))
        max_M = self.moment_graph(V, x, graph)

        return max_V, max_M
    
    def torsion_graph(self, torsion, graph):
        x = np.linspace(0, self.L, 1000)
        T = np.zeros(1000)

        for pos, t in torsion.items():
            T[idx(x, pos):] += t

        if graph:
            plt.figure(figsize=(10, 5))
            plt.plot(x, T, label="Torsion T(x)", color='blue')
            plt.show()

        return np.max(np.abs(T))

    def solve(self, graph = False):
        pos = self.supports[1] - self.supports[0]

        # Forces en x
        moment_x = 0

        for l, F in self.x_forces.items():
                moment_x += F * (l - self.supports[0])

        A_x = np.array([[1, 1], [0, pos]])
        b_x = np.array([sum(self.x_forces.values()), moment_x])

        self.supports_reaction_x = np.linalg.solve(A_x, b_x)

        self.max_V_x, self.max_M_x = self.shear_graph(self.supports_reaction_x, self.x_forces, graph)

        # Forces en y
        moment_y = 0

        for l, F in self.y_forces.items():
                moment_y += F * (l - self.supports[0])

        A_y = np.array([[1, 1], [0, pos]])
        b_y = np.array([sum(self.y_forces.values()), moment_y])

        self.supports_reaction_y = np.linalg.solve(A_y, b_y)

        self.max_V_y, self.max_M_y = self.shear_graph(self.supports_reaction_y, self.y_forces, graph)

        # Torsion
        self.max_T = self.torsion_graph(self.torsion, graph)

        return self.max_V_x, self.max_M_x, self.max_V_y, self.max_M_y, self.max_T

h_B = 0.05
h_C = 0.075
F_B = 1000
F_C = h_B/h_C*F_B

test = Poutre()

test.L = 0.95

test.add_force(F_B, 0.2)
test.add_force(F_C, 0.6, 90)
test.add_torsion(h_B*F_B, 0.2)
test.add_torsion(-h_C*F_C, 0.6)

print(test.solve(graph = True))
