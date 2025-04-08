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

        self.x_forces = {}
        self.y_forces = {}
        self.x_loads = {}
        self.y_loads = {}
        self.torsion = {}
        self.supports = {0 : 0, 1: 0.95}

    def add_force(self, F, l, angle=0):
        self.x_forces[l] = F*np.cos(np.deg2rad(angle))
        self.y_forces[l] = F*np.sin(np.deg2rad(angle))

    def add_load(self, F, l_i, l_f, angle=0):
        self.x_loads[(l_i, l_f)] = F*np.cos(np.deg2rad(angle))
        self.y_loads[(l_i, l_f)] = F*np.sin(np.deg2rad(angle))
        print(self.x_loads)
        print(self.y_loads)

    def add_torsion(self, T, l):
        self.torsion[l] = T

    def mod_support(self, i, l):
        self.supports[i] = l

    def moment_graph(self, V, graph):
        M = np.cumsum(V) * self.dx

        if graph:
            plt.figure(figsize=(10, 5))
            plt.plot(self.x, M, label="Moment de torsion M(x)", color='blue')
            plt.show()

        # return np.max(np.abs(M)), M
        return M

    def shear_graph(self, supports_reaction, forces, loads, graph):
        V = np.zeros(self.N)

        V[idx(self.x, self.supports[0]):] += supports_reaction[0]
        V[idx(self.x, self.supports[1]):] += supports_reaction[1]

        for pos, F in forces.items():
            V[idx(self.x, pos):] -= F

        for l, F in loads.items():
            for x_0 in self.x[idx(self.x, l[0]):idx(self.x, l[1])]:
                V[idx(self.x, x_0):] -= F / (l[1]-l[0]) *self.dx

        V[0] = 0
        V[-1] = 0

        if graph:
            plt.figure(figsize=(10, 5))
            plt.plot(self.x, V, label="Effort tranchant V(x)", color='blue')
            plt.show()

        max_V = np.max(np.abs(V))
        # max_M, M = self.moment_graph(V, graph)
        M = self.moment_graph(V, graph)

        # return max_V, V, max_M, M
        return V, M
    
    def torsion_graph(self, torsion, graph):
        T = np.zeros(self.N)

        for pos, t in torsion.items():
            T[idx(self.x, pos):] += t

        if graph:
            plt.figure(figsize=(10, 5))
            plt.plot(self.x, T, label="Torsion T(x)", color='blue')
            plt.show()

        # return np.max(np.abs(T))
        return T

    def solve(self, graph = False):
        pos = self.supports[1] - self.supports[0]

        # Forces en x
        moment_x = 0

        for l, F in self.x_forces.items():
                moment_x += F * (l - self.supports[0])

        for l, F in self.x_loads.items():
                moment_x += F * ((l[0]+l[1])/2 - self.supports[0])

        A_x = np.array([[1, 1], [0, pos]])
        b_x = np.array([sum(self.x_forces.values()) + sum(self.x_loads.values()), moment_x])

        self.supports_reaction_x = np.linalg.solve(A_x, b_x)

        # self.max_V_x, self.max_M_x, M_x = self.shear_graph(self.supports_reaction_x, self.x_forces, graph)
        self.V_x, self.M_x = self.shear_graph(self.supports_reaction_x, self.x_forces, self.x_loads, graph)

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
        self.V_y, self.M_y = self.shear_graph(self.supports_reaction_y, self.y_forces, self.y_loads, graph)

        self.V = np.sqrt(self.V_x**2 + self.V_y**2)
        self.M = np.sqrt(self.M_x**2 + self.M_y**2)

        plt.figure(figsize=(10, 5))
        plt.plot(self.x, self.V, label="Torsion T(x)", color='blue')
        plt.show()

        plt.figure(figsize=(10, 5))
        plt.plot(self.x, self.M, label="Torsion T(x)", color='blue')
        plt.show()

        # Torsion
        self.max_T = self.torsion_graph(self.torsion, False)

        # return self.max_V_x, self.max_M_x, self.max_V_y, self.max_M_y, self.max_T
        return self.V, self.M, self.max_T

h_B = 0.05
h_C = 0.075
F_B = 3000
F_C = h_B/h_C*F_B

print(F_C)

test = Poutre()

test.L = 0.95

test.add_force(F_B, 0.2)
test.add_force(F_C, 0.6, 90)
test.add_load(51.71, 0, 0.95, 90)
test.add_load(17.61, 0.2-0.032/2, 0.2+0.032/2, 90)
test.add_load(41.41, 0.6-0.032/2, 0.6+0.032/2, 90)
# test.add_load(1000, 0.3, 0.5, 90)
test.add_torsion(h_B*F_B, 0.2)
test.add_torsion(-h_C*F_C, 0.6)

test.solve(graph = True)
# print(test.solve(graph = True))
