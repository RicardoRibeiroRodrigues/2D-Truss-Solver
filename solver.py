from funcoesTermosol import *
import numpy as np
from math import degrees

class Node:
    def __init__(self, id, x, y) -> float:
        self.id = id
        self.x = x
        self.y = y

    def __str__(self) -> str:
        return f'Node({self.id}, x:{self.x}, y:{self.y})'

class Element:
    def __init__(self, n1, n2, young_modulus, area, id) -> float:
        self.n1 = n1
        self.n2 = n2
        self.y_modulos = young_modulus
        self.area = area
        self.id = id
        self.dist = self.calc_dist()
    
    def __str__(self) -> str:
        return f"Elemento: {self.n1}<->{self.n2}"

    def calc_dist(self) -> float:
        return ((self.n2.x - self.n1.x) **2 + (self.n2.y - self.n1.y) **2) ** 0.5
    
    def calc_c_s(self) -> tuple:
        c = (self.n2.x - self.n1.x) / self.dist 
        s = (self.n2.y - self.n1.y) / self.dist
        return c, s

    def calc_rigidity_matrix(self) -> np.array:
        c, s = self.calc_c_s()
        M = np.array([
                    [c**2, c*s, -c**2, -c*s],
                    [c*s, s**2, -c*s, -s**2],
                    [-c**2, -c*s, c**2, c*s],
                    [-c*s, -s**2, c*s, s**2]
                ])
        return ((self.y_modulos * self.area) / self.dist) * M

    def element_indexes(self, nn) -> np.array:
        indexes_list = []
        for i in range(self.id * 2 - 2, self.id * 2 + 2):
            num = i % (nn * 2)
            indexes_list.append(num)
        
        return np.array(indexes_list)

class Solver:

    def __init__(self, in_path, out_path) -> None:
        self.out_path = out_path
        variables = importa(in_path)
        self.nn = variables[0]
        self.N = variables[1]
        self.nm = variables[2]
        self.Inc = variables[3]
        self.nc = variables[4]
        self.F = variables[5]
        self.nr = variables[6]
        self.R = variables[7]
        self.enumerate_elements()

    def enumerate_elements(self):
        self.elements = []
        for i in range(self.nm):
            n1 = Node(self.Inc[i, 0], self.N[0, int(self.Inc[i, 0]-1)], self.N[1, int(self.Inc[i, 0]-1)])
            n2 = Node(self.Inc[i, 1], self.N[0, int(self.Inc[i, 1]-1)], self.N[1, int(self.Inc[i, 1]-1)])
            self.elements.append(Element(n1, n2, self.Inc[i, 2], self.Inc[i, 3], i+1))
    
    def calc_global_rigidity_matrix(self) -> np.array:
        K = np.zeros((self.nn * 2, self.nn * 2))
        for element in self.elements:
            K_e = element.calc_rigidity_matrix()
            element_indexes = element.element_indexes(self.nn)

            for i, g_i in enumerate(element_indexes):
                for j, g_j in enumerate(element_indexes):
                    K[g_i, g_j] += K_e[i, j]
        return K

    def apply_restraints(self, K) -> np.array:
        """
        Deletes the rows and columns of the rigidity matrix that correspond to the restrained nodes
        """
        for i, restrain in enumerate(self.R):
            K = np.delete(K, int(restrain - i), 1)
            K = np.delete(K, int(restrain - i), 0)
        return K

    def calc_displacement(self, k_g) -> np.array:
        """
        Compute the displacement of the nodes: U = Kg^-1 * F
        """
        F = self.F
        for i, restrain in enumerate(self.R):
            F = np.delete(F, int(restrain - i), 0)
        return np.linalg.solve(k_g, F)
    
    def calc_reactions(self, U, K) -> np.array:
        """
        Compute the reactions of the full F vector, including the reactions, that we wanted to discover: F = K * U
        """
        U_linha = np.zeros((self.nn * 2, 1))
        for i in range(self.nn * 2):
            if i not in self.R:
                U_linha[i] = U[i - len(self.R[self.R < i])]
        return np.matmul(K, U_linha)

            



if __name__ == "__main__":
    solver = Solver("entrada.xlsx", "saida.txt")
    # plota(solver.N, solver.Inc)
    K = solver.calc_global_rigidity_matrix()
    rest = solver.apply_restraints(K)
    U = solver.calc_displacement(rest)
    F = solver.calc_reactions(U, K)
    for i in range(len(solver.R)):
        print(f"R{i + 1} = {F[int(solver.R[i])]}")