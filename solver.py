from funcoesTermosol import *
import numpy as np
from math import degrees
from metodoJacob_Gauss import met_gauss

class Node:
    def __init__(self, id, x, y, free_degrees) -> float:
        self.id = id
        self.x = x
        self.y = y
        self.free_degrees = free_degrees

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

    def __calc_internal_tension(self, u) -> float:
        c, s = self.calc_c_s()
        T = np.array([-c, -s, c, s])
        return np.dot((self.y_modulos / self.dist) * T, u)

    def calc_internal_tension_and_force(self, u):
        internal_tension = self.__calc_internal_tension(u)
        return (internal_tension, internal_tension * self.area)
    
    def calc_internal_deformation(self, u):
        c, s = self.calc_c_s()
        T = np.array([-c, -s, c, s])
        return np.dot((1 / self.dist) * T, u)


    def degrees_of_freedom(self):
        return (*self.n1.free_degrees, *self.n2.free_degrees)

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
        self.nodes = []
        db = 0
        for i in range(self.nn):
            n = Node(i + 1, self.N[0, i], self.N[1, i], (db, db + 1))
            self.nodes.append(n)
            db += 2

        for i in range(self.nm):
            n1 = self.nodes[int(self.Inc[i, 0] - 1)]
            n2 = self.nodes[int(self.Inc[i, 1] - 1)]
            self.elements.append(Element(n1, n2, self.Inc[i, 2], self.Inc[i, 3], i+1))
    
    def calc_global_rigidity_matrix(self) -> np.array:
        K = np.zeros((self.nn * 2, self.nn * 2))
        for element in self.elements:
            K_e = element.calc_rigidity_matrix()
            # element_indexes = element.element_indexes(self.nn)
            element_indexes = element.degrees_of_freedom()

            for i, g_i in enumerate(element_indexes):
                for j, g_j in enumerate(element_indexes):
                    K[g_i, g_j] += K_e[i, j]
        return K

    def apply_restraints(self, K) -> np.array:
        """
        Deletes the rows and columns of the rigidity matrix that correspond to the restrained nodes
        """
        R = self.R.astype(int)
        KgR = np.delete(K, R, axis=0)
        KgR = np.delete(KgR, R, axis=1)
        return KgR

    def calc_displacement(self, k_g) -> np.array:
        """
        Compute the displacement of the nodes: U = Kg^-1 * F
        """
        F = self.F
        R = self.R.astype(int)
        F = np.delete(F, R)
        res, err = met_gauss(1e5, 1e-9, k_g, F)
        print(f"{err =}")
        # return np.linalg.solve(k_g, F)
        return res


    def full_displacement_matrix(self, u) -> np.array:
        """
        Add the displacement of the restrained nodes to the displacement matrix
        """
        u_full = np.zeros((self.nn * 2, 1))
        j = 0
        for i in range(self.nn * 2):
            if i in self.R:
                u_full[i] = 0
            else:
                u_full[i] = u[j]
                j += 1
        return u_full
    
    def calc_reactions(self, U, K) -> np.array:
        """
        Compute the reactions of the full F vector, including the reactions, that we wanted to discover: F = K * U
        """
        U_linha = np.zeros((self.nn * 2, 1))
        for i in range(self.nn * 2):
            if i not in self.R:
                U_linha[i] = U[i - len(self.R[self.R < i])]
        return np.matmul(K, U_linha)

    def calc_internal_tensions_forces(self, U) -> np.array:
        tensions = np.zeros(self.nm)
        forces = np.zeros(self.nm)
        for i, element in enumerate(self.elements):
            df = element.degrees_of_freedom()
            valores = np.array([U[df[0], 0], U[df[1], 0], U[df[2], 0], U[df[3], 0]])
            tensions[i], forces[i] = element.calc_internal_tension_and_force(valores)
        return tensions, forces

    def calc_internal_deformations(self, U) -> np.array:
        deformations = np.zeros(self.nm)
        for i, element in enumerate(self.elements):
            df = element.degrees_of_freedom()
            valores = np.array([U[df[0], 0], U[df[1], 0], U[df[2], 0], U[df[3], 0]])
            deformations[i] = element.calc_internal_deformation(valores)
        return deformations

    def plot_displacement(self, U):
        # New node matrix
        N = np.zeros((2, self.nn))
        for i in range(self.nn):
            N[0, i] = self.N[0, i] + U[2 * i]
            N[1, i] = self.N[1, i] + U[2 * i + 1]
        plota(N, self.Inc)


    def solve(self):
        plota(self.N,self.Inc)
        # Calcula a matriz de rigidez global
        K = self.calc_global_rigidity_matrix()
        # Aplica as condicoes de contorno
        rest = self.apply_restraints(K)
        # Calcula a matriz de deslocamentos
        U = self.calc_displacement(rest)
        full_u = self.full_displacement_matrix(U)
        # Calcula as reacoes
        F = self.calc_reactions(U, K)
        Ft = np.zeros((self.nr, 1))
        for i in range(self.nr):
            Ft[i, 0] = F[int(self.R[i]), 0]

        int_tensions, int_forces = self.calc_internal_tensions_forces(full_u)
        int_deformations = self.calc_internal_deformations(full_u)

        geraSaida(self.out_path, Ft, full_u, int_deformations, int_forces, int_tensions)
        self.plot_displacement(full_u)        



if __name__ == "__main__":
    solver = Solver("entrada2.xlsx", "saida.txt")
    solver.solve()