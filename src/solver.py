"""
Main module for the truss solver, contains all logic to solve the truss problem.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from .numeric_method import gauss_method
from .truss_data_io import import_data, generate_output, SolverData
from .node import Node
from .element import Element


class Solver:
    """
    Solver class for the truss.

    Solves the truss problem and generates the output file with the solve method.

    Just instantiate the class and call the solve method.
    """

    elements: list[Element]
    nodes: list[Node]

    def __init__(self, in_path) -> None:
        variables: SolverData = import_data(in_path)
        self._load_variables(variables)
        self._create_elements()

    def _load_variables(self, variables: SolverData) -> None:
        self.nodes_number = variables.nodes_number
        self.nodes_matrix = variables.nodes
        self.elements_number = variables.elements_number
        self.incidence = variables.incidence
        self.loads_number = variables.loads_number
        self.loads = variables.loads
        self.restraints_number = variables.restraints_number
        self.restraints = variables.restraints

    def _create_elements(self):
        self.elements = []
        self.nodes = []
        db = 0
        for i in range(self.nodes_number):
            node = Node(
                i + 1, self.nodes_matrix[0, i], self.nodes_matrix[1, i], (db, db + 1)
            )
            self.nodes.append(node)
            db += 2

        for i in range(self.elements_number):
            n1 = self.nodes[int(self.incidence[i, 0] - 1)]
            n2 = self.nodes[int(self.incidence[i, 1] - 1)]
            self.elements.append(
                Element(n1, n2, self.incidence[i, 2], self.incidence[i, 3])
            )

    def calc_global_rigidity_matrix(self) -> np.array:
        """
        Compute the global rigidity matrix of the truss.
        """
        global_rigidity_matrix = np.zeros(
            (self.nodes_number * 2, self.nodes_number * 2)
        )
        for element in self.elements:
            element_rigidity_matrix = element.calc_rigidity_matrix()
            element_indexes = element.degrees_of_freedom()

            for i, g_i in enumerate(element_indexes):
                for j, g_j in enumerate(element_indexes):
                    global_rigidity_matrix[g_i, g_j] += element_rigidity_matrix[i, j]
        return global_rigidity_matrix

    def apply_restraints(self, global_rigidity_matrix) -> np.array:
        """
        Deletes the rows and columns of the rigidity matrix that correspond to the restrained nodes
        """
        restraints = self.restraints.astype(int)
        updated_rigidity_matrix = np.delete(global_rigidity_matrix, restraints, axis=0)
        updated_rigidity_matrix = np.delete(updated_rigidity_matrix, restraints, axis=1)
        return updated_rigidity_matrix

    def calc_displacement(self, k_g) -> np.array:
        """
        Compute the displacement of the nodes: U = Kg^-1 * F
        """
        loads = self.loads
        restraints = self.restraints.astype(int)
        loads = np.delete(loads, restraints)
        displacement, _ = gauss_method(1e5, 1e-9, k_g, loads)
        return displacement

    def full_displacement_matrix(self, node_displacements) -> np.array:
        """
        Add the displacement of the restrained nodes to the displacement matrix
        """
        displacement_matrix = np.zeros((self.nodes_number * 2, 1))
        j = 0
        for i in range(self.nodes_number * 2):
            if i in self.restraints:
                displacement_matrix[i] = 0
            else:
                displacement_matrix[i] = node_displacements[j]
                j += 1
        return displacement_matrix

    def calc_reactions(self, displacement, stiffness_matrix) -> np.array:
        """
        Compute the reactions of the full F vector, including the reactions:
            -> F = K * U
        """
        updated_displacements = np.zeros((self.nodes_number * 2, 1))
        for i in range(self.nodes_number * 2):
            if i not in self.restraints:
                updated_displacements[i] = displacement[
                    i - len(self.restraints[self.restraints < i])
                ]
        return np.matmul(stiffness_matrix, updated_displacements)

    def _calc_internal_tensions_forces(self, displacements) -> np.array:
        tensions = np.zeros(self.elements_number)
        forces = np.zeros(self.elements_number)
        for i, element in enumerate(self.elements):
            df = element.degrees_of_freedom()
            valores = np.array(
                [
                    displacements[df[0], 0],
                    displacements[df[1], 0],
                    displacements[df[2], 0],
                    displacements[df[3], 0],
                ]
            )
            tensions[i], forces[i] = element.calc_internal_tension_and_force(valores)
        return tensions, forces

    def _calc_internal_deformations(self, displacements) -> np.array:
        deformations = np.zeros(self.elements_number)
        for i, element in enumerate(self.elements):
            df = element.degrees_of_freedom()
            valores = np.array(
                [
                    displacements[df[0], 0],
                    displacements[df[1], 0],
                    displacements[df[2], 0],
                    displacements[df[3], 0],
                ]
            )
            deformations[i] = element.calc_internal_deformation(valores)
        return deformations

    def plot_displacement(self, displacements, save_plot: bool):
        """
        Method to plot the displacement of the truss.
        """
        plt.style.use("seaborn-v0_8")
        _ = plt.figure()
        # Passa por todos os membros
        for element in self.elements:
            element.display()
            element.n1.display((0, 0))
            element.n2.display((0, 0))
            # pos deformacao
            element.n1.display(displacements[element.n1.free_degrees, 0], color="b")
            element.n1.display(displacements[element.n2.free_degrees, 0], color="b")
            plt.plot(
                (
                    element.n1.x + displacements[element.n1.node_id - 1],
                    element.n2.x + displacements[element.n2.node_id - 1],
                ),
                (
                    element.n1.y + displacements[element.n2.node_id],
                    element.n2.y + displacements[element.n2.node_id],
                ),
                "--",
                color="b",
                linewidth=3,
            )

        plt.xlabel("x [m]")
        plt.ylabel("y [m]")
        red_patch = mpatches.Patch(color="red", label="Antes da aplicação da carga")
        blue_patch = mpatches.Patch(color="blue", label="Depois da aplicação da carga")
        plt.legend(handles=[red_patch, blue_patch])
        plt.grid(True)
        plt.axis("equal")
        save_plot and plt.savefig("imgs/displacement.png", dpi=300)
        plt.show()

    def plot_internal_tensions(self, tensions, save_plot: bool):
        """
        Method to plot the internal tensions in each element.
        """
        plt.style.use("seaborn-v0_8")
        _ = plt.figure()
        # Scale the color of the element according to the tension
        for i, element in enumerate(self.elements):
            element.display(color="deepskyblue", linewidth=3)
            # Plot the point too
            element.n1.display((0, 0), color="black")
            element.n2.display((0, 0), color="black")
            ang_1 = (element.n2.y - element.n1.y) / element.dist
            ang_2 = (element.n2.x - element.n1.x) / element.dist

            plt.text(
                (element.n1.x + element.n2.x) / 2,
                (element.n1.y + element.n2.y) / 2,
                f"{tensions[i]:.2e} N",
                rotation=(180 * np.arctan(ang_1 / ang_2)) / np.pi,
                color="black",
                horizontalalignment="center",
                verticalalignment="center",
                fontsize=12,
                rotation_mode="anchor",
                transform_rotates_text=True,
            )
        plt.xlabel("x [m]")
        plt.ylabel("y [m]")
        plt.title("Tensões internas nos elementos")
        plt.grid(True)
        save_plot and plt.savefig("imgs/tensões_internas.png", dpi=300)
        plt.show()

    def solve(self, output_path: str, save_plots: bool = False) -> None:
        """
        Top level method to solve the truss problem.
        """
        # Calcula a matriz de rigidez global
        global_rigidity_matrix = self.calc_global_rigidity_matrix()
        # Aplica as condicoes de contorno
        rest = self.apply_restraints(global_rigidity_matrix)
        # Calcula a matriz de deslocamentos
        displacement_matrix = self.calc_displacement(rest)
        full_displacement = self.full_displacement_matrix(displacement_matrix)
        # Calcula as reacoes
        reactions = self.calc_reactions(displacement_matrix, global_rigidity_matrix)
        restraint_reactions = np.zeros((self.restraints_number, 1))
        for i in range(self.restraints_number):
            restraint_reactions[i, 0] = reactions[int(self.restraints[i]), 0]

        int_tensions, int_forces = self._calc_internal_tensions_forces(
            full_displacement
        )
        int_deformations = self._calc_internal_deformations(full_displacement)

        generate_output(
            output_path,
            restraint_reactions,
            full_displacement,
            int_deformations,
            int_forces,
            int_tensions,
        )
        self.plot_displacement(full_displacement, save_plots)
        self.plot_internal_tensions(int_tensions, save_plots)
