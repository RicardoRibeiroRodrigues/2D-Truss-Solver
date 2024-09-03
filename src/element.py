"""Element Class module."""

import matplotlib.pyplot as plt
import numpy as np
from .node import Node


class Element:
    """
    Define the element class.

    Attributes:
    n1: Node
        Initial node.
    n2: Node
        Final node.
    young_modulus: float
        Young modulus of the element.
    area: float
        Area of the element.
    element_id: int
        Element id.
    """

    current_id = 0

    def _get_id(self):
        Element.current_id += 1
        return Element.current_id

    def __init__(self, n1: Node, n2: Node, young_modulus: float, area: float) -> float:
        self.n1 = n1
        self.n2 = n2
        self.y_modulos = young_modulus
        self.area = area
        self.element_id = self._get_id()
        self.dist = self.calc_dist()

    def display(self, color="r", linewidth=3, line_style="--"):
        """
        Display the element in the plot.
        """
        plt.plot(
            (self.n1.x, self.n2.x),
            (self.n1.y, self.n2.y),
            line_style,
            color=color,
            linewidth=linewidth,
        )

    def __str__(self) -> str:
        return f"Elemento: {self.n1}<->{self.n2}"

    def calc_dist(self) -> float:
        """
        Euclidean distance between the nodes.
        """
        return ((self.n2.x - self.n1.x) ** 2 + (self.n2.y - self.n1.y) ** 2) ** 0.5

    def calc_cos_sin(self) -> tuple:
        """
        Method to calculate the cos and sin of the element.
        """
        cos = (self.n2.x - self.n1.x) / self.dist
        sin = (self.n2.y - self.n1.y) / self.dist
        return cos, sin

    def calc_rigidity_matrix(self) -> np.array:
        """
        Method to calculate the rigidity matrix of the element
        """
        cos, sin = self.calc_cos_sin()
        rigidity_matrix = np.array(
            [
                [cos**2, cos * sin, -(cos**2), -cos * sin],
                [cos * sin, sin**2, -cos * sin, -(sin**2)],
                [-(cos**2), -cos * sin, cos**2, cos * sin],
                [-cos * sin, -(sin**2), cos * sin, sin**2],
            ]
        )
        return ((self.y_modulos * self.area) / self.dist) * rigidity_matrix

    def __calc_internal_tension(self, displacement) -> float:
        """
        Calculate the internal tension of the element.
        """
        cos, sin = self.calc_cos_sin()
        tension_vector = np.array([-cos, -sin, cos, sin])
        return np.dot((self.y_modulos / self.dist) * tension_vector, displacement)

    def calc_internal_tension_and_force(self, u):
        """
        Calculate the internal tension and force of the element.
        """
        internal_tension = self.__calc_internal_tension(u)
        return (internal_tension, internal_tension * self.area)

    def calc_internal_deformation(self, u):
        """
        Calculate the internal deformation of the element.
        """
        cos, sin = self.calc_cos_sin()
        tension_vector = np.array([-cos, -sin, cos, sin])
        return np.dot((1 / self.dist) * tension_vector, u)

    def degrees_of_freedom(self):
        """
        Get the degrees of freedom of the element.
        """
        return (*self.n1.free_degrees, *self.n2.free_degrees)
