"""Node class for the nodes of the Truss."""

import matplotlib.pyplot as plt


class Node:
    """
    Node class for the nodes of the Truss.

    Attributes:
    id: int
        Node id.
    x: float
        x coordinate.
    y: float
        y coordinate.
    free_degrees: tuple
        Degrees of freedom of the node.
    """

    def __init__(self, node_id, x, y, free_degrees) -> float:
        self.node_id = node_id
        self.x = x
        self.y = y
        self.free_degrees = free_degrees

    def display(self, displacement, color="r"):
        """
        Display the node in the plot.
        """
        disp_x, disp_y = displacement
        plt.plot(self.x + disp_x, self.y + disp_y, "o", color=color)

    def __str__(self) -> str:
        return f"Node({self.node_id}, x:{self.x}, y:{self.y})"
