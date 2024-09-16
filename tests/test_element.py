from src import Node
from src import Element

def test_element_distance_calculation():
    node1 = Node(node_id=1, x=0, y=0, free_degrees=(0, 1))
    node2 = Node(node_id=2, x=3, y=4, free_degrees=(2, 3))  # Distance should be 5 (3-4-5 triangle)
    element = Element(n1=node1, n2=node2, young_modulus=210e9, area=1e-4)
    
    assert element.calc_dist() == 5.0  # Check if the calculated distance is correct
