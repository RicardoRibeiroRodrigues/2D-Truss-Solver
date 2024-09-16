from src import Node

def test_node_initialization():
    node = Node(node_id=1, x=5.0, y=10.0, free_degrees=(0, 1))
    assert node.node_id == 1
    assert node.x == 5.0
    assert node.y == 10.0
    assert node.free_degrees == (0, 1)
