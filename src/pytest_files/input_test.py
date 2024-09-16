import pytest
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from truss_data_io import *

def test_input_1():
    current_dir = os.path.dirname(__file__)
    data_dir = os.path.abspath(os.path.join(current_dir, '../../data'))
    excel_file_path = os.path.join(data_dir, 'entrada2.xlsx')

    solver_data = import_data(excel_file_path)

    assertion_object = SolverData(
        7,
        [[0., 0.144, 0.192, 0.288, 0.384, 0.432, 0.576],
        [0.,    0.072, 0.,    0.144, 0.,    0.072, 0.   ]],
        11,
        [[1.0000e+00, 2.0000e+00, 1.9314e+11, 5.2500e-06],
        [2.0000e+00, 4.0000e+00, 1.9314e+11, 5.2500e-06],
        [4.0000e+00, 6.0000e+00, 1.9314e+11, 5.2500e-06],
        [6.0000e+00, 7.0000e+00, 1.9314e+11, 5.2500e-06],
        [2.0000e+00, 3.0000e+00, 1.9314e+11, 5.2500e-06],
        [3.0000e+00, 4.0000e+00, 1.9314e+11, 5.2500e-06],
        [4.0000e+00, 5.0000e+00, 1.9314e+11, 5.2500e-06],
        [5.0000e+00, 6.0000e+00, 1.9314e+11, 5.2500e-06],
        [1.0000e+00, 3.0000e+00, 1.9314e+11, 5.2500e-06],
        [3.0000e+00, 5.0000e+00, 1.9314e+11, 5.2500e-06],
        [5.0000e+00, 7.0000e+00, 1.9314e+11, 5.2500e-06]],
        6,
        [[    0.],
        [    0.],
        [-1300.],
        [-1500.],
        [    0.],
        [    0.],
        [-1300.],
        [-1500.],
        [    0.],
        [    0.],
        [-1300.],
        [-1500.],
        [    0.],
        [    0.]],
        3,
        [[ 0.],
        [ 1.],
        [13.]]
    )

    assert assertion_object == solver_data
