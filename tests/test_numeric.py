import numpy as np
from src.numeric_method import gauss_method  

def test_gauss_solution():
    '''
    Test with known matrix and force vector
    '''
    A = 1e8 * np.array([[1.59, -0.4, -0.54], [-0.4, 1.7, 0.4], [-0.54, 0.4, 0.54]])
    B = np.array([0, 150, -100])
    max_iterations = 1000
    tol = 1e-6

    expected_solution = np.linalg.solve(A, B)
    with np.errstate(divide='ignore', invalid='ignore'):
        X, error = gauss_method(max_iterations, tol, A, B)

    assert np.allclose(X, expected_solution), "Calculated solution is not close to expected."
    assert error < tol, "The final error is greater than the tolerance."

    assert not np.isnan(X).any(), "The solution contains NaN values."
    assert not np.isinf(X).any(), "The solution contains infinite values."

def test_convergence():
    '''
    Verify that the method does not converge within the maximum number of iterations
    '''
    A = np.array([[2, 1], [5, 7]])
    B = np.array([11, 13])
    max_iterations = 1 
    tol = 1e-6

    X, error = gauss_method(max_iterations, tol, A, B)

    assert error > tol, "Unexpectedly small error for limited number of iterations."

