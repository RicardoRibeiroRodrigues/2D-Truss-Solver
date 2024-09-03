"""
Module to solve a linear system using the Gauss method.

Exports:
    gauss_method: function to solve a linear system using the Gauss method.
"""

import numpy as np


def gauss_method(max_iterations, tol, k_matrix, force_vector):
    """
    Gauss method to solve a linear system.

    Parameters
        max_iterations: maximum number of iterations
        tol: tolerance.
        k_matrix: rigidity matrix
        force_vector: vector of forces
    Returns
        x_m: vector of displacements
        err: error
    """
    error = 1e9
    ite = 0
    output_size = k_matrix.shape[0]
    x_m = np.zeros(output_size)
    x_m_ant = x_m.copy()
    while error > tol and ite <= max_iterations:
        for i in range(output_size):
            coef = k_matrix[i]
            x_m[i] = force_vector[i] / coef[i]
            for j in range(output_size):
                if j != i:
                    x_m[i] -= (x_m[j] * coef[j]) / coef[i]
        error = np.nanmax(
            abs((x_m[0:output_size] - x_m_ant[0:output_size]) / x_m[0:output_size])
        )
        ite += 1
        x_m_ant = x_m.copy()
    return x_m, error


if __name__ == "__main__":
    # Read the input file
    A = 1e8 * np.array([[1.59, -0.4, -0.54], [-0.4, 1.7, 0.4], [-0.54, 0.4, 0.54]])
    B = np.array([0, 150, -100])
    X, e = gauss_method(1e3, 1e-6, A, B)
    assert np.allclose(np.linalg.solve(A, B), X)
