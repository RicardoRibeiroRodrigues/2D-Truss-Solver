import numpy as np


def met_gauss(n_ite, tol, K, F):
    """
    n_ite: número máximo de iterações
    tol: tolerância
    K: matriz de rigidez
    F: vetor de forças
    """
    err = 1e9
    ite = 0
    n = K.shape[0]
    x_m = np.zeros(n)
    x_m_ant = x_m.copy()
    while err > tol and ite <= n_ite:
        for i in range(n):
            coef = K[i]
            x_m[i] = F[i] / coef[i]
            j = 0
            for j in range(n):
                if j != i:
                    x_m[i] -= (x_m[j] * coef[j]) / coef[i]
        err = np.nanmax(
            abs((x_m[0:n] - x_m_ant[0:n]) / x_m[0:n])
        )
        ite += 1
        x_m_ant = x_m.copy()
    print("Número de iterações: ", ite)
    return x_m, err


if __name__ ==  '__main__':
    # Read the input file
    A = 1e8 * np.array([[1.59, -0.4, -0.54], [-0.4, 1.7, 0.4], [-0.54, .4, .54]])
    B = np.array([0, 150, -100])
    X, e = met_gauss(1e3, 1e-6, A, B)    
    print(X)
    print(e)