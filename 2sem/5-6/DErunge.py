import math
import numpy as np
from NMclasses import DE
from NMclasses import Plots
import matplotlib.pyplot as plt

d2y = lambda x, y, dy: 2*x*dy +2*y - 4*x
yR = np.vectorize(lambda x: x + np.exp(x**2))

y0 = 1
dy0 = 1
a, b = 0, 1
xR = np.linspace(a, b, 100)


def RungeEps(y, epsMax):   
    eps = 0.1

    Eps = list()
    EpsReal = list()
    N = list()

    for i in range(1, epsMax):
        eps = math.pow(10, -i)
        n, d, k = 2, 1, 0
        DSp = DE(d2y, a, b, n, y0, dy0).RungeIV()
        while d > eps:
            D = []
            DSn = DE(d2y, a, b, 2*n, y0, dy0).RungeIV()
            for j in range(1, DSp.shape[1]):
                d = math.fabs(DSn[1][2*j] - DSp[1][j])/15
                D.append(d)
            DSp = DSn
            n *= 2
            d = max(D)
            k+=1
        # K[f'10E-{i}'] = n
        EpsReal.append(realDifference(y, DSn))
        N.append(k)       
        Eps.append(eps)
    return Eps, EpsReal, N

def realDifference(yReal, numericSolution):
    EpsTemp = list()
    for index in range(len(numericSolution)):
        EpsTemp.append(math.fabs(yReal(numericSolution[0][index]) - numericSolution[1][index]))
    EpsReal = max(EpsTemp)
    return EpsReal

if __name__ == '__main__':
    pass
    # Eps, EpsReal, N = RungeEps(yR, 13)
    # print(EpsReal)
    # Plots(1, [[Eps, N]], 'Количество итераций от заданной точности', 'Eps', 'K').build('logX')
    # Plots(2, [[N, EpsReal], [N, Eps]], 'Погрешность от итерации', 'K', 'Eps', ['Реальная погрешность', 'По правилу Рунге']).build('logY')
    # Plots(3, [[Eps, EpsReal]], 'Абсолютная погршность от заданной точности', 'Eps', 'EpsReal').build('loglog')
    # plt.show()