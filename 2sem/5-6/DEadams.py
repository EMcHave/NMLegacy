import math
import numpy as np
from NMclasses import DE
from NMclasses import Plots
import matplotlib.pyplot as plt
import DErunge

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
        n, d, k = 4, 1, 0
        DSp = DE(d2y, a, b, n, y0, dy0).RungeIV()
        while d > eps:
            D = []
            DSn = DE(d2y, a, b, 2*n, y0, dy0).AdamsIV()
            for j in range(1, DSp.shape[1]):
                d = math.fabs(DSn[1][2*j] - DSp[1][j])/15
                D.append(d)
            DSp = DSn
            n *= 2
            d = max(D)
            k+=1
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

    EpsA, EpsRealA, NA = RungeEps(yR, 13)
    EpsR, EpsRealR, NR = DErunge.RungeEps(yR, 13)
    Plots(1, [[EpsA, NA], [EpsA, NR]], 'Количество итераций от заданной точности', 'Eps', 'K',['Адамс', 'РК4']).build('logX')
    Plots(2, [[NA, EpsRealA], [NA, EpsA], [NA, EpsRealR], [NA, EpsR]], 'Погрешность от итерации', 'K', 'Eps', ['Реальная погрешность Адамс', 'По правилу Рунге Адамс', 'Реальная погрешность РК4' ,' По правилу Рунге РК4']).build('logY')
    Plots(3, [[EpsA, EpsRealA], [EpsA, EpsRealR]], 'Абсолютная погршность от заданной точности', 'Eps', 'EpsReal', ['Адамс', 'РК4']).build('loglog')
    plt.show()