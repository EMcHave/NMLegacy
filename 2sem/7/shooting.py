import math
import random
import statistics
import numpy as np
from NMclasses import DE
from NMclasses import Plots
import matplotlib.pyplot as plt
import pandas as pd

d2y = lambda x, y, dy: 2*x*dy +2*y - 4*x
yR = np.vectorize(lambda x: x + np.exp(x**2))

y0 = 1
dy0 = 1
a, b = 0, 1
xR = np.linspace(a, b, 100)

def tableMaker():
    de1 = DE(d2y, a, b, 20, y0, dy0)
    de2 = DE(d2y, a, b, 10, y0, dy0)
    sol10 = de1.Shooter(yR)
    sol5 = de2.Shooter(yR)
    result = {}
    Eps, Epsr = [], []
    for i in range(sol5.shape[1]):
        E = math.fabs(sol5[1][i]-sol10[1][2*i])/15
        Er = math.fabs(yR(sol10[0][2*i])-sol10[1][2*i])
        result[f"{sol10[0][2*i]}"] = [E, Er]
        Eps.append(E)
        Epsr.append(Er)
    table = pd.DataFrame(result, index = ['Рунге','Абсолютная'])
    table.to_excel(r"C:\Users\mchav\OneDrive\Учеба\Предметы\ЧМ\2 сем\7\Отчет\table.xlsx")

    return Eps, Epsr

def disturbance():
    result = {}
    de = DE(d2y, a, b, 5, y0, dy0)
    undisturbedSolution = de.Shooter(yR)
    for i in range(5):    
        r = random.random()
        de.y0 += r
        disturbedSolution = de.Shooter(yR)
        d = []
        for j in range(disturbedSolution.shape[1]): 
            d.append( math.fabs(disturbedSolution[1][j] - undisturbedSolution[1][j]) )
        de.y0 -=r
        result[f"{r}"] = d
    table = pd.DataFrame(result, index = ['Ошибка в точке 0', 'Ошибка в точке 1', 'Ошибка в точке 2', 'Ошибка в точке 3', 'Ошибка в точке 4', 'Ошибка в точке 5'])
    table.to_excel(r"C:\Users\mchav\OneDrive\Учеба\Предметы\ЧМ\2 сем\7\Отчет\tableDisturbance.xlsx")

if __name__ == '__main__':
    Eps, Epsr = tableMaker()
    disturbance()
    de = DE(d2y, a, b, 10, y0, dy0)
    sol = de.Shooter(yR)
     
    testSol = DE(d2y, a, b, 10, y0, 1.08).RungeIV()
    Plots(1, [[sol[0], Eps], [sol[0], Epsr]], 'Погрешность в узлах', 'x', 'Погрешность', ['По Рунге', 'Абсолютная']).build()
    plt.show()