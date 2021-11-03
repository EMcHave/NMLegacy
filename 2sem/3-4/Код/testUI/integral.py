import numpy as np 
import NMclasses as NM
import math 
import matplotlib.pyplot as plt
import pandas as pd
import time

from PyQt5 import *



f = lambda x: np.cos(x)
a = -3
b = 7
x = np.linspace(a, b, 500)

def integralEps(epsMax):   
    eps = 0.1

    K = dict()
    Eps = list()
    N = list()

    for i in range(1, epsMax):
        eps = math.pow(10, -i)
        n, d = 1, 1
        t_in = time.time()
        Ip = NM.Integtration(f, a, b, n).lobattoInt()
        while d > eps:
            In = NM.Integtration(f, a, b, 2*n).lobattoInt()
            d = math.fabs(In - Ip)/63
            Ip = In
            n *= 2
        t_out = time.time()
        ar = [n, t_out - t_in]
        K[f'10E-{i}'] = ar
        N.append(n)       
        Eps.append(eps)
    
    return K, Eps, N

if __name__ == '__main__':

    print(NM.Integtration(f, a, b, 2).lobattoInt())
    K, Eps, N = integralEps(16)

    table = pd.DataFrame(K, index = ['n, точек разб', 't, сек'])
    table.to_excel(r'C:\Users\mchav\OneDrive\Учеба\Предметы\ЧМ\2 сем\3-4\Код\table.xlsx')
    
    NM.Plots(1, [[x, f(x)]], 'График функции', 'x', 'f(x)')

    plt.show()
    # NM.Plots(2, [[Eps, N]], 'Зависимость кол-ва интервалов от точности', 'Eps', 'N')