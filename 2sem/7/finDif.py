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

x, y= DE(d2y, a, b, 10, y0, dy0).finiteDifferences(lambda x: -2*x, lambda x: -2*(x**0), lambda x: -4*x, yR)

plt.plot(x, y)
plt.plot(xR, yR(xR))
plt.show()