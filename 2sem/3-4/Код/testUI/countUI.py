import sys
import ctypes
import webbrowser
import os.path
import pandas as pd
from PyQt5.QtWidgets import (QWidget, QProgressBar, QTabWidget, QPushButton, QHBoxLayout, QVBoxLayout, QApplication, QLineEdit, QLabel, QMainWindow, QFileDialog) 
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

import integral
from NMclasses import interpolation
import numpy as np

from numpy import cos

class MplCanvas(FigureCanvas):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        self.axes.grid(True)
        super(MplCanvas, self).__init__(fig)


class Tabs(QTabWidget):
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        self.tab1 = QWidget()
        self.tab2 = QWidget()
        self.addTab(self.tab1,"Tab 1")
        self.tab1UI()
        # self.show()


    def tab1UI(self):
        self.countButton = QPushButton("Start count")
        self.countButton.clicked.connect(self.plot)
        self.folderButton = QPushButton('Open folder')
        self.folderButton.clicked.connect(self.openFolder)
            
        self.folderLabel = QLabel('Папка с графиком')
        self.paramA = QLineEdit()
        self.paramB = QLineEdit()
        vbox = QVBoxLayout()
        hbox1 = QHBoxLayout()
        hbox2 = QHBoxLayout()
        hbox1.addWidget(QLabel('a'))
        hbox1.addWidget(self.paramA) 
        hbox2.addWidget(QLabel('b'))
        hbox2.addWidget(self.paramB)
        vbox.addLayout(hbox1)
        vbox.addLayout(hbox2)
        vbox.addWidget(self.countButton)
        vbox.addWidget(QProgressBar(self))

        self.tab1.setLayout(vbox)
        
    def tab2UI(self):
        self.canvas = MplCanvas(self,  width=5, height=4, dpi=100)
        toolbar = NavigationToolbar(self.canvas, self)
        vbox = QVBoxLayout()
        vbox.addWidget(toolbar)
        vbox.addWidget(self.canvas)
        vbox.addWidget(self.folderButton)
        self.tab2.setLayout(vbox)
            

        # widget = QWidget()
        # widget.setLayout(vbox)
        # self.setCentralWidget(widget)



    def openFolder(self):
        # os.mkdir(r"C:\Users\mchav\OneDrive\Учеба\Предметы\ЧМ\2 сем\3-4\Код\testUI\plots")
        webbrowser.open(r"C:\Users\mchav\OneDrive\Учеба\Предметы\ЧМ\2 сем\3-4\Код\testUI\plots")


    def plot(self):

        self.addTab(self.tab2,"Tab 2")
        self.tab2UI()

        f = lambda x: cos(x)
        x, y, splines = interpolation(f, int(self.paramA.text()), int(self.paramB.text()), 30, 'Normal').splineInterpolation()
        data = pd.DataFrame(x, y)
        data.to_csv(r"C:\Users\mchav\OneDrive\Учеба\Предметы\ЧМ\2 сем\3-4\Код\testUI\plots\result.csv", sep = ';')

        self.canvas.axes.cla()
        self.canvas.axes.grid(True)
        self.canvas.axes.plot(x,y)
        # self.canvas.axes.savefig(r"C:\Users\mchav\OneDrive\Учеба\Предметы\ЧМ\2 сем\3-4\Код\testUI\plots\image.jpg")
        self.canvas.draw()
    
    def callDllFunction(self, buf):
        dllLibrary = ctypes.cdll.LoadLibrary(r'C:\Users\mchav\OneDrive\Учеба\Предметы\ЧМ\2 сем\3-4\Код\testUI\calc.dll') 
        c_length = ctypes.c_size_t(len(buf))
        temp_ptr = ctypes.c_void_p(buf.ctypes.data)
        returnCode = dllLibrary.calc(temp_ptr, c_length)
        
        return returnCode


class Window(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setFixedSize(400, 400)
        self.widget = Tabs()
        self.setCentralWidget(self.widget)
        self.show() 


if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = Window()
    sys.exit(app.exec())