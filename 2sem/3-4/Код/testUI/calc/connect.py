#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ctypes
import numpy as np

def callDllFunction(buf):
    ## Solution for Windows
    dllLibrary = ctypes.cdll.LoadLibrary(r'C:\Users\mchav\OneDrive\Учеба\Предметы\ЧМ\2 сем\3-4\Код\testUI\calc\calc.dll')
    
    ## Solution for macOS
    # dllLibrary = ctypes.cdll.LoadLibrary('main.dylib')
    
    ## Solution for Linux
    # dllLibrary = ctypes.cdll.LoadLibrary('main.so')
    
    c_length = ctypes.c_size_t(len(buf))
    temp_ptr = ctypes.c_void_p(buf.ctypes.data)
    returnCode = dllLibrary.calc(temp_ptr, c_length)
    
    return returnCode

if __name__ == '__main__':
    temp = np.empty(100, dtype=np.float64)
    temp.fill(20.0)
    print(temp)
    callDllFunction(temp)
    print(temp)
