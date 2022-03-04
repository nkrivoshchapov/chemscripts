from math import sin, cos
import numpy as np


DEG2RAD = 3.14159265358979/180
def rotation_matrix_2axis(ang):
    ang *= DEG2RAD
    return np.array([[cos(ang), 0, sin(ang)],
                     [0, 1, 0],
                     [-sin(ang), 0, cos(ang)]])

def rotation_matrix_1axis(ang):
    ang *= DEG2RAD
    return np.array([[1, 0, 0],
                     [0, cos(ang), -sin(ang)],
                     [0, sin(ang),  cos(ang)]])

def rotation_matrix_3axis(ang):
    ang *= DEG2RAD
    return np.array([[cos(ang), -sin(ang), 0],
                     [sin(ang), cos(ang), 0],
                     [0,0,1]])
