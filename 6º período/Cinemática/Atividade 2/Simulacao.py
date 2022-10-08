
#Código para solução do problema do mecanismo de 4 barras
#Autor: Kleber Junior

import numpy as np
import math as m
from scipy import linalg

#Entrada:
L1 = 300
L2 = 120
L3 = 180
L4 = 240
theta2 = 135
theta_2 = theta2
theta = np.linspace(0, 360, 361)

#Metodo de Análise trigonométrica

z = m.sqrt((L1**2) + (L2**2) - (2*L1*L2*m.cos(theta2*m.pi/180.0)))
gamma = (180*m.acos(((z**2) - (L3**2) - (L4**2))/(-2*L3*L4))/m.pi)
alpha = (180*m.acos(((z**2) + (L4**2) - (L3**2))/(2*z*L4))/m.pi)
beta = (180*m.acos(((z**2) + (L1**2) - (L2**2))/(2*z*L1))/m.pi)
theta4 = (180-(alpha+beta))
theta3 = (180-(alpha+beta+gamma))

print('Método de Análise Trigonométrica')
print("Theta 2 =",theta2)
print("Theta 3 =",theta3)
print("Theta 4 =",theta4)

#Metodo de Laços de vetores

K1 = L1/L2
K2 = L1/L4
K3 = ((L2**2) - (L3**2) + (L4**2) + (L1**2)) / (2*L2*L4)

A = np.cos(np.deg2rad(theta_2)) - K1 - (K2*np.cos(np.deg2rad(theta_2))) + K3
B = -2*np.sin(np.deg2rad(theta_2))
C = K1 - ((K2 + 1)*np.cos(np.deg2rad(theta_2))) + K3

theta_4v = 2*np.rad2deg(np.arctan((-B-np.sqrt(B**2 - 4*A*C))/(2*A)))
theta_3v = np.rad2deg(np.arcsin(((L4*np.sin(np.deg2rad(theta_4v))) - (L2*np.sin(np.deg2rad(theta_2))))/L3))

print('')
print('Método de Laço de Vetores')
print("Theta 3 =", theta_3v)
print("Theta 4 =", theta_4v)


#Metódo de Newton Raphson

theta_3 = 8   
theta_4 = 100
iteracoes_de_metodo = 5

for i in range(iteracoes_de_metodo):
    A_A = np.array([[-L3*np.sin(np.deg2rad(theta_3)), L4*np.sin(np.deg2rad(theta_4))],
                   [L3*np.cos(np.deg2rad(theta_3)), -L4*np.cos(np.deg2rad(theta_4))]])
    B_B = np.array([[-L2*np.cos(np.deg2rad(theta_2)) - L3*np.cos(np.deg2rad(theta_3)) + L4*np.cos(np.deg2rad(theta_4)) + L1],
                   [-L2*np.sin(np.deg2rad(theta_2)) - L3*np.sin(np.deg2rad(theta_3)) + L4*np.sin(np.deg2rad(theta_4))]])
    X = linalg.solve(A_A, B_B)
    AX = B
    theta_3 = theta_3 + np.rad2deg(float(X[0]))
    theta_4 = theta_4 + np.rad2deg(float(X[1]))

print('')
print('Método de Newton Raphson')
print("Theta 3 =", theta_3)
print("Theta 4 =", theta_4)


#Método Algébrico

a = L2 * np.ones(theta.shape)
b = L3 * np.ones(theta.shape)
c = L4 * np.ones(theta.shape)
d = L1 * np.ones(theta.shape)

A_x = a * np.cos(np.deg2rad(theta))
A_y = a * np.sin(np.deg2rad(theta))
S = ((a**2) - (b**2) + (c**2) - (d**2)) / (2*(A_x - d))
P = ((A_y**2)/(A_x - d)**2) + 1
Q = (2*A_y * (d - S)) / (A_x - d)
R = (d - S)**2 - (c**2)
B_y = (-Q + np.sqrt((Q**2) - 4*P*R)) / (2*P)
B_x = S - ((2*A_y*B_y) / (2*(A_x - d)))
B_y1 = (-Q - np.sqrt((Q**2) - 4*P*R)) / (2*P)
B_x1 = S - ((2*A_y*B_y1) / (2*(A_x - d)))
theta_3 = np.rad2deg(np.arctan((B_y[theta2]-A_y[theta2])/(B_x[theta2]-A_x[theta2])))
theta_4 = np.rad2deg(np.pi + np.arctan((B_y[theta2])/(B_x[theta2]-L1)))

print('')
print('Método Algébrico')
print("Theta 3 =", theta_3)
print("Theta 4 =", theta_4)