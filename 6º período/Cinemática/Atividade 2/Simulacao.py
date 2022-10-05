#Código para solução do problema do mecanismo de 4 barras
#Autor: Kleber Junior
import numpy as np
import math as m
import matplotlib.pyplot as plt

#Entrada:
L1 = 300
L2 = 120
L3 = 180
L4 = 240
theta2 = 1
theta = np.linspace(0, 360, 361)

#suposição para valores iniciais para o método de Newton Raphson
theta_3 = 8   #Método Geométrico =  15
theta_4 = 100    #Método Geométrico = 140
iteracoes_de_metodo = 500   #Número de iteracoes

z = m.sqrt((L1**2) + (L2**2) - (2*L1*L2*m.cos(theta2*m.pi/180.0)))
gamma = (180*m.acos(((z**2) - (L3**2) - (L4**2))/(-2*L3*L4))/m.pi)
alpha = (180*m.acos(((z**2) + (L4**2) - (L3**2))/(2*z*L4))/m.pi)
beta = (180*m.acos(((z**2) + (L1**2) - (L2**2))/(2*z*L1))/m.pi)
theta4 = (180-(alpha+beta))

#criar vetores de mesma dimensão do vetor de theta2 e com os valores dos links
a = L2 * np.ones(theta.shape)
b = L3 * np.ones(theta.shape)
c = L4 * np.ones(theta.shape)
d = L1 * np.ones(theta.shape)

