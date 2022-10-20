# -*- coding: utf-8 -*-
"""
Criado em 19/10/2022
Autor: Kleber Junior e Robson Junior
 
"""

import numpy as np
from matplotlib import pyplot as plt
import control as ct
import math as m

plt.close('all')

t0 = 0                           # Tempo inicial
tf = 1000                        # Tempo final
t = np.linspace(t0, tf, 1999)    # Instantes que desejo ter a solucao


#criando a função de tranferencia

print('------------------------------------------------------------------------')
print('\n Função de transferência \n')

s = ct.tf('s')                         #Cria uma variável 's' para permitir operações de álgebra para sistemas SISO
#G = (s + 3)/(s**2 + 6*s + 5)           #Definindo a função de transferencia
num=[1,3]
den=[1,6,5]
G=ct.TransferFunction(num, den)
G_s = ct.tf2ss(G)
print(G)

print('------------------------------------------------------------------------')
print('\n Representação em espaço de estados na forma Controlável\n')
[A,B] = ct.canonical_form(G_s, form='reachable')       # Controlável
print(A)

print('------------------------------------------------------------------------')
print('\n Representação em espaço de estados na forma Observável\n')
[A,B] = ct.canonical_form(G_s, form='observable')       # observavel
print(A)

print('------------------------------------------------------------------------')
print('\n Representação em espaço de estados na forma Diagonal\n')
[A,B] = ct.canonical_form(G_s, form='modal')       # Diagonal
print(A)
print('------------------------------------------------------------------------')

