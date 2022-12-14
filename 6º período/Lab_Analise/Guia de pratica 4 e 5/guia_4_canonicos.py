# -*- coding: utf-8 -*-
"""
Criado em 19/10/2022
Autor: Kleber Junior e Robson Junior
 
"""
<<<<<<< Updated upstream
import numpy as np
import control as ct
<<<<<<< Updated upstream
=======

plt.close('all')

t0 = 0                           # Tempo inicial
tf = 1000                        # Tempo final
t = np.linspace(t0, tf, 1999)    # Instantes que desejo ter a solucao

>>>>>>> Stashed changes
=======

import numpy as np #Importa a biblioteca numpy
import control as ct #Importa a biblioteca control 
>>>>>>> Stashed changes

#criando a função de tranferencia e printando variáveis do sistema
print('------------------------------------------------------------------------')
print('\n Função de transferência \n') 

s = ct.tf('s')                         #Cria uma variável 's' para permitir operações de álgebra para sistemas SISO
G = (s + 3)/(s**2 + 6*s + 5)           #Definindo a função de transferencia
G_s = ct.tf2ss(G) # convertendo a função de tranferencia para espaço de estados
print(G) #printa a função de transferencia

print('------------------------------------------------------------------------')
print('\n Representação em espaço de estados na forma Controlável\n')
[A,B] = ct.canonical_form(G_s, form='reachable')       # Forma Controlável
print(A)

print('------------------------------------------------------------------------')
print('\n Representação em espaço de estados na forma Observável\n')
[A,B] = ct.canonical_form(G_s, form='observable')       # Forma observavel
print(A)

print('------------------------------------------------------------------------')
print('\n Representação em espaço de estados na forma Diagonal\n')
[A,B] = ct.canonical_form(G_s, form='modal')       # Forma Diagonal
print(A)
print('------------------------------------------------------------------------')