# -*- coding: utf-8 -*-
"""
Criado em 19/10/2022
Autor: Kleber Junior e Robson Junior
 
"""
import numpy as np
import control as ct

#criando a função de tranferencia
print('------------------------------------------------------------------------')
print('\n Função de transferência \n')

s = ct.tf('s')                         #Cria uma variável 's' para permitir operações de álgebra para sistemas SISO
G = (s + 3)/(s**2 + 6*s + 5)           #Definindo a função de transferencia
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