import numpy as np # importando biblioteca numpy
import matplotlib.pyplot as plt # importando biblioteca para plotar as figuras
import control as ct  #importanto biblioteca control

plt.close('all') #comando para fechar todas janelas de plot

t0 = 0    # tempo inicial
tf = 15   # tempo final
t = np.linspace(t0,tf,10001) # instantes que desejo ter a solucao

K = 1.9187*10**-5
N = 2.4451
C = 51.161

T = t[1]-t[0] #período de amostragem

u= 0* np.ones(t.shape) #Sinal de controle necessário para levar o sistema para o 
#ponto de operação desejado em malha aberta.

