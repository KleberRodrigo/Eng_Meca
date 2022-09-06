
from cmath import exp
from re import T
from scipy.integrate import solve_ivp
from scipy.integrate import quad
from scipy.misc import derivative
import matplotlib.pyplot as plt
import numpy as np
import math as m

plt.close('all')
 
#Variáveis de simulação 
T  = 0 #tempo real 
Ti = 0 #segundos
Tf = 40 #segundos
step = 0.01 #segundos

#Variaveis do sistema 

I=0 # corrente inicial
V=0

#Eixos de plotagem
I_eixo = []
V_eixo = []
T_eixo = []

while(T<=Tf):
    #Calcula a corrente em função do tempo
    I = ((m.exp(-0.125 * T) *m.cos(0.3307*T))+(1.1339*(m.exp(-0.125*T)*m.sin(0.3307*T))))
    V = ( ((5*m.exp(-0.125*T))*((2.6459*m.sin(0.3307*T))-(m.cos(0.3307*T))))-5 )
    #Calcula a tensão em função do tempo
    #Preenche os Vetores
    I_eixo.append(I)
    V_eixo.append(V)
    #Tempo instantânio
    T_eixo.append(T)
    T += (Ti + step)

    #V_eixo.append(V)
    

plt.figure(1)
plt.subplot(2,1,1)
plt.plot(T_eixo, I_eixo, 'r',label='Corrente, I(t)')
plt.ylabel('$I(t)$', fontsize=12)
plt.xlabel('$t$', fontsize=12)
plt.grid(axis='both')
plt.legend(fontsize=10)

plt.subplot(2,1,2)
plt.plot(T_eixo, V_eixo, 'b',label='Tensão, I(t)')
plt.ylabel('$V(t)$', fontsize=12)
plt.xlabel('$t$', fontsize=12)
plt.grid(axis='both')
plt.legend(fontsize=10)
plt.show()