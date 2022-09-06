
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
Tf = 100 #segundos
step = 0.01 #segundos


#Eixos de plotagem
I_eixo = []
V_eixo = []
T_eixo = []

while(T<=Tf):
    #Preenche os Vetores
    I_eixo.append(I)
    V_eixo.append(V)
    T_eixo.append(T)

    #Armazena o valor Numérico 
    I = ((m.exp(-0.125 * T)*m.cos(0.3307))+(1.1339*(m.exp(-0.125*T)*m.sin(0.3307))))
    #Tempo instantânio
    T += (Ti + step)

plt.subplot(1,1,1)
plt.plot(T_eixo, I_eixo, 'r',label='Corrente, I(t)')
plt.ylabel('$u(t)$', fontsize=12)
plt.xlabel('$t$', fontsize=12)
plt.grid(axis='both')
plt.legend(fontsize=10)
plt.show()