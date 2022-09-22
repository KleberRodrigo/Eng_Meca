#Código para solução do problema de Secagem de grão
#Autor: Kleber Junior E Robson Resende

import numpy as np
import math as m
import matplotlib.pyplot as plt

t0=0
tf=100
eixo_t=[]
eixo_M=[]

temperatura=32
Rh=.19
mh=28
#k=(0.1579+(0.0001746*temperatura)-(0.0001413*Rh))
#n=(0.6545+(0.002425*temperatura)+(0.0007867*Rh))

k=(0.2958+(0.01215*temperatura)-(0.44565*Rh))
n=(0.13365+(0.009468*temperatura)+(0.0193653*Rh)-(0.000177431*(Rh**2)))

while(t0<tf):
    eixo_t.append(t0)
    eixo_M.append(mh)
    mh= (18.69*(m.exp(-0.42*(t0**0.567))))+9.31
    t0=t0+0.1

plt.figure(1)
plt.plot(eixo_t,eixo_M,'b',label='M')
plt.ylabel('$Humidade do grão$', fontsize=12)
plt.xlabel('$Tempo$', fontsize=12)
plt.grid(axis='both')
plt.legend(fontsize=10)
plt.title('Simulação da umidade em relação ao tempo')
plt.show()
