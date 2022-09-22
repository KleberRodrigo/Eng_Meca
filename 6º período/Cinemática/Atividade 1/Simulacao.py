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
theta2 = 0
eixo_z=[]
eixo_gamma=[]
eixo_alpha=[]
eixo_beta=[]
eixo_theta4=[]
eixo_theta2=[]

while(theta2<720):
    z = m.sqrt((L1**2) + (L2**2) - (2*L1*L2*m.cos(theta2*m.pi/180.0)))
    gamma = (180*m.acos(((z**2) - (L3**2) - (L4**2))/(-2*L3*L4))/m.pi)
    alpha = (180*m.acos(((z**2) + (L4**2) - (L3**2))/(2*z*L4))/m.pi)
    beta = (180*m.acos(((z**2) + (L1**2) - (L2**2))/(2*z*L1))/m.pi)
    theta4 = (180-(alpha+beta))
    eixo_z.append(z)
    eixo_gamma.append(gamma)
    eixo_alpha.append(alpha)
    eixo_beta.append(beta)
    eixo_theta4.append(theta4)
    eixo_theta2.append(theta2)
    theta2 = theta2 + 0.1
    
plt.figure(1)
plt.plot(eixo_theta2,eixo_z ,'b', label='z')
plt.plot(eixo_theta2,eixo_gamma ,'r', label='gamma')
plt.plot(eixo_theta2,eixo_alpha ,'y', label='alpha')
plt.plot(eixo_theta2,eixo_theta4 ,'g', label='theta4')

plt.ylabel('$Ângulo$', fontsize=12)
plt.xlabel('$Theta2$', fontsize=12)
plt.grid(axis='both')
plt.legend(fontsize=10)
plt.title('Ângulos em função de theta 2')
plt.show()
