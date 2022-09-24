from scipy.integrate import solve_ivp
from scipy.integrate import quad
from scipy.misc import derivative
import matplotlib.pyplot as plt
import numpy as np
import math as m

plt.close('all') # Fecha gráficos

#Derivada da função desejada:
def dhdt(t, h, A, beta, Fi):
    return ((Fi/A) - ((beta/A)*m.sqrt(h)))

#Erro:
def erro(t, setpoint, h):
    return (setpoint- h)

#Constantes para o controlador PI:
Kp = 1
Ki = 1

#Parâmetros que vem do próprio sistema:
area = 0.2
beta = 0.1

#Condições iniciais:
u = 0
h = 1

#Parâmetros de simulação:
t = 0.00
tf = 5.00
step = 0.01

#Vetores

h_eixo = []
u_eixo = []
t_eixo = []
integral = 0

#Altura desejada (ponto de equilíbrio):
setpoint = 0.5

while(t < tf):
    #Controle:
    u_eixo.append(u)
    h_eixo.append(h)
    t_eixo.append(t)
    #Proporcional:
    prop = erro(None, setpoint, h)
    #Integral
    integral = (prop*step) + integral
    if(prop==0):
        integral=0
    #Calculo do controle PI
    u = (Kp*prop)+(Ki*integral)
    if(u<0):
        u=0
    #Solução da EDO:
    sol = solve_ivp(dhdt, t_span=(t, t+step), y0=[h], method='RK45', t_eval=[t, t+step], args=(area, beta, u))
 
    h = sol.y[0][-1] + np.random.normal(0, 0.002)
    #Tempo instantânio
    t += step


#Plotagem do gráfico:
plt.figure(1)
plt.subplot(2,1,1)
#plt.ylim([0, 1])
plt.plot(t_eixo, h_eixo, 'b',label='Nível, h(t)')
plt.ylabel('$h(t)$', fontsize=12)
plt.xlabel('$t$', fontsize=12)
plt.grid(axis='both')
plt.legend(fontsize=10)
plt.subplot(2,1,2)
#plt.ylim([-0.02, 0.02])
plt.plot(t_eixo, u_eixo, 'r',label='Controle, u(t)')
plt.ylabel('$u(t)$', fontsize=12)
plt.xlabel('$t$', fontsize=12)
plt.grid(axis='both')
plt.legend(fontsize=10)
plt.show()