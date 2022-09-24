from scipy.integrate import solve_ivp
from scipy.integrate import quad
from scipy.misc import derivative
import matplotlib.pyplot as plt
import numpy as np
import math as m

#Derivada da função desejada:
def dhdt(t, h, A, beta, Fi):
    if(h < 0):
        h = 0
    return ((Fi/A) - ((beta/A)*m.sqrt(h)))

#Erro:
def erro(t, x_d, x=None):
    if(x is None):
        x = x_d
        x_d = t
    return (x_d-x)

#Constantes do controlador PI:
Kp = 50
Ki = 0.1
kd = 1

#Parâmetros do sistema:
area = 1
beta = .1

#Condições iniciais:
h0 = 1
u = 0

#Parâmetros de simulação:
t = 0
tf = 50
step = 0.01

#Looping:
h = h0
h_axis = []
h_axis.append(h)
u_axis = []
u_axis.append(u)
t_axis = []
t_axis.append(t)

#Altura desejada (ponto de equilíbrio):
h_D = 0.8

while(t < tf):
    #Controle:
    #Proporcional:
    prop = erro(h_D, h)
    #Integrativo:
    integration, er = quad(erro, t, t+step, args=(h_D, h))
    #Derivativo
    derivation = derivative(erro, h, dx=1e-6, args=(h_D, h))
    u = ((Kp*prop) + (Ki*integration) + (kd*derivation))
    if(u < 0):
        u = 0
    u_axis.append(u)

    #Solução da EDO:
    sol = solve_ivp(dhdt, t_span=(t, t+step), y0=[h], method='RK23', t_eval=[t, t+step], args=(area, beta, u))

    h = sol.y[0][-1] #+ np.random.normal(0, 0.002)
    h_axis.append(h)
    t += step
    t_axis.append(t)
    

#Plotagem do gráfico:
plt.figure(1)
plt.subplot(2,1,1)
#plt.ylim([0.74, 0.76])
plt.plot(t_axis, h_axis, 'k',label='Nível, h(t)')
plt.ylabel('$h(t)$', fontsize=12)
plt.xlabel('$t$', fontsize=12)
plt.grid(axis='both')
plt.legend(fontsize=10)
plt.subplot(2,1,2)
#plt.ylim([-0.02, 0.02])
plt.plot(t_axis, u_axis, 'b',label='Controle, u(t)')
plt.ylabel('$u(t)$', fontsize=12)
plt.xlabel('$t$', fontsize=12)
plt.grid(axis='both')
plt.legend(fontsize=10)
plt.show()