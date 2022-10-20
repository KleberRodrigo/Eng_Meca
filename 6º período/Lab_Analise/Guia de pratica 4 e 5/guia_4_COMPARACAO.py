""" Autor: Kleber Junior e RObson Junior 
    Criado em : 19/10/2022 """

import control as ct
import matplotlib.pyplot as plt
import numpy as np

t0=0
tf=10
time_step= np.linspace(t0, tf, 1001)

u = np.ones(time_step.shape)

x1_0 = 0
x2_0 = 0
vet_estado = [x1_0, x2_0]

#utilizando a forma controlável temos que 
# dx1/dt = x2(t)
# dx2/dt = -5x1(t) -6x2(t) + u(t)

def model_update(time_step, vet_estado, u, params):  
    dx1 = vet_estado[1]
    dx2 = -5*vet_estado[0] -6*vet_estado[1] + u 

    return [dx1, dx2]

def model_output(time_step, vet_estado, u, params):
    return vet_estado

SYSTEM = ct.NonlinearIOSystem(model_update, model_output,
    states=('x1', 'x2'), name='SYSTEM', inputs=('u'), outputs=('x1', 'x2'))

time_step, y = ct.input_output_response(SYSTEM, time_step, u, vet_estado)

est_1 = ((-1/4)*np.exp(- time_step)) + ((1/20)*np.exp(-5*time_step))+(1/5)
est_2 = ((1/4)*np.exp(- time_step)) - ((1/4)*np.exp(-5*time_step))

plt.figure(1)

plt.subplot(2,1,1)
plt.title('Comparação entre a forma computacional e a analítica', fontsize=10)
plt.plot(time_step, y[0], 'b', label='Espaço de estados', linewidth = 2)
plt.plot(time_step, est_1, 'r--', label='Matriz de transição de estados', linewidth = 2.5)
plt.ylabel('Estado de x1', fontsize=10)
plt.legend(fontsize=10)
plt.tick_params(labelsize=10)

plt.subplot(2,1,2)
plt.plot(time_step, y[1], 'b', label='Espaço de estados', linewidth = 2.0)
plt.plot(time_step, est_2, 'r--', label='Matriz de transição de estados', linewidth = 2.5)
plt.ylabel('Estado de x2', fontsize=10)
plt.xlabel('Tempo (s)', fontsize=10)
plt.legend(fontsize=10)
plt.tick_params(labelsize=10)

plt.show()