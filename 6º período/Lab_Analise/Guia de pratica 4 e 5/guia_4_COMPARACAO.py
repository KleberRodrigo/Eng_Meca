""" Autor: Kleber Junior e RObson Junior 
    Criado em : 19/10/2022 """

import control as ct #Importa a biblioteca control
import matplotlib.pyplot as plt #Importa a biblioeteca matplot
import numpy as np #Importa a biblioteca numpy 

t0=0 #Tempo inicial 
tf=10 #Tempo final de simulação
time_step= np.linspace(t0, tf, 1001) #Tempo de simulação, divide o tempo de simulação em 1000 espaços 

u = np.ones(time_step.shape) #Gera um vetor de valor unitario do tamanho do time_step

x1_0 = 0 #Variável de espaço 1
x2_0 = 0 #Variável de espaço 2
vet_estado = [x1_0, x2_0] # vetor de estados

#Define a equação diferencia nos estados
def model_update(time_step, vet_estado, u, params):
    dx1 = vet_estado[1] #Definição de x1
    dx2 = -5*vet_estado[0] -6*vet_estado[1] + u #Definição de x2
    return [dx1, dx2]

#definição da biblioteca
def model_output(time_step, vet_estado, u, params):
    return vet_estado #Retorna o vetor de estados

SYSTEM = ct.NonlinearIOSystem(model_update, model_output,
    states=('x1', 'x2'), name='SYSTEM', inputs=('u'), outputs=('x1', 'x2')) #Função para resolver a equação

time_step, y = ct.input_output_response(SYSTEM, time_step, u, vet_estado) #retorna o vetor de tempo e resposta em y
est_1 = ((-1/4)*np.exp(- time_step)) + ((1/20)*np.exp(-5*time_step))+(1/5) #Função para resolver a equação pelo formato de estados
est_2 = ((1/4)*np.exp(- time_step)) - ((1/4)*np.exp(-5*time_step)) #Função para resolver a equação pelo formato de equação de estados

plt.figure(1)#plotando figura
plt.subplot(2,1,1)#plotando subplot
plt.plot(time_step, y[0], 'b', label='Equação do espaço de estados', linewidth = 2) #Plotando o gráfico
plt.plot(time_step, est_1, 'r--', label='Resposta temporal',linewidth = 2.5)#Plotando a matriz de estados
plt.ylabel('Estado de x1', fontsize=10)#Legenda de x1
plt.legend(fontsize=10) #tamanho da legenda
plt.tick_params(labelsize=10)#Parametros do gráfico

plt.subplot(2,1,2) #plotando subfigura
plt.plot(time_step, y[1], 'b', label='Equação do espaço de estados', linewidth = 2.0) #Plotando e escrevendo o Titulo do plot
plt.plot(time_step, est_2, 'r--', label='Resposta temporal', linewidth = 2.5) #Plotando outro vetor da matriz
plt.ylabel('Estado de x2', fontsize=10) #Legenda do estado de x2
plt.xlabel('Tempo (s)', fontsize=10) #Legenda  do estado do eixo do tempo
plt.legend(fontsize=10) # define fontes do gráfico
plt.tick_params(labelsize=10) #define o tamanho da legenda
plt.show() # plota o gráfico na tela 