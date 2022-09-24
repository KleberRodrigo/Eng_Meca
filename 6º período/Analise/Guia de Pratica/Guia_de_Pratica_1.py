import numpy as np # importando biblioteca numpy
import matplotlib.pyplot as plt # importando biblioteca para plotar as figuras
import control as ct  #importanto biblioteca control

plt.close('all') #comando para fechar todas janelas de plot

#Criando um função para cálculo da equação diferencial do sistema------------------------------------
def model_update(t,x,u,params):
    # parametros do sistema
          # Resistencia da fornalha
    
    #Atenção: A vazão de entrada Q1 (Kg/S) é dada pela variável interna u, enquanto o estado do sistema, H,
    #é descrita na variável interna x. 
    T=50
    Rh=40
    Me=(((np.log(1-(Rh/100))))/((-K*(T+C))))**(1/N)
    k=(0.1579+(0.0001746*T))-(0.0001413*(Rh))
    n=(0.6545+(0.002425*T)+(0.0007867*(Rh)))
    dx = -k*(x-Me)

    return dx

def model_output(t,x,u,params):
    #Atenção os dados retornados nessa função estarão disponíveis para aquisição e tratamento externo. 
    #Para o caso em estudo, a saída do sistema, y, é o próprio estado, ou seja, y = H
    return x

#Definindo o sistema não linear do tanque, segundo biblioteca----------------------------------------
#Observe a notação adotada para os seguintes parâmetros: name, inputs, outputs e states. Esses
#parâmetros serão utilizados para realização da conexão dos blocos da malha fechada, quando  necessário.
SYSTEM = ct.NonlinearIOSystem(model_update, model_output , states=1, name='SYSTEM',inputs = ('u'), outputs = ('y'))

#Definindo os parâmetros de simulação----------------------------------------------------------------
# tempo    
t0 = 0    # tempo inicial
tf = 20  # tempo final
t = np.linspace(t0,tf,10001) # instantes que desejo ter a solucao

K = 1.9187*10**-5
N = 2.4451
C = 51.161
T = 50
Rh = 40

M0=28

#Sinal de controle necessário para levar o sistema para o 
#ponto de operação desejado em malha aberta.
#Condição inicial do sistema.

#Executando a simulação do sistema não linear em malha aberta----------------------------------------
#Observe que para simulação em malha aberta a função exige os seguintes parâmetros:
#Sistema a ser simualado, vetor com o tempo de simulação, vetor sinal de controle, condição
#inicial do sistema.
u=T

t, y = ct.input_output_response(SYSTEM, t, u, M0) 

#Plotando o resultado da simulação-------------------------------------------------------------------
plt.figure(1)
plt.plot(t,y,'b',label='T(t)')
plt.ylabel('T(t)[ºC]')
plt.xlabel('Tempo[s]')
plt.legend()
plt.title('Resposta temporal da umidade')
plt.show()
