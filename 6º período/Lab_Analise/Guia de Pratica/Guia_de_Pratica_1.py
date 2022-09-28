import numpy as np # importando biblioteca numpy
import matplotlib.pyplot as plt # importando biblioteca para plotar as figuras
import control as ct  #importanto biblioteca control


plt.close('all')

def model_update(t,x,u,params):
    
    Rh = 0.20 #umidade relativa
    T = 45 # temperatura de secagem
    C = 51.161 #constante 
    K = 1.9187 * (10**(-2)) #constante
    N = 2.4451 #constante
    
    k = 0.2958 + (0.01215*T) + (0.0193653*Rh) #constante
    Me = ((np.log(1-Rh))/(-K*(T+C)))**(1/N) #umidade de equilibrio
    
    dx =(x - Me)*(-k) #equação diferencial da umidade
    
    return dx          

def model_output(t,x,u,params):
   return x

SYSTEM = ct.NonlinearIOSystem(model_update, model_output , states=1, name='SYSTEM',inputs = ('u'), outputs = ('y'))
   
t0 = 0    # tempo inicial
tf = 1  # tempo final
t = np.linspace(t0,tf,10001) # instantes que desejo ter a solucao
M0=28 #condição inicial
u=0

t, y = ct.input_output_response(SYSTEM, t, u, M0) 

#Plotando o resultado da simulação-------------------------------------------------------------------
plt.figure(1)
plt.plot(t,y,'b',label='Umidade(t)')
''' #plot com variação da temperatura
def model_update1(t,x,u,params):
    
    Rh = 0.20
    T = 90
    C = 51.161
    K = 1.9187 * (10**(-2))
    N = 2.4451
    
    k = 0.2958 + (0.01215*T) + (0.0193653*Rh)
    Me = ((np.log(1-Rh))/(-K*(T+C)))**(1/N)
    
    dx =(x - Me)*(-k)
    
    return dx          

def model_output1(t,x,u,params):
   return x

SYSTEM = ct.NonlinearIOSystem(model_update1, model_output1 , states=1, name='SYSTEM',inputs = ('u'), outputs = ('y'))
   
M0=28
u=0

t, y = ct.input_output_response(SYSTEM, t, u, M0) 

#Plotando o resultado da simulação-------------------------------------------------------------------
plt.plot(t,y,'r',label='Umidade(t)')
'''
plt.ylabel('Umidade(t)[%]')
plt.xlabel('Tempo[h]')
plt.legend()
plt.show()