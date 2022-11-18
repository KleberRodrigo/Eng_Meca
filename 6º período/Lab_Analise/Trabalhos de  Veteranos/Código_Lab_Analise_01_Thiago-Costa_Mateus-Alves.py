#@Autor: Thiago Henrique de Faria Costa e Mateus Alves de Sales
#@Data: 28/11/2021

import numpy as np                      #\
from matplotlib import pyplot as plt    # \   importando as bibliotecas necessárias
import control as ct                    # /
import statistics as sta                #/

plt.close('all') #fechar todas as abas de gráficos

#declaração dos parâmetros utilizados
def model_uptade(t,T,Q,params):  
    V = params.get('V', 10)         #volume (m^3)
    c = params.get('c', 4500)       #calor específico a pressão constante (J/KgºC)
    F = params.get('F', 3)          #vazão volumétrica (m^3/s)
    T_0 = params.get('T_0', 20)     #temperatura ambiente
    h = params.get('h', 15)         #coeficiente de convecção natural (J/sm^2ºC)
    A = params.get('A', 31.4)       #superfície do tanque para a troca de calor por convecção (m^2)
    rho = params.get('rho', 1000)   #massa específica da solução (Kg/m^3)

    dT = (Q/(rho*V*c))+(h*A*(T_0-T)/(rho*V*c))+((F*(T_0-T)/V))
    return dT

s=ct.tf("s")
rho =1000
c = 4500
F = 3
V = 10
h = 15
A = 31.4
T_0 = 20

def model_output(t,T,Q,params): #retorna a saída de acordo com os parâmetros fornecidos em "model_update"
    return T

#Definição do sistema não linear
SYSTEM = ct.NonlinearIOSystem(model_uptade, model_output,
states=1, name='SYSTEM', inputs=('u'), outputs=('y'))

tempo_final = 500                                 #duração total da simulação
time = np.linspace(0, tempo_final, tempo_final+1) #definição de um vetor tempo com espaçamento de 1s
Q0 = 270009420 * np.ones(time.shape)              #sinal de controle necessário para levar o sistema ao ponto de operação desejado -> 40ºC
ref = 40 * np.ones(time.shape)                    #cria a referência para o ponto de operação



#simulação do sitema não linear em malha aberta
T,y = ct.input_output_response(SYSTEM, time, Q0, T_0) #calcula a resposta de saída do sistema

#Plotando a resposta temporal
plt.figure(1)
plt.subplot(2,1,1)
plt.plot(T,ref,'--k', label='ref')
plt.plot(T,y,'b', label='T(t)')
plt.ylabel('T(t)[C]')
plt.xlim(0,30)
plt.legend()
plt.ylim(0,100)
plt.title('Resposta temporal de tanques acoplados em malha aberta')
plt.grid()
plt.subplot(2,1,2)
plt.plot(time,Q0,'b', label='Q0')
plt.ylabel('Q(t)')
plt.legend()
plt.xlabel('Tempo[s]')
plt.xlim(0,30)
plt.ylim(0,1000000000)
plt.grid()
plt.show()

#Modelo de Frequência ---

G=1/((s+(F/V+h*A/(V*rho*c)))*(rho*c*V)) #equação obtida pela transformada de Laplace
t, ysp = ct.forced_response(G,T=time, U = Q0) # Validação da parametrização pelo degrau de subida pelo método dos três parâmetros
ysp += T_0 # Somando o referencial

#Plotando a função de transferência
plt.figure(1)
plt.subplot(2,1,1)
plt.plot(t,ref,'--k', label='ref')
plt.plot(t,ysp,'b', label='T(s)')
plt.ylabel('T(s)[C]')
plt.xlim(0,30)
plt.legend()
plt.ylim(0,100)
plt.title('Equação de transferência')
plt.grid()
plt.subplot(2,1,2)
plt.plot(time,Q0,'b', label='Q0')
plt.ylabel('Q(t)')
plt.legend()
plt.xlabel('Tempo[s]')
plt.xlim(0,30)
plt.ylim(0,1000000000)
plt.grid()
plt.show()

#Degraus

temperatura_transf = np.empty(len(time)) #alocando o vetor para a temperatura do liquido
temperatura_transf.fill(np.nan)          #limpando o vetor e definindo-o como NAN

Q2 = Q0[0] * np.ones(len(time))          #criando um novo vetor de Q0

Q2[0:125]= Q0[0]
Q2[125:250]=1.15*Q0[0]                   #degrau positivo -> +15%
Q2[250:375]= Q0[0]
Q2[375:500]=0.85*Q0[0]                   #degrau negativo -> -15%

temperatura_transf[0] = 20               #definindo o primeiro valor do vetor com o valor da temperatura ambiente

for k in range(len(time)-1):             #integração pelo método de Euler
    temperatura_transf[k+1] = temperatura_transf[k]+((F*(T_0 - temperatura_transf[k])/V)+ (Q2[k]/(rho*V*c)) + ((h*A)*(T_0 - temperatura_transf[k])/(rho*V*c)))

#Plotando os degraus
plt.figure(3)
plt.subplot(211)
plt.title('Degraus')
plt.plot(time,temperatura_transf,'b',label='$\\t$')
plt.plot(time,ref,'--k',label='ref')
plt.grid()
plt.ylabel('T(t)')
plt.legend()
plt.ylim(15,55)
plt.xlim(0,tempo_final)
plt.subplot(212)
plt.plot(time,Q2,'b',label='Q1')
plt.grid()
plt.ylabel('Q [J] ')
plt.xlabel('Tempo [s]')
plt.legend()
plt.ylim(0,1000000000)
plt.xlim(0,tempo_final)
plt.show()

#Degrau positivo
c_a1 = (temperatura_transf[126]-temperatura_transf[125])/(time[126]- time[125]) #coeficiente angular da curva tangente ao degrau de subida

eq_t1 = c_a1*time - c_a1*time[126] + temperatura_transf[126]                    #equação da reta tangente ao degrau de subida

Kp_dp = temperatura_transf[249]*np.ones(len(time))                              #crinado vetor do ponto em que o sistema se estabiliza após o primeiro degrau positivo.

K1 = (temperatura_transf[249]-40)/(Q2[249]-270009420)                           #calculando K do modelo de degrau positivo

tau_1 = 130-125                                                                 #tau relativo ao primeiro degrau positivo -> análise gráfica

#como não existe atraso, a variável temperatura_transf contida no método é igual a zero

GN_DP = K1/(tau_1*s+1) #modelo do degrau positivo

t, yps = ct.forced_response(GN_DP,T= time, U = Q2)
yps += T_0

#Degrau negativo
c_a2 = (temperatura_transf[376]-temperatura_transf[375])/(time[376]- time[375]) #coeficiente angular da curva tangente ao degrau de descida

eq_t2 = c_a2*time - c_a2*time[376] + temperatura_transf[376]                    #equação da reta tangente ao degrau de descida

Kp_dp2 = temperatura_transf[499]*np.ones(len(time))                             #crinado vetor do ponto em que o sistema se estabiliza após o primeiro degrau negativo.

K2 = (temperatura_transf[499]-temperatura_transf[374])/(Q2[499]-Q2[374])        #calculando K do modelo de degrau negativo

tau_2 = 380-375                                                                 #tau relativo ao primeiro degrau negativo -> análise gráfica

#como não existe atraso, a variável temperatura_transf contida no método é igual a zero

GN_DN = K2/(tau_2*s+1) #modelo do degrau negativo

t, ypn = ct.forced_response(GN_DN,T= time, U = Q2)
ypn += T_0

#Plotando os degraus positivos e negativos -> método de Ziegler
plt.figure(4)
plt.subplot(211)
plt.title('Degrau positivo')
plt.plot(time,temperatura_transf,'b',label='T')
plt.plot(t,yps,'g',label='T(s)')
plt.plot(time,eq_t1,'r',label='eq_t1')
plt.plot(time,ref,'--k',label='ref')
plt.plot(time,Kp_dp,'--y',label='K')
plt.grid()
plt.ylabel('T[C]')
plt.legend()
plt.ylim(35,45)
plt.xlim(120,140)
plt.subplot(212)
plt.plot(time,temperatura_transf,'b',label='T')
plt.plot(t,ypn,'g',label='T(s)')
plt.plot(time,eq_t2,'r',label='eq_t2')
plt.plot(time,ref,'--k',label='ref')
plt.plot(time,Kp_dp2,'--y',label='K')
plt.grid()
plt.xlabel('Tempo [s]')
plt.ylabel('T[C]')
plt.legend()
plt.ylim(35,45)
plt.xlim(370,390)
plt.show()

#Método de Miller
#degrau positivo

tau_3 = 0.63*tau_1

GM_DP = K1/(tau_3*s+1) #modelo do degrau positivo

t, yes = ct.forced_response(GM_DP,T= time, U = Q2)
yes += T_0

#degrau negativo

tau_4 = 0.63*tau_2

GM_DN = K2/(tau_4*s+1) #modelo do degrau negativo

t, yfs = ct.forced_response(GM_DN,T= time, U = Q2)
yfs = yfs + T_0

#Plotando os degraus relativos ao método de Miller
plt.figure(5)
plt.subplot(211)
plt.title('Degrau positivo Miller')
plt.plot(time,temperatura_transf,'b',label='T')
plt.plot(t,yes,'g',label='T(s)')
plt.plot(time,eq_t1,'r',label='eq_t1')
plt.plot(time,ref,'--k',label='ref')
plt.plot(time,Kp_dp,'--y',label='K')
plt.grid()
plt.ylabel('T[C]')
plt.legend()
plt.ylim(35,45)
plt.xlim(120,140)
plt.subplot(212)
plt.plot(time,temperatura_transf,'b',label='T')
plt.plot(t,yfs,'g',label='T(s)')
plt.plot(time,eq_t2,'r',label='eq_t2')
plt.plot(time,ref,'--k',label='ref')
plt.plot(time,Kp_dp2,'--y',label='K')
plt.grid()
plt.xlabel('Tempo [s]')
plt.ylabel('T[C]')
plt.legend()
plt.ylim(35,45)
plt.xlim(370,390)
plt.show()

#Plotando todos modelos juntos com a EDO
Q3 = 270009420*np.ones(len(time)) #criando o vetor para a quantidade de calor inicial

Q3[0:125]= Q0[0]
Q3[125:250]=1.07*Q0[0]
Q3[250:375]= Q0[0]
Q3[375:500]=0.93*Q0[0]

for k in range(len(time)-1): #integração pelo método de Euler 
    temperatura_transf[k+1] = temperatura_transf[k]+1*((F*(T_0 - temperatura_transf[k])/V)+ (Q3[k]/(rho*V*c)) + ((h*A)*(T_0 - temperatura_transf[k])/(rho*V*c))) #em que +1 é o tempo de amostragem, em segundos

t, yps = ct.forced_response(GN_DP,T= time, U = Q3)
yps += T_0

t, ypn = ct.forced_response(GN_DN,T= time, U = Q3)
ypn += T_0

t, yes = ct.forced_response(GM_DP,T= time, U = Q3)
yes += T_0

t, yfs = ct.forced_response(GM_DN,T= time, U = Q3)
yfs += T_0


#Plotando todos os gráficos juntos com a EDO
plt.figure(6)
plt.subplot(211)
plt.title('Todos os métodos')
plt.plot(time,temperatura_transf,'b',label='EDO') #modelo Edo  
plt.plot(t,yes,'r',label='Ziegler')               #modelo degrau positivo -> método de Ziegler
plt.plot(t,ypn,'y',label='Miller')                #modelo degrau negativo -> método de Ziegler
plt.plot(time,ref,'--k',label='ref')
plt.grid()
plt.xlabel('Tempo [s]')
plt.ylabel('T[C]')
plt.legend()
plt.ylim(25,50)
plt.xlim(0,tempo_final) 
plt.show()

#Comparando os sistemas pelo índice RMSE
RMSE = np.zeros(4)

real = temperatura_transf[125:500]                # T_0 no instante dos degraus

GN_DP = yps[125:500]                              # Erro da curva parametrizada pelo método de Ziegler -> degrau de subida
RMSE[0] = np.sqrt(sta.mean(((GN_DP-real)**2)))

GN_DN = ypn[125:500]                              # Erro da curva parametrizada pelo método de Ziegler -> degrau de descida 
RMSE[1] = np.sqrt(sta.mean(((GN_DN-real)**2)))

GM_DP = yes[125:500]                              # Erro da curva parametrizada pelo método de Miller -> degrau de subida
RMSE[2] = np.sqrt(sta.mean(((GM_DP-real)**2)))

GM_DN = yfs[125:500]                              # Erro da curva parametrizada pelo pelo método de Miller -> degrau de descida 
RMSE[3] = np.sqrt(sta.mean(((GM_DN-real)**2)))

print('Índice RMSE para o primeiro método de Ziegler degrau positivo', RMSE[0])
print('Índice RMSE para o primeiro método de Ziegler degrau negativo', RMSE[1])
print('Índice RMSE para o modelo de Miller degrau positivo          ', RMSE[2])
print('Índice RMSE para o modelo de Miller degrau negativo          ', RMSE[3])