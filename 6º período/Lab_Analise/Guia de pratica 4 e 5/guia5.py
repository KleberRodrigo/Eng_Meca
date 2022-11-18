""" Autor: Kleber Junior e RObson Junior 
    Criado em : 29/10/2022 """

import  numpy as np #Importa biblioteca Numpy
from    matplotlib import pyplot as plt  #Importa biblioteca pyplot 
import  control as ct   # Importa a biblioteca control 
import math as m  # Importa a biblioteca math 
plt.close('all') # fecha todos os graficos e plots 

def model_update(t,x,u,params): #Definição das equações diferenciais do espaço de estados
    x1 = x[0] # Espaço de estado x1
    x2 = x[1] # Espaço de estado x2
    return np.array([x2, K_1*(np.cos(x1)**2)*u - (K_2*np.sin(x1) + K_3*x2)]) # Retorna um vetor contendo a representação no espaço de estado

def model_output(t,x,u,params): #Função para definir o modelo da saída
    return x[0] #Retorna x[0]

def K_g(massa): #Criando função Para tratar as constante do sistema
    K_1 = (d_cm*rho_ar*C_a*L_a*L_1)/(2*massa*(((L_t**2)/12)+d_cm**2)) #Equação de K1
    K_2 = (g*d_cm)/ (((L_t**2)/12)+d_cm**2)#Equação de K2
    K_3 = (atrito*(d_cm**2))/(massa*(((L_t**2)/12)+d_cm**2))#Equação de K3
# Todas as equações saíram do guia de prática 5

    return K_1, K_2, K_3 # Retorna os Valores de K1, K2 e K3

L_a = 0.154 #Largura da placa móvel 
L_1 = 0.155 #Comprimento da placa abaixo do eixo de rotação
d_cm = 0.020 # Distancia do centro de massa da placa 
L_t = 0.270 #Comprimento da placa abaixo do eixo de rotação
rho_ar = 1.23 #Densidade do ar 
m_t = 0.100 # Massa Total da placa
C_a = 2.05#Coeficiente de arrasto
g = 9.81#Gravidade
atrito = 5#Coeficiente de atrito
ponto_inicial = 0#Ponto inicial da placa em graus 
pt_operacao = 50 #Ponto de operação da placa em graus
pt_operacao_rad = m.radians(pt_operacao) #Ponto de operação em radianos

K_1, K_2, K_3 = K_g(m_t) #Definindo as Variaveis das constantes

tmp_final = 24 # tempo final de Simulação
periodo = 0.001 #Periodo entre uma posição de vetor e outra 
time = np.arange(0, tmp_final, periodo)# Organiza o vetor com parametos das variaveis já definida
U_ent = (K_2 * m.sin(pt_operacao_rad))/(K_1 * (m.cos(pt_operacao_rad)**2)) #Definindo a equação do u de entrada
U_0 = U_ent * np.ones(time.shape) #Criando vetor do tamanho de time com os valores de u de entrada
ref = pt_operacao_rad * np.ones(time.shape) #Criando vetor do tamanho de time e com valores de 1

Degraus_U = U_ent * np.ones(len(time)) # Definindo vetor de degraus sequenciado
Degraus_U = np.array_split(Degraus_U, 3) #Dividindo o vetor em 3 partes

Degraus_U[0][:]=U_0[0] #Primeira Parte do Vetor
Degraus_U[1][:]= 1.20*U_0[0] #Segunda Parte do Vetor 
Degraus_U[2][:]= U_0[0] # Terceir parte do Vetor
Degraus_U_conc = np.concatenate([Degraus_U[0], Degraus_U[1], Degraus_U[2]]) # Junta todas as partes em um 

#função da python control para calcular EDO
FanPlate = ct.NonlinearIOSystem(model_update, model_output , states=2, name='FanPlate', inputs = ('u'), outputs = ('y')) #Passa parametro
t, y = ct.input_output_response(FanPlate, time, U_0, ponto_inicial)#Retorna resolvido

plt.figure('Sinal de Controle') #Plota a figura do gráfico
plt.subplot(2,1,1) #Divide em subplot
plt.plot(t,ref,'--k', label='referência(0,873 radianos)') # plotando o valor de referencia
plt.plot(t,y,'r', label='$\\theta(t)$') # Plotando o resultado de theta
plt.ylabel('$\\theta$(t)[rad]')#Define a legenda do eixo y
plt.xlim(0,5) #Define o limite do eixo x
plt.ylim(0,1.45) #Define o limite do eixo y
plt.legend() #Plota a legenda

plt.grid() #Ativa as linha de grade
plt.subplot(2,1,2) #Divide em subplot
plt.plot(time,U_0,'r', label='$u_0(t)$') #Plota a entrada
plt.ylabel('U(t)[m^2/s^2]') #Legenda do eixo y
plt.xlabel('Tempo[s]') #Legenda do eixo x
plt.legend() #Plota a Legenda
plt.xlim(0,20) # Define o limite do eixo x
plt.ylim(56,65)# Define o limite do eixo y
plt.grid() # Ativa as linha de grade
plt.show()#Plota o gráfico

t_degrau, y_degrau = ct.input_output_response(FanPlate, time, Degraus_U_conc, ponto_inicial) #passando o degrau como entrada
plt.figure('Resposta ao degrau de 20% do sinal de controle') #Plota a figura
plt.subplot(2,1,1)#Divide em subplot
plt.plot(t_degrau,ref,'--k', label='referência')#Passa o vetor de resposta no tempo
plt.plot(t_degrau,y_degrau,'r', label='u(t)')#Plota o resultado da EDO
plt.ylabel('$\\theta$(t)[rad]') # Define a legenda do eixo y
plt.xlim(0,20) #define o limite do eixo x
plt.ylim(0,1.45)#Define o limite do eixo y
plt.legend()#Plota a legenda

plt.grid() # Ativa as linhas de grade 
plt.subplot(2,1,2) #Divide em subplots 
plt.plot(time,Degraus_U_conc,'r', label='u(t)') # Imprime a entrada em degraus
plt.ylabel('U(t)[m^2/s^2]')#Legenda do eixo y
plt.xlabel('Tempo[s]')#LEgenda do eixo x
plt.legend() #Imprime Legenda
plt.xlim(0,20)#Impoe o limite de x
plt.ylim(60,75) # define o limite de y
plt.grid()#Define as linhas de grade
plt.show()#plota o grafico

y_n = y_degrau[8000:16000] # faz a normalização da curva
t_n = t_degrau[0:len(y_n)] # faz a normalização

theta_n = (y_n - min(y_n))/(y_degrau[11999] - min(y_n)) # equação de theta normalizado

ref_n = np.ones(t_n.shape) #vetor unitário
ref_tsp_2 = 1.02*np.ones(t_n.shape) #Vetor com 2% acima
ref_tsn_2 = 0.98*np.ones(t_n.shape) #Vetor com 2% a menos 
ref_tsp_5 = 1.05*np.ones(t_n.shape) #Vetor com 5% a mais
ref_tsn_5 = 0.95*np.ones(t_n.shape) #Vetor com 5% a menos

tp = 0.348 #Instante que acontece o pico
b = 0.5877 #y(tp)
Mp = b/1 #overshoot
tr = 0.1909# tempo de pico
ts = 2.498#tempode acomodação
zeta = 0.17465 #Coeficiente de amortecimento
wn = 9.16848 #frequencia natural do sistema

E_s = 1 + (np.e**((-zeta)*wn*t_n))/(np.sqrt(1-(zeta**2)))#equação da curva involtória superior
E_i = 1 - (np.e**((-zeta)*wn*t_n))/(np.sqrt(1-(zeta**2)))#equação da curva involtória inferior

plt.figure('Curvas envoltórias para a curva de resposta ao degrau unitário') #Plotando as curvas envoltorias
plt.plot(t_n,ref_n,'--k',label ='referência') #referencia
plt.plot(t_n,theta_n,'r',label ='$\\theta$(t)') #referencia
plt.plot(t_n,ref_tsp_2,'--g',label ='referência $\pm$ 2%') #referencia
plt.plot(t_n,ref_tsn_2,'--g') #referencia
plt.plot(t_n,ref_tsp_5,'--y',label='referência $\pm$ 5%') #referencia
plt.plot(t_n,ref_tsn_5,'--y') #referencia
plt.plot(t_n,E_s,'b', label = 'Curvas envoltórias') #curva envoltoria
plt.plot(t_n,E_i,'b') #curva envoltoria
plt.ylabel('$\\theta(t)$[rad]') #Legenda do eixo y
plt.xlabel('Tempo [s]')#Legenda do Eixo x
plt.xlim(0,3.5) #Limite do eixo x
plt.ylim(0,2)#Limite do eixo y
plt.legend()#Imprime Legenda
plt.grid()#Imprime linhas de grade
plt.show()#Plota o gráfico

Degraus_U_eq_transf = U_ent * np.ones(len(time)) #Definindo degraus
Degraus_U_eq_transf = np.array_split(Degraus_U_eq_transf, 8) # Vetor de degrau

Degraus_U_eq_transf[0][:]=U_ent #Degrau 1
Degraus_U_eq_transf[1][:]=1.10*U_ent#Degrau 2
Degraus_U_eq_transf[2][:]=0.90*U_ent#Degrau 3
Degraus_U_eq_transf[3][:]= 1.20*U_ent#Degrau 4
Degraus_U_eq_transf[4][:]= 0.80*U_ent#Degrau 5
Degraus_U_eq_transf[5][:]= U_ent#Degrau 6
Degraus_U_eq_transf[6][:]= 1.3*U_ent#Degrau 7
Degraus_U_eq_transf[7][:]= 0.70*U_ent#Degrau 8

Degraus_U_conc_eq_transf = np.concatenate([Degraus_U_eq_transf[0], Degraus_U_eq_transf[1], Degraus_U_eq_transf[2], Degraus_U_eq_transf[3], Degraus_U_eq_transf[4], Degraus_U_eq_transf[5], Degraus_U_eq_transf[6], Degraus_U_eq_transf[7]]) #Juntando todos os veotores em um

t_degrau2, y_degrau2 = ct.input_output_response(FanPlate, time, Degraus_U_conc_eq_transf, 0) #Passando parametro de entrada para a equação diferencial

s=ct.tf("s") #Cria variavel s para função de transferencia
K = abs((m.radians(50) - y_degrau[9999])/(1.20*U_ent- U_ent)) #Ganho estático
Eq_transferencia = K*(wn**2)/((s**2)+(2*s*wn*zeta)+(wn**2)) #Função de tranferencia

print() #Printa valor vazio
print(Eq_transferencia)#Printa a função de transferencia
print() #Printa valor vazio

timeTransfer, tempOut = ct.forced_response(Eq_transferencia,T=time, U=Degraus_U_conc_eq_transf-U_ent) # resposta forçada
timeTransfer, thetanonL = ct.input_output_response(FanPlate, time, Degraus_U_conc_eq_transf, [pt_operacao_rad,0])#Sistema real

plt.figure('Validação da resposta do modelo obtido')#Abre a Figura com legenda
plt.subplot(211)#divide em Subplot
plt.plot(timeTransfer,tempOut+(np.radians(50)),'r',label='Modelo') # plota a curva 1
plt.plot(t_degrau2,thetanonL,'y', label='Real')# plota a curva real
plt.plot(time,ref,'--k',label='ref(0.873rad)') # plota a referencia
plt.grid()#Ativa as linhas de grade
plt.ylabel('$\\theta$(t)[rad]') #Legenda do eixo y
plt.legend()#Imprime legenda 
plt.xlim(0,24) #Limite de x
plt.ylim(0.65,1.05) # Limite de y

plt.subplot(212) #subplot da entrada
plt.plot(timeTransfer,Degraus_U_conc_eq_transf,'r', label = '$Q(t)$')#Plota o vetor concatenado
plt.grid() # Ativa as linha de grade
plt.ylabel('${u(t)}[m^2/s^2]$')#Legenda de y
plt.xlabel('Tempo [s]')#Legenda de x
plt.legend()#Imprime legendas
plt.xlim(0,24) #Ajusta o limite de x
plt.ylim(40,80)#Ajusta o limite de y
plt.show() #Plota o grafico 