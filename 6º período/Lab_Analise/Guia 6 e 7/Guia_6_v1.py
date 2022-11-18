#Autores: Kleber Junior E Robson Junior
#Laboratório de Sistemas Lineares

import  numpy as np #importa a biblioteca numpy
from    matplotlib import pyplot as plt #Importa a biblioteca pyplot para os gráficos
import  control as ct #Importa a biblioteca control para fazer os calculos
import  math 

plt.close('all') #Fecha todos os gráficos abertos
def model_update(t,x,vetor_u,params): #Ajuste de parametros
    u = vetor_u[0]
    u_massa = vetor_u[1]
    MT = m_t*(1 + u_massa)

    K_1 = (d_cm*rho_ar*C_a*L_a*L_1)/(2*MT*(((L_t**2)/12)+d_cm**2)) #Constante passada do sistema 
    K_2 = (g*d_cm)/ (((L_t**2)/12)+d_cm**2) #Constante passada do sistema
    K_3 = (atrito*(d_cm**2))/(MT*(((L_t**2)/12)+d_cm**2)) #Constante passada do sistema
#Variavéis do espaço de estado
    x1 = x[0]
    x2 = x[1]
    return np.array([x2, K_1*(np.cos(x1)**2)*u - (K_2*np.sin(x1) + K_3*x2)], dtype=object)

def model_output(t,x,u,params): #Ajuste da Saída 
    return x[0]

def K_ganho(massa):  #Definindo função para o calculo de k1,k2 e k3
    K_1 = (d_cm*rho_ar*C_a*L_a*L_1)/(2*massa*(((L_t**2)/12)+d_cm**2)) #Equação de K1 segundo o guia de prática
    K_2 = (g*d_cm)/ (((L_t**2)/12)+d_cm**2)#Equação de K2 segundo o guia de prática
    K_3 = (atrito*(d_cm**2))/(massa*(((L_t**2)/12)+d_cm**2))#Equação de K3 segundo o guia de prática
    return K_1, K_2, K_3

def coeficiente_angular(valor_y, valor_x, j): #Calculo do coeficiente angular
    return (valor_y[j]-valor_y[j-1])/(valor_x[j]-valor_x[j-1])

def valor_K(vetor1, vetor2, ponto_de_operacao, variavel_de_controle, k): #Calculo do resultado de K
    return vetor1[k]-ponto_de_operacao, vetor2[k]-variavel_de_controle, (vetor1[k]-ponto_de_operacao)/(vetor2[k]-variavel_de_controle)

def equacao(vetor1, vetor2, c_angular, k): #Retorna a reta tangente da curva de ziegler
    return c_angular*vetor1 - c_angular*vetor1[k] + vetor2[k]


#Variáveis do sistema
L_a = 0.154
L_1 = 0.155
L_t = 0.270
d_cm = 0.020
rho_ar = 1.23
m_t = 0.005
C_a = 2.05
atrito = 5
g = 9.81
pt_operacao = 43
pt_operacao_rad = math.radians(pt_operacao)
pt_inicial = 0
#Fim variaveis do sistema

K_1, K_2, K_3 = K_ganho(m_t) #Chamando função do calculo das constantes
print(f'K1: ', K_1) #Imprime o valor de K1
print(f'K2: ', K_2)#Imprime o valor de K2
print(f'K3: ', K_3)#Imprime o valor de K3

tmp_final = 60 # Definindo tempo de simulação para 1 minuto
periodo = 0.001 #Definindo o time step
time = np.arange(0, tmp_final, periodo) # Cria um Vetor de 1 a 60
U_ent = (K_2 * math.sin(pt_operacao_rad))/(K_1 * (math.cos(pt_operacao_rad)**2)) #Equação de entrada
U0 = U_ent * np.ones(time.shape)# Criando um vetor de mesmo valor de U_ent
referencia = pt_operacao_rad * np.ones(time.shape)#Criando um vetor do ponto de operação
u_massa = np.zeros(len(time)) #Criando  um vetor de zeros

print(f'Sinal de entrada:',U_ent) #Printa o valor de U que foi calculado

Degraus_U = U_ent * np.ones(len(time)) #Criando vetor para alocar os degraus
Degraus_U = np.array_split(Degraus_U, 4) #Divide o vetor em 4 partes

Degraus_U[0][:]=U0[0] #Degrau inicial
Degraus_U[1][:]= 1.20*U0[0] #Degrau com 20% acima
Degraus_U[2][:]= U0[0]#Degrau com o mesmo valor do u de entrada
Degraus_U[3][:]= 0.80*U0[0]#Degrau com 20% menor

Degraus_U_conc = np.concatenate([Degraus_U[0], Degraus_U[1], Degraus_U[2], Degraus_U[3]]) #Função de concatenar Vetores

FanPlate = ct.NonlinearIOSystem(model_update, model_output , states=2, name='FanPlate', inputs = ('u', 'u_massa'), outputs = ('y')) #Função para 
t, y = ct.input_output_response(FanPlate, time, [U0, u_massa], pt_inicial) #Simulando sem os degraus

#Primeiro Plot sem Degrau
plt.figure(1)
plt.title('Resposta temporal do sistema em malha aberta sem degrau')
plt.subplot(2,1,1)
plt.plot(t,referencia,'--b', label='referencia')
plt.plot(t,y[0],'r', label='$\\theta$(t)')
plt.ylabel('$\\theta$(t)[rad]')
plt.xlim(0,10)
plt.ylim(0,1.5)
plt.legend()
plt.grid()
plt.subplot(2,1,2)
plt.plot(time,U0,'r', label='Q (t)')
plt.ylabel('u(t)[m^2/s^2]')
plt.xlabel('Tempo[s]')
plt.legend()
plt.xlim(0,10)
plt.ylim(2,2.15)
plt.grid()
plt.show()
#Fim do Plot sem degrau

#Resposta temporal do sistema em malha aberta com degrau
t_degrau, y_degrau = ct.input_output_response(FanPlate, time, [Degraus_U_conc, u_massa], pt_inicial)
y_max_deg = y_degrau[0][29499]
plt.figure(2)
plt.subplot(2,1,1)
plt.plot(t_degrau,referencia,'--b', label='referencia')
plt.plot(t_degrau,y_degrau[0],'r', label='T(t)')
plt.ylabel('$\\theta$(t)[rad]')
plt.xlim(0,tmp_final)
plt.ylim(0,1.5)
plt.legend()
plt.grid()
plt.subplot(2,1,2)
plt.plot(time,Degraus_U_conc,'r', label='Q (t)')
plt.ylabel('u(t)[m^2/s^2]')
plt.xlabel('Tempo[s]')
plt.legend()
plt.xlim(0,tmp_final)
plt.ylim(0,4)
plt.grid()
plt.show()
#Fim do plot com degrau 

#Criando Função de tranferencia para o metodo de ziegler
s = ct.TransferFunction.s #Define s para variavel 
#Reta tangente 
amp_degrau_pos = y_degrau[0][29999]#Amplitude do degrau positivo
ca_degrau_pos = coeficiente_angular(y_degrau[0], time, 15100)#Coeficiente angular 
equacao_degrau_positivo = equacao(time, y_degrau[0], ca_degrau_pos, 15100)#Equação da reta
referencia_2 = amp_degrau_pos * np.ones(time.shape)#Criando vetor para referencia
#Determinando os parâmetros de Ziegler 
tau = 16.02 - 15.01 #Variavel tau
delta_theta, delta_u, K1 = valor_K(y_degrau[0], Degraus_U_conc, pt_operacao_rad, U_ent, 29999)#Função de resolução
print(f'delta y: ',delta_theta) #Imprimindo o valor davariação de theta
print(f'delta u: ',delta_u)#Imprimindo a variação de u
print(f'K1: ', K1)#Imprimindo o valor de K1

#Modelo primeira order de Ziegler para o degrau positivo
modelo_degrau_positivo = K1/(tau*s+1)

t_degrau_positivo, y_degrau_positivo = ct.forced_response(modelo_degrau_positivo,T= time, U = Degraus_U_conc - U_ent) # Resolvendo a equação

#Inicio do plot do método de Ziegler para o degrau positivo - Gráfico Método de Ziegler-Nichols para degrau positivo
plt.figure(3)
plt.plot(time,y_degrau[0],'r',label='$\\theta$')
plt.plot(time,y_degrau_positivo + pt_operacao_rad,'g',label='$\\theta_{Ziegler-Nichols}$')
plt.plot(time,equacao_degrau_positivo,'g',label='')
plt.plot(time,referencia,'--b',label='referencia')
plt.plot(time,referencia_2,'--k',label='referencia')
plt.grid()
plt.ylabel('$\\theta$(t)[rad]')
plt.xlabel('Tempo[s]')
plt.legend()
plt.ylim(0.71,0.82)
plt.xlim(14.5,22)

#Definindo os limites do zoom, o quadrado pontilhado
x0 = 15.000
xf = 15.020
y0 = 0.7504
yf = 0.7508
plt.plot([x0,xf], [y0,y0], 'k--')
plt.plot([x0,xf], [yf,yf], 'k--')
plt.plot([x0,x0], [y0,yf], 'k--')
plt.plot([xf,xf], [y0,yf], 'k--')

#Inicio do zoom na plotagem 
a = plt.axes([0.45, 0.5, 0.15, 0.15]) # x-position, y-position, width, height
plt.xlim(x0,xf)
plt.ylim(y0,yf)
plt.grid()
plt.plot(time,y_degrau[0],'r',label='$g(s)$')
plt.plot(time,y_degrau_positivo + pt_operacao_rad,'g',label='$\\theta_{Ziegler-Nichols}$')
plt.plot(time,equacao_degrau_positivo,'g',label='referencia degrau')
plt.plot(time,referencia,'--b',label='referencia')
plt.show()
#Fim do zom na plotagem 

#Aproximação por Padé
Gi = K1/(tau*s+1) #Função de transferencia 

#Calculo pela 1ª ordem
n1 = 1 
N1 = ct.pade(0.15,n1) #Função da control para aproximação de padé
Gd1 = ct.TransferFunction(np.array(N1[0]),np.array(N1[1])) #Função de transferência.
Gr1 = Gi*Gd1
t1, y1 = ct.forced_response(Gr1,T= time, U = U0-U_ent)
y1 = y1 + pt_operacao_rad

#Calculo pela 3ª ordem
n3 = 3 
N3 = ct.pade(0.15,n3) #Função da control para aproximação de padé
Gd3 = ct.TransferFunction(np.array(N3[0]),np.array(N3[1])) #Função de transferência.
Gr3 = Gi*Gd3
t3, y3 = ct.forced_response(Gr1,T= time, U = U0-U_ent) 
y3 = y3 + pt_operacao_rad

#Calculo pela 5ª ordem
n5 = 5 
N5 = ct.pade(0.15,n5) #Função da control para aproximação de padé
Gd5 = ct.TransferFunction(np.array(N5[0]),np.array(N5[1])) #Função de transferência.
Gr5 = Gi*Gd5
t5, y5 = ct.forced_response(Gr1,T= time, U = U0-U_ent) 
y5 = y5 + pt_operacao_rad

#Calculo pela 9ª ordem
n9 = 9 
N9 = ct.pade(0.15,n9) #Função da control para aproximação de padé
Gd9 = ct.TransferFunction(np.array(N9[0]),np.array(N9[1])) #Função de transferência.
Gr9 = Gi*Gd9
t9, y9 = ct.forced_response(Gr1,T= time, U = U0-U_ent) 
y9 = y9 + pt_operacao_rad

#Espaço de estado
atraso_1 = ct.tf2io(Gd1, name='atraso_1', inputs='u', outputs='y1')#1ª ordem
atraso_3 = ct.tf2io(Gd3, name='atraso_3', inputs='u', outputs='y3')#3ª ordem
atraso_5 = ct.tf2io(Gd5, name='atraso_5', inputs='u', outputs='y5')#5ª ordem
atraso_9 = ct.tf2io(Gd9, name='atraso_9', inputs='u', outputs='y9')#9ª ordem

#construção da malha aberta interconectando sistemas
#1ª ordem
malha_aberta_01 = ct.InterconnectedSystem(
    (atraso_1, FanPlate), name='malha_aberta_01',
    connections = (('atraso_1.u',),('FanPlate.u','atraso_1.y1')),
    inplist = ('atraso_1.u'),
    outlist = ('FanPlate.y'))

X1 = np.zeros(n1+2) #para condições iniciais nulas
Ta, y1 = ct.input_output_response(malha_aberta_01, time, Degraus_U_conc, X1)#Chamando função i/o response

#3ª ordem
malha_aberta_03 = ct.InterconnectedSystem(
    (atraso_3, FanPlate), name='malha_aberta_03',
    connections = (('atraso_3.u',),('FanPlate.u','atraso_3.y3')),
    inplist = ('atraso_3.u'),
    outlist = ('FanPlate.y'))

X3 = np.zeros(n3+2) #para condições iniciais nulas
Ta, y3 = ct.input_output_response(malha_aberta_03, time, Degraus_U_conc, X3)

#5ª ordem
malha_aberta_05 = ct.InterconnectedSystem(
    (atraso_5, FanPlate), name='malha_aberta_05',
    connections = (('atraso_5.u',),('FanPlate.u','atraso_5.y5')),
    inplist = ('atraso_5.u'),
    outlist = ('FanPlate.y'))

X5 = np.zeros(n5+2) #para condições iniciais nulas
Ta, y5 = ct.input_output_response(malha_aberta_05, time, Degraus_U_conc, X5)

#9ª ordem
malha_aberta_09 = ct.InterconnectedSystem(
    (atraso_9, FanPlate), name='malha_aberta_09',
    connections = (('atraso_9.u',),('FanPlate.u','atraso_9.y9')),
    inplist = ('atraso_9.u'),
    outlist = ('FanPlate.y'))

X9 = np.zeros(n9+2) #para condições iniciais nulas
Ta, y9 = ct.input_output_response(malha_aberta_09, time, Degraus_U_conc, X9)

#Plotando resposta do sistema com atraso
plt.figure(4)
plt.plot(Ta,referencia,'--k',label='referencia')
plt.plot(Ta, y_degrau_positivo + pt_operacao_rad, 'black', label = 'Modelo - 3P')
plt.plot(Ta, y1, 'purple', label='Modelo - 1ª ordem') #atraso 1 
plt.plot(Ta, y3, 'b', label='Modelo - 2ª ordem') #atraso 3 
plt.plot(Ta, y5, 'r', label='Modelo - 5ª ordem') #atraso 5 
plt.plot(Ta, y9, 'g', label='Modelo - 9ª ordem') #atraso 9 
plt.legend()
plt.ylabel('$\\theta[rad]$')
plt.xlabel('Tempo[s]')
plt.ylim(0.74,0.82)
plt.xlim(14.5,21.5)
plt.grid()

#Definindo os limites do plot
x0_n = 14.95
xf_n = 15.20
y0_n = 0.7490
yf_n = 0.7510

#Plotando quadrado do zoom
plt.plot([x0_n,xf_n], [y0_n,y0_n], 'k--')
plt.plot([x0_n,xf_n], [yf_n,yf_n], 'k--')
plt.plot([x0_n,x0_n], [y0_n,yf_n], 'k--')
plt.plot([xf_n,xf_n], [y0_n,yf_n], 'k--')

#Fazendo o plot do Zoom
a = plt.axes([0.35, 0.5, 0.15, 0.15]) #posição e tamanho do grafico de zoom
plt.xlim(x0_n,xf_n)
plt.ylim(y0_n,yf_n)
plt.grid()
plt.plot(Ta,referencia,'--k',label='referencia')
plt.plot(Ta, y_degrau_positivo + pt_operacao_rad, 'black', label = 'Modelo - 3P') #3P
plt.plot(Ta, y1, 'purple', label='Modelo - 1ª ordem') #atraso 1ª ordem
plt.plot(Ta, y3, 'b', label='Modelo - 2ª ordem') #atraso 3ª ordem
plt.plot(Ta, y5, 'r', label='Modelo - 5ª ordem') #atraso 5ª ordem
plt.plot(Ta, y9, 'g', label='Modelo - 9ª ordem') #atraso 9ª ordem
plt.show()

"""Comparando o degrau positivo de Ziegler com a 5ª aproximação de Padé""" "Item 07"
#Gráfico de Comparação com Ziegler e a aproximação de 5ª ordem de padé
plt.figure(5)
plt.plot(Ta,referencia,'--b',label='referencia')
plt.plot(Ta, y_degrau_positivo + pt_operacao_rad, 'black', label = 'Modelo - 3P') #3P de Ziegler
plt.plot(Ta, y5, 'r', label='Modelo - 5ª ordem de Padé') #Modelo pela 5ª ordem de Padé 
plt.legend()
plt.ylabel('theta[rad]')
plt.xlabel('Tempo[s]')
plt.xlim(13.8, 26)
plt.ylim(0.675, 0.850)
plt.grid()
#Fim do gráfico
# Configuração dos limites do plot do zoom
x0_n = 14.950
xf_n = 15.200
y0_n = 0.7490
yf_n = 0.7540

#Plotando o quadrado do zoom
plt.plot([x0_n,xf_n], [y0_n,y0_n], 'k--')
plt.plot([x0_n,xf_n], [yf_n,yf_n], 'k--')
plt.plot([x0_n,x0_n], [y0_n,yf_n], 'k--')
plt.plot([xf_n,xf_n], [y0_n,yf_n], 'k--')

#plotando o gráfico com zoom 
a = plt.axes([0.35, 0.5, 0.15, 0.15]) #Configuração de posição e tamanho do zoom
plt.xlim(x0_n,xf_n) 
plt.ylim(y0_n,yf_n)
plt.grid()
plt.plot(Ta,referencia,'--b',label='referencia')
plt.plot(Ta, y_degrau_positivo + pt_operacao_rad, 'black', label = 'Modelo - 3P') #3P
plt.plot(Ta, y5, 'r', label='Modelo - 5ª ordem de Padé') #atraso 5ª ordem}
plt.show()

#Projeto de controlador P e PI pelo método de ziegler

theta_atraso = 0.15 #Atraso
num, dem = ct.pade(theta_atraso, 5) # Definindo as condiçoes para a função pade
Gd7 = ct.tf(num, dem) # Criando função de tranferencia
Atraso = ct.tf2io(Gd7,name ='atraso',inputs='u', outputs='y')

#Sintonização do controlador proporcional
Kc_P_ZN = tau / (K1 * theta_atraso) 
numerador_P_ZN = [Kc_P_ZN] 
denominador_P_ZN = 1
G_P_ZN = ct.tf(numerador_P_ZN,denominador_P_ZN)

#Sintonização controlador proporcional integral
Kc_PI_ZN = (0.9*tau)/(K1*theta_atraso) #Kc tabela 8.7
T_PI_ZN=(theta_atraso*10)/3 #Ti tabela 8.6
numerador_PI_ZN = [Kc_PI_ZN * T_PI_ZN, Kc_PI_ZN]
denominador_PI_ZN = [T_PI_ZN, 0]
G_PI_ZN = ct.tf(numerador_PI_ZN,denominador_PI_ZN)

#Método pela curva de CHR
#Sintonização do controlador Proporcional pela curva CHR
Kc_P_CHR = (0.3*tau) / (K1 * theta_atraso)
numerador_P_CHR = [Kc_P_CHR] 
denominador_P_CHR = 1
G_P_CHR = ct.tf(numerador_P_CHR,denominador_P_CHR)

#Sintonização do controlador proporcional integral pela curva CHR
Kc_PI_CHR = (0.6*tau)/(K1*theta_atraso) #Kc tabela 8.7
T_PI_CHR= 4*theta_atraso #Ti tabela 8.6
numerador_PI_CHR = [Kc_PI_CHR * T_PI_CHR, Kc_PI_CHR]
denominador_PI_CHR = [T_PI_CHR, 0]
G_PI_CHR = ct.tf(numerador_PI_CHR,denominador_PI_CHR)

time, y_degrau_G_P_ZN = ct.forced_response(G_P_ZN,T= time, U = Degraus_U_conc - U_ent) # Simulando a saída de um sistema linear
print(f"Controlador P -  Ziegler Nichols: {G_P_ZN}") #Imprimindo a Função de tranferencia do controlador P
print(f"Controlador PI - Ziegler Nichols: {G_PI_ZN}")#Imprimindo a Função de tranferencia do controlador PI
print(f"Controlador P -  CHR: {G_P_CHR}")#Imprimindo a Função de tranferencia do controlador P
print(f"Controlador PI - CHR: {G_PI_CHR}")#Imprimindo a Função de tranferencia do controlador PI

#Seguimento de referencia
tmp_final_2 = 160 #Tempo de Simulação 
T = 0.001#Time step 
time = np.arange(0, tmp_final_2, T) #Cria vetor espaçado do valor T indo de 0 à tmp_final_2 
u_massa = np.zeros(len(time)) #Cria um vetor de 0
U0 = U_ent * np.ones(time.shape) #Vetor de valores de entrada
referencia_3 = np.zeros(len(time))#Vetor de referencia

controlador_P_ZN = ct.tf2io(G_P_ZN, name='controlador_P_ZN', inputs='u', outputs='y') #controlador P de Ziegler
controlador_PI_ZN = ct.tf2io(G_PI_ZN, name='controlador_PI_ZN', inputs='u', outputs='y') #controlador PI de Ziegler
controlador_P_CHR = ct.tf2io(G_P_CHR, name='controlador_P_CHR', inputs='u', outputs='y') #controlador P do CHR
controlador_PI_CHR = ct.tf2io(G_PI_CHR, name='controlador_PI_CHR', inputs='u', outputs='y') #controlador PI do CHR


#Definindo os degraus do sistema para malha fechada
for k in range(len(time)):
    if time[k] < 15: # Ponto de equilíbrio
        referencia_3[k] = pt_operacao_rad
    elif time[k] < 30: # Degrau unitário superior
        referencia_3[k] = pt_operacao_rad + np.radians(3.45)
    elif time[k] < 45: # Condição de equilíbrio
        referencia_3[k] = pt_operacao_rad
    elif time[k] < 60: # Degrau unitário inferior
        referencia_3[k] = pt_operacao_rad - np.radians(3.45)
    elif time[k] < 75: # Metade do degrau unitário superior
        referencia_3[k] = pt_operacao_rad + np.radians(1.725)
    elif time[k] < 90: # Degrau unitário superior
        referencia_3[k] = pt_operacao_rad + np.radians(3.45)
    elif time[k] < 105: # Metade do degrau unitário inferior
        referencia_3[k] = pt_operacao_rad - np.radians(1.725)
    else: # Condição de equilíbrio
        referencia_3[k] = pt_operacao_rad

#controlador P - método da Curva de Reação - Ziegler-Nichols
FanPlate_P_ZN = ct.InterconnectedSystem(
    (controlador_P_ZN, Atraso, FanPlate), name='FanPlate_P_ZN',
    connections = (('controlador_P_ZN.u', '-FanPlate.y'),
                   ('atraso.u','controlador_P_ZN.y'),
                   ('FanPlate.u','atraso.y')),
    inplist = ('controlador_P_ZN.u','atraso.u', 'FanPlate.u_massa'),
    inputs = ('theta_referencia_P_ZN', 'U0', 'u_massa'),
    outlist = ('FanPlate.y', 'atraso.u'),
    outputs = ('y','u'))

#controladorador PI - método da Curva de Reação - Ziegler-Nichols
FanPlate_PI_ZN = ct.InterconnectedSystem(
    (controlador_PI_ZN, Atraso, FanPlate), name='FanPlate_PI_ZN',
    connections = (('controlador_PI_ZN.u', '-FanPlate.y'),
                   ('atraso.u','controlador_PI_ZN.y'),
                   ('FanPlate.u','atraso.y')),
    inplist = ('controlador_PI_ZN.u','atraso.u', 'FanPlate.u_massa'),
    inputs = ('theta_referencia_PI_ZN', 'U0', 'u_massa'),
    outlist = ('FanPlate.y', 'atraso.u'),
    outputs = ('y','u'))

#controladorador P - método CHR
FanPlate_P_CHR = ct.InterconnectedSystem(
    (controlador_P_CHR, Atraso, FanPlate), name='FanPlate_P_CHR',
    connections = (('controlador_P_CHR.u', '-FanPlate.y'),
                   ('atraso.u','controlador_P_CHR.y'),
                   ('FanPlate.u','atraso.y')),
    inplist = ('controlador_P_CHR.u','atraso.u', 'FanPlate.u_massa'),
    inputs = ('theta_referencia_P_CHR', 'U0', 'u_massa'),
    outlist = ('FanPlate.y', 'atraso.u'),
    outputs = ('y','u'))

#controladorador PI - método CHR
FanPlate_PI_CHR = ct.InterconnectedSystem(
    (controlador_PI_CHR, Atraso, FanPlate), name='FanPlate_PI_CHR',
    connections = (('controlador_PI_CHR.u', '-FanPlate.y'),
                   ('atraso.u','controlador_PI_CHR.y'),
                   ('FanPlate.u','atraso.y')),
    inplist = ('controlador_PI_CHR.u','atraso.u', 'FanPlate.u_massa'),
    inputs = ('theta_referencia_PI_CHR', 'U0', 'u_massa'),
    outlist = ('FanPlate.y', 'atraso.u'),
    outputs = ('y','u'))

pos_incial = np.radians(40) #Condição inicial do sistema
vel_inicial = 0 #Condição inicial do sistema
u_massa = np.zeros(len(time))#Vetor de entrada

t_P_PI, y_P_ZN = ct.input_output_response(FanPlate_P_ZN, time, [referencia_3, U0, u_massa], [0,0,0,0,0,0.9*pos_incial, vel_inicial])      #controlador P - método de Ziegler
t_P_PI, y_PI_ZN = ct.input_output_response(FanPlate_PI_ZN, time, [referencia_3, U0, u_massa], [0,0,0,0,0,0,0.9*pos_incial, vel_inicial])  #controlador PI do método de Ziegler
t_P_PI, y_P_CHR = ct.input_output_response(FanPlate_P_CHR, time, [referencia_3, U0, u_massa], [0,0,0,0,0,0.9*pos_incial, vel_inicial])    #controlador P - método CHR
t_P_PI, y_PI_CHR = ct.input_output_response(FanPlate_PI_CHR, time, [referencia_3, U0, u_massa], [0,0,0,0,0,0,0.9*pos_incial, vel_inicial])#controlador PI - método CHR


#Plot dos controladores 
plt.figure(6)
plt.subplot(2,1,1)
plt.plot(time, referencia_3,'--b', label='referencia')
plt.plot(time, y_P_ZN[0],'k', label= 'c1')  #controlador P - método da Curva de Reação - Ziegler-Nichols
plt.plot(time, y_PI_ZN[0],'g', label='c2')  #controlador PI - método da Curva de Reação - Ziegler-Nichols
plt.plot(time, y_P_CHR[0],'r', label='c3')  #controlador P - método CHR
plt.plot(time, y_PI_CHR[0],'purple', label='c4') #controlador PI - método CHR
plt.ylabel('$\\theta$[rad]')
plt.xlim(0,max(time))
plt.ylim(0.6, 0.88)
plt.legend(loc = 'upper right')
plt.grid()

#Definindo limites do plot do Zoom 
x0_n = 13.000
xf_n = 19.000
y0_n = 0.7425
yf_n = 0.86

#Plot da linha tracejada 
plt.plot([x0_n,xf_n], [y0_n,y0_n], 'k-')
plt.plot([x0_n,xf_n], [yf_n,yf_n], 'k--')
plt.plot([x0_n,x0_n], [y0_n,yf_n], 'k--')
plt.plot([xf_n,xf_n], [y0_n,yf_n], 'k--')

#Plot com zoom na curva
a = plt.axes([0.7, 0.57, 0.15, 0.13]) #Define a posição do grafico
plt.xlim(x0_n,xf_n)
plt.ylim(y0_n,yf_n)
plt.grid()
plt.plot(time, referencia_3,'--b', label='referencia')
plt.plot(time, y_P_ZN[0],'k', label='c1')         #controlador P - método da Curva de Ziegler
plt.plot(time, y_PI_ZN[0],'g', label='c2')        #controlador PI - método da Curva de Ziegler
plt.plot(time, y_P_CHR[0],'r', label='c3')        #controlador P - método CHR
plt.plot(time, y_PI_CHR[0],'purple', label='c4')  #controlador PI - método CHR
plt.subplot(2,1,2)
plt.plot(time,y_P_ZN[1,:],'b', label='u')                #controlador P - método da Curva de Ziegler
plt.plot(time, y_PI_ZN[1,:],'r', label='theta(t)')       #controlador PI - método da Curva de Ziegler
plt.plot(time, y_P_CHR[1,:],'g', label='theta(t)')       #controlador P - método CHR
plt.plot(time, y_PI_CHR[1,:],'purple', label='theta(t)') #controlador PI - método CHR
plt.ylabel('referencia')
plt.xlabel('Tempo [s]')
plt.xlim(0,max(time))
plt.legend(loc='upper right')
plt.grid()
plt.show()

#Alteração parametrica da massa
for k in range(len(time)): # Aplica a variação na massa do sistema entre 20s e 120s
    if time[k] >= 20 and time[k] < 120:
        u_massa[k] = 0.2

referencia_4 = pt_operacao_rad * np.ones(time.shape) # Vetor de referencia do ponto de operação
t1_, y1_ = ct.input_output_response(FanPlate_P_ZN, time, [referencia_4,U0,u_massa], [0,0,0,0,0,0.9*pos_incial,vel_inicial])     #controlador P - método da Curva de Ziegler
t2_, y2_ = ct.input_output_response(FanPlate_PI_ZN, time, [referencia_4,U0,u_massa], [0,0,0,0,0,0,0.9*pos_incial,vel_inicial])    #controlador PI - método da Curva Ziegler
t3_, y3_ = ct.input_output_response(FanPlate_P_CHR, time, [referencia_4,U0,u_massa], [0,0,0,0,0,0.9*pos_incial,vel_inicial])    #controlador P - método da Curva de CHR
t4_, y4_ = ct.input_output_response(FanPlate_PI_CHR, time, [referencia_4,U0,u_massa], [0,0,0,0,0,0,0.9*pos_incial,vel_inicial])   #controlador PI - método da Curva de CHR

#Plotando os Resultados dos controladores 
plt.figure(7)
plt.subplot(2,1,1)
plt.plot(t1_,referencia_4,'--b',label='referencia')
plt.plot(t1_,y1_[0,:],'r',label='c1') #controlador P - método da Curva de Ziegler
plt.plot(t2_,y2_[0,:],'k',label='c2') #controlador PI - método da Curva Ziegler
plt.plot(t3_,y3_[0,:],"c",label='c3')#controlador P - método da Curva de CHR
plt.plot(t4_,y4_[0,:],"g",label='c4')#controlador PI - método da Curva de CHR
plt.ylabel('$\\theta(t) [rad]$')
plt.ylim(0.57,0.87)
plt.xlim(0,tmp_final_2)
plt.legend(loc = 'lower right')
plt.grid()

#Definindo os limites do zoom
x0_n = 18
xf_n = 30
y0_n = 0.72
yf_n = 0.76

#Plotando a linha pontilhada 
plt.plot([x0_n,xf_n], [y0_n,y0_n], 'k--')
plt.plot([x0_n,xf_n], [yf_n,yf_n], 'k--')
plt.plot([x0_n,x0_n], [y0_n,yf_n], 'k--')
plt.plot([xf_n,xf_n], [y0_n,yf_n], 'k--')


a = plt.axes([0.17, 0.57, 0.15, 0.10]) # x-position, y-position, width, height
plt.xlim(x0_n,xf_n)
plt.ylim(y0_n,yf_n)
plt.grid()
plt.plot(t1_,referencia_4,'--b',label='referencia')
plt.plot(t1_,y1_[0,:],'r',label='c1')
plt.plot(t2_,y2_[0,:],'k',label='c2')
plt.plot(t3_,y3_[0,:],"c",label='c3')
plt.plot(t4_,y4_[0,:],"g",label='c4')

plt.subplot(2,1,2)
plt.plot(t1_,y1_[1,:],'r',label='c1')
plt.plot(t2_,y2_[1,:],'k',label='c2')
plt.plot(t3_,y3_[1,:],"c",label='c3')
plt.plot(t4_,y4_[1,:],"g",label='c4')
plt.ylabel('u(t) [m²/s²]')
plt.xlabel('Tempo [s]')
plt.ylim(0.5, 6.2)
plt.xlim(-0.5, tmp_final_2)
plt.legend()
plt.legend(prop={'size':10})
plt.grid()
plt.show()