#Laboratório de Análise de Sistemas Lineares
#@Autores: Thiago Henrique de Faria Costa e Mateus Alves de Sales
#@Data: 18/01/2022

import  numpy as np                      #\    
from    matplotlib import pyplot as plt   #\    
import  control as ct                       #- Importando as bibliotecas necessárias para desenvolvimento do código
import  statistics as sta                 #/ 
import  math                             #/

plt.close('all') #Fechando todas as abas de gráficos aberta

def model_update(t,x,u_array,params):
    u = u_array[0]
    u_massa = u_array[1]
    MT = m_t*(1 + u_massa)

    K_1 = (d_cm*rho_ar*C_a*L_a*L_1)/(2*MT*(((L_t**2)/12)+d_cm**2))
    K_2 = (g*d_cm)/ (((L_t**2)/12)+d_cm**2)
    K_3 = (atrito*(d_cm**2))/(MT*(((L_t**2)/12)+d_cm**2))

    # Variáveis de estado
    x1 = x[0] # Posição angular    \
    x2 = x[1] # Velocidade angular /  # Variáveis de estado
    return np.array([x2, K_1*(np.cos(x1)**2)*u - (K_2*np.sin(x1) + K_3*x2)], dtype=object) #Retorna a equação diferencial na forma de vertor


def model_output(t,x,u,params): # A saí­da do sistema (y) é o estado X[0], logo, y = theta.
    return x[0]


def K_ganho(massa):
    K_1 = (d_cm*rho_ar*C_a*L_a*L_1)/(2*massa*(((L_t**2)/12)+d_cm**2))
    K_2 = (g*d_cm)/ (((L_t**2)/12)+d_cm**2)
    K_3 = (atrito*(d_cm**2))/(massa*(((L_t**2)/12)+d_cm**2))

    return K_1, K_2, K_3

def coeficiente_angular(valor_y, valor_x, j): #Função que retorna o coeficiente angular da reta tangente a curva de Ziegler-Nichols
    return (valor_y[j]-valor_y[j-1])/(valor_x[j]-valor_x[j-1])

def valor_K(vetor1, vetor2, ponto_de_operacao, variavel_de_controle, k): #Função que retorna o valor do ganho 
    return vetor1[k]-ponto_de_operacao, vetor2[k]-variavel_de_controle, (vetor1[k]-ponto_de_operacao)/(vetor2[k]-variavel_de_controle)

def equacao(vetor1, vetor2, c_angular, k): #Função que retorna a equação da reta tangente a curva de Ziegler-Nichols
    return c_angular*vetor1 - c_angular*vetor1[k] + vetor2[k]


""" Definindo as varíaveis do sistema"""


L_a = 0.154                                                     # Largura da placa móvel de alumínio
L_1 = 0.155                                                     # Comprimento da placa abaixo do eixo de rotação.
L_t = 0.270                                                     # Comprimento da placa abaixo do eixo de rotação.
d_cm = 0.020                                                    # Distância do centro de massa da placa
rho_ar = 1.23                                                   # Densidade do ar
m_t = 0.005                                                     # Massa total da placa
C_a = 2.05                                                      # Coeficiente de arrasto
atrito = 5                                                      # Coeficiente de atrito viscoso
g = 9.81                                                        # Gravidade
pt_operacao = 43                                     # Ponto de operação
pt_operacao_rad = math.radians(pt_operacao)            # Ponto de operação em radianos
pt_inicial = 0                                               # Ponto inicial

"""Encontrando o valor das constantes K do sistema"""

K_1, K_2, K_3 = K_ganho(m_t)
print(f'K1: ', K_1)
print(f'K2: ', K_2)
print(f'K3: ', K_3)


tmp_final = 60   # Tempo final da análise (s)
periodo = 0.001    # Período de análise (s)
time = np.arange(0, tmp_final, periodo) # Criando um lista que vai de 0 até 60, com período 1
U_ent = (K_2 * math.sin(pt_operacao_rad))/(K_1 * (math.cos(pt_operacao_rad)**2))
U0 = U_ent * np.ones(time.shape) # Cria uma lista do tamanho da lista de tempo com todos os valores iguais ao valor de entrada
ref = pt_operacao_rad * np.ones(time.shape) # Cria uma lista do tamanho da lista de tempo com todos os valores iguais ao ponto de operação; é a linha de referência utilizada nos gráficos
u_massa = np.zeros(len(time))

print(f'Sinal de entrada:',U_ent)

""" Criando os degraus no sistema"""

Degraus_U = U_ent * np.ones(len(time))          # Criando um novo vetor para alocar os degraus
Degraus_U = np.array_split(Degraus_U, 4)            # Dividindo o vetor em 4 partes iguais

Degraus_U[0][:]=U0[0]                               # Definindo a primeira parte do vetor como U_0
Degraus_U[1][:]= 1.20*U0[0]                         # Definindo a segunda parte do vetor como +20% de U_0
Degraus_U[2][:]= U0[0]                              # Definindo a terceira parte do vetor como U_0   
Degraus_U[3][:]= 0.80*U0[0]                         # Definindo a segunda parte do vetor como +20% de U_0  


Degraus_U_conc = np.concatenate([Degraus_U[0], Degraus_U[1], Degraus_U[2], Degraus_U[3]]) # Unindo os 4 vetores que antes foram divididos

FanPlate = ct.NonlinearIOSystem(model_update, model_output , states=2, name='FanPlate', inputs = ('u', 'u_massa'), outputs = ('y'))

""" Simulação do sistema sem aplicação de degraus na entrada:"""

t, y = ct.input_output_response(FanPlate, time, [U0, u_massa], pt_inicial)

""" Plotando sistema sem aplicação de degraus na entrada:"""

plt.figure(1)
plt.subplot(2,1,1)
plt.plot(t,ref,'--k', label='ref')
plt.plot(t,y[0],'b', label='$\\theta$(t)')
plt.ylabel('$\\theta$(t)[rad]')
plt.xlim(0,10)
plt.ylim(0,1.5)
plt.legend()
# plt.title('Resposta temporal do sistema em malha aberta sem degrau')
plt.grid()
plt.subplot(2,1,2)
plt.plot(time,U0,'b', label='Q (t)')
plt.ylabel('u(t)[m^2/s^2]')
plt.xlabel('Tempo[s]')
plt.legend()
plt.xlim(0,tmp_final)
plt.ylim(1,3)
plt.grid()
plt.show()

""" Plotando sistema com aplicação de degraus na entrada:"""

t_degrau, y_degrau = ct.input_output_response(FanPlate, time, [Degraus_U_conc, u_massa], pt_inicial)
y_max_deg = y_degrau[0][29499]
plt.figure(2)
plt.subplot(2,1,1)
plt.plot(t_degrau,ref,'--k', label='ref')
plt.plot(t_degrau,y_degrau[0],'b', label='T(t)')
plt.ylabel('$\\theta$(t)[rad]')
plt.xlim(0,tmp_final)
plt.ylim(0,1.5)
plt.legend()
# plt.title('Resposta temporal do sistema em malha aberta com degrau')
plt.grid()
plt.subplot(2,1,2)
plt.plot(time,Degraus_U_conc,'b', label='Q (t)')
plt.ylabel('u(t)[m^2/s^2]')
plt.xlabel('Tempo[s]')
plt.legend()
plt.xlim(0,tmp_final)
plt.ylim(0,4)
plt.grid()
plt.show()

"""Determinando os parâmetros de Ziegler-Nichols" "Item 02"""

s = ct.TransferFunction.s # Cria um sistema de função de transferência.

"""Reta tangente ao degrau positivo"""

amp_degrau_pos = y_degrau[0][29999]                                           # Definindo a amplitude do degrau positivo
ca_degrau_pos = coeficiente_angular(y_degrau[0], time, 15100)                 # Definindo o coeficiente angular pela função
equacao_degrau_positivo = equacao(time, y_degrau[0], ca_degrau_pos, 15100)    # Definindo a equação da reta
ref_2 = amp_degrau_pos * np.ones(time.shape)                                   # Cria uma lista do tamanho da lista de tempo, que possui todos os valores iguais ao ponto de operação

"""Determinando os parâmetros de Ziegler através da análise gráfica para o degrau positivo, o valor de theta é de 0.011, o que pode ser considerado aproximadamente nulo:"""

tau = 16.02 - 15.01
delta_theta, delta_u, K1 = valor_K(y_degrau[0], Degraus_U_conc, pt_operacao_rad, U_ent, 29999)
print(f'delta y: ',delta_theta)
print(f'delta u: ',delta_u)
print(f'K1: ', K1)

"""Modelo de primeira ordem - Ziegler-Nichols - degrau positivo""" 

modelo_degrau_positivo = K1/(tau*s+1)  # Modelo para degrau positivo

t_degrau_positivo, y_degrau_positivo = ct.forced_response(modelo_degrau_positivo,T= time, U = Degraus_U_conc - U_ent) # Simulando a saída do sistema linear

"""Gráficos método de Ziegler-Nichols - degrau positivo""" 

plt.figure(3)
# plt.title("Gráfico Método de Ziegler-Nichols - degrau positivo")
plt.plot(time,y_degrau[0],'b',label='$\\theta$')
plt.plot(time,y_degrau_positivo + pt_operacao_rad,'g',label='$\\theta_{Ziegler-Nichols}$')
plt.plot(time,equacao_degrau_positivo,'r',label='')
plt.plot(time,ref,'--k',label='ref')
plt.plot(time,ref_2,'--c',label='ref')
plt.grid()
plt.ylabel('$\\theta$(t)[rad]')
plt.xlabel('Tempo[s]')
plt.legend()
plt.ylim(0.71,0.82)
plt.xlim(14.5,22)

x0 = 15.000
xf = 15.020
y0 = 0.7504
yf = 0.7508

# PLOT NORMAL

plt.plot([x0,xf], [y0,y0], 'c--')
plt.plot([x0,xf], [yf,yf], 'c--')
plt.plot([x0,x0], [y0,yf], 'c--')
plt.plot([xf,xf], [y0,yf], 'c--')

# PLOT COM ZOOM

a = plt.axes([0.35, 0.5, 0.15, 0.15]) # x-position, y-position, width, height
plt.xlim(x0,xf)
plt.ylim(y0,yf)
plt.grid()
plt.plot(time,y_degrau[0],'b',label='$g(s)$')
plt.plot(time,equacao_degrau_positivo,'r',label='ref degrau')
plt.plot(time,ref,'--k',label='ref')
plt.show()

"""Obtendo aproximações por Padé""" "Item 05"
Gi = K1/(tau*s+1)

#1ª ordem
n1 = 1 
N1 = ct.pade(0.15,n1) #aproximação de padé
Gd1 = ct.TransferFunction(np.array(N1[0]),np.array(N1[1])) #construção da função de transferência.
Gr1 = Gi*Gd1
t1, y1 = ct.forced_response(Gr1,T= time, U = U0-U_ent) 
y1 = y1 + pt_operacao_rad

#3ª ordem
n3 = 3 
N3 = ct.pade(0.15,n3) #aproximação de padé
Gd3 = ct.TransferFunction(np.array(N3[0]),np.array(N3[1])) #construção da função de transferência.
Gr3 = Gi*Gd3
t3, y3 = ct.forced_response(Gr1,T= time, U = U0-U_ent) 
y3 = y3 + pt_operacao_rad

#5ª ordem
n5 = 5 
N5 = ct.pade(0.15,n5) #aproximação de padé
Gd5 = ct.TransferFunction(np.array(N5[0]),np.array(N5[1])) #construção da função de transferência.
Gr5 = Gi*Gd5
t5, y5 = ct.forced_response(Gr1,T= time, U = U0-U_ent) 
y5 = y5 + pt_operacao_rad

#9ª ordem
n9 = 9 
N9 = ct.pade(0.15,n9) #aproximação de padé
Gd9 = ct.TransferFunction(np.array(N9[0]),np.array(N9[1])) #construção da função de transferência.
Gr9 = Gi*Gd9
t9, y9 = ct.forced_response(Gr1,T= time, U = U0-U_ent) 
y9 = y9 + pt_operacao_rad

#conversão das funções de transferência do atraso para um sistema em espaço de estado
atraso_1 = ct.tf2io(Gd1, name='atraso_1', inputs='u', outputs='y1')#1ª ordem
atraso_3 = ct.tf2io(Gd3, name='atraso_3', inputs='u', outputs='y3')#3ª ordem
atraso_5 = ct.tf2io(Gd5, name='atraso_5', inputs='u', outputs='y5')#5ª ordem
atraso_9 = ct.tf2io(Gd9, name='atraso_9', inputs='u', outputs='y9')#9ª ordem


#construção da malha aberta

#1ª ordem
malha_aberta_01 = ct.InterconnectedSystem(
    (atraso_1, FanPlate), name='malha_aberta_01',
    connections = (('atraso_1.u',),('FanPlate.u','atraso_1.y1')),
    inplist = ('atraso_1.u'),
    outlist = ('FanPlate.y'))

X1 = np.zeros(n1+2) #condições iniciais nulas
Ta, y1 = ct.input_output_response(malha_aberta_01, time, Degraus_U_conc, X1)


#3ª ordem
malha_aberta_03 = ct.InterconnectedSystem(
    (atraso_3, FanPlate), name='malha_aberta_03',
    connections = (('atraso_3.u',),('FanPlate.u','atraso_3.y3')),
    inplist = ('atraso_3.u'),
    outlist = ('FanPlate.y'))

X3 = np.zeros(n3+2) #condições iniciais nulas
Ta, y3 = ct.input_output_response(malha_aberta_03, time, Degraus_U_conc, X3)


#5ª ordem
malha_aberta_05 = ct.InterconnectedSystem(
    (atraso_5, FanPlate), name='malha_aberta_05',
    connections = (('atraso_5.u',),('FanPlate.u','atraso_5.y5')),
    inplist = ('atraso_5.u'),
    outlist = ('FanPlate.y'))

X5 = np.zeros(n5+2) #condições iniciais nulas
Ta, y5 = ct.input_output_response(malha_aberta_05, time, Degraus_U_conc, X5)


#9ª ordem
malha_aberta_09 = ct.InterconnectedSystem(
    (atraso_9, FanPlate), name='malha_aberta_09',
    connections = (('atraso_9.u',),('FanPlate.u','atraso_9.y9')),
    inplist = ('atraso_9.u'),
    outlist = ('FanPlate.y'))

X9 = np.zeros(n9+2) #condições iniciais nulas
Ta, y9 = ct.input_output_response(malha_aberta_09, time, Degraus_U_conc, X9)

"""plotagem da resposta temporal do sistema com atraso - Item 06"""

plt.figure(4)
plt.plot(Ta,ref,'--k',label='ref')
plt.plot(Ta, y_degrau_positivo + pt_operacao_rad, 'black', label = 'Modelo - 3P') #3P
plt.plot(Ta, y1, 'purple', label='Modelo - 1ª ordem') #atraso 1 
plt.plot(Ta, y3, 'g', label='Modelo - 2ª ordem') #atraso 3 
plt.plot(Ta, y5, 'r', label='Modelo - 5ª ordem') #atraso 5 
plt.plot(Ta, y9, 'b', label='Modelo - 9ª ordem') #atraso 9 
plt.legend()
plt.ylabel('$\\theta[rad]$')
plt.xlabel('Tempo[s]')
plt.ylim(0.70,0.81)
plt.xlim(14.5,21.5)
plt.grid()

x0_n = 14.95
xf_n = 15.20
y0_n = 0.7490
yf_n = 0.7510

# PLOT NORMAL

plt.plot([x0_n,xf_n], [y0_n,y0_n], 'c--')
plt.plot([x0_n,xf_n], [yf_n,yf_n], 'c--')
plt.plot([x0_n,x0_n], [y0_n,yf_n], 'c--')
plt.plot([xf_n,xf_n], [y0_n,yf_n], 'c--')

# PLOT COM ZOOM

a = plt.axes([0.35, 0.5, 0.15, 0.15]) # x-position, y-position, width, height
plt.xlim(x0_n,xf_n)
plt.ylim(y0_n,yf_n)
plt.grid()
plt.plot(Ta,ref,'--k',label='ref')
plt.plot(Ta, y_degrau_positivo + pt_operacao_rad, 'black', label = 'Modelo - 3P') #3P
plt.plot(Ta, y1, 'purple', label='Modelo - 1ª ordem') #atraso 1 
plt.plot(Ta, y3, 'g', label='Modelo - 2ª ordem') #atraso 3 
plt.plot(Ta, y5, 'r', label='Modelo - 5ª ordem') #atraso 5 
plt.plot(Ta, y9, 'b', label='Modelo - 9ª ordem') #atraso 9 
plt.show()

"""Comparando o degrau positivo de Ziegler com a 5ª aproximação de Padé""" "Item 07"

plt.figure(5)
plt.plot(Ta,ref,'--k',label='ref')
plt.plot(Ta, y_degrau_positivo + pt_operacao_rad, 'black', label = 'Modelo - 3P') #3P - ZN
plt.plot(Ta, y5, 'r', label='Modelo - 5ª ordem') #Modelo - 5ª ordem 
plt.legend()
plt.ylabel('theta[rad]')
plt.xlabel('Tempo[s]')
plt.xlim(13.8, 26)
plt.ylim(0.675, 0.850)
plt.grid()

x0_n = 14.950
xf_n = 15.200
y0_n = 0.7490
yf_n = 0.7540

# PLOT NORMAL

plt.plot([x0_n,xf_n], [y0_n,y0_n], 'c--')
plt.plot([x0_n,xf_n], [yf_n,yf_n], 'c--')
plt.plot([x0_n,x0_n], [y0_n,yf_n], 'c--')
plt.plot([xf_n,xf_n], [y0_n,yf_n], 'c--')

# PLOT COM ZOOM

a = plt.axes([0.35, 0.5, 0.15, 0.15]) # x-position, y-position, width, height
plt.xlim(x0_n,xf_n)
plt.ylim(y0_n,yf_n)
plt.grid()
plt.plot(Ta,ref,'--k',label='ref')
plt.plot(Ta, y_degrau_positivo + pt_operacao_rad, 'black', label = 'Modelo - 3P') #3P
plt.plot(Ta, y5, 'r', label='Modelo - 5ª ordem') #atraso 5
plt.show()


"""Projetando os controladores P e PI - Item 08"""

"""Método da curva de Ziegler - Nichols"""

theta_atraso = 0.15

num, dem = ct.pade(theta_atraso, 5) # Definindo as aproximações do atraso por pade de 5ª ordem
Gd7 = ct.tf(num, dem) # Função de transferência da aproximação do atraso de 5ª ordem
Atraso = ct.tf2io(Gd7,name ='atraso',inputs='u', outputs='y')

#Sintonização e determinação dos valores do controlador Proporcional
Kc_P_ZN = tau / (K1 * theta_atraso) 
numerador_P_ZN = [Kc_P_ZN] 
denominador_P_ZN = 1
G_P_ZN = ct.tf(numerador_P_ZN,denominador_P_ZN)

#Sintonização e determinação dos valores do controlador Proporcional Integral
Kc_PI_ZN = (0.9*tau)/(K1*theta_atraso) #Kc tabela 8.7
T_PI_ZN=(theta_atraso*10)/3 #Ti tabela 8.6
numerador_PI_ZN = [Kc_PI_ZN * T_PI_ZN, Kc_PI_ZN]
denominador_PI_ZN = [T_PI_ZN, 0]
G_PI_ZN = ct.tf(numerador_PI_ZN,denominador_PI_ZN)

"""Método da curva de CHR"""
#Sintonização e determinação dos valores do controlador Proporcional
Kc_P_CHR = (0.3*tau) / (K1 * theta_atraso)
numerador_P_CHR = [Kc_P_CHR] 
denominador_P_CHR = 1
G_P_CHR = ct.tf(numerador_P_CHR,denominador_P_CHR)

#Sintonização e determinação dos valores do controlador Proporcional Integral
Kc_PI_CHR = (0.6*tau)/(K1*theta_atraso) #Kc tabela 8.7
T_PI_CHR= 4*theta_atraso #Ti tabela 8.6
numerador_PI_CHR = [Kc_PI_CHR * T_PI_CHR, Kc_PI_CHR]
denominador_PI_CHR = [T_PI_CHR, 0]
G_PI_CHR = ct.tf(numerador_PI_CHR,denominador_PI_CHR)

time, y_degrau_G_P_ZN = ct.forced_response(G_P_ZN,T= time, U = Degraus_U_conc - U_ent) # Simulando a saída de um sistema linear
print(f"Controlador P -  Ziegler Nichols: {G_P_ZN}")
print(f"Controlador PI - Ziegler Nichols: {G_PI_ZN}")
print(f"Controlador P -  CHR: {G_P_CHR}")
print(f"Controlador PI - CHR: {G_PI_CHR}")

"""Seguimento de referência - Item 10-a"""

tmp_final_2 = 160 # Duração da simulação (s)
T = 0.001         # Período de amostragem (s)
time = np.arange(0, tmp_final_2, T) # Vetor do tempo de simulação espaçados de 0.01s
u_massa = np.zeros(len(time))
U0 = U_ent * np.ones(time.shape) # Cria uma lista do tamanho da lista de tempo com todos os valores iguais ao valor de entrada
ref_3 = np.zeros(len(time))

controlador_P_ZN = ct.tf2io(G_P_ZN, name='controlador_P_ZN', inputs='u', outputs='y') #controlador P da Curva de Reação Ziegler-Nichols
controlador_PI_ZN = ct.tf2io(G_PI_ZN, name='controlador_PI_ZN', inputs='u', outputs='y') #controlador PI da Curva de Reação Ziegler-Nichols
controlador_P_CHR = ct.tf2io(G_P_CHR, name='controlador_P_CHR', inputs='u', outputs='y') #controlador P do CHR
controlador_PI_CHR = ct.tf2io(G_PI_CHR, name='controlador_PI_CHR', inputs='u', outputs='y') #controlador PI do CHR

# saturacao = ct.NonlinearIOSystem(0, saida_saturacao, states=0, name='saturacao', inputs=('u'), outputs=('y5'))

"""Definindo os degraus do sistema para malha fechada"""

for k in range(len(time)):
    if time[k] < 15: # Estabilizando no ponto de equilíbrio
        ref_3[k] = pt_operacao_rad
    elif time[k] < 30: # Degrau unitário superior
        ref_3[k] = pt_operacao_rad + np.radians(3.45)
    elif time[k] < 45: # Condição de equilíbrio
        ref_3[k] = pt_operacao_rad
    elif time[k] < 60: # Degrau unitário inferior
        ref_3[k] = pt_operacao_rad - np.radians(3.45)
    elif time[k] < 75: # Metade do degrau unitário superior
        ref_3[k] = pt_operacao_rad + np.radians(1.725)
    elif time[k] < 90: # Degrau unitário superior
        ref_3[k] = pt_operacao_rad + np.radians(3.45)
    elif time[k] < 105: # Metade do degrau unitário inferior
        ref_3[k] = pt_operacao_rad - np.radians(1.725)
    else: # Condição de equilíbrio
        ref_3[k] = pt_operacao_rad

#controladorador P - método da Curva de Reação - Ziegler-Nichols
FanPlate_P_ZN = ct.InterconnectedSystem(
    (controlador_P_ZN, Atraso, FanPlate), name='FanPlate_P_ZN',
    connections = (('controlador_P_ZN.u', '-FanPlate.y'),
                   ('atraso.u','controlador_P_ZN.y'),
                   ('FanPlate.u','atraso.y')),
    inplist = ('controlador_P_ZN.u','atraso.u', 'FanPlate.u_massa'),
    inputs = ('theta_ref_P_ZN', 'U0', 'u_massa'),
    outlist = ('FanPlate.y', 'atraso.u'),
    outputs = ('y','u'))

#controladorador PI - método da Curva de Reação - Ziegler-Nichols
FanPlate_PI_ZN = ct.InterconnectedSystem(
    (controlador_PI_ZN, Atraso, FanPlate), name='FanPlate_PI_ZN',
    connections = (('controlador_PI_ZN.u', '-FanPlate.y'),
                   ('atraso.u','controlador_PI_ZN.y'),
                   ('FanPlate.u','atraso.y')),
    inplist = ('controlador_PI_ZN.u','atraso.u', 'FanPlate.u_massa'),
    inputs = ('theta_ref_PI_ZN', 'U0', 'u_massa'),
    outlist = ('FanPlate.y', 'atraso.u'),
    outputs = ('y','u'))

#controladorador P - método CHR
FanPlate_P_CHR = ct.InterconnectedSystem(
    (controlador_P_CHR, Atraso, FanPlate), name='FanPlate_P_CHR',
    connections = (('controlador_P_CHR.u', '-FanPlate.y'),
                   ('atraso.u','controlador_P_CHR.y'),
                   ('FanPlate.u','atraso.y')),
    inplist = ('controlador_P_CHR.u','atraso.u', 'FanPlate.u_massa'),
    inputs = ('theta_ref_P_CHR', 'U0', 'u_massa'),
    outlist = ('FanPlate.y', 'atraso.u'),
    outputs = ('y','u'))

#controladorador PI - método CHR
FanPlate_PI_CHR = ct.InterconnectedSystem(
    (controlador_PI_CHR, Atraso, FanPlate), name='FanPlate_PI_CHR',
    connections = (('controlador_PI_CHR.u', '-FanPlate.y'),
                   ('atraso.u','controlador_PI_CHR.y'),
                   ('FanPlate.u','atraso.y')),
    inplist = ('controlador_PI_CHR.u','atraso.u', 'FanPlate.u_massa'),
    inputs = ('theta_ref_PI_CHR', 'U0', 'u_massa'),
    outlist = ('FanPlate.y', 'atraso.u'),
    outputs = ('y','u'))

pos_incial = np.radians(40) #Condição inicial do sistema (posição).
vel_inicial = 0 #Condição inicial do sistema (velocidade).
u_massa = np.zeros(len(time))

t_P_PI, y_P_ZN = ct.input_output_response(FanPlate_P_ZN, time, [ref_3, U0, u_massa], [0,0,0,0,0,0.9*pos_incial, vel_inicial])      #controlador P - método da Curva de Reação - Ziegler-Nichols
t_P_PI, y_PI_ZN = ct.input_output_response(FanPlate_PI_ZN, time, [ref_3, U0, u_massa], [0,0,0,0,0,0,0.9*pos_incial, vel_inicial])    #controlador PI do método da Curva de Reação - Ziegler-Nichols
t_P_PI, y_P_CHR = ct.input_output_response(FanPlate_P_CHR, time, [ref_3, U0, u_massa], [0,0,0,0,0,0.9*pos_incial, vel_inicial])    #controlador P - método CHR
t_P_PI, y_PI_CHR = ct.input_output_response(FanPlate_PI_CHR, time, [ref_3, U0, u_massa], [0,0,0,0,0,0,0.9*pos_incial, vel_inicial])  #controlador PI - método CHR

plt.figure(6)
plt.subplot(2,1,1)
plt.plot(time, ref_3,'--k', label='ref')
plt.plot(time, y_P_ZN[0],'b', label= 'c1')  #controlador P - método da Curva de Reação - Ziegler-Nichols
plt.plot(time, y_PI_ZN[0],'r', label='c2')  #controlador PI - método da Curva de Reação - Ziegler-Nichols
plt.plot(time, y_P_CHR[0],'g', label='c3')  #controlador P - método CHR
plt.plot(time, y_PI_CHR[0],'purple', label='c4') #controlador PI - método CHR
plt.ylabel('$\\theta$[rad]')
plt.xlim(0,max(time))
plt.ylim(0.6, 0.88)
plt.legend(loc = 'upper right')
plt.grid()

x0_n = 13.000
xf_n = 19.000
y0_n = 0.7425
yf_n = 0.86

# PLOT NORMAL

plt.plot([x0_n,xf_n], [y0_n,y0_n], 'c--')
plt.plot([x0_n,xf_n], [yf_n,yf_n], 'c--')
plt.plot([x0_n,x0_n], [y0_n,yf_n], 'c--')
plt.plot([xf_n,xf_n], [y0_n,yf_n], 'c--')

# PLOT COM ZOOM

a = plt.axes([0.7, 0.57, 0.15, 0.13]) # x-position, y-position, width, height
plt.xlim(x0_n,xf_n)
plt.ylim(y0_n,yf_n)
plt.grid()
plt.plot(time, ref_3,'--k', label='ref')
plt.plot(time, y_P_ZN[0],'b', label='c1')         #controlador P - método da Curva de Reação - Ziegler-Nichols
plt.plot(time, y_PI_ZN[0],'r', label='c2')        #controlador PI - método da Curva de Reação - Ziegler-Nichols
plt.plot(time, y_P_CHR[0],'g', label='c3')        #controlador P - método CHR
plt.plot(time, y_PI_CHR[0],'purple', label='c4')  #controlador PI - método CHR
# plt.show()


plt.subplot(2,1,2)
plt.plot(time,y_P_ZN[1,:],'b', label='u')                #controlador P - método da Curva de Reação - Ziegler-Nichols
plt.plot(time, y_PI_ZN[1,:],'r', label='theta(t)')       #controlador PI - método da Curva de Reação - Ziegler-Nichols
plt.plot(time, y_P_CHR[1,:],'g', label='theta(t)')       #controlador P - método CHR
plt.plot(time, y_PI_CHR[1,:],'purple', label='theta(t)') #controlador PI - método CHR
plt.ylabel('ref')
plt.xlabel('Tempo [s]')
plt.xlim(0,max(time))
plt.legend(loc='upper right')
plt.grid()
plt.show()

"""Rejeição a pertubação - alteração paramétrica da massa - Item 10-b"""

for k in range(len(time)): # Aplica a variação na massa do sistema entre 20s e 120s
    if time[k] >= 20 and time[k] < 120:
        u_massa[k] = 0.2

ref_4 = pt_operacao_rad * np.ones(time.shape) # Cria uma lista do tamanho da lista de tempo, em que os valores são iguais ao ponto de operação

t1_, y1_ = ct.input_output_response(FanPlate_P_ZN, time, [ref_4,U0,u_massa], [0,0,0,0,0,0.9*pos_incial,vel_inicial])     #controlador P - método da Curva de Reação - Ziegler-Nichols
t2_, y2_ = ct.input_output_response(FanPlate_PI_ZN, time, [ref_4,U0,u_massa], [0,0,0,0,0,0,0.9*pos_incial,vel_inicial])    #controlador PI - método da Curva de Reação - Ziegler-Nichols
t3_, y3_ = ct.input_output_response(FanPlate_P_CHR, time, [ref_4,U0,u_massa], [0,0,0,0,0,0.9*pos_incial,vel_inicial])    #controlador P - método da Curva de CHR
t4_, y4_ = ct.input_output_response(FanPlate_PI_CHR, time, [ref_4,U0,u_massa], [0,0,0,0,0,0,0.9*pos_incial,vel_inicial])   #controlador PI - método da Curva de CHR

plt.figure(7)
plt.subplot(2,1,1)
plt.plot(t1_,ref_4,'--k',label='ref')
plt.plot(t1_,y1_[0,:],'mediumvioletred',label='c1')
plt.plot(t2_,y2_[0,:],'goldenrod',label='c2')
plt.plot(t3_,y3_[0,:],"-.b",label='c3')
plt.plot(t4_,y4_[0,:],"-.g",label='c4')
plt.ylabel('$\\theta(t) [rad]$')
plt.ylim(0.57,0.87)
plt.xlim(0,tmp_final_2)
plt.legend(loc = 'lower right')
plt.grid()
x0_n = 18
xf_n = 30
y0_n = 0.72
yf_n = 0.76

# PLOT NORMAL

plt.plot([x0_n,xf_n], [y0_n,y0_n], 'c--')
plt.plot([x0_n,xf_n], [yf_n,yf_n], 'c--')
plt.plot([x0_n,x0_n], [y0_n,yf_n], 'c--')
plt.plot([xf_n,xf_n], [y0_n,yf_n], 'c--')

# PLOT COM ZOOM

a = plt.axes([0.17, 0.57, 0.15, 0.10]) # x-position, y-position, width, height
plt.xlim(x0_n,xf_n)
plt.ylim(y0_n,yf_n)
plt.grid()
plt.plot(t1_,ref_4,'--k',label='ref')
plt.plot(t1_,y1_[0,:],'mediumvioletred',label='c1')
plt.plot(t2_,y2_[0,:],'goldenrod',label='c2')
plt.plot(t3_,y3_[0,:],"-.b",label='c3')
plt.plot(t4_,y4_[0,:],"-.g",label='c4')


plt.subplot(2,1,2)
plt.plot(t1_,y1_[1,:],'mediumvioletred',label='c1')
plt.plot(t2_,y2_[1,:],'goldenrod',label='c2')
plt.plot(t3_,y3_[1,:],"-.b",label='c3')
plt.plot(t4_,y4_[1,:],"-.g",label='c4')
plt.ylabel('u(t) [m²/s²]')
plt.xlabel('Tempo [s]')
plt.ylim(0.5, 6.2)
plt.xlim(-0.5, tmp_final_2)
plt.legend()
plt.legend(prop={'size':10})
plt.grid()
plt.show()


"""Avaliação da resposta da malha fechada - item 11"""

#Normalização do sistema   

# Normalizando a curva de resposta com o controlador 1 (P_zn)
yn_P_ZN = (y_P_ZN[0] - np.radians(43))/np.radians(3.45) # Normalização

# Normalizando a curva de resposta com o controlador 2 (PI_zn)
yn_PI_ZN = (y_PI_ZN[0] - np.radians(43))/np.radians(3.45) # Normalização

# Normalizando a curva de resposta com o controlador 3 (P_chr)
yn_P_CHR = (y_P_CHR[0] - np.radians(43))/np.radians(3.45) # Normalização

# Normalizando a curva de resposta com o controlador 4 (PI_chr)
yn_PI_CHR = (y_PI_CHR[0] - np.radians(43))/np.radians(3.45) # Normalização

# Vetores de referência para os gráficos normalizados
refn = np.ones(time.shape)              # Referência do ponto de equilíbrio
reftp = (1+ 0.02)*np.ones(time.shape)   # Referência da margem de +2%
reftn = (1- 0.02)*np.ones(time.shape)   # Referência da margem de -2%

plt.figure(8)
# Controlador 1 (P_zn)
plt.subplot(2,2,1)
plt.plot(time,yn_P_ZN,label='c1')
plt.plot(time, refn*yn_P_ZN[24999], "--r", label = 'Ref')
plt.plot(time, reftp*yn_P_ZN[24999], "--k", label = '$\pm 2\%$ ref')
plt.plot(time, reftn*yn_P_ZN[24999], "--k")
plt.ylabel('y(t)[rad]')
plt.xlim(15.0,19)
plt.ylim(0.0,1.40)
plt.legend()
plt.grid()

# Controlador 2 (PI_zn)
plt.subplot(2,2,2)
plt.plot(time,yn_PI_ZN,label='c2')
plt.plot(time, refn*yn_PI_ZN[24999], "--r", label = 'Ref')
plt.plot(time, reftp*yn_PI_ZN[24999], "--k", label = '$\pm 2\%$ ref')
plt.plot(time, reftn*yn_PI_ZN[24999], "--k")
plt.xlim(14.7,22)
plt.ylim(0.0,1.75)
plt.legend()
plt.grid()

# Controlador 3 (P_chr)
plt.subplot(2,2,3)
plt.plot(time,yn_P_CHR,label='c3')
plt.plot(time, refn*yn_P_CHR[24999], "--r", label = 'Ref')
plt.plot(time, reftp*yn_P_CHR[24999], "--k", label = '$\pm 2\%$ ref')
plt.plot(time, reftn*yn_P_CHR[24999], "--k")
plt.ylabel('y(t)[rad]')
plt.xlabel('Tempo [s]')
plt.xlim(15,17)
plt.ylim(0.0,0.9)
plt.legend()
plt.grid()

# Controlador 4 (PI_chr)
plt.subplot(2,2,4)
plt.plot(time,yn_PI_CHR,label='c4')
plt.plot(time, refn*yn_PI_CHR[24999], "--r", label = 'Ref')
plt.plot(time, reftp*yn_PI_CHR[24999], "--k", label = '$\pm 2\%$ ref')
plt.plot(time, reftn*yn_PI_CHR[24999], "--k")
plt.xlabel('Tempo [s]')
plt.xlim(15,18.5)
plt.ylim(0.0,1.35)
plt.legend()
plt.grid()
plt.show()


"""Performance da malha fechada - IAE, ITAE e RMSE - Item 12"""

# Para todos os métdos abaixo, os 4 primeiros valores dos vetores estão relacionados ao seguimento de referência e os 4 últimos estão relacionados à rejeição de perturbação.

# Calculando o erro em relação à referência 
erro_0 = np.abs(y_P_ZN[0] - ref_3)      # Calculo do erro do Controlador Proporcional - Zigler-Nichols - seguimento de referência
erro_1 = np.abs(y_PI_ZN[0] - ref_3)     # Calculo do erro do Controlador Proporcional Integral - Zigler-Nichols - seguimento de referência
erro_2 = np.abs(y_P_CHR[0] - ref_3)     # Calculo do erro do Controlador Proporcional - CHR - seguimento de referência
erro_3 = np.abs(y_PI_CHR[0] - ref_3)    # Calculo do erro do Controlador Proporcional Integral - seguimento de referências

erro_4 = np.abs(y1_[0] - ref_4)        # Calculo do erro do Controlador Proporcional - Zigler-Nichols - rejeição à perturbação
erro_5 = np.abs(y2_[0] - ref_4)        # Calculo do erro do Controlador Proporcional Integral - Zigler-Nichols - rejeição à perturbação
erro_6 = np.abs(y3_[0] - ref_4)        # Calculo do erro do Controlador Proporcional - CHR - rejeição à perturbação
erro_7 = np.abs(y4_[0] - ref_4)        # Calculo do erro do Controlador Proporcional Integral - CHR - rejeição à perturbação

# Cálculo pelo método IAE 
IAE = np.zeros(8) # Alocando um vetor de 8 posições para o erro pelo método de IAE

IAE[0] = T*np.sum(erro_0)               # Calculo do IAE do Controlador Proporcional - Zigler-Nichols - seguimento de referência
IAE[1] = T*np.sum(erro_1)               # Calculo do IAE do Controlador Proporcional Integral - Zigler-Nichols - seguimento de referência
IAE[2] = T*np.sum(erro_2)               # Calculo do IAE do Controlador Proporcional - CHR - seguimento de referência
IAE[3] = T*np.sum(erro_3)               # Calculo do IAE do Controlador Proporcional Integral - seguimento de referências

IAE[4] = T*np.sum(erro_4)               # Calculo do IAE do Controlador Proporcional - Zigler-Nichols - rejeição à perturbação
IAE[5] = T*np.sum(erro_5)               # Calculo do IAE do Controlador Proporcional Integral - Zigler-Nichols - rejeição à perturbação
IAE[6] = T*np.sum(erro_6)               # Calculo do IAE do Controlador Proporcional - CHR - rejeição à perturbação
IAE[7] = T*np.sum(erro_7)               # Calculo do IAE do Controlador Proporcional Integral - CHR - rejeição à perturbação


print(f"Seguimento de referência - IAE - ZN  - Proporcional:             {IAE[0]}")
print(f"Seguimento de referência - IAE - ZN  - Proporcional Intregratal: {IAE[1]}")
print(f"Seguimento de referência - IAE - CHR - Proporcional:             {IAE[2]}")
print(f"Seguimento de referência - IAE - CHR - Proporcional Integral:    {IAE[3]}")

print(f"Rejeição a pertubação -    IAE - ZN  - Proporcional:             {IAE[4]}")
print(f"Rejeição a pertubação -    IAE - ZN  - Proporcional Intregratal: {IAE[5]}")
print(f"Rejeição a pertubação -    IAE - CHR - Proporcional:             {IAE[6]}")
print(f"Rejeição a pertubação -    IAE - CHR - Proporcional Integral:    {IAE[7]}")


# Cálculo pelo método ITAE para degrais seccionados 

ITAE = np.zeros(8) # Alocando um vetor de 8 posições para o erro pelo método de ITAE

ITAE[0] = T*np.sum(time[0:20000]*(erro_0[0:20000] + erro_0[20000:40000] + erro_0[40000:60000] + erro_0[60000:80000] + erro_0[80000:100000] + erro_0[100000:120000] + erro_0[120000:140000] + erro_0[140000:160000])/8) # Calculo do ITAE Controlador Proporcional por Zigler-Nichols para a sequencia de degraus
ITAE[1] = T*np.sum(time[0:20000]*(erro_1[0:20000] + erro_1[20000:40000] + erro_1[40000:60000] + erro_1[60000:80000] + erro_1[80000:100000] + erro_1[100000:120000] + erro_1[120000:140000] + erro_1[140000:160000])/8) # Calculo do ITAE Controlador Proporcional Integral por Zigler-Nichols para a sequencia de degraus
ITAE[2] = T*np.sum(time[0:20000]*(erro_2[0:20000] + erro_2[20000:40000] + erro_2[40000:60000] + erro_2[60000:80000] + erro_2[80000:100000] + erro_2[100000:120000] + erro_2[120000:140000] + erro_2[140000:160000])/8) # Calculo do ITAE Controlador Proporcional por CHR para a sequencia de degraus
ITAE[3] = T*np.sum(time[0:20000]*(erro_3[0:20000] + erro_3[20000:40000] + erro_3[40000:60000] + erro_3[60000:80000] + erro_3[80000:100000] + erro_3[100000:120000] + erro_3[120000:140000] + erro_3[140000:160000])/8) # Calculo do ITAE Controlador Proporcional Integral por CHR para a sequencia de degraus

ITAE[4] = T*np.sum(time[0:20000]*(erro_4[0:20000] + erro_4[20000:40000] + erro_4[40000:60000] + erro_4[60000:80000] + erro_4[80000:100000] + erro_4[100000:120000] + erro_4[120000:140000] + erro_4[140000:160000])/8) # Calculo do ITAE Controlador Proporcional por Zigler-Nichols para a rejeição à perturbação;
ITAE[5] = T*np.sum(time[0:20000]*(erro_5[0:20000] + erro_5[20000:40000] + erro_5[40000:60000] + erro_5[60000:80000] + erro_5[80000:100000] + erro_5[100000:120000] + erro_5[120000:140000] + erro_6[140000:160000])/8) # Calculo do ITAE Controlador Proporcional Integral por Zigler-Nichols para a rejeição à perturbação;
ITAE[6] = T*np.sum(time[0:20000]*(erro_6[0:20000] + erro_6[20000:40000] + erro_6[40000:60000] + erro_6[60000:80000] + erro_6[80000:100000] + erro_6[100000:120000] + erro_6[120000:140000] + erro_7[140000:160000])/8) # Calculo do ITAE Controlador Proporcional por CHR para a rejeição à perturbação;
ITAE[7] = T*np.sum(time[0:20000]*(erro_7[0:20000] + erro_7[20000:40000] + erro_7[40000:60000] + erro_7[60000:80000] + erro_7[80000:100000] + erro_7[100000:120000] + erro_7[120000:140000] + erro_7[140000:160000])/8) # Calculo do ITAE Controlador Proporcional Integral por CHR para a rejeição à perturbação;

print(f"Seguimento de referência - ITAE - ZN  - Proporcional:             {ITAE[0]}")
print(f"Seguimento de referência - ITAE - ZN  - Proporcional Intregratal: {ITAE[1]}")
print(f"Seguimento de referência - ITAE - CHR - Proporcional:             {ITAE[2]}")
print(f"Seguimento de referência - ITAE - CHR - Proporcional Integral:    {ITAE[3]}")

print(f"Rejeição a pertubação -    ITAE - ZN  - Proporcional:             {ITAE[4]}")
print(f"Rejeição a pertubação -    ITAE - ZN  - Proporcional Intregratal: {ITAE[5]}")
print(f"Rejeição a pertubação -    ITAE - CHR - Proporcional:             {ITAE[6]}")
print(f"Rejeição a pertubação -    ITAE - CHR - Proporcional Integral:    {ITAE[7]}")


# Cálculo pelo método RMSE ----------------------------------------------------------------

RMSE = np.zeros(8) # Alocando um vetor de 8 posições para o erro pelo método de RMSE

RMSE[0] = np.sqrt((erro_0**2).mean()) # Calculo do RMSE do Controlador Proporcional - Zigler-Nichols - seguimento de referência
RMSE[1] = np.sqrt((erro_1**2).mean()) # Calculo do RMSE do Controlador Proporcional Integral - Zigler-Nichols - seguimento de referência
RMSE[2] = np.sqrt((erro_2**2).mean()) # Calculo do RMSE do Controlador Proporcional - CHR - seguimento de referência
RMSE[3] = np.sqrt((erro_3**2).mean()) # Calculo do RMSE do Controlador Proporcional Integral - CHR - seguimento de referência

RMSE[4] = np.sqrt((erro_4**2).mean()) # Calculo do RMSE do Controlador Proporcional - Zigler-Nichols - rejeição a pertubação
RMSE[5] = np.sqrt((erro_5**2).mean()) # Calculo do RMSE do Controlador Proporcional Integral - Zigler-Nichols - rejeição a pertubação
RMSE[6] = np.sqrt((erro_6**2).mean()) # Calculo do RMSE do Controlador Proporcional - CHR - rejeição a pertubação
RMSE[7] = np.sqrt((erro_7**2).mean()) # Calculo do RMSE do Controlador Proporcional Integral - CHR - rejeição a pertubação


print(f"Seguimento de referência - RMSE - ZN  - Proporcional:             {RMSE[0]}")
print(f"Seguimento de referência - RMSE - ZN  - Proporcional Intregratal: {RMSE[1]}")
print(f"Seguimento de referência - RMSE - CHR - Proporcional:             {RMSE[2]}")
print(f"Seguimento de referência - RMSE - CHR - Proporcional Integral:    {RMSE[3]}")

print(f"Rejeição a pertubação -    RMSE - ZN  - Proporcional:             {RMSE[4]}")
print(f"Rejeição a pertubação -    RMSE - ZN  - Proporcional Intregratal: {RMSE[5]}")
print(f"Rejeição a pertubação -    RMSE - CHR - Proporcional:             {RMSE[6]}")
print(f"Rejeição a pertubação -    RMSE - CHR - Proporcional Integral:    {RMSE[7]}")