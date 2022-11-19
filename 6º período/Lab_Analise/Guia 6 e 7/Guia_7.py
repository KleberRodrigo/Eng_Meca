#Autores: Kleber Junior E Robson Junior
#Laboratório de Sistemas Lineares

import  numpy as np #importa a biblioteca numpy
from    matplotlib import pyplot as plt #Importa a biblioteca pyplot para os gráficos
import  control as ct #Importa a biblioteca control para fazer os calculos
import  math #Importa a biblioteca math

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
    y1 = x[0]
    y2 = x[1]
    return y1,y2

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

def saturation_output(t,x,u,params): #função para realizar a saturação 
    if u < 0:
        return 0  # Quando u assume valores negativos a função retorna 0
    else:
        return u  # Quando u assume valores positivos, a função retorna o próprio u

Saturacao = ct.NonlinearIOSystem(0, saturation_output, states=0, name='saturacao', inputs = ('u'), outputs = ('y'))

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
#t, y = ct.input_output_response(FanPlate, time, [U0, u_massa], pt_inicial) #Simulando sem os degraus

t, y = ct.input_output_response(FanPlate, time, [U0,u_massa], [20,0])  # Função que retorna a resposta do sistema para o degrau


# Plot do sistema com aplicação de degraus na entrada
plt.figure(1)
plt.subplot(3,1,1)
plt.plot(t,referencia,'--k', label='ref')
plt.plot(t,y[0],'b', label='${\\theta(t)}$')
plt.ylabel('${\\theta(t)}$[$\\degree$C]')
#plt.xlim(0,tmp_final)
#plt.ylim(20,55)
plt.legend()
plt.grid()
""""
plt.subplot(3,1,2)
plt.plot(t,U0,'--k', label='ref')
plt.plot(t,y[1],'r', label='$Q_r(t)$')
plt.ylabel('$Q_r(t)$}[J]')
plt.xlim(0,tmp_final)
plt.ylim(20,100000)
plt.legend()
plt.grid()
"""
plt.subplot(3,1,3)
plt.plot(time,Degraus_U_conc,'c', label='${u(t)}$')
plt.ylabel('${u(t)[A]}$')
plt.plot(t,referencia,'--k', label='ref')
plt.xlabel('Tempo[s]')
plt.legend()
#plt.xlim(0,tmp_final)
#plt.ylim(4,10)
plt.grid()
plt.show()
