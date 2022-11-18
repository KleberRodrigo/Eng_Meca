#Laboratório de Análise de Sistemas Lineares
#@Autores: Thiago Henrique de Faria Costa e Mateus Alves de Sales
#@Data: 18/01/2022

import  numpy as np                      #\     
from    matplotlib import pyplot as plt   #\    
import  control as ct                      #- Importando as bibliotecas necessárias para desenvolvimento do código  
import  statistics as sta                 #/   
import  math                             #/          

plt.close('all') #Fechando todas as abas de gráficos aberta

def model_update(t,x,u,params):
    x1 = x[0] # Posição angular    \
    x2 = x[1] # Velocidade angular /  # Variáveis de estado
    
    return np.array([x2, K_1*(np.cos(x1)**2)*u - (K_2*np.sin(x1) + K_3*x2)]) #Retorna a equação diferencial na forma de vertor


def model_output(t,x,u,params): # A saí­da do sistema (y) é o estado X[0], logo, y = theta.
       
    return x[0]


def K_ganho(massa):
    K_1 = (d_cm*rho_ar*C_a*L_a*L_1)/(2*massa*(((L_t**2)/12)+d_cm**2))
    K_2 = (g*d_cm)/ (((L_t**2)/12)+d_cm**2)
    K_3 = (atrito*(d_cm**2))/(massa*(((L_t**2)/12)+d_cm**2))

    return K_1, K_2, K_3

""" Definindo as varíaveis do sistema"""

L_a = 0.154                                         # Largura da placa móvel de alumínio
L_1 = 0.155                                         # Comprimento da placa abaixo do eixo de rotação.
L_t = 0.270                                         # Comprimento da placa abaixo do eixo de rotação.
d_cm = 0.020                                        # Distância do centro de massa da placa
rho_ar = 1.23                                       # Densidade do ar
m_t = 0.100                                         # Massa total da placa
C_a = 2.05                                          # Coeficiente de arrasto
atrito = 5                                          # Coeficiente de atrito viscoso
g = 9.81                                            # Gravidade
pt_operacao = 50                                    # Ponto de operação
pt_operacao_rad = math.radians(pt_operacao)         # Ponto de operação em radianos
pt_inicial = 0                                      # Ponto inicial

"""Encontrando o valor das constantes K do sistema"""

K_1, K_2, K_3 = K_ganho(m_t)

tmp_final = 24   # Tempo final da análise (s)
periodo = 0.001  # Período de análise (s)
time = np.arange(0, tmp_final, periodo) # Criando um lista que vai de 0 até 24, com período 1
U_ent = (K_2 * math.sin(pt_operacao_rad))/(K_1 * (math.cos(pt_operacao_rad)**2))
U_0 = U_ent * np.ones(time.shape) # Cria uma lista do tamanho da lista de tempo com todos os valores iguais ao valor de entrada
ref = pt_operacao_rad * np.ones(time.shape) # Cria uma lista do tamanho da lista de tempo com todos os valores iguais ao ponto de operação; é a linha de referência utilizada nos gráficos


""" Criando os degraus no sistema"""

Degraus_U = U_ent * np.ones(len(time))       # Criando um novo vetor para alocar os degraus
Degraus_U = np.array_split(Degraus_U, 3)     # Dividindo o vetor em 3 partes iguais

Degraus_U[0][:]=U_0[0]                       # Definindo a primeira parte do vetor como U_0
Degraus_U[1][:]= 1.20*U_0[0]                 # Definindo a segunda parte do vetor como +20% de U_0
Degraus_U[2][:]= U_0[0]                      # Definindo a terceira parte do vetor como U_0   

Degraus_U_conc = np.concatenate([Degraus_U[0], Degraus_U[1], Degraus_U[2]]) # Unindo os 3 vetores que antes foram divididos

FanPlate = ct.NonlinearIOSystem(model_update, model_output , states=2, name='FanPlate', inputs = ('u'), outputs = ('y'))

""" Simulação do sistema sem aplicação de degraus na entrada:"""

t, y = ct.input_output_response(FanPlate, time, U_0, pt_inicial)

""" Plotando sistema sem aplicação de degraus na entrada:"""

plt.figure(1)
plt.subplot(2,1,1)
plt.plot(t,ref,'--k', label='ref(0,873 rad)')
plt.plot(t,y,'b', label='$\\theta$(t)')
plt.ylabel('$\\theta$(t)[rad]')
plt.xlim(0,10)
plt.ylim(0,1.5)
plt.legend()
# plt.title('Resposta temporal do sistema em malha aberta sem degrau') #título comentado para evitar redundâncias no relatório
plt.grid()
plt.subplot(2,1,2)
plt.plot(time,U_0,'b', label='$u_0$')
plt.ylabel('U(t)[m^2/s^2]')
plt.xlabel('Tempo[s]')
plt.legend()
plt.xlim(0,24)
plt.ylim(55,70)
plt.grid()
plt.show()


""" Plotando sistema com aplicação de degraus na entrada:"""

t_degrau, y_degrau = ct.input_output_response(FanPlate, time, Degraus_U_conc, pt_inicial)
plt.figure(2)
plt.subplot(2,1,1)
plt.plot(t_degrau,ref,'--k', label='ref')
plt.plot(t_degrau,y_degrau,'b', label='u(t)')
plt.ylabel('$\\theta$(t)[rad]')
plt.xlim(0,24)
plt.ylim(0,1.5)
plt.legend()
# plt.title('Resposta temporal do sistema em malha aberta com degrau')
plt.grid()
plt.subplot(2,1,2)
plt.plot(time,Degraus_U_conc,'b', label='u (t)')
plt.ylabel('U(t)[m^2/s^2]')
plt.xlabel('Tempo[s]')
plt.legend()
plt.xlim(0,24)
plt.ylim(60,75)
plt.grid()
plt.show()


""" Fazendo a normalização da curva obtida do degrau positivo analisado:"""

y_n = y_degrau[8000:16000]
t_n = t_degrau[0:len(y_n)] 

theta_n = (y_n - min(y_n))/(y_degrau[11999] - min(y_n))

""" Referências para acomadação: +- 2% e +- 5%:"""

ref_n = np.ones(t_n.shape)            # Referência do ponto de equilíbrio
ref_tsp_2 = 1.02*np.ones(t_n.shape)   # Referência +2%
ref_tsn_2 = 0.98*np.ones(t_n.shape)   # Referência -2%
ref_tsp_5 = 1.05*np.ones(t_n.shape)   # Referência +5%
ref_tsn_5 = 0.95*np.ones(t_n.shape)   # Referência -5%

""" Verificando o último gráfico
Definindo os principais parâmetros que caracterizam a resposta transitória:"""

tp = 0.348          # Instante de pico
b = 0.5877          # y(tp)
Mp = b/1            # Sobressinal máximo 
tr = 0.1909         # Tempo de subida
ts = 2.498          # Tempo de acomodação (2%)
zeta = 0.17465      # Constante de amortecimento
wn = 9.16848        # Frequência natural do sistema

""" Definindo as curvas envoltórias:"""

E_s = 1 + (np.e**((-zeta)*wn*t_n))/(np.sqrt(1-(zeta**2)))  # Curva envoltória superior
E_i = 1 - (np.e**((-zeta)*wn*t_n))/(np.sqrt(1-(zeta**2)))  # Curva envoltória inferior

""" Plotando Gráfico da resposta e das curvas envoltórias:"""
plt.figure(3)
plt.plot(t_n,ref_n,'--k',label ='ref')
plt.plot(t_n,theta_n,'b',label ='$\\theta$(t)')
plt.plot(t_n,ref_tsp_2,'--y',label ='ref $\pm$ 2%')
plt.plot(t_n,ref_tsn_2,'--y')
plt.plot(t_n,ref_tsp_5,'--g',label='ref $\pm$ 5%')
plt.plot(t_n,ref_tsn_5,'--g')
plt.plot(t_n,E_s,'r', label = 'Curvas envoltórias')
plt.plot(t_n,E_i,'r')
plt.ylabel('$\\theta(t)$ [rad]')
plt.xlabel('Tempo [s]')
plt.xlim(0,4)
plt.ylim(0,2)
plt.legend()
plt.grid()
plt.show()

""" Definindo novos degraus para sistema"""

Degraus_U_eq_transf = U_ent * np.ones(len(time))                #Criando um novo vetor para alocar os degraus
Degraus_U_eq_transf = np.array_split(Degraus_U_eq_transf, 8)    # Dividindo o vetor em 8 partes iguais

Degraus_U_eq_transf[0][:]=U_ent                                 # Definindo a primeira parte do vetor como U_ent
Degraus_U_eq_transf[1][:]=1.10*U_ent                            # Definindo a segunda parte do vetor como +10% de U_ent
Degraus_U_eq_transf[2][:]=0.90*U_ent                            # Definindo a terceira parte do vetor como -10% de U_ent
Degraus_U_eq_transf[3][:]= 1.20*U_ent                           # Definindo a quarta parte do vetor como +20% de U_ent
Degraus_U_eq_transf[4][:]= 0.80*U_ent                           # Definindo a quinta parte do vetor como -20% de U_ent
Degraus_U_eq_transf[5][:]= U_ent                                # Definindo a sexta parte do vetor como +30% de U_ent
Degraus_U_eq_transf[6][:]= 1.3*U_ent                            # Definindo a sétima parte do vetor como U_ent
Degraus_U_eq_transf[7][:]= 0.70*U_ent                           # Definindo a oitava parte do vetor como -30% de U_ent

Degraus_U_conc_eq_transf = np.concatenate([Degraus_U_eq_transf[0], Degraus_U_eq_transf[1], Degraus_U_eq_transf[2], Degraus_U_eq_transf[3], Degraus_U_eq_transf[4], Degraus_U_eq_transf[5], Degraus_U_eq_transf[6], Degraus_U_eq_transf[7]]) # Junção dos vetores que foram separados

t_degrau2, y_degrau2 = ct.input_output_response(FanPlate, time, Degraus_U_conc_eq_transf, 0)

""" Equação após transformada de Laplace - Modelo Linear"""

s=ct.tf("s") # Domínio da frequência

""" Determinando o valor do ganho estático - K"""

K = abs((math.radians(50) - y_degrau[9999])/(1.20*U_ent- U_ent))

" Equação geral para sistemas subamortecidos"

Eq_transferencia = K*(wn**2)/((s**2)+(2*s*wn*zeta)+(wn**2)) #Equação de transferência

timeTransfer, tempOut = ct.forced_response(Eq_transferencia,T=time, U=Degraus_U_conc_eq_transf-U_ent) 

timeTransfer, thetanonL = ct.input_output_response(FanPlate, time, Degraus_U_conc_eq_transf, [pt_operacao_rad,0]) # Sistema Real

""" Plotando Gráfico para comparação do sistema real e do modelo obtido:"""

plt.figure(4)
plt.subplot(211) 
#plt.title('Equação de Transfêrencia')
plt.plot(timeTransfer,tempOut+(np.radians(50)),'b',label='Modelo')
plt.plot(t_degrau2,thetanonL,'g', label='Real')
plt.plot(time,ref,'--k',label='ref(0.873rad)')
plt.grid()
plt.ylabel('$\\theta$(t)[rad]')
plt.legend()
plt.xlim(0,24)
plt.ylim(0.65,1.05)

plt.subplot(212) 
plt.plot(timeTransfer,Degraus_U_conc_eq_transf,'b', label = '$Q(t)$')
plt.grid()
plt.ylabel('${u(t)}[m^2/s^2]$')
plt.xlabel('Tempo [s]')
plt.legend()
plt.xlim(0,24)
plt.ylim(40,80)
plt.show()