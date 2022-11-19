# Relatório 03 - Laboratório de Análise de Sistema Lineares - Análise das Características da Malha Fechada
#@Alunos: Thiago Henrique de Faria Costa e Mateus Alves de Sales
#@Data: 29/01/2022

import  numpy as np
from    matplotlib import pyplot as plt
import  control as ct


plt.close('all')

'''Definindo as variáveis do sistema'''

pt_operacao = 40    # Ponto de operação (ºC)
V = 10              # Volume do reator químico (m^3)
c_p = 4500          # Capacidade térmica da solução (4500 J/c.m^3)
F = 0.5             # Vazao volumétrica (m^3/s)
h = 15              # Coeficiente de conveção natural (W/mC)
A_ext = 31.4        # Superfície do tanque para troca de calor por convecção (m^2)
Temp_inicial = 20   # Temperatura ambiente (ºC)
Qq = 7900           # Energia química necessária para catalização da reação (J/m^3) 
R = 10000           # Resistência térmica (ohm)
s=ct.tf("s")        # Frequência 

Qr = -((c_p*F*(Temp_inicial-pt_operacao))-Qq-(h*A_ext*(pt_operacao-Temp_inicial))) #Valor de Qr para o ponto de operação
U = Qr/R  #valor do sinal de controle necessário para levar o sistema ao ponto de operação

def model_update(t,x,u,params):
    Temp = x[0] 
    Qr = x[1] 
    Td = u[1]
    Ft = F*(1+Td)
    u = u[0] 
    dTheta = (1/(c_p*V))*(Qr+(c_p*Ft*(Temp_inicial-Temp))-Qq-(h*A_ext*(Temp-Temp_inicial)))
    dQr = (1/(12.5*h*np.pi))*(-Qr+(R*u)) 
    return dTheta,dQr # Retornando a EDO do sistema 

def model_output(t,x,u,params):   
    # Para o caso em estudo, a saída do sistema (y) é o estado X[0], ou seja, y = theta.
    y1 = x[0]
    y2 = x[1]

    return y1,y2

#Realizando a definição do sistema pela biblioteca control
system = ct.NonlinearIOSystem(model_update, model_output, states = 2, name = 'system', inputs = ('u','Td'), outputs = ('y1','y2') )

def coeficiente_angular(valor_y, valor_x, j, k): #função que retorna o coeficiente angular da reta tangente
    return abs(valor_y[j]-valor_y[k]) / abs(valor_x[j]-valor_x[k])

def equacao(vetor1, vetor2, c_angular, k): #função que retorna a equação da reta tangente
    return c_angular*vetor1 - c_angular*vetor1[k] + vetor2[k]

def valor_K(vetor1, vetor2, ponto_de_operacao, variavel_de_controle, k): #função que retorna o valor do ganho método Ziegler-Nichols
    return vetor1[k]-ponto_de_operacao, vetor2[k]-variavel_de_controle, abs(vetor1[k]-ponto_de_operacao)/abs(vetor2[k]-variavel_de_controle)

def saturation_output(t,x,u,params): #função para realizar a saturação 
    if u < 0:
        return 0  # Quando u assume valores negativos a função retorna 0
    else:
        return u  # Quando u assume valores positivos, a função retorna o próprio u

# Definindo a saturação da entrada
Saturacao = ct.NonlinearIOSystem(0, saturation_output, states=0, name='saturacao', inputs = ('u'), outputs = ('y'))


def criar_Degrau(tempo, intervalo, U = U, isDegr = 0): #função que cria degraus de adição ou de multiplicação no sinal desejado. O que define o tipo de degrau são os argumentos passados na inicialização da mesma.
    array = U * np.ones(len(tempo))
    newArray = []
    array = np.array_split(array, len(intervalo))  #Separação do array em n partes
    # Adicionando os degraus para cada intervalo
    if isDegr != 1:
        for index, value in enumerate(intervalo):
            array[index].fill(U * (1 + (value/10)))
            newArray = np.concatenate([newArray, array[index]]) # Concatenamento dos arrays que haviam sido separados
    else:
        for index, value in enumerate(intervalo):
            array[index].fill(U  + value)
            newArray = np.concatenate([newArray, array[index]]) # Concatenamento dos arrays que haviam sido separados

    return newArray

''' Variáveis de simulação do sistema (perído, tempo de simulação...) '''
 
tempo_final_validacao = 30000 # Tempo de execução para validação so sistema 
periodo = 0.1 
tempo_validacao = np.arange(0,tempo_final_validacao,periodo)  # Array do tempo de validação do sistema

X0 = 0                                      # Condição inicial do sistema
tempo_final = 20000                         # Tempo de execução do sistema
tempo = np.arange(0,tempo_final,periodo)    # Array do tempo de execução do sistema
degrau_u = criar_Degrau(tempo, [0, 2])      # Array preenchido com a entrada para o ponto de operação que não será alterada
Td = np.zeros(len(tempo))                   # Criando vetor Td com zeros, com o tamanho do tempo

ref_T = pt_operacao*np.ones(tempo.shape)    # Definindo o vetor de referência do ponto de operação
ref_Qr = Qr*np.ones(tempo.shape)            # Definindo o vetor de referência do ponto de operação para o tempo de validação
ref_U = U*np.ones(tempo.shape)

t, y = ct.input_output_response(system, tempo, [degrau_u,Td], [20,0])  # Função que retorna a resposta do sistema para o degrau


# Plot do sistema com aplicação de degraus na entrada
plt.figure(1)
plt.subplot(3,1,1)
plt.plot(t,ref_T,'--k', label='ref')
plt.plot(t,y[0],'b', label='${\\theta(t)}$')
plt.ylabel('${\\theta(t)}$[$\\degree$C]')
plt.xlim(0,tempo_final)
plt.ylim(20,55)
plt.legend()
plt.grid()

plt.subplot(3,1,2)
plt.plot(t,ref_Qr,'--k', label='ref')
plt.plot(t,y[1],'r', label='$Q_r(t)$')
plt.ylabel('$Q_r(t)$}[J]')
plt.xlim(0,tempo_final)
plt.ylim(20,100000)
plt.legend()
plt.grid()

plt.subplot(3,1,3)
plt.plot(tempo,degrau_u,'c', label='${u(t)}$')
plt.ylabel('${u(t)[A]}$')
plt.plot(t,ref_U,'--k', label='ref')
plt.xlabel('Tempo[s]')
plt.legend()
plt.xlim(0,tempo_final)
plt.ylim(4,10)
plt.grid()
plt.show()

coeficiente_angular_degrau = coeficiente_angular(y[0], tempo, 100500, 101000) # Calculando o coeficiente angular da reta tangente ao degrau de subida
equacao_degrau = equacao(tempo, y[0], coeficiente_angular_degrau, 100500) # Equação da reta tangente ao degrau de subida 

plt.figure(2)
plt.plot(tempo,y[0],'b',label='$\\theta(t)$')
plt.plot(tempo,equacao_degrau,'r',label='')
plt.plot(tempo,y[0][99999] * np.ones(tempo.shape),'--k',label='ref resposta degrau')
plt.plot(tempo,y[0][195000] * np.ones(tempo.shape),'--g',label='ref resposta')
plt.grid()
plt.ylabel('${\\theta(t)}$[$\\degree$C]')
plt.legend()
plt.ylim(35,55)
plt.xlim(9300,13000)

# Definindo as coordenadas do zoom

x0 = 9950
xf = 10200
y0 = 38
yf = 42

x0_ = 9950
xf_ = 10200
y0_ = 38
yf_ = 42

# PLOT NORMAL

plt.plot([x0_,xf_], [y0_,y0_], 'c--')
plt.plot([x0_,xf_], [yf_,yf_], 'c--')
plt.plot([x0_,x0_], [y0_,yf_], 'c--')
plt.plot([xf_,xf_], [y0_,yf_], 'c--')

# PLOT COM ZOOM

a = plt.axes([0.60, 0.55, 0.15, 0.15]) # x-position, y-position, width, height
plt.xlim(x0,xf)
plt.ylim(y0,yf)
plt.grid()
plt.plot(tempo,y[0],'b',label='$g(s)$')
plt.plot(tempo,equacao_degrau,'r',label='ref degrau')
plt.plot(tempo,y[0][49000] * np.ones(tempo.shape),'--k',label='ref resposta degrau')
plt.show()

degrau_u_validacao = criar_Degrau(tempo_validacao, [0, 1, 0, -1, 0]) # Criando os degraus para realizar a validação do modelo junto ao sistema (10)
Td = np.zeros(len(tempo_validacao)) # Criando um vetor de zeros do tamanho do tempo de valiação


tau = 10635 - 10022                                                              # Calculando o tau para o obter a função de transferência do sistema
theta, u, K1 = valor_K(y[0], degrau_u, pt_operacao, U, 195000)                   # Calculando o ganho estático do sistema 
num_pos = [K1]                                                                   # Numerador da função de transferência 
den_pos = [tau, 1]                                                               # Denominador da função de transferência 
G_pos = ct.tf(num_pos,den_pos)                                                   # Função de transferência do sistema
t, mz = ct.forced_response(G_pos,T= tempo_validacao, U = degrau_u_validacao - U) # Simulando a saída de um sistema linear
t_validacao, ms = ct.input_output_response(system, tempo_validacao, [degrau_u_validacao, Td], [pt_operacao,Qr]) # Resposta do sistema original para os degraus aplicados

print(f'Valor de Qr: ', Qr)
print(f'Valor de U:  ', U)
print(f'delta Y:   ', theta)
print(f'delta u:   ', u)
print(f'K      :   ', K1)
print(f'tau    : ', tau)
print(f'Função transferência modelo de primeira ordem: ', G_pos)

'''Validação do modelo de primeira ordem junto ao sistema (10)'''

plt.figure(3)
plt.subplot(211)
plt.plot(t_validacao, ms[0], 'r', label = 'Sistema não linear')
plt.plot(t_validacao, mz+pt_operacao, 'g', label = 'Modelo')
plt.plot(t_validacao, pt_operacao*np.ones(t_validacao.shape), '--k', label = 'Ref.')
plt.ylabel('$\\theta(t)$ [°C]')
plt.xlim(0,30000)
plt.grid()
plt.legend()


plt.subplot(212)
plt.plot(t_validacao, degrau_u_validacao, 'c',label = 'Sinal de controle')
plt.plot(t_validacao, U*np.ones(t_validacao.shape), '--k', label = 'Ref.')
plt.xlabel('Tempo [s]')
plt.ylabel('u(t) [A]')
plt.xlim(0,30000)
plt.grid()
plt.legend()
plt.show()

'''Controladores P e PI considerando o atraso de 8 segundos'''

atraso = 8

N = ct.pade(atraso,5) # Aproximação de padé de quinta ordem
Gd = ct.tf(N[0], N[1])
Atraso = ct.tf2io(Gd, name = 'atraso', inputs = 'u', outputs = 'y') # Criando bloco para o atraso
GdI = Gd * G_pos

#controlador P
Kp_chr = (0.3*tau)/(K1*atraso)          # Valor do Kp - controlador P
Gp_chr = ct.tf(Kp_chr,1)                # Função de transferência - controlador P - Método chr  

#controlador PI
Kpi_chr = ((0.6*tau)/(K1*atraso))/10    # Valor de Kp - controlador PI
Ti_chr = 4*atraso                       # Calculando o tempo integral - controlador PI
num = [Kpi_chr*Ti_chr,Kpi_chr]          # Numerador da função de transferência
den = [Ti_chr,0]                        # Denominador da função de transferência
Gpi_chr = ct.tf(num,den)                # Função de transferência - controlador PI - Método chr      

print(f'Função de transferência controlador P:  ', Gp_chr)
print(f'Função de transferência controlador PI: ', Gpi_chr)

controlador_P = ct.tf2io(Gp_chr, name='controlador_P', inputs='u', outputs='y')    # Controlador P
controlador_PI = ct.tf2io(Gpi_chr, name='controlador_PI', inputs='u', outputs='y') # Controlador PI

# Controlador P
SysP = ct.InterconnectedSystem(
    (controlador_P,Atraso,Saturacao,system), name='sysdelayed1',
    connections = (('controlador_P.u','-system.y1'),('saturacao.u','controlador_P.y'),
                    ('atraso.u','saturacao.y'),('system.u','atraso.y')),
    inplist = ('controlador_P.u','system.Td','saturacao.u'),
    inputs = ('ref','Td','u0'),
    outlist = ('system.y1','system.y2', 'atraso.u'),
    outputs = ('y1','y2','u'))

# Controlador PI
SysPI = ct.InterconnectedSystem(
    (controlador_PI,Atraso,Saturacao,system), name='sysdelayed2',
    connections = (('controlador_PI.u','-system.y1'),('saturacao.u','controlador_PI.y'),
                    ('atraso.u','saturacao.y'),('system.u','atraso.y')),
    inplist = ('controlador_PI.u','system.Td','saturacao.u'),
    inputs = ('ref','Td','u0'),
    outlist = ('system.y1','system.y2', 'atraso.u'),
    outputs = ('y1','y2','u'))

'''Determinando as funções de ganho de malha, função de sensitividade e sensitividade complementar'''

# Ganho de malha

Lp =  Gp_chr * G_pos     # Ganho de malha controlador P
Lpi = Gpi_chr * G_pos    # Ganho de malha controlador PI
print(f'Função ganho de malha controlador P : ', Lp)
print(f'Função ganho de malha controlador PI: ', Lpi)

# Função sensitividade

Sp = 1/(1+Lp)       # Função sensitividade controlador P
Spi = 1/(1+Lpi)     # Função sensitividade controlador PI
print(f'Função de sensitividade controlador P : ', Sp)
print(f'Função de sensitividade controlador PI : ', Spi)
print(f'S(s)G(s) P: ', Sp*G_pos)
print(f'S(s)G(s) PI: ', Spi*G_pos)

# Função sensitividade complementar

Cp = Lp/(1+Lp)      # Função sensitividade complementar controlador P
Cpi = Lpi/(1+Lpi)   # Função sensitividade complementar controlador PI
print(f'Função de sensitividade complementar controlador P : ', Cp)
print(f'Função de sensitividade complementar controlador PI : ', Cpi)

'''Validação da resposta do sistema em malha fechada para um impulso, um degrau e pertubação terporária na vazão volumétrica'''

tempo_validacao_mf_final = 7000  # Tempo final da simulação de validação da malha fechada
tempo_validacao_mf = np.arange(0, tempo_validacao_mf_final, periodo)   # Vetor tempo da simulação de validação da malha fechada
ref = pt_operacao*np.ones(tempo_validacao_mf.shape)    # Vetor referência com o ponto de operação
Td = np.zeros(len(tempo_validacao_mf))   # Vetor de zeros para o tempo de validação da malha fechada

# Rejeição a pertubação para um impulso

u_impulso = U*np.ones(tempo_validacao_mf.shape)

for k in range(len(tempo_validacao_mf)):
    if tempo_validacao_mf[k] >= 2000 and tempo_validacao_mf[k] <= 2010: # Criando um impulso de 1s
        u_impulso[k] = 4.0*U
    else:
        u_impulso[k] = U

# Função da biblioteca que retorna a resposta da malha fechada para impulso (Controlador P) 
t2, y1 = ct.input_output_response(SysP, tempo_validacao_mf, [ref, Td, u_impulso], [0,0,0,0,0,pt_operacao,Qr]) 
# Função da biblioteca que retorna a resposta da malha fechada para impulso (Controlador PI) 
t2, y2 = ct.input_output_response(SysPI, tempo_validacao_mf, [ref, Td, u_impulso], [0,0,0,0,0,0,pt_operacao,Qr]) 


plt.figure(4)
plt.subplot(311)
plt.plot(t2, y1[0], 'r', label = 'Controlador P')
plt.plot(t2, y2[0], 'b', label = 'Controlador PI')
plt.plot(t2, pt_operacao*np.ones(t2.shape), '--k', label = 'Ref.' )
plt.ylabel('$\\theta (t)$ [°C]')
plt.xlim(0,tempo_validacao_mf_final)
plt.grid()
plt.legend()

plt.subplot(312)
plt.plot(t2, y1[1], 'r', label = 'Controlador P')
plt.plot(t2, y2[1], 'b', label = 'Controlador PI')
plt.plot(t2, Qr*np.ones(t2.shape), '--k', label = 'Ref.' )
plt.ylabel('Qr(t) [J]')
plt.xlim(0,tempo_validacao_mf_final)
plt.grid()
plt.legend()

plt.subplot(313)
plt.plot(t2, u_impulso, 'b', label = 'Impulso')
plt.plot(t2, U*np.ones(t2.shape), '--k', label = 'Ref.')
plt.ylabel('u [A]')
plt.xlabel('Tempo [s]')
plt.xlim(0,tempo_validacao_mf_final)
plt.grid()
plt.legend()
plt.show()

# Rejeição a pertubação para uma entrada degrau com amplitude de 0.25

u_degrau = criar_Degrau(tempo_validacao_mf, [0, 2.5, 0]) # Cria um degrau de multiplicação com amplitude de 0.25 (é passado um valor de 2.5 pois dentro da função esse valor é dividido por 10)

# Função da biblioteca que retorna a resposta da malha fechada para degrau e degrau (Controlador P) 
t4, y3 = ct.input_output_response(SysP, tempo_validacao_mf, [ref, Td, u_degrau], [0,0,0,0,0,pt_operacao,Qr]) 
# Função da biblioteca que retorna a resposta da malha fechada para degrau e degrau (Controlador PI) 
t4, y4 = ct.input_output_response(SysPI, tempo_validacao_mf, [ref, Td, u_degrau], [0,0,0,0,0,0,pt_operacao,Qr]) 

plt.figure(5)
plt.subplot(311)
plt.plot(t4, y3[0], 'r', label = 'Controlador P')
plt.plot(t4, y4[0], 'b', label = 'Controlador PI')
plt.plot(t4, pt_operacao*np.ones(t4.shape), '--k', label = 'Ref.' )
plt.ylabel('$\\theta (t)$ [°C]')
plt.xlim(0,tempo_validacao_mf_final)
plt.grid()
plt.legend()

plt.subplot(312)
plt.plot(t4, y3[1], 'r', label = 'Controlador P')
plt.plot(t4, y4[1], 'b', label = 'Controlador PI')
plt.plot(t4, Qr*np.ones(t4.shape), '--k', label = 'Ref.' )
plt.ylabel('Qr (t) [J]')
plt.xlim(0,tempo_validacao_mf_final)
plt.grid()
plt.legend()

plt.subplot(313)
plt.plot(t4, u_degrau, 'b', label = 'Degrau')
plt.plot(t4, U*np.ones(t4.shape), '--k', label = 'Ref.')
plt.ylabel('u [A]')
plt.xlabel('Tempo [s]')
plt.xlim(0,tempo_validacao_mf_final)
plt.grid()
plt.legend()
plt.show()

# Rejeição a perturbação para uma variação temporária na vazão volumétrica na ordem de 30%

F_volumetrica = criar_Degrau(tempo_validacao_mf, [0,3,0],F) # Cria um degrau de multiplicação da ordem de 30% (é passado um valor de 3 pois dentro da função esse valor é dividido por 10)

# Aplica-se os efeitos do degrau na variação volumétrica para os valores de Qr e U:
Qr_vol = np.ones(len(tempo_validacao_mf))
for j in range(len(tempo_validacao_mf)):
    Qr_vol[j] = -((c_p*F_volumetrica[j]*(Temp_inicial-pt_operacao))-Qq-(h*A_ext*(pt_operacao-Temp_inicial)))

u_vol = np.ones(len(tempo_validacao_mf))
for i in range(len(tempo_validacao_mf)):
    u_vol[i] = Qr_vol[i]/R


# Função da biblioteca que retorna a resposta da malha fechada para variação da vazão volumétrica (Controlador P) 
t6, y5 = ct.input_output_response(SysP, tempo_validacao_mf, [ref, Td, u_vol], [0,0,0,0,0,pt_operacao,Qr])   
# Função da biblioteca que retorna a resposta da malha fechada para variação da vazão volumétrica (Controlador PI) 
t6, y6 = ct.input_output_response(SysPI, tempo_validacao_mf, [ref, Td, u_vol], [0,0,0,0,0,0,pt_operacao,Qr]) 

plt.figure(6)
plt.subplot(311)
plt.plot(t6, y5[0,:], 'r', label = 'Controlador P')
plt.plot(t6, y6[0,:], 'b', label = 'Controlador PI')
plt.plot(t6, pt_operacao*np.ones(t6.shape), '--k', label = 'Ref.' )
plt.ylabel('$\\theta (t)$ [°C]')
plt.xlim(0,tempo_validacao_mf_final)
plt.grid()
plt.legend()

plt.subplot(312)
plt.plot(t6, y5[1,:], 'r', label = 'Controlador P')
plt.plot(t6, y6[1,:], 'b', label = 'Controlador PI')
plt.plot(t6, Qr*np.ones(t6.shape), '--k', label = 'Ref.' )
plt.ylabel('Qr (t)[J]')
plt.xlim(0,tempo_validacao_mf_final)
plt.grid()
plt.legend()

plt.subplot(313)
plt.plot(t6, F_volumetrica, 'b', label = 'Pertubação Volumétrica')
plt.ylabel('F [$m^3/s$]')
plt.xlabel('Tempo [s]')
plt.xlim(0,tempo_validacao_mf_final)
plt.grid()
plt.legend()
plt.show()