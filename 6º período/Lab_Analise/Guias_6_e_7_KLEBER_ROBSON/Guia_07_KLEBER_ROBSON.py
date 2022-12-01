#Autores: Kleber Junior E Robson Junior
#Laboratório de Sistemas Lineares

import numpy as np # importando biblioteca numpy
import matplotlib.pyplot as plt # importando biblioteca para plotar os graficos
import control as ct  #importando biblioteca Control

plt.close('all') #fecha todas janelas de plot

#Variaveis passado do sistema:
La = .154
L1 = .155
Lt = .270
d = .02
mt = .005
rho = 1.23
c = 2.05
mi = 5
g = 9.81

#Constantes
K1 = ((d*rho*c*La*L1)/(2*mt*(((Lt**2)/12)+(d**2))))
K2 = ((g*d)/(((Lt**2)/12)+(d**2)))
K3 = ((mi*d**2)/(mt*(((Lt**2)/12)+(d**2))))
#Dinâmica
X0 = [0, 0] #condições iniciais do sistema
def model_update(t, x, u, params): #Definindo a edo do espaço de estado
    
    x1 = x[0] # variavel que relaciona a posicao
    x2 = x[1] # variavel que relaciona a velocidade
    M = u[1] # variavel que relaciona a massa variável
    u = u[0] # variavel que relaciona o sinal de controle
    #Constantes do modelo:
    Constante_1 = ((d*rho*c*La*L1)/(2*M*(((Lt**2)/12)+(d**2))))
    Constante_2 = ((g*d)/(((Lt**2)/12)+(d**2)))
    Constante_3 = ((mi*d**2)/(M*(((Lt**2)/12)+(d**2))))

    #Retorna as equações diferenciais 
    return [x2, ((Constante_1*(np.cos(x1)**2)*u) - ((Constante_2*np.sin(x1)) + (Constante_3*x2)))]

def model_output(t, x, u, params): #ajustando o parametro de saida 
    return x

#Resolvendo o sistema
fanplate = ct.NonlinearIOSystem(model_update, model_output, states=2, name='fanplate', inputs=('u', 'u_m'), outputs=('x1','x2'))

pop = 43 #Ponto de operação escolhido pela dupla
u_eq = (((K2/K1)*np.sin(np.radians(pop)))/(np.cos(np.radians(pop))**2)) #Calculo do u que leva o sistema para o desejado
Ad = 5 #amplitude do degrau
ud = (((K2/K1)*np.sin(np.radians(pop+Ad)))/(np.cos(np.radians(pop+Ad))**2)) #sintonia do controlador
s = ct.tf('s') #Criando variavel s para função de tranferencia
#Duração do atraso:
delay = 4
#Aproximação de padé:
o_p = 5 #Aproximação de quinta ordem de padé
N = ct.pade(delay, o_p) #Função para resolução pelo método de padé
#Função de transferênciaquie relaciona o atraso
Gd = ct.TransferFunction(np.array(N[0]), np.array(N[1]))#Função de tranferencia
print(f'\n Gd: \t {Gd}')#imprimindo a função de tranferencia do atraso
atraso = ct.tf2io(Gd, name='atraso', inputs='u', outputs='y')#transformando a função de tranferencia em um sistema de entradas e saídas:

#Parâmetros para os controladores
k = (np.radians(Ad)/(ud-u_eq))
theta = delay
tau = 1.01
#Controlador P - Ziegler Nichols:
Kc_p = (tau/(k*theta))
ft_p = ct.tf(Kc_p, 1)
sys_p = ct.tf2io(ft_p, name='cp', inputs='u', outputs='y')
omega_p = 7.3732
#Controlador PI - Ziegler Nichols
Kc_pi = ((0.6*tau)/(k*theta)) # Ganho 
Ti_pi = (4*theta) #Constante de tempo
ft_pi1 = ct.tf(Kc_pi, 1) # PI - Ziegler Nichols da Parte Proporcional
ft_pi2 = (Kc_pi/(Ti_pi*s)) # PI - Ziegler Nichols da Parte Integral
sys_pi1 = ct.tf2io(ft_pi1, name='cpi1', inputs='u', outputs='y')
sys_pi2 = ct.tf2io(ft_pi2, name='cpi2', inputs='u', outputs='y')
omega_pi = 5.9132
#Condições iniciais
X0 = np.hstack((np.zeros(o_p), np.radians(pop), 0))
top = 200#Tempo para alcançar o ponto de operação
t = np.arange(0, top+1600, .01)#vetor de tempo:
u0 = u_eq*np.ones(t.shape)#Sinal de equilíbrio:
mc = mt*np.ones(t.shape)#Massa constante
#Sequencia de degrau pré estabelecido
r0 = [0, 1, 0, -1, 0.5, 1, -0.5, 0]
tdeg = 200 #duração de cada degrau
rs = []
rs.append(np.radians(pop)*np.ones(top*100))
for i in range(len(r0)):
    ri = np.radians(pop+(Ad*r0[i]))
    rs.append(ri*np.ones(tdeg*100))
r = np.hstack(rs)

#Função interconectar sistemas
closed_loop1 = ct.InterconnectedSystem(
    (fanplate, atraso, sys_p), name='closed_loop1', 
    connections=(('cp.u', '-fanplate.x1'), ('atraso.u', 'cp.y'), ('fanplate.u', 'atraso.y')), 
    inplist=('cp.u', 'atraso.u', 'fanplate.u_m'),
    inputs=('xref', 'u0', 'u_m'),
    outlist=('fanplate.x1', 'fanplate.x2', 'fanplate.u'),
    outputs=('x1', 'x2', 'u')
)

#Simulação da malha fechada
t, xout = ct.input_output_response(closed_loop1, t, [r, u0, mc], X0)
cp = xout[0]

#Pot da simulação
plt.figure(1)
plt.rcParams['figure.figsize'] = (12, 5)
plt.subplot(2, 1, 1)
plt.plot(t[(top*100):]-top, np.degrees(r[(top*100):]), 'k', linestyle='dashed', label='Ref')
plt.plot(t[(top*100):]-top, np.degrees(cp[(top*100):]),'r',label='$C_P$')
plt.ylabel('$\\theta$ [°]')
plt.legend(loc='lower right')
# plt.title('Controlador Proporcional')
plt.grid()

plt.subplot(2, 1, 2)
plt.plot(t[(top*100):]-top, xout[2][(top*100):], 'r', label='$u_{(t)}$')
plt.ylabel('$u_{(t)}$')
plt.xlabel('Tempo [s]')
plt.grid()
plt.show()

########################################################################################

#Conexão da malha fechada:
closed_loop2 = ct.InterconnectedSystem(
    (fanplate, atraso, sys_pi1, sys_pi2), name='closed_loop2', 
    connections=(('cpi1.u', '-fanplate.x1'), ('atraso.u', 'cpi1.y'), ('cpi2.u', '-fanplate.x1'), ('atraso.u', 'cpi2.y'), ('fanplate.u', 'atraso.y')), 
    inplist=('cpi1.u', 'cpi2.u', 'atraso.u', 'fanplate.u_m'),
    inputs=('xref1', 'xref2', 'u0', 'u_m'),
    outlist=('fanplate.x1', 'fanplate.x2', 'fanplate.u'),
    outputs=('x1', 'x2', 'u')
)

# Simulando com o proporcional integral
t, xout = ct.input_output_response(closed_loop2, t, [r, r, u0, mc], X0)
cpi = xout[0]

#Plot dos resultados
plt.figure(2)
plt.rcParams['figure.figsize'] = (12, 5)
plt.subplot(2, 1, 1)
plt.plot(t[(top*100):]-top, np.degrees(r[(top*100):]), 'k', linestyle='dashed', label='Ref')
plt.plot(t[(top*100):]-top, np.degrees(cpi[(top*100):]),'r',label='$C_{PI}$')
plt.ylabel('$\\theta$ [°]')
plt.legend(loc='lower right')
#plt.title('Controlador Proporcional Integral')
plt.grid()

plt.subplot(2, 1, 2)
plt.plot(t[(top*100):]-top, xout[2][(top*100):], 'r', label='$u_{(t)}$')
plt.ylabel('$u_{(t)}$')
plt.xlabel('Tempo [s]')
plt.grid()
plt.show()

############################################################################################
#Definindo o Ganho de malha:
Gc = ft_p #ganho do controlador

#Aproximação para modelo de primeira ordem
G = ((k*(np.exp(-theta*1j*omega_p)))/((tau*s)+1))
L = Gc*G
print('\nGanho de malha L(s):')
print(L)

#Sensitividade:
S = (1/(1+L))
print('\nSensitividade S(s):')
print(S)
print('\n')

#Sensitividade complementar:
C = (L/(1+L))
print('\nSensitividade Complementar C(s):')
print(C)
print('\n')

############################################################################################3
#Ganho de malha:
Gc = (ft_pi1+ft_pi2) #ganho do controlador
#Via aproximação para modelo de primeira ordem, vide guia 5:
G = ((k*(np.exp(-theta*1j*omega_p)))/((tau*s)+1))
L = Gc*G
print('\nGanho de malha L(s):')
print(L)

#Sensitividade:
S = (1/(1+L))
print('\nSensitividade S(s):')
print(S)

#Sensitividade complementar:
C = (L/(1+L))
print('\nSensitividade Complementar C(s):')
print(C)

#################################################################
# Parâmetros de tempos:
timp = 0.5 #duração do impulso
tdeg = 200 #duração do degrau

#Array de tempo:
t = np.arange(0, ((6*tdeg)+(timp)), .01)

#Criação dos sinais:
Aimp = 20 #amplitude do impulso
Adeg = 3 #amplitude do degrau
rop = (np.radians(pop)*np.ones(int(tdeg*100))) #sinal que aplica o degrau
rimp = (np.radians(pop+Aimp)*np.ones(int(timp*100))) #imulso
rdeg = (np.radians(pop+Adeg)*np.ones(int(tdeg*100))) #degrau
r = np.concatenate((rop, rimp, rop, rdeg, rop, rop, rop)) #array de referência
#Sinal de equilíbrio:
u0 = u_eq*np.ones(t.shape)

#Vetor de massa:
M = np.concatenate(((mt*np.ones(int(((4*tdeg)+(timp))*100))), (1.3*mt*np.ones(int(tdeg*100))), (mt*np.ones(int(tdeg*100)))))

###################################################################################

# Simulando Malha Fechada:
X0 = np.zeros(7)
t, xout = ct.input_output_response(closed_loop1, t, [r, u0, M], X0)
cp = xout[0]

#Plotando o resultado da simulação-------------------------------------------------------------------
plt.figure(3)
plt.rcParams['figure.figsize'] = (12, 5)
plt.subplot(3, 1, 1)
plt.plot(t, np.degrees(r), 'k', linestyle='dashed', label='$Ref$')
plt.plot(t, np.degrees(cp), 'r', label='$P_{ZN}$')
plt.ylabel('$\\theta$ [°]')
plt.legend(loc='lower right')
#plt.ylim(18, 35)
#plt.title('Rejeição a Perturbação: Controlador Proporcional')
plt.grid()

plt.subplot(3, 1, 2)
plt.plot(t, M, 'b', label='$m_{(t)}$')
plt.ylabel('m [Kg]')
plt.grid()

plt.subplot(3, 1, 3)
plt.plot(t, xout[2], 'g', label='$u_{(t)}$')
plt.ylabel('$u_{(t)}$')
plt.xlabel('Tempo [s]')
plt.grid()
plt.show()

##########################################################################

# Simulando Malha Fechada:
X0 = np.zeros(7)
t, xout = ct.input_output_response(closed_loop2, t, [r, r, u0, M], X0)
cpi = xout[0]

#Plotando o resultado da simulação-------------------------------------------------------------------
plt.figure(4)
plt.rcParams['figure.figsize'] = (12, 5)
plt.subplot(3, 1, 1)
plt.plot(t, np.degrees(r), 'k', linestyle='dashed', label='$Ref$')
plt.plot(t, np.degrees(cpi), 'r', label='$PI_{CHR}$')
plt.ylabel('$\\theta$ [°]')
plt.legend(loc='lower right')
#plt.ylim(18, 35)
#plt.title('Rejeição a Perturbação: Controlador Proporcional Integral')
plt.grid()

plt.subplot(3, 1, 2)
plt.plot(t, M, 'b', label='$m_{(t)}$')
plt.ylabel('m [Kg]')
plt.grid()

plt.subplot(3, 1, 3)
plt.plot(t, xout[2], 'g', label='$u_{(t)}$')
plt.ylabel('$u_{(t)}$')
plt.xlabel('Tempo [s]')
plt.grid()
plt.show()

############################################################################################