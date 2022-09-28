import numpy as np               # importando biblioteca numpy
import matplotlib.pyplot as plt  # importando biblioteca para plot
import control as ct             # importanto biblioteca control

plt.close('all')                 # Fecha todas as janelas de plot

# Equação diferencial que modela o sistema=====================================
def model_update(t, x, u, params):
    
    if(x < 0): x=0 
    return ((-(x**(0.75)) + 0.53025*(x**(0.0037))*R*u)/13.76023)


def model_output(t, x, u, params):
    return x


# Função para uma equação não linear, a entrada será o sinal de controle "u"
# E a saída a lista "y"
SYSTEM = ct.NonlinearIOSystem(
    model_update, model_output, states=1, name='SYSTEM', inputs=('u'), outputs=('y'))

#===========Definindo os parâmetros de simulação===============================
# Variaveis relacionadas ao tempo
t0 = 0                           # Tempo inicial
tf = 2299                        # Tempo final
t = np.linspace(t0, tf, 2300)    # Instantes que desejo ter a solucao

# Tempo de amostragem
T = t[1]-t[0]

# Valor Banco de resistencias ohms
R = 100

# Condição inicial do sistema.
X0 = 30 

# Pontos de Equilíbrio para temperaturas de 115, 120, 121, 124, e 125°C
u_68 = (68**(3/4))/(0.53025*68**(0.0315/8.5)*100)  
u_69 = (69**(3/4))/(0.53025*69**(0.0315/8.5)*100)  
u_70 = (70**(3/4))/(0.53025*70**(0.0315/8.5)*100)
u_71 = (71**(3/4))/(0.53025*71**(0.0315/8.5)*100)
u_72 = (72**(3/4))/(0.53025*72**(0.0315/8.5)*100)
u_73 = (73**(3/4))/(0.53025*73**(0.0315/8.5)*100)


# Sinal de Controle
u = u_70 * np.ones(t.shape)

# Laço que aplica os degraus no sinal de controle
for i in range(len(u)):
    if (i > 400 and i <= 700):
        u[i] = u_68

    elif (i > 700 and i <= 1000):
        u[i] = u_69

    elif (i > 1000 and i <= 1300):
        u[i] = u_70

    elif (i > 1300 and i <= 1600):
        u[i] = u_71
    
    elif (i > 1600 and i <= 1900):
        u[i] = u_72

    elif (i > 1900 and i <= 2200):
        u[i] = u_73

# Executando a simulação do sistema não linear em malha aberta=================
t, y = ct.input_output_response(SYSTEM, t, u, X0)

# Inclinações máximas em cada degrau
inc_max1 = 0
inc_max2 = 0

# Instantes onde cada inclinação máxima ocorre
inc_t1 = 0
inc_t2 = 0

# Laço que calcula pontos de maior inclinação para cada degrau
for i in range(len(y)):
    if(i > 400 and i <= 700 and (y[i]-y[i-1]) < inc_max1):
        inc_t1 = i
        inc_max1 = y[i]-y[i-1]

    elif(i > 700 and i <= 1000 and (y[i]-y[i-1]) > inc_max2):
        inc_t2 = i
        inc_max2 = y[i]-y[i-1]

# Retas tangentes aos pontos de maior inclinação
K1 = inc_max1 * t
K2 = inc_max2 * t

#Variáveis de controle
j = 0
k = 0
l = 0
m = 0
n = 0
o = 0

e_sz = np.nan * np.ones(t.shape)     #Curva de Subida por Ziegler theta 6.2, AC 72.8
e_sm = np.nan * np.ones(t.shape)     #Curva de Subida por Miller theta 6.2, AB 46
e_dz = np.nan * np.ones(t.shape)     #Curva de Descida por Ziegler theta 15.7, AC 85.6
e_dm = np.nan * np.ones(t.shape)     #Curva de Descida por Miller theta 15.7, AB 54.1


for i in range(len(t)):
    if(i > 400 and i <= 700):
        e_sz[i] = 70 + 5*((u[i]-u[400])/(u_73-u_70))*(1 - np.exp((-j + 6.2)/72.8))        #Descida 
        e_sm[i] = 70 + 5*((u[i]-u[400])/(u_73-u_70))*(1 - np.exp((-j + 6.2)/46))
        e_dz[i] = 70 + 5*((u[i]-u[400])/(u_73-u_70))*(1 - np.exp((-j + 15.7)/85.6))
        e_dm[i] = 70 + 5*((u[i]-u[400])/(u_73-u_70))*(1 - np.exp((-j + 15.7)/54.1))
        j = j+1
    
    if(i > 700 and i <= 1000):
        e_sz[i] = e_sz[700] + 5*((u[i]-u[700])/(u_73-u_70))*(1 - np.exp((-k + 6.2)/72.8))        #Descida 
        e_sm[i] = e_sm[700] + 5*((u[i]-u[700])/(u_73-u_70))*(1 - np.exp((-k + 6.2)/46))
        e_dz[i] = e_dz[700] + 5*((u[i]-u[700])/(u_73-u_70))*(1 - np.exp((-k + 15.7)/85.6))
        e_dm[i] = e_dm[700] + 5*((u[i]-u[700])/(u_73-u_70))*(1 - np.exp((-k + 15.7)/54.1))
        k = k+1
    
    if(i > 1000 and i <= 1300):
        e_sz[i] = e_sz[1000] + 5*((u[i]-u[1000])/(u_73-u_70))*(1 - np.exp((-l + 6.2)/72.8))        #Descida 
        e_sm[i] = e_sm[1000] + 5*((u[i]-u[1000])/(u_73-u_70))*(1 - np.exp((-l + 6.2)/46))
        e_dz[i] = e_dz[1000] + 5*((u[i]-u[1000])/(u_73-u_70))*(1 - np.exp((-l + 15.7)/85.6))
        e_dm[i] = e_dm[1000] + 5*((u[i]-u[1000])/(u_73-u_70))*(1 - np.exp((-l + 15.7)/54.1))
        l = l+1
        
    if(i > 1300 and i <= 1600):
        e_sz[i] = e_sz[1300] + 5*((u[i]-u[1300])/(u_73-u_70))*(1 - np.exp((-m + 6.2)/72.8))        #Descida 
        e_sm[i] = e_sm[1300] + 5*((u[i]-u[1300])/(u_73-u_70))*(1 - np.exp((-m + 6.2)/46))
        e_dz[i] = e_dz[1300] + 5*((u[i]-u[1300])/(u_73-u_70))*(1 - np.exp((-m + 15.7)/85.6))
        e_dm[i] = e_dm[1300] + 5*((u[i]-u[1300])/(u_73-u_70))*(1 - np.exp((-m + 15.7)/54.1))
        m = m+1    
    
    if(i > 1600 and i <= 1900):
        e_sz[i] = e_sz[1600] + 5*((u[i]-u[1600])/(u_73-u_70))*(1 - np.exp((-n + 6.2)/72.8))        #Descida 
        e_sm[i] = e_sm[1600] + 5*((u[i]-u[1600])/(u_73-u_70))*(1 - np.exp((-n + 6.2)/46))
        e_dz[i] = e_dz[1600] + 5*((u[i]-u[1600])/(u_73-u_70))*(1 - np.exp((-n + 15.7)/85.6))
        e_dm[i] = e_dm[1600] + 5*((u[i]-u[1600])/(u_73-u_70))*(1 - np.exp((-n + 15.7)/54.1))
        n = n+1  
    
    if(i > 1900 and i <= 2200):
        e_sz[i] = e_sz[1900] + 5*((u[i]-u[1900])/(u_73-u_70))*(1 - np.exp((-o + 6.2)/72.8))        #Descida 
        e_sm[i] = e_sm[1900] + 5*((u[i]-u[1900])/(u_73-u_70))*(1 - np.exp((-o + 6.2)/46))
        e_dz[i] = e_dz[1900] + 5*((u[i]-u[1900])/(u_73-u_70))*(1 - np.exp((-o + 15.7)/85.6))
        e_dm[i] = e_dm[1900] + 5*((u[i]-u[1900])/(u_73-u_70))*(1 - np.exp((-o + 15.7)/54.1))
        o = o+1  

teste = e_sz[401:2201]

#       Raiz do erro quadrático médio para cada aproximacao
RMSE_e_sz = np.sqrt(np.square(np.subtract(y[401:2201], e_sz[401:2201])).mean())   #  Ziegler curva de subida
RMSE_e_dz = np.sqrt(np.square(np.subtract(y[401:2201], e_dz[401:2201])).mean())   #  Ziegler curva de descida
RMSE_e_sm = np.sqrt(np.square(np.subtract(y[401:2201], e_sm[401:2201])).mean())   #  Miller curva de subida
RMSE_e_dm = np.sqrt(np.square(np.subtract(y[401:2201], e_dm[401:2201])).mean())   #  Miller curva de descida

print(inc_t1)
print(inc_t2)
print(inc_max1)
print(inc_max2)

print("Ziegler curva de subida: ", RMSE_e_sz)
print("Ziegler curva de descida: ", RMSE_e_dz)
print("Miller curva de subida: ", RMSE_e_sm)
print("Miller curva de descida: ", RMSE_e_dm)

#==============Temperatura do Ar e Sinal de Controle===========================
plt.figure(1)
#plt.subplot(2,1,1)

plt.plot(t,e_sz,'r', alpha=1, label='Ziegler degrau positivo')
plt.plot(t,e_dz,'g',label='Ziegler degrau negativo')
plt.plot(t,e_sm,'b',label='Miller degrau positivo')
plt.plot(t,e_dm,'m',label='Miller degrau negativo')
plt.plot(t,y,'k',label='T_aq(t)')

plt.title('Comparação dos modelos obtidos', fontsize=15)
plt.ylabel('Temperatura(°C)', fontsize=15)
plt.xlabel('Tempo(s)', fontsize=15)
plt.legend(loc='lower center' ,fontsize=15)
plt.tick_params(labelsize=13)
plt.grid()

#plt.subplot(2,1,2)
#plt.plot(t,u,'b',label='u(t)')

#plt.ylabel('Sinal de controle', fontsize=15)
#plt.legend(fontsize=15)
#plt.xlabel('Tempo (s)', fontsize=15)
#plt.xlim(350,tf)
#plt.grid()
#plt.tick_params(labelsize=13)
plt.show()

'''
plt.figure(1)
#==================Primeiro Degrau=============================================
plt.subplot(1,2,1)

#plt.vlines(x=400, ymin=114, ymax=121, colors = 'b', linestyles='dashed', label='Degrau no Sinal')
#plt.vlines(x=700, ymin=114, ymax=121, colors = 'b', linestyles='dashed')

plt.hlines(y=120, xmin=400, xmax=700, colors = 'c', linestyles='solid', label='Referência')
plt.hlines(y=115, xmin=350, xmax=750, colors = 'r', linestyles='dashed', label='K')
plt.hlines(y=(-5*0.632+120), xmin=350, xmax=750, colors = 'b', linestyles='dashed', label='0,63K')

plt.plot(t,K1+144.28,'g-.', label='Tangente')
plt.plot(t,y,'k',label='Temperatura °C')

#plt.plot(t,e_dm, 'r--' ,label='Miller')
#plt.plot(t,e_dz, 'b--' ,label='Ziegler-Nichols')
plt.xlim(400,600)
plt.ylim(114,121)

plt.scatter([416], [120],label='A')
plt.scatter([469.8], [116.83],label='B')
plt.scatter([501], [115],label='C')

plt.title('Resposta ao primeiro degrau', fontsize=15)
plt.ylabel('Temperatura do Grão (°C)', fontsize=15)
plt.xlabel('Tempo (s)', fontsize=15)
plt.legend(loc='upper right', fontsize=15)
plt.tick_params(labelsize=13)
plt.grid()

#==================Segundo Degrau==============================================
plt.subplot(1,2,2)

#plt.vlines(x=700, ymin=114, ymax=121, colors = 'b', linestyles='dashed', label='Degrau no Sinal')
#plt.vlines(x=1000, ymin=114, ymax=121, colors = 'b', linestyles='dashed')

plt.hlines(y=115.02, xmin=650, xmax=1050, colors = 'c', linestyles='solid', label='Referência')
plt.hlines(y=120, xmin=650, xmax=1050, colors = 'r', linestyles='dashed', label='K')
plt.hlines(y=(+5*0.632+115), xmin=650, xmax=1050, colors = 'b', linestyles='dashed', label='0,63K')

plt.plot(t,K2+66.5,'g-.', label='Tangente')
plt.plot(t,y,'k',label='Temperatura °C')
#plt.plot(t,e_sm, 'r-.' ,label='Miller')
#plt.plot(t,e_sm, 'y-.' ,label='Ziegler-Nichols')

plt.xlim(700,900)
plt.ylim(114,121)

plt.scatter([706.2], [115],label='A')
plt.scatter([752], [118.16],label='B')
plt.scatter([779], [120],label='C')

plt.title('Resposta ao segundo degrau', fontsize=15)
#plt.ylabel('Temperatura do Grão (°C)', fontsize=15)
plt.xlabel('Tempo (s)', fontsize=15)
plt.legend(loc='lower right', fontsize=15)
plt.tick_params(labelsize=13)
plt.grid()
'''
#==================Terceiro Degrau=============================================
'''plt.subplot(2,2,3)
plt.vlines(x=1000, ymin=119, ymax=126, colors = 'b', linestyles='dashed', label='Degrau no Sinal')
plt.vlines(x=1300, ymin=119, ymax=126, colors = 'b', linestyles='dashed')
plt.plot(t,y,'k',label='Temperatura °C')
plt.plot(t,e_dz, 'r-.' ,label='Miller')
plt.plot(t,e_dz, 'y-.' ,label='Ziegler-Nichols')
plt.xlim(950,1350)
plt.ylim(119,126)
plt.legend(loc='lower right')
plt.grid()

#==================Quarto Degrau=============================================
plt.subplot(2,2,4)
plt.vlines(x=1300, ymin=119, ymax=126, colors = 'b', linestyles='dashed', label='Degrau no Sinal')
plt.vlines(x=1600, ymin=119, ymax=126, colors = 'b', linestyles='dashed')
plt.plot(t,y,'k',label='Temperatura °C')
plt.plot(t,e_dmm, 'r-.' ,label='Miller')
plt.plot(t,e_dmz, 'y-.' ,label='Ziegler-Nichols')
plt.xlim(1250,1650)
plt.ylim(119,126)
plt.legend(loc='upper right')
plt.grid()'''

plt.show()
    
