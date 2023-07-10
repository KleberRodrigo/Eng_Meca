from scipy import signal
import matplotlib.pyplot as plt 
import numpy as np
# import scienceplots

# # Configurações de plot:
# plt.style.use([
#     'grid',
#     'retro'
# ])
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['font.size'] = 10

R = 1000
L = 0.0076
C = 230e-9

Fr = 3807   #Frequência de Ressonância
Fs = 21612  #Frequência de corte superior
Fi = 671    #Frequência de corte inferior

#Colocar  a função de transferência aqui.
Vout = [(R/L), 0]
Vin =[1, (R/L), (1/(L*C))]
#------------------------------------__-#

sys = signal.TransferFunction(Vout,Vin)

w = (np.arange(10,1e6,1)*2*np.pi)

w,mag,phase = signal.bode(sys,w)

plt.rcParams['figure.figsize'] = (12, 8)
plt.subplots(sharex=True)
plt.subplots_adjust(hspace=0.03)

plt.subplot(2,1,1)
plt.semilogx((w/(np.pi*2)),mag)
plt.plot((Fr,Fr),((np.min(mag)+2),(np.max(mag))+2),'green',linestyle='dotted',label='Frequência de ressonância')
plt.plot((Fi,Fi),((np.min(mag)+2),(np.max(mag))+2),'gray',linestyle='dotted',label = 'Limite inferior e superior' )
plt.plot((Fs,Fs),((np.min(mag)+2),(np.max(mag))+2),'gray',linestyle='dotted')
plt.grid(True,which="both")
plt.ylabel('$Magnitude(db)$', fontsize=12)
plt.xticks([])
plt.legend()

plt.subplot(2,1,2)
plt.semilogx((w/(np.pi*2)),phase)
plt.plot((Fr,Fr),((np.min(phase)+5),(np.max(phase))+5),'green',linestyle='dotted')
plt.plot((Fi,Fi),((np.min(phase)+5),(np.max(phase))+5),'gray',linestyle='dotted')
plt.plot((Fs,Fs),((np.min(phase)+5),(np.max(phase))+5),'gray',linestyle='dotted')
plt.grid(True,which="both")
plt.xlabel('$Frequência(Hz)$', fontsize=12)
plt.ylabel('$Fase(º)$', fontsize=12)


plt.savefig('Bode_passa_faixa.eps', dpi=600, transparent=True, bbox_inches='tight')
plt.show()

