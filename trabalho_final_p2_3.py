import numpy as np
from scipy.special import erf, erfinv
import matplotlib.pyplot as plt
import math

SNRdb= np.arange(1,36,1) #vetor com os valores de SNR média para o plot
SNR = np.power(10,SNRdb/10)
sigmaI=0.04 #turbulence strength
sigmaX=0.1 
media=1

#Funções Q e Q^-1
# Q(f) = 0.5 - 0.5 * erf(x/sqrt(2)) -> Função Q(x)
def qfunc(x):
  return 0.5 - 0.5 * erf(x/np.sqrt(2))

# Q^-1(f) = sqrt(2) * erf^-1(1-2*x) -> Função Q(x) inversa
def qfuncinv(x):
  return np.sqrt(2) * erfinv(1 - 2*x)

# Função de cálculo da eficiência espectral
def S(x,N, Po, sigmax):
    resultado = 0
    for j in range(1,N+1):

        if j == 1:
            I = np.sqrt((1/2*x))*qfuncinv(Po)
        
        else:
          M = np.power(2,j)
          I = (1/np.sin(np.pi/M)) * np.sqrt(1/(2*x))*qfuncinv((np.log2(M)/2)*Po)
        
        xj = (np.log(I) + 2*np.power(sigmax,2))/(2*sigmax)
        resultado += qfunc(xj)
        print(xj)
        
    return resultado/2 

# Plotando as curvas de eficiência espectral
# do modelo adaptativo proposto
Ns = np.array([3, 5])
Pos = np.array([1e-2, 1e-3])
y = []
legendas  = []
cont = 0

fig1 = plt.figure(figsize=(9, 6))

for N in Ns:
  for Po in Pos:
    y.append(S(SNR,N,Po,sigmaX))
    plt.plot(SNRdb,y[cont])
    legendas.append(f"N = {N} e Po = {Po}")
    cont += 1
  
# Plotando a curva da capacidade do canal
W = 1 # Considerando a banda unitária
Cup = (W/2)*(np.log2(SNR/np.e) - (4*sigmaX**2)/np.log(2))
plt.plot(SNRdb, Cup, '-')
legendas.append("Cap. do Canal")

# Plotando a eficiência espectral do BPSK
y_bpsk1 = []
y_bpsk2 = []

for i in SNRdb:
  if i <= 5.2:
    y_bpsk1.append(0)
    y_bpsk2.append(0)
  elif i > 5.2 and i <= 8:
     y_bpsk1.append(0.5)
     y_bpsk2.append(0)
  else:
    y_bpsk1.append(0.5)
    y_bpsk2.append(0.5)

plt.plot(SNRdb, y_bpsk1, ':')
legendas.append("BPSK, Po = 0.01")
plt.plot(SNRdb, y_bpsk2, ':')
legendas.append("BPSK, Po = 0.001")


plt.legend(legendas)
plt.xlabel('SNR [dB]', fontsize=16)
plt.ylabel('Eficiência Espectral [bit/s/Hz]', fontsize=16)
plt.xlim(0,35)
plt.ylim(0,3)
plt.grid()
plt.show()
