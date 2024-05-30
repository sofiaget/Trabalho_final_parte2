import numpy as np
from scipy.special import erf, erfinv
import matplotlib.pyplot as plt
import math

# Define the SNR values for the plot in dB, extending up to 50 dB
SNRdb = np.arange(1, 51, 1)
SNR = np.power(10, SNRdb / 10)

# Turbulence strengths
sigmaI = 0.43
sigmaX = 0.3
media = 1

# Q-function and its inverse
def qfunc(x):
    return 0.5 - 0.5 * erf(x / np.sqrt(2))

def qfuncinv(x):
    return np.sqrt(2) * erfinv(1 - 2 * x)

# Function to calculate spectral efficiency
def S(x, N, Po, sigmaX):
    resultado = 0
    for j in range(1, N + 1):
        if j == 1:
            I = np.sqrt((1 / 2 * x)) * qfuncinv(Po)
        else:
            M = np.power(2, j)
            I = (1 / np.sin(np.pi / M)) * np.sqrt(1 / (2 * x)) * qfuncinv((np.log2(M) / 2) * Po)

        xj = (np.log(I) + 2 * np.power(sigmaX, 2)) / (2 * sigmaX)
        resultado += qfunc(xj)

    return resultado / 2

# Plotting the curves for the adaptive PSK model
Ns = np.array([3, 5])
Pos = np.array([1e-2, 1e-3])
y = []
legendas = []
cont = 0

fig2 = plt.figure(figsize=(9, 6))

for N in Ns:
    for Po in Pos:
        y.append(S(SNR, N, Po, sigmaX))
        plt.plot(SNRdb, y[cont], marker='o' if N == 3 else 's', linestyle='-', label=f"Adaptive PSK, N={N}, $P_o={Po}$")
        cont += 1

# Plotting the channel capacity curve
W = 1
Cup = (W / 2) * (np.log2(SNR / np.e) - (4 * sigmaX**2) / np.log(2))
plt.plot(SNRdb, Cup, '-', label="High SNR Upper Bound, $C_{up}$")

# Plotting the efficiency of non-adaptive BPSK
y_bpsk1 = []
y_bpsk2 = []

for i in SNRdb:
    if i <= 11:
        y_bpsk1.append(0)
        y_bpsk2.append(0)
    elif i > 11 and i <= 16:
        y_bpsk1.append(0.5)
        y_bpsk2.append(0)
    else:
        y_bpsk1.append(0.5)
        y_bpsk2.append(0.5)

plt.plot(SNRdb, y_bpsk1, ':', label="Non-adaptive BPSK, $P_o=10^{-2}$")
plt.plot(SNRdb, y_bpsk2, ':', label="Non-adaptive BPSK, $P_o=10^{-3}$")

plt.legend()
plt.xlabel('Average Electrical SNR, $\overline{\gamma}$, (dB)', fontsize=16)
plt.ylabel('Spectral Efficiency (bit/s/Hz)', fontsize=16)
plt.xlim(0, 50)
plt.ylim(0, 4)
plt.grid()
plt.show()
