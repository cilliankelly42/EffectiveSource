import numpy as np 
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt

m=0
Trad = 346.751
Omega_r = 0.0181202 
Omega_phi = 0.024575 
N = 41# This needs to change depending on length of orbit data
data= {}

for i in np.arange(N):
    data[i] = np.loadtxt("/home/cillian/.effective_source/data/datam%d_ecc03/loop_%d.txt" %(m,i)) #The columns in data are as follows: r, theta, tp (= t), Re(eff_src), Im(eff_src),

t = np.zeros(len(data))
re_source, im_source = np.zeros(len(data)), np.zeros(len(data))
for i in np.arange(len(data)):
    t[i] = data[i][2]
    re_source[i] = data[i][3]
    im_source[i] = data[i][4]

eff_source = re_source + 1j * im_source

dt = t[1] - t[0]
zero_padding = 1
frequency = np.fft.fftshift(np.fft.fftfreq(n = N * zero_padding, d = dt))
omega = 2 * np.pi * frequency

# Fourier transforms: 
# 1: Real 
fft_re_source = np.fft.fftshift(np.fft.fft(re_source, n = N * zero_padding)) * dt

# 2: Imag
fft_im_source = np.fft.fftshift(np.fft.fft(im_source, n = N * zero_padding)) * dt

#Times by e^i m Omega_phi t and compute n-modes individually

# 3: Complex
fft_eff_source = np.fft.fftshift(np.fft.fft(eff_source, n = N * zero_padding)) * dt

# print(fft_re_source + 1j * fft_im_source - fft_eff_source)
# print(eff_source - re_source - 1j * im_source)

# print(np.argmax(fft_eff_source))
# print(fft_eff_source[264])

# plt.plot(omega, fft_im_source.imag)
# plt.plot(omega, np.abs(fft_eff_source))
number_of_modes = 20
fourier_modes = np.zeros(number_of_modes)
n = np.zeros(number_of_modes)
for i in np.arange(number_of_modes):
    n[i] = i
    fourier_modes[i] = ((1 / Trad) * np.trapezoid(eff_source * np.exp(1j * i * Omega_r * t), t)) 

"""
integral = 0
for i in np.arange(len(eff_source)):
    integral += eff_source[i] * np.exp(1j * t[i]) * dt
print(1/(2 * np.pi) * integral, "\n")
print(1/(2 * np.pi) * np.trapezoid(eff_source * np.exp(1j * t), t))
"""

# plt.scatter(n, np.abs(fourier_modes))
plt.yscale("log")
plt.plot(t, np.abs(eff_source))
plt.xlabel("t")
plt.ylabel("$S_{eff}^{m=0}$")
#plt.plot(omega, fft_re_source**2)
# plt.plot(omega, np.sqrt(fft_im_source**2 + fft_re_source**2))
#plt.plot(omega, np.abs(fft_eff_source))
# plt.xlim(-0.6,0.6)
# plt.vlines(x = 0, ymin = min(fft_im_source), ymax = max(fft_im_source))
plt.savefig("/home/cillian/.effective_source/visual/Fourier/images/eff_source_m%d_e06.png" %m)
plt.show()

"""
# plt.plot(t, np.abs(eff_source))
plt.plot(t, eff_source.imag)
plt.plot(t, eff_source.real)
# plt.vlines(x=[280.476, 2 * 280.476], ymin = min(eff_source), ymax = max(eff_source), color = "red")
plt.xlabel("$\omega$")
plt.ylabel("$|F(S_{eff}^{m=%d})$|" %m)
plt.title("FFT of $S_{eff}^{m=%d}$" %m)
plt.xlim(200,650)
# plt.xlim(-0.05,0.05)
#plt.savefig("/home/cillian/.my_spec_proj/Fourier/images/FFT_m%d.png" %m)
plt.show()
"""