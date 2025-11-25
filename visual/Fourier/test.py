import numpy as np 
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt

N = 5000

def box(x):
    if -1 <= x <= 1:
        return 1
    else:
        return 0

t = np.linspace(-2, 2, N)
dt = t[1] - t[0]

f = np.zeros(len(t))
 
for i in np.arange(len(t)):
    f[i] = box(t[i])

fig = plt.figure()
ax = fig.add_subplot(111)
fft = np.fft.fftshift(np.fft.fft(f, n = N*10)) * dt
freq = np.fft.fftshift(np.fft.fftfreq(N*10, dt))

ax.plot(freq, np.abs(fft))
ax.set_xlim(-5,5)
plt.show()