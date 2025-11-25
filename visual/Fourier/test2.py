import numpy as np
import matplotlib.pyplot as plt

# Sampled time domain
N = 4000
tmin, tmax = -2, 2
t = np.linspace(tmin, tmax, N, endpoint=False)
dt = t[1] - t[0]

pad_factor = 20
# Rectangular pulse: 1 on [-1,1]
def box(t):
    if np.abs(t) <= 1:
        return 1
    else:
        return 0

x = np.zeros(len(t))
for i in np.arange(len(t)):
    x[i] = box(t[i])

# Continuous FT approximation
X = np.fft.fftshift(np.fft.fft(x, N*5)) * dt

# Frequency axis in Hz (cycles per time unit)
f = np.fft.fftshift(np.fft.fftfreq(N*5, d=dt))
omega = 2 * np.pi * f

plt.plot(omega, np.abs(X))
plt.xlabel("f (Hz) -- exactly Mathematica’s frequency")
plt.ylabel("|X(f)|")
plt.xlim(t[0],t[-1])
plt.title("NumPy FFT matching Mathematica's {0, -2π} FourierTransform")
plt.grid(True)
plt.show()
