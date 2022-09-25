import numpy as np #alias the library as np
import matplotlib.pyplot as plt

# SineWave: y(t) = A*sin(2*pi*f*t + phi)

A = 1      # Amplitude
f = 5      # frequency
phi = 0    # starting position
sr = 100   # num of data points(samples) in 1 second

# numpy arange method returns evenly spaced values w/in
# a given interval
# arguments: startTime, stopTime, step
time = np.arange(0, 10, 1/sr)
y = A * np.sin(2*np.pi*f*time + phi)

plt.figure(figsize=(100,4))
plt.plot(time, y)
plt.title("Sine Wave")
plt.xlabel("Time(s)")
plt.ylabel("Amplitude")
plt.show()
