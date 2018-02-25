import matplotlib.pyplot as plt
import numpy as np


# definitions

# worm-like chain
def WLC(f, p, L, S, x0):
    return (L * (1 - 0.5 * (np.sqrt(kBT / (f * p))) + f / S)) / 1000 + x0


# L_app
def NewP(L, p, alpha, i):
    C = 8 * (1 - np.cos(alpha / 4))
    ro = i / L
    return p / ((1 + p * i * C * ro) ** 2)


# constants

kBT = 4.114  # (pN nm) - Boltzmann factor
Lc = 4092  # contour length (bp)
L = Lc * 0.34  # contour length (nm)
p = 50  # persistence length (nm)
S = 1000  # stretch modulus (pN)
x0 = 0  # offset (nm)
alpha = 110  # opening angle (degrees)
i = 1  # number of dimers

force = []
extension = []

for f in range(1, 175):
    f /= 10
    force.append(f)
    extension.append(WLC(f, p, L, S, x0))

files = []
files.append('180109_data_031_72.fit')
files.append('180109_data_046_72.fit')
files.append('180109_data_068_72.fit')


for file_all in files:

    # read file
    f = open(file_all, 'r')
    data_lines = f.readlines()[1:]
    f.close()

    data_force = []
    data_extension = []

    for x in data_lines:
        data_force.append(float(x.split()[0]))
        data_extension.append(float(x.split()[2]))

    plt.scatter(data_extension, data_force, label=str(file_all), s=10)

plt.plot(extension, force, label="WLC", zorder=100, color="black")
plt.xlim(0, 1.5)
plt.ylim(-1, 20)
plt.xlabel('Extension ($\mu$m)')
plt.ylabel('Force (pN)')
plt.legend(frameon=False)
plt.tick_params(direction='in', axis="both", bottom="on", top="on", left="on", right="on")
plt.show()