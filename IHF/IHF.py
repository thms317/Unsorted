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
alpha = 120  # opening angle (degrees)
i = 1  # number of dimers

force = []
extension = []
extension_kulic = []

for f in range(1, 175):
    f /= 10
    force.append(f)
    extension.append(WLC(f, p, L, S, x0))
    extension_kulic.append(WLC(f, NewP(L, p, alpha, i), L, S, x0))

plt.plot(extension, force, label="DNA only ($p = 50$ nm)")
plt.plot(extension_kulic, force, label="DNA + IHF ($p_{app} = 32.4$ nm)")
plt.legend(frameon=False)
plt.tick_params(direction='in', axis="both", bottom="on", top="on", left="on", right="on")
plt.xlabel('Extension ($\mu$m)')
plt.ylabel('Force (pN)')
plt.xlim(0, 1.5)
plt.ylim(-1, 20)
plt.show()

