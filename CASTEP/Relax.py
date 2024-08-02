import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as cons
#
eVperA = cons.e/1e-10
forces = np.loadtxt('Forces.dat')
fig, ax = plt.subplots()
ax.plot(forces*eVperA)
ax.set_xlabel("Iteration")
ax.set_ylabel("Force in N")
fig.savefig("Relax.png")
plt.show()