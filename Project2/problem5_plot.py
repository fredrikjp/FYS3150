import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

with open("problem5.txt") as f:
    n = []
    t = []
    for line in f.readlines():
        n.append(float(line.split()[0]))
        t.append(float(line.split()[1]))
    f.close()


plt.plot(n, t, "*")
plt.title("Runtime of matrix size N")
plt.xlabel("N")
plt.ylabel("time[s]")
plt.savefig("runtime.pdf")
plt.show()

linreg = stats.linregress(np.log(n), np.log(t))
slope = linreg[0]
intercept = linreg[1]
data_fit = slope*np.log(n) + intercept

plt.plot(np.log(n), data_fit)
plt.plot(np.log(n), np.log(t), "*")
plt.xlabel("ln(N)")
plt.ylabel("ln(time [s])")
plt.legend([f"Fitted line with slope {slope:.2f}", "Data points"])
plt.savefig("ln_runtime.pdf")
plt.show()

