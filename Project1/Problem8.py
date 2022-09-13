import pyarma as pa
import matplotlib.pyplot as plt
import numpy as np

def exact(x):
    return 1 - (1 - np.exp(-10))*x - np.exp(-10*x)

def error(filename):
    with open(filename) as f:
        x = []
        u = []
        i = 0
        for line in f.readlines():
            x.append(float(line.split()[0]))
            u.append(float(line.split()[1]))
            i += 1
        f.close()
    x = np.array(x[1:-1])
    v = np.array(u[1:-1])
    u = exact(x)
    abs_err = np.log10(np.abs(u - v))
    rel_err =  np.log10(np.abs((u - v)/u))
    max_rel_err = np.max(np.abs((u - v)/u))
    return x, abs_err, rel_err, max_rel_err

x10, abs_err10, rel_err10 , max_rel_err10= error("Problem7_n10.txt")
x100, abs_err100, rel_err100, max_rel_err100 = error("Problem7_n100.txt")
x1000, abs_err1000, rel_err1000, max_rel_err1000 = error("Problem7_n1000.txt")
x10000, abs_err10000, rel_err10000, max_rel_err10000 = error("Problem7_n10000.txt")

print(f" n=10       max_rel_err={max_rel_err10} \n n=100      max_rel_err={max_rel_err100} \n n=1000     max_rel_err={max_rel_err1000} \n n=10000    max_rel_err={max_rel_err10000} \n")

plt.plot(x10, abs_err10, "--")
plt.plot(x100, abs_err100, "--")
plt.plot(x1000, abs_err1000, "--")
plt.plot(x10000, abs_err10000, "--")
plt.xlabel("x")
plt.ylabel("y")
plt.legend(["n=10", "n=100", "n=1000", "n=10000"])
plt.title("Log10 of absolute error")
plt.savefig("abs_err.pdf")
plt.show()

plt.plot(x10, rel_err10, "--")
plt.plot(x100, rel_err100, "--")
plt.plot(x1000, rel_err1000, "--")
plt.plot(x10000, rel_err10000, "--")
plt.xlabel("x")
plt.ylabel("y")
plt.legend(["n=10", "n=100", "n=1000", "n=10000"])
plt.title("Log10 of relative error")
plt.savefig("rel_err.pdf")
plt.show()
