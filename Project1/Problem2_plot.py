import pyarma as pa
import matplotlib.pyplot as plt

with open("x_u(x).txt") as f:
    x = []
    u = []
    for line in f.readlines():
        x.append(float(line.split()[0]))
        u.append(float(line.split()[1]))
    f.close()


plt.plot(x,u)
plt.title("u(x)")
plt.xlabel("x")
plt.ylabel("y")
plt.savefig("u(x)_plot.pdf")
plt.show()

