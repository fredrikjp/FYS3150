import pyarma as pa
import matplotlib.pyplot as plt

u = pa.mat()
x = pa.mat()
u.load("u.bin")
x.load("x.bin")

plt.plot(x,u)
plt.title("u(x)")
plt.xlabel("x")
plt.ylabel("y")
plt.savefig("u(x)_plot.pdf")
plt.show()

