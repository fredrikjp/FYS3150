import pyarma as pa
import matplotlib.pyplot as plt

def plot_from_file(filename, line_style = "k-"):    
    with open(filename) as f:
        x = []
        u = []
        for line in f.readlines():
            x.append(float(line.split()[0]))
            u.append(float(line.split()[1]))
        f.close()

    plt.plot(x, u, line_style)

for i in range(1, 5):
    n = 10**i 
    plot_from_file("x_u(x).txt", )
    plot_from_file(f"Problem7_n{n}.txt", "y--")
    plt.title(f"n={n}")
    plt.legend(["u", "v"])
    plt.xlabel("x")
    plt.ylabel("y")
    plt.savefig(f"Problem7_n{n}.pdf")
    plt.show()

