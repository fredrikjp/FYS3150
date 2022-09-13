import pyarma as pa
import matplotlib.pyplot as plt
import numpy as np

def plot_from_file(filename):    
    with open(filename) as f:
        x = []
        u = []
        for line in f.readlines():
            x.append(float(line.split()[0]))
            u.append(float(line.split()[1]))
        f.close()
    plt.plot(np.log10(x), u, "o")

plot_from_file("Thomas_duration.txt")
plot_from_file("Special_duration.txt")
plt.title("Algorithm duration")
plt.legend(["General_alg", "Special_alg"])
plt.xlabel("log10(n)")
plt.ylabel("time[s]")
plt.savefig("Alg_duration.pdf")
plt.show()

 


