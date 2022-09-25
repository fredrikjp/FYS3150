import matplotlib.pyplot as plt
import numpy as np

for a in range(1,3):
    n = int(10**a)
    with open(f"problem6_n{n}.txt") as f:
        x = []
        jacobi = []
        anal = []
        for line in f.readlines():
            x.append(float(line.split()[0]))
            jacobi.append(float(line.split()[1]))
            anal.append(float(line.split()[2]))
        f.close()

    for i in range(3):
        #Eigenvector 3 for n=100 has opposite signs between Jacobi and analytic solution
        # such that swiching signs of one of them, makes it easier to compare:
        if i == 2 and a == 2: 
            b=np.array(jacobi[i*(n + 2):(i + 1)*(n + 2)])
            plt.plot(x[:(n + 2)], -b)
        else:
            plt.plot(x[:(n + 2)], jacobi[i*(n + 2):(i + 1)*(n + 2)])
        plt.plot(x[:(n + 2)], anal[i*(n + 2):(i + 1)*(n + 2)], "--")
    plt.title(f"Eigen vectors with n = {10**a} steps and L = length of beam")
    plt.xlabel("x[1/L]")
    plt.ylabel("u[1/L]")
    plt.legend(["E1 Jacobi", "E1 Analytic", "E2 Jacobi", "E2 Analytic", "E3 Jacobi", "E3 Analytic"])
    plt.savefig(f"eigenvec_n{n}.pdf")
    plt.show()

