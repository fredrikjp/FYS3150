import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import os


work_path = os.path.dirname(__file__) 
fig_path = os.path.join(work_path, "fig/") #path to directory storing figures

plt.rcParams["text.usetex"] =True

def get_data(filename):
    with open(filename) as f:
        e = []
        m = []
        mean_e = []
        mean_m = []
        CV = []
        chi = []
        boolean = len(f.readline().split())==4
        f.seek(0)
        if boolean:
            for line in f.readlines():
                mean_e.append(float(line.split()[0]))
                mean_m.append(float(line.split()[1]))
                CV.append(float(line.split()[2]))
                chi.append(float(line.split()[3]))
            f.close()
            return np.array(mean_e), np.array(mean_m), np.array(CV), np.array(chi)
        else:
            for line in f.readlines():
                e.append(float(line.split()[0]))
                m.append(float(line.split()[1]))
            f.close()
            return np.array(e), np.array(m)

def analytic(T):
    kB = 1 #.380649e-23
    beta = 1/(kB*T)
    N = 4
    mean_e = 4*(np.exp(-4*beta)-np.exp(4*beta))/(N*(np.exp(-4*beta)+np.exp(4*beta)+6))
    mean_m = 4*(np.exp(4*beta)+2)/(4*(np.exp(-4*beta)+np.exp(4*beta)+6))
    CV = 1/(N*kB*T**2)*32*(3*np.exp(-4*beta)+3*np.exp(4*beta)+2)/(np.exp(-4*beta)+np.exp(4*beta)+6)**2
    chi = 1/(N*kB*T)*16*(np.exp(-4*beta)+3*np.exp(4*beta)+3)/(np.exp(-4*beta)+np.exp(4*beta)+6)**2
    return mean_e, mean_m, CV, chi

# 2X2 ising model
""" 
T = 1
ana_mean_e, ana_mean_m, ana_CV, ana_chi = analytic(T)
mean_e, mean_m, CV, chi = get_data("problem4.txt")
n = len(mean_e)
ana_mean_e = np.ones(n)*(ana_mean_e)
ana_mean_m = np.ones(n)*(ana_mean_m)
ana_CV = np.ones(n)*(ana_CV)
ana_chi = np.ones(n)*(ana_chi)


x = np.linspace(0, n, n)

plt.title(f"Temperature T = {T} J/$k_B$")
plt.plot(x, mean_e, label=f"Mean $\\epsilon$")
plt.plot(x, ana_mean_e, "--",label=f"Analytic mean $\\epsilon$")
plt.legend()
plt.xlabel("Monte Carlo cycles")
plt.ylabel("Energy [J]")
plt.savefig("e(cycles).pdf")
plt.show()

plt.title(f"Temperature T = {T} J/k$_B$")
plt.plot(x, mean_m, label="Mean m")
plt.plot(x, ana_mean_m, "--", label=f"Analytic mean m")
plt.legend()
plt.xlabel("Monte Carlo cycles")
plt.ylabel("Magnetization [unitless]")
plt.savefig("m(cycles).pdf")
plt.show()

plt.title(f"Temperature T = {T} J/k$_B$")
plt.plot(x, CV, label=f"$C_V$")
plt.plot(x, ana_CV, "--", label=f"Analytic $C_V$")
plt.legend()
plt.xlabel("Monte Carlo cycles")
plt.ylabel("Heat capacity [1/J]")
plt.savefig("CV(cycles).pdf")
plt.show()

plt.title(f"Temperature T = {T} J/k$_B$")
plt.plot(x, chi, label=f"$\\chi$")
plt.plot(x, ana_chi, "--", label=f"Analytic $\\chi$")
plt.legend()
plt.xlabel("Monte Carlo cycles")
plt.ylabel("Susceptibility [unitless]")
plt.savefig("chi(cycles).pdf")
plt.show()
"""


# 20X20 ising model
"""
mean_e1, mean_m1, CV1, chi1 = get_data("T1_ordered.txt")
mean_e2, mean_m2, CV2, chi2 = get_data("T1_unordered.txt")
mean_e3, mean_m3, CV3, chi3 = get_data("T2_ordered.txt")
mean_e4, mean_m4, CV4, chi4 = get_data("T2_unordered.txt")

n = len(mean_e1)
x = np.linspace(0, n, n)


plt.title("Mean $\\epsilon$")
plt.plot(x, mean_e1, label=f" T=1 J/k$_B$ | ordered")
plt.plot(x, mean_e2, label=f" T=1 J/k$_B$ | unordered")
plt.plot(x, mean_e3, label=f" T=2.4 J/k$_B$ | ordered")
plt.plot(x, mean_e4, label=f" T=2.4 J/k$_B$ | unordered")
plt.legend()
plt.xlabel("Monte Carlo cycles")
plt.ylabel("Energy [J]")
plt.savefig("e(cycles,ordered).pdf")
plt.show()

plt.title("Mean m")
plt.plot(x, mean_m1, label=f" T=1 J/k$_B$ | ordered")
plt.plot(x, mean_m2, label=f" T=1 J/k$_B$ | unordered")
plt.plot(x, mean_m3, label=f" T=2.4 J/k$_B$ | ordered")
plt.plot(x, mean_m4, label=f" T=2.4 J/k$_B$ | unordered")
plt.legend()
plt.xlabel("Monte Carlo cycles")
plt.ylabel("Magnetization [unitless]")
plt.savefig("m(cycles,ordered).pdf")
plt.show()
"""

# 20X20 ising model probability distribution
"""
e1, m1 = get_data("T1_e.txt")
e2, m2 = get_data("T2_e.txt")

bins1 = len(np.unique(e1)) + 1
bins2 = len(np.unique(e2)) + 1

plt.title(f"p($\\epsilon$;T=1[J/k$_B$])")
plt.hist(e1, bins=bins1, density=True)
plt.xlabel(f"$\\epsilon$[J]")
plt.ylabel("Probability density")
plt.savefig("T1_pdf.pdf")
plt.show()

plt.title("p($\\epsilon$;T=2.4[J/k$_B$])")
plt.hist(e2, bins=bins2, density=True)
plt.xlabel("$\\epsilon$[J]")
plt.ylabel("Probability density")
plt.savefig("T2_pdf.pdf")
plt.show()
"""

# 20X20 ising model for differnet temperatures

"""
mean_e1, mean_m1, CV1, chi1 = get_data("L40_model.txt")
mean_e2, mean_m2, CV2, chi2 = get_data("L60_model.txt")
mean_e3, mean_m3, CV3, chi3 = get_data("L80_model.txt")
mean_e4, mean_m4, CV4, chi4 = get_data("L100_model.txt")

n = len(mean_e1)
T = np.linspace(2.1, 2.4, n)

plt.title("Mean $\epsilon$")
plt.plot(T, mean_e1, "--", label=f"L=40")
plt.plot(T, mean_e2, "--", label=f"L=60")
plt.plot(T, mean_e3, "--", label=f"L=80")
plt.plot(T, mean_e4, "--", label=f"L=100")
plt.legend()
plt.xlabel("Temperature[$J/k_B$]")
plt.ylabel("Energy [J]")
plt.savefig("e(T).pdf")
plt.show()

plt.title("Mean $m$")
plt.plot(T, mean_m1, "--", label=f"L=40")
plt.plot(T, mean_m2, "--", label=f"L=60")
plt.plot(T, mean_m3, "--", label=f"L=80")
plt.plot(T, mean_m4, "--", label=f"L=100")
plt.legend()
plt.xlabel("Temperature[$J/k_B$]")
plt.ylabel("Magnetization [unitless]")
plt.savefig("m(T).pdf")
plt.show()

plt.title("Heat capacity $C_V$")
plt.plot(T, CV1, "--", label=f"L=40")
plt.plot(T, CV2, "--", label=f"L=60")
plt.plot(T, CV3, "--", label=f"L=80")
plt.plot(T, CV4, "--", label=f"L=100")
plt.legend()
plt.xlabel("Temperature[$J/k_B$]")
plt.ylabel("$C_V$ [1/J]")
plt.savefig("CV(T).pdf")
plt.show()

plt.title("Susceptibility $\chi$")
plt.plot(T, chi1, "--", label=f"L=40")
plt.plot(T, chi2, "--", label=f"L=60")
plt.plot(T, chi3, "--", label=f"L=80")
plt.plot(T, chi4, "--", label=f"L=100")
plt.legend()
plt.xlabel("Temperature[$J/k_B$]")
plt.ylabel("$\chi$ [unitless]")
plt.savefig("chi(T).pdf")
plt.show()
"""

# Estimating Tc(L=infinity)

"""
L = np.array([40, 60, 80, 100])
Tc = np.array([, , , ])

res = stats.linregress(L, Tc*L)
plt.plot(L, res.intercept + res.slope*L, "--", label = f"Linear least-squares slope {res.slope:.4f}$\pm${res.stderr:.4f}")
plt.plot(L, L*Tc, "*", label = "L*$T_c$")
plt.legend()
plt.xlabel("L")
plt.ylabel("$T_c[J/k_B]$")
plt.show()
"""


