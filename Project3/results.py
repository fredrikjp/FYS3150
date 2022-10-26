import matplotlib.pyplot as plt
import numpy as np
import pyarma as pa
import os

work_path = os.path.dirname(__file__) 
fig_path = os.path.join(work_path, "fig/") #path to directory storing figures


def get_data(filename):
    R = []
    V = []
    with open(filename) as f:
        line = f.readline()
        while len(line.strip()) != 0:
            R.append(list(map(float, line.split())))
            line = f.readline()

        line = f.readline()
        while len(line.strip()) != 0:
            V.append(list(map(float, line.split())))
            line = f.readline()

        line = f.readline()
        t = list(map(float, line.split())) 
        f.close()
    return np.array(t), np.array(R), np.array(V)

t, R, V = get_data("P1_t50.txt")
z = R[2]
plt.plot(t, z)
plt.xlabel(r"t[$\mu s$]")
plt.ylabel(r"z[$\mu m$]")
plt.savefig(fig_path+"z.pdf")
plt.show()

t, R, V = get_data("P2_t50_PIfalse.txt")
x1, y1 = R[0], R[1]
x2, y2 = R[3], R[4]

plt.plot(x1, y1)
plt.plot(x2, y2)
plt.title("Without particle interactions")
plt.xlabel(r"x[$\mu m$]")
plt.ylabel(r"y[$\mu m$]")
plt.legend(["Particle 1", "Particle 2"])
plt.axis("equal")
plt.savefig(fig_path+"P2_t50_PIfalse.pdf")
plt.show()

z1, z2 = R[2], R[5]
vx1, vx2 = V[0], V[3]
vz1, vz2 = V[2], V[5]
plt.plot(x1, vx1, "--")

plt.plot(x2, vx2)
plt.title("Without particle interactions")
plt.xlabel(r"x[$\mu m$]")
plt.ylabel(r"dx/dt[$\mu m/\mu s$]")
plt.legend(["Particle 1", "Particle 2"])
plt.savefig(fig_path+"PhaseSpace_x_PIfalse.pdf")
plt.show()

plt.plot(z1, vz1, "--")
plt.plot(z2, vz2)
plt.title("Without particle interactions")
plt.xlabel(r"z[$\mu m$]")
plt.ylabel(r"dz/dt[$\mu m/\mu s}$]")
plt.legend(["Particle 1", "Particle 2"])
plt.savefig(fig_path+"PhaseSpace_z_PIfalse.pdf")
plt.show()

ax = plt.axes(projection="3d")
ax.plot3D(x1, y1, z1, label = "Particle 1")
ax.plot3D(x2, y2, z2, label = "Particle 2")
ax.set_xlabel(r"x[$\mu m$]")
ax.set_ylabel(r"y[$\mu m$]")
ax.set_zlabel(r"z[$\mu m$]")
ax.set_title("Without particle interactions")
ax.legend()
plt.savefig(fig_path+"3D_PIfalse.pdf")
plt.show()

t, R, V = get_data("P2_t50_PItrue.txt")
x1, y1 = R[0], R[1]
x2, y2 = R[3], R[4]

plt.plot(x1, y1)
plt.plot(x2, y2)
plt.title("With particle interactions")
plt.xlabel(r"x[$\mu m$]")
plt.ylabel(r"y[$\mu m$]")
plt.legend(["Particle 1", "Particle 2"])
plt.axis("equal")
plt.savefig(fig_path+"P2_t50_PItrue.pdf")
plt.show()

z1, z2 = R[2], R[5]
vx1, vx2 = V[0], V[3]
vz1, vz2 = V[2], V[5]

plt.plot(x1, vx1, "--")
plt.plot(x2, vx2)
plt.title("With particle interactions")
plt.xlabel(r"x[$\mu m$]")
plt.ylabel(r"dx/dt[$\mu m/\mu s$]")
plt.legend(["Particle 1", "Particle 2"])
plt.savefig(fig_path+"PhaseSpace_x_PItrue.pdf")
plt.show()

plt.plot(z1, vz1, "--")
plt.plot(z2, vz2)
plt.title("With particle interactions")
plt.xlabel(r"z[$\mu m$]")
plt.ylabel(r"dz/dt[$\mu m/\mu s}$]")
plt.legend(["Particle 1", "Particle 2"])
plt.savefig(fig_path+"PhaseSpace_z_PItrue.pdf")
plt.show()

ax = plt.axes(projection="3d")
ax.plot3D(x1, y1, z1, label = "Particle 1")
ax.plot3D(x2, y2, z2, label = "Particle 2")
ax.set_xlabel(r"x[$\mu m$]")
ax.set_ylabel(r"y[$\mu m$]")
ax.set_zlabel(r"z[$\mu m$]")
ax.set_title("With particle interactions")
ax.legend()
plt.savefig(fig_path+"3D_PItrue.pdf")
plt.show()

# Analtical solution for Paricle 1 equations of motion 
def analytic(t):
    v0 = 25 
    x0 = 20
    z0 = 20
    
    B_0 = 9.65e+1
    V_0 = 2.41e+6
    d = 500
    q = 1 
    m = 40.078 

    w_0 = q*B_0/m 
    w_z = np.sqrt(2*q*V_0/(m*d**2))
    w_minus = (w_0 - np.sqrt(w_0**2-2*w_z**2))/2 
    w_plus = (w_0 + np.sqrt(w_0**2-2*w_z**2))/2 

    A_plus = (v0 + w_minus*x0)/(w_minus - w_plus)
    A_minus = -(v0 + w_plus*x0)/(w_minus - w_plus)
    f = A_plus * np.exp(-1j*(w_plus*t)) + A_minus * np.exp(-1j*(w_minus*t))
    x = np.real(f)
    y = np.imag(f)
    z = z0 * np.cos(w_z*t)
    R = np.array([x, y, z])
    return R 

def relative_error(R_analytic, R_approximation):
    rel_err = np.linalg.norm((R_analytic - R_approximation), axis = 0)/np.linalg.norm(R_analytic, axis=0)
    return rel_err 

# Error convergence rate for simulations of 50 microseconds
def r_err(R_analytic_list, R_approximation_list):
    r_err = 0
    for k in range(1, len(R_analytic_list)):
        delta = np.max(np.linalg.norm(R_analytic_list[k] - R_approximation_list[k], axis = 0))
        delta_ = np.max(np.linalg.norm(R_analytic_list[k-1] - R_approximation_list[k-1], axis = 0))
        h = 50/(len(R_analytic_list[k][0])-1)
        h_ = 50/(len(R_analytic_list[k-1][0])-1)
        
        r_err += 1/3*np.log(delta/delta_)/np.log(h/h_)
    return r_err


t0, R0, V = get_data("P1_RK4_t50_n4000.txt")
R_0 = analytic(t0)
err0 = relative_error(R_0, R0)


t1, R1, V = get_data("P1_RK4_t50_n8000.txt")
R_1 = analytic(t1)

err1 = relative_error(R_1, R1)

t2, R2, V = get_data("P1_RK4_t50_n16000.txt")
R_2 = analytic(t2)

err2 = relative_error(R_2, R2)


t3, R3, V = get_data("P1_RK4_t50_n32000.txt")
R_3 = analytic(t3)

err3 = relative_error(R_3, R3)

plt.plot(t0, err0)
plt.plot(t1, err1)
plt.plot(t2, err2)
plt.plot(t3, err3)
plt.title("RK4")
plt.xlabel(r"t[$\mu$s]")
plt.ylabel("Realative error")
plt.legend(["n=4000", "n=8000", "n=16000", "n=32000"])
plt.savefig(fig_path+"RE_RK4.pdf")
plt.show()

plt.plot(t2, err2)
plt.plot(t3, err3)
plt.title("RK4")
plt.xlabel(r"t[$\mu$s]")
plt.ylabel("Realative error")
plt.legend(["n=16000", "n=32000"])
plt.savefig(fig_path+"RE_RK4_2.pdf")
plt.show()

RK4_r_err = r_err([R_0, R_1, R_2, R_3], [R0, R1, R2, R3])
print(f"RK4 error convergence rate = {RK4_r_err:.3f}")

t0, R0, V = get_data("P1_FE_t50_n4000.txt")
R_0 = analytic(t0)

err0 = relative_error(R_0, R0)

t1, R1, V = get_data("P1_FE_t50_n8000.txt")
R_1 = analytic(t1)

err1 = relative_error(R_1, R1)

t2, R2, V = get_data("P1_FE_t50_n16000.txt")
R_2 = analytic(t2)

err2 = relative_error(R_2, R2)

t3, R3, V = get_data("P1_FE_t50_n32000.txt")
R_3 = analytic(t3)

err3 = relative_error(R_3, R3)

plt.plot(t0, err0)
plt.plot(t1, err1)
plt.plot(t2, err2)
plt.plot(t3, err3)
plt.title("Forward Euler")
plt.xlabel(r"t[$\mu$s]")
plt.ylabel("Realative error")
plt.legend(["n=4000", "n=8000", "n=16000", "n=32000"])
plt.savefig(fig_path+"RE_FE.pdf")
plt.show()

FE_r_err = r_err([R_0, R_1, R_2, R_3], [R0, R1, R2, R3])
print(f"Forward Euler error convergence rate = {FE_r_err:.3f}")

w, particcles_left, f = get_data("Particles_left.txt")
f = f.ravel()

plt.plot(w, particcles_left[0], ".-", label = f"f={f[0]}") 
plt.plot(w, particcles_left[1], ".-", label = f"f={f[1]}")
plt.plot(w, particcles_left[2], ".-", label = f"f={f[2]}")
plt.legend()
plt.title("Without particle interactions")
plt.xlabel("Frequency [MHz]")
plt.ylabel("Particles in trap")
plt.savefig(fig_path+"Particles_left.pdf")
plt.show()


w, particcles_left, f = get_data("Particles_left_Coulomb_off.txt")
f = f.ravel()

plt.plot(w, particcles_left[0], ".-", label = f"f={f[0]}") 
plt.plot(w, particcles_left[1], ".-", label = f"f={f[1]}")
plt.plot(w, particcles_left[2], ".-", label = f"f={f[2]}")
plt.legend()
plt.title("Without particle interactions")
plt.xlabel("Frequency [MHz]")
plt.ylabel("Particles in trap")
plt.savefig(fig_path+"Particles_left_Coulomb_off.pdf")
plt.show()

w, particcles_left, f = get_data("Particles_left_Coulomb_on.txt")
f = f.ravel()

plt.plot(w, particcles_left[0], ".-", label = f"f={f[0]}") 
plt.plot(w, particcles_left[1], ".-", label = f"f={f[1]}")
plt.plot(w, particcles_left[2], ".-", label = f"f={f[2]}")
plt.legend()
plt.title("With particle interactions")
plt.xlabel("Frequency [MHz]")
plt.ylabel("Particles in trap")
plt.savefig(fig_path+"Particles_left_Coulomb_on.pdf")
plt.show()


