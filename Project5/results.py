import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib
from matplotlib.animation import FuncAnimation

work_path = os.path.dirname(__file__) 
fig_path = os.path.join(work_path, "fig/") #path to directory storing figures

plt.rcParams["text.usetex"] =True

def get_data(filename):
    with open(filename) as f:
        U = []
        t = f.readline().split() 
        positions = f.readline().split() 
        for line in f.readlines():
            cx_row = []
            for cx in line.split():
                # Create string of complex number compatible with complex()
                cx = cx.strip("()").replace(",", "+").replace("+-","-")+"j"
                cx_row.append(complex(cx))
            U.append(cx_row)
        f.close()
        return np.array(t, dtype="float"), np.array(positions, dtype="float"), np.array(U)

# Problem 7  
#"""

t, positions, U0 = get_data("problem7a.txt")
U1 = get_data("problem7b.txt")[2]

N = len(t)
tot_prob_deviation0 = np.zeros(N)
tot_prob_deviation1 = np.zeros(N)
for i in range(N):
    tot_prob_deviation0[i] = 1 - U0[:,i] @ np.conjugate(U0[:,i]).T
    tot_prob_deviation1[i] = 1 - U1[:,i] @ np.conjugate(U1[:,i]).T

plt.figure()
plt.plot(t, tot_prob_deviation0)
plt.title(r"$v_0=0$")
plt.xlabel("t")
plt.ylabel("probability deviation from 1")
plt.savefig("problem7a.pdf")

plt.figure()
plt.plot(t, tot_prob_deviation1)
plt.title(r"$v_0=10^{10}$")
plt.xlabel("t")
plt.ylabel("probability deviation from 1")
plt.savefig("problem7b.pdf")
plt.show()
#"""

# Problem 8
#"""
t, positions, U = get_data("problem8.txt")
u0 = U[:,0]
u1 = U[:, np.where(t==0.001)[0][0]]
u2 = U[:, np.where(t==0.002)[0][0]]

N = int(np.sqrt(len(u0)))
U0 = (np.conjugate(u0)*u0).reshape(N, N)
U1 = (np.conjugate(u1)*u1).reshape(N, N)
U2 = (np.conjugate(u2)*u2).reshape(N, N)

u0_Re, u0_Im = np.real(u0).reshape(N, N), np.imag(u0).reshape(N, N)
u1_Re, u1_Im = np.real(u1).reshape(N, N), np.imag(u1).reshape(N, N)
u2_Re, u2_Im = np.real(u2).reshape(N, N), np.imag(u2).reshape(N, N)


x_values = positions
y_values = positions
x, y = np.meshgrid(x_values, y_values)


extent = np.min(positions), np.max(positions), np.min(positions), np.max(positions)


plt.figure()
plt.title("p(t = 0)")
plt.imshow(np.real(U0), extent = extent, vmax = np.max(np.real(U0)))
plt.colorbar()
plt.xlabel("x")
plt.ylabel("y")
plt.savefig("p_heatmap_0.pdf")

plt.figure()
plt.title("p(t = 0.001)")
plt.imshow(np.real(U1), extent = extent, vmax = np.max(np.real(U1)))
plt.colorbar()
plt.xlabel("x")
plt.ylabel("y")
plt.savefig("p_heatmap_1.pdf")

plt.figure()
plt.title("p(t = 0.002)")
plt.imshow(np.real(U2), extent = extent, vmax = np.max(np.real(U2)))
plt.colorbar()
plt.xlabel("x")
plt.ylabel("y")
plt.savefig("p_heatmap_2.pdf")
plt.show()


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 6))
fig.suptitle("u(t = 0)", fontsize = 15)
im1 = ax1.imshow(u0_Re, extent = extent, vmax = np.max(u0_Re))
im2 = ax2.imshow(u0_Im, extent = extent, vmax = np.max(u0_Im))
divider1 = make_axes_locatable(ax1)
divider2 = make_axes_locatable(ax2)
cax = divider1.append_axes('bottom', size='5%', pad=0.5) 
fig.colorbar(im1, cax=cax, orientation='horizontal', label="Real(u)")
cax = divider2.append_axes('bottom', size='5%', pad=0.5) 
fig.colorbar(im2, cax=cax, orientation='horizontal', label="Imag(u)")
ax1.set(xlabel = "x", ylabel = "y")
ax2.set(xlabel = "x", ylabel = "y")
plt.subplots_adjust(left=0, right=1, top=0.9, bottom=0.1, wspace=-0.1)
plt.savefig("u_heatmap_0.pdf")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 6))
fig.suptitle("u(t = 0.001)", fontsize = 15)
im1 = ax1.imshow(u1_Re, extent = extent, vmax = np.max(u1_Re))
im2 = ax2.imshow(u1_Im, extent = extent, vmax = np.max(u1_Im))
divider1 = make_axes_locatable(ax1)
divider2 = make_axes_locatable(ax2)
cax = divider1.append_axes('bottom', size='5%', pad=0.5) 
fig.colorbar(im1, cax=cax, orientation='horizontal', label="Real(u)")
cax = divider2.append_axes('bottom', size='5%', pad=0.5) 
fig.colorbar(im2, cax=cax, orientation='horizontal', label="Imag(u)")
ax1.set(xlabel = "x", ylabel = "y")
ax2.set(xlabel = "x", ylabel = "y")
plt.subplots_adjust(left=0, right=1, top=0.9, bottom=0.1, wspace=-0.1)
plt.savefig("u_heatmap_1.pdf")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 6))
fig.suptitle("u(t = 0.002)", fontsize = 15)
im1 = ax1.imshow(u2_Re, extent = extent, vmax = np.max(u2_Re))
im2 = ax2.imshow(u2_Im, extent = extent, vmax = np.max(u2_Im))
divider1 = make_axes_locatable(ax1)
divider2 = make_axes_locatable(ax2)
cax = divider1.append_axes('bottom', size='5%', pad=0.5) 
fig.colorbar(im1, cax=cax, orientation='horizontal', label="Real(u)")
cax = divider2.append_axes('bottom', size='5%', pad=0.5) 
fig.colorbar(im2, cax=cax, orientation='horizontal', label="Imag(u)")
ax1.set(xlabel = "x", ylabel = "y")
ax2.set(xlabel = "x", ylabel = "y")
plt.subplots_adjust(left=0, right=1, top=0.9, bottom=0.1, wspace=-0.1)
plt.savefig("u_heatmap_2.pdf")

plt.show()
#"""

# Problem 9
#"""
def probability_detector(U, n, m, j):
    """ Adding the probability at x[j] over time """
    p = np.zeros(m)
    for i in range(n):
        P_i = np.real(np.conjugate(U[:,i])*U[:,i]).reshape(m, m)
        p += P_i[:, j]
        if any(elem > 1e-8 for elem in P_i[:, -1]):
            print("Warning! Reflection of back wall")
    p = p / np.linalg.norm(p)
    return p


t, positions, U_double = get_data("problem8.txt")
t, positions, U_single = get_data("problem9_single.txt")
t, positions, U_triple = get_data("problem9_triple.txt")

n = len(t)
m = int(np.sqrt(len(U_double[:,0])))
j = np.where(positions == 0.8)[0][0]

p_double = probability_detector(U_double, n, m, j)
p_single = probability_detector(U_single, n, m, j)
p_triple = probability_detector(U_triple, n, m, j)

plt.figure()
plt.plot(positions, p_double)
plt.title("Double slit")
plt.xlabel("y")
plt.ylabel(r"p(y$|$x=0.8;t=0.002)")
plt.savefig("p_double.pdf")

plt.figure()
plt.plot(positions, p_single)
plt.title("Single slit")
plt.xlabel("y")
plt.ylabel(r"p(y$|$x=0.8;t=0.002)")
plt.savefig("p_single.pdf")

plt.figure()
plt.plot(positions, p_triple)
plt.title("Triple slit")
plt.xlabel("y")
plt.ylabel(r"p(y$|$x=0.8;t=0.002)")
plt.savefig("p_triple.pdf")

plt.show()
#"""

# Problem X
#"""

t, positions, U = get_data("problem8.txt")
extent = np.min(positions), np.max(positions), np.min(positions), np.max(positions)
N = int(np.sqrt(len(U[:,0])))

dt = t[1]

# Fill z_data_list with U
z_data_list = []
for i in range(len(t)):
    z_data = np.real(np.conjugate(U[:,i])*U[:,i]).reshape(N, N)
    z_data_list.append(z_data)


#
# Now the list z_data_list contains a series of "frames" of z(x,y,t), 
# where each frame can be plotted as a 2D image using imshow. Let's
# animate it!
#

# Some settings
fontsize = 12
t_min = t[0]

# Create figure
fig = plt.figure()
ax = plt.gca()

# Create a colour scale normalization according to the max z value in the first frame
norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(z_data_list[0]))

# Plot the first frame
img = ax.imshow(z_data_list[0], extent=extent, cmap=plt.get_cmap("viridis"), norm=norm)

# Axis labels
plt.xlabel("x", fontsize=fontsize)
plt.ylabel("y", fontsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)

# Add a colourbar
cbar = fig.colorbar(img, ax=ax)
cbar.set_label("p(x,y;t)", fontsize=fontsize)
cbar.ax.tick_params(labelsize=fontsize)

# Add a text element showing the time
time_txt = plt.text(0.95, 0.95, "t = {:.3e}".format(t_min), color="white", 
                    horizontalalignment="right", verticalalignment="top", fontsize=fontsize)

# Function that takes care of updating the z data and other things for each frame
def animation(i):
    # Normalize the colour scale to the current frame?
    norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(z_data_list[i]))
    img.set_norm(norm)

    # Update z data
    img.set_data(z_data_list[i])

    # Update the time label
    current_time = t_min + i * dt
    time_txt.set_text("t = {:.3e}".format(current_time))

    return img

# Use matplotlib.animation.FuncAnimation to put it all together
anim = FuncAnimation(fig, animation, interval=1, frames=np.arange(0, len(z_data_list), 2), repeat=False, blit=0)

# Run the animation!
plt.show()

# # Save the animation
anim.save('animation_8.gif', writer="writegif", bitrate=-1, fps=20)

#"""
