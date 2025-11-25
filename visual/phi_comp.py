import numpy as np 
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt

"""
The columns in data are as follows:

r, theta, tp (= t), Re(eff_src), Im(eff_src),

The src the effective source (box of puncture)
"""

m=0
data = {}
data_phi ={}
for i in np.arange(301):
    data[i] = np.loadtxt("/home/cillian/.effective_source/data0/loop_%d.txt" %i)
    data_phi[i] = np.loadtxt("/home/cillian/.effective_source/data_phi/loop_%d.txt" %i)

time = np.zeros(len(data))
re_source = np.zeros(len(data))
im_source = np.zeros(len(data))
for i in np.arange(len(data)):
    time[i] = data_phi[i][2]
    re_source[i] = data_phi[i][3]
    im_source[i] = data_phi[i][4]

re_source_wrong = np.zeros(len(data))
im_source_wrong = np.zeros(len(data))
for i in np.arange(len(data)):
    time[i] = data[i][2]
    re_source_wrong[i] = data[i][3]
    im_source_wrong[i] = data[i][4]

mod_source = np.sqrt(re_source**2 + im_source**2)

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(time, mod_source)
#ax1.plot(time, re_source)
#ax1.plot(time, im_source)

ax1.vlines(x=[140.238,280.476], ymin = min(mod_source), ymax = max(mod_source), color = "red")
ax1.set_ylabel("$|S_{eff}^{m=%d}|$" %m, fontsize=15)
ax1.set_xlabel("t")
ax1.set_title(r"$|S_{eff}^{m=%d}|$ as a function of time for $(r,\theta) = (%d,\pi/2)$" %(m,int(data[0][0])))
ax1.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

plt.savefig("Eff_source_m%d_r%d_t" %(m, int(data[0][0])))
plt.show()

"""
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
#ax2.plot(time, np.sqrt(re_source_wrong**2 + im_source_wrong**2))
ax2.plot(time, re_source_wrong)
#ax2.plot(time, im_source_wrong)
ax2.set_title("No phi evolution")
ax2.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
fig2.savefig("./phi_ev_comparison/no_phi_evolution_601")

"""

"""
# Create dictionaries to store data at different times
r = {}
theta = {}
time = np.zeros(len(data))
Re_effsource = {}
Im_effsource = {}
for i in np.arange(len(data)):
    r[i] = data[i][:,0]
    theta[i] = data[i][:,1]
    time[i] = data[i][0,2]
    Re_effsource[i] = data[i][:,3]
    Im_effsource[i] = data[i][:,4]


print(r[0])
# Create solution grid
r_grid = {}
theta_grid = {}
Re_effsource_grid = {}
Im_effsource_grid = {}
for i in np.arange(len(data)):
    r_grid[i], theta_grid[i] = np.meshgrid(np.unique(r[i]), np.unique(theta[i]))
    Re_effsource_grid[i] = Re_effsource[i].reshape(r_grid[i].shape) 
    Im_effsource_grid[i] = Im_effsource[i].reshape(r_grid[i].shape) 

for i in np.arange(len(data)):
"""
"""
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.plot_surface(r_grid[i], theta_grid[i], Re_effsource_grid[i])
    ax.set_xlabel("$r$",fontsize = 15)
    ax.set_ylabel(r"$\theta$", fontsize = 15)
    ax.set_zlabel("Re :$S_{eff}^{m = 2}$", fontsize = 10)
    ax.set_title("t=%f" %time[i])
    ax.view_init(azim=-45)

    #plt.show()
    """

"""
#print(r_grid[5][0,r_index])
print(theta_grid[5][theta_index, 0] - np.pi/2)

Im_effsource_t = np.zeros(len(data))
Re_effsource_t = np.zeros(len(data))
for i in np.arange(len(data)):
    Re_effsource_t[i] = Re_effsource_grid[i][theta_index, r_index]
    Im_effsource_t[i] = Im_effsource_grid[i][theta_index, r_index]

fig1 = plt.figure()
#fig2 = plt.figure()
ax1 = fig1.add_subplot(111)
#ax2 = fig2.add_subplot(111)
ax1.plot(time, np.sqrt(Re_effsource_t**2 + Im_effsource_t**2))
#plt.show()
print()
"""
"""
ax2.plot(time, Re_effsource_t)
ax1.set_xlabel("t", fontsize = 15)
ax1.set_ylabel("$Im(S_{eff}^{m = 2})$",fontsize = 15)
ax2.set_xlabel("t", fontsize = 15)
ax2.set_ylabel("$Re(S_{eff}^{m = 2})$",fontsize = 15)
#ax1.set_xlim(100,200)
plt.show()
"""
"""
r = data[150][:,0] # + 10
theta = data[150][:,1]# + np.pi/2
t = data[150][:,2]
ReSrc = data[150][:, 3]
ImSrc = data[150][:, 4]
print(t[0])

eff_src = np.sqrt(ReSrc**2 + ImSrc**2)

r, theta = np.meshgrid(np.unique(theta), np.unique(r))
ReSrc = ReSrc.reshape(r.shape) 
ImSrc = ImSrc.reshape(r.shape) 
eff_src = eff_src.reshape(r.shape)

#print(ReSrc)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot_surface(r, theta, ReSrc)
ax.set_ylabel("$\delta r$",fontsize = 15)
ax.set_xlabel(r"$\delta \theta$", fontsize = 15)
ax.set_zlabel("Re :$S_{eff}^{m = 2}$", fontsize = 10)
ax.set_title("t=0")
ax.view_init(azim=-45)

#plt.savefig("./phi_0")
plt.show()
"""
