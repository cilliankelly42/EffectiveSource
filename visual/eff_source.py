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
data= {}
for i in np.arange(601):
    data[i] = np.loadtxt("/home/cillian/.effective_source/datam0/loop_%d.txt" %i)

"""
for i in np.arange(len(datam1)):
    print(datam0[i][2] - datam1[i][2])

time = np.zeros(len(datam0))
re_source_m0 = np.zeros(len(datam0))
im_source_m0 = np.zeros(len(datam0))
for i in np.arange(len(datam0)):
    time[i] = datam0[i][2]
    re_source_m0[i] = datam0[i][3]
    im_source_m0[i] = datam0[i][4]

re_source_m1 = np.zeros(len(datam1))
im_source_m1 = np.zeros(len(datam1))
for i in np.arange(len(datam1)):
    re_source_m1[i] = datam1[i][3]
    im_source_m1[i] = datam1[i][4]

mod_m0 = np.sqrt(re_source_m0**2 + im_source_m0**2)
mod_m1 = np.sqrt(re_source_m1**2 + im_source_m1**2)
#print(mod_m0 - mod_m1)

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(time, mod_m0)
ax1.plot(time, mod_m1)
plt.show()
"""

time = np.zeros(len(data))
re_source = np.zeros(len(data))
im_source = np.zeros(len(data))
for i in np.arange(len(data)):
    time[i] = data[i][2]
    re_source[i] = data[i][3]
    im_source[i] = data[i][4]

mod_source = np.sqrt(re_source**2 + im_source**2)

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(time, mod_source)
#ax1.plot(time, re_source)
#ax1.plot(time, im_source)
#ax1.vlines(x=[140.238,280.476,420.714], ymin = min(mod_source), ymax = max(mod_source), color = "red")

ax1.set_ylabel("$Im(S_{eff}^{m=%d})$" %m, fontsize=15)
ax1.set_xlabel("t")
#ax1.set_title(r"$|S_{eff}^{m=%d}|$ as a function of time for $(r,\theta) = (%d,%d)$" %(m,int(data[0][0]),int(data[0][1])))
ax1.set_title(r"$(r,\theta) = (%d,%d)$" %(int(data[0][0]),int(data[0][1])))
ax1.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

print(data[2][2] - data[1][2])
#plt.savefig("Im_Eff_source_m%d_r%d_theta%d" %(m, int(data[0][0]), int(data[0][1])))
plt.show()