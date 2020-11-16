# Script to plot chemical shifts of imino protons and nitrogens from BMRB chemical shift data
# Input "nuc_sparse.txt" and outputs .eps-plot

import numpy as np
import matplotlib

matplotlib.use('TKAgg',warn=False, force=True) # matplotlib.use('Agg')
from matplotlib import pyplot as plt
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)

# 2JTP 15417 1 G N G C s H1 12.0

inpf = open("../lists/imino_nuc_sparse_20.txt","r")
lines = inpf.readlines()

i = 0

# UG
ug_H3 = []
ug_N3 = []

# UA
ua_H3 = []
ua_N3 = []

# GC
gc_H1 = []
gc_N1 = []

# GU
gu_H1 = []
gu_N1 = []

while i < (len(lines)-1):

    current_line = lines[i]
    split_curr = current_line.split()

    # 1269 2QH3 7404 2 G G A C s H1 N1 12.0 146.0

    num, c_pdb_id, c_bmrbid, c_pos, c_base, c_basebf, c_basea, c_base_part, c_feat, c_H, c_N, c_Hcs, c_Ncs = split_curr


    if (c_base == "U" and c_base_part == "G") and (c_feat == "s" and c_H == "H3"):

        ug_H3.append(float(c_Hcs))
        ug_N3.append(float(c_Ncs))

    elif (c_base == "U" and c_base_part == "A") and (c_feat == "s" and c_H == "H3"):

        ua_H3.append(float(c_Hcs))
        ua_N3.append(float(c_Ncs))

    elif (c_base == "G" and c_base_part == "C") and (c_feat == "s" and c_H == "H1"):

        gc_H1.append(float(c_Hcs))
        gc_N1.append(float(c_Ncs))

    elif (c_base == "G" and c_base_part == "U") and (c_feat == "s" and c_H == "H1"):

        gu_H1.append(float(c_Hcs))
        gu_N1.append(float(c_Ncs))

    i+=1

inpf.close()


# Average and std U-G
avg_ug_H3 = round(np.average(np.asarray(ug_H3)),3) # number is number of significant figures
avg_ug_N3 = round(np.average(np.asarray(ug_N3)),3)

std_ug_H3 = round(np.std(np.asarray(ug_H3)),3)
std_ug_N3 = round(np.std(np.asarray(ug_N3)),3)

# Average and std U-A

avg_ua_H3 = round(np.average(np.asarray(ua_H3)),3) # number is number of significant figures
avg_ua_N3 = round(np.average(np.asarray(ua_N3)),3)

std_ua_H3 = round(np.std(np.asarray(ua_H3)),3)
std_ua_N3 = round(np.std(np.asarray(ua_N3)),3)

# Average and std G-C
avg_gc_H1 = round(np.average(np.asarray(gc_H1)),3) # number is number of significant figures
avg_gc_N1 = round(np.average(np.asarray(gc_N1)),3)

std_gc_H1 = round(np.std(np.asarray(gc_H1)),3)
std_gc_N1 = round(np.std(np.asarray(gc_N1)),3)

# Average and std G-U
avg_gu_H1 = round(np.average(np.asarray(gu_H1)),3) # number is number of significant figures
avg_gu_N1 = round(np.average(np.asarray(gu_N1)),3)

std_gu_H1 = round(np.std(np.asarray(gu_H1)),3)
std_gu_N1 = round(np.std(np.asarray(gu_N1)),3)


# plot
fig = plt.figure()
params = {'legend.fontsize': 20,
          'legend.handlelength': 2}
plt.rcParams.update(params)
plt.plot(ug_H3, ug_N3, "mo", markersize=6, zorder=1, label="UG n="+str(len(ug_H3)))
plt.errorbar(avg_ug_H3,avg_ug_N3,std_ug_N3,std_ug_H3,"kx",markersize=6)
plt.plot(ua_H3, ua_N3, "bo", markersize=6, zorder=1, label="UA n="+str(len(ua_H3)))
plt.errorbar(avg_ua_H3,avg_ua_N3,std_ua_N3,std_ua_H3,"kx",markersize=6)
plt.plot(gc_H1, gc_N1, "co", markersize=6, zorder=1, label="GC n="+str(len(gc_H1)))
plt.errorbar(avg_gc_H1,avg_gc_N1,std_gc_N1,std_gc_H1,"kx",markersize=6)
plt.plot(gu_H1, gu_N1, "yo", markersize=6, zorder=1, label="GU n="+str(len(gu_H1)))
plt.errorbar(avg_gu_H1,avg_gu_N1,std_gu_N1,std_gu_H1,"kx",markersize=6)




# Adding ellipses with 1 sd
# GU
gu_horizontal_radius = std_gu_H1
gu_vertical_radius = std_gu_N1
gu_x_center = avg_gu_H1
gu_y_center = avg_gu_N1

gu_x = np.linspace(0,15,1000)
gu_y = np.linspace(135,150,1000)

gu_X,gu_Y = np.meshgrid(gu_x,gu_y)

gu_eqn = (((gu_X-gu_x_center)**2)/((gu_horizontal_radius)**2))+ (((gu_Y-gu_y_center)**2)/((gu_vertical_radius)**2))
Z = 1

plt.contour(gu_X,gu_Y,gu_eqn,[Z], colors="k", linewidths=1, linestyles="dashed")

# GC
gc_horizontal_radius = std_gc_H1
gc_vertical_radius = std_gc_N1
gc_x_center = avg_gc_H1
gc_y_center = avg_gc_N1

gc_x = np.linspace(0,15,1000)
gc_y = np.linspace(135,155,1000)

gc_X,gc_Y = np.meshgrid(gc_x,gc_y)

gc_eqn = (((gc_X-gc_x_center)**2)/((gc_horizontal_radius)**2))+ (((gc_Y-gc_y_center)**2)/((gc_vertical_radius)**2))
Z = 1

plt.contour(gc_X,gc_Y,gc_eqn,[Z], colors="k", linewidths=1, linestyles="dashed")

# UG
ug_horizontal_radius = std_ug_H3
ug_vertical_radius = std_ug_N3
ug_x_center = avg_ug_H3
ug_y_center = avg_ug_N3

ug_x = np.linspace(0,15,1000)
ug_y = np.linspace(150,165,1000)

ug_X,ug_Y = np.meshgrid(ug_x,ug_y)

ug_eqn = (((ug_X-ug_x_center)**2)/((ug_horizontal_radius)**2))+ (((ug_Y-ug_y_center)**2)/((ug_vertical_radius)**2))
Z = 1

plt.contour(ug_X,ug_Y,ug_eqn,[Z], colors="k", linewidths=1, linestyles="dashed")

# UA
ua_horizontal_radius = std_ua_H3
ua_vertical_radius = std_ua_N3
ua_x_center = avg_ua_H3
ua_y_center = avg_ua_N3

ua_x = np.linspace(0,15,1000)
ua_y = np.linspace(150,165,1000)

ua_X,ua_Y = np.meshgrid(ua_x,ua_y)

ua_eqn = (((ua_X-ua_x_center)**2)/((ua_horizontal_radius)**2))+ (((ua_Y-ua_y_center)**2)/((ua_vertical_radius)**2))
Z = 1

plt.contour(ua_X,ua_Y,ua_eqn,[Z], colors="k", linewidths=1, linestyles="dashed")


plt.title("Imino chemical shifts", fontsize=20)
plt.xlabel("1H (ppm)", fontsize=20)
plt.ylabel("15N (ppm)", fontsize=20)
plt.legend()

plt.xlim(15.0,10.0)
plt.ylim(170.0,135.0)
plt.show()

fig.savefig("bmrb_imino_20.eps",format='eps')


