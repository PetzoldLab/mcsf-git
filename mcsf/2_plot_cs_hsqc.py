#!/usr/local/bin/python
# coding: utf-8

# **************************************************************************** #
# ||                    ~~~  2_plot_cs_hsqc.py  ~~~                         || #
# ||                                                                        || #
# ||           Script to plot chemical shifts of imino protons              || #
# ||             and nitrogens from BMRB chemical shift data.               || #
# ||                                                                        || #
# ||                     Written by Hampus Karlsson May-18                  || #
# ||                Updated by Noah Hopkins and Magdalena Riad              || #
# ||                           for the publication                          || #
# ||      "Mutate-and-Chemical-Shift-Fingerprint (MCSF) to characterize     || #
# ||              excited states in RNA using NMR spectroscopy"             || #
# ||                                                                        || #
# ||                           Katja Petzold Group                          || #
# ||                                May 2021                                || #
# **************************************************************************** #

# Script to create 2D plots using a list of chemical shifts.
#
# Input: "out_ALL_(X1[Y1]_A1-B1)_(X2[Y2]_A2-B2)_ ... _(Xn[Yn]_An-Bn).txt"
#
# Outputs a plot in .eps or .pdf format.
#
# Note: This script is not fully automated and needs to be adjusted/configured
# according to ones needs. Currently it is set up using one of the MCSF
# publication figures as an example.

# **********--------------------------------------------------------********** #
# |                           ~~ [1] IMPORTS ~~                              | #
# **********--------------------------------------------------------********** #

# Import built-ins
import os.path

# Import third party modules / packages
import numpy as np
import matplotlib
from platform import python_version

if float(python_version()[0]) >= 3.0:
    from matplotlib import pyplot as plt
else:
    matplotlib.use('TKAgg', warn=False, force=True)
    from matplotlib import pyplot as plt


# **********--------------------------------------------------------********** #
# |                    ~~ [2] OPTIONS AND INPUT DATA ~~                      | #
# **********--------------------------------------------------------********** #

INPUT_FILE = 'out/out_ALL_(A[U]_H3-N3)_(G[U]_H3-N3)_(U[G]_H1-N1)_(C[G]_H1-N1).txt'
# PATH TO INPUT FILE
# Input file is a list of shifts and related data. Can be retrieved by running
# script '1_extract_from_pair_list.py'.

SAVE_FIGURE = False
# SAVE OPTION
# Set to True to save the figure.

FIGURE_FORMAT = 'pdf'
# CHOOSE FORMAT
# Choose to save the figure in .eps or .pdf format.
# Set to either 'eps' or 'pdf'.

PAIR_TYPE = 'all'
# CHOOSE PAIR TYPE TO PLOT.
# Choose 'all' to plot all base-pairs.
# Choose 'brackets' to plot only base-pairs with brackets in the
# dot-bracket notation.

# Configure plot ticks and label sizes
plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=12)


# NOTE: Additional settings and configuration need to be made below
#       depending on your needs.


# **********--------------------------------------------------------********** #
# |                              ~~ [5] RUN ~~                               | #
# **********--------------------------------------------------------********** #

inpf = open(os.path.normpath(INPUT_FILE), "r")
lines = inpf.readlines()

i = 0

# GU
gu_H3 = []
gu_N3 = []

# AU
au_H3 = []
au_N3 = []

# CG
cg_H1 = []
cg_N1 = []

# UG
ug_H1 = []
ug_N1 = []

for i, line in enumerate(lines):
    current_line = line
    split_curr = current_line.split()
    
    # Columns
    # 0:   idx
    # 1:   pdb_id
    # 2:   entry_id
    # 3:   shift_list_id
    # 4:   assembly_id
    # 5:   entity_id
    # 6:   entity_assembly_id
    # 7:   residue_id
    # 8:   nt_type
    # 9:   nt_dbn
    # 10:  residue_id_bp_parner
    # 11:  nt_type_bp_partner
    # 12:  nt_dbn_bp_partner
    # 13:  nuc_type_h
    # 14:  cs_data_h
    # 15:  cs_err_data_h
    # 16:  ambiguity_code_h
    # 17:  nuc_type_n
    # 18:  cs_data_n
    # 19:  cs_err_data_n
    # 20:  ambiguity_code_n
    
    # Example row:
    # 1 1JO7 4816 1 1 1 1 28 U ) 4 A ( H3 13.648 0.002 1 N3 165.283 0.0 1
    (idx, pdb_id, bmrbid, shift_list_id, assembid, entityid,
     entity_assid, resid, nt_type, nt_dbn, resid_bp_part, nt_type_part,
     nt_dbn_part, c_nuctype_H, c_Hcs, c_Hcs_err, c_amb_h, c_nuctype_n,
     c_Ncs, c_Ncs_err, c_amb_n) = split_curr
    
    # Adjust/configure according to needs:
    if PAIR_TYPE == "all":
        if (nt_type == "U" and nt_type_part == "G") and c_nuctype_H == "H3":

            gu_H3.append(float(c_Hcs))
            gu_N3.append(float(c_Ncs))

        elif (nt_type == "U" and nt_type_part == "A") and c_nuctype_H == "H3":

            au_H3.append(float(c_Hcs))
            au_N3.append(float(c_Ncs))

        elif (nt_type == "G" and nt_type_part == "C") and c_nuctype_H == "H1":

            cg_H1.append(float(c_Hcs))
            cg_N1.append(float(c_Ncs))

        elif (nt_type == "G" and nt_type_part == "U") and c_nuctype_H == "H1":

            ug_H1.append(float(c_Hcs))
            ug_N1.append(float(c_Ncs))

    elif PAIR_TYPE == "brackets":
        if (nt_type == "U" and nt_type_part == "G") and (nt_dbn == "(" or nt_dbn == ")") and (nt_dbn_part == "(" or nt_dbn_part == ")") and c_nuctype_H == "H3":

            gu_H3.append(float(c_Hcs))
            gu_N3.append(float(c_Ncs))

        elif (nt_type == "U" and nt_type_part == "A") and (nt_dbn == "(" or nt_dbn == ")") and (nt_dbn_part == "(" or nt_dbn_part == ")") and c_nuctype_H == "H3":

            au_H3.append(float(c_Hcs))
            au_N3.append(float(c_Ncs))

        elif (nt_type == "G" and nt_type_part == "C") and (nt_dbn == "(" or nt_dbn == ")") and (nt_dbn_part == "(" or nt_dbn_part == ")") and c_nuctype_H == "H1":

            cg_H1.append(float(c_Hcs))
            cg_N1.append(float(c_Ncs))

        elif (nt_type == "G" and nt_type_part == "U") and (nt_dbn == "(" or nt_dbn == ")") and (nt_dbn_part == "(" or nt_dbn_part == ")") and c_nuctype_H == "H1":

            ug_H1.append(float(c_Hcs))
            ug_N1.append(float(c_Ncs))

    i += 1

inpf.close()

# Adjust/configure according to needs:

# Average and std U-G
avg_ug_H3 = round(np.average(np.asarray(gu_H3)), 3)
avg_ug_N3 = round(np.average(np.asarray(gu_N3)), 3)

std_ug_H3 = round(np.std(np.asarray(gu_H3)), 3)
std_ug_N3 = round(np.std(np.asarray(gu_N3)), 3)

# Average and std U-A
avg_ua_H3 = round(np.average(np.asarray(au_H3)), 3)
avg_ua_N3 = round(np.average(np.asarray(au_N3)), 3)

std_ua_H3 = round(np.std(np.asarray(au_H3)), 3)
std_ua_N3 = round(np.std(np.asarray(au_N3)), 3)

# Average and std G-C
avg_gc_H1 = round(np.average(np.asarray(cg_H1)), 3)
avg_gc_N1 = round(np.average(np.asarray(cg_N1)), 3)

std_gc_H1 = round(np.std(np.asarray(cg_H1)), 3)
std_gc_N1 = round(np.std(np.asarray(cg_N1)), 3)

# Average and std G-U
avg_gu_H1 = round(np.average(np.asarray(ug_H1)), 3)
avg_gu_N1 = round(np.average(np.asarray(ug_N1)), 3)

std_gu_H1 = round(np.std(np.asarray(ug_H1)), 3)
std_gu_N1 = round(np.std(np.asarray(ug_N1)), 3)

# Plot
fig = plt.figure()
params = {'legend.fontsize': 14,
          'legend.handlelength': 1}
plt.rcParams.update(params)
plt.plot(gu_H3, gu_N3, "mo", markersize=6, zorder=1, label="GU n=" + str(len(gu_H3)))
plt.errorbar(avg_ug_H3, avg_ug_N3, std_ug_N3, std_ug_H3, "kx", markersize=6)
plt.plot(au_H3, au_N3, "bo", markersize=6, zorder=1, label="AU n=" + str(len(au_H3)))
plt.errorbar(avg_ua_H3, avg_ua_N3, std_ua_N3, std_ua_H3, "kx", markersize=6)
plt.plot(cg_H1, cg_N1, "ro", markersize=6, zorder=1, label="CG n=" + str(len(cg_H1)))
plt.errorbar(avg_gc_H1, avg_gc_N1, std_gc_N1, std_gc_H1, "kx", markersize=6)
plt.plot(ug_H1, ug_N1, "co", markersize=6, zorder=1, label="UG n=" + str(len(ug_H1)))
plt.errorbar(avg_gu_H1, avg_gu_N1, std_gu_N1, std_gu_H1, "kx", markersize=6)

# Adding ellipses with 1 sd
# GU
gu_horizontal_radius = std_gu_H1
gu_vertical_radius = std_gu_N1
gu_x_center = avg_gu_H1
gu_y_center = avg_gu_N1

gu_x = np.linspace(0, 15, 1000)
gu_y = np.linspace(135, 150, 1000)

gu_X, gu_Y = np.meshgrid(gu_x, gu_y)

gu_eqn = (((gu_X - gu_x_center) ** 2) / (gu_horizontal_radius ** 2)) + (
            ((gu_Y - gu_y_center) ** 2) / (gu_vertical_radius ** 2))
Z = 1

plt.contour(gu_X, gu_Y, gu_eqn, [Z], colors="k", linewidths=1, linestyles="dashed")

# GC
gc_horizontal_radius = std_gc_H1
gc_vertical_radius = std_gc_N1
gc_x_center = avg_gc_H1
gc_y_center = avg_gc_N1

gc_x = np.linspace(0, 15, 1000)
gc_y = np.linspace(135, 155, 1000)

gc_X, gc_Y = np.meshgrid(gc_x, gc_y)

gc_eqn = (((gc_X - gc_x_center) ** 2) / (gc_horizontal_radius ** 2)) + (((gc_Y - gc_y_center) ** 2) / (
        gc_vertical_radius ** 2))
Z = 1

plt.contour(gc_X, gc_Y, gc_eqn, [Z], colors="k", linewidths=1, linestyles="dashed")

# UG
ug_horizontal_radius = std_ug_H3
ug_vertical_radius = std_ug_N3
ug_x_center = avg_ug_H3
ug_y_center = avg_ug_N3

ug_x = np.linspace(0, 15, 1000)
ug_y = np.linspace(150, 165, 1000)

ug_X, ug_Y = np.meshgrid(ug_x, ug_y)

ug_eqn = (((ug_X - ug_x_center) ** 2) / (ug_horizontal_radius ** 2)) + (
            ((ug_Y - ug_y_center) ** 2) / (ug_vertical_radius ** 2))
Z = 1

plt.contour(ug_X, ug_Y, ug_eqn, [Z], colors="k", linewidths=1, linestyles="dashed")

# UA
ua_horizontal_radius = std_ua_H3
ua_vertical_radius = std_ua_N3
ua_x_center = avg_ua_H3
ua_y_center = avg_ua_N3

ua_x = np.linspace(0, 15, 1000)
ua_y = np.linspace(150, 165, 1000)

ua_X, ua_Y = np.meshgrid(ua_x, ua_y)

ua_eqn = (((ua_X - ua_x_center) ** 2) / (ua_horizontal_radius ** 2)) + (
            ((ua_Y - ua_y_center) ** 2) / (ua_vertical_radius ** 2))
Z = 1

plt.contour(ua_X, ua_Y, ua_eqn, [Z], colors="k", linewidths=1, linestyles="dashed")

plt.title("Imino chemical shifts", fontsize=16)
plt.xlabel("1H (ppm)", fontsize=16)
plt.ylabel("15N (ppm)", fontsize=16)
plt.legend()

plt.xlim(15.0, 9.0)
plt.ylim(170.0, 135.0)
plt.show()

if SAVE_FIGURE:
    if FIGURE_FORMAT == "eps":
        fig.savefig("bmrb_imino_all.eps", format='eps')
    elif FIGURE_FORMAT == "pdf":
        fig.savefig("bmrb_imino_all.pdf", format='pdf')
