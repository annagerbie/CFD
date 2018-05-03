#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 11:10:13 2018

@author: millanna
"""
import csv
# import os
import numpy as np
import matplotlib.pyplot as plt

# path = 'CFD_colorado_results'
# text_files = (['CFD_colorado_results/'
#               + str(f) for f in os.listdir(path) if f.endswith('.csv')])
# print(text_files)
text_files = ['CFD_colorado_slave2.csv']
aif = 0.314
bins = 1  # Unidirectional wind
numturbs = 37
poss_directions = [270, 280, 290, 300, 310, 320, 330, 340, 350]
poss_ws = [5., 10., 15., 20., 25.]
poss_ust = [0, 1, 2, 3, 4, 5]
error = [[0.] * numturbs for i in range(len(poss_directions))]
error_ct = [[0.] * numturbs for i in range(len(poss_directions))]
under = [[0.] * numturbs for i in range(len(poss_directions))]
direc = [[0.] * numturbs for i in range(len(poss_directions))]
howmuch = [[0.] * numturbs for i in range(len(poss_directions))]
ave_dev = [[0.] * numturbs for i in range(len(poss_directions))]

error2 = [[0.] * numturbs for i in range(len(poss_ws))]
error2_ct = [[0.] * numturbs for i in range(len(poss_ws))]
under2 = [[0.] * numturbs for i in range(len(poss_ws))]
direc2 = [[0.] * numturbs for i in range(len(poss_ws))]
howmuch2 = [[0.] * numturbs for i in range(len(poss_ws))]
ave_dev2 = [[0.] * numturbs for i in range(len(poss_ws))]

error3 = [[0.] * numturbs for i in range(len(poss_ust))]
error3_ct = [[0.] * numturbs for i in range(len(poss_ust))]
under3 = [[0.] * numturbs for i in range(len(poss_ust))]
direc3 = [[0.] * numturbs for i in range(len(poss_ust))]
howmuch3 = [[0.] * numturbs for i in range(len(poss_ust))]
ave_dev3 = [[0.] * numturbs for i in range(len(poss_ust))]

RSME = [[0.] * numturbs for i in range(len(poss_directions))]
# by wind direction
RSME2 = [[0.] * numturbs for i in range(len(poss_ws))]
# by wind speed
RSME3 = [[0.] * numturbs for i in range(len(poss_ust))]
# by number of upstream turbines
num_usturbines = [0] * 10
num_usturbines += [2, 3, 2, 2, 2, 1] + [0] * 6
num_usturbines += [3, 1, 2, 3, 3, 2, 3, 1] + [0] * 6 + [1]
# wind_conditions = []
all_error = []
for thisfile in text_files:
    print(thisfile)
    with open(thisfile, 'r', newline='') as outfile:
        info = csv.reader(outfile, delimiter=',', quotechar='|')

        for ii, row in enumerate(info):  # subset for testing!
                turbines_realsp = [float(temp) for temp in row[7:(7 + 37)]]
                # windspeeds measured from nacelles
                # print(row[-37:])
                cum_ws = [float(each) for each in row[-37:]]
                # print(row[5])
                # print(type(row[5]))
                # print(float(row[5]))
                # wind_conditions.append((float(row[5]), float(row[6])))
                direction_coord = poss_directions.index(int(float(row[5])))
                for i in range(numturbs):
                    this_error = turbines_realsp[i] - cum_ws[i] / (1. - aif)
                    ###this_error = turbines_realsp[i] - cum_ws[i]
                    all_error.append(this_error)
                    # print('real speed: ',turbines_realsp[i])
                    # print('sim speed: ',cum_ws[i])
                    # print('sim speed/(1-a): ',cum_ws[i] / (1 - aif))
                    coordi = pow(turbines_realsp[i] - cum_ws[i] / (1. - aif),
                                 2)
                    ###coordi = pow(turbines_realsp[i] - cum_ws[i],2)
                    if turbines_realsp[i] - cum_ws[i] / (1. - aif) > 0.:
                    ###if turbines_realsp[i] - cum_ws[i] > 0.:
                        # underpredicting
                        under[direction_coord][i] += 1.
                    by_much = turbines_realsp[i] - cum_ws[i] / (1. - aif)
                    ###by_much = turbines_realsp[i] - cum_ws[i]
                    howmuch[direction_coord][i] += by_much
                    error[direction_coord][i] += coordi
                    error_ct[direction_coord][i] += 1.

                windsp_coord = 0

                while poss_ws[windsp_coord] < float(row[6]):
                    windsp_coord += 1
                for i in range(numturbs):
                    coordi = pow(turbines_realsp[i] - cum_ws[i] / (1. - aif),
                                 2)
                    ###coordi = pow(turbines_realsp[i] - cum_ws[i], 2)
                    if turbines_realsp[i] - cum_ws[i] / (1. - aif) > 0.:
                    ###if turbines_realsp[i] - cum_ws[i] > 0.:
                        # underpredicting
                        under2[windsp_coord][i] += 1.
                    by_much = turbines_realsp[i] - cum_ws[i] / (1. - aif)
                    ###by_much = turbines_realsp[i] - cum_ws[i]
                    howmuch2[windsp_coord][i] += by_much
                    error2[windsp_coord][i] += coordi
                    error2_ct[windsp_coord][i] += 1.

                for i in range(numturbs):
                    ust_coord = poss_ust.index(num_usturbines[i])
                    coordi = pow(turbines_realsp[i] - cum_ws[i] / (1. - aif),
                                 2)
                    ###coordi = pow(turbines_realsp[i] - cum_ws[i], 2)
                    if turbines_realsp[i] - cum_ws[i] / (1. - aif) > 0.:
                    ###if turbines_realsp[i] - cum_ws[i] > 0.:
                        # underpredicting
                        under3[ust_coord][i] += 1.
                    by_much = turbines_realsp[i] - cum_ws[i] / (1. - aif)
                    ###by_much = turbines_realsp[i] - cum_ws[i]
                    howmuch3[ust_coord][i] += by_much
                    error3[ust_coord][i] += coordi
                    error3_ct[ust_coord][i] += 1.
    print('last ii: ', ii)

output = 'on'
if output == 'on':
    tot_dir = [0.] * len(error)
    tot_ws = [0.] * len(error2)
    tot_ust = [0.] * len(error3)
    for i, j in enumerate(error):
        try:
            tot_dir[i] = np.sqrt(sum(j) / sum(error_ct[i]))
            direc[i] = sum(under[i]) / sum(error_ct[i])
            # percent under predicting
            ave_dev[i] = sum(howmuch[i]) / sum(error_ct[i])
        except ZeroDivisionError:
            tot_dir[i] = 'N/A'
            direc[i] = 'N/A'
            # percent under predicting
            ave_dev[i] = 'N/A'
        for a, b in enumerate(j):
            try:
                RSME[i][a] = np.sqrt(b / error_ct[i][a])
            except ZeroDivisionError:
                RSME[i][a] = 'N/A'

    for i, j in enumerate(error2):
        try:
            tot_ws[i] = np.sqrt(sum(j) / sum(error2_ct[i]))
            direc2[i] = sum(under2[i]) / sum(error2_ct[i])
            # percent under predicting
            ave_dev2[i] = sum(howmuch2[i]) / sum(error2_ct[i])
        except ZeroDivisionError:
            tot_ws[i] = 'N/A'
            direc2[i] = 'N/A'
            # percent under predicting
            ave_dev2[i] = 'N/A'
        for a, b in enumerate(j):
            try:
                RSME2[i][a] = np.sqrt(b / error2_ct[i][a])
            except ZeroDivisionError:
                RSME2[i][a] = 'N/A'

    for i, j in enumerate(error3):
        try:
            tot_ust[i] = np.sqrt(sum(j) / sum(error3_ct[i]))
            direc3[i] = sum(under3[i]) / sum(error3_ct[i])
            # percent under predicting
            ave_dev3[i] = sum(howmuch3[i]) / sum(error3_ct[i])
        except ZeroDivisionError:
            tot_ust[i] = 'N/A'
            direc3[i] = 'N/A'
            # percent under predicting
            ave_dev3[i] = 'N/A'
        for a, b in enumerate(j):
            try:
                RSME3[i][a] = np.sqrt(b / error3_ct[i][a])
            except ZeroDivisionError:
                RSME3[i][a] = 'N/A'

    with open('RMSE_CFD_colorado_byturb.csv', 'w+',
              newline='') as outfile:
        print('opening byturb_rotor')
        write = csv.writer(outfile)
        blank = ['']
        blank.extend([i for i in range(37)])
        write.writerow(blank)
        for k, i in enumerate(RSME):
            start_dir = [poss_directions[k]]
            start_dir.extend(i)
            write.writerow(start_dir)
        write.writerow([])
        for k, i in enumerate(RSME2):
            start_ws = [poss_ws[k]]
            start_ws.extend(i)
            write.writerow(start_ws)
        write.writerow([])
        for k, i in enumerate(RSME3):
            start_ust = [poss_ust[k]]
            start_ust.extend(i)
            write.writerow(start_ust)
    with open('RMSE_CFD_colorado.csv', 'w+', newline='') as outfile:
        print('opening _rotor')
        write = csv.writer(outfile)
        write.writerow(['wind directions']+poss_directions)
        write.writerow(['RMSE']+tot_dir)
        write.writerow(['percent underpredicted']+direc)
        write.writerow(['average error']+ave_dev)
        write.writerow([])
        write.writerow(['wind speeds']+poss_ws)
        write.writerow(['RMSE']+tot_ws)
        write.writerow(['percent underpredicted']+direc2)
        write.writerow(['average error']+ave_dev2)
        write.writerow([])
        write.writerow(['number of upstream turbines']+poss_ust)
        write.writerow(['RMSE']+tot_ust)
        write.writerow(['percent underpredicted']+direc3)
        write.writerow(['average error']+ave_dev3)
        # print('The final layout with ' + str(initial_num_
        #       + ' turbines has a score of: ' + str(score))
# make a histogram of all errors
# plt.figure()
f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
ax1.hist(all_error)
ax1.set_title('CFD model error frequency plot')
ax1.set_ylabel('Frequency')
plt.xlabel('error (m/s)')
plt.xlim(-20, 15)
# all_error.sort()
# for i in range(5):
#     index = i * len(all_error) / 4.
#     quartiles.append(all_error[index])
outbox = ax2.boxplot(all_error, vert=False)
# outliers = outbox['caps']
# outliers = [float(i) for i in outliers]
# print(len(outliers))
plt.savefig('cfd_error_hist.png')
indices = []
# outlier_directions = []
for ct, i in enumerate(all_error):
    if i < -5 or i > 7:
        indices.append(int(ct / 37))
print(len(indices))
print(len(set(indices)))
outlier_conditions = []
outlier_ct = []
with open('CFD_colorado_slave2.csv', newline='') as csvfile:
    info = csv.reader(csvfile, delimiter=',', quotechar='|')
    for i, row in enumerate(info):
        if i in indices:
            condition = (float(row[5]), float(row[6]))
            if condition not in outlier_conditions:
                outlier_conditions.append(condition)
                outlier_ct.append(1)
            else:
                index = outlier_conditions.index(condition)
                outlier_ct[index] += 1
# statistical significance
n = 133644.
mean_CFD = np.mean(all_error)
std_CFD = np.std(all_error)
mean_jensen = 1.3474757059066387
mean_jensen_nwp = 1.3274146920382719
std_jensen = 2.0496650238745571
std_jensen_nwp = 2.0536242366793984
 
var_CFD = pow(std_CFD,2) / (n - 1.)
var_jensen = pow(std_jensen, 2) / (n - 1.)
var_jensennwp = pow(std_jensen_nwp, 2) / (n - 1.)
p = (mean_CFD - mean_jensen) / np.sqrt(var_CFD + var_jensen)
p2 = (mean_CFD - mean_jensen_nwp) / np.sqrt(var_CFD + var_jensennwp)
print('done!')
