# -*- coding: utf-8 -*-
"""
Created on Sun Apr 22 19:23:41 2018

@author: Annalise
"""

import csv

text_files = ['CFD_colorado_slave2.csv']
poss_directions = [270, 280, 290, 300, 310, 320, 330, 340, 350]
# poss_directions = [270, 290, 310, 330, 350]
poss_ws = [5., 10., 15., 20., 25.]
# for 10 degree increments
# prob_ct = [[0, 0, 0] for ii in range(poss_directions)]
# for 20 degree increments
prob_ct = [[0, 0] for ii in range(9)]
for thisfile in text_files:
    print(thisfile)
    with open(thisfile, 'r', newline='') as outfile:
        info = csv.reader(outfile, delimiter=',', quotechar='|')

        for ii, row in enumerate(info):  # subset for testing!
            direction_coord = poss_directions.index(int(float(row[5])))
#            except ValueError:
#                # for 20 degree increments
#                direction_coord = poss_directions.index(int(float(row[5]) - 10))
            windsp_coord = 0
            while poss_ws[windsp_coord] < float(row[6]):
                windsp_coord += 1
            if windsp_coord < 3 and windsp_coord > 0:
                prob_ct[direction_coord][windsp_coord - 1] += 1. / 2870.
print(prob_ct)
print(sum([sum(i) for i in prob_ct]))
