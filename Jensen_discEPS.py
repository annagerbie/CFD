# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 11:21:39 2016

@author: Annalise
"""

#import math
import random
#import os
import matplotlib.pyplot as plt
import os
global nemo
import numpy as np
from time import time


init_step = 8. #number of starting mesh intervals
mesh_width = 200. #meters - width of mesh
#area = 2000.0              #Length of one side of solution area, in m

site_x = 1999.             #Length of one side of solution area, in m  
site_y = 1999.              #Length of one side of solution area, in m  
XLength = site_x * 1.2             #Length of one side of solution area, in m
YLength = site_y * 1.2               #Length of one side of solution area, in m
U0 = [7.5, 12.5, 17.5]      #mean wind speed, in m/s
#z = 60                     #hub height, in m
#Rr = 20                    #Rotor Radius, in m
Ct = 0.88                   #Thrust Coefficient
#aif = 0.1909830056250526    #Axial induction factor based on Cp = 0.5
#aif = 0.25
aif = 1./3. #for testing
# wind_cases = [[1.0]]
#wind_cases = [(0.004, 0.009, 0.011), (0.004, 0.009, 0.011), (0.004, 0.009, 0.011), (0.004, 0.009, 0.011), (0.004, 0.009, 0.011), (0.004, 0.009, 0.011), (0.004, 0.009, 0.011), (0.004, 0.009, 0.011), (0.004, 0.009, 0.011), (0.004, 0.009, 0.011),
#              (0.004, 0.009, 0.011), (0.004, 0.009, 0.011), (0.004, 0.009, 0.011), (0.004, 0.009, 0.011), (0.004, 0.009, 0.011), (0.004, 0.009, 0.011), (0.004, 0.009, 0.011), (0.004, 0.009, 0.011), (0.004, 0.009, 0.011), (0.004, 0.009, 0.011),
#              (0.004, 0.009, 0.011), (0.004, 0.009, 0.011), (0.004, 0.009, 0.011), (0.004, 0.009, 0.011), (0.004, 0.009, 0.011), (0.004, 0.009, 0.011), (0.004, 0.009, 0.011), (0.004, 0.011, 0.012), (0.004, 0.013, 0.015), (0.004, 0.016, 0.017), 
#              (0.004, 0.015, 0.031), (0.004, 0.019, 0.037), (0.004, 0.015, 0.031), (0.004, 0.016, 0.017), (0.004, 0.013, 0.015), (0.004, 0.011, 0.012)]
        
minstep = 3.0              #minimum step size to break pattern search
#CHECK = 0
#CHECK2 = 0
#xzeroval = 0.0
#xareaval = 0.0
#yzeroval = 0.0
#yareaval = 0.0
prelim = 0.0
final = 0.0
transflag = 0.0
ro = 1.225                  #air density in kg/m3
        
num_pops = 5               #Number of Turbines the program pops
max_pop_tries = 1000        #Maximum number of placement attempts before exiting
tolerance = 0.000000001
random_vec = []             #For selecting turbine order in Pattern Search
hstep = 45.0                #Step size for hub height search, in m
rstep = 25.0                #Step size for rotor radius search, in m
hstepmin = 1.0              #Minimum step size before hub height search exits, in m
rstepmin = 1.0              #Minimum step size before rotor radius search exits, in m
rradmin = 19.0              #Minimum rotor radius, in m
rradmax = 67.0              #Maimum rotor radius, in m
hubmin = 38.0               #Minimum hub height, in m
hubmax = 135.0              #Maximum hub height, in m
Uref = 11.5                 #Reference wind speed for wind shear calculation, in m
Zref = 80.0                 #Reference height for wind shear calculation, in m
depth = 200.0               #water depth, 200m. User-defined in future
WCOE = 0.1                  #wholdsale cost of electricity, in $/kWh
yrs = 20.0

output = 'on'
nwp = False
#ma = 0.33
Cp = 16./27. #for testing
#Cp = 0.34
#Cf = 0.4
Cf = 1.
Cutin = 3.0         #cut-in wind speed

shore = 'on'
#onshore/offshore specific
if shore == 'on':
    z0 = 0.05                 #Surface roughness, in m ###changed to match NREL
    
elif shore == 'off':
    z0 = 0.0005                 #Surface roughness, in m
degree10 = False
speeds3 = True
directions  = []
#degrees = [0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 210., 220., 230., 240., 250., 260., 270., 280., 290., 300., 310., 320., 330., 340., 350.]
if degree10:
    degrees = [270, 280, 290, 300, 310, 320, 330, 340, 350]
    for l in range(0, len(degrees)):
        radians = (degrees[l] / 180.) * np.pi
        directions.append(radians)
    #prob = [1./len(directions)] * len(directions)   #probability of wind from each direction
    #probu0 = [1./len(U0)] * len(U0) #probability of each ambiant wind speed
    probwui = [[0.019063545150501682, 0.011036789297658863, 0.001003344481605351],
               [0.02073578595317727, 0.015719063545150506, 0.00033444816053511704],
               [0.05785953177257504, 0.017056856187290977, 0],
               [0.06923076923076892, 0.019732441471571917, 0.001003344481605351],
               [0.05551839464882924, 0.025752508361204032, 0.0026755852842809363],
               [0.12374581939799251, 0.0521739130434781, 0.00869565217391304],
               [0.12809364548494911, 0.11070234113712306, 0.011036789297658863],
               [0.07525083612040097, 0.06688963210702312, 0.012374581939799333],
               [0.06755852842809335, 0.023745819397993327, 0.0030100334448160534]]
elif speeds3:
    degrees = [270, 290, 310, 330, 350]
    for l in range(0, len(degrees)):
        radians = (degrees[l] / 180.) * np.pi
        directions.append(radians)
    #prob = [1./len(directions)] * len(directions)   #probability of wind from each direction
    #probu0 = [1./len(U0)] * len(U0) #probability of each ambiant wind speed
    probwui = [[0.03979933110367888, 0.026755852842809385, 0.0013377926421404682],
               [0.12709030100334373, 0.03678929765886285, 0.001003344481605351],
               [0.1792642140468237, 0.07792642140468188, 0.01137123745819398],
               [0.2033444816053529, 0.17759197324414805, 0.02341137123745821],
               [0.06755852842809335, 0.023745819397993327, 0.0030100334448160534]]
else:
    U0 = [7.5, 12.5]
    degrees = [270, 280, 290, 300, 310, 320, 330, 340, 350]
    for l in range(0, len(degrees)):
        radians = (degrees[l] / 180.) * np.pi
        directions.append(radians)
    probwui = [[0.01986062717770035, 0.011498257839721245],
               [0.021602787456446, 0.016376306620209048],
               [0.06027874564459945, 0.01777003484320557],
               [0.07212543554006969, 0.02055749128919861],
               [0.05783972125435554, 0.026829268292682954],
               [0.12891986062717695, 0.05435540069686424],
               [0.13344947735191573, 0.11533101045296097],
               [0.0783972125435539, 0.06968641114982582],
               [0.07038327526132407, 0.024738675958188173]]
#condition = raw_input("Input Condition NEUTRAL, STABLE, or UNSTABLE: ") #Input wind condition
condition = 'NEUTRAL'

if condition == "NEUTRAL":
    APow = 0.08835          #Neutral Conditions: WDC(h) = APow*h^BPow
    BPow = -0.1521 
    if shore == 'off':
        alphah = 0.11           #Power Law Exponent - Neutral (averaged over seasons)
    elif shore == 'on':
        alphah = 0.15567      
    
elif condition == "STABLE":
    APow = 0.07535
    BPow = -0.1496
    alphah = 0.14
    
elif condition == "UNSTABLE":
    APow = 0.09759
    BPow = -0.1352
    alphah = 0.08
 
layout_type = 'random'

if layout_type == 'test':
    YLocation = [[0.], [200.], [200.]]
    ZLocation = [[0.], [0.], [0.]]
    initial_num = 3
    spacing = 300.
    XLocation = [[0.], [-spacing / 2.], [spacing / 2.]]
    for i in range(len(XLocation)):
        xfinish1 = [XLocation[i][0]]
        yfinish1 = [YLocation[i][0]]
        for k in range(1, len(directions)):
            theta = directions[k]
            newy = (xfinish1[0] * np.sin(theta)) + (yfinish1[0] * np.cos(theta))            #find preliminary new x-location given step size  
            yfinish1.append(newy)
            newx = (xfinish1[0] * np.cos(theta)) - (yfinish1[0] * np.sin(theta))    
            xfinish1.append(newx)
        YLocation[i] = yfinish1
        XLocation[i] = xfinish1
    
if layout_type == 'grid':
    initExtent = 0.95
    rows = 4
    cols = 4
    xpos = np.linspace(-initExtent*(site_x - 40.),initExtent*(site_x - 40.),cols)
    ypos = np.linspace(-initExtent*(site_y - 40.),initExtent*(site_y - 40.),rows)
    #print('xlocations')
    #print(xpos)
    YLocation = []
    XLocation = []
    ZLocation = []
    for i in range(rows):
        for j in range(cols):
            xfinish1 = [xpos[j]]
            yfinish1 = [ypos[i]]
            for k in range(1, len(directions)):
                theta = directions[k]
                newy = (xpos[j] * np.sin(theta)) + (ypos[i] * np.cos(theta))            #find preliminary new x-location given step size  
                yfinish1.append(newy)
                newx = (xpos[j] * np.cos(theta)) - (ypos[i] * np.sin(theta))    
                xfinish1.append(newx)
            YLocation.append(yfinish1)
            XLocation.append(xfinish1)
            # # some starting noise sometimes helps
            # mx.append(Constant(xpos[j]+5.*np.random.randn()))
            # my.append(Constant(ypos[i]+5.*np.random.randn()))
            ZLocation.append(0.)
    #title = 'Grid Layout'
    initial_num = 16
                
elif layout_type == 'offset':
    rows = 4
    cols = 4
    xpos = np.linspace(-initExtent*(site_x - 40.),initExtent*(site_x - 40.),cols)
    ypos = np.linspace(-initExtent*(site_y - 40.),initExtent*(site_y - 40.),rows)
    offset = ypos[1] - ypos[0]
    xpos2 = []
    for i in ypos:
        xpos2.append(i + offset/2)
    #print('xlocations')
    #print(xpos)
    YLocation = []
    XLocation = []
    for i in range(rows):
        for j in range(cols):
            if i % 2 != 0:
                xfinish1 = [xpos[j]]
            else:
                xfinish1 = [xpos2[j]]
            yfinish1 = [ypos[i]]
            for k in range(1, len(directions)):
                theta = directions[k]
                newy = (xfinish1[0] * np.sin(theta)) + (yfinish1[0] * np.cos(theta))            #find preliminary new x-location given step size  
                yfinish1.append(newy)
                newx = (xfinish1[0] * np.cos(theta)) - (yfinish1[0] * np.sin(theta))    
                xfinish1.append(newx)
            YLocation.append(yfinish1)
            XLocation.append(xfinish1)
            # # some starting noise sometimes helps
            # mx.append(Constant(xpos[j]+5.*np.random.randn()))
            # my.append(Constant(ypos[i]+5.*np.random.randn()))
            ZLocation.append(0.)
    #title = 'Offset Grid Layout'
    initial_num = 16

elif layout_type == 'colorado':
    colorado_x = [-6680.980411, -6549.425159, -6380.075027, -6197.64205, -6015.935896, -5834.229741, -3668.29226, -3799.847528, -3435.70836, -3230.016971, -6651.907427, -6455.664786, -6265.236739, -6075.535515, -5810.971353, -5660.518655, -5505.705008, -5284.750318, -5082.693066, -4877.728514, -4743.992777, -4621.159409, -4735.270881, -4497.599216, -3722.077287, -3590.522019, -3462.600873, -3332.499253, -3274.353278, -3082.471558, -2533.718907, -2334.568937, -2151.409109, -1912.283778, -1732.758073, -1599.022325, -1551.778718]
    colorado_y = [-16892.41457, -16682.25616, -16486.55309, -16286.40222, -16109.60229, -15913.89922, -15251.17746, -15518.04528, -15084.38507, -14958.7348, -14737.4569, -14570.66451, -14410.54381, -14241.52752, -13962.42826, -13814.539, -13652.19441, -13447.59575, -13315.27378, -13175.16818, -12957.22612, -12707.03754, -11972.03907, -11891.97872, -11993.16611, -11794.12719, -11599.53607, -11399.3852, -11183.66704, -10937.92625, -10207.37558, -10036.1354, -9870.454957, -9691.431125, -9536.870177, -9323.375918, -9012.030123]
    #colorado_x = [6680.980411, 6549.425159, 6380.075027, 6197.64205, 6015.935896, 5834.229741, 3668.29226, 3799.847528, 3435.70836, 3230.016971, 6651.907427, 6455.664786, 6265.236739, 6075.535515, 5810.971353, 5660.518655, 5505.705008, 5284.750318, 5082.693066, 4877.728514, 4743.992777, 4621.159409, 4735.270881, 4497.599216, 3722.077287, 3590.522019, 3462.600873, 3332.499253, 3274.353278, 3082.471558, 2533.718907, 2334.568937, 2151.409109, 1912.283778, 1732.758073, 1599.022325, 1551.778718]
    #colorado_y = [16892.41457, 16682.25616, 16486.55309, 16286.40222, 16109.60229, 15913.89922, 15251.17746, 15518.04528, 15084.38507, 14958.7348, 14737.4569, 14570.66451, 14410.54381, 14241.52752, 13962.42826, 13814.539, 13652.19441, 13447.59575, 13315.27378, 13175.16818, 12957.22612, 12707.03754, 11972.03907, 11891.97872, 11993.16611, 11794.12719, 11599.53607, 11399.3852, 11183.66704, 10937.92625, 10207.37558, 10036.1354, 9870.454957, 9691.431125, 9536.870177, 9323.375918, 9012.030123]
    initial_num = 37  
    ZLocation = np.array([0.] * len(colorado_x))
    YLocation = []
    XLocation = []
    for j in range(len(colorado_x)):
        xfinish1 = [colorado_x[j]]
        yfinish1 = [colorado_y[j]]
        for k in range(1, len(directions)):
            theta = directions[k]
            newy = (xfinish1[0] * np.sin(theta)) + (yfinish1[0] * np.cos(theta))            #find preliminary new x-location given step size  
            yfinish1.append(newy)
            newx = (xfinish1[0] * np.cos(theta)) - (yfinish1[0] * np.sin(theta))    
            xfinish1.append(newx)
        YLocation.append(yfinish1)
        XLocation.append(xfinish1)
    
elif layout_type == 'linear':
    '''My adds'''
    ypos = np.linspace(-initExtent*(site_x - 40.),initExtent*(site_x - 40.),10)
    xpos = np.array([0.] * 10)
    ZLocation = np.array([0.] * 10)
    YLocation = []
    XLocation = []
    for j in range(10):
        xfinish1 = [xpos[j]]
        yfinish1 = [ypos[j]]
        for k in range(1, len(directions)):
            theta = directions[k]
            newy = (xfinish1[0] * np.sin(theta)) + (yfinish1[0] * np.cos(theta))            #find preliminary new x-location given step size  
            yfinish1.append(newy)
            newx = (xfinish1[0] * np.cos(theta)) - (yfinish1[0] * np.sin(theta))    
            xfinish1.append(newx)
        YLocation.append(yfinish1)
        XLocation.append(xfinish1)
    #title = 'Inline Layout'
    initial_num = 10
    
elif layout_type == 'optimized':
    xpos = [1304.2130491853543, 303.16195377176825, 1882.7665062612382, 734.1365122392898, 36.582152829805295, 1078.6666766839614, 1779.34556196998, 417.79552464946823, 1986.6903503041713, 133.48666173166612, 1430.1977444645, 954.8374274350235, 839.2773010244406, 1186.314920647951, 518.9482768461385, 1681.3296876680367, 621.3034610992612, 1567.283691418906]
    ypos = [204.01944774702838, 270.02673301978496, 276.97611783621744, 305.2706459401636, 1065.3030786448571, 104.67022680498235, 610.6008734856377, 450.6469581295796, 32.94310546319856, 1290.6284216353858, 413.28638998844144, 611.4810784612052, 117.69510957163021, 489.9943326680768, 225.4370943493467, 858.4028566999443, 505.9650483803546, 1212.01937488221]
    #title = 'Optimized Layout'
    ZLocation = np.array([0.] * len(xpos))
    YLocation = []
    XLocation = []
    for j in range(len(xpos)):
        xfinish1 = [xpos[j]]
        yfinish1 = [ypos[j]]
        for k in range(1, len(directions)):
            theta = directions[k]
            newy = (xfinish1[0] * np.sin(theta)) + (yfinish1[0] * np.cos(theta))            #find preliminary new x-location given step size  
            yfinish1.append(newy)
            newx = (xfinish1[0] * np.cos(theta)) - (yfinish1[0] * np.sin(theta))    
            xfinish1.append(newx)
        YLocation.append(yfinish1)
        XLocation.append(xfinish1)
    initial_num = len(xpos)
        
# Store Information on Turbines
class Turbine(object):
    hood = []
    Prated = 0.0
    Capital_Cost = 0.0
    O_M = 0.0
    Substation_Cost = 0.0
    Leasing_Cost = 0.0
    Mooring_Cost = 0.0
    Installation_Cost = 977620.0
    Cost = 0.0
    Stopped = 0
    def __init__(self, XLocation, YLocation, ZLocation, HubHeight, RotorRad, alpha, usturbinesrec, usturbines, dsturbinesrec, dsturbines, wakewidth, distance, percent, ui, windspeeds, xcoords, zcoords, Power, Area):
        self.XLocation = XLocation
        self.YLocation = YLocation
        self.ZLocation = ZLocation
        self.HubHeight = HubHeight
        self.ZHub = ZLocation + HubHeight
        self.RotorRad = RotorRad
        self.alpha = alpha
        self.usturbinesrec = usturbinesrec
        self.usturbines = usturbines
        self.dsturbinesrec = dsturbinesrec
        self.dsturbines = dsturbines
        self.wakewidth = wakewidth
        self.distance = distance
        self.percent = percent
        self.ui = ui
        self.windspeeds = windspeeds
        self.xcoords = xcoords
        self.zcoords = zcoords
        self.Power = Power
        self.Area = Area

        
#######################################################################################################



###################################################################################################################
def Translation_X(step_size, i, directions):
    transflag = 0
    xstart = turbines[i].XLocation[0]
    ystart = turbines[i].YLocation[0]
    #print(xstart, step_size)
    xfinish = xstart + step_size            #find preliminary new x-location given step size
    if xfinish >= 0 and xfinish <= site_x:    #if this new x-location is not out of bounds, translate it
        xfinish1 = [xfinish]
        yfinish1 = [ystart]
        for j in range(1, len(directions)):
            theta = directions[j]
            newy = (xfinish * np.sin(theta)) + (ystart * np.cos(theta))
            yfinish1.append(newy)
            newx = (xfinish * np.cos(theta)) - (ystart * np.sin(theta)) #find preliminary new x-location given step size  
            xfinish1.append(newx)
        turbines[i].XLocation = xfinish1
        turbines[i].YLocation = yfinish1
        return transflag
        
    else:
        transflag = 1
        return transflag


################################################################################################################        
def Translation_Y(step_size, i, directions):
    transflag = 0
    ystart = turbines[i].YLocation[0]
    xstart = turbines[i].XLocation[0]
    yfinish = ystart + step_size            #find preliminary new x-location given step size
    if yfinish >= 0 and yfinish <= site_y:    #if this new x-location is not out of bounds, translate it  
        xfinish1 = [xstart]
        yfinish1 = [yfinish]
        for j in range(1, len(directions)):
            theta = directions[j]
            newy = (xstart * np.sin(theta)) + (yfinish * np.cos(theta))            #find preliminary new x-location given step size  
            yfinish1.append(newy)
            newx = (xstart * np.cos(theta)) - (yfinish * np.sin(theta))    
            xfinish1.append(newx)
        turbines[i].YLocation = yfinish1
        turbines[i].XLocation = xfinish1
        return transflag
    else:
        transflag = 1
        return transflag  
        
#################################################################################################################    
def Initial_Layout():
    
    for n in range(0, initial_num):
        reset = 0        
        checkx = 0
        while checkx == 0 and reset < 50000:
            xmove = 0
            ymove = 0
            CHECK2 = 0
            x_opts = int(site_x / mesh_width) + 1
            y_opts = int(site_y / mesh_width) + 1
            xmove = int(random.uniform(0, x_opts)) * mesh_width
            ymove = int(random.uniform(0, y_opts)) * mesh_width
            #print('xmove = ', xmove)
            Translation_X(xmove, n, directions)
            #print('ymove = ', ymove)
            Translation_Y(ymove, n, directions)
            #print(turbines[n].YLocation)
            
            CHECK2 = Check_Interference(n)
            
            if CHECK2 != 1:
                checkx = 1                  #If there is no interference and the turbine can be placed, then choose new corresponding z-coordinate
                
            else:
                Translation_X(-xmove, n, directions)  #If there is interference, move turbine back to origin
                Translation_Y(-ymove, n, directions)
                reset += 1
                
        if reset == 50000:
            for l in range(0, initial_num):
                turbines[l].XLocation[0] = 0.0
                turbines[l].YLocation[0] = 0.0
                for j in range(1, len(directions)):
                    turbines[l].XLocation[j] = 0.
                    turbines[l].YLocation[j] = 0.
            return reset
    print('XLocation')
    print([i.XLocation[0] for i in turbines])
    print('YLocation')
    print([i.YLocation[0] for i in turbines])
    return reset

    
###################################################################################################################    
def Check_Interference(n):
    CHECK2 = 0
    checkx = 0
    x = turbines[n].XLocation[0]
    y = turbines[n].YLocation[0]
    for k in range(0, initial_num):
        if k != n:
            xnew = turbines[k].XLocation[0]
            ynew = turbines[k].YLocation[0]
            checkx = x - xnew
            checky = y - ynew
            checkrad = np.sqrt(checkx ** 2.0 + checky ** 2.0)
            if checkrad < 200:
                if checkrad > 0.:
                    print('failed interference: ',checkrad)
                CHECK2 = 1
    return CHECK2

#Code Check for Initial Layout
#XLocation, YLocation = Initial_Layout()
#print(XLocation)
#for n in range(0, initial_num):
    #print(turbines[n].XLocation) #Store new coordinates in turbine class
#checkcheck = []
#for k in range(0, initial_num):
#    for n in range(0, initial_num):
#        if n > k:
#            checkrad = []
#            checkx = XLocation[k] - XLocation[n]
#            checky = YLocation[k] - YLocation[n]
#            checkrad.append(np.sqrt(checkx ** 2 + checky ** 2))
#    checkcheck.append(min(checkrad))
#print (checkcheck)  a
###############################################################################
def check_layout(mx,my):
    coords = [(i,j) for i,j in zip(mx,my)] #zip into tuple
    has_been = True
    check = 0
    num = 0
    with open('layouts_E.txt', 'r') as layouts_file:
        #print(sum(1 for i in layouts_file))
        line = layouts_file.readline()
        while line:
            num += 1
            line = line.replace('[','').replace(']','')
            #print(line)
            #line = line.split(',')
            #print(line)
            old_coords = [float(x.strip()) for x in line.split(',')]
            old_coords = [(old_coords[i], old_coords[i + initial_num]) for i in range(initial_num)]
            for i in coords:
                #print(i)
                space_tot = [np.sqrt((j[0] - i[0]) ** 2 + (j[1] - i[1]) ** 2) for j in old_coords] #list of distance between each point and existing point
                if all(i >= 0.005 for i in space_tot): #new point
                    check += 1
                    #print('pt not in this layout')
                    break
            line = layouts_file.readline()
        layouts_file.close()
    if abs(check - num) < 0.00005:
        has_been = False
    if has_been == False:
        with open('layouts_E.txt', 'a+') as layouts_file:
            new_layout = mx + my
            #print('new layout: ',new_layout)
            layouts_file.write(str(new_layout)+'\n')
    else:
        pass
        #print('layout previously checked')
        
    return has_been
################################################################################################################
def Compute_Wake(initial_num, z0, U0, Zref, alphah, ro, aif):
    #lastj = 12
    for j in range(0, initial_num): #downstream turbine
        upstrm = []   
        distanceus = []
        wakewidthds = []
        for direction in range(0, len(directions)):
            In_Wakek = []
            khub = turbines[j].XLocation[direction]                 #turbine k's X Location
            krad = turbines[j].RotorRad
            kleft = khub - krad                 #define left most point of rotor swept area of turbine k
            kright = khub + krad                #define right most point of rotor swept area of turbine k
            wake_d = []
            disty = []
        
            for i in range(0, initial_num): #upstream turbine
                if j != i:
                    hubheight = turbines[i].HubHeight 
                    alpha = 0.5 / (np.log(hubheight / z0))
                    turbines[i].alpha = alpha
                    y = turbines[i].YLocation[direction]
                    x = turbines[i].XLocation[direction]
                    dis = YLength - y                          #Sets dis to max downstream distance
                    Rr = turbines[i].RotorRad
                    r1 = (alpha * dis) + Rr #Calculate maximum ds wake radius
                    space1 = x + r1
                    space2 = x - r1
                    
                    if (kleft >= space2 and kleft <= space1) or (kright >= space2 and kright <= space1):
                        #if either rotor swept area points are within maximum rectangular wake
                        if turbines[i].YLocation[direction] < turbines[j].YLocation[direction]: 
                            Y = turbines[j].YLocation[direction]
                            dist = Y - y                  #distance between turbines
                            #print(dist) #code check
                            wake_rad = (alpha * dist) + Rr  #define radius of triangular wake                        
                            kz = turbines[j].HubHeight   #kz is the z coordinate of the rotor k's hub
                            jz = turbines[i].HubHeight   #jz is the z coordinate of the wake's center
                            cd = np.sqrt((x - khub) ** 2.0 + (jz - kz) ** 2.0)   #distance between the centerline of wake and rotor hub
                            if cd < (wake_rad + krad):  #if distance between centers is less than the sum of the two radii, the rotor swept area is in the wake
                                In_Wakek.append(i)
                                wake_d.append(wake_rad * 2.0)
                                disty.append(dist)
            upstrm.append(In_Wakek) 
            distanceus.append(disty)
            wakewidthds.append(wake_d)
        turbines[j].usturbines = [this for this in upstrm]
        #print('turbine j: ',j)
        #print(turbines[j].usturbines)
        #print('turbine: ',lastj)
        #print(turbines[lastj].usturbines)
        #lastj = j*1
        turbines[j].wakewidth = [this for this in wakewidthds]
        turbines[j].distance = [this for this in distanceus]
        
        #print('usturbines for ', i, ': ',turbines[i].usturbines)
        #print('wakewidth for ', j, ': ',turbines[j].wakewidth)
        #print('distances for ', j, ': ',turbines[j].distance)
        
    for i in range(0, initial_num):  
        dsds = []
        for d in range(0, direction):
            dsone = []
            for j in range(0, initial_num):    
                if j != i:
                    if i in turbines[j].usturbines[d]:
                        dsone.append(j)            
            dsds.append(dsone)
        turbines[i].dsturbines = dsds
        #print('dsturbines for ', j, ': ',turbines[j].dsturbines)

    #code check           
    #print(turbines[7].dsturbines[0])
    #print(turbines[7].dsturbinesrec)
    #print(turbines[7].usturbines[0])
    #print(turbines[7].usturbinesrec)
    #print(turbines[7].wakewidth)
    #print(turbines[7].distance)
    #print(turbines[4].dsturbines[0])
    #print(turbines[4].dsturbinesrec)
    #print(turbines[4].usturbines[0])
    #print(turbines[4].usturbinesrec)
    #print(turbines[4].wakewidth)
    #print(turbines[4].distance)
                   
    #Now that we know which turbines are downstream of others, calculate the percentage of the rotor swept area that is within the wake                        
    for i in range(0, initial_num):
        complete_percent = []
        for wd in range(0, len(directions)):
            parpercent = []
            overlap_flag = 0
            kz = turbines[i].HubHeight   #turbine k's hub height (z-location)
            kx = turbines[i].XLocation[wd]   #turbine k's x location
            #ky = turbines[i].YLocation[wd]   #turbine k's y location
            krad = turbines[i].RotorRad  #turbine k's rotor radius
            
            if len(turbines[i].usturbines[wd]) == 1:        #if the turbine has one upstream turbine
                j = turbines[i].usturbines[wd][0]
                jz = turbines[j].HubHeight                  #z coordinate of the wake's center
                jx = turbines[j].XLocation[wd]                  #x coordinate of the wake's center
                #jy = turbines[j].YLocation[wd]                  #y coordinate of the wake's center
                jwakerad = (turbines[i].wakewidth[wd][0]) / 2.0   #radius of wake width
                dist = turbines[i].distance[wd][0]
                cd = np.sqrt(((jx-kx) ** 2.0) + ((jz - kz) ** 2.0))   #distance between centerline of wake and rotor hub           
                int1_den = 2.0 * cd * krad
    
                if cd + krad <= jwakerad:                         #if dsturbine is completely in usturbine wake, overlap = 100%
                    parpercent.append(1.0)
                    #print('works')
                
                elif cd + jwakerad <= krad:                 #if the wake is fully encompassed by the rotor diameter    
                    wakearea = np.pi * (jwakerad ** 2.0)
                    percentwake = wakearea / (np.pi * (krad ** 2.0))
                    parpercent.append(percentwake)
    
                else:
                    integrand1 = ((cd ** 2.0) + (krad ** 2.0) - (jwakerad ** 2.0)) / int1_den
                    #print(integrand1)
                    int2_den = 2.0 * cd * jwakerad                
                    integrand2 = ((cd ** 2.0) + (jwakerad ** 2.0) - (krad ** 2.0)) / int2_den
                    #print(integrand2) 
                    q = (krad ** 2.0) * (np.arccos(integrand1)) 
                    b = (jwakerad ** 2.0) * (np.arccos(integrand2))
                    c = 0.5 * np.sqrt((-cd + krad + jwakerad) * (cd + krad - jwakerad) * (cd - krad + jwakerad) * (cd + krad + jwakerad))
                    AOverlap = q + b - c
                    RSA = ((np.pi) * (krad ** 2.0))
                    z = AOverlap / RSA
                    parpercent.append(z)      #percentage of RSA that has wake interaction 
                
           
            elif len(turbines[i].usturbines[wd]) == 2:            #if the turbine has two upstream turbines
                first = turbines[i].usturbines[wd][0]
                second = turbines[i].usturbines[wd][1]
                firstx = turbines[first].XLocation[wd]
                firstz = turbines[first].HubHeight
                firstrad = turbines[i].wakewidth[wd][0] / 2.0
                secondx = turbines[second].XLocation[wd]
                secondz = turbines[second].HubHeight
                secondrad = turbines[i].wakewidth[wd][1] / 2.0
                cd = np.sqrt(((firstx - secondx) ** 2.0) + ((firstz - secondz) ** 2.0))   #distance between the centerline of wake and rotor hub
    
                if cd > (firstrad + secondrad):     #if wakes do not overlap at all within the rotor swept area
                    #m = []
                    overlap_flag = 1
                    for q in range(0, len(turbines[i].usturbines[wd])):
                        j = turbines[i].usturbines[wd][q]
                        jz = turbines[j].HubHeight             #z coordinate of the wake's center
                        jx = turbines[j].XLocation[wd]             #x location of the wake's center
                        #jy = turbines[j].YLocation[wd]             #y location of the wake's center
                        jwakerad = (turbines[i].wakewidth[wd][q]) / 2.0
                        dist = turbines[i].distance[wd][q]
                        cd = np.sqrt(((jx - kx) ** 2.0) + ((jz - kz) ** 2.0))     #distance between the centerline of wake and rotor hub
                           
    
                        if cd + krad <= jwakerad:
                            parpercent.append(1.0)
    
                        elif cd + jwakerad <= krad:           #if the wake is fully encompassed by the rotor diameter
                            wakearea = np.pi * (jwakerad ** 2.0)
                            percentwake = wakearea / (np.pi * (krad ** 2.0))
                            parpercent.append(percentwake)
                            
                        else:
                            integrand1 = ((cd ** 2.0) + (krad ** 2.0) - (jwakerad ** 2.0)) / (2.0 * cd * krad)
                            integrand2 = ((cd ** 2.0) + (jwakerad ** 2.0) - (krad ** 2.0)) / (2.0 * cd * jwakerad)             
                            d = (krad ** 2.0) * (np.arccos(integrand1)) 
                            b = (jwakerad ** 2.0) * (np.arccos(integrand2))
                            c = 0.5 * np.sqrt((-cd + krad + jwakerad) * (cd + krad - jwakerad) * (cd - krad + jwakerad) * (cd + krad + jwakerad))
                            AOverlap = d + b - c
                            #AOverlap = ((krad ** 2.0) * np.arccos(integrand1)) + ((jwakerad ** 2.0) * np.arccos(integrand2)) - 0.5 * np.sqrt(abs((-cd + krad + jwakerad) * (cd + krad - jwakerad) * (cd - krad + jwakerad) * (cd + krad + jwakerad)))
                            RSA = np.pi * (krad ** 2.0)
                            z = AOverlap / RSA
                            parpercent.append(z)      #percentage of RSA that has wake interaction
                    
            if len(turbines[i].usturbines[wd]) >= 2 and overlap_flag != 1:      #if there are at least 2 upstream turbines whose wakes overlap, discretize the RSA and evaluate each point
                Discretize_RSA(i, wd)
                #print(i)
            complete_percent.append(parpercent)
        turbines[i].percent = complete_percent
#Code Check
#Compute_Wake(initial_num, z0, U0, Zref, alphah, ro, aif)                   
                  
    #calculate wind speed for each downstream turbine based on downstream distance
    sorted_order = [(jj,ii.YLocation[0]) for jj,ii in enumerate(turbines)]
    sorted_order = sorted(sorted_order, key=lambda x: x[1])
    sorted_order = [ii[0] for ii in sorted_order]
    #print('sorted order')
    #print(sorted_order)
    wdsp_list = [[0.] * len(U0) for ii in range(initial_num)]
    windspeeds = [wdsp_list for ii in range(len(turbines[0].XLocation))]
    # windspeeds = [[] for i in range(len(turbines[0].XLocation))]
    for wd in range(len(turbines[0].XLocation)):
        wdsp_byturb = []
        analysis_order = [(i, turbines[i].YLocation[wd]) for i in range(initial_num)]
        analysis_order.sort(key=lambda x: x[1])
        analysis_order = [i[0] for i in analysis_order]
        for k in analysis_order:
            wdsp = []
            for u0i in range(0, len(U0)):
                if len(turbines[k].usturbines[wd]) == 0:        #if turbine has no upstream turbines, 
                       #INCORPORATE POWER LAW
                    hubheight = turbines[k].HubHeight
                    Uz = U0[u0i] * ((hubheight / Zref) ** alphah)         #corrects wind speed for hub height
                    wdsp.append(Uz)
                   
                elif len(turbines[k].usturbines[wd]) == 1:        #if turbine has 1 upstream turbine
                    total = 0.0
                    #USturb = turbines[k].usturbines[wd][0]
                    #USht = turbines[USturb].HubHeight
                    x = turbines[k].distance[wd][0]
                    hubheight = turbines[k].HubHeight
                    alpha = (0.5 / np.log(hubheight / z0))
                    Rr = turbines[k].RotorRad
                    
                    #Grady Model 
                    r1 = Rr * np.sqrt((1-aif) / (1 - 2*aif))
                    EWU = U0[u0i] * (1 - (2*aif)/((1+alpha*(x/r1))**(2)))
                    Uz = EWU * ((hubheight / Zref) ** alphah)
                    #print(turbines[k].percent[wd][0])
                    portion = Uz * turbines[k].percent[wd][0]
                    remainder = U0[u0i] * (1.0 - turbines[k].percent[wd][0]) * ((hubheight / Zref) ** alphah)
                    total = portion + remainder                 #weighted average of windspeeds
                    wdsp.append(total)   
                    
                elif len(turbines[k].usturbines[wd]) == 2 and len(turbines[k].percent[wd]) != 0:      #if the turbine has two upstream turbines whose wakes do not overlap
                    portion = 0.0
                    total = 0.0
                    for j in range(0, len(turbines[k].usturbines[wd])):
                        x = turbines[k].distance[wd][j]
                        #USturb = turbines[k].usturbines[wd][j]
                        hubheight = turbines[k].HubHeight
                        alpha = 0.5 / np.log(hubheight / z0)
                        Rr = turbines[k].RotorRad
                        r1 = Rr * np.sqrt((1 - aif) / (1 - 2 * aif))
                        EWU = U0[u0i] * (1 - (2*aif)/((1+alpha*(x/r1))**(2)))
                        Uz = EWU * ((hubheight / Zref) ** alphah)
                        portion += Uz * turbines[k].percent[wd][j]
                    remainder = U0[u0i] * (1.0 - turbines[k].percent[wd][0] - turbines[k].percent[wd][1]) * ((hubheight / Zref) ** alphah)
                    total = portion + remainder                 #weighted average of windspeeds
                    wdsp.append(total)
                 
                elif len(turbines[k].usturbines[wd]) >= 2 and len(turbines[k].percent[wd]) == 0:      #turbine has at least two upstream turbines whos wakes overlap
                    coordWS = []
                    usturbcoord = [[] for i in range(len(turbines[k].xcoords[wd]))]
                    for i in range(0, len(turbines[k].xcoords[wd])):        #xcoords created in Discretize_RSA
                        decWS = []
                        xval = turbines[k].xcoords[wd][i]
                        zval = turbines[k].zcoords[wd][i]
                        khub = turbines[k].HubHeight
                        #alpha = 0.5 / np.log(zval / z0)
                        Rr = turbines[k].RotorRad
                        #r1 = Rr * np.sqrt((1.0 - aif) / (1.0 - 2.0 * aif))
                        for j in range(0, len(turbines[k].usturbines[wd])):
                            x = turbines[k].distance[wd][j]
                            US = turbines[k].usturbines[wd][j]
                            r2 = turbines[k].wakewidth[wd][j] / 2.0
                            #print('r2: ',r2)
                            xc = turbines[US].XLocation[wd]         #'c' for centerline of wake
                            #yc = turbines[US].YLocation[wd]
                            zhubc = turbines[US].HubHeight
                            xturb = xval * 1.
                            #yturb = turbines[k].YLocation[winddir]
                            zhubturb = zval * 1.
                            rt2 = abs(zhubturb - zhubc)        #height of the triangular portion of the chord area in z
                            rt1 = abs(xturb - xc)              #height of the triangluar portion of the chord area in x
                            space = np.sqrt((rt2 ** 2) + (rt1 ** 2))      #distance between wake center and discritized point
                            if space <= r2:        #if point is within wake
                                Rr = turbines[k].RotorRad
                                alpha = 0.5 / np.log(zval / z0)
                                #Grady's a
                                r1 = Rr * np.sqrt((1 - aif) / (1 - 2 * aif))
                                Uz = U0[u0i] * (1 - (2*aif)/((1+alpha*(x/r1))**(2)))
                                decWS.append(Uz)
                                usturbcoord[i].append(US)
                        #print('i: ',i)
                        #print('usturbcoord: ',usturbcoord)        
                        coordui = 0.0        
                        if len(decWS) != 0:
        
                            if len(decWS) == 1:         #if the point only has one wake acting on it
                                coordui = decWS[0] * ((zval / Zref) ** alphah)
                                coordWS.append(coordui)
                                
                            elif len(decWS) > 1:          #if the pint has more than one wake acting on it
                                tally = 0.0
                                for l in range(0, len(decWS)):
                                    u = decWS[l]
                                    tally += ((1.0 - (u / U0[u0i])) ** 2.0)
                                #print('should be pos: ',np.sqrt(tally),': turbine ', k)
                                coordui = U0[u0i] * (1 - (np.sqrt(tally))) * ((zval / Zref) ** alphah)
                                coordWS.append(coordui)
                                
                        else:               #if the point has no wakes acting on it
                            Uz = U0[u0i] * ((zval / Zref) ** alphah)
                            coordui = Uz
                            coordWS.append(coordui)
                            #print('error if grid: turbine ', k)                            
                            
                    #Sum discretized wind speeds
                    #nested wake provision
                    if set([usturbcoord[0] == i for i in usturbcoord]) == {True} and nwp == True:
                        #find index of closest upstream wake
                        usindex = turbines[k].usturbines[wd][turbines[k].distance[wd].index(min(turbines[k].distance[wd]))]
                        #print("analyzing turbine: ",k)
                        #print('only reducing speed from turbine: ',usindex)
                        x = min(turbines[k].distance[wd])
                        #print(x)
                        hubheight = turbines[k].HubHeight
                        alpha = (0.5 / np.log(hubheight / z0))
                        Rr = turbines[k].RotorRad
                        
                        #Grady Model 
                        r1 = Rr * np.sqrt((1.-aif) / (1. - 2.*aif))
                        EWU = windspeeds[wd][usindex][u0i] * (1 - (2.*aif)/((1.+alpha*(x/r1))**(2)))
                        #print(EWU)
                        wdsp.append(EWU * ((hubheight / Zref) ** alphah))
                        
                    else: #no nested wake
                        tally2 = 0.0
                        percentage = 1.0 / float(len(coordWS))
                        #print(coordWS)
                        for f in range(0, len(coordWS)):
                            tally2 += percentage * coordWS[f]
            
                        d = len(coordWS)
                        wdsp.append(tally2)
            wdsp_byturb.append(wdsp)
        order_wdsp = list(zip(analysis_order, wdsp_byturb))
        order_wdsp.sort(key=lambda x: x[0])
        wdsp_byturb = [ii[1] for ii in order_wdsp]
        windspeeds[wd] = wdsp_byturb
        turbines[k].ui = wdsp
        turbines[k].windspeeds = wdsp
    # print(windspeeds)
    #calculate power developed for each turbine
    for i in range(0, initial_num):
        pwr = []
        rorad = turbines[i].RotorRad
        Area = (rorad ** 2.0) * np.pi
        for wd in range(0, len(turbines[i].XLocation)):
            for ws in range(len(U0)):
                #power curve consideration
                #incorporating power curve suggested by Pat, June 10th
                if windspeeds[wd][i][ws] < Uref and windspeeds[wd][i][ws] >= Cutin: #Calculate power for effective windspeeds between 3 and 11 m/s
                    temp1 = 0.5 * ro * Area * (windspeeds[wd][i][ws] ** 3.0) * Cp * Cf/1000.
                    p1 = temp1 * probwui[wd][ws]
                    pwr.append(p1)
        
                elif windspeeds[wd][i][ws] < Cutin:        #cut in speed is 3 m/s
                    pwr.append(0.0)
        
                elif windspeeds[wd][i][ws] >= Uref:      #constant for 11.5 m/s rated power and above
                    temp1 = 0.5 * ro * Area * (Uref ** 3.0) * Cp * Cf/1000.
                    p1 = temp1 * probwui[wd][ws]
                    pwr.append(p1)
                
                #temp1 = 0.5 * ro * Area * (windspeeds[wd][i][ws] ** 3.0) * Cp * Cf/1000.
                #p1 = temp1 * probwui[wd][ws]
                #pwr.append(p1)
        turbines[i].Power = [this_power for this_power in pwr]
################################################################################################################
def Compute_Cost(initial_num, ro, yrs, WCOE, condition, depth):
    cost = 0.0
    Capital_Cost = 0
    Cabling_Cost = 0
    Mooring_Cost = 0
    O_M = 0
    Substation_Cost = 0
    #Installation_Cost = 977620.0
    Leasing_Cost = 0
    
    for i in range (0, initial_num):
        rorad = turbines[i].RotorRad
        Area = (rorad ** 2.0) * np.pi
        turbines[i].Area = Area             #define Area for each turbine
        cf = 0.4    #capcity factor
        Prated = 0.5 * ro * Area * (11.5 ** 3.0) * 0.5      #calculate rated power in Watts
        turbines[i].Prated = Prated
        Capital_Cost = Prated * 1.48
        turbines[i].Capital_Cost = Capital_Cost
        O_M = (Prated * 133.0 / 1000.0) * yrs
        turbines[i].O_M = O_M
        #Substation_Cost = 0.26284 * Prated     #From Myhr
        #Substation_Cost = 0.662269 * Prated    #From JEDI
        Substation_Cost = 0.02 * Prated
        turbines[i].Substation_Cost = Substation_Cost
        Leasing_Cost = (Prated / 1000.0) * 8760.0 * cf * WCOE * (8.0 * 0.02 + (yrs - 8.0) * 0.04)
        turbines[i].Leasing_Cost = Leasing_Cost     #(8.0 * 0.02 + (yrs - 8.0) * 0.04) is r, r = 0.02 for first 8 years, 0.04 for rest of life of project
        Mooring_Cost = 4.0 * (140148.0 + 274.0 * depth)
        turbines[i].Mooring_Cost = Mooring_Cost
        #Installation_Cost = turbines[i].Installation_Cost
        turbines[i].Cost = turbines[i].Capital_Cost + turbines[i].O_M + turbines[i].Substation_Cost + turbines[i].Leasing_Cost + turbines[i].Mooring_Cost + turbines[i].Installation_Cost   #cost related to Prated
        cost += turbines[i].Cost

    #REAL CABLE COST
    d_t, networks = calcicl(initial_num)
    #d_t = 2. #DUMMY FOR TESTING CODE
    d_s = 16.0 #assume 16 km offshore (10 mi)
    Cabling_Cost = d_t * 312000 + d_s * 492000
    cost += Cabling_Cost      #66210000.0 is from substation cost
    cost += 2000000.0
    #print('d_t:',d_t)    
    return cost
 
#Compute_Cost(initial_num, ro, yrs, WCOE, condition)   
 
###############################################################################
def Eval_Objective(initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE,
                   condition, aif):
    Compute_Wake(initial_num, z0, U0, Zref, alphah, ro, aif)
    #cost = Compute_Cost(initial_num, ro, yrs, WCOE, condition, depth)
    #cost = 0.0
    global nemo
    nemo += 1

    objective = 0.
    # tot_power = 0.
    ws_combos = len(probwui)*len(probwui[0])
    for j in range(0, initial_num):
        for ws in range(ws_combos):
            objective -= turbines[j].Power[ws]
            #objective -= turbines[j].Power[ws]
            #windspeeds += turbines[i].ui
    # objective = 50. / tot_power
    #print(objective)
    return objective
###############################################################################
def Clear_Vectors():
    for i in range(0, initial_num):
        turbines[i].percent = [[] for i in range(len(directions))]
        turbines[i].distance = [[] for i in range(len(directions))]
        turbines[i].dsturbines = [[] for i in range(len(directions))]
        turbines[i].wakewidth = [[] for i in range(len(directions))]
        turbines[i].usturbines = [[] for i in range(len(directions))]
        #turbines[i].dsturbinesrec[winddir] = []
        #turbines[i].usturbinesrec[winddir] = []
        turbines[i].windspeeds = [[] for i in range(len(directions))]
###############################################################################
def Rand_Vector(initial_num):
    random_vec = []
    for i in range(0, initial_num):    
        random_vec.append(i)
        
    #shuffle elements by randomly exchanging each with one other
    for i in range(0, len(random_vec)):
        r = random.randint(0, len(random_vec)-1)  #select random value in vector
        temp = random_vec[i]
        random_vec[i] = random_vec[r]
        random_vec[r] = temp
    return random_vec
###############################################################################
def Pattern_Search(init_step, minstep, random_vec, initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, num_pops, max_pop_tries, hstep, i, hstepmin, rradmin, rradmax, rstep, hubmin, hubmax, rstepmin, aif, directions):
    eval_ct = 0
    with open('layouts_E.txt', 'w') as layouts_file: #clear layoutfile
        layouts_file.close()
        
    Clear_Vectors()
    #nomove = Eval_Objective(initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, aif)

    for h in range(0, 1):
        nomove = Eval_Objective(initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, aif)
        eval_ct += 1
        #for ii,jj in enumerate(turbines):
        #    print('turbine: ',ii)
        #    print(jj.XLocation)
        #    print(jj.YLocation)
        #    print(jj.ui)
        #    print(jj.usturbines)
            
        #print('objective eval: ', nomove)
        #print('h= ', h)
        step2 = init_step * mesh_width
        while step2 >= minstep:
            print(step2)
            #print_graph()
            random_vec = Rand_Vector(initial_num)    #creates a randomly ordered vector of turbines
            for j in range(0, len(random_vec)):
                # cter >= 32:
                #    break
                #cter += 1
                i = random_vec[j]
                turbines[i].Stopped = 0
                #print('Turbine ', i,' is being tested.', nomove)
                #print(Eval_Objective(initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, aif))
                flag = 0
                innerflag = 0
                transflag = 0
                ###comment out when finished debugging!!!
                #step2 = minstep - 2.
                #break
                #develop preliminary objective for comparison purposes
                
                Clear_Vectors()
                # print('The nomove value for turbine ', i, ' is ', nomove)
                #print('step size: ',step2)
                #('stepped into while loop')
                if innerflag == 0 and flag == 0:
                    # move 1 was just unsucessfully attempted
                    transflag = Translation_Y(step2, i, directions)
                    CHECK2 = 0
                    if transflag == 1:      #if the translation moved the turbine out of bounds, go to next translation
                        innerflag = 1       #move2 was attempted
                        # print('turbine not moved up. out of bounds')
                    else:
                        CHECK2 = Check_Interference(i)
                        if CHECK2 == 1:         #if interference occurs, move the turbine back, go to next translation
                            Translation_Y(-step2, i, directions)
                            innerflag = 1
                            # print('turbine not moved up. interferenec')
                        #elif check_layout([i.XLocation[0] for i in turbines],[i.YLocation[0] for i in turbines]) == True:
                        #        Translation_Y(-step2, i, directions)
                        #        innerflag = 1
                        #        print('turbine not moved up. Layout already attempted')
                        else:       #if there is no interference, evaluate and store
                            move2 = Eval_Objective(initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, aif)
                            eval_ct += 1
                            Clear_Vectors()
                            if move2 >= nomove:       #if evaluation is worse than initial, move back, go to next translation
                                Translation_Y(-step2, i, directions)
                                innerflag = 1
                                # print('turbine not moved up. worse eval: ', move2)
                            else:       #evaluation is better, keep move, go to next turbine
                                flag = 1
                                nomove = move2 * 1.
                                # print('turbine ', i, ' moved up.', move2)
                                #print(nomove)
                                #HubHeight_Search(hstep, i, hstepmin, XLocation, YLocation, initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, rradmin, rradmax, rstep, hubmin, hubmax, rstepmin, aif)
                            
                if innerflag == 1 and flag == 0:        #move 2 was just unsucessfully attempted
                    transflag = Translation_X(-step2, i, directions)
                    CHECK2 = 0
                    if transflag == 1:      #if the translation moved the turbine out of bounds, go to next translation
                        innerflag = 2       #move3 was attempted
                        # print('turbine not moved left. out of bounds')
                    else:
                        CHECK2 = Check_Interference(i)
                        if CHECK2 == 1:         #if interference occurs, move the turbine back, go to next translation
                            Translation_X(step2, i, directions)
                            innerflag = 2
                            # print('turbine not moved left. interferenec')
                        #elif check_layout([i.XLocation[0] for i in turbines],[i.YLocation[0] for i in turbines]) == True:
                        #    Translation_X(step2, i, directions)
                        #    innerflag = 2
                        #    print('turbine not moved up. Layout already attempted')
                        else:       #if there is no interference, evaluate and store
                            move3 = Eval_Objective(initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, aif)
                            eval_ct += 1
                            Clear_Vectors()
                            if move3 >= nomove:       #if evaluation is worse than initial, move back, go to next translation
                                Translation_X(step2, i, directions)
                                innerflag = 2
                                # print('turbine not moved left. worse eval: ', move3)
                            else:       #evaluation is better, keep move, go to next turbine
                                flag = 1
                                nomove = move3 * 1.
                                # print('turbine ', i, ' moved left.', move3)
                                #print(nomove)
                                #HubHeight_Search(hstep, i, hstepmin, XLocation, YLocation, initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, rradmin, rradmax, rstep, hubmin, hubmax, rstepmin, aif)
                            
                if innerflag == 2 and flag == 0:        #move 3 was just unsucessfully attempted
                    transflag = Translation_Y(-step2, i, directions)
                    CHECK2 = 0
                    if transflag == 1:      #if the translation moved the turbine out of bounds, go to next translation
                        innerflag = 3       #move3 was attempted
                        # print('turbine not moved down. out of bounds')
                    else:
                        CHECK2 = Check_Interference(i)
                        if CHECK2 == 1:         #if interference occurs, move the turbine back, go to next translation
                            Translation_Y(step2, i, directions)
                            innerflag = 3
                            # print('turbine not moved down. interference')
                        #elif check_layout([i.XLocation[0] for i in turbines],[i.YLocation[0] for i in turbines]) == True:
                        #    Translation_Y(step2, i, directions)
                        #    innerflag = 3
                        #    print('turbine not moved up. Layout already attempted')
                        else:       #if there is no interference, evaluate and store
                            move4 = Eval_Objective(initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, aif)
                            eval_ct += 1
                            Clear_Vectors()
                            if move4 >= nomove:       #if evaluation is worse than initial, move back, go to next translation
                                Translation_Y(step2, i, directions)
                                innerflag = 3
                                # print('turbine not moved down. worse eval: ', move4)
                            else:       #evaluation is better, keep move, go to next turbine
                                flag = 1
                                nomove = move4 * 1.
                                #print(nomove)
                                # print('turbine ', i, ' moved down.', move4)
                                #HubHeight_Search(hstep, i, hstepmin, XLocation, YLocation, initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, rradmin, rradmax, rstep, hubmin, hubmax, rstepmin, aif)
                if innerflag == 3 and flag == 0:
                    transflag = Translation_X(step2, i, directions)     #move the turbine one step right
                    CHECK2 = 0
                    if transflag == 1:          #if the translation moved the turbine out of bounds, go to next translation
                        innerflag = 4           #signifies move 1 was attempted
                        # print('Turbine not moved right. out of bounds')
                
                    else:       #if there is the turbine is in bounds, evaluate and store
                        CHECK2 = Check_Interference(i)
                        if CHECK2 == 1:         #if interference occurs, move the turbine back , go to next translation
                            Translation_X(-step2, i, directions)
                            innerflag = 4
                            # print('turbine not moved right, interference')
                        #elif check_layout([i.XLocation[0] for i in turbines],[i.YLocation[0] for i in turbines]) == True:
                        #    Translation_X(-step2, i, directions)
                        #    innerflag = 4
                        #    print('turbine not moved up. Layout already attempted')
                        else:       #if there is no interference, evaluate and store
                            move1 = Eval_Objective(initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, aif)
                            eval_ct += 1
                            Clear_Vectors()
                            if move1 >= nomove:             #if evaluation is worse than initial, move back, go to next translation
                                Translation_X(-step2, i, directions)
                                innerflag = 4
                                # print('Turbine not moved right. worse eval: ', move1)
                            else:
                                flag = 1           #signifies movement was kept
                                nomove = move1 * 1.
                                # print('turbine ', i, ' moved right.', move1)
                                #print(nomove)
                            #HubHeight_Search(hstep, i, hstepmin, XLocation, YLocation, initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, rradmin, rradmax, rstep, hubmin, hubmax, rstepmin, aif)
            
                if innerflag == 4 and flag == 0:        #translation at this step size has resulted in no moves for this turbine
                    turbines[i].Stopped = 1
                    #HubHeight_Search(hstep, i, hstepmin, XLocation, YLocation, initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, rradmin, rradmax, rstep, hubmin, hubmax, rstepmin, aif)
    
            exit_css = 0        #exit current step size
            for i in range(0, initial_num):
                exit_css += turbines[i].Stopped
            #print(exit_css)

            if exit_css == initial_num:
                #break   
                #print_graph()
                #all turbines have stopped moving at this step size, halving step size.
                #find worst performing turbine and randomly assign elsewhere
                #print("No moves at step size ", step2, " are possible. Popping weakest turbine.")
                for b in range(0, num_pops):
                    min_power = 5000000.     #initialized to first turbine power output
                    random_vec2 = Rand_Vector(initial_num)    #creates a randomly ordered vector of turbines
                    for j in range(0, initial_num):
                        randorder = random_vec2[j]
                        Power = sum(turbines[randorder].Power)                        
                        if Power < min_power:
                            min_power = Power
                            min_turb = randorder
                                
                    # print('The weakest turbine is turbine ', min_turb, ' with power currently at ', min_power)
                    #start_eval = Eval_Objective(initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, aif)

                    Clear_Vectors()
                    initialx = turbines[min_turb].XLocation[0]
                    initialy = turbines[min_turb].YLocation[0]
                        
                    checkx = 0
                    k = 0
                    flag = 0
                    while flag == 0 and k < max_pop_tries:
                        checkx = 0
                        while checkx != 1:      #will try random locations until one has no interference
                            CHECK2 = 0
                            x_opts = int(site_x / mesh_width) + 1
                            y_opts = int(site_y / mesh_width) + 1
                            turbines[min_turb].XLocation[0] = int(random.uniform(0, x_opts)) * mesh_width
                            turbines[min_turb].YLocation[0] = int(random.uniform(0, y_opts)) * mesh_width
                            CHECK2 = Check_Interference(min_turb)
                            if CHECK2 != 1:         #No interference
                                checkx = 1          #place turbine and exit poping loop
                                for j in range(1, len(directions)):
                                    theta = directions[j]
                                    turbines[min_turb].XLocation[j] = (turbines[min_turb].XLocation[0]* np.cos(theta)) - (turbines[min_turb].YLocation[0] * np.sin(theta))
                                    turbines[min_turb].YLocation[j] = (turbines[min_turb].XLocation[0] * np.sin(theta)) + (turbines[min_turb].YLocation[0] * np.cos(theta))
                                # print('Turbine ', min_turb, ' has moved to a new location.')
                            else:
                                turbines[min_turb].XLocation[0] = initialx
                                turbines[min_turb].YLocation[0] = initialy
                                #print('Turbine cannot be relocated without interference, trying agian.')
                            

                        new_eval = Eval_Objective(initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, aif)
                        eval_ct += 1
                        Clear_Vectors()
                        if new_eval < nomove:
                            flag = 1
                            nomove = new_eval * 1. #keep eval
                            #HubHeight_Search(hstep, i, hstepmin, XLocation, YLocation, initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, rradmin, rradmax, rstep, hubmin, hubmax, rstepmin, aif)
                            #print('Move has improved the evaluation: ',new_eval)
                            #print(nomove)
                        else:
                            turbines[min_turb].XLocation[0] = initialx
                            turbines[min_turb].YLocation[0] = initialy
                            for j in range(1, len(directions)):
                                    theta = directions[j]
                                    turbines[min_turb].XLocation[j] = (initialx* np.cos(theta)) - (initialy * np.sin(theta))
                                    turbines[min_turb].YLocation[j] = (initialx * np.sin(theta)) + (initialy * np.cos(theta))
                            #print('Move did not improve evaluation. Trying new moves.')
                        k += 1
                        
                    #print('pops till acceptance= ', k)
                                
                if init_step > 1. and init_step < 2.: #make sure it tries one mesh distance
                    init_step = int(1)
                else:
                    init_step = int(init_step / 2.)
                step2 = init_step * mesh_width
                #print('step sized is being reduced to: ', step2)
                        
            #else:
                #print('Turbines were not all placed at step size ', step2, '. Repeating.')
    return eval_ct
#################################################################################################################
def HubHeight_Search(hstep, i, hstepmin, XLocation, YLocation, initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, rradmin, rradmax, rstep, hubmin, hubmax, rstepmin, aif):            
    step = hstep * 1.
    #zlocal = turbines[i].ZLocation
    while step > hstepmin:
        zhubheight = turbines[i].HubHeight
        hflag = 0       #hub has not been moved
        startval = Eval_Objective(initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, aif)

        #print('starting objective: ', startval)
        Clear_Vectors()
        #First Move
        if hflag == 0:
            if (zhubheight + step) <= hubmax:
                turbines[i].OldHubHeight = zhubheight       #store current hub height prior to any changes
                oldrr = turbines[i].RotorRad 
                turbines[i].HubHeight = zhubheight + step       #increase hub height by one step size
                RotorRadius_Search(i, rradmin, rradmax, rstep, initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, rstepmin, aif)                
                tempspeed = turbines[i].ui
                Obj1 = Eval_Objective(initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, aif)

                Clear_Vectors()
                if Obj1 >= startval:        #if change in hub height worsens evaluation, change back to old hub height, go to next direction
                    turbines[i].HubHeight = turbines[i].OldHubHeight
                    turbines[i].RotorRad = oldrr
                    turbines[i].ui = tempspeed
                    #print('new hh worse, back to ', turbines[i].HubHeign, ' and ', turbines[i].RotorRad)
                    
                else:
                    hflag = 1       #hub has been moved
                    #print('HubHeight up to ', turbines[i].HubHeign, ', new objective: ', Obj1)
        
        #Second Move            
        if hflag == 0:      #if move 1 was not accepted
            if (zhubheight - step) >= hubmin:
                turbines[i].OldHubHeight = zhubheight       #store current hub height prior to any changes
                oldrr = turbines[i].RotorRad 
                turbines[i].HubHeight = zhubheight - step       #decrease hub height by one step size
                RotorRadius_Search(i, rradmin, rradmax, rstep, initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, rstepmin, aif)
                tempspeed = turbines[i].ui
                Obj2 = Eval_Objective(initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, aif)

                Clear_Vectors()
                if Obj2 >= startval:        #if change in hub height worsens evaluation, change back to old hub height, go to next direction
                    turbines[i].HubHeight = turbines[i].OldHubHeight
                    turbines[i].RotorRad = oldrr
                    turbines[i].ui = tempspeed
                    #print('new hh worse, back to ', turbines[i].HubHeigt, ' and ', turbines[i].RotorRad)
                    
                else:
                    hflag = 1       #hub has been moved
                    #print('HubHeight down to ', turbines[i].HubHeigt, ', new objective: ', Obj2)
                    
        #Third Move
        if hflag == 0:          #Second Move not accepted
            if (zhubheight + (step / 2.0)) <= hubmax:
                turbines[i].OldHubHeight = zhubheight       #store current hub height prior to any changes
                oldrr = turbines[i].RotorRad 
                turbines[i].HubHeight = zhubheight + (step / 2.0)       #increase hub height by one half-step size
                RotorRadius_Search(i, rradmin, rradmax, rstep, initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, rstepmin, aif)                
                tempspeed = turbines[i].ui
                Obj3 = Eval_Objective(initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, aif)
                
                Clear_Vectors()
                if Obj3 >= startval:        #if change in hub height worsens evaluation, change back to old hub height, go to next direction
                    turbines[i].HubHeight = turbines[i].OldHubHeight
                    turbines[i].RotorRad = oldrr
                    turbines[i].ui = tempspeed
                    #print('new hh worse, back to ', turbines[i].HubHeigt, ' and ', turbines[i].RotorRad)
                    
                else:
                    hflag = 1       #hub has been moved
                    #print('HubHeight up half to ', turbines[i].HubHeigt, ', new objective: ', Obj3)
                    
        #Fourth Move            
        if hflag == 0:      #if move 3 was not accepted
            if (zhubheight - (step / 2.0)) >= hubmin:
                turbines[i].OldHubHeight = zhubheight       #store current hub height prior to any changes
                oldrr = turbines[i].RotorRad 
                turbines[i].HubHeight = zhubheight - (step / 2.0)       #decrease hub height by one half-step size
                RotorRadius_Search(i, rradmin, rradmax, rstep, initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, rstepmin, aif)
                tempspeed = turbines[i].ui
                Obj4 = Eval_Objective(initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, aif)

                Clear_Vectors()
                if Obj4 >= startval:        #if change in hub height worsens evaluation, change back to old hub height, go to next direction
                    turbines[i].HubHeight = turbines[i].OldHubHeight
                    turbines[i].RotorRad = oldrr
                    turbines[i].ui = tempspeed
                    #print('new hh worse, back to ', turbines[i].HubHeigt, ' and ', turbines[i].RotorRad)
                    
                else:
                    hflag = 1       #hub has been moved    a
                    #print('HubHeight down to ', turbines[i].HubHeigt, ', new objective: ', Obj4)
                    
        if hflag == 0:      #if no moves were accepted, half step size
            #print('No changes in height were selected for the step size 1 = ', step, '. Halving step size.')
            step = step / 2.0
            
        #else:
            #print('Hub height has not remained constant through search, continuing at step size ', step, '.')        
#####################################################################################################################
def Discretize_RSA(i, winddir):
    XCenter = turbines[i].XLocation[winddir]
    ZCenter = turbines[i].HubHeight
    rad = turbines[i].RotorRad
    #turbines[i].xcoords[winddir] = []
    #turbines[i].zcoords[winddir] = []
    xcoords = [XCenter]
    zcoords = [ZCenter]
    #center row
    #center point
    for j in range(1, 5):
        xcoords.append(XCenter + (j * (rad / 4.0)))
        xcoords.append(XCenter - (j * (rad / 4.0)))
        zcoords.append(ZCenter)
        zcoords.append(ZCenter)
    '''    
    #next rows
    #+ in Z-Direction
    #center Point
    xcoords.append(XCenter)
    zcoords.append(ZCenter + (rad / 4.0))
    for j in range(1, 4):
        xcoords.append(XCenter + (j * (rad / 4.0)))
        xcoords.append(XCenter - (j * (rad / 4.0)))
        zcoords.append(ZCenter + (rad / 4.0))
        zcoords.append(ZCenter + (rad / 4.0))
        
    #- in Z-Direction
    #center Point
    xcoords.append(XCenter)
    zcoords.append(ZCenter - (rad / 4.0))
    for j in range(1, 4):
        xcoords.append(XCenter + (j * (rad / 4.0)))
        xcoords.append(XCenter - (j * (rad / 4.0)))
        zcoords.append(ZCenter - (rad / 4.0))
        zcoords.append(ZCenter - (rad / 4.0))
        
    #next rows
    #+ in Z-Direction
    #center Point
    xcoords.append(XCenter)
    zcoords.append(ZCenter + (rad / 2.0))
    for j in range(1, 4):
        xcoords.append(XCenter + (j * (rad / 4.0)))
        xcoords.append(XCenter - (j * (rad / 4.0)))
        zcoords.append(ZCenter + (rad / 2.0))
        zcoords.append(ZCenter + (rad / 2.0))
        
    #- in Z-Direction
    #center Point
    xcoords.append(XCenter)
    zcoords.append(ZCenter - (rad / 2.0))
    for j in range(1, 4):
        xcoords.append(XCenter + (j * (rad / 4.0)))
        xcoords.append(XCenter - (j * (rad / 4.0)))
        zcoords.append(ZCenter - (rad / 2.0))
        zcoords.append(ZCenter - (rad / 2.0))
        
    #next rows
    #+ in Z-Direction
    #center Point
    xcoords.append(XCenter)
    zcoords.append(ZCenter + (rad * (3.0 / 4.0)))
    for j in range(1, 3):
        xcoords.append(XCenter + (j * (rad / 4.0)))
        xcoords.append(XCenter - (j * (rad / 4.0)))
        zcoords.append(ZCenter + (rad * (3.0 / 4.0)))
        zcoords.append(ZCenter + (rad * (3.0 / 4.0)))
        
    #- in Z-Direction
    #center Point
    xcoords.append(XCenter)
    zcoords.append(ZCenter - (rad * (3.0 / 4.0)))
    for j in range(1, 3):
        xcoords.append(XCenter + (j * (rad / 4.0)))
        xcoords.append(XCenter - (j * (rad / 4.0)))
        zcoords.append(ZCenter - (rad * (3.0 / 4.0)))
        zcoords.append(ZCenter - (rad * (3.0 / 4.0)))
        
    #last points: Top Center & Bottom Center
    xcoords.append(XCenter)
    zcoords.append(ZCenter + rad)
    xcoords.append(XCenter)
    zcoords.append(ZCenter - rad)
    '''
    #print('xcoords: ',xcoords)
    dummyx = turbines[i].xcoords
    dummyy = turbines[i].zcoords
    dummyx[winddir] = [this for this in xcoords]
    dummyy[winddir] = [this for this in zcoords]
    turbines[i].xcoords = [this for this in dummyx]
    turbines[i].zcoords = [this for this in dummyy]
    #if i >= 13:
    #    print(turbines[13].xcoords)
#############################################################################################################
def RotorRadius_Search(i, rradmin, rradmax, rstep, initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, rstepmin, aif):
    #now that an optimal hub height has been found, perform sub-level EPS for rotor radius
    #Clear_Vectors()
    
    hubh = turbines[i].HubHeight
    #rrad = turbines[i].RotorRad         #initialized to 20m
    minrange = (hubh * (2.0 / 3.0)) / 2.0   #the minimum rotor diameter should be 2/3 of the hub height
    maxrange = hubh / 2.0       #the maximum rotor diameter should be 1x the hub height
    
    if minrange < rradmin:      #set min as larger of 19m OR ratio parameter
        minrange = rradmin * 1.
        #print('minimum range for rotor radius is ', minrange)
        
    if maxrange > rradmax:      #set max as smaller of 67 OR ratio parameter
        maxrange = rradmax * 1.
        #print('maximum range for rotor radius is ', maxrange)
        
    ranger = maxrange - minrange
    centerrad = minrange + (ranger / 2.0)       #center value of range of rotor radii
    #currentr = turbines[i].RotorRad

    turbines[i].RotorRad = centerrad        #set rotor radius to cener valid value
    
    step2 = rstep * 1.
    while step2 > rstepmin:
        rflag = 0               #radius has not been changed
        startval = Eval_Objective(initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, aif)
        #print('starting objective: ' startval)
        Clear_Vectors()
        
        #First Move
        if rflag == 0:
            if (turbines[i].RotorRad + step2) <= maxrange:
                turbines[i].OldRotorRad = turbines[i].RotorRad  #store current rotor radius prior to any changes
                turbines[i].RotorRad += step2   #increase rotor radius by one step size
                tempspeed = turbines[i].ui
                Obj1 = Eval_Objective(initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, aif)

                Clear_Vectors()
                if Obj1 >= startval:        #if change in rotor radius worsens evaluation, change back to old rotor radius, go to next direction
                    turbines[i].RotorRad = turbines[i].OldRotorRad
                    turbines[i].ui = tempspeed

                else:
                    rflag = 1       #move was accepted
                    #print('rr up to ', turbines[i].RotorRad, ', new objective: ' Obj1)
                    
        #Second Move
        if rflag == 0:      #if First Move was not accepted
            if (turbines[i].RotorRad - step2) >= minrange:
                turbines[i].OldRotorRad = turbines[i].RotorRad  #store current rotor radius prior to any changes
                turbines[i].RotorRad -= step2   #decrease rotor radius by one step size
                tempspeed = turbines[i].ui
                Obj2 = Eval_Objective(initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, aif)

                Clear_Vectors()
                if Obj2 >= startval:        #if change in rotor radius worsens evaluation, change back to old rotor radius, go to next direction
                    turbines[i].RotorRad = turbines[i].OldRotorRad
                    turbines[i].ui = tempspeed

                else:
                    rflag = 1       #move was accepted
                    #print('rr down to ', turbines[i].RotorRad, ', new objective: ' Obj2)
                    
        #Third Move   a
        if rflag == 0:      #if Second Move was not accepted
            if (turbines[i].RotorRad + (step2 / 2.0)) <= maxrange:
                turbines[i].OldRotorRad = turbines[i].RotorRad  #store current rotor radius prior to any changes
                turbines[i].RotorRad += (step2 / 2.0)   #increase rotor radius by one half-step size
                tempspeed = turbines[i].ui
                Obj3 = Eval_Objective(initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, aif)

                Clear_Vectors()
                if Obj3 >= startval:        #if change in rotor radius worsens evaluation, change back to old rotor radius, go to next direction
                    turbines[i].RotorRad = turbines[i].OldRotorRad
                    turbines[i].ui = tempspeed

                else:
                    rflag = 1       #move was accepted
                    #print('rr up to ', turbines[i].RotorRad, ', new objective: ' Obj3)
                    
        #Fourth Move
        if rflag == 0:      #if Third Move was not accepted
            if (turbines[i].RotorRad - (step2 / 2.0)) >= minrange:
                turbines[i].OldRotorRad = turbines[i].RotorRad  #store current rotor radius prior to any changes
                turbines[i].RotorRad -= (step2 / 2.0)   #decrease rotor radius by one half-step size
                tempspeed = turbines[i].ui
                Obj4 = Eval_Objective(initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, aif)

                Clear_Vectors()
                if Obj4 >= startval:        #if change in rotor radius worsens evaluation, change back to old rotor radius, go to next direction
                    turbines[i].RotorRad = turbines[i].OldRotorRad
                    turbines[i].ui = tempspeed

                else:
                    rflag = 1       #move was accepted
                    #print('rr down to ', turbines[i].RotorRad, ', new objective: ' Obj4)
                    
        if rflag == 0:      #if this flag is still zero at this point, no moves were taken at this step size, halve step size
            #print('No Changes in rotor radius were selected for the step size l = ', step2'. Halving step size.')
            step2 = step2 / 2.0
            
        #else:
            #print('Rotor Radius has not remained contant through search, continuing at step size ', step2, '.')  
###########################################################################################
def print_graph():   
    redcx = []
    redcy = []
    yellowcx = []
    yellowcy = []
    greencx = []
    greency = []
    bluecx = []
    bluecy = []
    redtrx = []
    redtry = []
    yellowtrx = []
    yellowtry = []
    greentrx = []
    greentry = []
    bluetrx = []
    bluetry = []
    redrhx = []
    redrhy = []
    yellowrhx = []
    yellowrhy = []
    greenrhx = []
    greenrhy = []
    bluerhx = []
    bluerhy = []
    redsqx = []
    redsqy = []
    yellowsqx = []
    yellowsqy = []
    greensqx = []
    greensqy = []
    bluesqx = []
    bluesqy = []
    noplotx = []
    noploty = []
    
    for i in range(0, initial_num):
        if turbines[i].HubHeight <= 60 and turbines[i].RotorRad <= 30:
            redcx.append(turbines[i].XLocation[0])
            redcy.append(turbines[i].YLocation[0])

        elif turbines[i].HubHeight <= 60 and turbines[i].RotorRad <= 40:
            yellowcx.append(turbines[i].XLocation[0])
            yellowcy.append(turbines[i].YLocation[0])

        elif turbines[i].HubHeight <= 60 and turbines[i].RotorRad <= 60:
            greencx.append(turbines[i].XLocation[0])
            greency.append(turbines[i].YLocation[0])

        elif turbines[i].HubHeight <= 60 and turbines[i].RotorRad > 60:
            bluecx.append(turbines[i].XLocation[0])
            bluecy.append(turbines[i].YLocation[0])

        elif turbines[i].HubHeight <= 80 and turbines[i].RotorRad <= 30:
            redtrx.append(turbines[i].XLocation[0])
            redtry.append(turbines[i].YLocation[0])

        elif turbines[i].HubHeight <= 80 and turbines[i].RotorRad <= 40:
            yellowtrx.append(turbines[i].XLocation[0])
            yellowtry.append(turbines[i].YLocation[0])

        elif turbines[i].HubHeight <= 80 and turbines[i].RotorRad <= 60:
            greentrx.append(turbines[i].XLocation[0])
            greentry.append(turbines[i].YLocation[0])

        elif turbines[i].HubHeight <= 80 and turbines[i].RotorRad > 60:
            bluetrx.append(turbines[i].XLocation[0])
            bluetry.append(turbines[i].YLocation[0])

        elif turbines[i].HubHeight <= 120 and turbines[i].RotorRad <= 30:
            redrhx.append(turbines[i].XLocation[0])
            redrhy.append(turbines[i].YLocation[0])

        elif turbines[i].HubHeight <= 120 and turbines[i].RotorRad <= 40:
            yellowrhx.append(turbines[i].XLocation[0])
            yellowrhy.append(turbines[i].YLocation[0])

        elif turbines[i].HubHeight <= 120 and turbines[i].RotorRad <= 60:
            greenrhx.append(turbines[i].XLocation[0])
            greenrhy.append(turbines[i].YLocation[0])

        elif turbines[i].HubHeight <= 120 and turbines[i].RotorRad > 60:
            bluerhx.append(turbines[i].XLocation[0])
            bluerhy.append(turbines[i].YLocation[0])

        elif turbines[i].HubHeight > 120 and turbines[i].RotorRad <= 30:
            redsqx.append(turbines[i].XLocation[0])
            redsqy.append(turbines[i].YLocation[0])
  
        elif turbines[i].HubHeight > 120 and turbines[i].RotorRad <= 40:
            yellowsqx.append(turbines[i].XLocation[0])
            yellowsqy.append(turbines[i].YLocation[0])

        elif turbines[i].HubHeight > 120 and turbines[i].RotorRad <= 60:
            greensqx.append(turbines[i].XLocation[0])
            greensqy.append(turbines[i].YLocation[0])

        elif turbines[i].HubHeight > 120 and turbines[i].RotorRad > 60:
            bluesqx.append(turbines[i].XLocation[0])
            bluesqy.append(turbines[i].YLocation[0])

        else:
            noplotx.append(turbines[i].XLocation[0])
            noploty.append(turbines[i].YLocation[0])

        #print('unplotted x ', noplotx)
        #print('unplotted y', noploty)

    #initial_num += 1
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.scatter(redcx, redcy, s=10, c='r', marker="o")
    ax1.scatter(yellowcx, yellowcy, s=10, c='y', marker="o")
    ax1.scatter(greencx, greency, s=10, c='g', marker="o")
    ax1.scatter(bluecy, bluecy, s=10, c='b', marker="o")
    ax1.scatter(redtrx, redtry, s=10, c='r', marker="^")
    ax1.scatter(yellowtrx, yellowtry, s=10, c='y', marker="^")
    ax1.scatter(greentrx, greentry, s=10, c='g', marker="^")
    ax1.scatter(bluetrx, bluetry, s=10, c='b', marker="^")
    ax1.scatter(redrhx, redrhy, s=10, c='r', marker="d")
    ax1.scatter(yellowrhx, yellowrhy, s=10, c='y', marker="d")
    ax1.scatter(greenrhx, greenrhy, s=10, c='g', marker="d")
    ax1.scatter(bluerhy, bluerhy, s=10, c='b', marker="d")
    ax1.scatter(redsqx, redsqy, s=10, c='r', marker="s")
    ax1.scatter(yellowsqx, yellowsqy, s=10, c='y', marker="s")
    ax1.scatter(greensqx, greensqy, s=10, c='g', marker="s")
    ax1.scatter(bluesqx, bluesqy, s=10, c='b', marker="s")
    #plt.axis([0, site_x, 0, site_y])
    
    for i in range(len(turbines)):
        ax1.annotate(i, (turbines[i].XLocation[0],turbines[i].YLocation[0]))
    plt.ylabel('Position (m)')
    plt.xlabel('Position (m)')
    plt.title(str('Optimization of ' + str(initial_num) + ' Turbines'))
    xmax = int(np.ceil(site_x / mesh_width))
    ymax = int(np.ceil(site_y / mesh_width))
    plt.xticks([i*mesh_width for i in range(xmax)])
    plt.yticks([i*mesh_width for i in range(ymax)])
    plt.show()
    #plt.savefig(str(str(initial_num) + 'turbinesWithDisc.png'), bbox_inches='tight')
#############################################################################################################
#calculating distance between turbines
def length(j, k):
    x1 = turbines[j].XLocation[0]
    x2 = turbines[k].XLocation[0]
    y1 = turbines[j].YLocation[0]
    y2 = turbines[k].YLocation[0]
    
    term1 = (x1 - x2) ** 2
    term2 = (y1 - y2) ** 2
    length = np.sqrt(term1 + term2)
    return length   
##############################################################################################################
# Find Closest Turbines
def closest(j, hood_size):
    dist = []
    turb = [n for n in range(0,initial_num)]
    turb.remove(j)
    for l in range(0, initial_num):
        if l != j:
            space = length(j, l)
            dist.append(space)

    current_trb = [j] * (initial_num - 1)
    Dist_Order = list(zip(current_trb, turb, dist))               #Create tuple to identify turbine number
    Dist_Order = sorted(Dist_Order, key=lambda x: x[2])
    closest = []
    for i in range(0, hood_size):
        closest.append(Dist_Order[i])
    
    turbines[j].hood = closest
        
    return closest
#####################################################################################################
def calcicl(initial_num):
    if initial_num > 10:
        hood_size = 10
    else:
        hood_size = initial_num - 1
    for i in range(0, initial_num):
        closest(i, hood_size)
    
    #create list of closest turbines and distances
    all_opt = []
    rand_vec = Rand_Vector(initial_num)
    for i in range(0, initial_num):
        rtrb = rand_vec[i]
        for j in range(0, hood_size):
            all_opt.append(turbines[rtrb].hood[j])
            
    all_opt = sorted(all_opt, key=lambda x: x[2])   #sorted list of near-by connections
           
    networks = []
    
    
    networks.append([all_opt[0]])
    first = all_opt[0][0]
    second = all_opt[0][1]
    #print(len(all_opt))
    deletables = []
    for k in range(0, len(all_opt)):
        if first == all_opt[k][0] and second == all_opt[k][1]:
            deletables.append(k) #eliminate looping options in mesh
        if first == all_opt[k][1] and second == all_opt[k][0]:
            deletables.append(k) #eliminate looping options in mesh
    deletables.sort()
    for i in range(1, len(deletables) + 1):
            j = len(deletables) - i
            delete = deletables[j]
            del all_opt[delete]
    
        
    #print(len(all_opt))
    while len(all_opt) > 0:
        #network.append(all_opt[0])
        sets = []
        first = all_opt[0][0]
        second = all_opt[0][1]
        #Create List of all values in a mesh
        for j in range(0, len(networks)):
            m = []
            n = len(networks[j])
            for i in range(0, n):
                m.append(networks[j][i][0])
                m.append(networks[j][i][1])
            b = list(set(m))
            sets.append(b)
        
        #print('sets')
        #print(sets)
        #print('networks')
        #print(networks)
        maxcount = 0
        maxsets = []
        #maxsetscheck = []
        for j in range(0, len(networks)):
            counta = 0
            countb = 0
            if first in sets[j] or second in sets[j]:
                counta = 1
                if j not in maxsets:
                    maxsets.append(j)
            for i in range(0, len(networks)):
                if i != j:
                    if first in sets[i] or second in sets[i]:
                        countb = 1
                        if i not in maxsets:
                            maxsets.append(i)
            if (counta + countb) > maxcount:
                maxcount = counta + countb
        
        deletables = []            
        if maxcount == 0: #form a new mesh
            networks.append([all_opt[0]])
            #print('count was 0')
            for k in range(0, len(all_opt)):
                if first == all_opt[k][0] and second == all_opt[k][1]:
                    deletables.append(k) #eliminate looping options in mesh
                if first == all_opt[k][1] and second == all_opt[k][0]:
                    deletables.append(k) #eliminate looping options in mesh
            
        elif maxcount == 1: #if one of the coordinates is already in a mesh, add to mesh
            a = maxsets[0]
            networks[a].append(all_opt[0])
            #print('count was 1')
            if first not in sets[a]:
                #print('enter part 1')
                f = len(sets[a])
                for h in range(0, f):
                    third = sets[a][h]
                    for k in range(0, len(all_opt)):
                        if first == all_opt[k][0] and third == all_opt[k][1]:
                            deletables.append(k) #eliminate looping options in mesh
                        if first == all_opt[k][1] and third == all_opt[k][0]:
                            deletables.append(k) #eliminate looping options in mesh
            if second not in sets[a]:
                #print('enter part 2')
                f = len(sets[a])
                for h in range(0, f):
                    third = sets[a][h]
                    for k in range(0, len(all_opt)):
                        if second == all_opt[k][0] and third == all_opt[k][1]:
                            deletables.append(k) #eliminate looping options in mesh
                        if second == all_opt[k][1] and third == all_opt[k][0]:
                            deletables.append(k) #eliminate looping options in mesh
        
        elif maxcount == 2: #if one number is in one mesh, and other is in another mesh, combine meshes and delete second mesh
            a = maxsets[0]
            b = maxsets[1]
            c = len(networks[b])
            for each in range(0, c):    #combine meshes
                networks[a].append(networks[b][each])
            networks[a].append(all_opt[0])
            del networks[b]
    
            #print('count was 2')
            c = len(sets[a])
            d = len(sets[b])
            for each in range(0, c): #eliminate looping options between meshes
                third = sets[a][each]
                for every in range(0, d):   
                    fourth = sets[b][every]
                    for k in range(0, len(all_opt)):
                        if third == all_opt[k][0] and fourth == all_opt[k][1]:
                            deletables.append(k) #eliminate looping options in mesh
                        if third == all_opt[k][1] and fourth == all_opt[k][0]:
                            deletables.append(k) #eliminate looping options in mesh
    
        #else: #something went wrong
            #print('oh Shit')
            
        #print('deletables')
        #print(deletables)
        #print('all_opt length: ', len(all_opt))
        deletables.sort()
        for i in range(1, len(deletables) + 1):
            j = len(deletables) - i
            delete = deletables[j]
            del all_opt[delete]
        #print('all_opt length: ', len(all_opt))
        
    #print(networks[0])
    icl = 0.
    a = len(networks[0])
    for i in range(0, a):
        icl += networks[0][i][2]
    icl = icl/1000. #m to km
    #print(icl, 'm of interior cable needed')
    #cableout = open('cablein.txt', 'w')
    #cableout.write(str(str(icl) + '\n'))
    #cableout.write(str('connections\n'))
    #cableout.write(str(networks))
    return icl, networks
########################################################################################################################
def plot_cables(networks):    
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    a = len(networks[0])
    for i in range(0, a):
        j = networks[0][i][0]
        k = networks[0][i][1]
        xnew = [turbines[j].XLocation[0], turbines[k].XLocation[0]]
        ynew = [turbines[j].YLocation[0], turbines[k].YLocation[0]]
        ax1.plot(xnew, ynew, c='k')
        ax1.scatter(xnew, ynew, s=10, c='r', marker="o")
        #plt.axis([0, 2000, 0, 2000])
        plt.ylabel('Position (m)')
        plt.xlabel('Position (m)')
        plt.title(str('Cable Layout for ' + str(initial_num) + ' Turbines'))
##########################################################################################################################
hardcode = False
rand_start = True
inf_start = False

if inf_start == True:
    inf = open('input_num.txt', 'r')     #Initial Number of Turbines 
    num_lines = sum(1 for line in inf)
    #num_lines = int(raw_input('how many turbine number options? '))
    inf.close()
    inf = open('input_num.txt', 'r')
else:
    num_lines = 1
    
if output == 'on':
    excel = open('Jensen_out/Case2, BreifDatapt.txt', 'w') 
    excel.write(str('number of turbines, final evaluation, power output, efficiency, number of evaluations made, input file line number \n'))
#i = 1
initial_num = 30
for line in range(0, num_lines):
    nemo = 0
    if inf_start == True:
        g = inf.readline()
        initial_num = int(g.rstrip())
        print(initial_num)
    if output == 'on':
        #initial_num = int(initial_num_file.readline())
        trial_num = 1
        while os.path.isfile('Jensen_out/Case2 ' + str(initial_num) + 'Test ' + str(trial_num) + '.txt'):
            trial_num += 1
        data_out = open('Jensen_out/Case2 ' + str(initial_num) + 'Test ' + str(trial_num) + '.txt', 'w')
    if rand_start == True:
        XLocation = [[0.0] * len(directions)] * initial_num
        YLocation = [[0.0] * len(directions)] * initial_num
    ZLocation = [0.0] * initial_num
    HubHeight = [80.0] * initial_num
    OldHubHeight = [80.0] * initial_num
    ZHub = []
    RotorRad = [40.0] * initial_num
    if layout_type == 'colorado':
        RotorRad = [77.0 / 2.] * initial_num
    OldRotorRad = [40.0] * initial_num
    alpha = [0.0] * initial_num
    usturbinesrec = []
    usturbines = [[]] * len(directions)
    dsturbinesrec = []
    dsturbines = [[]] * len(directions)
    wakewidth = [[]] * len(directions)
    distance = [[]] * len(directions)
    percent = [[]]* len(directions)
    ui = [10.0] * len(directions)     #changes in Compute_Wake function
    windspeeds = [0.] * len(directions)
    xcoords = [[]] * len(directions)
    zcoords = [[]] * len(directions)
    Power = [0.0] * len(directions)
    Area = 0

    for i in range(0, initial_num):
        ZHub.append(HubHeight[i] + ZLocation[i])

    #initialize turbines         
    turbines = [Turbine(XLocation[i], YLocation[i], ZLocation[i], HubHeight[i], RotorRad[i], alpha[i], usturbinesrec, usturbines, dsturbinesrec, dsturbines, wakewidth, distance, percent, ui, windspeeds, xcoords, zcoords, Power, Area) for i in range(0, initial_num)]        
    if rand_start == True:
        #set turbine starting locations       
        reset = 50000.    
        while reset >= 50000.:        
            reset = Initial_Layout()
    print_graph()        
    if output == 'on':    
        #print("Initial X Locations: ", XLocation)
        #print('Initial Y Locations: ', YLocation)
        XLocation = []
        YLocation = []
        ZLocation = []
        HubHeight = []
        RotorRad = []
        for l in range(0, initial_num):
            XLocation.append(turbines[l].XLocation[0])
            YLocation.append(turbines[l].YLocation[0])
            ZLocation.append(turbines[l].ZLocation)
            HubHeight.append(turbines[l].HubHeight)
            RotorRad.append(turbines[l].RotorRad)
    
        data_out.write(str('Initial X Locations: ' + str(XLocation) + '\n'))
        data_out.write(str('Initial Y Locations: ' + str(YLocation) + '\n'))
        data_out.write(str('Initial Z Locations: ' + str(ZLocation) + '\n'))
        data_out.write(str('Initial hub heights: ' + str(HubHeight) + '\n'))
        data_out.write(str('Initial rotor radii: ' + str(RotorRad) + '\n'))

    #Evaluate Initial Layout
    score = Eval_Objective(initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, aif)       
    if output == 'on':
        data_out.write(str('The initial layout has a score of: ' + str(score) + '\n'))
    if hardcode == False:
        start_time = time()
        evals = Pattern_Search(init_step, minstep, random_vec, initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, num_pops, max_pop_tries, hstep, i, hstepmin, rradmin, rradmax, rstep, hubmin, hubmax, rstepmin, aif, directions)
        print('evaluation conducted ' + str(evals) + ' times.')
        total_time = time() - start_time
        print_graph()
        print('time to run EPS: ',total_time)
        if output == 'on':
            XLocation = []
            YLocation = []
            ZLocation = []
            HubHeight = []
            RotorRad = []
            for l in range(0, initial_num):
                XLocation.append(turbines[l].XLocation[0])
                YLocation.append(turbines[l].YLocation[0])
                ZLocation.append(turbines[l].ZLocation)
                HubHeight.append(turbines[l].HubHeight)
                RotorRad.append(turbines[l].RotorRad)
            
            #write out final locations
            data_out.write(str('Final X Locations: ' + str(XLocation) + '\n'))
            data_out.write(str('Final Y Locations: ' + str(YLocation) + '\n'))
            data_out.write(str('Final Z Locations: ' + str(ZLocation) + '\n'))
            data_out.write(str('Final hub heights: ' + str(HubHeight) + '\n'))
            data_out.write(str('Final rotor radii: ' + str(RotorRad) + '\n'))
        #Evaluate Final Layout
        score = Eval_Objective(initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, aif)
        if output == 'on':
            data_out.write(str('The final layout has a score of: ' + str(score) + '\n'))
    print('The final layout with ', initial_num, ' turbines has a score of: ', score)

    '''
    for i in turbines:
        print(i.ui[0])
    
    for i in turbines:
        print(i.Power[0])
    '''
    #print('Farm Output Power is ', total, 'Watts')
    if output == 'on':
        for i in range(0, initial_num):
            data_out.write(str('The effective windspeed for turbine ' + str(i) + ' is ' + str(turbines[i].ui) + '\n'))
            data_out.write(str('Turbine ' + str(i) + ' power: ' + str(turbines[i].Power) + '\n'))
     
        data_out.write(str('The total power generated was: ' + str(score) + ' kW. \n'))    
    
        #rorad = turbines[k].RotorRad
        #Area = np.pi * (rorad ** 2.0)
        #temp1 = 0.3 * (12 ** 3.0)
        #possible += temp1
        #efficiency = (total / possible) * 100
        #print('Farm efficiency is ', efficiency, '%.')
        data_out.write(str('The Objective was calculated ' + str(nemo) + ' times. \n'))
        data_out.write(str('The pattern serach ran in '+str(total_time)+' s\n'))
    
        #print_graph()
        #plot_cables(networks)
    
        data_out.close() 
        
#initial_num += 1
#plt.plot(XLocation, YLocation, 'r^')
#plt.axis([0, 2000, 0, 2000])
#plt.ylabel('Position (m)')
#plt.xlabel('Position (m)')
#plt.title(str('Optimization of ' + str(initial_num) + ' Turbines'))
#plt.savefig(str(str(initial_num) + 'turbinesTest.png'), bbox_inches='tight')
#plt.clf()
#plt.show()
    
#initial_num += 1
if output == 'on':
    excel.close()
if inf_start == True:
    inf.close()
#data_out.close()
#initial_num_file.close()
#print('files are closed')
#for i in range(0, initial_num):
#    print([i, turbines[i].XLocation, turbines[i].YLocation])
