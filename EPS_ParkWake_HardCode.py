# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 11:21:39 2016

@author: Annalise
"""

import math
import random
import os
import matplotlib.pyplot as plt
global nemo
import numpy as np
nemo = 1

init_step = 400.0           #Step Size for translation, in m

LengthX = 50000.               #Length of one side of solution area, in m
LengthY = 50000.  #15. * 80.         
#LengthX = 2000.0               #Length of one side of solution area, in m
#LengthY = 2000.0               #Length of one side of solution area, in m
U0 = 8.0                   #mean wind speed, in m/s
#z = 60                     #hub height, in m
#Rr = 20                    #Rotor Radius, in m
Ct = 0.88                   #Thrust Coefficient

aif = 0.1909830056250526    #Axial induction factor based on Cp = 0.5
        
output = 'off'    
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
Cutin = 3.0                 #Cut in speed from power curve
Zref = 80.0                 #Reference height for wind shear calculation, in m
depth = 200.0               #water depth, 200m. User-defined in future
WCOE = 0.1                  #wholdsale cost of electricity, in $/kWh
yrs = 20.0

#Cp = 0.5
ma = 0.33
Cp = 4*float(ma)*(1.-float(ma))**2
Cf = 1.0


shore = 'on'
#onshore/offshore specific
if shore == 'on':
    z0 = 0.04                 #Surface roughness, in m ###changed to match NREL
    
elif shore == 'off':
    z0 = 0.0005                 #Surface roughness, in m

#condition = raw_input("Input Condition NEUTRAL, STABLE, or UNSTABLE: ") #Input wind condition
condition = 'NEUTRAL'

if condition == "NEUTRAL":
    APow = 0.08835          #Neutral Conditions: WDC(h) = APow*h^BPow
    BPow = -0.1521      
    if shore == 'off':
        alphah = 0.11           #Power Law Exponent - Neutral (averaged over seasons)
    elif shore == 'on':
        alphah = 0.15567           #Power Law Exponent - Neutral (averaged over seasons)
    
elif condition == "STABLE":
    APow = 0.07535
    BPow = -0.1496
    alphah = 0.14
    
elif condition == "UNSTABLE":
    APow = 0.09759
    BPow = -0.1352
    alphah = 0.08
    
    
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
def Translation_X(step_size, i):
    transflag = 0
    xstart = turbines[i].XLocation
    #print(xstart, step_size)
    xfinish = xstart + step_size            #find preliminary new x-location given step size
    if xfinish >= 0 and xfinish <= LengthX:    #if this new x-location is not out of bounds, translate it
        XLocation = xfinish   
        turbines[i].XLocation = XLocation
        return transflag
    else:
        transflag = 1
        return transflag


################################################################################################################        
def Translation_Y(step_size, i):
    transflag = 0
    ystart = turbines[i].YLocation
    yfinish = ystart + step_size            #find preliminary new x-location given step size
    if yfinish >= 0 and yfinish <= LengthY:    #if this new x-location is not out of bounds, translate it
        YLocation = yfinish   
        turbines[i].YLocation = YLocation
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

            random_X = random.uniform(0, LengthX)
            xmove = random_X
            #print('xmove = ', xmove)
            Translation_X(xmove, n)
            #print(turbines[n].XLocation)
            random_Y = random.uniform(0, LengthY)
            ymove = random_Y
            #print('ymove = ', ymove)
            Translation_Y(ymove, n)
            #print(turbines[n].YLocation)
            
            CHECK2 = Check_Interference(n)
            
            if CHECK2 != 1:
                checkx = 1                  #If there is no interference and the turbine can be placed, then choose new corresponding z-coordinate
                
            else:
                Translation_X(-xmove, n)  #If there is interference, move turbine back to origin
                Translation_Y(-ymove, n)
                reset += 1
                
        if reset == 5000:
            for l in range(0, initial_num):
                turbines[l].XLocation = 0.0
                turbines[l].YLocation = 0.0
            return reset

    return reset

    
###################################################################################################################    
def Check_Interference(n):
    CHECK2 = 0
    checkx = 0
    x = turbines[n].XLocation
    y = turbines[n].YLocation
    for k in range(0, initial_num):
        if k != n:
            xnew = turbines[k].XLocation
            ynew = turbines[k].YLocation
            checkx = x - xnew
            checky = y - ynew
            checkrad = math.sqrt(checkx ** 2.0 + checky ** 2.0)
            if checkrad < 200:
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
#            checkrad.append(math.sqrt(checkx ** 2 + checky ** 2))
#    checkcheck.append(min(checkrad))
#print (checkcheck)  a



################################################################################################################    
def Compute_Wake(initial_num, z0, U0, Zref, alphah, ro, aif):
    turb = [n for n in range(0,initial_num)]
    XLocation = []
    YLocation = []
    for l in range(0, initial_num):
        XLocation.append(turbines[l].XLocation)
        YLocation.append(turbines[l].YLocation)
    YOrder = list(zip(turb, YLocation))               #Create tuple to identify turbine number
    YOrder = sorted(YOrder, key=lambda x: x[1])          #Sort by upstream YLocation
    #print(YOrder) #code check
    #print('yorder: ',YOrder) 
    In_Wakej = []
    In_Wakek = []
    wake_d = []
    disty = []   
    
    #Check within maximum wake width for downstream turbines
    for i in range(0, initial_num):         
        y = YOrder[i][1]                        #Go through Y locations in order
        j = YOrder[i][0]
        x = turbines[j].XLocation              #Find corresponding x location
        hubheight = turbines[j].HubHeight 
        #Calculate and set the wake decay constant for each turbine based on its current height
        #turbines[j].alpha = APow * (hubheight ** BPow)      
        alpha = 0.5 / (math.log(hubheight / z0))
        turbines[j].alpha = alpha
        dis = LengthY - y                          #Sets dis to max downstream distance

        Rr = turbines[j].RotorRad
        r1 = (alpha * dis) + Rr                 #Calculate maximum ds wake radius
        if alpha * dis < 0:
            print(j)
            print(alpha)
            print(dis)
            print(y)
        space1 = x + r1
        space2 = x - r1
        #print(j) #code check
        
        for l in range(0, initial_num):
            k = YOrder[l][0]
            #print(k) #code check
            khub = turbines[k].XLocation                 #turbine k's X Location
            krad = turbines[k].RotorRad
            kleft = khub - krad                 #define left most point of rotor swept area of turbine k
            kright = khub + krad                #define right most point of rotor swept area of turbine k
            
            if j != k and turbines[j].YLocation < turbines[k].YLocation:
                #if k >= 12:
                #    print('pre on turbine: ', k)
                if (kleft >= space2 and kleft <= space1) or (kright >= space2 and kright <= space1):
                    #if either rotor swept area points are within maximum rectangular wake
                    if turbines[j].YLocation < turbines[k].YLocation:      #and turbine is actually downstream
                        #if k >= 12:
                        #    print('post on turbine: ', k, j)
                        dsy = turbines[k].YLocation              #location of downstream turbine
                        #print(dsy) #code check
                        dist = dsy - y                  #distance between turbines
                        #print(dist) #code check
                        wake_rad = (alpha * dist) + Rr  #define radius of triangular wake                        
                        #print(j, k)     #Code Check

                        kz = turbines[k].HubHeight   #kz is the z coordinate of the rotor k's hub
                        jz = turbines[j].HubHeight   #jz is the z coordinate of the wake's center
                        
                        cd = math.sqrt((x - khub) ** 2.0 + (jz - kz) ** 2.0)   #distance between the centerline of wake and rotor hub
                        if cd < (wake_rad + krad):  #if distance between centers is less than the sum of the two radii, the rotor swept area is in the wake
                            In_Wakej.append(j)
                            In_Wakek.append(k)
                            wake_d.append(wake_rad * 2.0)
                            disty.append(dist)
                elif k == 12 and j == 8:
                    
                    print('kleft: ', kleft)
                    print('space1: ', space1)
                    print('space2: ',space2)
                    print(turbines[j].YLocation < turbines[k].YLocation)
                    
    #CodeCheck
    #print(In_Wakej)
    #print(In_Wakek)
    #print(wake_d)
    #print(disty)
    
    q = set(In_Wakej)
    unique_us = list(q)
    q = set(In_Wakek)
    unique_ds = list(q)

    
    #downstream turbines
    for w in range(0, len(unique_us)):          #each unique value in upstream turbine list
        z = unique_us[w]
        q = []
        for l in range(0, len(In_Wakej)):       #find indicies associate with each unique value
            val = In_Wakej[l]                   
            if z == val:                        
                ds = In_Wakek[l]
                q.append(ds)
        turbines[z].dsturbines = q

    #upstream turbines
    for w in range(0, len(unique_ds)):
        z = unique_ds[w]
        q = []
        b = []
        c = []
        for l in range(0, len(In_Wakek)):
            val = In_Wakek[l]
            if z == val:
                us = In_Wakej[l]
                q.append(us)
                b.append(wake_d[l])
                c.append(disty[l])
        print(q)
        turbines[z].usturbines = q
        turbines[z].wakewidth = b
        turbines[z].distance = c          


#code check           
#print(turbines[7].dsturbines)
#print(turbines[7].dsturbinesrec)
#print(turbines[7].usturbines)
#print(turbines[7].usturbinesrec)
#print(turbines[7].wakewidth)
#print(turbines[7].distance)
#print(turbines[4].dsturbines)
#print(turbines[4].dsturbinesrec)
#print(turbines[4].usturbines)
#print(turbines[4].usturbinesrec)
#print(turbines[4].wakewidth)
#print(turbines[4].distance)
                   
    #Now that we know which turbines are downstream of others, calculate the percentage of the rotor swept area that is within the wake                        
    for i in range(0, initial_num):
        overlap_flag = 0
        k = YOrder[i][0]             #k represents the first ordered turbine, then through the rest in order
        kz = turbines[k].HubHeight   #turbine k's hub height (z-location)
        kx = turbines[k].XLocation   #turbine k's x location
        #ky = turbines[k].YLocation   #turbine k's y location
        krad = turbines[k].RotorRad  #turbine k's rotor radius
        
        if len(turbines[k].usturbines) == 1:        #if the turbine has one upstream turbine
            j = turbines[k].usturbines[0]
            jz = turbines[j].HubHeight                  #z coordinate of the wake's center
            jx = turbines[j].XLocation                  #x coordinate of the wake's center
            #jy = turbines[j].YLocation                  #y coordinate of the wake's center
            jwakerad = (turbines[k].wakewidth[0]) / 2.0   #radius of wake width
            dist = turbines[k].distance[0]
            cd = math.sqrt(((jx-kx) ** 2.0) + ((jz - kz) ** 2.0))   #distance between centerline of wake and rotor hub           
            int1_den = 2.0 * cd * krad

            if cd + krad <= jwakerad:                         #if dsturbine is completely in usturbine wake, overlap = 100%
                turbines[k].percent = [1.0]
                #print('works')
            
            elif cd + jwakerad <= krad:                 #if the wake is fully encompassed by the rotor diameter    
                wakearea = math.pi * (jwakerad ** 2.0)
                percentwake = wakearea / (math.pi * (krad ** 2.0))
                turbines[k].percent = [percentwake]

            else:
                integrand1 = ((cd ** 2.0) + (krad ** 2.0) - (jwakerad ** 2.0)) / int1_den
                #print(integrand1)
                int2_den = 2.0 * cd * jwakerad                
                integrand2 = ((cd ** 2.0) + (jwakerad ** 2.0) - (krad ** 2.0)) / int2_den
                #print(integrand2) 
                q = (krad ** 2.0) * (math.acos(integrand1)) 
                b = (jwakerad ** 2.0) * (math.acos(integrand2))
                c = 0.5 * math.sqrt((-cd + krad + jwakerad) * (cd + krad - jwakerad) * (cd - krad + jwakerad) * (cd + krad + jwakerad))
                AOverlap = q + b - c
                RSA = ((math.pi) * (krad ** 2.0))
                z = AOverlap / RSA
                turbines[k].percent = [(z)]      #percentage of RSA that has wake interaction 
            
       
        elif len(turbines[k].usturbines) == 2:            #if the turbine has two upstream turbines
            first = turbines[k].usturbines[0]
            second = turbines[k].usturbines[1]
            firstx = turbines[first].XLocation
            firstz = turbines[first].HubHeight
            firstrad = turbines[k].wakewidth[0] / 2.0
            secondx = turbines[second].XLocation
            secondz = turbines[second].HubHeight
            secondrad = turbines[k].wakewidth[1] / 2.0
            cd = math.sqrt(((firstx - secondx) ** 2.0) + ((firstz - secondz) ** 2.0))   #distance between the centerline of wake and rotor hub

            if cd > (firstrad + secondrad):     #if wakes do not overlap at all within the rotor swept area
                m = []
                overlap_flag = 1
                for q in range(0, len(turbines[k].usturbines)):
                    j = turbines[k].usturbines[q]
                    jz = turbines[j].HubHeight             #z coordinate of the wake's center
                    jx = turbines[j].XLocation             #x location of the wake's center
                    #jy = turbines[j].YLocation             #y location of the wake's center
                    jwakerad = (turbines[k].wakewidth[q]) / 2.0
                    dist = turbines[k].distance[q]
                    cd = math.sqrt(((jx - kx) ** 2.0) + ((jz - kz) ** 2.0))     #distance between the centerline of wake and rotor hub
                       

                    if cd + krad <= jwakerad:
                        m.append(1.0)

                    elif cd + jwakerad <= krad:           #if the wake is fully encompassed by the rotor diameter
                        wakearea = math.pi * (jwakerad ** 2.0)
                        percentwake = wakearea / (math.pi * (krad ** 2.0))
                        m.append(percentwake)
                        
                    else:
                        integrand1 = ((cd ** 2.0) + (krad ** 2.0) - (jwakerad ** 2.0)) / (2.0 * cd * krad)
                        integrand2 = ((cd ** 2.0) + (jwakerad ** 2.0) - (krad ** 2.0)) / (2.0 * cd * jwakerad)             
                        d = (krad ** 2.0) * (math.acos(integrand1)) 
                        b = (jwakerad ** 2.0) * (math.acos(integrand2))
                        c = 0.5 * math.sqrt((-cd + krad + jwakerad) * (cd + krad - jwakerad) * (cd - krad + jwakerad) * (cd + krad + jwakerad))
                        AOverlap = d + b - c
                        #AOverlap = ((krad ** 2.0) * math.acos(integrand1)) + ((jwakerad ** 2.0) * math.acos(integrand2)) - 0.5 * math.sqrt(abs((-cd + krad + jwakerad) * (cd + krad - jwakerad) * (cd - krad + jwakerad) * (cd + krad + jwakerad)))
                        RSA = math.pi * (krad ** 2.0)
                        z = AOverlap / RSA
                        m.append(z)      #percentage of RSA that has wake interaction
                turbines[k].percent = m
            else:
                turbines[k].percent = []
                
        if len(turbines[k].usturbines) >= 2 and overlap_flag != 1:      #if there are at least 2 upstream turbines whose wakes overlap, discretize the RSA and evaluate each point
            Discretize_RSA(k)

#Code Check
#Compute_Wake(initial_num, z0, U0, Zref, alphah, ro, aif)                   
                  
    #calculate wind speed for each downstream turbine based on downstream distance
    for i in range(0, initial_num):
        k = YOrder[i][0]                        #for each turbine in increasing y-order
        hubheight = turbines[k].HubHeight
        #temp = (0.5 / math.log(hubheight / z0))
        #turbines[k].alpha = temp
        #alpha = turbines[k].alpha
        wdsp = []
        if len(turbines[k].usturbines) == 0:        #if turbine has no upstream turbines, 
               #INCORPORATE POWER LAW
            hubheight = turbines[k].HubHeight
            Uz = U0 * ((hubheight / Zref) ** alphah)      #corrects wind speed for hub height
            wdsp.append(Uz)
            turbines[k].ui = Uz
           
        elif len(turbines[k].usturbines) == 1:        #if turbine has 1 upstream turbine
            total = 0.0
            #USturb = turbines[k].usturbines[0]
            #USht = turbines[USturb].HubHeight
            x = turbines[k].distance[0]
            hubheight = turbines[k].HubHeight
            temp = (0.5 / math.log(hubheight / z0))
            #turbines[k].alpha = temp
            alpha = temp
            Rr = turbines[k].RotorRad
            
            #Grady Model 
            r1 = Rr * math.sqrt((1-aif) / (1 - 2*aif))
            EWU = U0 * (1 - (2*aif)/((1+alpha*(x/r1))**(2)))
            Uz = EWU * ((hubheight / Zref) ** alphah)
            portion = Uz * turbines[k].percent[0]
            remainder = U0 * (1.0 - turbines[k].percent[0]) * ((hubheight / Zref) ** alphah)
            total = portion + remainder                 #weighted average of windspeeds
            wdsp.append(total)
            turbines[k].ui = total
            '''
            usturb = turbines[k].usturbines[0]
            EWUT = turbines[usturb].ui
            
            #INCORPORATE POWER LAW
            Uz = EWUT * ((hubheight / Zref) ** alphah)      #corrects wind speed for height of turbine
            EWUT = Uz       #set turbine's current effective wind speed to power law-related wind speed
            ueff = EWUT * (1.0 - (2.0 / 3.0) * (Rr / (Rr + alpha * x)) ** 2.0)
            portion = ueff * turbines[k].percent[0]
            Uz = U0 * ((hubheight / Zref) ** alphah)
            remainder = Uz * (1.0 - turbines[k].percent[0])
            total = portion + remainder                 #weighted average of windspeeds
            wdsp.append(total)
            turbines[k].ui = total
            '''           
            
        elif len(turbines[k].usturbines) == 2 and len(turbines[k].percent) != 0:      #if the turbine has two upstream turbines whose wakes do not overlap
            portion = 0.0
            total = 0.0
            for j in range(0, len(turbines[k].usturbines)):
                x = turbines[k].distance[j]
                #USturb = turbines[k].usturbines[j]
                hubheight = turbines[k].HubHeight
                alpha = 0.5 / math.log(hubheight / z0)
                Rr = turbines[k].RotorRad
                r1 = Rr * math.sqrt((1 - aif) / (1 - 2 * aif))
                EWU = U0 * (1 - (2*aif)/((1+alpha*(x/r1))**(2)))
                Uz = EWU * ((hubheight / Zref) ** alphah)
                portion += Uz * turbines[k].percent[j]
            remainder = U0 * (1.0 - turbines[k].percent[0] - turbines[k].percent[1]) * ((hubheight / Zref) ** alphah)
            total = portion + remainder                 #weighted average of windspeeds
            wdsp.append(total)
            turbines[k].ui = total
                
            '''
                usturb = turbines[k].usturbines[j]
                EWUT = turbines[usturb].ui
                #INCORPORATE POWER LAW
                Uz = EWUT * ((hubheight / Zref) ** alphah)
                EWUT = Uz       #set turbine's current effective wind speed to power law-related wind speed
                ueff = EWUT * (1.0 - (2.0/3.0) * ((Rr / (Rr + alpha * x)) ** 2.0))
                portion = portion + ueff * turbines[k].percent[j]
            Uz = U0 * ((hubheight / Zref) ** alphah)
            remainder = Uz * (1.0 - (turbines[k].percent[0] - turbines[k].percent[1]))
            total = portion + remainder
            wdsp.append(total)
            turbines[k].ui = total
            '''
        elif len(turbines[k].usturbines) >= 2 and len(turbines[k].percent) == 0:      #turbine has at least two upstream turbines whos wakes overlap
            coordWS = []
            for i in range(0, len(turbines[k].xcoords)):        #xcoords created in Discretize_RSA
                decWS = []
                xval = turbines[k].xcoords[i]
                zval = turbines[k].zcoords[i]
                khub = turbines[k].HubHeight
                #alpha = 0.5 / math.log(zval / z0)
                Rr = turbines[k].RotorRad
                #r1 = Rr * math.sqrt((1.0 - aif) / (1.0 - 2.0 * aif))
                for j in range(0, len(turbines[k].usturbines)):
                    x = turbines[k].distance[j]
                    
                    US = turbines[k].usturbines[j]
                    r2 = turbines[k].wakewidth[j] / 2.0
                    xc = turbines[US].XLocation         #'c' for centerline of wake
                    #yc = turbines[US].YLocation
                    zhubc = turbines[US].HubHeight
                    xturb = xval * 1.
                    #yturb = turbines[k].YLocation
                    zhubturb = zval * 1.
                    rt2 = abs(zhubturb - zhubc)        #height of the triangular portion of the chord area in z
                    rt1 = abs(xturb - xc)              #height of the triangluar portion of the chord area in x
                    space = math.sqrt((rt2 ** 2) + (rt1 ** 2))      #distance between wake center and discritized point
                 
                    if space <= r2:        #if point is within wake
                        Rr = turbines[k].RotorRad
                        alpha = 0.5 / math.log(zval / z0)
                        #Grady's a
                        r1 = Rr * math.sqrt((1 - aif) / (1 - 2 * aif))
                        Uz = U0 * (1 - (2*aif)/((1 + alpha*(x/r1))**(2)))
                        decWS.append(Uz)
                        
                        
                coordui = 0.0        
                if len(decWS) != 0:

                    if len(decWS) == 1:         #if the point only has one wake acting on it
                        coordui = (decWS[0]) * ((zval / Zref) ** alphah)
                        coordWS.append(coordui)
                        
                    elif len(decWS) > 1:          #if the pint has more than one wake acting on it
                        tally = 0.0
                        for l in range(0, len(decWS)):
                            u = decWS[l]
                            tally += ((1.0 - (u / U0)) ** 2.0)
                            
                        '''
                        if len(decWS) >= 3 and i == 48:
                            print(len(decWS))
                            print('turbine: ',k, ', pt: ', i)
                            print(tally)
                        '''
                        coordui = (U0 - (U0 * math.sqrt(tally))) * ((zval / Zref) ** alphah)
                        coordWS.append(coordui)
                        
                else:               #if the point has no wakes acting on it
                    Uz = U0 * ((zval / Zref) ** alphah)
                    coordui = Uz
                    coordWS.append(coordui)
                    '''
                        EWUT = turbines[US].ui
                        #INCORPORATE POWER LAW
                        Uz = EWUT * ((zval / Zref) ** alphah)
                        EWUT = Uz       #set turbine's current effective wind speed to power law-related wind speed
                        ueff = EWUT * (1.0 - (2.0 / 3.0) * ((Rr/ (Rr + alpha * x)) ** 2.0))
                        decWS.append(ueff)
                      
                coordui = 0.0        
                if len(decWS) != 0:
                    Uz = U0 * ((zval / Zref) ** alphah)
                    if len(decWS) == 1:         #if the point only has one wake acting on it
                        coordui = decWS[0]
                        coordWS.append(coordui)
                        
                    elif len(decWS) > 1:          #if the pint has more than one wake acting on it
                        tally = 0.0
                        for l in range(0, len(decWS)):
                            u = decWS[l]
                            tally += ((1.0 - (u / Uz)) ** 2.0)
                            
                        coordui = Uz - (Uz * math.sqrt(tally))
                        coordWS.append(coordui)
                        
                else:               #if the point has no wakes acting on it
                #INCORPORATE POWER LAW
                    hubheight = turbines[k].HubHeight
                    Uz = U0 * ((zval / Zref) ** alphah)
                    coordui = Uz
                    coordWS.append(coordui)
                '''
            #Sum discretized wind speeds
            tally2 = 0.0
            percentage = 1.0 / 49.0
            for f in range(0, len(coordWS)):
                tally2 += percentage * coordWS[f]

            d = len(coordWS)
            wdsp.append(tally2)
            turbines[k].ui = tally2
        turbines[k].windspeeds = wdsp
            
    #calculate power developed for each turbine
    for i in range(0, initial_num):
        temp1 = 0.0
        rorad = turbines[i].RotorRad
        Area = (rorad ** 2.0) * math.pi
        '''
        temp1 += 0.5 * ro * Area * (turbines[i].ui ** 3.0) * 0.5
        turbines[i].Power = temp1
        '''
        #incorporating power curve suggested by Pat, June 10th
        if turbines[i].ui < Uref and turbines[i].ui >= Cutin: #Calculate power for effective windspeeds between 3 and 11 m/s
            temp1 += 0.5 * ro * Area * (turbines[i].ui ** 3.0) * Cp * Cf / 1000.
            turbines[i].Power = temp1

        if turbines[i].ui < Cutin:        #cut in speed is 3 m/s
            turbines[i].Power = 0

        if turbines[i].ui >= Uref:      #constant for 11.5 m/s rated power and above
            temp1 += 0.5 * ro * Area * (Uref ** 3.0) * Cp * Cf / 1000.
            turbines[i].Power = temp1

    return YOrder
################################################################################################################    
def Compute_Wake_mod(initial_num, z0, U0, Zref, alphah, ro, aif):
    turb = [n for n in range(0,initial_num)]
    XLocation = []
    YLocation = []
    for l in range(0, initial_num):
        XLocation.append(turbines[l].XLocation)
        YLocation.append(turbines[l].YLocation)
    YOrder = list(zip(turb, YLocation))               #Create tuple to identify turbine number
    YOrder = sorted(YOrder, key=lambda x: x[1])          #Sort by upstream YLocation
    #print(YOrder) #code check
    #print('yorder: ',YOrder) 
    In_Wakej = []
    In_Wakek = []
    wake_d = []
    disty = []   
    
    #Check within maximum wake width for downstream turbines
    for i in range(0, initial_num):         
        y = YOrder[i][1]                        #Go through Y locations in order
        j = YOrder[i][0]
        x = turbines[j].XLocation              #Find corresponding x location
        hubheight = turbines[j].HubHeight 
        #Calculate and set the wake decay constant for each turbine based on its current height
        #turbines[j].alpha = APow * (hubheight ** BPow)      
        alpha = 0.5 / (math.log(hubheight / z0))
        turbines[j].alpha = alpha
        dis = LengthY - y                          #Sets dis to max downstream distance

        Rr = turbines[j].RotorRad
        r1 = (alpha * dis) + Rr                 #Calculate maximum ds wake radius
        if alpha * dis < 0:
            print('error: turbine beyond farm space')
            print(j)
            print(alpha)
            print(dis)
            print(y)
        space1 = x + r1
        space2 = x - r1
        #print(j) #code check
        
        for l in range(0, initial_num):
            k = YOrder[l][0]
            #print(k) #code check
            khub = turbines[k].XLocation                 #turbine k's X Location
            krad = turbines[k].RotorRad
            kleft = khub - krad                 #define left most point of rotor swept area of turbine k
            kright = khub + krad                #define right most point of rotor swept area of turbine k
            
            if j != k and turbines[j].YLocation < turbines[k].YLocation:
                #if k >= 12:
                #    print('pre on turbine: ', k)
                if (kleft >= space2 and kleft <= space1) or (kright >= space2 and kright <= space1):
                    #if either rotor swept area points are within maximum rectangular wake
                    if turbines[j].YLocation < turbines[k].YLocation:      #and turbine is actually downstream
                        #if k >= 12:
                        #    print('post on turbine: ', k, j)
                        dsy = turbines[k].YLocation              #location of downstream turbine
                        #print(dsy) #code check
                        dist = dsy - y                  #distance between turbines
                        #print(dist) #code check
                        wake_rad = (alpha * dist) + Rr  #define radius of triangular wake                        
                        #print(j, k)     #Code Check

                        kz = turbines[k].HubHeight   #kz is the z coordinate of the rotor k's hub
                        jz = turbines[j].HubHeight   #jz is the z coordinate of the wake's center
                        
                        cd = math.sqrt((x - khub) ** 2.0 + (jz - kz) ** 2.0)   #distance between the centerline of wake and rotor hub
                        if cd < (wake_rad + krad):  #if distance between centers is less than the sum of the two radii, the rotor swept area is in the wake
                            In_Wakej.append(j)
                            In_Wakek.append(k)
                            wake_d.append(wake_rad * 2.0)
                            disty.append(dist)
                '''
                elif k == 12 and j == 8:
                    
                    print('kleft: ', kleft)
                    print('space1: ', space1)
                    print('space2: ',space2)
                    print(turbines[j].YLocation < turbines[k].YLocation)
                '''    
    #CodeCheck
    #print(In_Wakej)
    #print(In_Wakek)
    #print(wake_d)
    #print(disty)
    
    q = set(In_Wakej)
    unique_us = list(q)
    q = set(In_Wakek)
    unique_ds = list(q)

    
    #downstream turbines
    for w in range(0, len(unique_us)):          #each unique value in upstream turbine list
        z = unique_us[w]
        q = []
        for l in range(0, len(In_Wakej)):       #find indicies associate with each unique value
            val = In_Wakej[l]                   
            if z == val:                        
                ds = In_Wakek[l]
                q.append(ds)
        turbines[z].dsturbines = q

    #upstream turbines
    for w in range(0, len(unique_ds)):
        z = unique_ds[w]
        q = []
        b = []
        c = []
        for l in range(0, len(In_Wakek)):
            val = In_Wakek[l]
            if z == val:
                us = In_Wakej[l]
                q.append(us)
                b.append(wake_d[l])
                c.append(disty[l])
        print(q)
        turbines[z].usturbines = q
        turbines[z].wakewidth = b
        turbines[z].distance = c          


#code check           
#print(turbines[7].dsturbines)
#print(turbines[7].dsturbinesrec)
#print(turbines[7].usturbines)
#print(turbines[7].usturbinesrec)
#print(turbines[7].wakewidth)
#print(turbines[7].distance)
#print(turbines[4].dsturbines)
#print(turbines[4].dsturbinesrec)
#print(turbines[4].usturbines)
#print(turbines[4].usturbinesrec)
#print(turbines[4].wakewidth)
#print(turbines[4].distance)
                   
    #Now that we know which turbines are downstream of others, calculate the percentage of the rotor swept area that is within the wake                        
    for i in range(0, initial_num):
        overlap_flag = 0
        k = YOrder[i][0]             #k represents the first ordered turbine, then through the rest in order
        kz = turbines[k].HubHeight   #turbine k's hub height (z-location)
        kx = turbines[k].XLocation   #turbine k's x location
        #ky = turbines[k].YLocation   #turbine k's y location
        krad = turbines[k].RotorRad  #turbine k's rotor radius
        
        if len(turbines[k].usturbines) == 1:        #if the turbine has one upstream turbine
            j = turbines[k].usturbines[0]
            jz = turbines[j].HubHeight                  #z coordinate of the wake's center
            jx = turbines[j].XLocation                  #x coordinate of the wake's center
            #jy = turbines[j].YLocation                  #y coordinate of the wake's center
            jwakerad = (turbines[k].wakewidth[0]) / 2.0   #radius of wake width
            dist = turbines[k].distance[0]
            cd = math.sqrt(((jx-kx) ** 2.0) + ((jz - kz) ** 2.0))   #distance between centerline of wake and rotor hub           
            int1_den = 2.0 * cd * krad

            if cd + krad <= jwakerad:                         #if dsturbine is completely in usturbine wake, overlap = 100%
                turbines[k].percent = [1.0]
                #print('works')
            
            elif cd + jwakerad <= krad:                 #if the wake is fully encompassed by the rotor diameter    
                wakearea = math.pi * (jwakerad ** 2.0)
                percentwake = wakearea / (math.pi * (krad ** 2.0))
                turbines[k].percent = [percentwake]

            else:
                integrand1 = ((cd ** 2.0) + (krad ** 2.0) - (jwakerad ** 2.0)) / int1_den
                #print(integrand1)
                int2_den = 2.0 * cd * jwakerad                
                integrand2 = ((cd ** 2.0) + (jwakerad ** 2.0) - (krad ** 2.0)) / int2_den
                #print(integrand2) 
                q = (krad ** 2.0) * (math.acos(integrand1)) 
                b = (jwakerad ** 2.0) * (math.acos(integrand2))
                c = 0.5 * math.sqrt((-cd + krad + jwakerad) * (cd + krad - jwakerad) * (cd - krad + jwakerad) * (cd + krad + jwakerad))
                AOverlap = q + b - c
                RSA = ((math.pi) * (krad ** 2.0))
                z = AOverlap / RSA
                turbines[k].percent = [(z)]      #percentage of RSA that has wake interaction 
            
       
        elif len(turbines[k].usturbines) == 2:            #if the turbine has two upstream turbines
            first = turbines[k].usturbines[0]
            second = turbines[k].usturbines[1]
            firstx = turbines[first].XLocation
            firstz = turbines[first].HubHeight
            firstrad = turbines[k].wakewidth[0] / 2.0
            secondx = turbines[second].XLocation
            secondz = turbines[second].HubHeight
            secondrad = turbines[k].wakewidth[1] / 2.0
            cd = math.sqrt(((firstx - secondx) ** 2.0) + ((firstz - secondz) ** 2.0))   #distance between the centerline of wake and rotor hub

            if cd > (firstrad + secondrad):     #if wakes do not overlap at all within the rotor swept area
                m = []
                overlap_flag = 1
                for q in range(0, len(turbines[k].usturbines)):
                    j = turbines[k].usturbines[q]
                    jz = turbines[j].HubHeight             #z coordinate of the wake's center
                    jx = turbines[j].XLocation             #x location of the wake's center
                    #jy = turbines[j].YLocation             #y location of the wake's center
                    jwakerad = (turbines[k].wakewidth[q]) / 2.0
                    dist = turbines[k].distance[q]
                    cd = math.sqrt(((jx - kx) ** 2.0) + ((jz - kz) ** 2.0))     #distance between the centerline of wake and rotor hub
                       

                    if cd + krad <= jwakerad:
                        m.append(1.0)

                    elif cd + jwakerad <= krad:           #if the wake is fully encompassed by the rotor diameter
                        wakearea = math.pi * (jwakerad ** 2.0)
                        percentwake = wakearea / (math.pi * (krad ** 2.0))
                        m.append(percentwake)
                        
                    else:
                        integrand1 = ((cd ** 2.0) + (krad ** 2.0) - (jwakerad ** 2.0)) / (2.0 * cd * krad)
                        integrand2 = ((cd ** 2.0) + (jwakerad ** 2.0) - (krad ** 2.0)) / (2.0 * cd * jwakerad)             
                        d = (krad ** 2.0) * (math.acos(integrand1)) 
                        b = (jwakerad ** 2.0) * (math.acos(integrand2))
                        c = 0.5 * math.sqrt((-cd + krad + jwakerad) * (cd + krad - jwakerad) * (cd - krad + jwakerad) * (cd + krad + jwakerad))
                        AOverlap = d + b - c
                        #AOverlap = ((krad ** 2.0) * math.acos(integrand1)) + ((jwakerad ** 2.0) * math.acos(integrand2)) - 0.5 * math.sqrt(abs((-cd + krad + jwakerad) * (cd + krad - jwakerad) * (cd - krad + jwakerad) * (cd + krad + jwakerad)))
                        RSA = math.pi * (krad ** 2.0)
                        z = AOverlap / RSA
                        m.append(z)      #percentage of RSA that has wake interaction
                turbines[k].percent = m
            else:
                turbines[k].percent = []
                
        if len(turbines[k].usturbines) >= 2 and overlap_flag != 1:      #if there are at least 2 upstream turbines whose wakes overlap, discretize the RSA and evaluate each point
            Discretize_RSA(k)

#Code Check
#Compute_Wake(initial_num, z0, U0, Zref, alphah, ro, aif)                   
                  
    #calculate wind speed for each downstream turbine based on downstream distance
    for i in range(0, initial_num):
        k = YOrder[i][0]                        #for each turbine in increasing y-order
        hubheight = turbines[k].HubHeight
        #temp = (0.5 / math.log(hubheight / z0))
        #turbines[k].alpha = temp
        #alpha = turbines[k].alpha
        wdsp = []
        if len(turbines[k].usturbines) == 0:        #if turbine has no upstream turbines, 
               #INCORPORATE POWER LAW
            hubheight = turbines[k].HubHeight
            Uz = U0 * ((hubheight / Zref) ** alphah)      #corrects wind speed for hub height
            wdsp.append(Uz)
            turbines[k].ui = Uz
           
        elif len(turbines[k].usturbines) == 1:        #if turbine has 1 upstream turbine
            total = 0.0
            #USturb = turbines[k].usturbines[0]
            #USht = turbines[USturb].HubHeight
            x = turbines[k].distance[0]
            hubheight = turbines[k].HubHeight
            temp = (0.5 / math.log(hubheight / z0))
            #turbines[k].alpha = temp
            alpha = temp
            Rr = turbines[k].RotorRad
            
            #Grady Model 
            r1 = Rr * math.sqrt((1-aif) / (1 - 2*aif))
            EWU = U0 * (1 - (2*aif)/((1+alpha*(x/r1))**(2)))
            Uz = EWU * ((hubheight / Zref) ** alphah)
            portion = Uz * turbines[k].percent[0]
            remainder = U0 * (1.0 - turbines[k].percent[0]) * ((hubheight / Zref) ** alphah)
            total = portion + remainder                 #weighted average of windspeeds
            wdsp.append(total)
            turbines[k].ui = total
            '''
            usturb = turbines[k].usturbines[0]
            EWUT = turbines[usturb].ui
            
            #INCORPORATE POWER LAW
            Uz = EWUT * ((hubheight / Zref) ** alphah)      #corrects wind speed for height of turbine
            EWUT = Uz       #set turbine's current effective wind speed to power law-related wind speed
            ueff = EWUT * (1.0 - (2.0 / 3.0) * (Rr / (Rr + alpha * x)) ** 2.0)
            portion = ueff * turbines[k].percent[0]
            Uz = U0 * ((hubheight / Zref) ** alphah)
            remainder = Uz * (1.0 - turbines[k].percent[0])
            total = portion + remainder                 #weighted average of windspeeds
            wdsp.append(total)
            turbines[k].ui = total
            '''           
            
        elif len(turbines[k].usturbines) == 2 and len(turbines[k].percent) != 0:      #if the turbine has two upstream turbines whose wakes do not overlap
            portion = 0.0
            total = 0.0
            for j in range(0, len(turbines[k].usturbines)):
                x = turbines[k].distance[j]
                #USturb = turbines[k].usturbines[j]
                hubheight = turbines[k].HubHeight
                alpha = 0.5 / math.log(hubheight / z0)
                Rr = turbines[k].RotorRad
                r1 = Rr * math.sqrt((1 - aif) / (1 - 2 * aif))
                EWU = U0 * (1 - (2*aif)/((1+alpha*(x/r1))**(2)))
                Uz = EWU * ((hubheight / Zref) ** alphah)
                portion += Uz * turbines[k].percent[j]
            remainder = U0 * (1.0 - turbines[k].percent[0] - turbines[k].percent[1]) * ((hubheight / Zref) ** alphah)
            total = portion + remainder                 #weighted average of windspeeds
            wdsp.append(total)
            turbines[k].ui = total
                
            '''
                usturb = turbines[k].usturbines[j]
                EWUT = turbines[usturb].ui
                #INCORPORATE POWER LAW
                Uz = EWUT * ((hubheight / Zref) ** alphah)
                EWUT = Uz       #set turbine's current effective wind speed to power law-related wind speed
                ueff = EWUT * (1.0 - (2.0/3.0) * ((Rr / (Rr + alpha * x)) ** 2.0))
                portion = portion + ueff * turbines[k].percent[j]
            Uz = U0 * ((hubheight / Zref) ** alphah)
            remainder = Uz * (1.0 - (turbines[k].percent[0] - turbines[k].percent[1]))
            total = portion + remainder
            wdsp.append(total)
            turbines[k].ui = total
            '''
        elif len(turbines[k].usturbines) >= 2 and len(turbines[k].percent) == 0:      #turbine has at least two upstream turbines whos wakes overlap
            coordWS = []
            for i in range(0, len(turbines[k].xcoords)):        #xcoords created in Discretize_RSA
                decWS = []
                xval = turbines[k].xcoords[i]
                zval = turbines[k].zcoords[i]
                khub = turbines[k].HubHeight
                #alpha = 0.5 / math.log(zval / z0)
                Rr = turbines[k].RotorRad
                #r1 = Rr * math.sqrt((1.0 - aif) / (1.0 - 2.0 * aif))
                for j in range(0, len(turbines[k].usturbines)):
                    x = turbines[k].distance[j]
                    
                    US = turbines[k].usturbines[j]
                    r2 = turbines[k].wakewidth[j] / 2.0
                    xc = turbines[US].XLocation         #'c' for centerline of wake
                    #yc = turbines[US].YLocation
                    zhubc = turbines[US].HubHeight
                    xturb = xval * 1.
                    #yturb = turbines[k].YLocation
                    zhubturb = zval * 1.
                    rt2 = abs(zhubturb - zhubc)        #height of the triangular portion of the chord area in z
                    rt1 = abs(xturb - xc)              #height of the triangluar portion of the chord area in x
                    space = math.sqrt((rt2 ** 2) + (rt1 ** 2))      #distance between wake center and discritized point
                 
                    if space <= r2:        #if point is within wake
                        check = 0 ###check to see if ws deduction has been acounted for
                        for each in turbines[US].dsturbines: ### go through each turbine in the downstream of the upstream turbine
                            if each not in turbines[k].usturbines: ### if there is no turbine between these two, go ahead and cound it
                                check += 0
                            '''WORK IN PROGRESS
                            elif 
                                ### if there is a turbine between them, but we're at a point where wakes don't overlap
                                check += 0
                            else: ### double counting the wakes of us turbines
                                check += 2
                            '''
                        
                        Rr = turbines[k].RotorRad
                        alpha = 0.5 / math.log(zval / z0)
                        #Grady's a
                        r1 = Rr * math.sqrt((1 - aif) / (1 - 2 * aif))
                        Uz = U0 * (1 - (2*aif)/((1 + alpha*(x/r1))**(2)))
                        decWS.append(Uz)
                        
                        
                coordui = 0.0        
                if len(decWS) != 0:

                    if len(decWS) == 1:         #if the point only has one wake acting on it
                        coordui = (decWS[0]) * ((zval / Zref) ** alphah)
                        coordWS.append(coordui)
                        
                    elif len(decWS) > 1:          #if the pint has more than one wake acting on it
                        tally = 0.0
                        for l in range(0, len(decWS)):
                            u = decWS[l]
                            tally += ((1.0 - (u / U0)) ** 2.0)
                            
                        '''
                        if len(decWS) >= 3 and i == 48:
                            print(len(decWS))
                            print('turbine: ',k, ', pt: ', i)
                            print(tally)
                        '''
                        coordui = (U0 - (U0 * math.sqrt(tally))) * ((zval / Zref) ** alphah)
                        coordWS.append(coordui)
                        
                else:               #if the point has no wakes acting on it
                    Uz = U0 * ((zval / Zref) ** alphah)
                    coordui = Uz
                    coordWS.append(coordui)
                    '''
                        EWUT = turbines[US].ui
                        #INCORPORATE POWER LAW
                        Uz = EWUT * ((zval / Zref) ** alphah)
                        EWUT = Uz       #set turbine's current effective wind speed to power law-related wind speed
                        ueff = EWUT * (1.0 - (2.0 / 3.0) * ((Rr/ (Rr + alpha * x)) ** 2.0))
                        decWS.append(ueff)
                      
                coordui = 0.0        
                if len(decWS) != 0:
                    Uz = U0 * ((zval / Zref) ** alphah)
                    if len(decWS) == 1:         #if the point only has one wake acting on it
                        coordui = decWS[0]
                        coordWS.append(coordui)
                        
                    elif len(decWS) > 1:          #if the pint has more than one wake acting on it
                        tally = 0.0
                        for l in range(0, len(decWS)):
                            u = decWS[l]
                            tally += ((1.0 - (u / Uz)) ** 2.0)
                            
                        coordui = Uz - (Uz * math.sqrt(tally))
                        coordWS.append(coordui)
                        
                else:               #if the point has no wakes acting on it
                #INCORPORATE POWER LAW
                    hubheight = turbines[k].HubHeight
                    Uz = U0 * ((zval / Zref) ** alphah)
                    coordui = Uz
                    coordWS.append(coordui)
                '''
            #Sum discretized wind speeds
            tally2 = 0.0
            percentage = 1.0 / 49.0
            for f in range(0, len(coordWS)):
                tally2 += percentage * coordWS[f]

            d = len(coordWS)
            wdsp.append(tally2)
            turbines[k].ui = tally2
        turbines[k].windspeeds = wdsp
            
    #calculate power developed for each turbine
    for i in range(0, initial_num):
        temp1 = 0.0
        rorad = turbines[i].RotorRad
        Area = (rorad ** 2.0) * math.pi
        '''
        temp1 += 0.5 * ro * Area * (turbines[i].ui ** 3.0) * 0.5
        turbines[i].Power = temp1
        '''
        #incorporating power curve suggested by Pat, June 10th
        if turbines[i].ui < Uref and turbines[i].ui >= Cutin: #Calculate power for effective windspeeds between 3 and 11 m/s
            temp1 += 0.5 * ro * Area * (turbines[i].ui ** 3.0) * Cp * Cf / 1000.
            turbines[i].Power = temp1

        if turbines[i].ui < Cutin:        #cut in speed is 3 m/s
            turbines[i].Power = 0

        if turbines[i].ui >= Uref:      #constant for 11.5 m/s rated power and above
            temp1 += 0.5 * ro * Area * (Uref ** 3.0) * Cp * Cf / 1000.
            turbines[i].Power = temp1

    return YOrder
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
        Area = (rorad ** 2.0) * math.pi
        turbines[i].Area = Area             #define Area for each turbine
        #cf = 0.4    #capcity factor
        Prated = 0.5 * ro * Area * (Uref ** 3.0) * Cp/1000000      #calculate rated power in Watts
        turbines[i].Prated = Prated
        Capital_Cost = Prated * 1480000
        turbines[i].Capital_Cost = Capital_Cost
        O_M = (Prated * 133000.0) * yrs
        turbines[i].O_M = O_M
        #Substation_Cost = 0.26284 * Prated     #From Myhr
        #Substation_Cost = 0.662269 * Prated    #From JEDI
        Substation_Cost = 20000 * Prated
        turbines[i].Substation_Cost = Substation_Cost
        Leasing_Cost = (Prated) * 8760.0 * Cf * WCOE * (8.0 * 0.02 + (yrs - 8.0) * 0.04)
        turbines[i].Leasing_Cost = Leasing_Cost     #(8.0 * 0.02 + (yrs - 8.0) * 0.04) is r, r = 0.02 for first 8 years, 0.04 for rest of life of project
        Mooring_Cost = 4.0 * (140148.0 + 274.0 * depth)
        turbines[i].Mooring_Cost = Mooring_Cost
        #Installation_Cost = turbines[i].Installation_Cost
        turbines[i].Cost = turbines[i].Capital_Cost + turbines[i].O_M + turbines[i].Substation_Cost + turbines[i].Leasing_Cost + turbines[i].Mooring_Cost + turbines[i].Installation_Cost   #cost related to Prated
        cost += turbines[i].Cost


    d_t, networks = calcicl(initial_num)
    d_s = 32.0 #assume 32 km offshore (20 mi)
    Cabling_Cost = d_t * 307000 + d_s * 492000
    cost += Cabling_Cost    
    cost += 2000000.0 #66210000.0 or 2000000 from substation cost

        
    return cost
 
#Compute_Cost(initial_num, ro, yrs, WCOE, condition)   
 
################################################################################################################
def Eval_Objective(initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, aif):
    Compute_Wake(initial_num, z0, U0, Zref, alphah, ro, aif)
    cost = Compute_Cost(initial_num, ro, yrs, WCOE, condition, depth)
    #cost = 0.0
    global nemo
    nemo += 1
    Ptot = 0.0
    windspeeds = 0.0
    for i in range(0, initial_num):
        Ptot = Ptot + turbines[i].Power
        windspeeds += turbines[i].ui

    objective = cost / Ptot
    
    return objective
    
    
#####################################################################################################################
def Clear_Vectors():
    for i in range(0, initial_num):
        turbines[i].percent = []
        turbines[i].distance = []
        turbines[i].dsturbines = []
        turbines[i].wakewidth = []
        turbines[i].usturbines = []
        turbines[i].dsturbinesrec = []
        turbines[i].usturbinesrec = []
        turbines[i].windspeeds = []



#######################################################################################################################
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


######################################################################################################################

def Pattern_Search(init_step, minstep, random_vec, initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, num_pops, max_pop_tries, hstep, i, hstepmin, rradmin, rradmax, rstep, hubmin, hubmax, rstepmin, aif):
    
    Clear_Vectors()
        
    for h in range(0, 1):
        #objective = Eval_Objective(initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, aif)

        #print('objective eval: ', objective)
        #print('h= ', h)
        step2 = init_step
        while step2 >= minstep:
            random_vec = Rand_Vector(initial_num)    #creates a randomly ordered vector of turbines
            for j in range(0, len(random_vec)):
                i = random_vec[j]
                turbines[i].Stopped = 0
                #print('Turbine ', i,' is being tested.')
                flag = 0
                innerflag = 0
                transflag = 0
                nomove = Eval_Objective(initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, aif)

                #develop preliminary objective for comparison purposes
                
                Clear_Vectors()
                #print('The nomove value for turbine ', i, ' is ', nomove)
                
                #('stepped into while loop')
                if turbines[i].YLocation >= (LengthY / 2.0):
                    
                    if innerflag == 0 and flag == 0:        #move 1 was just unsucessfully attempted
                        transflag = Translation_Y(step2, i)
                        CHECK2 = 0
                        if transflag == 1:      #if the translation moved the turbine out of bounds, go to next translation
                            innerflag = 1       #move2 was attempted
                            #print('turbine not moved up.')
                        else:
                            CHECK2 = Check_Interference(i)
                            if CHECK2 == 1:         #if interference occurs, move the turbine back, go to next translation
                                Translation_Y(-step2, i)
                                innerflag = 1
                                #print('turbine not moved up.')
                            else:       #if there is no interference, evaluate and store
                                move2 = Eval_Objective(initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, aif)

                                Clear_Vectors()
                                if move2 >= nomove:       #if evaluation is worse than initial, move back, go to next translation
                                    Translation_Y(-step2, i)
                                    innerflag = 1
                                    #print('turbine not moved up.')
                                else:       #evaluation is better, keep move, go to next turbine
                                    flag = 1
                                    #print('turbine ', i, ' moved up.', move2)
                                    HubHeight_Search(hstep, i, hstepmin, XLocation, YLocation, initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, rradmin, rradmax, rstep, hubmin, hubmax, rstepmin, aif)
                                
                    if innerflag == 1 and flag == 0:        #move 2 was just unsucessfully attempted
                        transflag = Translation_X(-step2, i)
                        CHECK2 = 0
                        if transflag == 1:      #if the translation moved the turbine out of bounds, go to next translation
                            innerflag = 2       #move3 was attempted
                            #print('turbine not left.')
                        else:
                            CHECK2 = Check_Interference(i)
                            if CHECK2 == 1:         #if interference occurs, move the turbine back, go to next translation
                                Translation_X(step2, i)
                                innerflag = 2
                                #print('turbine not left.')
                            else:       #if there is no interference, evaluate and store
                                move3 = Eval_Objective(initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, aif)

                                Clear_Vectors()
                                if move3 >= nomove:       #if evaluation is worse than initial, move back, go to next translation
                                    Translation_X(step2, i)
                                    innerflag = 2
                                    #print('turbine not moved left.')
                                else:       #evaluation is better, keep move, go to next turbine
                                    flag = 1
                                    #print('turbine ', i, ' moved left.', move3)
                                    HubHeight_Search(hstep, i, hstepmin, XLocation, YLocation, initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, rradmin, rradmax, rstep, hubmin, hubmax, rstepmin, aif)
                                
                    if innerflag == 2 and flag == 0:        #move 3 was just unsucessfully attempted
                        transflag = Translation_Y(-step2, i)
                        CHECK2 = 0
                        if transflag == 1:      #if the translation moved the turbine out of bounds, go to next translation
                            innerflag = 3       #move3 was attempted
                            #print('turbine not moved down.')
                        else:
                            CHECK2 = Check_Interference(i)
                            if CHECK2 == 1:         #if interference occurs, move the turbine back, go to next translation
                                Translation_Y(step2, i)
                                innerflag = 3
                                #print('turbine not moved down.')
                            else:       #if there is no interference, evaluate and store
                                move4 = Eval_Objective(initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, aif)

                                Clear_Vectors()
                                if move4 >= nomove:       #if evaluation is worse than initial, move back, go to next translation
                                    Translation_Y(step2, i)
                                    innerflag = 3
                                    #print('turbine not moved down.')
                                else:       #evaluation is better, keep move, go to next turbine
                                    flag = 1
                                    #print('turbine ', i, ' moved down.', move4)
                                    HubHeight_Search(hstep, i, hstepmin, XLocation, YLocation, initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, rradmin, rradmax, rstep, hubmin, hubmax, rstepmin, aif)
                    if innerflag == 3 and flag == 0:
                        transflag = Translation_X(step2, i)     #move the turbine one step right
                        CHECK2 = 0
                        if transflag == 1:          #if the translation moved the turbine out of bounds, go to next translation
                            innerflag = 4           #signifies move 1 was attempted
                            #print('Turbine not moved right.')
                    
                        else:       #if there is the turbine is in bounds, evaluate and store
                            CHECK2 = Check_Interference(i)
                            if CHECK2 == 1:         #if interference occurs, move the turbine back , go to next translation
                                Translation_X(-step2, i)
                                innerflag = 4
                                #print('turbine not moved right')
                            else:       #if there is no interference, evaluate and store
                                move1 = Eval_Objective(initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, aif)

                                Clear_Vectors()
                                if move1 >= nomove:             #if evaluation is worse than initial, move back, go to next translation
                                    Translation_X(-step2, i)
                                    innerflag = 4
                                    #print('Turbine not moved right.')
                                else:
                                    flag = 1           #signifies movement was kept
                                #print('turbine ', i, ' moved right.', move1)
                                HubHeight_Search(hstep, i, hstepmin, XLocation, YLocation, initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, rradmin, rradmax, rstep, hubmin, hubmax, rstepmin, aif)
                
                    if innerflag == 4 and flag == 0:        #translation at this step size has resulted in no moves for this turbine
                        turbines[i].Stopped = 1
                        HubHeight_Search(hstep, i, hstepmin, XLocation, YLocation, initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, rradmin, rradmax, rstep, hubmin, hubmax, rstepmin, aif)
                                
                                
                elif turbines[i].YLocation < (LengthY / 2.0):
                    
                    if innerflag == 0 and flag == 0:        #move 1 was just unsucessfully attempted
                        transflag = Translation_Y(-step2, i)
                        CHECK2 = 0
                        if transflag == 1:      #if the translation moved the turbine out of bounds, go to next translation
                            innerflag = 1       #move2 was attempted
                            #print('turbine not moved up.')
                        else:
                            CHECK2 = Check_Interference(i)
                            if CHECK2 == 1:         #if interference occurs, move the turbine back, go to next translation
                                Translation_Y(step2, i)
                                innerflag = 1
                                #print('turbine not moved up.')
                            else:       #if there is no interference, evaluate and store
                                move2 = Eval_Objective(initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, aif)

                                Clear_Vectors()
                                if move2 >= nomove:       #if evaluation is worse than initial, move back, go to next translation
                                    Translation_Y(step2, i)
                                    innerflag = 1
                                    #print('turbine not moved up.')
                                else:       #evaluation is better, keep move, go to next turbine
                                    flag = 1
                                    #print('turbine ', i, ' moved up.', move2)
                                    HubHeight_Search(hstep, i, hstepmin, XLocation, YLocation, initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, rradmin, rradmax, rstep, hubmin, hubmax, rstepmin, aif)
                                
                    if innerflag == 1 and flag == 0:        #move 2 was just unsucessfully attempted
                        transflag = Translation_X(-step2, i)
                        CHECK2 = 0
                        if transflag == 1:      #if the translation moved the turbine out of bounds, go to next translation
                            innerflag = 2       #move3 was attempted
                            #print('turbine not left.')
                        else:
                            CHECK2 = Check_Interference(i)
                            if CHECK2 == 1:         #if interference occurs, move the turbine back, go to next translation
                                Translation_X(step2, i)
                                innerflag = 2
                                #print('turbine not left.')
                            else:       #if there is no interference, evaluate and store
                                move3 = Eval_Objective(initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, aif)

                                Clear_Vectors()
                                if move3 >= nomove:       #if evaluation is worse than initial, move back, go to next translation
                                    Translation_X(step2, i)
                                    innerflag = 2
                                    #print('turbine not moved left.')
                                else:       #evaluation is better, keep move, go to next turbine
                                    flag = 1
                                    #print('turbine ', i, ' moved left.', move3)
                                    HubHeight_Search(hstep, i, hstepmin, XLocation, YLocation, initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, rradmin, rradmax, rstep, hubmin, hubmax, rstepmin, aif)
                                
                    if innerflag == 2 and flag == 0:        #move 3 was just unsucessfully attempted
                        transflag = Translation_Y(step2, i)
                        CHECK2 = 0
                        if transflag == 1:      #if the translation moved the turbine out of bounds, go to next translation
                            innerflag = 3       #move3 was attempted
                            #print('turbine not moved down.')
                        else:
                            CHECK2 = Check_Interference(i)
                            if CHECK2 == 1:         #if interference occurs, move the turbine back, go to next translation
                                Translation_Y(-step2, i)
                                innerflag = 3
                                #print('turbine not moved down.')
                            else:       #if there is no interference, evaluate and store
                                move4 = Eval_Objective(initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, aif)

                                Clear_Vectors()
                                if move4 >= nomove:       #if evaluation is worse than initial, move back, go to next translation
                                    Translation_Y(-step2, i)
                                    innerflag = 3
                                    #print('turbine not moved down.')
                                else:       #evaluation is better, keep move, go to next turbine
                                    flag = 1
                                    #print('turbine ', i, ' moved down.', move4)
                                    HubHeight_Search(hstep, i, hstepmin, XLocation, YLocation, initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, rradmin, rradmax, rstep, hubmin, hubmax, rstepmin, aif)
                    if innerflag == 3 and flag == 0:
                        transflag = Translation_X(step2, i)     #move the turbine one step right
                        CHECK2 = 0
                        if transflag == 1:          #if the translation moved the turbine out of bounds, go to next translation
                            innerflag = 4           #signifies move 1 was attempted
                            #print('Turbine not moved right.')
                    
                        else:       #if there is the turbine is in bounds, evaluate and store
                            CHECK2 = Check_Interference(i)
                            if CHECK2 == 1:         #if interference occurs, move the turbine back , go to next translation
                                Translation_X(-step2, i)
                                innerflag = 4
                                #print('turbine not moved right')
                            else:       #if there is no interference, evaluate and store
                                move1 = Eval_Objective(initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, aif)

                                Clear_Vectors()
                                if move1 >= nomove:             #if evaluation is worse than initial, move back, go to next translation
                                    Translation_X(-step2, i)
                                    innerflag = 4
                                    #print('Turbine not moved right.')
                                else:
                                    flag = 1           #signifies movement was kept
                                #print('turbine ', i, ' moved right.', move1)
                                HubHeight_Search(hstep, i, hstepmin, XLocation, YLocation, initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, rradmin, rradmax, rstep, hubmin, hubmax, rstepmin, aif)
                             
                     
                    if innerflag == 4 and flag == 0:        #translation at this step size has resulted in no moves for this turbine
                        turbines[i].Stopped = 1
                        HubHeight_Search(hstep, i, hstepmin, XLocation, YLocation, initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, rradmin, rradmax, rstep, hubmin, hubmax, rstepmin, aif)
                    
                        #for every turbine in random order
                
                        #check to see if all turbies have not been moved at this step
                        #print_graph()
            exit_css = 0        #exit current step size
            for i in range(0, initial_num):
                exit_css += turbines[i].Stopped
            #print(exit_css)

            if exit_css == initial_num:
                #all turbines have stopped moving at this step size, halving step size.
                #find worst performing turbine and randomly assign elsewhere
                for b in range(0, num_pops):
                    #print("No moves at step size ", step2, " are possible. Popping weakest turbine.")
                    min_power = 5000000.     #initialized to first turbine power output
                    random_vec2 = Rand_Vector(initial_num)    #creates a randomly ordered vector of turbines
                    for j in range(0, initial_num):
                        randorder = random_vec2[j]
                        if turbines[randorder].Power < min_power:
                            min_power = turbines[randorder].Power
                            min_turb = randorder
                                
                    #print('The weakest turbine is turbine ', min_turb, ' with power currently at ', min_power)
                    start_eval = Eval_Objective(initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, aif)

                    Clear_Vectors()
                    initialx = turbines[min_turb].XLocation
                    initialy = turbines[min_turb].YLocation
                    checkx = 0
                    k = 0
                    flag = 0
                    while flag == 0 and k < max_pop_tries:
                        checkx = 0
                        while checkx != 1:      #will try random locations until one has no interference
                            CHECK2 = 0
                            random_X = random.uniform(0, LengthX)
                            xmove = random_X
                            turbines[min_turb].XLocation = xmove
                            random_Y = random.uniform(0, LengthY)
                            ymove = random_Y
                            turbines[min_turb].YLocation = ymove
                                
                            CHECK2 = Check_Interference(min_turb)
                            if CHECK2 != 1:         #No interference
                                checkx = 1          #place turbine and exit poping loop
                                #print('Turbine ', min_turb, ' has moved to a new location.')
                            else:
                                turbines[min_turb].XLocation = initialx
                                turbines[min_turb].YLocation = initialy
                                #print('Turbine cannot be relocated without interference, trying agian.')
                            

                        new_eval = Eval_Objective(initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, aif)

                        Clear_Vectors()
                        if new_eval < start_eval:
                            flag = 1
                            HubHeight_Search(hstep, i, hstepmin, XLocation, YLocation, initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, rradmin, rradmax, rstep, hubmin, hubmax, rstepmin, aif)
                            #print('Move has improved the evaluation. Continuing pattern serach.')
                        else:
                            turbines[min_turb].XLocation = initialx
                            turbines[min_turb].YLocation = initialy
                            #print('Move did not improve evaluation. Trying new moves.')
                        k += 1
                        
                    #print('pops till acceptance= ', k)
                                
                step2 = step2 / 2.0
                #print('step sized is being reduced to: ', step2)
                        
            #else:
                #print('Turbines were not all placed at step size ', step2, '. Repeating.')
                        
############################################################################################################################

#################################################################################################################
def HubHeight_Search(hstep, i, hstepmin, XLocation, YLocation, initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, rradmin, rradmax, rstep, hubmin, hubmax, rstepmin, aif):            
    step = hstep
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
def Discretize_RSA(i):
    XCenter = turbines[i].XLocation
    ZCenter = turbines[i].HubHeight
    rad = turbines[i].RotorRad
    turbines[i].xcoords = []
    turbines[i].zcoords = []
    xcoords = []
    zcoords = []
    
    #center row
    #center point
    xcoords.append(XCenter)
    zcoords.append(ZCenter)
    
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
    turbines[i].xcoords = xcoords
    turbines[i].zcoords = zcoords

#############################################################################################################
def RotorRadius_Search(i, rradmin, rradmax, rstep, initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, rstepmin, aif):
    #now that an optimal hub height has been found, perform sub-level EPS for rotor radius
    #Clear_Vectors()
    
    hubh = turbines[i].HubHeight
    #rrad = turbines[i].RotorRad         #initialized to 20m
    minrange = (hubh * (2.0 / 3.0)) / 2.0   #the minimum rotor diameter should be 2/3 of the hub height
    maxrange = hubh / 2.0       #the maximum rotor diameter should be 1x the hub height
    
    if minrange < rradmin:      #set min as larger of 19m OR ratio parameter
        minrange = rradmin
        #print('minimum range for rotor radius is ', minrange)
        
    if maxrange > rradmax:      #set max as smaller of 67 OR ratio parameter
        maxrange = rradmax
        #print('maximum range for rotor radius is ', maxrange)
        
    ranger = maxrange - minrange
    centerrad = minrange + (ranger / 2.0)       #center value of range of rotor radii
    #currentr = turbines[i].RotorRad

    turbines[i].RotorRad = centerrad        #set rotor radius to cener valid value
    
    step2 = rstep
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
            redcx.append(turbines[i].XLocation)
            redcy.append(turbines[i].YLocation)

        elif turbines[i].HubHeight <= 60 and turbines[i].RotorRad <= 40:
            yellowcx.append(turbines[i].XLocation)
            yellowcy.append(turbines[i].YLocation)

        elif turbines[i].HubHeight <= 60 and turbines[i].RotorRad <= 60:
            greencx.append(turbines[i].XLocation)
            greency.append(turbines[i].YLocation)

        elif turbines[i].HubHeight <= 60 and turbines[i].RotorRad > 60:
            bluecx.append(turbines[i].XLocation)
            bluecy.append(turbines[i].YLocation)

        elif turbines[i].HubHeight <= 80 and turbines[i].RotorRad <= 30:
            redtrx.append(turbines[i].XLocation)
            redtry.append(turbines[i].YLocation)

        elif turbines[i].HubHeight <= 80 and turbines[i].RotorRad <= 40:
            yellowtrx.append(turbines[i].XLocation)
            yellowtry.append(turbines[i].YLocation)

        elif turbines[i].HubHeight <= 80 and turbines[i].RotorRad <= 60:
            greentrx.append(turbines[i].XLocation)
            greentry.append(turbines[i].YLocation)

        elif turbines[i].HubHeight <= 80 and turbines[i].RotorRad > 60:
            bluetrx.append(turbines[i].XLocation)
            bluetry.append(turbines[i].YLocation)

        elif turbines[i].HubHeight <= 120 and turbines[i].RotorRad <= 30:
            redrhx.append(turbines[i].XLocation)
            redrhy.append(turbines[i].YLocation)

        elif turbines[i].HubHeight <= 120 and turbines[i].RotorRad <= 40:
            yellowrhx.append(turbines[i].XLocation)
            yellowrhy.append(turbines[i].YLocation)

        elif turbines[i].HubHeight <= 120 and turbines[i].RotorRad <= 60:
            greenrhx.append(turbines[i].XLocation)
            greenrhy.append(turbines[i].YLocation)

        elif turbines[i].HubHeight <= 120 and turbines[i].RotorRad > 60:
            bluerhx.append(turbines[i].XLocation)
            bluerhy.append(turbines[i].YLocation)

        elif turbines[i].HubHeight > 120 and turbines[i].RotorRad <= 30:
            redsqx.append(turbines[i].XLocation)
            redsqy.append(turbines[i].YLocation)
  
        elif turbines[i].HubHeight > 120 and turbines[i].RotorRad <= 40:
            yellowsqx.append(turbines[i].XLocation)
            yellowsqy.append(turbines[i].YLocation)

        elif turbines[i].HubHeight > 120 and turbines[i].RotorRad <= 60:
            greensqx.append(turbines[i].XLocation)
            greensqy.append(turbines[i].YLocation)

        elif turbines[i].HubHeight > 120 and turbines[i].RotorRad > 60:
            bluesqx.append(turbines[i].XLocation)
            bluesqy.append(turbines[i].YLocation)

        else:
            noplotx.append(turbines[i].XLocation)
            noploty.append(turbines[i].YLocation)

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
    #plt.axis([0, 2000, 0, 2000])
    plt.ylabel('Position (m)')
    plt.xlabel('Position (m)')
    plt.title(str('Optimization of ' + str(initial_num) + ' Turbines'))
    plt.savefig(str(str(initial_num) + 'turbinesWithDisc.png'), bbox_inches='tight')
        
#################################################################################################################

#############################################################################################################
#calculating distance between turbines
def length(j, k):
    x1 = turbines[j].XLocation
    x2 = turbines[k].XLocation
    y1 = turbines[j].YLocation
    y2 = turbines[k].YLocation
    
    term1 = (x1 - x2) ** 2
    term2 = (y1 - y2) ** 2
    length = math.sqrt(term1 + term2)
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
    if initial_num < 11:
        hood_size = initial_num - 1
    else:
        hood_size = 10
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
        xnew = [turbines[j].XLocation, turbines[k].XLocation]
        ynew = [turbines[j].YLocation, turbines[k].YLocation]
        ax1.plot(xnew, ynew, c='k')
        ax1.scatter(xnew, ynew, s=10, c='r', marker="o")
        plt.axis([0, 2000, 0, 2000])
        plt.ylabel('Position (m)')
        plt.xlabel('Position (m)')
        plt.title(str('Cable Layout for ' + str(initial_num) + ' Turbines'))
        
##########################################################################################################################
rows = 4
cols = 4
initExtent = 0.95
radius = 40.
siteX = 15. * 80. ###must be much smaller than LocationX
siteY = 15. * 80. ###must be much smaller than LocationY

xpos = np.linspace(-initExtent*(siteX - radius),initExtent*(siteX - radius),cols)
ypos = np.linspace(-initExtent*(siteY - radius),initExtent*(siteY - radius),rows)
mx = []
my = []
for i in range(rows):
    for j in range(cols):
        mx.append(xpos[j] - xpos[0])
        my.append(ypos[i] - ypos[0])

layouts_x = [mx]
layouts_y = [my]

mx = []
my = []
for i in range(rows):
    for j in range(cols):
        mx.append(xpos[j])
        my.append(ypos[i])
        
layouts_x.append(mx)
layouts_y.append(my)

#optimized layout
#layouts_x.append([1304.2130491853543, 303.16195377176825, 1882.7665062612382, 734.1365122392898, 36.582152829805295, 1078.6666766839614, 1779.34556196998, 417.79552464946823, 1986.6903503041713, 133.48666173166612, 1430.1977444645, 954.8374274350235, 839.2773010244406, 1186.314920647951, 518.9482768461385, 1681.3296876680367, 621.3034610992612, 1567.283691418906])
#layouts_y.append([204.01944774702838, 270.02673301978496, 276.97611783621744, 305.2706459401636, 1065.3030786448571, 104.67022680498235, 610.6008734856377, 450.6469581295796, 32.94310546319856, 1290.6284216353858, 413.28638998844144, 611.4810784612052, 117.69510957163021, 489.9943326680768, 225.4370943493467, 858.4028566999443, 505.9650483803546, 1212.01937488221])

#line layout
#layouts_x.append([20., 20., 20., 20., 20., 20., 20., 20., 20., 20.])
#layouts_y.append([0., 200., 400., 600., 800., 1000., 1200., 1400., 1600., 1800.])
#layouts_y.append([0., 2., 4., 6., 8., 10., 12., 14., 16., 18.])
if output == 'on':
    enum = 1
    while os.path.isfile('3D_Park BreifDatapt' + str(enum) + '.txt'):
        enum +=1
    excel = open('3D_Park BreifDatapt' + str(enum) + '.txt', 'w') 
    
    excel.write(str('number of turbines, final evaluation, power output, efficiency, number of evaluations made, input file line number \n'))
#i = 1
#initial_num = 15
for line in range(0, len(layouts_x)):
    initial_num = len(layouts_x[line])
    #initial_num = int(initial_num_file.readline())
    if output == 'on':
        trialnum = 1
        while os.path.isfile('3D_Park_' + str(initial_num) + 'Test ' + str(trialnum) + '.txt'):
            trialnum += 1
        data_out = open('3D_Park_' + str(initial_num) + 'Test ' + str(trialnum) + '.txt', 'w')
    
    XLocation = layouts_x[line]
    YLocation = layouts_y[line]
    ZLocation = [0.0] * initial_num
    HubHeight = [80.0] * initial_num
    OldHubHeight = [80.0] * initial_num
    ZHub = []
    RotorRad = [40.0] * initial_num
    OldRotorRad = [40.0] * initial_num
    alpha = [0.0] * initial_num
    usturbinesrec = []
    usturbines = []
    dsturbinesrec = []
    dsturbines = []
    wakewidth = []
    distance = []
    percent = []
    ui = 12.0     #changes in Compute_Wake function
    windspeeds = []
    xcoords = []
    zcoords = []
    Power = 0.0
    Area = 0

    for i in range(0, initial_num):
        ZHub.append(HubHeight[i] + ZLocation[i])

    #initialize turbines         
    turbines = [Turbine(XLocation[i], YLocation[i], ZLocation[i], HubHeight[i], RotorRad[i], alpha[i], usturbinesrec, usturbines, dsturbinesrec, dsturbines, wakewidth, distance, percent, ui, windspeeds, xcoords, zcoords, Power, Area) for i in range(0, initial_num)]        
     #Evaluate Initial Layout
    score = Eval_Objective(initial_num, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, aif)       
    if output == 'on':
        data_out.write(str('X Locations: ' + str(XLocation) + '\n'))
        data_out.write(str('Y Locations: ' + str(YLocation) + '\n'))
        data_out.write(str('Z Locations: ' + str(ZLocation) + '\n'))
        data_out.write(str('hub heights: ' + str(HubHeight) + '\n'))
        data_out.write(str('rotor radii: ' + str(RotorRad) + '\n'))
    
       
        data_out.write(str('The layout has a score of: ' + str(score) + '\n'))

    total = 0.0
    possible = 0.0
    windspeeds = 0.0

    for i in range(0, initial_num):
        total += turbines[i].Power
        Area = math.pi * ((turbines[i].RotorRad) ** 2.)
        hubheight = turbines[i].HubHeight
        Uiii = U0 * ((hubheight / Zref) ** alphah)
        if Uiii < Cutin:
            possible += 0.
        elif Uiii < Uref:
            possible += 0.5 * (Uiii ** 3) * ro * Area * Cp * Cf / 1000.
        else:
            possible += 0.5 * (Uref ** 3) * ro * Area * Cp * Cf / 1000.
    efficiency = total/possible * 100
    print('efficiency = ', efficiency)
    print('total power = ', total)
    #print('Farm Output Power is ', total, 'Watts')
    if output == 'on':
        for i in range(0, initial_num):
            data_out.write(str('The effective windspeed for turbine ' + str(i) + ' is ' + str(turbines[i].ui) + '\n'))
            data_out.write(str('Turbine ' + str(i) + ' power: ' + str(turbines[i].Power) + '\n'))
     
        data_out.write(str('The total power generated was: ' + str(total) + ' kW. \n'))    

        #rorad = turbines[k].RotorRad
        #Area = math.pi * (rorad ** 2.0)
        #temp1 = 0.3 * (12 ** 3.0)
        #possible += temp1
        #efficiency = (total / possible) * 100
        #print('Farm efficiency is ', efficiency, '%.')
        data_out.write(str('The efficiency of the layout was: ' + str(efficiency) + '%. \n'))
        data_out.write(str('The Objective was calculated ' + str(nemo) + ' times. \n'))
        
        excel.write(str(str(initial_num) + ',' + str(score) + ',' + str(total) + ',' + str(efficiency) + ',' + str(nemo) + ',' + str(line + 1) + '\n'))
    
    #print_graph()
    #plot_cables(networks)
    sum_jensen = 0.
    for i in turbines:
        sum_jensen += i.windspeeds[0] ** 3
    #    print('turbine ',i,' windspeed: ', i.windspeeds[0])
    print('sum_jensen: ',sum_jensen)
    if output == 'on':
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
wdsps = [3.32213844322,
3.28122913386,
3.26099178756,
3.32506440617,
3.75566888684,
3.77181817163,
3.77094739129,
3.78139110206,
4.314573208,
4.30851591476,
4.34337267494,
4.30917666,
5.22184407986,
5.30417856035,
5.18963032343,
5.27707050801]
sum_CFD = 0.
for i in wdsps:
    sum_CFD += i ** 3
print('sum_CFD: ',sum_CFD)

#initial_num += 1
if output == 'on':
    excel.close()

#data_out.close()
#initial_num_file.close()
#print('files are closed')
#for i in range(0, initial_num):
#    print([i, turbines[i].XLocation, turbines[i].YLocation])
