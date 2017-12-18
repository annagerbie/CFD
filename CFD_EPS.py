
__author__ = "Ryan King <ryan.king@nrel.gov>"
__date__ = "2017-11-21"
__copyright__ = "Copyright (c) NREL 2017. All rights reserved."
__license__  = "Apache License, Version 2.0"


from dolfin import *
import numpy as np
import random as random
from scipy import integrate
from scipy import stats
import time
import sys
import csv
from math import sqrt

sys.getrecursionlimit()
sys.setrecursionlimit(10000)

set_log_level(INFO) ###edits what gets printed out
# mpirun -np 4 python windse4dir.py 2>&1 | tee log
'''
CRITICAL  = 50, // errors that may lead to data corruption and suchlike
ERROR     = 40, // things that go booms
WARNING   = 30, // things that may go boom later
INFO      = 20, // information of general interest
PROGRESS  = 16, // what's happening (broadly)
TRACE     = 13, // what's happening (in detail)
DBG       = 10  // sundry

'''
air_density = 1.225 #kg/m^3
init_step = 400. #initial step size
minstep = 3. #minimum step size
num_pops = 5
max_pop_tries = 1000
# Print log messages only from the root process in parallel
parameters["std_out_all_processes"] = False; ### Print log messages only from the root process in parallel
parameters['form_compiler']['cpp_optimize_flags'] = '-O3 -fno-math-errno -march=native'        
parameters["form_compiler"]["optimize"]     = True
parameters["form_compiler"]["cpp_optimize"] = True
parameters['form_compiler']['representation'] = 'uflacs'
parameters['form_compiler']['quadrature_degree'] = 12

#actual turbine diam in meters
RD = 80.

# mesh parameters
degree = 2
Lx = 60.*RD ###size of mesh (NOT FARM)
Ly = 60.*RD ###size of mesh (NOT FARM)
nx = 50 ###number of points in mesh
ny = 50 ###number of points in mesh

#WTG parameters
numturbs = 16

inflowVel=8 ### ambient wind speed?
#number of inflow direction bins
bins = 1 ### Unidirectional wind
WTGexp = 8. ### gamma for smoothing kernel (Eq. 12)
radius = RD/2.
thickness = RD/20.
numRefine = 2
A=RD # weird for 2D  ### WHAT's this actually? - area b/c it's got no z?
HH=80    
initExtent=.95 ### create buffer on grid layout
mlDenom=2. ### mixing length denom - not sure how this is arrived at


#site/refinement constraints
site_x = 15.*RD ###size of farm
site_y = 15.*RD ###size of farm

restart = False ###seeded layout (must go into layout function to change seed)
randStart = False ###random layout
gridStart = True ###gridded layout - must have 16 turbines
optimize = False ###
linearStart = False ###turbines in a row

if optimize == True:
    from dolfin_adjoint import *


# inflow is always from left in simulation
dirs = np.linspace(pi/2., 2*pi+pi/2., bins, endpoint = False)
weights = np.ones(bins)/bins ### assume even weightings


# actuator disk distribution for normalization constant
def WTGdist(x,y): ###smoothing kernel
    return np.exp(-((x/thickness)**WTGexp + (y/radius)**WTGexp))

def refine_mesh(mesh, site_x, site_y, HH):
    #refines the mesh around the site boundaries
    h = mesh.hmin()
    
    cell_f = CellFunction('bool', mesh, False) ### create a mesh of the same size, where every cell is set to 'False'
    for cell in cells(mesh):
        if (cell.midpoint()[0]**2 + cell.midpoint()[1]**2 < site_x**2 + site_y**2 + h) :
            cell_f[cell] = True
            ###if the midpoint of the cell is within the circumradius of the farm, change the cell's value to true

    mesh = refine(mesh, cell_f) ###refine the cells within the circumradius of the farm

    return mesh

def createLayout(numturbs):
    mx=[]
    my=[]
    mz=[]

    ###if you've specified a random starting array
    if randStart == True:
        for i in range(numturbs):
            ### select a random x location, within the farm bounds, and such that the entire rotor will be within the farm bounds
            mx.append(Constant(np.random.uniform(low=-(site_x - radius),high=(site_x - radius))))
            ### select a random y location, within the farm bounds, and such that the entire rotor will be within the farm bounds
            my.append(Constant(np.random.uniform(low=-(site_y - radius), high=(site_y - radius))))
            ### single hub height, z = HH
            mz.append(Constant(HH))
    elif gridStart == True:
        ###if you've specified a grided starting array 
        ###Can only accept 16 turbines
        if numturbs == 16:
            rows = 4
            cols = 4
            xpos = np.linspace(-initExtent*(site_x - radius),initExtent*(site_x - radius),cols)
            ypos = np.linspace(-initExtent*(site_y - radius),initExtent*(site_y - radius),rows)
            #print('xlocations')
            #print(xpos)
            for i in range(rows):
                for j in range(cols):
                    mx.append(Constant(xpos[j]))
                    my.append(Constant(ypos[i]))
                    # # some starting noise sometimes helps
                    # mx.append(Constant(xpos[j]+5.*np.random.randn()))
                    # my.append(Constant(ypos[i]+5.*np.random.randn()))
                    mz.append(Constant(HH))
                    
    
    elif linearStart == True:
        '''My adds'''
        mx = np.linspace(-initExtent*(site_x - radius),initExtent*(site_x - radius),numturbs)
        my = np.array([0] * numturbs)
        mz = np.array([HH] * numturbs)
        
    elif offgridStart == True:
        pass
    elif restart == True:
        # fixed layout here
        m_temp = [Constant(-113.961988283),Constant(-386.535837904),Constant(-512.116113959),Constant(-237.354391531),Constant(638.697968355),Constant(13.6826901448),Constant(386.535838424),Constant(-113.961987466),Constant(13.6826875361),Constant(-638.697971072),Constant(-887.942379804),Constant(-813.542880381),Constant(813.542880031),Constant(-887.942379852),Constant(237.354391629),Constant(-512.116113931),Constant(-237.3543916),Constant(512.116113865),Constant(-813.542880345),Constant(887.942379783),Constant(887.942379753),Constant(813.542880265),Constant(-13.6826884631),Constant(638.697970038),Constant(-386.535837846),Constant(113.961988218),Constant(-638.697970958),Constant(-13.6826879195),Constant(512.116113711),Constant(237.354391612),Constant(113.961988),Constant(386.535838129)]
        ###Splits a listinto x and y coordinates, function is shown further down
        mx,my = splitSolution(m_temp,numturbs)
        for i in range(numturbs):
            mz.append(Constant(HH))

    return mx, my, mz

def createRotatedTurbineForce(mx,my,ma,A,beta,numturbs,alpha,V):
    #mx, my = x, y coordinates of each turbine
    #ma = vector of value 0.33 for each turbine - axial induction factor
    # A = RD at beginning of code. Maybe area since we're using 2D...
    # beta = integral over actuator disk area in x and y: used to normalize
    x=SpatialCoordinate(mesh)
    WTGbase = project(Expression(("1.0","0.0"),degree=2),V)
    tf = Function(V)

    for i in range(numturbs):
        #rotation
        xrot = cos(alpha)*mx[i] - sin(alpha)*my[i]
        yrot = sin(alpha)*mx[i] + cos(alpha)*my[i]

        ### force imported on the flow by each wind turbine = 0.5 * rho * An*c't*smoothing kernal / beta - no multiplication by magnitude because we assume turbine is turned into wind
        ### why no windspeed included? why no rho included?
        tf = tf + 0.5*4.*A*ma[i]/(1.-ma[i])/beta*exp(-(((x[0] - xrot)/thickness)**WTGexp + ((x[1] - yrot)/radius)**WTGexp))*WTGbase.copy(deepcopy=True)

    return tf

def rotatedPowerFunctional(alpha,A,beta,mx,my,ma,u,numturbs,V):
    ### used only in optimization
    # functional for dolfin-adjoint
    x=SpatialCoordinate(mesh)
    J=Functional(0.)
    for i in range(numturbs):
        #rotation
        xrot = cos(alpha)*mx[i] - sin(alpha)*my[i]
        yrot = sin(alpha)*mx[i] + cos(alpha)*my[i]

        J = J + 0.5*4.*np.pi*radius**2*ma[i]/(1.-ma[i])/beta*Functional(exp(-(((x[0] - xrot)/thickness)**WTGexp + ((x[1] - yrot)/radius)**WTGexp))*u[0]**3.*dx) ### /beta added by Annalise 12/13/17

    return J

def rotatedPowerFunction(alpha,A,beta,mx,my,ma,up,numturbs,V):
    ### used in developing final power
    #emulating an actual power curve
    x=SpatialCoordinate(mesh)
    J = []
    for i in range(numturbs):
        #rotation
        xrot = cos(alpha)*mx[i] - sin(alpha)*my[i] ### -5 added by Annalise 12/14 to understand effects of smoothing kernal on usturbine
        yrot = sin(alpha)*mx[i] + cos(alpha)*my[i]
        print(up.sub(0)(xrot,yrot)[0])
        #J = 0.5*np.pi*radius**2*4*float(ma[i])*(1.-float(ma[i]))**2*up.sub(0)(xrot,yrot)[0]**3 
        J.append(0.5*air_density*np.pi*radius**2*4*float(ma[i])/(1.-float(ma[i]))*up.sub(0)(xrot,yrot)[0]**3)
        ### up.sub(0)(xrot,yrot)[0] --> up == u and p combined --> sub(0) == just u --> (xrot, yrot) == position of interest (center pt) --> [0] == x-velocity

    return J


def createControl(mx,my,numturbs):
        m = []
        for i in range(numturbs):
            m.append(Control(mx[i]))
            m.append(Control(my[i]))
            
        return m

def createBounds(mx,my,numturbs):
        ub=[]
        lb=[]
        for i in range(numturbs):
            lb.append(Constant(-(site_x - radius)))
            lb.append(Constant(-(site_y - radius)))
            ub.append(Constant((site_x - radius)))
            ub.append(Constant((site_y - radius)))
            
        bounds = [lb,ub]
        return bounds

def createAxialBounds(ma,numturbs):
        ub=[]
        lb=[]
        for i in range(numturbs):
            lb.append(Constant(0.))
            ub.append(Constant(0.75))
            
        bounds = [lb,ub]
        return bounds

def splitSolution(m_opt,numturbs):
    mx_opt = []
    my_opt = []
    j=0
    for i in range(numturbs):
        mx_opt.append(m_opt[j])
        j+=1
        my_opt.append(m_opt[j])
        j+=1
        
    return mx_opt,my_opt

def main(tf):
    nu = Constant(.00005) ### kinematic viscosity
    f = Constant((0.,0.))
    up_next = Function(VQ) ### up_next becomes tuple of vector and finite elements?
    u_next,p_next = split(up_next) ### split vector (wind speed) and finite (pressure) elements?
    v,q = TestFunctions(VQ)
    class InitialConditions(Expression): ###inherits from Expression class in fenics
        def __init__(self,**kwargs):
            random.seed(2) ### WHY IS THIS HERE?
        def eval(self, values, x):
            values[0] = inflowVel
            values[1] = 0.0
            values[2] = 0.0
        def value_shape(self):
            return (3,)  
    #boundary conditions
    class NoSlipBoundary(SubDomain):
        def inside(self, x, on_boundary): ### windspeed has w = 0 at top and bottom
            return near(x[1]**2 - (Ly/2.)**2, 0.) and on_boundary

    class InflowBoundary(SubDomain): ### windspeeed has u = inflow velocity and w = 0 at locations of inflow
        def inside(self, x, on_boundary):
            return near(x[0],-(Lx/2.)) and on_boundary

    class PeriodicBoundary(SubDomain):

        def inside(self, x, on_boundary):
            # return True if on left or bottom boundary AND NOT on one of the two slave edges
            return bool((near(x[0], -(Lx/2.)) or near(x[1], -(Ly/2.))) and 
                    (not (near(x[0], (Lx/2.)) or near(x[1], (Ly/2.)))) and on_boundary)
                          
        def map(self, x, y):
            if near(x[0], (Lx/2.)) and near(x[1], (Ly/2.)):
                y[0] = x[0] - 2*(Lx/2.)
                y[1] = x[1] - 2*(Ly/2.)
            elif near(x[0], (Lx/2.)):
                y[0] = x[0] - 2*(Lx/2.)
                y[1] = x[1]
            else: # near(x[1], (Ly/2.)):
                y[0] = x[0]
                y[1] = x[1] - 2*(Ly/2.)

    u0 = InitialConditions(degree=2)  
    up_next.interpolate(u0)

    lmix = radius/mlDenom ### mixing lenth
    ### mean rate of strain tensor ... so confused
    S = sqrt(2.*inner(0.5*(grad(u_next)+grad(u_next).T),0.5*(grad(u_next)+grad(u_next).T)))
    nu_T=lmix**2.*S ### eddie viscosity

    F = inner(grad(u_next)*u_next, v)*dx + (nu+nu_T)*inner(grad(u_next), grad(v))*dx - inner(div(v),p_next)*dx - inner(div(u_next),q)*dx - inner(f,v)*dx + inner(tf*u_next[0]**2,v)*dx 
    ### u_next[0] ** 2 because it's not added to the force on f_AD calc?
    ### I think this is ok without density --> http://libelemental.org/featured/2011/pdesys.html
    
    # lateral BC
    # bc1 = DirichletBC(VQ.sub(0), Constant((8.0,0.0)), NoSlipBoundary())
    bc1a = DirichletBC(VQ.sub(0).sub(1), Constant(0.0), NoSlipBoundary())

    # inflow BC
    bc2 = DirichletBC(VQ.sub(0), Constant((inflowVel,0.0)), InflowBoundary())
    # bc2a = DirichletBC(VQ.sub(0).sub(0), Constant(8.), InflowBoundary())

    # outflow pressure BC is implicitly 0

    bc = [bc1a,bc2]

    solve(F == 0, up_next, bc, solver_parameters={"newton_solver":{"absolute_tolerance": 1e-8}})
    u_next,p_next = split(up_next)

    if optimize == False:
        nu_T_out=project(nu_T, Q)
        lStr= 'nu_t.pvd'
        file = File(lStr)
        file << nu_T_out
    return u_next, up_next
###############################################################################
def Eval_Objective(mx,my,ma,A,B,numturb):
    J = 0.
    cumulative_power = [0.] * numturbs
    for i in range(bins): ###for each wind direction
        tf_rot= createRotatedTurbineForce(mx,my,ma,A,B,numturbs,dirs[i],V) ### calculate force imparted by turbines
        u_rot, up_rot = main(tf_rot)  ### RANS solver
        power_dev = rotatedPowerFunction(alpha,A,beta,mx,my,ma,up_rot,numturbs,V)
        J = J + weights[i]*sum(power_dev)
        cumulative_power = [i+j for i,j in zip(power_dev,cumulative_power)]
    
    return J, cumulative_power
#######################################################################################################################
def Rand_Vector(numturbs):
    random_vec = []
    for i in range(0, numturbs):    
        random_vec.append(i)
        
    #shuffle elements by randomly exchanging each with one other
    for i in range(0, len(random_vec)):
        r = random.randint(0, len(random_vec)-1)  #select random value in vector
        temp = random_vec[i]
        random_vec[i] = random_vec[r]
        random_vec[r] = temp
    return random_vec
###################################################################################################################
def Translation_X(step_size, i, mx):
    transflag = 0
    #print(xstart, step_size)
    xfinish = mx[i] + step_size            #find preliminary new x-location given step size
    if xfinish >= 0 and xfinish <= site_x:    #if this new x-location is not out of bounds, translate it
        mx[i] = xfinish
        return mx, transflag
    else:
        transflag = 1
        return mx, transflag
################################################################################################################        
def Translation_Y(step_size, i, my):
    transflag = 0
    yfinish = my[i] + step_size            #find preliminary new x-location given step size
    if yfinish >= 0 and yfinish <= site_y:    #if this new x-location is not out of bounds, translate it
        my[i] = yfinish   
        return my, transflag
    else:
        transflag = 1
        return my, transflag  
###################################################################################################################    
def Check_Interference(mx, my, n):
    CHECK2 = 0
    x = mx[n]
    y = my[n]
    for k in range(0, numturbs):
        if k != n:
            xnew = mx[k]
            ynew = my[k]
            checkx = x - xnew
            checky = y - ynew
            checkrad = sqrt(checkx ** 2.0 + checky ** 2.0)
            if checkrad < 200:
                CHECK2 = 1
    return CHECK2
###############################################################################
def EPS(mx, my, mz, ma):
    Stopped = [0.] * numturbs
    for h in range(0, 1):
        step2 = init_step
        while step2 >= minstep:
            random_vec = Rand_Vector(numturbs)    #creates a randomly ordered vector of turbines
            for j in range(0, len(random_vec)):
                i = random_vec[j]
                Stopped[i] == 0
                print('Turbine ', i,' is being tested.')
                flag = 0
                innerflag = 0
                transflag = 0
                nomove, turb_power = Eval_Objective(mx,my,ma,A,B,numturbs)
                #develop preliminary objective for comparison purposes                
                print('stepped into while loop')
                if my[i] >= 0.:
                    if innerflag == 0 and flag == 0:        #move 1 was just unsucessfully attempted
                        transflag = Translation_Y(step2, i, my)
                        CHECK2 = 0
                        if transflag == 1:      #if the translation moved the turbine out of bounds, go to next translation
                            innerflag = 1       #move2 was attempted
                            print('turbine not moved up.')
                        else:
                            CHECK2 = Check_Interference(mx, my, i)
                            if CHECK2 == 1:         #if interference occurs, move the turbine back, go to next translation
                                Translation_Y(-step2, i, my)
                                innerflag = 1
                                print('turbine not moved up.')
                            else:       #if there is no interference, evaluate and store
                                move2, turb_power2 = Eval_Objective(mx,my,ma,A,B,numturbs)

                                if move2 >= nomove:       #if evaluation is worse than initial, move back, go to next translation
                                    Translation_Y(-step2, i, my)
                                    innerflag = 1
                                    print('turbine not moved up.')
                                else:       #evaluation is better, keep move, go to next turbine
                                    flag = 1
                                    turb_power = [i for i in turb_power2]
                                    objective = move2 * 1.
                                    print('turbine ', i, ' moved up.', move2)
                                    #HubHeight_Search(hstep, i, hstepmin, XLocation, YLocation, numturbs, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, rradmin, rradmax, rstep, hubmin, hubmax, rstepmin, aif)
                                
                    if innerflag == 1 and flag == 0:        #move 2 was just unsucessfully attempted
                        transflag = Translation_X(-step2, i, mx)
                        CHECK2 = 0
                        if transflag == 1:      #if the translation moved the turbine out of bounds, go to next translation
                            innerflag = 2       #move3 was attempted
                            print('turbine not left.')
                        else:
                            CHECK2 = Check_Interference(i)
                            if CHECK2 == 1:         #if interference occurs, move the turbine back, go to next translation
                                Translation_X(step2, i, mx)
                                innerflag = 2
                                print('turbine not left.')
                            else:       #if there is no interference, evaluate and store
                                move3, turb_power2 = Eval_Objective(mx,my,ma,A,B,numturbs)

                                if move3 >= nomove:       #if evaluation is worse than initial, move back, go to next translation
                                    Translation_X(step2, i, mx)
                                    innerflag = 2
                                    print('turbine not moved left.')
                                else:       #evaluation is better, keep move, go to next turbine
                                    flag = 1
                                    turb_power = [i for i in turb_power2]
                                    objective = move3 * 1.
                                    print('turbine ', i, ' moved left.', move3)
                                    #HubHeight_Search(hstep, i, hstepmin, XLocation, YLocation, numturbs, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, rradmin, rradmax, rstep, hubmin, hubmax, rstepmin, aif)
                                
                    if innerflag == 2 and flag == 0:        #move 3 was just unsucessfully attempted
                        transflag = Translation_Y(-step2, i, my)
                        CHECK2 = 0
                        if transflag == 1:      #if the translation moved the turbine out of bounds, go to next translation
                            innerflag = 3       #move3 was attempted
                            print('turbine not moved down.')
                        else:
                            CHECK2 = Check_Interference(i)
                            if CHECK2 == 1:         #if interference occurs, move the turbine back, go to next translation
                                Translation_Y(step2, i, my)
                                innerflag = 3
                                print('turbine not moved down.')
                            else:       #if there is no interference, evaluate and store
                                move4, turb_power2 = Eval_Objective(mx,my,ma,A,B,numturbs)

                                if move4 >= nomove:       #if evaluation is worse than initial, move back, go to next translation
                                    Translation_Y(step2, i, my)
                                    innerflag = 3
                                    print('turbine not moved down.')
                                else:       #evaluation is better, keep move, go to next turbine
                                    flag = 1
                                    turb_power = [i for i in turb_power2]
                                    objective = move4 * 1.
                                    print('turbine ', i, ' moved down.', move4)
                                    #HubHeight_Search(hstep, i, hstepmin, XLocation, YLocation, numturbs, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, rradmin, rradmax, rstep, hubmin, hubmax, rstepmin, aif)
                    if innerflag == 3 and flag == 0:
                        transflag = Translation_X(step2, i, mx)     #move the turbine one step right
                        CHECK2 = 0
                        if transflag == 1:          #if the translation moved the turbine out of bounds, go to next translation
                            innerflag = 4           #signifies move 1 was attempted
                            print('Turbine not moved right.')
                    
                        else:       #if there is the turbine is in bounds, evaluate and store
                            CHECK2 = Check_Interference(i)
                            if CHECK2 == 1:         #if interference occurs, move the turbine back , go to next translation
                                Translation_X(-step2, i, mx)
                                innerflag = 4
                                print('turbine not moved right')
                            else:       #if there is no interference, evaluate and store
                                move1, turb_power2 = Eval_Objective(mx,my,ma,A,B,numturbs)

                                if move1 >= nomove:             #if evaluation is worse than initial, move back, go to next translation
                                    Translation_X(-step2, i, mx)
                                    innerflag = 4
                                    print('Turbine not moved right.')
                                else:
                                    flag = 1           #signifies movement was kept
                                    turb_power = [i for i in turb_power2]
                                    objective = move1 * 1.
                                print('turbine ', i, ' moved right.', move1)
                                #HubHeight_Search(hstep, i, hstepmin, XLocation, YLocation, numturbs, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, rradmin, rradmax, rstep, hubmin, hubmax, rstepmin, aif)
                
                    if innerflag == 4 and flag == 0:        #translation at this step size has resulted in no moves for this turbine
                        Stopped[i] = 1
                        objective = nomove * 1.
                        #HubHeight_Search(hstep, i, hstepmin, XLocation, YLocation, numturbs, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, rradmin, rradmax, rstep, hubmin, hubmax, rstepmin, aif)
                                
                                
                elif my[i] < 0.:
                    
                    if innerflag == 0 and flag == 0:        #move 1 was just unsucessfully attempted
                        transflag = Translation_Y(-step2, i, my)
                        CHECK2 = 0
                        if transflag == 1:      #if the translation moved the turbine out of bounds, go to next translation
                            innerflag = 1       #move2 was attempted
                            print('turbine not moved up.')
                        else:
                            CHECK2 = Check_Interference(i)
                            if CHECK2 == 1:         #if interference occurs, move the turbine back, go to next translation
                                Translation_Y(step2, i, my)
                                innerflag = 1
                                print('turbine not moved up.')
                            else:       #if there is no interference, evaluate and store
                                move2, turb_power2 = Eval_Objective(mx,my,ma,A,B,numturbs)

                                if move2 >= nomove:       #if evaluation is worse than initial, move back, go to next translation
                                    Translation_Y(step2, i, my)
                                    innerflag = 1
                                    print('turbine not moved up.')
                                else:       #evaluation is better, keep move, go to next turbine
                                    flag = 1
                                    turb_power = [i for i in turb_power2]
                                    objective = move2 * 1.
                                    print('turbine ', i, ' moved up.', move2)
                                    #HubHeight_Search(hstep, i, hstepmin, XLocation, YLocation, numturbs, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, rradmin, rradmax, rstep, hubmin, hubmax, rstepmin, aif)
                                
                    if innerflag == 1 and flag == 0:        #move 2 was just unsucessfully attempted
                        transflag = Translation_X(-step2, i, mx)
                        CHECK2 = 0
                        if transflag == 1:      #if the translation moved the turbine out of bounds, go to next translation
                            innerflag = 2       #move3 was attempted
                            print('turbine not left.')
                        else:
                            CHECK2 = Check_Interference(i)
                            if CHECK2 == 1:         #if interference occurs, move the turbine back, go to next translation
                                Translation_X(step2, i, mx)
                                innerflag = 2
                                print('turbine not left.')
                            else:       #if there is no interference, evaluate and store
                                move3, turb_power2 = Eval_Objective(mx,my,ma,A,B,numturbs)

                                if move3 >= nomove:       #if evaluation is worse than initial, move back, go to next translation
                                    Translation_X(step2, i, mx)
                                    innerflag = 2
                                    print('turbine not moved left.')
                                else:       #evaluation is better, keep move, go to next turbine
                                    flag = 1
                                    turb_power = [i for i in turb_power2]
                                    objective = move3 * 1.
                                    print('turbine ', i, ' moved left.', move3)
                                    #HubHeight_Search(hstep, i, hstepmin, XLocation, YLocation, numturbs, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, rradmin, rradmax, rstep, hubmin, hubmax, rstepmin, aif)
                                
                    if innerflag == 2 and flag == 0:        #move 3 was just unsucessfully attempted
                        transflag = Translation_Y(step2, i, my)
                        CHECK2 = 0
                        if transflag == 1:      #if the translation moved the turbine out of bounds, go to next translation
                            innerflag = 3       #move3 was attempted
                            print('turbine not moved down.')
                        else:
                            CHECK2 = Check_Interference(i)
                            if CHECK2 == 1:         #if interference occurs, move the turbine back, go to next translation
                                Translation_Y(-step2, i, my)
                                innerflag = 3
                                print('turbine not moved down.')
                            else:       #if there is no interference, evaluate and store
                                move4, turb_power2 = Eval_Objective(mx,my,ma,A,B,numturbs)

                                if move4 >= nomove:       #if evaluation is worse than initial, move back, go to next translation
                                    Translation_Y(-step2, i, my)
                                    innerflag = 3
                                    print('turbine not moved down.')
                                else:       #evaluation is better, keep move, go to next turbine
                                    flag = 1
                                    turb_power = [i for i in turb_power2]
                                    objective = move4 * 1.
                                    print('turbine ', i, ' moved down.', move4)
                                    #HubHeight_Search(hstep, i, hstepmin, XLocation, YLocation, numturbs, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, rradmin, rradmax, rstep, hubmin, hubmax, rstepmin, aif)
                    if innerflag == 3 and flag == 0:
                        transflag = Translation_X(step2, i, mx)     #move the turbine one step right
                        CHECK2 = 0
                        if transflag == 1:          #if the translation moved the turbine out of bounds, go to next translation
                            innerflag = 4           #signifies move 1 was attempted
                            print('Turbine not moved right.')
                    
                        else:       #if there is the turbine is in bounds, evaluate and store
                            CHECK2 = Check_Interference(i)
                            if CHECK2 == 1:         #if interference occurs, move the turbine back , go to next translation
                                Translation_X(-step2, i, mx)
                                innerflag = 4
                                print('turbine not moved right')
                            else:       #if there is no interference, evaluate and store
                                move1, turb_power2 = Eval_Objective(mx,my,ma,A,B,numturbs)

                                if move1 >= nomove:             #if evaluation is worse than initial, move back, go to next translation
                                    Translation_X(-step2, i, mx)
                                    innerflag = 4
                                    print('Turbine not moved right.')
                                else:
                                    flag = 1           #signifies movement was kept
                                    turb_power = [i for i in turb_power2]
                                    objective = move1 * 1.
                                print('turbine ', i, ' moved right.', move1)
                                #HubHeight_Search(hstep, i, hstepmin, XLocation, YLocation, numturbs, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, rradmin, rradmax, rstep, hubmin, hubmax, rstepmin, aif)
                             
                     
                    if innerflag == 4 and flag == 0:        #translation at this step size has resulted in no moves for this turbine
                        Stopped[i] = 1
                        objective = nomove * 1.
                        #HubHeight_Search(hstep, i, hstepmin, XLocation, YLocation, numturbs, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, rradmin, rradmax, rstep, hubmin, hubmax, rstepmin, aif)
                    
                        #for every turbine in random order
                
                        #check to see if all turbies have not been moved at this step
                        #print_graph()
            exit_css = sum(Stopped)        #exit current step size
            if exit_css == numturbs:
                #all turbines have stopped moving at this step size, halving step size.
                #find worst performing turbine and randomly assign elsewhere
                for b in range(0, num_pops):
                    print("No moves at step size ", step2, " are possible. Popping weakest turbine.")
                    min_power = 5000000.     #initialized to first turbine power output
                    random_vec2 = Rand_Vector(numturbs)    #creates a randomly ordered vector of turbines
                    for j in range(0, numturbs):
                        randorder = random_vec2[j]
                        if turb_power[randorder] < min_power:
                            min_power = turb_power[randorder]
                            min_turb = randorder
                                
                    print('The weakest turbine is turbine ', min_turb, ' with power currently at ', min_power)
                    #start_eval, turb_power = Eval_Objective(mx,my,ma,A,B,numturb)

                    initialx = mx[min_turb]
                    initialy = my[min_turb]
                    checkx = 0
                    k = 0
                    flag = 0
                    while flag == 0 and k < max_pop_tries:
                        checkx = 0
                        while checkx != 1:      #will try random locations until one has no interference
                            CHECK2 = 0
                            mx[min_turb] = random.uniform(0, site_x)
                            my[min_turb] = random.uniform(0, site_y)
                                
                            CHECK2 = Check_Interference(min_turb)
                            if CHECK2 != 1:         #No interference
                                checkx = 1          #place turbine and exit poping loop
                                print('Turbine ', min_turb, ' has moved to a new location.')
                            else:
                                mx[min_turb] = initialx
                                my[min_turb] = initialy
                                print('Turbine cannot be relocated without interference, trying agian.')
                            
                        new_eval, turb_power2 = Eval_Objective(mx,my,ma,A,B,numturbs)
                        if new_eval < objective:
                            flag = 1
                            turb_power = [i for i in turb_power2]
                            #HubHeight_Search(hstep, i, hstepmin, XLocation, YLocation, numturbs, z0, U0, Zref, alphah, ro, yrs, WCOE, condition, rradmin, rradmax, rstep, hubmin, hubmax, rstepmin, aif)
                            print('Move has improved the evaluation. Continuing pattern serach.')
                        else:
                            mx[min_turb] = initialx
                            my[min_turb] = initialy
                            print('Move did not improve evaluation. Trying new moves.')
                        k += 1
                        
                    #print('pops till acceptance= ', k)
                                
                step2 = step2 / 2.0
                print('step sized is being reduced to: ', step2)
                        
            #else:
                #print('Turbines were not all placed at step size ', step2, '. Repeating.')
    return mx, my, mz
###############################################################################
'''things done each time code runs'''
#which inflow angle to plot
alpha = dirs[0]

#domain centered on (0,0)
mesh = RectangleMesh(Point(-Lx/2., -Ly/2.), Point(Lx/2., Ly/2.), nx, ny)

h = mesh.hmin()

'''refine mesh twice in circumradius of farm'''
for nums in range(numRefine):
    # print 'refining mesh'
    mesh=refine_mesh(mesh, site_x, site_y, HH)
    h = mesh.hmin()

'''WHAT'S HAPPENING HERE?! - somehow setting up the mesh to store the values we need?'''
# function spaces, mixed function space syntax not backwards compatible
V = VectorElement('Lagrange', mesh.ufl_cell(), 2) 
Q = FiniteElement('Lagrange', mesh.ufl_cell(), 1)
VQ = FunctionSpace(mesh, MixedElement([V,Q]))   #NSE equations
V = VectorFunctionSpace(mesh, 'Lagrange', 2)
Q = FunctionSpace(mesh, 'Lagrange', 1)



###############################################################################
if __name__ == "__main__":
    
    '''create layout after previous specification of random, gridded, or seeded'''
    mx,my,mz = createLayout(numturbs)
    ### array of 0.33?? Why?
    ma=[Constant(mm) for mm in 0.33*np.ones(numturbs)] ###axial induction factor
    ### calculate the double integral of the smoothing kernel for 
    beta = integrate.dblquad(WTGdist,-3*radius,3*radius,lambda x: -3*radius,lambda x: 3*radius)

    B=beta[0] ### B=estimate of integral, without estimate of error

    if optimize == True: 
        print('optimizing - line 345')
        mx_opt,my_opt,mz_opt = EPS(mx, my, mz, ma)
        '''
        # power functional
        J = Functional(0.) ### what is this??? - data structure type
        for i in range(bins): ###for each wind direction
            tf_rot= createRotatedTurbineForce(mx,my,ma,A,B,numturbs,dirs[i],V) ### calculate force imparted by turbines
            u_rot, p_rot = main(tf_rot)  ### RANS solver
            J = J + weights[i]*rotatedPowerFunctional(dirs[i],A,B,mx,my,ma,u_rot,numturbs,V)

        # position control variables
        m=createControl(mx,my,numturbs)
        bounds = createBounds(mx,my,numturbs)

        rf = ReducedFunctional(J,m)

        def iter_cb(m):
            if MPI.rank(mpi_comm_world()) == 0:
                print("m = ")
                for mm in m:
                    print("Constant("+ str(mm)+ "),")
        ### DOES THIS NOT USE ADJOINT? SEEMS LIKE BFGS -- see last page of journal article
        
        m_opt = maximize(rf, method="L-BFGS-B", options = {"disp": True}, bounds = bounds, callback = iter_cb)
        # m_opt = maximize(rf, method="SLSQP", options = {"disp": True}, bounds = bounds, callback = iter_cb)
        mx_opt,my_opt=splitSolution(m_opt,numturbs)
        mz_opt=mz
        '''

    # otherwise go strait to final plotting routine, plot last inflow direction
    else:
        mx_opt = mx
        my_opt = my
        mz_opt = mz

    tf_rot_opt= createRotatedTurbineForce(mx_opt,my_opt,ma,A,B,numturbs,dirs[0],V)
    tf_rot_opt_out = project(tf_rot_opt, V)
    # plot(tf_rot_opt_out)
    tfStr='tf.pvd' 
    file1 = File(tfStr)
    file1 << tf_rot_opt_out

    u_rot_opt,up_rot_opt  = main(tf_rot_opt) ### RANS analysis
    #print('mx: ',mx)
    #print('my: ',my)
    for i in range(bins):
        
        # report power curve scalar function instead of functional for adjoint derivatives
        if i ==0:
            Jfunc = weights[i]*rotatedPowerFunction(dirs[i],A,B,mx_opt,my_opt,ma,up_rot_opt,numturbs,V)
        else:
            Jfunc = Jfunc + weights[i]*rotatedPowerFunction(dirs[i],A,B,mx_opt,my_opt,ma,up_rot_opt,numturbs,V)
    # power = assemble(Jfunc)
    print('power output: ')
    print(Jfunc)
    data = [Jfunc]
    output = 'off'
    if output == 'on':
        with open('CFD_output.csv', "a+") as csv_file:
            writer = csv.writer(csv_file, delimiter=',')
            for line in data:
                writer.writerow(line)

    u_rot_opt_out=project(u_rot_opt, V)

    uStr='u.pvd'
    file2 = File(uStr)
    file2 << u_rot_opt_out
