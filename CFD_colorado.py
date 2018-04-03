
__author__ = "Annalise Miller <millanna@oregonstate.edu>"
__coauthor__ = "Ryan King <ryan.king@nrel.gov>"
__date__ = "2017-12-1"

from dolfin import *
import numpy as np
import random as random
from scipy import integrate
from scipy import stats
from time import time
from datetime import datetime
import sys
import csv
import matplotlib.pyplot as plt
import os
# import multiprocessing as mp
# mpirun -np 4 python CFD_colorado.py #2>&1 | tee log
plt.switch_backend('agg')

sys.getrecursionlimit()
sys.setrecursionlimit(10000)

set_log_level(ERROR)  # edits what gets printed out
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
air_density = 1.225  # kg/m^3
init_step = 400.  # initial step size
# init_step = 200. #changed for testing
minstep = 3.  # minimum step size
num_pops = 5
max_pop_tries = 1000
# Print log messages only from the root process in parallel
parameters["std_out_all_processes"] = False
# Print log messages only from the root process in parallel
parameters['form_compiler']['cpp_optimize_flags'] = '-O3 -fno-math-errno -march=native' 
parameters["form_compiler"]["optimize"] = True
parameters["form_compiler"]["cpp_optimize"] = True
parameters['form_compiler']['representation'] = 'uflacs'
parameters['form_compiler']['quadrature_degree'] = 12
parameters['linear_algebra_backend'] = "PETSc"
# parameters["mesh_partitioner"] = "SCOTCH"
# solver.parameters.linear_solver = "gmres"

# comm = mpi_comm_world()
# mpiRank = MPI.rank(comm)
# print('test')
# print('MPI PROCESS RANK ', mpiRank)

# actual turbine diam in meters
RD = 80.

# mesh parameters
degree = 2
Lx = 300. * RD  # size of mesh (NOT FARM)
Ly = 300. * RD  # size of mesh (NOT FARM)
numx = 100  # number of points in mesh
numy = 100  # number of points in mesh

# number of inflow direction bins
bins = 1  # Unidirectional wind

WTGexp = 8.  # gamma for smoothing kernel (Eq. 12)
radius = RD/2.
thickness = RD/20.
numRefine = 2
A = RD  # weird for 2D  # WHAT's this actually? - area b/c it's got no z?
HH = 80
initExtent = 0.95  # create buffer on grid layout
mlDenom = 2.  # mixing length denom - not sure how this is arrived at
rad2 = 4. * RD

# Ct = 0.75
# Ct = 8./9. # limit for testing
Ct = 0.861  # GE 1.5sle from Sandia
# Cp = 0.34
Cp = 16./27.  # betz limit for testing
# site/refinement constraints
site_x = 7.*RD  # size of farm
site_y = 7.*RD  # size of farm

adaptive_meshing = True

# only one of these should be true to work
restart = False  # seeded layout (must go into layout function to change seed)
randStart = False  # random layout
gridStart = False  # gridded layout - must have 16 turbines
linearStart = False  # turbines in a row
offgridStart = False
test_turb = False
old_good = False  # previous good layout
coloradoStart = True
if old_good:
    site_x = 30. * RD
    sity_y = 30. * RD
    Lx = 120. * RD
    Ly = 120. * RD
    numturbs = 18
elif offgridStart or gridStart:
    numturbs = 16
elif linearStart:
    numturbs = 10
elif coloradoStart:
    RD = 77.
    radius = RD / 2.
    site_x = 50. * RD
    sity_y = 50. * RD
    Lx = 200. * RD
    Ly = 200. * RD
    numturbs = 37
else:
    numturbs = 10

heat_output = False
print_mesh = True
optimize = False
checkpts = False  # creates plots of turbine movement throughout optimization
'''only for adjoint method'''
# if optimize == True:
#     from dolfin_adjoint import *
tot_evals = 0
# inflow is always from left in simulation


def WTGdist(x, y):  # smoothing kernel
    # actuator disk distribution for normalization constant
    return np.exp(-((x/thickness)**WTGexp + (y/radius)**WTGexp))


def refine_mesh(mesh, site_x, site_y, refine_where, mx, my, mz, ma, rad2):
    # refines the mesh around the site boundaries
    h = mesh.hmin()
    # cell_ct = 0.
    cell_f = CellFunction('bool', mesh, False)
    # create a mesh of the same size, where every cell is set to 'False'
    if refine_where == 'farm':
        for cell in cells(mesh):
            # cell_ct += 1
            if (cell.midpoint()[0]**2 + cell.midpoint()[1]**2
                < site_x**2 + site_y**2 + h):
                cell_f[cell] = True
                # if the midpoint of the cell is within the circumradius of the
                # farm, change the cell's value to true
    else:
        for cell in cells(mesh):  # cycle through each cell
            # cell_ct += 1
            for i in range(numturbs):  # cycle through each turbine
                if (pow(cell.midpoint()[0] - mx[i], 2)
                    + pow(cell.midpoint()[1] - my[i], 2)
                    < pow(rad2, 2) + h):
                    cell_f[cell] = True
                    # if the midpoint of the cell is within a radius of a
                    # turbine
    mesh = refine(mesh, cell_f)
    # refine the cells within the circumradius of the farm
    # print('mesh cell count: ',cell_ct)
    return mesh


def createLayout(numturbs):
    mx = []
    my = []
    mz = []

    if test_turb:
        my = [0., 200., 200.]
        spacing = 300.
        mz = [HH] * numturbs
        mx = [0., -spacing / 2., spacing / 2.]
        # mx = [0.] * numturbs
        # my = [i*spacing for i in range(numturbs)]

    # if you've specified a random starting array
    if randStart:
        for i in range(numturbs):
            # select a random x location, within the farm bounds,
            # and such that the entire rotor will be within the farm bounds
            # mx.append(Constant(np.random.uniform(low=-(site_x - radius),
            #                                      high=(site_x - radius))))
            mx.append(np.random.uniform(low=-(site_x - radius),
                                        high=(site_x - radius)))
            # select a random y location, within the farm bounds,
            # and such that the entire rotor will be within the farm bounds
            # my.append(Constant(np.random.uniform(low=-(site_y - radius),
            #                                      high=(site_y - radius))))
            my.append(np.random.uniform(low=-(site_y - radius),
                                        high=(site_y - radius)))
            # single hub height, z = HH
            # mz.append(Constant(HH))
            mz.append(HH)
    elif gridStart:
        # if you've specified a grided starting array
        # Can only accept 16 turbines
        if numturbs == 16:
            rows = 4
            # rows = 1
            cols = 4
            xpos = np.linspace(-initExtent*(site_x - radius),
                               initExtent*(site_x - radius),
                               cols)
            ypos = np.linspace(-initExtent*(site_y - radius),
                               initExtent*(site_y - radius),
                               rows)
            for i in range(rows):
                for j in range(cols):
                    # mx.append(Constant(xpos[j]))
                    # my.append(Constant(ypos[i]))
                    mx.append(xpos[j])
                    my.append(ypos[i])
                    # some starting noise sometimes helps
                    # mx.append(Constant(xpos[j]+5.*np.random.randn()))
                    # my.append(Constant(ypos[i]+5.*np.random.randn()))
                    # mz.append(Constant(HH))
                    mz.append(HH)

    elif linearStart:
        '''My adds'''
        my = np.linspace(-initExtent*(site_x - radius),
                         initExtent*(site_x - radius),
                         numturbs)
        mx = np.array([0] * numturbs)
        mz = np.array([HH] * numturbs)

    elif offgridStart:
        rows = 4
        cols = 4
        xpos = np.linspace(-initExtent*(site_x - radius),
                           initExtent*(site_x - radius),
                           cols)
        ypos = np.linspace(-initExtent*(site_y - radius),
                           initExtent*(site_y - radius),
                           rows)
        offset = ypos[1] - ypos[0]
        xpos2 = []
        for i in ypos:
            xpos2.append(i + offset/2)
        # print('xlocations')
        # print(xpos)
        for i in range(rows):
            for j in range(cols):
                if i % 2 == 0:
                    mx.append(xpos[j])
                else:
                    mx.append(xpos2[j])
                my.append(ypos[i])
                mz.append(HH)

    elif old_good:
        mx = [1304.2130491853543, 303.16195377176825, 1882.7665062612382,
              734.1365122392898, 36.582152829805295, 1078.6666766839614,
              1779.34556196998, 417.79552464946823, 1986.6903503041713,
              133.48666173166612, 1430.1977444645, 954.8374274350235,
              839.2773010244406, 1186.314920647951, 518.9482768461385,
              1681.3296876680367, 621.3034610992612, 1567.283691418906]
        my = [204.01944774702838, 270.02673301978496, 276.97611783621744,
              305.2706459401636, 1065.3030786448571, 104.67022680498235,
              610.6008734856377, 450.6469581295796, 32.94310546319856,
              1290.6284216353858, 413.28638998844144, 611.4810784612052,
              117.69510957163021, 489.9943326680768, 225.4370943493467,
              858.4028566999443, 505.9650483803546, 1212.01937488221]
        mz = [HH] * numturbs

    elif coloradoStart:
        mx = [2564.600847, 2433.045595, 2263.695463, 2081.262486, 1899.556332,
              1717.850177, -448.0873045, -316.5320365, -680.6712045,
              -886.3625935, 2535.527863, 2339.285222, 2148.857175, 1959.155951,
              1694.591789, 1544.139091, 1389.325444, 1168.370754, 966.3135015,
              761.3489495, 627.6132125, 504.7798445, 618.8913165, 381.2196515,
              -394.3022775, -525.8575455, -653.7786915, -783.8803115,
              -842.0262865, -1033.908007, -1582.660658, -1781.810628,
              -1964.970456, -2204.095787, -2383.621492, -2517.35724,
              -2564.600847]
        my = [3940.192224, 3730.033814, 3534.330744, 3334.179874, 3157.379944,
              2961.676874, 2298.955114, 2565.822934, 2132.162724, 2006.512454,
              1785.234554, 1618.442164, 1458.321464, 1289.305174, 1010.205914,
              862.3166535, 699.9720635, 495.3734035, 363.0514335, 222.9458335,
              5.0037735, -245.1848065, -980.1832765, -1060.243627,
              -959.0562365, -1158.095157, -1352.686277, -1552.837147,
              -1768.555307, -2014.296097, -2744.846767, -2916.086947,
              -3081.76739, -3260.791222, -3415.35217, -3628.846429,
              -3940.192224]
        mz = [HH] * len(mx)

    elif restart:
        # fixed layout here
        m_temp = [Constant(-113.961988283), Constant(-386.535837904),
                  Constant(-512.116113959), Constant(-237.354391531),
                  Constant(638.697968355), Constant(13.6826901448),
                  Constant(386.535838424), Constant(-113.961987466),
                  Constant(13.6826875361), Constant(-638.697971072),
                  Constant(-887.942379804), Constant(-813.542880381),
                  Constant(813.542880031), Constant(-887.942379852),
                  Constant(237.354391629), Constant(-512.116113931),
                  Constant(-237.3543916), Constant(512.116113865),
                  Constant(-813.542880345), Constant(887.942379783),
                  Constant(887.942379753), Constant(813.542880265),
                  Constant(-13.6826884631), Constant(638.697970038),
                  Constant(-386.535837846), Constant(113.961988218),
                  Constant(-638.697970958), Constant(-13.6826879195),
                  Constant(512.116113711), Constant(237.354391612),
                  Constant(113.961988), Constant(386.535838129)]
        # Splits a listinto x and y coordinates, function is shown further down
        mx, my = splitSolution(m_temp, numturbs)
        for i in range(numturbs):
            mz.append(Constant(HH))

    return mx, my, mz


def createRotatedTurbineForce(mx, my, ma, A, beta, numturbs, alpha, V, mesh):
    # mx, my = x, y coordinates of each turbine
    # ma = vector of value 0.33 for each turbine - axial induction factor
    # A = RD at beginning of code. Maybe area since we're using 2D...
    # beta = integral over actuator disk area in x and y: used to normalize
    x = SpatialCoordinate(mesh)
    WTGbase = project(Expression(("1.0", "0.0"), degree=2), V)
    tf = Function(V)
    # print(tf)
    if checkpts:
        check_it = project(tf, V)
        n = ([check_it(cos(alpha)*mx[i] - sin(alpha)*my[i],
                       sin(alpha)*mx[i] + cos(alpha)*my[i])
             for i in range(numturbs)])
        nx = [cos(alpha) * mx[i] - sin(alpha) * my[i] for i in range(numturbs)]
        ny = [sin(alpha) * mx[i] + cos(alpha) * my[i] for i in range(numturbs)]
        fig, ax = plt.subplots()
        ax.scatter(nx, ny)
        for i, txt in enumerate(n):
            ax.annotate(txt, (nx[i], ny[i]))
        plt.savefig('initial_tf.png', bbox_inches='tight')
    for i in range(numturbs):
        # rotation
        xrot = cos(alpha)*mx[i] - sin(alpha)*my[i]
        yrot = sin(alpha)*mx[i] + cos(alpha)*my[i]
        tf = tf + (0.5*A*Ct/((1.-ma[i]) ** 2)/beta
                   * exp(-(((x[0] - xrot)/thickness)**WTGexp
                           + ((x[1] - yrot)/radius)**WTGexp))
                   * WTGbase.copy(deepcopy=True))
    if checkpts:
        check_it = project(tf, V)
        n = ([check_it(cos(alpha)*mx[i] - sin(alpha)*my[i],
                       sin(alpha)*mx[i] + cos(alpha)*my[i])
             for i in range(numturbs)])
        fig, ax = plt.subplots()
        ax.scatter(nx, ny)
        for i, txt in enumerate(n):
            ax.annotate(txt, (nx[i], ny[i]))
        plt.savefig('final_tf.png', bbox_inches='tight')
    return tf


def rotatedPowerFunctional(alpha, A, beta, mx, my, ma, u, numturbs, V):
    # used only in optimization
    # functional for dolfin-adjoint
    x = SpatialCoordinate(mesh)
    J = Functional(0.)
    for i in range(numturbs):
        # rotation
        xrot = cos(alpha)*mx[i] - sin(alpha)*my[i]
        yrot = sin(alpha)*mx[i] + cos(alpha)*my[i]

        J = J + (0.5 * np.pi * radius**2 * Cp / ((1. - ma[i]) ** 3) / beta
                 * Functional(
                         exp(-(((x[0] - xrot) / thickness)**WTGexp
                               + ((x[1] - yrot)/radius)**WTGexp))
                         * u[0]**3. * dx))
    # /beta added by Annalise 12/13/17
    return J


def rotatedPowerFunction(alpha, A, beta, mx, my, ma, up, numturbs,
                         V, mesh, heat=False):
    # used in developing final power
    # emulating an actual power curve
    x = SpatialCoordinate(mesh)
    J = []
    if checkpts:
        nx = [cos(alpha)*mx[i] - sin(alpha)*my[i] for i in range(numturbs)]
        ny = [sin(alpha)*mx[i] + cos(alpha)*my[i] for i in range(numturbs)]
        fig, ax = plt.subplots()
        ax.scatter(nx, ny)
        n = [up.sub(0)(nx[i], ny[i])[0] for i in range(numturbs)]
        for i, txt in enumerate(n):
            ax.annotate(txt, (nx[i], ny[i]))
        plt.savefig('windspeeds.png', bbox_inches='tight')
    for i in range(numturbs):
        # rotation
        xrot = cos(alpha) * mx[i] - sin(alpha) * my[i]
        # -5 added by Annalise 12/14 to understand effects of
        # smoothing kernal on usturbine
        yrot = sin(alpha)*mx[i] + cos(alpha)*my[i]

        # print(up.sub(0)(xrot, yrot)[0])

        J.append(up.sub(0)(xrot, yrot)[0])
        # except:
        #     print(xrot, yrot)
        # up.sub(0)(xrot,yrot)[0] --> up == u
        # p combined --> sub(0) == just u --> (xrot, yrot)
        # == position of interest (center pt) --> [0] == x-velocity
    if heat:
        heat_out = [[]]
        outvals = 500
        nx = [cos(alpha) * mx[i] - sin(alpha) * my[i] for i in range(numturbs)]
        ny = [sin(alpha) * mx[i] + cos(alpha) * my[i] for i in range(numturbs)]
        interval_x = (max(nx) - min(nx)) * 2. / outvals
        if interval_x > 0.01:
            x_start = min(nx) - (max(nx) - min(nx)) * 0.5 + interval_x / 2.
            x1 = min(nx) - (max(nx) - min(nx)) * 0.5
            x2 = max(nx) + (max(nx) - min(nx)) * 0.5
        else:
            x_start = min(nx) - 100. + 200. / (2. * outvals)
            interval_x = 200. / outvals
            x1 = min(nx) - 100.
            x2 = max(nx) + 100.
        interval_y = (max(ny) - min(ny)) * 2. / outvals
        if interval_y > 0.01:
            y_start = min(ny) - (max(ny) - min(ny)) * 0.5 + interval_y / 2.
            y1 = min(ny) - (max(ny) - min(ny)) * 0.5
            y2 = max(ny) + (max(ny) - min(ny)) * 0.5
        else:
            y_start = min(ny) - 100. + 200. / (2. * outvals)
            interval_y = 200. / outvals
            y1 = min(ny) - 100.
            y2 = max(ny) + 100.
        spacing_outx = [i * interval_x + x_start for i in range(outvals)]
        spacing_outy_sub = [i * interval_y + y_start for i in range(outvals)]
        spacing_outy = ([spacing_outy_sub[-i]
                        for i in range(1, len(spacing_outy_sub) + 1)])
        for j in spacing_outy:
            heat_out[0].append([up.sub(0)(i, j)[0] for i in spacing_outx])
        heat_out.append([x1, x2, y1, y2])
        # print(heat_out)
        return J, heat_out
    else:
        return J


def createControl(mx, my, numturbs):
        m = []
        for i in range(numturbs):
            m.append(Control(mx[i]))
            m.append(Control(my[i]))

        return m


def createBounds(mx, my, numturbs):
        ub = []
        lb = []
        for i in range(numturbs):
            lb.append(Constant(-(site_x - radius)))
            lb.append(Constant(-(site_y - radius)))
            ub.append(Constant((site_x - radius)))
            ub.append(Constant((site_y - radius)))

        bounds = [lb, ub]
        return bounds


def createAxialBounds(ma, numturbs):
        ub = []
        lb = []
        for i in range(numturbs):
            lb.append(Constant(0.))
            ub.append(Constant(0.75))

        bounds = [lb, ub]
        return bounds


def splitSolution(m_opt, numturbs):
    mx_opt = []
    my_opt = []
    j = 0
    for i in range(numturbs):
        mx_opt.append(m_opt[j])
        j += 1
        my_opt.append(m_opt[j])
        j += 1

    return mx_opt, my_opt


def main(tf, wind_case, VQ):
    nu = Constant(.00005)
    # kinematic viscosity
    f = Constant((0., 0.))
    up_next = Function(VQ)
    # up_next becomes tuple of vector and finite elements for wind speed
    # and pressure
    u_next, p_next = split(up_next)
    # split vector (wind speed) and finite (pressure) elements?
    v, q = TestFunctions(VQ)

    class InitialConditions(Expression):
        # inherits from Expression class in fenics
        def __init__(self, **kwargs):
            random.seed(2)

        def eval(self, values, x):
            values[0] = wind_cases[wind_case][1]
            values[1] = 0.0
            values[2] = 0.0

        def value_shape(self):
            return (3,)
    # boundary conditions

    class NoSlipBoundary(SubDomain):
        def inside(self, x, on_boundary):
            # windspeed has w = 0 at top and bottom
            return near(x[1]**2 - (Ly/2.)**2, 0.) and on_boundary

    class InflowBoundary(SubDomain):
        # windspeeed has u = inflow velocity and w = 0 at locations of inflow
        def inside(self, x, on_boundary):
            return near(x[0], -(Lx/2.)) and on_boundary

    class PeriodicBoundary(SubDomain):

        def inside(self, x, on_boundary):
            # return True if on left or bottom boundary
            # AND NOT on one of the two slave edges
            return bool((near(x[0], -(Lx/2.)) or near(x[1], -(Ly/2.)))
                        and (not (near(x[0], (Lx/2.)) or near(x[1], (Ly/2.))))
                        and on_boundary)

        def map(self, x, y):
            if near(x[0], (Lx/2.)) and near(x[1], (Ly/2.)):
                y[0] = x[0] - 2*(Lx/2.)
                y[1] = x[1] - 2*(Ly/2.)
            elif near(x[0], (Lx/2.)):
                y[0] = x[0] - 2*(Lx/2.)
                y[1] = x[1]
            else:  # near(x[1], (Ly/2.)):
                y[0] = x[0]
                y[1] = x[1] - 2*(Ly/2.)

    u0 = InitialConditions(degree=2)
    up_next.interpolate(u0)

    lmix = radius/mlDenom  # mixing lenth
    # mean rate of strain tensor ... so confused
    S = sqrt(2.*inner(0.5*(grad(u_next) + grad(u_next).T),
                      0.5*(grad(u_next)+grad(u_next).T)))
    nu_T = lmix**2. * S  # eddie viscosity

    F = (inner(grad(u_next) * u_next, v) * dx
         + (nu+nu_T) * inner(grad(u_next), grad(v)) * dx
         - inner(div(v), p_next) * dx
         - inner(div(u_next), q) * dx
         - inner(f, v) * dx
         + inner(tf * u_next[0]**2, v) * dx)
    # u_next[0] ** 2 because it's not added to the force on f_AD calc?
    # I think this is ok without density
    # --> http://libelemental.org/featured/2011/pdesys.html

    # lateral BC
    # bc1 = DirichletBC(VQ.sub(0), Constant((8.0,0.0)), NoSlipBoundary())
    bc1a = DirichletBC(VQ.sub(0).sub(1), Constant(0.0), NoSlipBoundary())

    # inflow BC
    bc2 = DirichletBC(VQ.sub(0), Constant((wind_cases[wind_case][1], 0.0)),
                      InflowBoundary())
    # bc2a = DirichletBC(VQ.sub(0).sub(0), Constant(8.), InflowBoundary())

    # outflow pressure BC is implicitly 0

    bc = [bc1a, bc2]
    J = derivative(F, up_next)
    problem = NonlinearVariationalProblem(F, up_next, bc, J=J)
    solver = NonlinearVariationalSolver(problem)
    prm = solver.parameters
    list_linear_solver_methods()
    solver.nonlinear_variational_solver = 'newton_solver'
    prm["newton_solver"]["absolute_tolerance"] = 1E-3  # 1E-8
    prm["newton_solver"]["relative_tolerance"] = 1E-3  # 1E-7
    prm["newton_solver"]["maximum_iterations"] = 120
    prm["newton_solver"]["relaxation_parameter"] = 1.0
    prm["newton_solver"]["linear_solver"] = 'mumps'
    # Added for server
#    prm["newton_solver"]["krylov_solver"]["relative_tolerance"]=1e-3
#    prm["newton_solver"]["krylov_solver"]["maximum_iterations"]=1000
#    prm["newton_solver"]["linear_solver"] = 'gmres'
#    prm["newton_solver"]["preconditioner"] = 'ilu'
#    solver.solve()

    try:
        solver.solve()
    except:
        print('solver error')
        return 'fuck', 'fuck'
    # solve(F == 0, up_next, bc,
    #       solver_parameters={"newton_solver":{"absolute_tolerance": 1e-8},
    #       "newton_solver":{"maximum_iterations": 120}})
    u_next, p_next = split(up_next)
    # print('u_next: ', up_next.sub(0)(-Lx/2 + 0.000005,-Ly/2 + 0.000005)[0])
    # if optimize == False:
    #     nu_T_out=project(nu_T, Q)
    #     lStr= 'nu_t.pvd'
    #     file = File(lStr)
    #     file << nu_T_out
    return u_next, up_next
###############################################################################


def Eval_Objective(mx, my, ma, A, B, numturb, heat=False):
    global tot_evals
    tot_evals += 1
    start = time()
    J = 0.
    cumulative_ws = [0.] * numturbs
    for i in range(bins * len(windsp)):  # for each wind direction
        V, Q, VQ, mesh = create_mesh(mx, my, mz, ma, rad2)
        tf_rot = createRotatedTurbineForce(mx, my, ma, A, B,
                                           numturbs, wind_cases[i][0], V, mesh)
        # calculate force imparted by turbines
        u_rot, up_rot = main(tf_rot, i, VQ)  # RANS solver
        if u_rot == 'fuck':
            if heat:
                return 'error', ['error'] * numturbs, 'error'
            else:
                print('error')
                return 'error', ['error']*numturbs
        if heat and i == 0:  # only calc heat for wind from left
            power_dev, heat_out = rotatedPowerFunction(wind_cases[i][0], A,
                                                       beta, mx, my, ma,
                                                       up_rot, numturbs, V,
                                                       mesh, True)
        else:
            wind_sp = rotatedPowerFunction(wind_cases[i][0], A, beta, mx, my,
                                           ma, up_rot, numturbs, V, mesh)
        J = J - weights[i]*sum(wind_sp)
        cumulative_ws = ([k*weights[i] + j
                          for k, j in zip(wind_sp, cumulative_ws)])
        print('time to evaluate bin: ', (time() - start)/60.)
        start = time()
    print('new evaluation: ', J)
    if heat:
        return J, cumulative_ws, heat_out
    else:
        return J, cumulative_ws
###############################################################################


def Rand_Vector(numturbs):
    random_vec = []
    for i in range(0, numturbs):
        random_vec.append(i)
    # shuffle elements by randomly exchanging each with one other
    for i in range(0, len(random_vec)):
        r = random.randint(0, len(random_vec)-1)
        # select random value in vector
        temp = random_vec[i]
        random_vec[i] = random_vec[r]
        random_vec[r] = temp
    return random_vec
###############################################################################


def Translation_X(step_size, i, mx):
    transflag = 0
    # print(xstart, step_size)
    xfinish = mx[i] + step_size
    # find preliminary new x-location given step size
    if xfinish >= -site_x and xfinish <= site_x:
        # if this new x-location is not out of bounds, translate it
        mx[i] = xfinish
        return mx, transflag
    else:
        transflag = 1
        return mx, transflag
###############################################################################


def Translation_Y(step_size, i, my):
    transflag = 0
    yfinish = my[i] + step_size
    # find preliminary new x-location given step size
    if yfinish >= -site_y and yfinish <= site_y:
        # if this new x-location is not out of bounds, translate it
        my[i] = yfinish
        return my, transflag
    else:
        transflag = 1
        return my, transflag
###############################################################################


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
    # clear last layouts file
    with open('layouts.txt', 'w') as layouts_file:
        layouts_file.close()

    Stopped = [0.] * numturbs
    # print('type mx: ',type(mx))
    # print('type my: ',type(my))
    enum = 1
    while os.path.isdir('trial_'+str(enum)):
        enum += 1
    directory = 'trial_'+str(enum)
    objective, turb_power = Eval_Objective(mx, my, ma, A, B, numturbs)
    # objective = nomove * 1.
    for h in range(0, 1):
        step2 = init_step
        while step2 >= minstep:
            random_vec = Rand_Vector(numturbs)
            # creates a randomly ordered vector of turbines
            for j in range(0, len(random_vec)):
                i = random_vec[j]
                # print('type(mx[i]): ',type(mx[i]))
                # print('mx[i]: ',mx[i])
                Stopped[i] == 0
                print('Turbine ', i, ' is being tested.')
                flag = 0
                innerflag = 0
                transflag = 0

                if innerflag == 0 and flag == 0:
                    # move 1 was just unsucessfully attempted
                    # print('initial my: ', my)
                    my, transflag = Translation_Y(step2, i, my)
                    # print('new my: ', my)
                    CHECK2 = 0
                    if transflag == 1:
                        # if the translation moved the turbine out of bounds,
                        # go to next translation
                        innerflag = 1  # move2 was attempted
                        print('turbine not moved up. Outside bounds')
                    else:
                        CHECK2 = Check_Interference(mx, my, i)
                        if CHECK2 == 1:
                            # if interference occurs, move the turbine back,
                            # go to next translation
                            my, transflag = Translation_Y(-step2, i, my)
                            innerflag = 1
                            print('turbine not moved up. Interference')
                        elif check_layout(mx, my):
                            my, transflag = Translation_Y(-step2, i, my)
                            innerflag = 1
                            print('turbine not moved up. ' +
                                  'Layout already attempted')
                        else:
                            # if there is no interference, evaluate and store
                            move2, turb_power2 = Eval_Objective(mx, my, ma,
                                                                A, B, numturbs)

                            if move2 >= objective:
                                # if evaluation is worse than initial,
                                # move back, go to next translation
                                my, transflag = Translation_Y(-step2, i, my)
                                innerflag = 1
                                print('turbine not moved up. Worse Evaluation')
                            else:
                                # evaluation is better, keep move,
                                # go to next turbine
                                flag = 1
                                turb_power = [w for w in turb_power2]
                                objective = move2 * 1.
                                print('turbine ', i, ' moved up.', move2)
                                print(my)
                                print_graph(mx, my, i, directory)

                if innerflag == 1 and flag == 0:
                    # move 2 was just unsucessfully attempted
                    mx, transflag = Translation_X(-step2, i, mx)
                    CHECK2 = 0
                    if transflag == 1:
                        # if the translation moved the turbine out of bounds,
                        # go to next translation
                        innerflag = 2  # move3 was attempted
                        print('turbine not left. Outside Bounds')
                    else:
                        CHECK2 = Check_Interference(mx, my, i)
                        if CHECK2 == 1:
                            # if interference occurs, move the turbine back,
                            # go to next translation
                            mx, transflag = Translation_X(step2, i, mx)
                            innerflag = 2
                            print('turbine not left. Interference')
                        elif check_layout(mx, my):
                            mx, transflag = Translation_X(step2, i, my)
                            innerflag = 2
                            print('turbine not moved left. ' +
                                  'Layout already attempted')
                        else:
                            # if there is no interference, evaluate and store
                            move3, turb_power2 = Eval_Objective(mx, my, ma,
                                                                A, B, numturbs)

                            if move3 >= objective:
                                # if evaluation is worse than initial,
                                # move back, go to next translation
                                mx, transflag = Translation_X(step2, i, mx)
                                innerflag = 2
                                print('turbine not moved left. Worse Eval')
                            else:
                                # evaluation is better, keep move,
                                # go to next turbine
                                flag = 1
                                turb_power = [w for w in turb_power2]
                                objective = move3 * 1.
                                print('turbine ', i, ' moved left.', move3)
                                print(mx)
                                print_graph(mx, my, i, directory)

                if innerflag == 2 and flag == 0:
                    # move 3 was just unsucessfully attempted
                    my, transflag = Translation_Y(-step2, i, my)
                    CHECK2 = 0
                    if transflag == 1:
                        # if the translation moved the turbine out of bounds,
                        # go to next translation
                        innerflag = 3  # move3 was attempted
                        print('turbine not moved down. Outside Bounds')
                    else:
                        CHECK2 = Check_Interference(mx, my, i)
                        if CHECK2 == 1:
                            # if interference occurs, move the turbine back,
                            # go to next translation
                            my, transflag = Translation_Y(step2, i, my)
                            innerflag = 3
                            print('turbine not moved down. Interference')
                        elif check_layout(mx, my):
                            my, transflag = Translation_Y(step2, i, my)
                            innerflag = 3
                            print('turbine not moved down.' +
                                  'Layout already attempted')
                        else:
                            # if there is no interference, evaluate and store
                            move4, turb_power2 = Eval_Objective(mx, my, ma,
                                                                A, B, numturbs)

                            if move4 >= objective:
                                # if evaluation is worse than initial,
                                # move back, go to next translation
                                my, transflag = Translation_Y(step2, i, my)
                                innerflag = 3
                                print('turbine not moved down.' +
                                      'Worse Evaluation')
                            else:
                                # evaluation is better, keep move,
                                # go to next turbine
                                flag = 1
                                turb_power = [w for w in turb_power2]
                                objective = move4 * 1.
                                print('turbine ', i, ' moved down.', move4)
                                print(my)
                                print_graph(mx, my, i, directory)

                if innerflag == 3 and flag == 0:
                    mx, transflag = Translation_X(step2, i, mx)
                    # move the turbine one step right
                    CHECK2 = 0
                    if transflag == 1:
                        # if the translation moved the turbine out of bounds,
                        # go to next translation
                        innerflag = 4  # signifies move 1 was attempted
                        print('Turbine not moved right. Outside Bounds')

                    else:
                        # if there is the turbine is in bounds,
                        # evaluate and store
                        CHECK2 = Check_Interference(mx, my, i)
                        if CHECK2 == 1:
                            # if interference occurs, move the turbine back,
                            # go to next translation
                            mx, transflag = Translation_X(-step2, i, mx)
                            innerflag = 4
                            print('turbine not moved right. Interference')
                        elif check_layout(mx, my):
                            mx, transflag = Translation_X(-step2, i, my)
                            innerflag = 4
                            print('turbine not moved left.' +
                                  'Layout already attempted')
                        else:
                            # if there is no interference, evaluate and store
                            move1, turb_power2 = Eval_Objective(mx, my, ma,
                                                                A, B, numturbs)

                            if move1 >= objective:
                                # if evaluation is worse than initial,
                                # move back, go to next translation
                                mx, transflag = Translation_X(-step2, i, mx)
                                innerflag = 4
                                print('Turbine not moved right.' +
                                      'Worse Evaluation')
                            else:
                                flag = 1  # signifies movement was kept
                                turb_power = [w for w in turb_power2]
                                objective = move1 * 1.
                            print('turbine ', i, ' moved right.', move1)
                            print(mx)
                            print_graph(mx, my, i, directory)

                if innerflag == 4 and flag == 0:
                    # translation at this step size has resulted in no moves
                    # for this turbine
                    Stopped[i] = 1

            exit_css = sum(Stopped)  # exit current step size
            if exit_css == numturbs:
                # all turbines have stopped moving at this step size,
                # halving step size.
                # find worst performing turbine and randomly assign elsewhere
                for b in range(0, num_pops):
                    print("No moves at step size " + str(step2) +
                          " are possible. Popping weakest turbine.")
                    min_power = 5000000.
                    # initialized to first turbine power output
                    random_vec2 = Rand_Vector(numturbs)
                    # creates a randomly ordered vector of turbines
                    for j in range(0, numturbs):
                        randorder = random_vec2[j]
                        if turb_power[randorder] < min_power:
                            min_power = turb_power[randorder]
                            min_turb = randorder

                    print('The weakest turbine is turbine ' + str(min_turb) +
                          ' with power currently at ' + str(min_power))
                    # start_eval, turb_power = Eval_Objective(mx, my, ma,
                    #                                         A,B,numturb)

                    initialx = mx[min_turb]
                    initialy = my[min_turb]
                    checkx = 0
                    k = 0
                    flag = 0
                    while flag == 0 and k < max_pop_tries:
                        checkx = 0
                        while checkx != 1:
                            # will try random locations until
                            # one has no interference
                            CHECK2 = 0
                            mx[min_turb] = random.uniform(0, site_x)
                            my[min_turb] = random.uniform(0, site_y)

                            CHECK2 = Check_Interference(mx, my, min_turb)
                            if CHECK2 != 1:  # No interference
                                checkx = 1
                                # place turbine and exit poping loop
                                print('Turbine ' + str(min_turb) +
                                      ' has moved to a new location.')
                            else:
                                mx[min_turb] = initialx
                                my[min_turb] = initialy
                                print('Turbine cannot be relocated without '
                                      + 'interference, trying agian.')

                        new_eval, turb_power2 = Eval_Objective(mx, my, ma,
                                                               A, B, numturbs)
                        if new_eval < objective:
                            flag = 1
                            turb_power = [w for w in turb_power2]
                            objective = new_eval * 1.
                            # HubHeight_Search(hstep, i, hstepmin, XLocation,
                            #                  YLocation, numturbs, z0, U0,
                            #                  Zref, alphah, ro, yrs, WCOE,
                            #                  condition, rradmin, rradmax,
                            #                  rstep, hubmin, hubmax, rstepmin,
                            #                  aif)
                            print('Move has improved the evaluation.' +
                                  'Continuing pattern serach.')
                            print_graph(mx, my, min_turb, directory)
                        else:
                            mx[min_turb] = initialx
                            my[min_turb] = initialy
                            print('Move did not improve evaluation.' +
                                  'Trying new moves.')
                        k += 1

                    # print('pops till acceptance= ', k)

                step2 = step2 / 2.0
                print('step sized is being reduced to: ', step2)

            # else:
                # print('Turbines were not all placed at step size '
                #       + str(step2) + '. Repeating.')
    return mx, my, mz
###############################################################################


def check_layout(mx, my):
    coords = [(i, j) for i, j in zip(mx, my)]  # zip into tuple
    has_been = True
    check = 0
    num = 0
    with open('layouts.txt', 'a+') as layouts_file:
        # print(sum(1 for i in layouts_file))
        for line in layouts_file:
            num += 1
            print(line)
            line = line.strip('[]')
            # line = line.split(',')
            # print(line)
            old_coords = [float(x.strip()) for x in line.split(',')]
            old_coords = ([(old_coords[i], old_coords[i + numturbs])
                          for i in range(numturbs)])
            for i in coords:
                if i not in old_coords:
                    check += 1
                    break
        if abs(check - num) < 0.00005:
            has_been = False
        if not has_been:
            new_layout = mx + my
            layouts_file.write(str(new_layout))
        else:
            print('layout previously checked')
        layouts_file.close()
    return has_been
###############################################################################


def print_graph(all_x, all_y, i, directory):
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.scatter(all_x, all_y, s=30, c='k', marker=(3, 2))
    ax1.scatter(all_x[i], all_y[i], s=30, c='r', marker=(3, 2))

    plt.ylabel('Position (m)')
    plt.xlabel('Position (m)')
    plt.title(str('Optimization of ' + str(numturbs) + ' Turbines'))
    if not os.path.exists(directory):
        os.makedirs(directory)
    plt.savefig(directory + str('/step_turb_' + str(i) + '.png'),
                bbox_inches='tight')
###############################################################################


'''things done each time code runs'''
# which inflow angle to plot
eval_num = 0
# domain centered on (0,0)


def create_mesh(mx, my, mz, ma, rad2):
    # print('Ly ',Ly)
    # print('Lx ', Lx)
    mesh = RectangleMesh(Point(-Lx/2., -Ly/2.),
                         Point(Lx/2., Ly/2.), numx, numy)
    dummy = mesh.coordinates()
    # print('min x: ',min([i[0] for i in dummy]))
    # print('max x: ',max([i[0] for i in dummy]))
    # print('min y: ',min([i[1] for i in dummy]))
    # print('max y: ',max([i[1] for i in dummy]))
    h = mesh.hmin()
    '''refine mesh twice in circumradius of farm'''
    for nums in range(numRefine):
        # print 'refining mesh'
        mesh = refine_mesh(mesh, site_x, site_y,
                           'farm', mx, my, mz, ma, rad2)
        h = mesh.hmin()
    if print_mesh:
        mesh1 = []
        for each in mesh.coordinates():
            if abs(each[0]) < site_x and abs(each[1]) < site_y:
                mesh1.append((float(each[0]), float(each[1])))
    if adaptive_meshing:
        mesh = refine_mesh(mesh, site_x, site_y,
                           'turbines', mx, my, mz, ma, rad2)
        h = mesh.hmin()
    if print_mesh:
        meshx2 = []
        meshy2 = []
        for each in mesh.coordinates():
            if (abs(each[0]) < site_x
                and abs(each[1]) < site_y
                and ((each[0], each[1])
                not in mesh1)):

                meshx2.append(float(each[0]))
                meshy2.append(float(each[1]))
        plt.figure()
        # print(meshx2)
        plt.scatter([iii[0] for iii in mesh1],
                    [iii[1] for iii in mesh1],
                    s=1,
                    c='k')
        plt.scatter(meshx2, meshy2, s=1, c='r')
        # plt.scatter(mx,my,color = 'r', marker='*')
        plt.savefig('mesh_vis_colorado.png')
    print('mesh size: ', len(mesh.coordinates()))
    h = mesh.hmin()

    '''WHAT'S HAPPENING HERE?!
    - somehow setting up the mesh to store the values we need?'''
    # function spaces, mixed function space syntax not backwards compatible
    V = VectorElement('Lagrange', mesh.ufl_cell(), 2)
    Q = FiniteElement('Lagrange', mesh.ufl_cell(), 1)
    VQ = FunctionSpace(mesh, MixedElement([V, Q]))  # NSE equations
    V = VectorFunctionSpace(mesh, 'Lagrange', 2)
    Q = FunctionSpace(mesh, 'Lagrange', 1)
    return V, Q, VQ, mesh


###############################################################################
if __name__ == "__main__":

    '''create layout after previous specification of random,
    gridded, or seeded'''
#    processes = []
#    print('cores available: ',mp.cpu_count())
#    cores = 2.
#    print('starting multi-core tests on', cores, 'cores')
#    print('running...')
#    for core in range(cores):
#            processes.append(mp.Process(target=refine_waypoints,
#                                        args=(generate_waypoints(num_points),
#                                              point_step,
#                                              instance,
#                                              core,
#                                              score_limit,
#                                              num_runs)))

    mx, my, mz = createLayout(numturbs)
    # array of 0.33?? Why?
    aif = 0.314
    ma = [Constant(mm) for mm in aif * np.ones(numturbs)]
    # axial induction factor for GE 1.5sle
    # calculate the double integral of the smoothing kernel for
    beta = integrate.dblquad(WTGdist,
                             -3*radius, 3*radius,
                             lambda x: -3*radius,
                             lambda x: 3*radius)

    B = beta[0]  # B=estimate of integral, without estimate of error
    # otherwise go strait to final plotting routine, plot last inflow direction
    mx_opt = [i for i in mx]
    my_opt = [i for i in my]
    mz_opt = [i for i in mz]
    poss_directions = [270, 280, 290, 300, 310, 320, 330, 340, 350]
    poss_ws = [0., 5., 10., 15., 20., 25.]
    poss_ust = [0, 1, 2, 3, 4, 5]
    error = [[0.] * numturbs for i in range(len(poss_directions))]
    error_ct = [[0.] * numturbs for i in range(len(poss_directions))]
    error2 = [[0.] * numturbs for i in range(len(poss_ws))]
    error2_ct = [[0.] * numturbs for i in range(len(poss_ws))]
    error3 = [[0.] * numturbs for i in range(len(poss_ust))]
    error3_ct = [[0.] * numturbs for i in range(len(poss_ust))]
    RSME = [[0.] * numturbs for i in range(len(poss_directions))]
    # by wind direction
    RSME2 = [[0.] * numturbs for i in range(len(poss_ws))]
    # by wind speed
    RSME3 = [[0.] * numturbs for i in range(len(poss_ust))]
    # by number of upstream turbines
    num_usturbines = ([0] * 10 + [2, 3, 2, 2, 2, 1]
                      + [0] * 6 + [3, 1, 2, 3, 3, 2, 3, 1] + [0] * 6 + [1])
    with open('colorado_wind_data_CLEANED.csv', newline='') as csvfile:
        info = csv.reader(csvfile, delimiter=',', quotechar='|')
        enum = 1
        while os.path.isfile('CFD_colorado_results/RMSE_CFD_colorado_ws' +
                             str(enum) + '.csv'):
            enum += 1
        with open('CFD_colorado_results/RMSE_CFD_colorado_ws' +
                  str(enum) + '.csv', 'w+', newline='') as outfile:
            full_write = csv.writer(outfile)
            # num_lines = 0.
            counter = 0
            beggining = 144 + (enum - 1 - 24) * 4
            ending = 144 + (enum - 24) * 4
            for i, row in enumerate(info):  # subset for testing!
                # if counter >= 10:
                #     break
                # counter += 1
                if i >= beggining and i < ending:
                    print('test number: ' + str(i))
                    windsp = [float(row[6])]  # one wind speed
                    dirs = [float(row[5]) / 180. * np.pi - 3*np.pi/2]
                    # dirs = np.linspace(0, 2*pi, bins, endpoint = False)
                    weights = np.ones(bins * len(windsp))/(bins * len(windsp))
                    # assume even weightings
                    wind_cases = []
                    for i in dirs:
                        for j in windsp:
                            wind_cases.append((i, j))
                            # tuple of wind speed + wind diretion
                    print('wind_cases: ', wind_cases)
                    if heat_output:
                        Jfunc, cum_power, heat_out = Eval_Objective(mx_opt,
                                                                    my_opt, ma,
                                                                    A, B,
                                                                    numturbs,
                                                                    True)
                        if heat_out != 'error':
                            plt.figure()
                            plt.imshow(heat_out[0], cmap='hot',
                                       interpolation='nearest',
                                       extent=heat_out[1])
                            plt.colorbar()
                            nx = ([cos(dirs[0])*mx[i] - sin(dirs[0])*my[i]
                                  for i in range(numturbs)])
                            ny = ([sin(dirs[0])*mx[i] + cos(dirs[0])*my[i]
                                  for i in range(numturbs)])
                            for i, j in zip(nx, ny):
                                plt.plot([i, i],
                                         [j - radius, j + radius],
                                         '-k')
                            now = datetime.now()
                            plt.title('Cross-Field Wind Speed')
                            plt.xlabel('Crosswind Distance (m)')
                            plt.ylabel('Downwind Distance (m)')
                            plt.savefig('heat_out_'+str(now.month)+'_'
                                        + str(now.day)+'_' + str(now.hour)
                                        + '_' + str(now.minute) + '.png',
                                        bbox_inches='tight')
                    else:
                        Jfunc, cum_ws = Eval_Objective(mx_opt, my_opt,
                                                       ma, A, B, numturbs)

                    turbines_realsp = [float(temp) for temp in row[7:]]
                    # windspeeds measured from nacelles
                    output_row = [float(this) for this in row] + cum_ws
                    full_write.writerow(output_row)
                    # print(row[5])
                    rmse = sum([pow(turbines_realsp[iiii] - cum_ws[iiii],
                                    2) for iiii in range(numturbs)])/numturbs
                    rmse = np.sqrt(rmse)
                    print('RMSE: ' + str(rmse))
                    direction_coord = poss_directions.index(int(row[5]))
    output = 'off'
    '''
    if output == 'on':
        tot_dir = [0.] * len(error)
        tot_ws = [0.] * len(error2)
        tot_ust = [0.] * len(error3)
        for i, j in enumerate(error):
            try:
                tot_dir[i] = np.sqrt(sum(j) / sum(error_ct[i]))
            except:
                tot_dir[i] = 'N/A'
            for a, b in enumerate(j):
                try:
                    RSME[i][a] = np.sqrt(b / error_ct[i][a])
                except:
                    RSME[i][a] = 'N/A'
        for i, j in enumerate(error2):
            try:
                tot_ws[i] = np.sqrt(sum(j) / sum(error2_ct[i]))
            except:
                tot_ws[i] = 'N/A'
            for a, b in enumerate(j):
                try:
                    RSME2[i][a] = np.sqrt(b / error2_ct[i][a])
                except:
                    RSME2[i][a] = 'N/A'
        for i, j in enumerate(error3):
            try:
                tot_ust[i] = np.sqrt(sum(j) / sum(error3_ct[i]))
            except:
                tot_ust[i] = 'N/A'
            for a, b in enumerate(j):
                try:
                    RSME3[i][a] = np.sqrt(b / error3_ct[i][a])
                except:
                    RSME3[i][a] = 'N/A'
        with open('RMSE_CFD_colorado_byturb_rotor.csv',
                  'w+', newline='') as outfile:
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
        with open('RMSE_CFD_colorado_rotor.csv',
                  'w+', newline='') as outfile:
            write = csv.writer(outfile)
            write. writerow(poss_directions)
            write.writerow(tot_dir)
            write.writerow([])
            write. writerow(poss_ws)
            write.writerow(tot_ws)
            write.writerow([])
            write. writerow(poss_ust)
            write.writerow(tot_ust)
            # print('The final layout with '
            #       + str(initial_num) + ' turbines has a score of: '
            #       + str(score))
    '''
    print('done!')
    '''
    u_rot_opt_out=project(u_rot_opt, V)

    uStr='u.pvd'
    file2 = File(uStr)
    file2 << u_rot_opt_out
    '''
