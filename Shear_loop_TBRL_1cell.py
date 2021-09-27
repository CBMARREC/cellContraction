__author__ = 'Georgios Grekas (grekas.g@gmail.com)'
results_path = 'results/'

from init_parameters import *
# import linecache

# --------- IMPORTANT -------------------------------------

# Body starts with an initial_shear, will be tilted so that the top advances in small steps of length steps_length until 
# it reaches max_shear.
initial_shear = 0.05 # this is tan theta (theata is the angle between the lateral side at the beginning and the same side at the end)
max_shear = 0.30 # TOP side will be displaced a (max_shearx100)% from the initial position

step_lenght = 0.01

print_Niter = 1

resolution = 7

uc0 = -0.0 # uniform radial displacement, - gives contraction, + gives expansion 
###u0 can also be a vector, i.e. [a, b] for some a and b,  or an expression.

my_k = 10 # Value of the constant k

my_model = 'PNIPAAm_particles' # 'linear_spring' 

# choose if the centers are free to move or not. Cells can move only if one has called the method problem.set_k
solver_t = 'free_centers' # 'fixed_centers' 

# --------- Create domain and mesh START ---------------------

# Le cuadrado
bottom =  CompiledSubDomain("near(x[1], side) && on_boundary", side = 0)
top = CompiledSubDomain("near(x[1], side) && on_boundary", side = 1)#

left = CompiledSubDomain("near(x[0], side) && on_boundary", side = 0)#
right = CompiledSubDomain("near(x[0], side) && on_boundary", side = 1)

# inside = CompiledSubDomain("near(x[0], side) && on_boundary", side = 1)

# The cell

#c_x, c_y, rho = 0.5, 0.5, 0.3 # center and radius of a big cirlce (The ECM)  NO EN MI CASO
c1_x, c1_y, i_rho1 = 0.5, 0.5, 0.24    # center and radius of a cell
#c2_x, c2_y, i_rho2 = -4, 0., 2.   # center and radius of the second cell


cell_domain = CircularDomain(i_rho1, c1_x, c1_y, uc0) # creates a circle with radius rho and center at (c_x, c_y), 0 indicates that
# the outer boundary is fixed
# domain =CircularDomain(rho, c_x, c_y) # here the outer boundary is free


# Define Dirichlet boundary (x = 0 or x = 1)
b = Expression(("0.0",
                "0.0"),pr = initial_shear, degree = 2)
t = Expression(("1*pr",
                "0.0"),
                pr = initial_shear, degree = 2)

l = Expression(("pr*x[1]",
                "0.0"),
                pr = initial_shear, degree = 2)

r = Expression(("pr*x[1]",
                "0.0"),
                pr = initial_shear, degree = 2)


domain = RectangularDomain(0, 0, 1, 1, [[b,bottom], [l,left], [r,right], [t, top]], isMeshUniform = False)

domain.remove_subdomain(cell_domain)
#isMeshUnifor default value is "False". When "True" the area of all triangles in the original mesh are equal. Otherwise triangles ares are as better fit.
# !!! When there are irregular shapes inside (ex: circles for cells) this value has to be False so triangles can addapt to the irregular form.

domain.create_mesh(resolution) # create a mesh over the domain for the given resolution
mesh = domain.get_mesh()
# --------- Create domain and mesh  END ---------------------



# -----------------Define problem type START --------------------------------
# ---------- i.e. mEnergyMinimization or mNonlinearVariational --------------

el_order = 1 # the elements order 
# problem = mEnergyMinimizationProblem('1.0/96.0 * (5*Ic**3 - 9*Ic**2 - \
#                                     12*Ic*J**2 + 12*J**2 + 8) + exp( 80*(0.22 - J) )',
#                                      domain, el_order=el_order)

# problem.set_k(my_k, cell_model = my_model) # accounts for cell response  

# -----------------Define problem type END --------------------------------


print('---------------------------------')
print('------------- PASO 1 ------------')
print('---------------------------------')

# ---------------- Define solver type--------------------------------------


problem = mEnergyMinimizationProblem('1.0/96.0 * (5*Ic**3 - 9*Ic**2 - \
                                    12*Ic*J**2 + 12*J**2 + 8) + exp( 80*(0.22 - J) )',
                                     domain, el_order = el_order)
 
problem.set_k(my_k, cell_model = my_model) # accounts for cell response  


solver = m_ncgSolver(problem, res_path = results_path, solver = solver_t, print_Niter = print_Niter) # let save_run = False, this has to be in the while loop every time

# -----------------Define problem type END --------------------------------


solver.initialization('Polyconvex') # initialize the displacement vector using a polyconvex functional
solver.add_disturbances() # If the program finds a SADdle point or local minimum, disturbances will take the solution out of it and the solution can find the global min.
u = solver.solve()

try:
	solver.save_all_function_info()
except:
	print('hdf5 is not supported')

# ------------- Open necessary files ---------------------------

energies = np.load("results/vtkFiles/saved_functions/energy_i.npy", "r") #opens the file where energies are stored

energyVsShear = open("results/vtkFiles/saved_functions/EnergyVsShear.csv", "w") 

last_line1st = energies[-2] # last energy
energyVsShear.write(str(initial_shear)+','+str(last_line1st)+'\n')

# -------------Ends files management -----------------

# ---------------- solving Ends -----------------------------------------

count = 1

# ------------------- loop begins (shear) ---------------------


while (initial_shear < max_shear):

    u_0 = u
    count = count + 1
    shear = initial_shear + step_lenght
    initial_shear = shear


    print('---------------------------------')
    print '------------- PASO', count, '------------'
    print '        Shear', shear
    print('---------------------------------')

    # Define Dirichlet boundary (y = 0 or y = 1)
    b = Expression(("0.0",
                    "0.0"),pr = initial_shear, degree = 2)
    t = Expression(("1*pr",
                    "0.0"),
                    pr = initial_shear, degree = 2)

    l = Expression(("pr*x[1]",
                    "0.0"),
                    pr = initial_shear, degree = 2)

    r = Expression(("pr*x[1]",
                    "0.0"),
                    pr = initial_shear, degree = 2)

    domain = RectangularDomain(0, 0, 1, 1, [ [b,bottom], [l,left], [r,right], [t, top]], isMeshUniform = False)

    domain.remove_subdomain( CircularDomain(i_rho1, c1_x, c1_y, uc0) )

    domain.create_mesh(resolution) # create a mesh over the domain for the given resolution


    problem = mEnergyMinimizationProblem('1.0/96.0 * (5*Ic**3 - 9*Ic**2 - \
                                        12*Ic*J**2 + 12*J**2 + 8) + exp( 80*(0.22 - J) )',
                                         domain, el_order = el_order)

    problem.set_k(my_k, cell_model = my_model) # accounts for cell response  (you can exclude this line)


    solver = m_ncgSolver(problem, res_path = results_path, solver = solver_t, print_Niter = print_Niter) # let save_run = False

    # initialize the displacement vector using a the previous result
    solver.init_from_function(u_0)

    u = solver.solve()

    try:
	    solver.save_all_function_info()
    except:
	    print('hdf5 is not supported')

    energies = np.load("results/vtkFiles/saved_functions/energy_i.npy", "r") #opens the file where energies are stored

    last_line1st = energies[-2] # last energy
    energyVsShear.write(str(initial_shear)+','+str(last_line1st)+'\n')


# ------------------- loop ends -----------------------


energyVsShear.close()

solver.plot_results()
try:
	solver.save_all_function_info()
except:
	print('hdf5 is not supported')

