# this file will handle running simulations.

#Ideas:
'''
    1. Accept txtfile or json file input data to determine what sim we are going to run
    2. create functions to output data in a pretty way
        2a. use paraview functions to make super cool paraview visualizations
        2b. create an option to run viz off completely
    3. 

'''

import sys
from dolfin import *
# import fenicstools
#import mshr
import time
import ufl
import csv
# import meshio
import numpy as np
import os
import subprocess
from fenicstools import Probes
from mpi4py import MPI
from fiber_tools import *
from fiber_sim_tools import *
from nodal_tools import *
from scipy import linalg

'''
useful terminal functions: 

To copy using the scp command (transfers files): 
scp -r alex@wccms-a32264.oden.utexas.edu:'/home/alex/3D_TFM_Cell_Model/surface_displacement_bc/output' '/Users/alexkhang/Desktop/FEniCS_output'

To mount:
sudo sshfs -o allow_other,default_permissions alex@wccms-a32264.oden.utexas.edu:/home/alex/3D_TFM_Cell_Model ~/Desktop/system76

To unmount: 
sudo umount ~/Desktop/system76
'''

# list_linear_solver_methods()
# list_krylov_solver_preconditioners()
# quit()

# info(NonlinearVariationalSolver.default_parameters(), 1)
# quit()

start_time = time.time()
comm = MPI.comm_world
rank = comm.Get_rank()
max_rank = comm.Get_size()-1

# this data can be inputed through a txt file
mesh_dir = "../meshes/contracting_cell/"
intput = ""
mesh_file = "process.xdmf"
domain_file = "process_domains.xdmf"
boundary_file = "process_boundaries.xdmf"
output_dir = "../output/contracting_cell_bc" 

# don't touch this yet, idk what it means.
# get defaults for iterative solvers
parameters['linear_algebra_backend'] = 'PETSc'
parameters['form_compiler']['representation'] = 'uflacs'
parameters['form_compiler']['optimize'] = True
parameters['form_compiler']['cpp_optimize'] = True
parameters['form_compiler']['quadrature_degree'] = 3
parameters['krylov_solver']['absolute_tolerance' ]= 1E-8
parameters['krylov_solver']['relative_tolerance'] = 1E-4
parameters['krylov_solver']['maximum_iterations'] = 100000


# parameters["form_compiler"]["quadrature_degree"] = 2
ffc_options = {"optimize": True, \
               "eliminate_zeros": True, \
               "precompute_basis_const": True, \
               "precompute_ip_const": True}
set_log_level(20)

# mu_g = 1
# mu_c = 50
# mu_n = 5000
mu_g = 108e-12
mu_c = 750e-12
mu_n = mu_c * 10
# mu_c = 15e-12
# mu_n = mu_c * 10
mu_sf = 390e-12
nu_g = 0.49
nu_c = 0.49
nu_n = nu_c

def run(ffc_options):
# def run(output_dir,ffc_options):
    '''
    Define the solution, and external fields.
    Run solver to compute the solution and save it periodically.
    '''
    fm = 0
    # chunk = int(fm/100)
    chunk = 50
    dt = 1./chunk
    freq_checkout = 2 

    # loads in data for fiber orientation update algorithm 
    points_1 = np.loadtxt(mesh_dir + 'predicted_normal_vertices.txt',delimiter=',')
    points_2 = np.loadtxt(mesh_dir + 'predicted_endo1_vertices.txt',delimiter=',')
    # centroid = np.mean(points_1,axis=0)
    # stretch = 0.95
    # # F = np.array([[1/np.sqrt(stretch), 0, 0],[0, 1/np.sqrt(stretch), 0,],[0, 0, stretch]])
    # F = np.array([[1, 0, 0],[0, 1, 0,],[0, 0, stretch]])

    # points_2 = np.transpose(np.matmul(F,np.transpose(points_1-centroid))) + centroid


    displacements = points_2 - points_1

    # print(points_1)
    # print(points_2)
    # print(np.max(displacements,axis=0))
    # print(displacements)
    # quit()

    faces = np.loadtxt(mesh_dir + 'faces.txt',delimiter=',')
    faces = faces - 1
    cytoplasm_points = np.loadtxt(mesh_dir + 'cytoplasm_vertices.txt',delimiter=' ')


    # vectors = np.loadtxt('../meshes/bird_cell/vectors_level_set.txt',delimiter=',')
    # points = np.loadtxt('../meshes/bird_cell/points_level_set.txt',delimiter=',')


    mesh, domains, boundaries, V0, V, V2 = setup(mesh_dir,mesh_file,domain_file,boundary_file)
    # V2 = TensorFunctionSpace(mesh, "CG", 1, shape=(3, 3))

    # V2 = TensorFunctionSpace(mesh, "CG", 1, shape=(3, 3))
    u = Function(V)

    B = Expression(('0.','0.','0.'),t=0.,element=V.ufl_element())

    # no traction boundary condition
    T = Constant((0.,0.,0.))

    # assignment of contractile strength, fiber expression levels, and fiber orientation
    f = ContractileStrength(t=0, element=V0.ufl_element())
    # phi = ExpressionLevel(t=1, points_2=points_2, element=V0.ufl_element())
    phi = 1
    # m0 = FiberOrientation(points=points,vectors=vectors)
    m0 = as_vector([0,0,1])

    # create XDMF files for output
    u_all_file = XDMFFile(output_dir + '/solution_all.xdmf')
    material_file = XDMFFile(output_dir + '/domains.xdmf')
    material_file.write(domains)
    m_all_file = XDMFFile(output_dir + '/deformed_fiber_all.xdmf')
    I4_file = XDMFFile(output_dir + '/I4.xdmf')
    Heavy_file = XDMFFile(output_dir + '/Heavy.xdmf')
    Tp_file = XDMFFile(output_dir + '/Tp.xdmf')
    Ta_file = XDMFFile(output_dir + '/Ta.xdmf')
    Tsf_file = XDMFFile(output_dir + '/Tsf.xdmf')
    m0_all_file = XDMFFile(output_dir + '/undeformed_fiber_all.xdmf')
    phi_file = XDMFFile(output_dir + '/phi.xdmf')
        
    # within this, we can add certain info/things
    for n in range(chunk):

        t = n*dt

        if rank == 0:
            print('chunk number = %d'%(n))
            sys.stdout.flush()

        t1 = time.time()
        u, m, I4, Heavyside, Tp, Ta, Tsf, F = solver_simple(n,chunk,u,m0,mesh,domains,V0,V, points_1, faces, displacements, boundaries, B,T,phi,f,mu_g,mu_c,mu_n,mu_sf,nu_g,nu_c,nu_n,ffc_options)
        t2 = time.time()

        if rank == 0:
            print('solve time (s) = ', t2-t1)
            sys.stdout.flush()

        # F = project(F, V2, solver_type = 'gmres', preconditioner_type = 'hypre_amg')
    
        F = project(F, V2, solver_type = 'cg', preconditioner_type = 'amg')


        # # project
        # proj_m = project(m,V)
        # proj_I4 = project(I4,V0)
        # proj_H = project(Heavyside,V0)
        # proj_m0 = project(m0,V)
        # proj_phi = project(phi,V0)
        # # Tp = project(Tp, V2, solver_type = 'cg', preconditioner_type = 'amg')
        # # Ta = project(Ta, V2, solver_type = 'cg', preconditioner_type = 'amg')
        # # Tsf = project(Tsf, V2, solver_type = 'cg', preconditioner_type = 'amg')
        # Tp = project(Tp,V0)
        # Ta = project(Ta,V0)
        # Tsf = project(Tsf,V0)

        # # rename 
        # u.rename('u','displacement')
        # proj_m.rename('m','deformed_fiber_orientation')
        # proj_I4.rename('I4','fiber_stretch')
        # proj_H.rename('H','heaviside')
        # proj_m0.rename('m0','initial_fiber_orientation')
        # proj_phi.rename('phi','fiber_expression_level')
        # Tp.rename('Tp','passive_fiber_stress')
        # Ta.rename('Ta','active_fiber_stress')
        # Tsf.rename('Tsf','total_fiber_stress')

        # # Write Outputs
        # u_all_file.write(u,t)
        # m0_all_file.write(proj_m0,t)
        # m_all_file.write(proj_m,t)
        # I4_file.write(proj_I4,t)
        # Heavy_file.write(proj_H,t)
        # Tp_file.write(Tp,t)
        # Ta_file.write(Ta,t)
        # Tsf_file.write(Tsf,t)
        # phi_file.write(proj_phi,t)

    # evaluates function across all processes and assigns evaluation to each process
    probes = Probes(cytoplasm_points,V2) # defines points to be 'probed' or queried 
    probes(F) # evaluates function at probes locations
    F_nodes = probes.array() # concatenantes function evaluations

    if rank == 0:
        # u_nodes = [u(x) for x in points_1]
        F_nodes = F_nodes
        numData = F_nodes.shape
    else: 
        # while not comm.Iprobe(source = 0):
        #     time.sleep(1)
        numData = None
        # time.sleep(1)

    # broadcast numData and allocate array on other ranks:
    numData = comm.bcast(numData, root=0)
    if rank != 0:
        F_nodes = np.empty(numData, dtype='d') 

    comm.Bcast(F_nodes, root=0) # broadcast the array from rank 0 to all others

    lambda3 = np.zeros((len(cytoplasm_points),1))
    e3 = np.zeros((len(cytoplasm_points),3))

    for i in range(len(cytoplasm_points)):

        F_arr = np.zeros((3,3))
        F_arr[0,0] = F_nodes[i,0]; F_arr[0,1] = F_nodes[i,1]; F_arr[0,2] = F_nodes[i,2]
        F_arr[1,0] = F_nodes[i,3]; F_arr[1,1] = F_nodes[i,4]; F_arr[1,2] = F_nodes[i,5]
        F_arr[2,0] = F_nodes[i,6]; F_arr[2,1] = F_nodes[i,7]; F_arr[2,2] = F_nodes[i,8]

        R, U = linalg.polar(F_arr)
        U_w, U_v = np.linalg.eig(U)
        U_arg0 = np.argmax(U_w)
        U_arg2 = np.argmin(U_w)
        U_arg1 = np.argmin(np.abs(np.median(U_w) - U_w))

        U_lambda_1 = U_w[U_arg0]; U_lambda_2 = U_w[U_arg1]; U_lambda_3 = U_w[U_arg2]
        e1_x = U_v[0,U_arg0];    e1_y = U_v[1,U_arg0];    e1_z = U_v[2,U_arg0] 
        e2_x = U_v[0,U_arg1];    e2_y = U_v[1,U_arg1];    e2_z = U_v[2,U_arg1] 
        e3_x = U_v[0,U_arg2];    e3_y = U_v[1,U_arg2];    e3_z = U_v[2,U_arg2] 
        lambda3[i] = U_lambda_3
        e3[i,0] = e3_x
        e3[i,1] = e3_y
        e3[i,2] = e3_z

    np.savetxt(mesh_dir + 'lambda3.txt',lambda3)
    np.savetxt(mesh_dir + 'e3.txt', e3, delimiter=',')

    # lambda3 = np.zeros((len(cytoplasm_points),1))
    # e3 = np.zeros((len(cytoplasm_points),3))
    # for i in range(len(cytoplasm_points)):
    #     F_arr = F(cytoplasm_points[i,:]).reshape((3,3))
    #     R, U = linalg.polar(F_arr)
    #     U_w, U_v = np.linalg.eig(U)
    #     U_arg0 = np.argmax(U_w)
    #     U_arg2 = np.argmin(U_w)
    #     U_arg1 = np.argmin(np.abs(np.median(U_w) - U_w))

    #     U_lambda_1 = U_w[U_arg0]; U_lambda_2 = U_w[U_arg1]; U_lambda_3 = U_w[U_arg2]
    #     e1_x = U_v[0,U_arg0];    e1_y = U_v[1,U_arg0];    e1_z = U_v[2,U_arg0] 
    #     e2_x = U_v[0,U_arg1];    e2_y = U_v[1,U_arg1];    e2_z = U_v[2,U_arg1] 
    #     e3_x = U_v[0,U_arg2];    e3_y = U_v[1,U_arg2];    e3_z = U_v[2,U_arg2] 
    #     lambda3[i] = U_lambda_3
    #     e3[i,0] = e3_x
    #     e3[i,1] = e3_y
    #     e3[i,2] = e3_z

    # np.savetxt(mesh_dir + 'lambda3.txt',lambda3)
    # np.savetxt(mesh_dir + 'e3.txt', e3, delimiter=',')


    with open(os.path.join(output_dir,"sim_params.txt"), "w+") as f:
        f.write("Mesh: {:s}\n".format(mesh_file))
        # f.write("No. Elements: {:d}\n".format(int(MPI.sum(MPI.comm_world, mesh.num_cells()))))
        f.write("No. Processors: {:f}\n".format(int(comm.Get_size())))
        f.write("No. chunks: {:f}\n".format(chunk))
        f.write("Gel shear modulus (Pa): {:f}\n".format(mu_g))
        f.write("Cyotplasm shear modulus (Pa): {:f}\n".format(mu_c))
        f.write("Nucleus shear modulus (Pa): {:f}\n".format(mu_n))
        f.write("Stress fiber shear modulus (Pa): {:f}\n".format(mu_sf))
        f.write("Gel poisson ratio: {:f}\n".format(nu_g))
        f.write("Cytoplasm poisson ratio: {:f}\n".format(nu_c))
        f.write("Nucleus poisson ratio: {:f}\n".format(nu_n))
        f.write("contraction strength (Pa) = {:f}\n".format(fm))
        # f.write("sim type = {:s}\n".format(sim))
        f.write("Total Time (min.) = {:f}\n".format((time.time() - start_time)/60))

if __name__ == '__main__':
    # run(output_dir,ffc_options)
    run(ffc_options)
    
print("--- %s minutes ---" % ((time.time() - start_time)/60))

