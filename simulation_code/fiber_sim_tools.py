from dolfin import *
#import mshr
import ufl
from ufl.core.ufl_type import update_global_expr_attributes
import meshio
import numpy as np
from fiber_tools import *
from nodal_tools import *

def setup(mesh_dir,mesh_file,domain_file,boundary_file):
    '''
    Create simple cell geometry, define different material domains,
    initiate finite element function spaces for scalar and vector variables.
    '''
    # mesh = Mesh()
    # #with XDMFFile('e10_dense_pruned.xdmf') as infile:
    # with XDMFFile(mesh_dir + mesh_file) as infile:
    #     infile.read(mesh)
    # domains = MeshFunction('size_t',mesh,mesh.topology().dim())
    # with XDMFFile(mesh_dir + mesh_file) as infile:
    #     infile.read(domains)

    '''
    100: gel
    200: cytoplasm
    300: nucleus
    '''


    # Gel Volume Mesh
    mesh = Mesh()
    with XDMFFile(mesh_dir + mesh_file) as infile:
        infile.read(mesh)

    mvc = MeshValueCollection("size_t", mesh, 2)
    with XDMFFile(mesh_dir + domain_file) as infile:
        infile.read(mvc, "domains") 
    domains = cpp.mesh.MeshFunctionSizet(mesh, mvc)

    mvc = MeshValueCollection("size_t", mesh, 2)
    with XDMFFile(mesh_dir + boundary_file) as infile:
        infile.read(mvc, "boundaries") 
    boundaries = cpp.mesh.MeshFunctionSizet(mesh, mvc)

    # # Initialization Mesh
    # mesh_init = Mesh()
    # with XDMFFile(mesh_dir + mesh_file) as infile:
    #     infile.read(mesh_init)


    # degree of the finite element space
    degree = 2
    V2 = TensorFunctionSpace(mesh, 'P', 1, shape=(3, 3))
    V = VectorFunctionSpace(mesh, 'P', degree)
    V0 = FunctionSpace(mesh, 'P', degree)
    return mesh, domains, boundaries, V0, V, V2


def solver_simple(i,chunk,u, m0, mesh, domains, V0, V, surface_nodes, surface_faces, displacements, boundaries, B, T, phi, f, mu_g, mu_c, mu_n, mu_sf, nu_g, nu_c, nu_n, ffc_options):
    '''
    Solve the boundary value problem with body force B and boundary traction T
    and active fiber contraction f
    u is the solution to be computed
    '''
    # boundary of the full domain where Dirichlet condition is prescribed
    # def boundary(x, on_boundary):
    #     return on_boundary
    # u_D = Expression(('0.','0.','0.'), degree=2)
    # bc = DirichletBC(V, u_D, boundary)

    # Boundary Conditions
    midpoints = get_midpoints(surface_nodes, surface_faces)
    midpoint_disp = get_midpoint_disp(displacements, surface_faces)
    face_map = get_face_mapping(midpoints, mesh, boundaries, 1)
    face2disp = dict(zip(face_map, midpoint_disp))

    zero = Constant((0.0, 0.0, 0.0))
    bf = BoundaryFunc(mesh, face2disp, 0)
    bf.scalar = (i+1)/chunk  

    outer_bc = DirichletBC(V, zero, boundaries, 101)
    inner_bc = DirichletBC(V, bf, boundaries, 1)
    bcs = [inner_bc, outer_bc]


    # Define variational problem
    du = TrialFunction(V)
    v = TestFunction(V)

    d = u.geometric_dimension()
    I = Identity(d)
    F = I+grad(u)
    C = F.T*F
    Ic = tr(C)
    J = det(F)
    C_bar = C/J**(2./3)
    Ic_bar = tr(C_bar)

    # prescribe material properites for different regions
    # E_g = 1.
    # E_c = 10.
    # E_n = 100.
    # nu_g = 0.49
    # nu_c = 0.2
    # nu_n = 0.4999
    # nu_g = 0.2
    # nu_c = 0.2
    # nu_n = 0.2
    # mu_g = E_g/2/(1+nu_g) 
    # mu_c = E_c/2/(1+nu_c) 
    # mu_n = E_n/2/(1+nu_n) 
    # K_g = E_g/3/(1-2*nu_g)
    # K_c = E_c/3/(1-2*nu_c)
    # K_n = E_n/3/(1-2*nu_n)
    # K_g = mu_g*1000
    # K_c = mu_c*1000
    # K_n = mu_n*1000



    # mu_g = 108
    # # mu_g = 10
    # # mu_c = 5
    # mu_c = 108
    # # mu_n = 15
    # mu_n = mu_g * 10
    # mu_sf = 390
    # # mu_sf = 20
    # # nu_g = 0.4999
    # # nu_c = 0.4999
    # # nu_n = 0.4999
    # nu_g = 0.2
    # nu_c = 0.2
    # nu_n = 0.4999

    K_g = 2*mu_g*(1+nu_g)/(3*(1-2*nu_g))
    K_c = 2*mu_c*(1+nu_c)/(3*(1-2*nu_c))
    K_n = 2*mu_n*(1+nu_n)/(3*(1-2*nu_n))

    dx = Measure('dx',domain=mesh,subdomain_data=domains)
    psi_g = mu_g/2*(Ic_bar-3) + K_g/2*(ln(J))**2
    psi_c = mu_c/2*(Ic_bar-3) + K_c/2*(ln(J))**2
    psi_n = mu_n/2*(Ic_bar-3) + K_n/2*(ln(J))**2
 
    # assemble the total potential energy
    Pi = psi_g*dx(100)+psi_c*dx(200)+psi_n*dx(300) - dot(B,u)*dx('everywhere') - dot(T,u)*ds

    I4 = sqrt(dot(m0,C*m0))

    Heavyside = conditional(I4 > 1, 1, 0)

    m = F*m0/I4 # fiber orientation (deformed)

    # calcualtes passive cauchy stress in the sf
    Tp = phi*Heavyside*2*I4/J*mu_sf*(I4-1)*as_matrix([[m[0]*m[0], m[0]*m[1], m[0]*m[2]],
            [m[1]*m[0], m[1]*m[1], m[1]*m[2]],
            [m[2]*m[0], m[2]*m[1], m[2]*m[2]]
            ])
    Tp_prime = Tp - tr(Tp)/3*as_matrix([[1, 0, 0], # deviatoric stress
            [0, 1, 0],
            [0, 0, 1]
            ])
    Tp_vm = sqrt((3/2)*inner(Tp_prime,Tp_prime)) # Von Mises stress

    # calcualtes active cauchy stress in the sf
    Ta = phi*f*I4/J*as_matrix([[m[0]*m[0], m[0]*m[1], m[0]*m[2]],
            [m[1]*m[0], m[1]*m[1], m[1]*m[2]],
            [m[2]*m[0], m[2]*m[1], m[2]*m[2]]
            ])
    Ta_prime= Ta - tr(Ta)/3*as_matrix([[1, 0, 0], # deviatoric stress
            [0, 1, 0],
            [0, 0, 1]
            ])
    Ta_vm = sqrt((3/2)*inner(Ta_prime,Ta_prime)) # Von Mises stress

    # calcualtes total cauchy stress in the sf
    Tsf = Tp + Ta
    Tsf_prime = Tsf - tr(Tsf)/3*as_matrix([[1, 0, 0], # deviatoric stress
            [0, 1, 0],
            [0, 0, 1]
            ])
    Tsf_vm = sqrt((3/2)*inner(Tsf_prime,Tsf_prime)) # Von Mises stress

    # could try to split f into 2 components to put sf density into the model
    # for now, assume that everything is uniformly distributed
    # still assuming fiber only concentrates along mean direction.

    # pilo-kirchoff stress in the sf
    Psf = J*Tsf*inv(F.T)

    # take Gateaux derivative of Pi
    # A = derivative(Pi, u, v) + inner(Psf,grad(v))*dx(200)
    A = derivative(Pi, u, v) #excludes fibers
    # calculate Jacobian
    J = derivative(A, u, du)

    # Compute solution
    solve(A == 0, u, bcs, J=J, 
        # solver_parameters={"newton_solver":{"relative_tolerance": 1e-8,
        #     "absolute_tolerance": 1e-8}},
        solver_parameters={"newton_solver":{
            "relative_tolerance": 1e-8,
            "absolute_tolerance": 1e-8,
            "linear_solver" : 'gmres',
            'preconditioner': 'hypre_amg'}},
            # "linear_solver" : 'cg'}},
        form_compiler_parameters=ffc_options)

    # return u, B, m
    return u, m, I4, Heavyside, Tp_vm, Ta_vm, Tsf_vm, F

def solver(u, m0, mesh, domains, V0, V, B, T, phi, f, mu_g, mu_c, mu_n, mu_sf, nu_g, nu_c, nu_n, ffc_options):
    '''
    Solve the boundary value problem with body force B and boundary traction T
    and active fiber contraction f
    u is the solution to be computed
    '''
    # boundary of the full domain where Dirichlet condition is prescribed
    def boundary(x, on_boundary):
        return on_boundary
    u_D = Expression(('0.','0.','0.'), degree=2)
    bc = DirichletBC(V, u_D, boundary)

    # Define variational problem
    du = TrialFunction(V)
    v = TestFunction(V)

    d = u.geometric_dimension()
    I = Identity(d)
    F = I+grad(u)
    C = F.T*F
    Ic = tr(C)
    J = det(F)
    C_bar = C/J**(2./3)
    Ic_bar = tr(C_bar)

    # prescribe material properites for different regions
    # E_g = 1.
    # E_c = 10.
    # E_n = 100.
    # nu_g = 0.49
    # nu_c = 0.2
    # nu_n = 0.4999
    # nu_g = 0.2
    # nu_c = 0.2
    # nu_n = 0.2
    # mu_g = E_g/2/(1+nu_g) 
    # mu_c = E_c/2/(1+nu_c) 
    # mu_n = E_n/2/(1+nu_n) 
    # K_g = E_g/3/(1-2*nu_g)
    # K_c = E_c/3/(1-2*nu_c)
    # K_n = E_n/3/(1-2*nu_n)
    # K_g = mu_g*1000
    # K_c = mu_c*1000
    # K_n = mu_n*1000



    # mu_g = 108
    # # mu_g = 10
    # # mu_c = 5
    # mu_c = 108
    # # mu_n = 15
    # mu_n = mu_g * 10
    # mu_sf = 390
    # # mu_sf = 20
    # # nu_g = 0.4999
    # # nu_c = 0.4999
    # # nu_n = 0.4999
    # nu_g = 0.2
    # nu_c = 0.2
    # nu_n = 0.4999

    K_g = 2*mu_g*(1+nu_g)/(3*(1-2*nu_g))
    K_c = 2*mu_c*(1+nu_c)/(3*(1-2*nu_c))
    K_n = 2*mu_n*(1+nu_n)/(3*(1-2*nu_n))

    dx = Measure('dx',domain=mesh,subdomain_data=domains)
    psi_g = mu_g/2*(Ic_bar-3) + K_g/2*(ln(J))**2
    psi_c = mu_c/2*(Ic_bar-3) + K_c/2*(ln(J))**2
    psi_n = mu_n/2*(Ic_bar-3) + K_n/2*(ln(J))**2
 
    # assemble the total potential energy
    Pi = psi_g*dx(100)+psi_c*dx(200)+psi_n*dx(300) - dot(B,u)*dx('everywhere') - dot(T,u)*ds

    I4 = sqrt(dot(m0,C*m0))

    Heavyside = conditional(I4 > 1, 1, 0)

    m = F*m0/I4 # fiber orientation (deformed)

    # calcualtes passive cauchy stress in the sf
    Tp = phi*Heavyside*2*I4/J*mu_sf*(I4-1)*as_matrix([[m[0]*m[0], m[0]*m[1], m[0]*m[2]],
            [m[1]*m[0], m[1]*m[1], m[1]*m[2]],
            [m[2]*m[0], m[2]*m[1], m[2]*m[2]]
            ])
    Tp_prime = Tp - tr(Tp)/3*as_matrix([[1, 0, 0], # deviatoric stress
            [0, 1, 0],
            [0, 0, 1]
            ])
    Tp_vm = sqrt((3/2)*inner(Tp_prime,Tp_prime)) # Von Mises stress

    # calcualtes active cauchy stress in the sf
    Ta = phi*f*I4/J*as_matrix([[m[0]*m[0], m[0]*m[1], m[0]*m[2]],
            [m[1]*m[0], m[1]*m[1], m[1]*m[2]],
            [m[2]*m[0], m[2]*m[1], m[2]*m[2]]
            ])
    Ta_prime= Ta - tr(Ta)/3*as_matrix([[1, 0, 0], # deviatoric stress
            [0, 1, 0],
            [0, 0, 1]
            ])
    Ta_vm = sqrt((3/2)*inner(Ta_prime,Ta_prime)) # Von Mises stress

    # calcualtes total cauchy stress in the sf
    Tsf = Tp + Ta
    Tsf_prime = Tsf - tr(Tsf)/3*as_matrix([[1, 0, 0], # deviatoric stress
            [0, 1, 0],
            [0, 0, 1]
            ])
    Tsf_vm = sqrt((3/2)*inner(Tsf_prime,Tsf_prime)) # Von Mises stress

    # could try to split f into 2 components to put sf density into the model
    # for now, assume that everything is uniformly distributed
    # still assuming fiber only concentrates along mean direction.

    # pilo-kirchoff stress in the sf
    Psf = J*Tsf*inv(F.T)

    # take Gateaux derivative of Pi
    A = derivative(Pi, u, v) + inner(Psf,grad(v))*dx(200)
    # calculate Jacobian
    J = derivative(A, u, du)

    # Compute solution
    solve(A == 0, u, bc, J=J, 
        # solver_parameters={"newton_solver":{"relative_tolerance": 1e-8,
        #     "absolute_tolerance": 1e-8}},
        solver_parameters={"newton_solver":{"relative_tolerance": 1e-8,
            "absolute_tolerance": 1e-8,
            "linear_solver" : 'gmres',
            'preconditioner': 'hypre_amg'}},
        form_compiler_parameters=ffc_options)

    # return u, B, m
    return u, m, I4, Heavyside, Tp_vm, Ta_vm, Tsf_vm