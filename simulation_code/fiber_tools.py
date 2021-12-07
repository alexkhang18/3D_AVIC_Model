# This is a file that contains modularized elements of a stress fiber (sf) contraction simulation.

from dolfin import *
#import mshr
import ufl
import meshio
import numpy as np
from scipy import linalg

'''
	This class overloads the dolfin "Expression" function.
	It defines the initial orientation of the sf distribution.
'''
class FiberOrientation(UserExpression):
        '''
            User-defined stress fiber orientation in the cell
            in the undeformed configuration
        '''
        def __init__(self, F, *args, **kwargs):
            self.F = F
            super().__init__(*args, **kwargs)

        #current expression: surface displacement based update algorithm 
        def eval(self, value, x):
            F_arr = self.F(x).reshape((3,3))
            
            R, U = linalg.polar(F_arr)
            U_w, U_v = np.linalg.eig(U)
            U_arg0 = np.argmax(U_w)
            U_arg2 = np.argmin(U_w)
            U_arg1 = np.argmin(np.abs(np.median(U_w) - U_w))

            U_lambda_1 = U_w[U_arg0]; U_lambda_2 = U_w[U_arg1]; U_lambda_3 = U_w[U_arg2]
            e1_x = U_v[0,U_arg0];    e1_y = U_v[1,U_arg0];    e1_z = U_v[2,U_arg0] 
            e2_x = U_v[0,U_arg1];    e2_y = U_v[1,U_arg1];    e2_z = U_v[2,U_arg1] 
            e3_x = U_v[0,U_arg2];    e3_y = U_v[1,U_arg2];    e3_z = U_v[2,U_arg2] 

            value[0] = e3_x
            value[1] = e3_y
            value[2] = e3_z
            print(value)

        def value_shape(self):
            return (3,)

class ExpressionLevel(UserExpression):
        '''
            User-defined stress fiber orientation in the cell
            in the undeformed configuration
        '''
        def __init__(self, F, *args, **kwargs):
            self.F = F
            super().__init__(*args, **kwargs)

        #current expression: surface displacement based update algorithm 
        def eval(self, value, x):
            F_arr = self.F(as_vector([x[0],x[1],x[2]])).reshape((3,3))
            
            R, U = linalg.polar(F_arr)
            U_w, U_v = np.linalg.eig(U)
            U_arg0 = np.argmax(U_w)
            U_arg2 = np.argmin(U_w)
            U_arg1 = np.argmin(np.abs(np.median(U_w) - U_w))

            U_lambda_1 = U_w[U_arg0]; U_lambda_2 = U_w[U_arg1]; U_lambda_3 = U_w[U_arg2]
            e1_x = U_v[0,U_arg0];    e1_y = U_v[1,U_arg0];    e1_z = U_v[2,U_arg0] 
            e2_x = U_v[0,U_arg1];    e2_y = U_v[1,U_arg1];    e2_z = U_v[2,U_arg1] 
            e3_x = U_v[0,U_arg2];    e3_y = U_v[1,U_arg2];    e3_z = U_v[2,U_arg2] 

            value[0] = 1/U_lambda_3

        def value_shape(self):
            return (1,)

# this represents the scalar contractile strength field that is to be applied to the simulated fibers in the Cytoplasm mesh           
class ContractileStrength(UserExpression):

    # overloading the __init__ function in order to allow time to be tracked.
    def __init__(self, t,**kwargs):

        '''
            User defined contraction strength f
        '''
        self.t = t
        super().__init__(**kwargs)

    
    def eval(self, value, x):
        
        # define contractile strength spatially
        # value[0] = self.t*5
        value[0] = self.t

    def value_shape(self):
        return (1,)

