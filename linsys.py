from decimal import Decimal, getcontext
from copy import deepcopy

from vector import Vector
from plane import Plane

getcontext().prec = 30


class LinearSystem(object):

    ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG = 'All planes in the system should live in the same dimension'
    NO_SOLUTIONS_MSG = 'No solutions'
    INF_SOLUTIONS_MSG = 'Infinitely many solutions'

    def __init__(self, planes):
        try:
            d = planes[0].dimension
            for p in planes:
                assert p.dimension == d

            self.planes = planes
            self.dimension = d

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)


    def swap_rows(self, row1, row2):
        self[row1], self[row2] = self[row2], self[row1]


    def multiply_coefficient_and_row(self, coefficient, row):
        
        n = self[row].normal_vector
        k = self[row].constant_term
        
        new_normal_vector = n.scalar_multiply(coefficient)
        new_constant_term = k*coefficient
        
        self[row] = Plane(normal_vector=new_normal_vector, constant_term=new_constant_term)


    def add_multiple_times_row_to_row(self, coefficient, row_to_add, row_to_be_added_to):
        
        n1 = self[row_to_add].normal_vector
        n2 = self[row_to_be_added_to].normal_vector
        k1 = self[row_to_add].constant_term
        k2 = self[row_to_be_added_to].constant_term
        
        new_normal_vector = n1.scalar_multiply(coefficient).plus(n2)
        new_constant_term = (k1*coefficient) + k2
        
        self[row_to_be_added_to] = Plane(normal_vector=new_normal_vector, constant_term=new_constant_term)
        
        
    def compute_triangular_form(self):
        system = deepcopy(self)
        
        m = len(system.planes)
        n = system.dimension
        j = 0 #the jth variable out of total n variables
        
        for i in range(m):
            while j < n: 
                c = MyDecimal(system.planes[i].normal_vector.coordinates[j])
                if c.is_near_zero():
                    swap_succeeded = system.swap_with_row_below_with_nonzero_coeff(row, col)
                    if not swap_succeeded:
                        j +=1
                        continue
                    
                system.clear_coefficients_below(i,j)
                j += 1
                break
        
        return system
    
    
    def compute_rref(self):
        tf = self.compute_triangular_form()
        
        num_equations = len(tf)
        pivot_indices = tf.indices_of_first_nonzero_terms_in_each_row()
        
        for i in range(num_equations)[::-1]:
            j = pivot_indices[i]
            if j < 0:
                continue
            tf.scale_row_to_make_coefficients_equal_one(i, j)
            tf.clear_coefficients_above(i,j)
        
        return tf
    
    def scale_row_to_make_coefficients_equal_one(self, row, col):
        n = self.planes[row].normal_vector.coordinates
        beta = Decimal(1.0)/n[col]
        self.multiply_coefficient_and_row(beta,row)
        
    def clear_coefficients_above(self, row, col):
        for k in range(row)[::-1]:
            n = self.planes[k].normal_vector.coordinates
            alpha = -(n[col])
            self.add_multiple_times_row_to_row(alpha, row, k)
        
    def swap_with_row_below_with_nonzero_coeff(self,row, col):
            num_equations = len(self.planes)
            for k in range(row+1, num_equations):
                coefficient = MyDecimal(self.planes[k].normal_vector.coordinates[col])
                if not coefficient.is_near_zero():
                    self.swap_rows(row, k)
                    return True
                
            return False
        
    def clear_coefficients_below(self, row, col):
        num_equations = len(self.planes)
        beta = MyDecimal(self.planes[row].normal_vector.coordinates[col])
        for k in range(row+1, num_equations):
            n = self.planes[k].normal_vector.coordinates
            gamma = n[col]
            alpha = -gamma/beta
            self.add_multiple_times_row_to_row(alpha, row, k)
            
            
    def compute_solution(self):
        try:
            return self.do_gaussian_elimination_and_parameterize_solution()
        
        except Exception as e:
            if str(e) == self.NO_SOLUTIONS_MSG:
                return str(e)
            
            else: 
                raise e
                
    def do_gaussian_elimination_and_parameterize_solution(self):
        '''
        Tries to find the vector that runs through a system of equations
        
        '''
        rref = self.compute_rref()
        
        rref.raise_exception_if_contradictory_equation()
        direction_vectors = rref.extract_direction_vectors_for_parameterization()
        basepoint = rref.extract_basepoint_for_parameterization()
        
        
        return Parameterization(basepoint,direction_vectors)
    
    def extract_direction_vectors_for_parameterization(self):
        '''
        
        '''
        num_variables = self.dimension
        pivot_indices = self.indices_of_first_nonzero_terms_in_each_row()
        #returns the inidices that have missing variables
        free_variables_indices = set(range(num_variables)) - set(pivot_indices)
        
        direction_vectors = []

        for free_var in free_variable_indices:
            #create empty coord list for each direction vector
            vector_coords = [0] * num_variables
            #for that given free var, its value is 1 x=x, y=y, z=z
            #when parameterized we remove the var and have a vector of values 
            vector_coords[free_var] = 1
            
            for i, p in enumerate(self.planes):
                pivot_var = pivot_indices[i]
                if pivot_var < 0:
                    break
                
                #set the coordinate to the negative coefficient of the free variable
                vector_coords[pivot_var] = -p.normal_vector[free_var]
            direction_vectors.append(Vector(vector_coords))
            
        return direction_vectors
            
    def extract_basepoint_for_parameterization(self):
        '''
        '''
        num_variables = self.dimension
        pivot_indices = self.indices_of_first_nonzero_terms_in_each_row()
        
        basepoint_coords = [0]*num_variables
        
        for i, p in enumerate(self.planes):
            pivot_var = pivot_indices[i]
            if pivot_var < 0:
                break
                
            basepoint_coords[pivot_var] = p.constant_term
            
        return Vector(basepoint_coords)
    
    
    def raise_exception_if_contradictory_equation():
        '''
        Identifies systems where one or more equation 
        results in 0 = k
        '''
        #make sure each of the planes has a non_zero first coeff
        for p in self.planes:
            try:
                p.first_non_zero_index(p.normal_vector)
                
            #if there are no non_zero elements
            #check the constant to make sure it is also non_zero
            except Exception as e:
                if str(e) == 'No nonzero elements found':
                    constant_term = MyDecimal(p.constant_term)
                    
                    #if constant is non_zero then we no the system is
                    #contradictory and no solutions exist
                    if not constant_term.is_near_zero():
                        raise Exception(self.NO_SOLUTIONS_MSG)
                        
                else:
                    raise e
                    
                    
    def raise_exception_if_too_few_pivots(self):
        '''
        Checks to make sure that if we have a system of n equations that 
        for each row in the system our pivot point is one index to the right
        of the row above. 
        '''
        #this returns a list of the index of the first non-zero element
        #if we are in rref it should be something like [0,1,2,3,..n]
        #for a system of n equations
        pivot_indices = self.indices_of_first_nonzero_terms_in_each_row()
        num_pivots = sum([1 if index >= 0 else 0 for index in pivot_indices])
        num_variables = self.dimension
        
        #If we have fewer equations/pivot points then variables than
        #our system has infinite solutions
        if num_pivots < num_variables:
            raise Exception(self.INF_SOLUTIONS_MSG)
            
        
            
    def indices_of_first_nonzero_terms_in_each_row(self):
        num_equations = len(self.planes)
        num_variables = self.dimension

        indices = [-1] * num_equations

        for i,p in enumerate(self.planes):
            try:
                indices[i] = p.first_nonzero_index(p.normal_vector.coordinates)
            except Exception as e:
                if str(e) == Plane.NO_NONZERO_ELTS_FOUND_MSG:
                    continue
                else:
                    raise e

        return indices


    def __len__(self):
        return len(self.planes)


    def __getitem__(self, i):
        return self.planes[i]


    def __setitem__(self, i, x):
        try:
            assert x.dimension == self.dimension
            self.planes[i] = x

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)


    def __str__(self):
        ret = 'Linear System:\n'
        temp = ['Equation {}: {}'.format(i+1,p) for i,p in enumerate(self.planes)]
        ret += '\n'.join(temp)
        return ret


class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps


# p0 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
# p1 = Plane(normal_vector=Vector(['0','1','0']), constant_term='2')
# p2 = Plane(normal_vector=Vector(['1','1','-1']), constant_term='3')
# p3 = Plane(normal_vector=Vector(['1','0','-2']), constant_term='2')

# s = LinearSystem([p0,p1,p2,p3])

# print(s.indices_of_first_nonzero_terms_in_each_row())
# print('{},{},{},{}'.format(s[0],s[1],s[2],s[3]))
# print(len(s))
# print(s)

# s[0] = p1
# print(s)

# print(MyDecimal('1e-9').is_near_zero())
# print(MyDecimal('1e-11').is_near_zero())