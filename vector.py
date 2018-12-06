#imports
import math
from decimal import Decimal, getcontext

getcontext().prec = 30

#Vector class

#notes on python magic methods and how to implement your own are here
#https://rszalski.github.io/magicmethods/
class Vector(object):
    
    CANNOT_NORM_ZERO_VEC_MSG = 'Cannot normalize the zero vector'
    NO_UNIQUE_PARALLEL_COMPONENT_MSG = 'No unique parallel component vectors'
    NO_UNIQUE_ORTHOG_COMPONENT_MSG = 'No unique orthogonal component vectors'
    
    def __init__(self, coordinates):
        try:
            if not coordinates:
                raise ValueError
            self.coordinates = tuple([Decimal(x) for x in coordinates])
            self.dimension = len(self.coordinates)
            
        except ValueError:
            raise ValueError('Coordinates cannot be empty')
            
        except TypeError:
            raise TypeError('Coordinates must be an iterable')
            
    
    def __str__(self):
        return 'Vector: {}'.format(self.coordinates)
    
    def __eq__(self, v):
        return self.coordinates == v.coordinates
    
    def plus(self, v):
        #list comprehension with zip!
        new_coords = [x + y for x,y in zip(self.coordinates, v.coordinates)]
        return Vector(new_coords)
    
    def minus(self, v):
        #list comprehension with zip!
        new_coords = [x - y for x,y in zip(self.coordinates, v.coordinates)]
        return Vector(new_coords)
    
    #think about how numpy does broadcasting across vectors 
    def scalar_multiply(self, c):
        '''
        Parameters
        ==========
        c: float,int 
    
        Returns
        =======
        Vector
        
        '''
        new_coords = [c*x for x in self.coordinates]
        return Vector(new_coords)
    
    def magnitude(self):
        
        coords_squared = [x**2 for x in self.coordinates]
        return math.sqrt(sum(coords_squared))
    
    def normalize(self):
        try:
            magn = self.magnitude()
            return self.scalar_multiply(Decimal('1.0')/magn)
        
        #if vector has magnitude 0
        except ZeroDivisionError:
            raise Exception(self.CANNOT_NORM_ZERO_VEC_MSG)
            
    def compute_inner_product(self, v):
        ip = sum([x*y for x,y in zip(self.coordinates, v.coordinates)])
        return ip
    
    def compute_angle_radians(self, v):
        ip = self.inner_product(v)
        magnitude_self =  self.magnitude()
        magnitude_v = v.magnitude()
        angle_rad = math.acos((ip/(magnitude_self*magnitude_v)))
        return angle_rad
    
    def compute_angle_degrees(self, v):
        degrees_per_radian = 180./math.pi
        return degrees_per_radian * self.angle_radians(v)
    
    
    def compute_angle_with(self, v, in_degrees=False):
        try: 
            u1 = self.normalize()
            u2 = v.normalize()
            angle_in_radians = math.acos(u1.compute_inner_product(u2))
            
            if in_degrees:
                degrees_per_radian = 180./math.pi
                return angle_in_radians*degrees_per_radian
            
        except Exception as e:
            if str(e) == self.CANNOT_NORM_ZERO_VEC_MSG:
                raise Exception('Cannot compute angle with zero vector')
                
            else: 
                raise e

    def compute_orthogonal_component(self, basis):
        
        try:
            parellel_component = compute_parallel_component(basis)
            return self.minus(parallel_component)
        
        except Exception as e:
            if str(e) == self.CANNOT_NORM_ZERO_VEC_MSG:
                raise Exception(self.NO_UNIQUE_ORTHOG_COMPONENT_MSG)
            
            else: 
                raise e
        
    def compute_parallel_component(self, basis):
        
        try:
            unit_vector = basis.normalize()
            scalar_weight = self.inner_product(unit_vector)
            return unit_vector.scalar_multiply(scalar_weight)
        
        except Exception as e:
            if str(e) == self.CANNOT_NORM_ZERO_VEC_MSG:
                raise Exception(self.NO_UNIQUE_PARALLEL_COMPONENT_MSG)
            
            else: 
                raise e
                
    def is_orthogonal_to(self, v, tolerance=1e-10):
        return abs(self.inner_product(v) < tolerance)
    
    def is_parallel_to(self, v):
        return ( 
                    self.is_zero() or 
                    v.is_zero() or 
                    self.compute_angle_with(v) == 0 or 
                    self.compute_angle_with(v) == math.pi 
               )
                
    def is_zero(self, tolerance=1e-10):
        return self.magnitude() < tolerance 
    
    def compute_cross_product(self, v):
        
        x_1, y_1, z_1 = self.coordinates
        x_2, y_2, z_2 = v.coordinates
        cross_coords = [ (y_1*z_2 - y_2*z_1), -(x_1*z_2 - x_2*z-1), (x_1*y_2 - x_2*y_1) ]
        return Vector(new_coords)
    
    def compute_parallelogram_area(self, v):
        cross_product = self.compute_cross_product(v)
        return cross_product.magnitude()
    
    def compute_triangle_area(self, v):
        return self.compute_parallelogram_area(v)/ Decimal('2.0')   