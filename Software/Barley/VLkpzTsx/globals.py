# Global variables for stoma simulations

# Directory and File details
initial_mesh = 'baseline-stress.mdxm'
out_dir = '/home/cdurney/stomata/data/symplastic/ml_plot/'
outdirSuffix = ''

# Initialization variables
rod_thickness = 3.0
wall_thickness = 1.0
## Guard Cell
E1_E3_gc = 40 
E2_gc = 75 
## Subsidiary Cell
E1_E3_sc = 40 
E2_sc = 75 
## Initial Pressure
p_gc_i = 0.5 
p_sc_i = 0.2 

poisson_coef = 0.3

# pressure incrementer
num_points = 10 
gc_start = 1.0 
gc_end = 10.0 
sc_start = 0.2 
sc_end = 0.2 



