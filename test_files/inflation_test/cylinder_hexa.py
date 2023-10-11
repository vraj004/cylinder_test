"""
Author: Liam Murray, murrayla@student.unimelb.edu.au
Descrption: test openCMISS-iron implementation via hexahedral element 
                cylinder inflation.
Input: runtime_files/
                    cylinder_hexa.ele
                    cylinder_hexa.nodesList
        Files contain element and node data for hexahedral cylinder.
Output: vtk_files/cylinder_hexa.vtk
        vtk files of deformation under inflation
"""

import numpy as np
import meshio
from opencmiss.iron import iron
import cmfe
import os

# +==+ ^\_/^ +==+ ^\_/^ +==+ 
# Parameter Setup
# +==+ ^\_/^ +==+ ^\_/^ +==+ 

# Runtime required parameters
DIM = 3
XI_N = 3
N_N_EL = 27
QUAD_ORDER = 4
X, Y, Z = (1, 2, 3)
PRESSURE_TEST = True
LOADSTEPS = 10
INNER_RAD = 0.375
C_VALS = [2.0, 6.0]
RUNTIME_PATH = "/home/jovyan/work/docker-iron/test_files/inflation_test/runtime_files/"

# Unique user number identifiers
(
    coord_n, basis_n, region_n, mesh_n, decomp_n,
    geo_field_n, dep_field_n, mat_field_n,
    eqs_field_n, eqs_set_n, problem_n, def_field_n,
    pre_field_n
) = range(1, 14)

# +==+ ^\_/^ +==+ ^\_/^ +==+
# Node and Element setup from input files
# +==+ ^\_/^ +==+ ^\_/^ +==+ 

def nodes(test_name):
    n_idx = []
    n_xyz = []
    with open(RUNTIME_PATH + test_name + ".nodes", 'r') as n_file:
        for i, line in enumerate(n_file):
            if i == 0: continue
            line = line.strip().split('\t')
            n_idx.append(int(line[0]))
            n_xyz.append(line[1:])

    return np.array(n_xyz).astype(float), n_idx, i

def elems(test_name):
    e_idx = []
    e_map = []
    with open(RUNTIME_PATH + test_name + ".ele", 'r') as e_file:
        for i, line in enumerate(e_file):
            if i == 0: continue
            line = line.strip().split('\t')
            e_idx.append(i)
            e_map.append(line[3:])

    return np.array(e_map).astype(int), e_idx, i

# +==+ ^\_/^ +==+ ^\_/^ +==+ 
# Main function for safe operation of inflation test
# +==+ ^\_/^ +==+ ^\_/^ +==+ 

def main(test_name):

    # +============+  
    # Node and Element infrastructure
    # +============+ 

    n_np_xyz, n_idx, n_n = nodes(test_name)
    e_np_map, e_idx, e_n = elems(test_name)

    # +============+  
    # Base infrastructure
    # +============+  
   
    cmfe_coord = cmfe.coordinate_setup(coord_n, DIM)
    cmfe_basis = cmfe.basis_setup(basis_n, XI_N)
    cmfe_region = cmfe.region_setup(region_n, cmfe_coord)
    _, cmfe_mesh = cmfe.mesh_setup(
        cmfe_region, n_n, mesh_n, DIM, 
        e_n, cmfe_basis, e_idx, e_np_map
    )
    cmfe_decomp = cmfe.decomposition_setup(cmfe_mesh, decomp_n)

    # +============+  
    # Field and Equation infrastructure
    # +============+  

    cmfe_geo_field = cmfe.geometric_setup(geo_field_n, cmfe_region, cmfe_decomp, n_n, n_idx, n_np_xyz)
    cmfe_mat_field = cmfe.material_setup(mat_field_n, cmfe_decomp, cmfe_geo_field, cmfe_region, C_VALS)
    cmfe_dep_field = cmfe.dependent_setup(dep_field_n, cmfe_region, cmfe_decomp, cmfe_geo_field)
    _, cmfe_eqs_set, _ = cmfe.equations_setup(
        eqs_field_n, eqs_set_n, cmfe_region, 
        cmfe_geo_field, dep_field_n, cmfe_dep_field, 
        mat_field_n, cmfe_mat_field
    )

    # +============+ 
    # Problem and Solution infrastructure
    # +============+ 
    
    cmfe_problem, _, cmfe_solver_eqs = cmfe.problem_solver_setup(problem_n, cmfe_eqs_set, LOADSTEPS)
    _, _ = cmfe.boundary_conditions_setup(cmfe_solver_eqs, cmfe_dep_field, n_n, n_np_xyz)

    # +============+ 
    # Solve
    # +============+  

    try:
        cmfe_problem.Solve()
    except:
        print("Failed")

    # +============+ 
    # Export and Complete
    # +============+ 

    # Wrap it up
    cmfe_problem.Destroy()
    cmfe_coord.Destroy()
    cmfe_region.Destroy()
    cmfe_basis.Destroy()
    iron.Finalise()

# +==+ ^\_/^ +==+ ^\_/^ +==+ 
# Run check
# +==+ ^\_/^ +==+ ^\_/^ +==+ 

if __name__ == '__main__':
    test_name = "cylinder_hexa"
    main(test_name)