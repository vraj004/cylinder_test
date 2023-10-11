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
RUNTIME_PATH = "/home/jovyan/work/docker-iron/test_files/inflation_test/runtime_files/"

# Unique user number identifiers
(
    coord_n, basis_n, region_n, mesh_n, decomp_n,
    geo_field_n, dep_field_n, mat_field_n,
    eqs_field_n, eqs_n, problem_n, def_field_n,
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
    # Load node and element data
    n_np_xyz, n_idx, n_n = nodes(test_name)
    e_np_map, e_idx, e_n = elems(test_name)

    # +==+ computational infrastructure
    
    # Coordiante system
    cmfe_coord = cmfe.coordinate_setup(coord_n, DIM)
    # Basis {LAGRANGE_HERMITE_TP & QUADRATIC_LAGRANGE}
    cmfe_basis = cmfe.basis_setup(basis_n, XI_N)
    # Region
    cmfe_region = cmfe.region_setup(region_n, cmfe_coord)
    # Mesh (?)
    cmfe_mesh_elem, cmfe_mesh = cmfe.mesh_setup(
        cmfe_region, 
        n_n, 
        mesh_n, 
        DIM, 
        e_n, 
        cmfe_basis, 
        e_idx, 
        e_np_map
    )
    # Decomposition {CALCULATED}
    cmfe_decomp = cmfe.decomposition_setup(cmfe_mesh, decomp_n)

    # +==+ field infrastructure
    # Geometric Field Setup {VALUES}
    cmfe_geo_field = cmfe.geometric_setup(geo_field_n, cmfe_region, cmfe_decomp, n_n, n_idx, n_np_xyz)
    mat_field = cmfe.material_setup(decomp, geo_field, region)

    eqs_set_field, eqs_set = cmfe.equations_setup(region, geo_field, mat_field)
    dep_field, geo_field, eqs_set, eqs = cmfe.dependent_setup(region, decomp, geo_field, eqs_set)
    problem, ctrl_loop = cmfe.problem_setup()
    problem, nl_solver, li_solver, solver, solver_eqs = cmfe.solver_setup(problem, eqs_set)
    bcs, solver_eqs = cmfe.boundary_conditions_setup(solver_eqs, dep_field, n_n, n_np_xyz)

    try:
        problem.Solve()
    except:
        print("Failed")

    cmfe.output(region, decomp, dep_field, mesh, geo_field, n_n, e_n, mesh_e)

    # Wrap it up
    problem.Destroy()
    coord_sys.Destroy()
    region.Destroy()
    basis.Destroy()
    iron.Finalise()

# +==+ ^\_/^ +==+ ^\_/^ +==+ 
# Run check
# +==+ ^\_/^ +==+ ^\_/^ +==+ 

if __name__ == '__main__':
    test_name = "cylinder_hexa"
    main(test_name)