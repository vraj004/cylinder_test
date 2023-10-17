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
from opencmiss.iron import iron
import cmfe

# +==+ ^\_/^ +==+ ^\_/^ +==+ 
# Parameter Setup
# +==+ ^\_/^ +==+ ^\_/^ +==+ 

# Runtime required parameters
DIM = 3
XI_N = 3
N_N_EL = 27
QUAD_ORDER = 4
X, Y, Z, P = (1, 2, 3, 4)
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
    n_np_xyz = np.array(n_xyz).astype(float)
    return n_np_xyz, n_idx, i

def elems(test_name):
    e_idx = []
    e_map = []
    with open(RUNTIME_PATH + test_name + ".ele", 'r') as e_file:
        for i, line in enumerate(e_file):
            if i == 0: continue
            line = line.strip().split('\t')
            e_idx.append(i)
            e_map.append(line[3:])
    e_np_map = np.array(e_map).astype(int)
    return e_np_map, e_idx, i

# +==+ ^\_/^ +==+ ^\_/^ +==+ 
# Main function for safe operation of inflation test
# +==+ ^\_/^ +==+ ^\_/^ +==+ 

def main(test_name):

    # +============+  
    # Base infrastructure
    # +============+  
   
    cmfe_coord = cmfe.coordinate_setup(coord_n, DIM)
    cmfe_basis = cmfe.basis_setup(basis_n, XI_N)
    cmfe_region = cmfe.region_setup(region_n, cmfe_coord)

    # +============+  
    # Nodes and Element infrastructure
    # +============+

    n_np_xyz, n_idx, n_n = nodes(test_name)
    e_np_map, e_idx, e_n = elems(test_name)

    # +============+
    # Mesh infrastructure
    # +============+ 

    cmfe_node = iron.Nodes()
    cmfe_mesh = iron.Mesh()
    cmfe_node.CreateStart(cmfe_region, n_n)
    cmfe_node.CreateFinish()
    cmfe_mesh.CreateStart(mesh_n, cmfe_region, DIM)
    cmfe_mesh.NumberOfElementsSet(e_n)
    cmfe_mesh.NumberOfComponentsSet(1) 
    cmfe_mesh_e = iron.MeshElements()
    cmfe_mesh_e.CreateStart(cmfe_mesh, 1, cmfe_basis)
    for i in range(e_n):
        nodesList = list(
            map(int,e_np_map[i][:])
        )
        cmfe_mesh_e.NodesSet(e_idx[i], nodesList)
    cmfe_mesh_e.CreateFinish()
    cmfe_mesh.CreateFinish()
    print('+==+ ^\_/^ MESH COMPLETE')
    
    # +============+  
    # Decomposition and Geometry infrastructure
    # +============+  
    
    cmfe_decomp = cmfe.decomposition_setup(cmfe_mesh, decomp_n)
    cmfe_geo_field = cmfe.geometric_setup(geo_field_n, cmfe_region, cmfe_decomp, n_idx, n_np_xyz)
    comp_nodes_n = iron.ComputationalNumberOfNodesGet()
    for i, idx in enumerate(n_idx):
        cmfe_geo_field.ParameterSetUpdateNodeDP(
            iron.FieldVariableTypes.U, 
            iron.FieldParameterSetTypes.VALUES,
            1, 1, idx, X, n_np_xyz[i, 0]
        )
        cmfe_geo_field.ParameterSetUpdateNodeDP(
            iron.FieldVariableTypes.U, 
            iron.FieldParameterSetTypes.VALUES,
            1, 1, idx, Y, n_np_xyz[i, 1]
        )
        cmfe_geo_field.ParameterSetUpdateNodeDP(
            iron.FieldVariableTypes.U, 
            iron.FieldParameterSetTypes.VALUES,
            1, 1, idx, Z, n_np_xyz[i, 2]
        )
    cmfe_geo_field.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
    cmfe_geo_field.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
    print('+==+ ^\_/^ GEOMETRIC FIELD COMPLETE')

    # +============+  
    # Material and Dependent infrastructure
    # +============+  

    cmfe_mat_field = cmfe.material_setup(mat_field_n, cmfe_decomp, cmfe_geo_field, cmfe_region, C_VALS)
    cmfe_dep_field = cmfe.dependent_setup(dep_field_n, cmfe_region, cmfe_decomp, cmfe_geo_field)

    # +============+  
    # Equation infrastructure
    # +============+ 

    cmfe_eqs_set_field = iron.Field()
    cmfe_eqs_set = iron.EquationsSet()
    cmfe_eqs_set_specs = [
        iron.ProblemClasses.ELASTICITY,
        iron.ProblemTypes.FINITE_ELASTICITY,
        iron.EquationsSetSubtypes.MOONEY_RIVLIN
    ]
    cmfe_eqs_set.CreateStart(
        eqs_set_n, cmfe_region, cmfe_geo_field, 
        cmfe_eqs_set_specs, eqs_field_n, cmfe_eqs_set_field
    )
    cmfe_eqs_set.CreateFinish()
    cmfe_eqs_set.DependentCreateStart(dep_field_n, cmfe_dep_field)
    cmfe_eqs_set.DependentCreateFinish()
    cmfe_eqs_set.MaterialsCreateStart(mat_field_n, cmfe_mat_field)
    cmfe_eqs_set.MaterialsCreateFinish()
    cmfe.equations_setup(cmfe_eqs_set)

    # +============+ 
    # Solve
    # +============+  

    p = 0.0

    for _ in range(0,1):
        pre_inc = [1500/6] * 6
        tolerances = [1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5]
        its = 10

        for i in range(0, len(pre_inc)):
            inc = pre_inc[i]
            tol = tolerances[i]

            # +============+ 
            # Problem and Solution infrastructure
            # +============+ 
            
            cmfe_problem, cmfe_solver, cmfe_solver_eqs = cmfe.problem_solver_setup(problem_n, cmfe_eqs_set, its)
            cmfe.boundary_conditions_setup(cmfe_solver_eqs, cmfe_dep_field, n_n, n_np_xyz, inc)

            cmfe_problem.Solve()
            cmfe_problem.Finalise()
            cmfe_solver_eqs.Finalise()

    # +============+ 
    # Deformed Fields and Export infrastructure
    # +============+ 

    cmfe_def_field = cmfe.deformed_setup(def_field_n, cmfe_region, cmfe_decomp, cmfe_dep_field)
    cmfe_pre_field = cmfe.pressure_setup(cmfe_region, cmfe_decomp, pre_field_n)
    cmfe_dep_field.ParametersToFieldParametersComponentCopy(
        iron.FieldVariableTypes.U,
        iron.FieldParameterSetTypes.VALUES, X,
        cmfe_def_field, iron.FieldVariableTypes.U,
        iron.FieldParameterSetTypes.VALUES, X
    )
    cmfe_dep_field.ParametersToFieldParametersComponentCopy(
        iron.FieldVariableTypes.U,
        iron.FieldParameterSetTypes.VALUES, Y,
        cmfe_def_field, iron.FieldVariableTypes.U,
        iron.FieldParameterSetTypes.VALUES, Y
    )
    cmfe_dep_field.ParametersToFieldParametersComponentCopy(
        iron.FieldVariableTypes.U,
        iron.FieldParameterSetTypes.VALUES, Z,
        cmfe_def_field, iron.FieldVariableTypes.U,
        iron.FieldParameterSetTypes.VALUES, Z
    )
    cmfe_dep_field.ParametersToFieldParametersComponentCopy(
        iron.FieldVariableTypes.U,
        iron.FieldParameterSetTypes.VALUES,
        P,
        cmfe_pre_field, 
        iron.FieldVariableTypes.U,
        iron.FieldParameterSetTypes.VALUES, 
        1
    )
    cmfe.vtk_output(cmfe_mesh, n_n, cmfe_geo_field, cmfe_dep_field, e_np_map, cmfe_mesh_e, RUNTIME_PATH)

    # +============+ 
    # Wrap it up
    # +============+ 

    cmfe_problem.Destroy()
    cmfe_coord.Destroy()
    cmfe_region.Destroy()
    cmfe_basis.Destroy()
    iron.Finalise()

# +==+ ^\_/^ +==+ ^\_/^ +==+ 
# Run check
# +==+ ^\_/^ +==+ ^\_/^ +==+ 

if __name__ == '__main__':
    test_name = "cylinder_hexa_test"
    main(test_name)