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
import argparse
from opencmiss.iron import iron
import cmfe_pressure as cmfe

# +==+ ^\_/^ +==+ ^\_/^ +==+ 
# Parameter Setup
# +==+ ^\_/^ +==+ ^\_/^ +==+ 

# Runtime required parameters
DIM = 3
XI_N = 3
N_N_EL = 27
QUAD_ORDER = 4
X, Y, Z, P = (1, 2, 3, 4)
PRESSURE_NODE = True
LOADSTEPS = 1
INNER_RAD = 0.375
C_VALS = [2, 6]
PRESSURE_VAL = 0.0
TRANSLATION_VAL = 0.0
RUNTIME_PATH = "/home/jovyan/work/vraj004/cylinder_test/test_files/inflation_test/runtime_files/"

# Iron Numbering for 27?
#   *  z = 0           z = 0.5         z = 1          
#   *  6--7 --8     * 15--16--17    x 24--25--26
#   *  |      |     *  |      |     x  |      |
#   *  3  4   5     * 12  13  14    x 21  22  23
#   *  |      |     *  |      |     x  |      |
#   *  0--1 --2     *  9--10--11    x 18--19--20

# Gmsh Numbering for Hexa-27
#   *  z = 0           z = 0.5         z = 1    
#   *  3--13--2     * 15--24--14    *  7--19--6      
#   *  |      |     *  |      |     *  |      |       
#   *  9  20  11    * 22  26  23    * 17  25  18     
#   *  |      |     *  |      |     *  |      |     
#   *  0-- 8--1     * 10--21--12    *  4--16--5       

# VTK Numbering for Hexa-27
#   *  z = 0           z = 0.5         z = 1    
#   *  3--10--2     * 19--23--18    x  7--14--6
#   *  |      |     *  |      |     x  |      |
#   * 11  24  9     * 20  26  21    x 15  25  13
#   *  |      |     *  |      |     x  |      |
#   *  0-- 8--1     * 16--22--17    x  4--12--5 

GMSH2VTK = [
    0, 1, 2, 3, 4, 5, 6, 
    7, 8, 11, 13, 9, 16, 18,
    19, 17, 10, 12, 14, 15, 22,
    23, 21, 24, 20, 25, 26
]
GMSH2IRON = [
    0, 8, 1, 9, 20, 11, 3, 
    13, 2, 10, 21, 12, 22, 
    26, 23, 15, 24, 14, 4, 16, 
    5, 17, 25, 18, 7, 19, 6
]
IRON_VTK = [
    0, 2, 8, 6, 18, 20, 26, 
    24, 1, 5, 7, 3, 19, 13, 
    25, 21, 9, 11, 17, 15, 12,
    14, 10, 16, 4, 22, 13
]

# # Bottom Left -> Anti-Clockwise {VTK Standard}
# VERTICES = {
#     "IRON": [0, 2, 8, 6, 18, 20, 26, 24],
#     "GMSH": [0, 1, 2, 3, 4, 5, 6, 7],
#     "VTK": [0, 1, 2, 3, 4, 5, 6, 7]
# }
# # Bottom Middle -> Anti-Clockwise {VTK Standard}
# MIDEDGE = {
#     "IRON": [1, 5, 7, 3, 19, 23, 25, 21, 9, 11, 17, 15],
#     "GMSH": [8, 11, 13, 9, 16, 18, 19, 17, 10, 12, 14, 15],
#     "VTK": [8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
# }
# # Middle Left -> Across, Down, Up, Bottom, Top, Middle {VTK Standard}
# FACEINTERIOR = {
#     "IRON": [12, 14, 10, 16, 4, 22, 13],
#     "GMSH": [22, 23, 21, 24, 20, 25, 26],
#     "VTK": [20, 21, 22, 23, 24, 25, 26]
# }

# Iron Numbering for 27?
#   *  z = 0       x  z = 1         
#   *  2------3    x  6------7
#   *  |      |    x  |      |
#   *  |      |    x  |      |
#   *  |      |    x  |      |
#   *  0------1    x  4------5

# Gmsh Numbering for Hexa-27
#   *  z = 0       x  z = 1         
#   *  3------2    x  7------6
#   *  |      |    x  |      |
#   *  |      |    x  |      |
#   *  |      |    x  |      |
#   *  0------1    x  4------5   

# VTK Numbering for Hexa-27
#   *  z = 0       x  z = 1         
#   *  3------2    x  7------6
#   *  |      |    x  |      |
#   *  |      |    x  |      |
#   *  |      |    x  |      |
#   *  0------1    x  4------5

GMSH2IRON_LIN = [
    0, 1, 3, 2,
    4, 5, 7, 6
]
IRON2VTK_LIN = [
    0, 1, 3, 2,
    4, 5, 7, 6
]

# Unique user number identifiers
(
    coord_n, quad_basis_n, lin_basis_n, 
    region_n, mesh_n, decomp_n,
    geo_field_n, dep_field_n, mat_field_n,
    eqs_field_n, eqs_set_n, problem_n, 
    def_field_n, pre_field_n
) = range(1, 15)

# +==+ ^\_/^ +==+ ^\_/^ +==+
# Node and Element setup from input files
# +==+ ^\_/^ +==+ ^\_/^ +==+ 

# def cvtNodeNumbering(prior, after):
#     vtk = np.array(VERTICES["VTK"] + MIDEDGE["VTK"] + FACEINTERIOR["VTK"])
#     input_nodes = np.array(VERTICES[prior] + MIDEDGE[prior] + FACEINTERIOR[prior])
#     output_nodes = np.array(VERTICES[after] + MIDEDGE[after] + FACEINTERIOR[after])
#     cvtInput2VTK = []
#     cvtVTK2Output = []
#     for node in vtk:
#         cvtInput2VTK.append(np.where(input_nodes==node)[0][0])
#     for node in output_nodes:
#         cvtVTK2Output.append(np.where(input_nodes==node)[0][0])
#     return cvtVTK2Output

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

def main(test_name, test_type, pressure_test):

    # +============+  
    # Base infrastructure
    # +============+  
    
    # +==+ cmfe coordainte system

    cmfe_coord = cmfe.coordinate_setup(coord_n, DIM)
    print('+==+ COORDINATE SYSTEM COMPLETE')
    # +==+ cmfe basis system
    cmfe_quad_basis = cmfe.basis_setup(quad_basis_n, XI_N, "quadratic")
    print('+==+ QUAD BASIS SYSTEM COMPLETE')
    # +==+ cmfe region system
    cmfe_region = cmfe.region_setup(region_n, cmfe_coord)
    print('+==+ REGION COMPLETE')

    # +============+  
    # Nodes and Element infrastructure
    # +============+

    n_np_xyz, n_idx, n_n = nodes(test_name)
    e_np_map, e_idx, e_n = elems(test_name)

    # +============+  
    # Pressure Setup
    # +============+
    if pressure_test:
        cmfe_lin_basis = cmfe.basis_setup(lin_basis_n, XI_N, "linear")
        print('+==+ LIN BASIS SYSTEM COMPLETE')

    # +============+
    # Mesh infrastructure
    # +============+ 

    # +==+ cmfe nodes
    cmfe_node = iron.Nodes()
    cmfe_node.CreateStart(cmfe_region, n_n)
    cmfe_node.CreateFinish()
    # +==+ cmfe mesh
    cmfe_mesh = iron.Mesh()
    cmfe_mesh.CreateStart(mesh_n, cmfe_region, DIM)
    cmfe_mesh.NumberOfElementsSet(e_n)
    if pressure_test:
        cmfe_mesh.NumberOfComponentsSet(2) 
    else:
        cmfe_mesh.NumberOfComponentsSet(1) 
    # +==+ cmfe quadratic mesh elements
    cmfe_mesh_e_quad = iron.MeshElements()
    cmfe_mesh_e_quad.CreateStart(cmfe_mesh, 1, cmfe_quad_basis)
    # += allocating nodes to elements
    print('+= ... begin quadratic mesh allocation')
    for i in range(e_n):
        node_allocs = list(
            map(int,[e_np_map[i][idx] for idx in GMSH2IRON])
        )
        cmfe_mesh_e_quad.NodesSet(e_idx[i], node_allocs)
    # +=
    cmfe_mesh_e_quad.CreateFinish()
    if pressure_test:
        # +==+ cmfe linear mesh elements
        cmfe_mesh_e_lin = iron.MeshElements()
        cmfe_mesh_e_lin.CreateStart(cmfe_mesh, 2, cmfe_lin_basis)
        # += allocating nodes to elements
        print('+= ... begin linear mesh allocation')
        for i in range(e_n):
            node_allocs = list(
                map(int,[e_np_map[i][idx] for idx in GMSH2IRON_LIN])
            )
            cmfe_mesh_e_lin.NodesSet(e_idx[i], node_allocs)
        # +=
        cmfe_mesh_e_lin.CreateFinish()
    cmfe_mesh.CreateFinish()
    print('+==+ MESH ALLOCATION COMPLETE')
    
    # +============+  
    # Decomposition and Geometry infrastructure
    # +============+  
    
    # +==+ cmfe decomposition field
    cmfe_decomp = cmfe.decomposition_setup(cmfe_mesh, decomp_n)
    print('+==+ DECOMPOSITION FIELD COMPLETE')
    # +==+ cmfe geometric field
    cmfe_geo_field = cmfe.geometric_setup(geo_field_n, cmfe_region, cmfe_decomp, n_idx, n_np_xyz)
    # += setup geometric field with undeformed coordiantes
    print('+= ... begin undeformed mesh setup')
    for i, idx in enumerate(n_idx):
        for j in [X, Y, Z]:
            cmfe_geo_field.ParameterSetUpdateNodeDP(
                iron.FieldVariableTypes.U, 
                iron.FieldParameterSetTypes.VALUES,
                1, 1, idx, j, n_np_xyz[i, j-1]
            )
    # += 
    cmfe_geo_field.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
    cmfe_geo_field.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
    print('+==+ GEOMETRIC FIELD COMPLETE')

    # +============+  
    # Material and Dependent infrastructure
    # +============+  

    # +==+ cmfe material field
    cmfe_mat_field = cmfe.material_setup(mat_field_n, cmfe_decomp, cmfe_geo_field, cmfe_region, C_VALS)
    print('+==+ MATERIAL FIELD COMPLETE')
    # +==+ cmfe dependent field
    cmfe_dep_field = cmfe.dependent_setup(dep_field_n, cmfe_region, cmfe_decomp, cmfe_geo_field, pressure_test)
    for i in [X, Y, Z]:
        iron.Field.ParametersToFieldParametersComponentCopy(
            cmfe_geo_field, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, i,
            cmfe_dep_field, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, i
        )
    iron.Field.ComponentValuesInitialiseDP(
        cmfe_dep_field, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, P, -8.0
    )
    print('+==+ DEPENDENT FIELD COMPLETE')

    # +============+  
    # Equation infrastructure
    # +============+ 

    # +==+ cmfe equation set field 
    cmfe_eqs_set_field = iron.Field()
    cmfe_eqs_set_specs = [
        iron.ProblemClasses.ELASTICITY,
        iron.ProblemTypes.FINITE_ELASTICITY,
        iron.EquationsSetSubtypes.MOONEY_RIVLIN
    ]
    # +==+ cmfe equation set 
    cmfe_eqs_set = iron.EquationsSet()
    cmfe_eqs_set.CreateStart(
        eqs_set_n, cmfe_region, cmfe_geo_field, cmfe_eqs_set_specs, eqs_field_n, cmfe_eqs_set_field
    )
    cmfe_eqs_set.CreateFinish()
    cmfe_eqs_set.DependentCreateStart(dep_field_n, cmfe_dep_field)
    cmfe_eqs_set.DependentCreateFinish()
    cmfe_eqs_set.MaterialsCreateStart(mat_field_n, cmfe_mat_field)
    cmfe_eqs_set.MaterialsCreateFinish()
    # +==+ cmfe equations
    cmfe.equations_setup(cmfe_eqs_set)
    print('+==+ EQUATION FIELD COMPLETE')

    # +============+ 
    # Solve
    # +============+  

    # +==+ Export field information so far
    fields = iron.Fields()
    fields.CreateRegion(cmfe_region)
    fields.NodesExport("Output", "FORTRAN")
    fields.ElementsExport("Output", "FORTRAN")
    fields.Finalise()
    # += iterations through increments for solution
    # pressure_inc = [1.0 * 1/25, 1.0 * 1/25, 1.0 * 1/25, 2.0 * 1/25, 2.0 * 1/25, 2.0 * 1/25, 1.0 * 1/25]
    pressure_inc = [PRESSURE_VAL] * LOADSTEPS
    translation_inc = [TRANSLATION_VAL/LOADSTEPS] * len(pressure_inc)
    print('+= ... begin solver')
    for i in range(0, len(pressure_inc), 1):
        curr_p = pressure_inc[i]
        curr_t = translation_inc[i]
        # +============+ 
        # Problem and Solution infrastructure
        # +============+ 
        
        # +==+ cmfe problem and solver field
        cmfe_problem, cmfe_solver, cmfe_solver_eqs = cmfe.problem_solver_setup(
            problem_n, cmfe_eqs_set, LOADSTEPS
        )
        # +==+ cmfe boundary conditions
        cmfe.boundary_conditions_setup(cmfe_solver_eqs, cmfe_dep_field, n_n, n_np_xyz, curr_p, curr_t)
        # += solver for current iterations
        print("+===============================================================+")
        print(f'+= increment {i} =+')
        print("+===============================================================+")
        cmfe_problem.Solve()
        cmfe_problem.Finalise()
        cmfe_solver_eqs.Finalise()
    print('+==+ SOLVER COMPLETE')

    # +============+ 
    # Deformed Fields and Export infrastructure
    # +============+ 

    # +==+ cmfe deformed field
    cmfe_def_field = cmfe.deformed_setup(def_field_n, cmfe_region, cmfe_decomp, cmfe_dep_field)
    print('+==+ DEFORMED FIELD COMPLETE')
    # +==+ cmfe pressure field
    cmfe_pre_field = cmfe.pressure_setup(cmfe_region, cmfe_decomp, pre_field_n, pressure_test)
    print('+==+ PRESSURE FIELD COMPLETE')
    # += setup deformed field with new values
    print('+= ... begin deformed mesh setup')
    for i in [X, Y, Z]:
        cmfe_dep_field.ParametersToFieldParametersComponentCopy(
            iron.FieldVariableTypes.U,
            iron.FieldParameterSetTypes.VALUES, i,
            cmfe_def_field, iron.FieldVariableTypes.U,
            iron.FieldParameterSetTypes.VALUES, i
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
    print('+==+ DEPENDENT FIELD COMPLETE')

    # cmfe & meshio output
    cmfe.vtk_output(
        cmfe_mesh, n_n, cmfe_geo_field, cmfe_dep_field, e_np_map, 
        cmfe_mesh_e_quad, RUNTIME_PATH, test_name, None, test_type
    )
    print('+==+ EXPORT COMPLETE')

    # +============+ 
    # Wrap it up
    # +============+ 

    # results = {}
    # elementNumber = 1
    # p = cmfe_dep_field.ParameterSetGetElement(
    # iron.FieldVariableTypes.U,
    # iron.FieldParameterSetTypes.VALUES,elementNumber,4)
    # results['Hydrostatic pressure'] = p
    # print("Hydrostatic pressure")
    # print(p)

    # cmfe_problem.Destroy()
    cmfe_coord.Destroy()
    cmfe_region.Destroy()
    cmfe_quad_basis.Destroy()
    iron.Finalise()

# +==+ ^\_/^ +==+ ^\_/^ +==+ 
# Run check
# +==+ ^\_/^ +==+ ^\_/^ +==+ 

if __name__ == '__main__':
    
    # +==+ Intake Arguments
    argparser = argparse.ArgumentParser("A python program to run inflation tests in OpenCMISS Iron.")
    argparser.add_argument("test_file")
    argparser.add_argument("test_type")
    argparser.add_argument("pressure_test")
    args = argparser.parse_args()
    test_name = args.test_file
    type_name = args.test_type
    pressure_test = args.pressure_test

    if pressure_test == "False":
        main(test_name, type_name, False)
    else:
        main(test_name, type_name, True)