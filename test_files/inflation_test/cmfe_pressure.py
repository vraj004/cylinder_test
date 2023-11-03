"""
Author: Liam Murray, murrayla@student.unimelb.edu.au
Descrption: cmiss setup funtions for finite element implementation of 
                cylinder inflation
"""

import numpy as np
from opencmiss.iron import iron
import meshio

# Runtime required parameters
X, Y, Z, P = (1, 2, 3, 4)

IRON_VTK = [
    0, 2, 8, 6, 
    18, 20, 26, 24,
    1, 5, 7, 3, 
    19, 13, 25, 21,
    9, 11, 17, 15, 
    12, 14, 10, 16,
    4, 22, 13
]

GMSH2VTK = [
    0, 1, 2, 3, 4, 5, 6, 7,
    8, 11, 13, 9, 16, 18, 19, 17,
    10, 12, 14, 15, 22, 23, 21, 24,
    20, 25, 26
]

GMSH2IRON_LIN = [
    0, 1, 3, 2,
    4, 5, 7, 6
]
IRON2VTK_LIN = [
    0, 1, 3, 2,
    4, 5, 7, 6
]


EXTENSION_TEST = False
PRESSURE_TEST = True
INNER_RAD = 1
OUTER_RAD = 1.5
TOL_VEC = 1e-7

# +==+ ^\_/^ +==+ ^\_/^ +==+ 
# Coordinate System:
#   Simply defines the coordinate system that is used for 
#       solving the finite element problem (i.e. 3D).
# +==+ ^\_/^ +==+ ^\_/^ +==+ 

def coordinate_setup(coord_n, dim):
    coord_sys = iron.CoordinateSystem()
    coord_sys.CreateStart(coord_n)
    coord_sys.DimensionSet(dim)
    coord_sys.CreateFinish()
    return coord_sys

# +==+ ^\_/^ +==+ ^\_/^ +==+ 
# Basis:
#   Denotes the interpolation of the problem geometry. 
#       The type of basis here is LAGRANGE_HERMITE_TP which
#       is describing the form of the interpolation functions. This 
#       inidicates the sort of continuity expected between each node.  
#       The Interpolation Scheme is denoted here as QUADRATIC_LAGRANGE, 
#       meaning there are 3 nodes along each direction i.e.,
#       -x <- n1 ---- n2 ---- n3 -> x
# +==+ ^\_/^ +==+ ^\_/^ +==+ 

def basis_setup(basis_n, xi_n, type):

    if type == "quadratic":
        basis = iron.Basis()
        basis.CreateStart(basis_n)
        basis.TypeSet(iron.BasisTypes.LAGRANGE_HERMITE_TP)
        basis.NumberOfXiSet(xi_n)
        basis.InterpolationXiSet(
            [iron.BasisInterpolationSpecifications.QUADRATIC_LAGRANGE] * xi_n
        )
        basis.QuadratureNumberOfGaussXiSet([3]*xi_n)
        basis.CreateFinish()

    elif type == "linear":
        basis = iron.Basis()
        basis.CreateStart(basis_n)
        basis.TypeSet(iron.BasisTypes.LAGRANGE_HERMITE_TP)
        basis.NumberOfXiSet(xi_n)
        basis.InterpolationXiSet(
            [iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE] * xi_n
        )
        basis.QuadratureNumberOfGaussXiSet([3]*xi_n)
        basis.CreateFinish()

    return basis

# +==+ ^\_/^ +==+ ^\_/^ +==+ 
# Region: 
#   The region is where a certain physics is occuring, we can always
#       specify multiple regions to couple different operations and
#       seperate them in our computations.
# +==+ ^\_/^ +==+ ^\_/^ +==+ 

def region_setup(region_n, coord_sys):
    region = iron.Region()
    region.CreateStart(region_n, iron.WorldRegion)
    region.CoordinateSystemSet(coord_sys)
    region.LabelSet("Region")
    region.CreateFinish()
    return region

# +==+ ^\_/^ +==+ ^\_/^ +==+ 
# Decomposition:
#   This step is implementing parallel processing by
#       breaking the mesh into different components for
#       calculation. The Domains are each a parallel stream
#       which is set by "comp_nodes_n"; the available 
#       quantity of nodes. The parameter "CALCULATED" indicates
#       that OpenCMISS will calculate the decomposition with
#       a graph partitioning approach.
# +==+ ^\_/^ +==+ ^\_/^ +==+ 

def decomposition_setup(mesh, decomp_n):
    comp_nodes_n = iron.ComputationalNumberOfNodesGet()
    decomp = iron.Decomposition()
    decomp.CreateStart(decomp_n, mesh)
    decomp.TypeSet(iron.DecompositionTypes.CALCULATED)
    decomp.NumberOfDomainsSet(comp_nodes_n)
    decomp.CreateFinish()
    return decomp

# +==+ ^\_/^ +==+ ^\_/^ +==+ 
# Geometric Field:
#   Once a decomposition into domains has been set up, these domains
#       can have parts of the geometry allocated to them. This is the
#       geometric field. Variables in this field and others can have
#       labels associated with them. Field values are by default 0
#       so after creating the field we can then allocate them values. 
# +==+ ^\_/^ +==+ ^\_/^ +==+ 

def geometric_setup(geo_field_n, region, decomp, n_idx, n_xyz):
    geo_field = iron.Field()
    geo_field.CreateStart(geo_field_n, region) 
    geo_field.MeshDecompositionSet(decomp)
    geo_field.TypeSet(iron.FieldTypes.GEOMETRIC)
    geo_field.VariableLabelSet(iron.FieldVariableTypes.U, "Geometry")
    geo_field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,X,1)
    geo_field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,Y,1) 
    geo_field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,Z,1) 
    geo_field.CreateFinish()
    return geo_field

# +==+ ^\_/^ +==+ ^\_/^ +==+ 
# Material Field:
#   The next field is the material field. This field allows
#       us to dictate the constitutive properties of the 
#       material being tested. Here we intend to use a Mooney-Rivlin
#       material so this requires two constitutive components c1 and c2.
#       Here we tell the componenets to be ELEMENT_BASED which indicates
#       that the parameters can be different between elements but not
#       between nodes. 
# +==+ ^\_/^ +==+ ^\_/^ +==+ 

def material_setup(mat_field_n, decomp, geo_field, region, c):
    mat_field = iron.Field()
    mat_field.CreateStart(mat_field_n, region)
    mat_field.TypeSet(iron.FieldTypes.MATERIAL)
    mat_field.MeshDecompositionSet(decomp)
    mat_field.GeometricFieldSet(geo_field)
    mat_field.VariableLabelSet(iron.FieldVariableTypes.U, "Material")
    mat_field.NumberOfComponentsSet(iron.FieldVariableTypes.U, len(c))
    mat_field.ScalingTypeSet(iron.FieldScalingTypes.ARITHMETIC_MEAN)
    for i, c_val in enumerate(c):
        mat_field.ComponentInterpolationSet(
            iron.FieldVariableTypes.U, i+1, 
            iron.FieldInterpolationTypes.GAUSS_POINT_BASED
        )
    mat_field.CreateFinish()
    for i, c_val in enumerate(c):
        mat_field.ComponentValuesInitialiseDP(
            iron.FieldVariableTypes.U, 
            iron.FieldParameterSetTypes.VALUES, 
            i+1, c_val
        )
    return mat_field

# +==+ ^\_/^ +==+ ^\_/^ +==+ 
# Dependent Field:
#   To hold the values of our nodes after deformation we require a 
#       dependent field. To intialise this before solving we first
#       set the values to be that of the undeformed geometry from
#       the geometric field. And add a fourth component which is 
#       indicating the hydrostatic pressure. Here, we denote 
#       ELEMENT_BASED for pressure, meaning it is per element. 
# +==+ ^\_/^ +==+ ^\_/^ +==+ 

def dependent_setup(dep_field_n, region, decomp, geo_field, pressure_test):
    dep_field = iron.Field()
    dep_field.CreateStart(dep_field_n, region)
    dep_field.MeshDecompositionSet(decomp)
    dep_field.TypeSet(iron.FieldTypes.GEOMETRIC_GENERAL)
    dep_field.GeometricFieldSet(geo_field)
    dep_field.DependentTypeSet(iron.FieldDependentTypes.DEPENDENT)
    dep_field.VariableLabelSet(iron.FieldVariableTypes.U, "Dependent")
    dep_field.NumberOfVariablesSet(2)
    dep_field.NumberOfComponentsSet(iron.FieldVariableTypes.U, 4)
    dep_field.NumberOfComponentsSet(iron.FieldVariableTypes.DELUDELN, 4)

    if pressure_test:
        dep_field.ComponentMeshComponentSet(iron.FieldVariableTypes.U, X, 1)
        dep_field.ComponentMeshComponentSet(iron.FieldVariableTypes.U, Y, 1) 
        dep_field.ComponentMeshComponentSet(iron.FieldVariableTypes.U, Z, 1) 
        dep_field.ComponentMeshComponentSet(iron.FieldVariableTypes.U, P, 2)
        dep_field.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN, P, 2)
        dep_field.ComponentInterpolationSet(
            iron.FieldVariableTypes.U, 4, iron.FieldInterpolationTypes.NODE_BASED
        )
        dep_field.ComponentInterpolationSet(
            iron.FieldVariableTypes.DELUDELN, 4, iron.FieldInterpolationTypes.NODE_BASED
        )
        dep_field.ScalingTypeSet(iron.FieldScalingTypes.ARITHMETIC_MEAN) # ??

    else:
        dep_field.ComponentInterpolationSet(
            iron.FieldVariableTypes.U, 4, iron.FieldInterpolationTypes.ELEMENT_BASED
        )
        dep_field.ComponentInterpolationSet(
            iron.FieldVariableTypes.DELUDELN, 4, iron.FieldInterpolationTypes.ELEMENT_BASED
        )

    dep_field.CreateFinish()
    
    return dep_field

# +==+ ^\_/^ +==+ ^\_/^ +==+ 
# Equations Field:
#    This field allows us to indicate what type of equations are used 
#       to solve the problem set. ELASTICITY indicates that the geometry
#       is being tested for an elastic case; FINITE_ELASTICITY is the 
#       the type of the FEM equation being solved. We then classify the 
#       constitutive equation for the material properties defined
#       earlier, MOONEY_RIVLIN. We can then explicitly link the other
#       fields and define output and sparsity types for solving.
# +==+ ^\_/^ +==+ ^\_/^ +==+ 

def equations_setup(eqs_set):
    eqs = iron.Equations()
    eqs_set.EquationsCreateStart(eqs)
    eqs.SparsityTypeSet(iron.EquationsSparsityTypes.SPARSE)
    eqs.OutputTypeSet(iron.EquationsOutputTypes.NONE)
    eqs_set.EquationsCreateFinish()
    return 

# +==+ ^\_/^ +==+ ^\_/^ +==+ 
# Problem Solver:
#    
# +==+ ^\_/^ +==+ ^\_/^ +==+ 

def problem_solver_setup(problem_n, eqs_set, loadsteps):
    problem = iron.Problem()
    problems_specs = (
        [
            iron.ProblemClasses.ELASTICITY,
            iron.ProblemTypes.FINITE_ELASTICITY,
            iron.ProblemSubtypes.NONE
        ]
    )
    problem.CreateStart(problem_n, problems_specs)
    problem.CreateFinish()
    problem.ControlLoopCreateStart()
    ctrl_loop = iron.ControlLoop()
    problem.ControlLoopGet(
        [iron.ControlLoopIdentifiers.NODE], ctrl_loop
    )
    ctrl_loop.MaximumIterationsSet(loadsteps)
    problem.ControlLoopCreateFinish()
    non_solver = iron.Solver()
    lin_solver = iron.Solver()
    problem.SolversCreateStart()
    problem.SolverGet(
        [iron.ControlLoopIdentifiers.NODE], 1, non_solver
    )
    non_solver.OutputTypeSet(iron.SolverOutputTypes.PROGRESS)
    non_solver.NewtonJacobianCalculationTypeSet(iron.JacobianCalculationTypes.EQUATIONS)
    # non_solver.NewtonAbsoluteToleranceSet(1e-3)
    # non_solver.NewtonSolutionToleranceSet(1e-2)
    # non_solver.NewtonConvergenceTestTypeSet(iron.NewtonConvergenceTypes.PETSC_DEFAULT)
    non_solver.NewtonLinearSolverGet(lin_solver)
    # non_solver.NewtonLineSearchTypeSet(iron.NewtonLineSearchTypes.QUADRATIC)
    lin_solver.LinearTypeSet(iron.LinearSolverTypes.DIRECT)
    # lin_solver.LibraryTypeSet(iron.SolverLibraries.MUMPS)
    problem.SolversCreateFinish()
    solver = iron.Solver()
    solver_eqs = iron.SolverEquations()
    problem.SolverEquationsCreateStart()
    problem.SolverGet(
        [iron.ControlLoopIdentifiers.NODE], 1, solver
    )
    solver.SolverEquationsGet(solver_eqs)
    solver_eqs.SparsityTypeSet(iron.SolverEquationsSparsityTypes.SPARSE)
    _ = solver_eqs.EquationsSetAdd(eqs_set)
    problem.SolverEquationsCreateFinish()
    return problem, solver, solver_eqs

def boundary_conditions_setup(solver_eqs, dep_field, n_n, n_np_xyz, pre, tra):
    bcs = iron.BoundaryConditions()
    solver_eqs.BoundaryConditionsCreateStart(bcs)
    min_z = np.min(n_np_xyz[:, 2])
    max_z = np.max(n_np_xyz[:, 2])
    bounds = [0, 1]
    # for i in range(0, n_n, 1):
#         x = n_np_xyz[i, 0]
#         y = n_np_xyz[i, 1]
#         z = n_np_xyz[i, 2]
#         # +==+ EXTENSION OPTION
#         if EXTENSION_TEST:
#             if n_np_xyz[i, 2] == max_z: 
#                 bcs.AddNode(
#                     dep_field, 
#                     iron.FieldVariableTypes.U,
#                     1, 1, i+1, Z,
#                     iron.BoundaryConditionsTypes.FIXED,
#                     tra
#                 )
#                 # bcs.AddNode(dep_field, 
#                 #     iron.FieldVariableTypes.DELUDELN, 1, 1, i+1, Z,
#                 #     iron.BoundaryConditionsTypes.PRESSURE_INCREMENTED, 
#                 #     pre
#                 # )
#             else: 
#                 if n_np_xyz[i, 2] == min_z: 
#                     bcs.AddNode(
#                     dep_field, 
#                     iron.FieldVariableTypes.U,
#                     1, 1, i+1, Z,
#                     iron.BoundaryConditionsTypes.FIXED,
#                     0.0
#                 )
#                 for j in [X, Y]:
#                     bcs.AddNode(
#                         dep_field, 
#                         iron.FieldVariableTypes.U,
#                         1, 1, i+1, j,
#                         iron.BoundaryConditionsTypes.FIXED,
#                         0.0
#                     )
#         # +==+ PRESSURE OPTION
#         if PRESSURE_TEST:
#             # vec_norm = np.linalg.norm([n_np_xyz[i, 0], n_np_xyz[i, 1]])
#             vec = np.array([x, y])
#             vec_norm = np.sqrt(vec.dot(vec))
#             conditions = [
#                 ((x >= bounds[0] and x <= bounds[1]) and y == bounds[0]),
#                 ((x >= bounds[0] and x <= bounds[1]) and y == bounds[1]),
#                 ((y >= bounds[0] and y <= bounds[1]) and x == bounds[0]),
#                 ((y >= bounds[0] and y <= bounds[1]) and x == bounds[1]),
#             ]
# #             if np.any(conditions):
# #                 bcs.AddNode(dep_field, 
# #                     iron.FieldVariableTypes.DELUDELN, 1, 1, i+1, Z,
# #                     iron.BoundaryConditionsTypes.PRESSURE_INCREMENTED, 
# #                     pre
# #                 )
# #                 print(i, conditions, x, y, z)
#             # if i == 102: #x == 0 and (y == 0.25): # or y == 0.5 or y == 0.75):
#             #     bcs.AddNode(dep_field, 
#             #         iron.FieldVariableTypes.DELUDELN, 1, 1, i+1, Z,
#             #         iron.BoundaryConditionsTypes.PRESSURE_INCREMENTED, 
#             #         pre
#             #     )
#             # if z == min_z or x == 1:
#             #     for j in [X, Y, Z]:
#             #         bcs.AddNode(
#             #             dep_field, 
#             #             iron.FieldVariableTypes.U,
#             #             1, 1, i+1, j,
#             #             iron.BoundaryConditionsTypes.FIXED,
#             #             0.0
#             #         )
#             # if z == max_z:
#             #     for j in [Z]:
#             #         bcs.AddNode(
#             #             dep_field, 
#             #             iron.FieldVariableTypes.U,
#             #             1, 1, i+1, j,
#             #             iron.BoundaryConditionsTypes.FIXED,
#             #             0.0
#             #         )
# #             if np.isclose(vec_norm, OUTER_RAD, 1e-5) and np.isclose(np.min(n_np_xyz[:, 2]), n_np_xyz[i, 2], 1e-3): 
# #                 for j in [Z]:
# #                     bcs.AddNode(
# #                         dep_field, 
# #                         iron.FieldVariableTypes.U,
# #                         1, 1, i+1, j,
# #                         iron.BoundaryConditionsTypes.FIXED,
# #                         0.0
# #                     )
# #             if np.isclose(vec_norm, INNER_RAD, 1e-5) or np.isclose(vec_norm, OUTER_RAD, 1e-5):
# #                 bcs.AddNode(
# #                         dep_field, 
# #                         iron.FieldVariableTypes.U,
# #                         1, 1, i+1, X,
# #                         iron.BoundaryConditionsTypes.FIXED,
# #                         0.0005 * n_np_xyz[i, 0]
# #                     )
# #                 bcs.AddNode(
# #                         dep_field, 
# #                         iron.FieldVariableTypes.U,
# #                         1, 1, i+1, Y,
# #                         iron.BoundaryConditionsTypes.FIXED,
# #                         0.0005 * n_np_xyz[i, 1]
# #                     )
#                 # bcs.AddNode(dep_field, 
#                 #     iron.FieldVariableTypes.DELUDELN, 1, 1, i+1, Z,
#                 #     iron.BoundaryConditionsTypes.PRESSURE_INCREMENTED, 
#                 #     pre
#                 # )

    solver_eqs.BoundaryConditionsCreateFinish()
    return 

# +==+ ^\_/^ +==+ ^\_/^ +==+ 
# Deformed Field:
#    
# +==+ ^\_/^ +==+ ^\_/^ +==+ 

def deformed_setup(def_field_n, region, decomp, dep_field):
    def_field = iron.Field()
    def_field.CreateStart(def_field_n, region)
    def_field.MeshDecompositionSet(decomp)
    def_field.TypeSet(iron.FieldTypes.GEOMETRIC)
    def_field.VariableLabelSet(iron.FieldVariableTypes.U, "DeformedGeometry")
    def_field.ComponentMeshComponentSet(iron.FieldVariableTypes.U, X, 1)
    def_field.ComponentMeshComponentSet(iron.FieldVariableTypes.U, Y, 1)
    def_field.ComponentMeshComponentSet(iron.FieldVariableTypes.U, Z, 1)
    def_field.CreateFinish()
    return dep_field

# +==+ ^\_/^ +==+ ^\_/^ +==+ 
# Pressure Field:
#    
# +==+ ^\_/^ +==+ ^\_/^ +==+ 

def pressure_setup(region, decomp, pre_field_n, pressure_test):
    pre_field = iron.Field()
    pre_field.CreateStart(pre_field_n, region)
    pre_field.MeshDecompositionSet(decomp)
    pre_field.VariableLabelSet(iron.FieldVariableTypes.U, "Pressure")

    if pressure_test:
        pre_field.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 1, 2)
        pre_field.ComponentInterpolationSet(
            iron.FieldVariableTypes.U, 
            1, 
            iron.FieldInterpolationTypes.NODE_BASED
        )
        pre_field.NumberOfComponentsSet(iron.FieldVariableTypes.U, 1)

    else:
        pre_field.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 1, 1)
        pre_field.ComponentInterpolationSet(
            iron.FieldVariableTypes.U, 
            1, 
            iron.FieldInterpolationTypes.ELEMENT_BASED
        )
        pre_field.NumberOfComponentsSet(iron.FieldVariableTypes.U, 1)

    pre_field.CreateFinish()
    return pre_field

# +==+ ^\_/^ +==+ ^\_/^ +==+ 
# Meshio VTK Output:
#    
# +==+ ^\_/^ +==+ ^\_/^ +==+ 

def vtk_output(mesh, n_n, geo_field, dep_field, e_np_map, mesh_e, runtime_path, test_name, convention, test_type):
    # +============+ 
    # Store nodes Before & After deformation
    # +============+ 
    meshNodes = iron.MeshNodes()
    mesh.NodesGet(1,meshNodes)
    n_list = []
    n_bef = []
    n_aft = []
    for i in range(0, n_n, 1):
        n_bef.append(
            [
                geo_field.ParameterSetGetNodeDP(
                    iron.FieldVariableTypes.U,
                    iron.FieldParameterSetTypes.VALUES,
                    1, 1, i+1, X
                ),
                geo_field.ParameterSetGetNodeDP(
                    iron.FieldVariableTypes.U,
                    iron.FieldParameterSetTypes.VALUES,
                    1, 1, i+1, Y
                ),
                geo_field.ParameterSetGetNodeDP(
                    iron.FieldVariableTypes.U,
                    iron.FieldParameterSetTypes.VALUES,
                    1, 1, i+1, Z
                )
            ]
        )
        n_aft.append(
            [
                dep_field.ParameterSetGetNodeDP(
                    iron.FieldVariableTypes.U,
                    iron.FieldParameterSetTypes.VALUES,
                    1, 1, i+1, X
                ),
                dep_field.ParameterSetGetNodeDP(
                    iron.FieldVariableTypes.U,
                    iron.FieldParameterSetTypes.VALUES,
                    1, 1, i+1, Y
                ),
                dep_field.ParameterSetGetNodeDP(
                    iron.FieldVariableTypes.U,
                    iron.FieldParameterSetTypes.VALUES,
                    1, 1, i+1, Z
                )
            ]
        )
        n_list.append(
            [
                i+1, 
                n_bef[i][0], n_bef[i][1], n_bef[i][2], 
                n_aft[i][0], n_aft[i][1], n_aft[i][2]
            ]
        )
    # +============+ 
    # Store elements
    # +============+ 
    e_list = e_np_map
    # +============+ 
    # Store data files
    # +============+
    node_file = open(runtime_path + 'output_mesh.node', 'w')
    node_file.writelines([str(line) + "\n" for line in n_list])
    node_file.close()
    elem_file=open(runtime_path + 'output_mesh.ele','w')
    elem_file.writelines([str(line) + "\n" for line in e_list])
    elem_file.close()
    bef_def = np.array(n_bef)
    aft_def = np.array(n_aft)
    np.save(runtime_path + 'output_mesh_before_def.npy',bef_def)
    np.save(runtime_path + 'output_mesh_after_def.npy',aft_def)
    # +============+ 
    # VTK export
    # +============+
    e_list_gmsh = np.array(e_list)[:,:] - 1
    e_list_vtk = e_list_gmsh[:, GMSH2VTK]
    meshio.write_points_cells(
        "vtk_files/" + test_name + test_type + ".vtk", 
        bef_def, 
        [("hexahedron27", e_list_vtk)] + [("hexahedron27", e_list_vtk)], 
        {"deformed": aft_def-bef_def}
    )
    return 