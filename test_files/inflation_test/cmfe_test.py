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
import meshio

# +==+ ^\_/^ +==+ ^\_/^ +==+ 
# Node numbering convension
# +==+ ^\_/^ +==+ ^\_/^ +==+ 

#         VTK            Gmsh
#               top
#   *  7--14--6     *  7--19--6
#   *  |      |     *  |      |
#   * 15  25  13    * 17  25  18
#   *  |      |     *  |      |
#   *  4--12--5     *  4--16--5
#   *
#   *           middle
#   * 19--23--18    * 15--24--14
#   *  |      |     *  |      |
#   * 20  26  21    * 22  26  23
#   *  |      |     *  |      |
#   * 16--22--17    * 10--21--12
#   *
#   *           bottom
#   *  3--10--2     *  3--13--2
#   *  |      |     *  |      |
#   * 11  24  9     *  9  20  11
#   *  |      |     *  |      |
#   *  0-- 8--1     *  0-- 8--1

# +==+ ^\_/^ +==+ ^\_/^ +==+ 
# Parameter Setup
# +==+ ^\_/^ +==+ ^\_/^ +==+ 

DIM = 3
XI_N = 3
N_N_EL = 27
QUAD_ORDER = 4
X, Y, Z, P = (1, 2, 3, 4)
PRESSURE_TEST = True
LOADSTEPS = 5
INNER_RAD = 0.375
C_VALS = [1, 0.2]
RUNTIME_PATH = "/home/jovyan/work/docker-iron/test_files/inflation_test/runtime_files/"
GMSH2VTK = [
    0, 1, 2, 3, 4, 5, 6, 7,
    8, 11, 13, 9, 16, 18, 19, 17,
    10, 12, 14, 15, 22, 23, 21, 24,
    20, 25, 26
]
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
# Meshio VTK Output:
#    
# +==+ ^\_/^ +==+ ^\_/^ +==+ 

def vtk_output(mesh, n_n, geo_field, dep_field, e_np_map, mesh_e, runtime_path, test_name):
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
        runtime_path + test_name + ".vtk", 
        bef_def, 
        [("hexahedron27", e_list_vtk)] + [("hexahedron27", e_list_vtk)], 
        {"deformed": aft_def}
    )
    return

# +==+ ^\_/^ +==+ ^\_/^ +==+ 
# Main function for safe operation of inflation test
# +==+ ^\_/^ +==+ ^\_/^ +==+ 

def main(test_name):

    # +==+ ^\_/^ +==+ ^\_/^ +==+ 
    # Coordinate System:
    #   Simply defines the coordinate system that is used for 
    #       solving the finite element problem (i.e. 3D).
    # +==+ ^\_/^ +==+ ^\_/^ +==+ 

    cmfe_coord = iron.CoordinateSystem()
    cmfe_coord.CreateStart(coord_n)
    cmfe_coord.DimensionSet(DIM)
    cmfe_coord.CreateFinish()
    print('+==+ COORDINATE SYSTEM COMPLETE')
    
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

    cmfe_basis = iron.Basis()
    cmfe_basis.CreateStart(basis_n)
    cmfe_basis.TypeSet(iron.BasisTypes.LAGRANGE_HERMITE_TP)
    cmfe_basis.NumberOfXiSet(XI_N)
    cmfe_basis.InterpolationXiSet(
        [iron.BasisInterpolationSpecifications.QUADRATIC_LAGRANGE] * XI_N
    )
    cmfe_basis.QuadratureNumberOfGaussXiSet([3]*XI_N)
    cmfe_basis.CreateFinish()
    print('+==+ BASIS SYSTEM COMPLETE')

    # +==+ ^\_/^ +==+ ^\_/^ +==+ 
    # Region: 
    #   The region is where a certain physics is occuring, we can always
    #       specify multiple regions to couple different operations and
    #       seperate them in our computations.
    # +==+ ^\_/^ +==+ ^\_/^ +==+ 

    cmfe_region = iron.Region()
    cmfe_region.CreateStart(region_n, iron.WorldRegion)
    cmfe_region.CoordinateSystemSet(cmfe_coord)
    cmfe_region.LabelSet("Region")
    cmfe_region.CreateFinish()
    print('+==+ REGION COMPLETE')

    # +==+ ^\_/^ +==+ ^\_/^ +==+ 
    # Mesh: 
    #   
    # +==+ ^\_/^ +==+ ^\_/^ +==+ 

    # +==+ loading data from functions
    n_np_xyz, n_idx, n_n = nodes(test_name)
    e_np_map, e_idx, e_n = elems(test_name)
    # +==+ cmfe nodes
    cmfe_node = iron.Nodes()
    cmfe_node.CreateStart(cmfe_region, n_n)
    cmfe_node.CreateFinish()
    # +==+ cmfe mesh
    cmfe_mesh = iron.Mesh()
    cmfe_mesh.CreateStart(mesh_n, cmfe_region, DIM)
    cmfe_mesh.NumberOfElementsSet(e_n)
    cmfe_mesh.NumberOfComponentsSet(1) 
    # +==+ cmfe mesh elements
    cmfe_mesh_e = iron.MeshElements()
    cmfe_mesh_e.CreateStart(cmfe_mesh, 1, cmfe_basis)
    # += allocating nodes to elements
    print('+= ... begin mesh allocation')
    for i in range(e_n):
        nodesList = list(
            map(int,e_np_map[i][:])
        )
        cmfe_mesh_e.NodesSet(e_idx[i], nodesList)
    # +=
    cmfe_mesh_e.CreateFinish()
    cmfe_mesh.CreateFinish()
    print('+==+ MESH ALLOCATION COMPLETE')
    
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

    comp_nodes_n = iron.ComputationalNumberOfNodesGet()
    cmfe_decomp = iron.Decomposition()
    cmfe_decomp.CreateStart(decomp_n, cmfe_mesh)
    cmfe_decomp.TypeSet(iron.DecompositionTypes.CALCULATED)
    cmfe_decomp.NumberOfDomainsSet(comp_nodes_n)
    cmfe_decomp.CreateFinish()
    print('+==+ DECOMPOSITION FIELD COMPLETE')

    # +==+ ^\_/^ +==+ ^\_/^ +==+ 
    # Geometric Field:
    #   Once a decomposition into domains has been set up, these domains
    #       can have parts of the geometry allocated to them. This is the
    #       geometric field. Variables in this field and others can have
    #       labels associated with them. Field values are by default 0
    #       so after creating the field we can then allocate them values. 
    # +==+ ^\_/^ +==+ ^\_/^ +==+ 

    cmfe_geo_field = iron.Field()
    cmfe_geo_field.CreateStart(geo_field_n, cmfe_region) 
    cmfe_geo_field.MeshDecompositionSet(cmfe_decomp)
    cmfe_geo_field.TypeSet(iron.FieldTypes.GEOMETRIC)
    cmfe_geo_field.VariableLabelSet(iron.FieldVariableTypes.U, "Geometry")
    for i in [X, Y, Z]:
        cmfe_geo_field.ComponentMeshComponentSet(iron.FieldVariableTypes.U, i, 1)
    cmfe_geo_field.CreateFinish()
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

    cmfe_mat_field = iron.Field()
    cmfe_mat_field.CreateStart(mat_field_n, cmfe_region)
    cmfe_mat_field.TypeSet(iron.FieldTypes.MATERIAL)
    cmfe_mat_field.MeshDecompositionSet(cmfe_decomp)
    cmfe_mat_field.GeometricFieldSet(cmfe_geo_field)
    cmfe_mat_field.VariableLabelSet(iron.FieldVariableTypes.U, "Material")
    cmfe_mat_field.NumberOfComponentsSet(iron.FieldVariableTypes.U, len(C_VALS))
    cmfe_mat_field.ScalingTypeSet(iron.FieldScalingTypes.ARITHMETIC_MEAN)
    for i, c in enumerate(C_VALS):
        cmfe_mat_field.ComponentInterpolationSet(
            iron.FieldVariableTypes.U, i+1, 
            iron.FieldInterpolationTypes.GAUSS_POINT_BASED
        )
    cmfe_mat_field.CreateFinish()
    for i, c in enumerate(C_VALS):
        cmfe_mat_field.ComponentValuesInitialiseDP(
            iron.FieldVariableTypes.U, 
            iron.FieldParameterSetTypes.VALUES, 
            i+1, c
        )
    print('+==+ MATERIAL FIELD COMPLETE')

    # +==+ ^\_/^ +==+ ^\_/^ +==+ 
    # Dependent Field:
    #   To hold the values of our nodes after deformation we require a 
    #       dependent field. To intialise this before solving we first
    #       set the values to be that of the undeformed geometry from
    #       the geometric field. And add a fourth component which is 
    #       indicating the hydrostatic pressure. Here, we denote 
    #       ELEMENT_BASED for pressure, meaning it is per element. 
    # +==+ ^\_/^ +==+ ^\_/^ +==+ 

    cmfe_dep_field = iron.Field()
    cmfe_dep_field.CreateStart(dep_field_n, cmfe_region)
    cmfe_dep_field.MeshDecompositionSet(cmfe_decomp)
    cmfe_dep_field.TypeSet(iron.FieldTypes.GEOMETRIC_GENERAL)
    cmfe_dep_field.GeometricFieldSet(cmfe_geo_field)
    cmfe_dep_field.DependentTypeSet(iron.FieldDependentTypes.DEPENDENT)
    cmfe_dep_field.VariableLabelSet(iron.FieldVariableTypes.U, "Dependent")
    cmfe_dep_field.NumberOfVariablesSet(2)
    cmfe_dep_field.NumberOfComponentsSet(iron.FieldVariableTypes.U, 4)
    cmfe_dep_field.NumberOfComponentsSet(iron.FieldVariableTypes.DELUDELN, 4)
    cmfe_dep_field.ComponentInterpolationSet(
        iron.FieldVariableTypes.U, 4,
        iron.FieldInterpolationTypes.ELEMENT_BASED)
    cmfe_dep_field.ComponentInterpolationSet(
        iron.FieldVariableTypes.DELUDELN, 4,
        iron.FieldInterpolationTypes.ELEMENT_BASED)
    cmfe_dep_field.CreateFinish()
    for i in [X, Y, Z]:
        iron.Field.ParametersToFieldParametersComponentCopy(
            cmfe_geo_field, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, i,
            cmfe_dep_field, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, i
        )
    iron.Field.ComponentValuesInitialiseDP(
        cmfe_dep_field, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, P, 0
    )
    print('+==+ DEPENDENT FIELD COMPLETE')

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
    cmfe_eqs = iron.Equations()
    cmfe_eqs_set.EquationsCreateStart(cmfe_eqs)
    cmfe_eqs.SparsityTypeSet(iron.EquationsSparsityTypes.SPARSE)
    cmfe_eqs.OutputTypeSet(iron.EquationsOutputTypes.NONE)
    cmfe_eqs_set.EquationsCreateFinish()
    print('+==+ EQUATION FIELD COMPLETE')

    # +==+ ^\_/^ +==+ ^\_/^ +==+ 
    # Problem Solver:
    #    
    # +==+ ^\_/^ +==+ ^\_/^ +==+ 
        
    # +==+ Export field information so far
    fields = iron.Fields()
    fields.CreateRegion(cmfe_region)
    fields.NodesExport("Output", "FORTRAN")
    fields.ElementsExport("Output", "FORTRAN")
    fields.Finalise()

    # += iterations through increments for solution
    pre_inc = [15000/LOADSTEPS] * LOADSTEPS
    print('+= ... begin solver')
    for i, inc in enumerate(range(0, len(pre_inc))):

        # +============+ 
        # Problem and Solution infrastructure
        # +============+ 
        
        # +==+ cmfe problem and solver field
        cmfe_problem = iron.Problem()
        cmfe_problems_specs = (
            [
                iron.ProblemClasses.ELASTICITY,
                iron.ProblemTypes.FINITE_ELASTICITY,
                iron.ProblemSubtypes.NONE
            ]
        )
        cmfe_problem.CreateStart(problem_n, cmfe_problems_specs)
        cmfe_problem.CreateFinish()
        cmfe_problem.ControlLoopCreateStart()
        cmfe_ctrl_loop = iron.ControlLoop()
        cmfe_problem.ControlLoopGet(
            [iron.ControlLoopIdentifiers.NODE], cmfe_ctrl_loop
        )
        cmfe_ctrl_loop.MaximumIterationsSet(LOADSTEPS)
        cmfe_problem.ControlLoopCreateFinish()
        cmfe_non_solver = iron.Solver()
        cmfe_lin_solver = iron.Solver()
        cmfe_problem.SolversCreateStart()
        cmfe_problem.SolverGet(
            [iron.ControlLoopIdentifiers.NODE], 1, cmfe_non_solver
        )
        cmfe_non_solver.OutputTypeSet(iron.SolverOutputTypes.PROGRESS)
        cmfe_non_solver.NewtonJacobianCalculationTypeSet(iron.JacobianCalculationTypes.EQUATIONS)
        # non_solver.NewtonAbsoluteToleranceSet(1e-3)
        # non_solver.NewtonSolutionToleranceSet(1e-2)
        # non_solver.NewtonConvergenceTestTypeSet(iron.NewtonConvergenceTypes.PETSC_DEFAULT)
        cmfe_non_solver.NewtonLinearSolverGet(cmfe_lin_solver)
        # non_solver.NewtonLineSearchTypeSet(iron.NewtonLineSearchTypes.QUADRATIC)
        cmfe_lin_solver.LinearTypeSet(iron.LinearSolverTypes.DIRECT)
        # lin_solver.LibraryTypeSet(iron.SolverLibraries.MUMPS)
        cmfe_problem.SolversCreateFinish()
        cmfe_solver = iron.Solver()
        cmfe_solver_eqs = iron.SolverEquations()
        cmfe_problem.SolverEquationsCreateStart()
        cmfe_problem.SolverGet(
            [iron.ControlLoopIdentifiers.NODE], 1, cmfe_solver
        )
        cmfe_solver.SolverEquationsGet(cmfe_solver_eqs)
        cmfe_solver_eqs.SparsityTypeSet(iron.SolverEquationsSparsityTypes.SPARSE)
        _ = cmfe_solver_eqs.EquationsSetAdd(cmfe_eqs_set)
        cmfe_problem.SolverEquationsCreateFinish()

        cmfe_bcs = iron.BoundaryConditions()
        cmfe_solver_eqs.BoundaryConditionsCreateStart(cmfe_bcs)
        min_z = np.min(n_np_xyz[:, 2])
        max_z = np.max(n_np_xyz[:, 2])
        for i in range(0, n_n, 1):
            if n_np_xyz[i, 0] == max_z: 
                cmfe_bcs.AddNode(
                    cmfe_dep_field, 
                    iron.FieldVariableTypes.U,
                    1, 1, i+1, X,
                    iron.BoundaryConditionsTypes.FIXED,
                    0.02
                )
            else: 
                if n_np_xyz[i, 0] == max_z: 
                    cmfe_bcs.AddNode(
                    cmfe_dep_field, 
                    iron.FieldVariableTypes.U,
                    1, 1, i+1, X,
                    iron.BoundaryConditionsTypes.FIXED,
                    0.0
                )
                for j in [Y, Z]:
                    cmfe_bcs.AddNode(
                    cmfe_dep_field, 
                    iron.FieldVariableTypes.U,
                    1, 1, i+1, j,
                    iron.BoundaryConditionsTypes.FIXED,
                    0.0
                )
                # bcs.AddNode(dep_field, 
                #             iron.FieldVariableTypes.DELUDELN, 1,
                #             iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 
                #             i+1, Z,
                #             iron.BoundaryConditionsTypes.PRESSURE_INCREMENTED, 
                #             pre
                # )
        cmfe_solver_eqs.BoundaryConditionsCreateFinish()

        # += solver for current iterations
        print("+===============================================================+")
        print(f'+= ... begin increment {i}')
        print("+===============================================================+")
        cmfe_problem.Solve()
        cmfe_problem.Finalise()
        cmfe_solver_eqs.Finalise()
    print('+==+ SOLVER COMPLETE')

    # +==+ ^\_/^ +==+ ^\_/^ +==+ 
    # Deformed Field:
    #    
    # +==+ ^\_/^ +==+ ^\_/^ +==+ 

    cmfe_def_field = iron.Field()
    cmfe_def_field.CreateStart(def_field_n, cmfe_region)
    cmfe_def_field.MeshDecompositionSet(cmfe_decomp)
    cmfe_def_field.TypeSet(iron.FieldTypes.GEOMETRIC)
    cmfe_def_field.VariableLabelSet(iron.FieldVariableTypes.U, "DeformedGeometry")
    cmfe_def_field.ComponentMeshComponentSet(iron.FieldVariableTypes.U, X, 1)
    cmfe_def_field.ComponentMeshComponentSet(iron.FieldVariableTypes.U, Y, 1)
    cmfe_def_field.ComponentMeshComponentSet(iron.FieldVariableTypes.U, Z, 1)
    cmfe_def_field.CreateFinish()
    print('+==+ DEFORMED FIELD COMPLETE')

    # +==+ ^\_/^ +==+ ^\_/^ +==+ 
    # Pressure Field:
    #    
    # +==+ ^\_/^ +==+ ^\_/^ +==+ 

    cmfe_pre_field = iron.Field()
    cmfe_pre_field.CreateStart(pre_field_n, cmfe_region)
    cmfe_pre_field.MeshDecompositionSet(cmfe_decomp)
    cmfe_pre_field.VariableLabelSet(iron.FieldVariableTypes.U, "Pressure")
    cmfe_pre_field.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 1, 1)
    cmfe_pre_field.ComponentInterpolationSet(
        iron.FieldVariableTypes.U, 
        1, 
        iron.FieldInterpolationTypes.ELEMENT_BASED
    )
    cmfe_pre_field.NumberOfComponentsSet(iron.FieldVariableTypes.U, 1)
    cmfe_pre_field.CreateFinish()
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
    vtk_output(cmfe_mesh, n_n, cmfe_geo_field, cmfe_dep_field, e_np_map, cmfe_mesh_e, RUNTIME_PATH, test_name)
    print('+==+ EXPORT COMPLETE')

    # +============+ 
    # Wrap it up
    # +============+ 

    # cmfe_problem.Destroy()
    cmfe_coord.Destroy()
    cmfe_region.Destroy()
    cmfe_basis.Destroy()
    iron.Finalise()

# +==+ ^\_/^ +==+ ^\_/^ +==+ 
# Run check
# +==+ ^\_/^ +==+ ^\_/^ +==+ 

if __name__ == '__main__':
    test_name = "hexa_test"
    main(test_name)