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
import os

DIM = 3
N_N_EL = 27
QUAD_ORDER = 4
X, Y, Z = (1, 2, 3)
PRESSURE_TEST = True
LOADSTEPS = 10
INNER_RAD = 0.375

(
    coor_n, basi_n, regi_n, mesh_n, deco_n,
    geo_field_n, dep_field_n, mat_field_n,
    eqs_field_n, eqs_n, prob_n, def_field_n,
    pre_field_n
) = range(1, 14)

def coordinate_setup():
    coord_sys = iron.CoordinateSystem()
    coord_sys.CreateStart(coor_n)
    coord_sys.DimensionSet(DIM)
    coord_sys.CreateFinish()

    return coord_sys

def node_setup(test_name):
    n_idx = []
    n_xyz = []
    with open("/home/jovyan/work/test_files/runtime_files/" + test_name + ".nodes", 'r') as n_file:
        for i, line in enumerate(n_file):
            if i == 0: continue
            line = line.strip().split('\t')
            n_idx.append(int(line[0]))
            n_xyz.append(line[1:])

    return np.array(n_xyz).astype(float), n_idx, i

def element_setup(test_name):
    e_idx = []
    e_map = []
    with open("/home/jovyan/work/test_files/runtime_files/" + test_name + ".ele", 'r') as e_file:
        for i, line in enumerate(e_file):
            if i == 0: continue
            line = line.strip().split('\t')
            e_idx.append(i)
            e_map.append(line[3:])

    return e_map, e_idx, i

def basis_setup():
    xi_n = 3
    basis = iron.Basis()
    basis.CreateStart(basi_n)
    basis.TypeSet(iron.BasisTypes.LAGRANGE_HERMITE_TP)
    basis.NumberOfXiSet(xi_n)
    basis.InterpolationXiSet(
        [iron.BasisInterpolationSpecifications.QUADRATIC_LAGRANGE] * xi_n
    )
    basis.QuadratureNumberOfGaussXiSet([3]*xi_n)
    basis.CreateFinish()

    return basis

def region_setup(coord_sys):
    region = iron.Region()
    region.CreateStart(regi_n, iron.WorldRegion)
    region.CoordinateSystemSet(coord_sys)
    region.LabelSet("Region")
    region.CreateFinish()

    return region

def mesh_setup(region, n_n, e_n, basis, e_idx, e_np_map):
    nodesList = iron.Nodes()
    nodesList.CreateStart(region, n_n)
    nodesList.CreateFinish()
    #
    mesh = iron.Mesh()
    mesh.CreateStart(mesh_n, region, DIM)
    mesh.NumberOfElementsSet(e_n)
    mesh.NumberOfComponentsSet(1)
    #
    mesh_e = iron.MeshElements()
    mesh_e.CreateStart(mesh, 1, basis)

    for i in range(e_n):
        
        nodesList = list(
            map(int,e_np_map[i][:])
        )
        mesh_e.NodesSet(e_idx[i], nodesList)

    mesh_e.CreateFinish()
    mesh.CreateFinish()

    return mesh_e, mesh

def decomposition_setup(mesh):
    comp_nodes_n = iron.ComputationalNumberOfNodesGet()
    decomp = iron.Decomposition()
    decomp.CreateStart(deco_n, mesh)
    decomp.TypeSet(iron.DecompositionTypes.CALCULATED)
    decomp.NumberOfDomainsSet(comp_nodes_n)
    decomp.CreateFinish()

    return decomp

def geometric_setup(region, decomp, n_n, n_idx, n_xyz):
    geo_field = iron.Field()
    geo_field.CreateStart(geo_field_n, region) 
    geo_field.MeshDecompositionSet(decomp)
    geo_field.VariableLabelSet(iron.FieldVariableTypes.U, "Geometry")
    geo_field.ScalingTypeSet(iron.FieldScalingTypes.ARITHMETIC_MEAN)
    geo_field.CreateFinish()
    #
    geo_field.ParameterSetUpdateStart(
        iron.FieldVariableTypes.U, 
        iron.FieldParameterSetTypes.VALUES
    )

    # computationalNodeNumber = iron.ComputationalNodeNumberGet()
    for idx in range(n_n):
        n_id = n_idx[idx]
        n_x, n_y, n_z = (n_xyz[idx][0], n_xyz[idx][1], n_xyz[idx][2])

        geo_field.ParameterSetUpdateNodeDP(
            iron.FieldVariableTypes.U, 
            iron.FieldParameterSetTypes.VALUES,
            1, 
            1, 
            n_id, 
            1, 
            n_x
        )
        geo_field.ParameterSetUpdateNodeDP(
            iron.FieldVariableTypes.U, 
            iron.FieldParameterSetTypes.VALUES,
            1, 
            1, 
            n_id, 
            2, 
            n_y
        )
        geo_field.ParameterSetUpdateNodeDP(
            iron.FieldVariableTypes.U, 
            iron.FieldParameterSetTypes.VALUES,
            1, 
            1, 
            n_id, 
            3, 
            n_z
        )
    
    geo_field.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
    
    return geo_field

def material_setup(decomp, geo_field, region):
    mat_field = iron.Field()
    mat_field.CreateStart(mat_field_n, region)
    mat_field.TypeSet(iron.FieldTypes.MATERIAL)
    mat_field.MeshDecompositionSet(decomp)
    mat_field.GeometricFieldSet(geo_field)
    mat_field.VariableLabelSet(iron.FieldVariableTypes.U, "Material")
    mat_field.NumberOfComponentsSet(iron.FieldVariableTypes.U, 2)
    mat_field.ScalingTypeSet(iron.FieldScalingTypes.ARITHMETIC_MEAN)
    # Loop through components
    for component in [1, 2]:
        mat_field.ComponentInterpolationSet(
        iron.FieldVariableTypes.U, component,
        iron.FieldInterpolationTypes.ELEMENT_BASED)
    mat_field.CreateFinish()
    mat_field.ComponentValuesInitialiseDP(
        iron.FieldVariableTypes.U, 
        iron.FieldParameterSetTypes.VALUES, 
        1, 
        2.0
    )
    mat_field.ComponentValuesInitialiseDP(
        iron.FieldVariableTypes.U, 
        iron.FieldParameterSetTypes.VALUES, 
        2, 
        6.0
    )

    return mat_field

def equations_setup(region, geo_field, mat_field):
    eqs_set_field = iron.Field()
    eqs_set = iron.EquationsSet()
    #
    eqs_set_specs = [
        iron.ProblemClasses.ELASTICITY,
        iron.ProblemTypes.FINITE_ELASTICITY,
        iron.EquationsSetSubtypes.MOONEY_RIVLIN
    ]
    #
    eqs_set.CreateStart(
        eqs_n, 
        region, 
        geo_field, 
        eqs_set_specs,
        eqs_field_n, 
        eqs_set_field
    )

    eqs_set.MaterialsCreateStart(mat_field_n, mat_field)
    eqs_set.MaterialsCreateFinish()
    eqs_set.CreateFinish()

    return eqs_set_field, eqs_set

def dependent_setup(region, decomp, geo_field, eqs_set):
    dep_field = iron.Field()
    eqs_set.DependentCreateStart(dep_field_n, dep_field)
    dep_field.VariableLabelSet(iron.FieldVariableTypes.U, "Dependent")

    for i in [1, 2, 3]:
        dep_field.ComponentMeshComponentSet(iron.FieldVariableTypes.U, i, 1)
        dep_field.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN, i, 1)
    dep_field.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 4, 1)
    dep_field.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN, 4, 1)

    dep_field.ComponentInterpolationSet(iron.FieldVariableTypes.U, 4,
                                                iron.FieldInterpolationTypes.NODE_BASED)
    dep_field.ComponentInterpolationSet(iron.FieldVariableTypes.DELUDELN, 4,
                                                iron.FieldInterpolationTypes.NODE_BASED)
    dep_field.ScalingTypeSet(iron.FieldScalingTypes.ARITHMETIC_MEAN)

    eqs_set.DependentCreateFinish()

    for i in [1, 2, 3]:
        iron.Field.ParametersToFieldParametersComponentCopy(geo_field, iron.FieldVariableTypes.U,
                                                             iron.FieldParameterSetTypes.VALUES, i, dep_field,
                                                             iron.FieldVariableTypes.U,
                                                             iron.FieldParameterSetTypes.VALUES, i)

    # Set hydrostatic pressure initial guess.
    dep_field.ComponentValuesInitialise(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 4,
                                             0.0)

    dep_field.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
    dep_field.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
    
    # if PRESSURE_TEST:
    #     dep_field.ComponentInterpolationSet(
    #         iron.FieldVariableTypes.U, 4, iron.FieldInterpolationTypes.GAUSS_POINT_BASED) # CHECK
    #     dep_field.ComponentMeshComponentSet(
    #         iron.FieldVariableTypes.U, 4, 1)
    #     # DELUDELN variable (forces)
    #     dep_field.ComponentInterpolationSet(
    #         iron.FieldVariableTypes.DELUDELN, 4,
    #         iron.FieldInterpolationTypes.NODE_BASED)
    #     dep_field.ComponentMeshComponentSet(
    #         iron.FieldVariableTypes.DELUDELN, 4, 1)
    # else:
    #     # Set the hydrostatic pressure to be constant within each element.
    #     dep_field.ComponentInterpolationSet(
    #         iron.FieldVariableTypes.U, 4,
    #         iron.FieldInterpolationTypes.ELEMENT_BASED)
    #     dep_field.ComponentInterpolationSet(
    #         iron.FieldVariableTypes.DELUDELN, 4,
    #         iron.FieldInterpolationTypes.ELEMENT_BASED)
    # # dep_field.CreateFinish()
    # eqs_set.DependentCreateFinish()
    # #
    # iron.Field.ParametersToFieldParametersComponentCopy(
    #     geo_field,
    #     iron.FieldVariableTypes.U, 
    #     iron.FieldParameterSetTypes.VALUES, 
    #     1,
    #     dep_field, 
    #     iron.FieldVariableTypes.U, 
    #     iron.FieldParameterSetTypes.VALUES, 
    #     1
    # ) 
    # iron.Field.ParametersToFieldParametersComponentCopy(
    #     geo_field, 
    #     iron.FieldVariableTypes.U, 
    #     iron.FieldParameterSetTypes.VALUES, 
    #     2,
    #     dep_field, 
    #     iron.FieldVariableTypes.U, 
    #     iron.FieldParameterSetTypes.VALUES, 
    #     2
    # ) 
    # iron.Field.ParametersToFieldParametersComponentCopy(
    #     geo_field, 
    #     iron.FieldVariableTypes.U, 
    #     iron.FieldParameterSetTypes.VALUES, 
    #     3, 
    #     dep_field, 
    #     iron.FieldVariableTypes.U, 
    #     iron.FieldParameterSetTypes.VALUES, 
    #     3
    # ) 
    # iron.Field.ComponentValuesInitialiseDP(
    #     dep_field, 
    #     iron.FieldVariableTypes.U, 
    #     iron.FieldParameterSetTypes.VALUES, 
    #     4, 
    #     0.0
    # ) 

    eqs = iron.Equations()
    eqs_set.EquationsCreateStart(eqs)
    eqs.SparsityTypeSet(iron.EquationsSparsityTypes.SPARSE)
    eqs.OutputTypeSet(iron.EquationsOutputTypes.NONE)
    eqs_set.EquationsCreateFinish()

    return dep_field, geo_field, eqs, eqs_set

def problem_setup():
    problem = iron.Problem()
    problems_specs = (
        [
            iron.ProblemClasses.ELASTICITY,
            iron.ProblemTypes.FINITE_ELASTICITY,
            iron.ProblemSubtypes.NONE
        ]
    )
    problem.CreateStart(prob_n, problems_specs)
    problem.CreateFinish()
    #
    problem.ControlLoopCreateStart()
    ctrl_loop = iron.ControlLoop()
    problem.ControlLoopGet([iron.ControlLoopIdentifiers.NODE], ctrl_loop)
    ctrl_loop.MaximumIterationsSet(LOADSTEPS)
    problem.ControlLoopCreateFinish()

    return problem, ctrl_loop

def solver_setup(problem, eqs_set):
    nl_solver = iron.Solver()
    li_solver = iron.Solver()
    problem.SolversCreateStart()
    problem.SolverGet(
        [iron.ControlLoopIdentifiers.NODE], 
        1, 
        nl_solver
    )
    nl_solver.OutputTypeSet(iron.SolverOutputTypes.PROGRESS)
    nl_solver.NewtonJacobianCalculationTypeSet(
        iron.JacobianCalculationTypes.EQUATIONS
    )
    nl_solver.NewtonLinearSolverGet(li_solver)
    li_solver.LinearTypeSet(iron.LinearSolverTypes.DIRECT)
    problem.SolversCreateFinish()
    # 
    solver = iron.Solver()
    solver_eqs = iron.SolverEquations()
    problem.SolverEquationsCreateStart()
    problem.SolverGet([iron.ControlLoopIdentifiers.NODE], 1, solver)
    solver.SolverEquationsGet(solver_eqs)
    solver_eqs.SparsityTypeSet(iron.SolverEquationsSparsityTypes.SPARSE)
    _ = solver_eqs.EquationsSetAdd(eqs_set)
    problem.SolverEquationsCreateFinish()

    return problem, nl_solver, li_solver, solver, solver_eqs

def boundary_conditions_setup(solver_eqs, dep_field, n_n, n_np_xyz):
    bcs = iron.BoundaryConditions()
    solver_eqs.BoundaryConditionsCreateStart(bcs)
    
    for i in range(0, n_n, 1):

        if (n_np_xyz[i, 2] == np.min(n_np_xyz[:, 2])) or (n_np_xyz[i, 2] == np.max(n_np_xyz[:, 2])):
            # bcs.AddNode(dep_field, iron.FieldVariableTypes.U,1,1,i,X,iron.BoundaryConditionsTypes.FIXED,0.0)
            # bcs.AddNode(dep_field, iron.FieldVariableTypes.U,1,1,i,Y,iron.BoundaryConditionsTypes.FIXED,0.0)
            bcs.AddNode(
                dep_field, 
                iron.FieldVariableTypes.U,
                1,
                1,
                i,
                Z,
                iron.BoundaryConditionsTypes.FIXED,
                0.0
            )

        vecNorm = np.linalg.norm([n_np_xyz[i, 0], n_np_xyz[i, 1]])
        if abs(vecNorm - INNER_RAD) < 1e-5:
            bcs.AddNode(
                dep_field, 
                iron.FieldVariableTypes.U,
                1,
                1,
                i,
                X,
                iron.BoundaryConditionsTypes.PRESSURE,
                1.5
            )
            bcs.AddNode(
                dep_field, 
                iron.FieldVariableTypes.U,
                1,
                1,
                i,
                Y,
                iron.BoundaryConditionsTypes.PRESSURE,
                1.5
            )
            bcs.AddNode(
                dep_field, 
                iron.FieldVariableTypes.U,
                1,
                1,
                i,
                Z,
                iron.BoundaryConditionsTypes.PRESSURE,
                1.5
            )

    solver_eqs.BoundaryConditionsCreateFinish()

    return bcs, solver_eqs

def output(region, decomp, dep_field, mesh, geo_field, n_n, e_n, mesh_e):

    def_field = iron.Field()
    def_field.CreateStart(def_field_n, region)
    def_field.MeshDecompositionSet(decomp)
    def_field.TypeSet(iron.FieldTypes.GEOMETRIC)
    def_field.VariableLabelSet(iron.FieldVariableTypes.U, "DeformedGeometry")
    for component in [1, 2, 3]:
        def_field.ComponentMeshComponentSet(
                iron.FieldVariableTypes.U, component, 1)
    def_field.CreateFinish()

    pre_field = iron.Field()
    pre_field.CreateStart(pre_field_n, region)
    pre_field.MeshDecompositionSet(decomp)
    pre_field.VariableLabelSet(iron.FieldVariableTypes.U, "Pressure")
    pre_field.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 1, 1)
    pre_field.ComponentInterpolationSet(
        iron.FieldVariableTypes.U, 
        1, 
        iron.FieldInterpolationTypes.ELEMENT_BASED
    )
    pre_field.NumberOfComponentsSet(iron.FieldVariableTypes.U, 1)
    pre_field.CreateFinish()

    for component in [1, 2, 3]:
        dep_field.ParametersToFieldParametersComponentCopy(
            iron.FieldVariableTypes.U,
            iron.FieldParameterSetTypes.VALUES, component,
            def_field, iron.FieldVariableTypes.U,
            iron.FieldParameterSetTypes.VALUES, component)

    dep_field.ParametersToFieldParametersComponentCopy(
        iron.FieldVariableTypes.U,
        iron.FieldParameterSetTypes.VALUES,
        4,
        pre_field, 
        iron.FieldVariableTypes.U,
        iron.FieldParameterSetTypes.VALUES, 
        1
    )

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
                    1,
                    1,
                    i+1,
                    1
                ),
                geo_field.ParameterSetGetNodeDP(
                    iron.FieldVariableTypes.U,
                    iron.FieldParameterSetTypes.VALUES,
                    1,
                    1,
                    i+1,
                    2
                ),
                geo_field.ParameterSetGetNodeDP(
                    iron.FieldVariableTypes.U,
                    iron.FieldParameterSetTypes.VALUES,
                    1,
                    1,
                    i+1,
                    3
                )
            ]
        )
        n_aft.append(
            [
                dep_field.ParameterSetGetNodeDP(
                    iron.FieldVariableTypes.U,
                    iron.FieldParameterSetTypes.VALUES,
                    1,
                    1,
                    i+1,
                    1
                ),
                dep_field.ParameterSetGetNodeDP(
                    iron.FieldVariableTypes.U,
                    iron.FieldParameterSetTypes.VALUES,
                    1,
                    1,
                    i+1,
                    2
                ),
                dep_field.ParameterSetGetNodeDP(
                    iron.FieldVariableTypes.U,
                    iron.FieldParameterSetTypes.VALUES,
                    1,
                    1,
                    i+1,
                    3
                )
            ]
        )
        n_list.append(
            [
                i+1, 
                n_bef[i][0], 
                n_bef[i][1],
                n_bef[i][2], 
                n_aft[i][0], 
                n_aft[i][1], 
                n_aft[i][2]
            ]
        )

    n_file = open('/home/jovyan/work/test_files/runtime_files/output_mesh.node', 'w')
    n_file.writelines([str(line) + "\n" for line in n_list])
    n_file.close()

    e_list = []
    for j in range(0, e_n, 1):
        elem_n = mesh_e.NodesGet(1+j, 27)
        e_list[j].append(elem_n)

    elem_file=open('/home/jovyan/work/test_files/runtime_files/output_mesh.ele','w')
    elem_file.writelines([str(line) + "\n" for line in e_list])
    elem_file.close()

    before_def = np.array(n_bef)
    np.save('/home/jovyan/work/test_files/runtime_files/output_mesh_before_def.npy',before_def)
    after_def = np.array(n_aft)
    np.save('/home/jovyan/work/test_files/runtime_files/output_mesh_after_def.npy',after_def)

    points = before_def
    point_data = {"deformed": after_def}

    cells_array = np.array(e_list)[:,:] - 1
    cells = [("tetra10", cells_array)] + [("tetra10", cells_array)]

    meshio.write_points_cells("/home/jovyan/work/test_files/runtime_files/output_mesh.vtk",points, cells, point_data)

    meshio.write_points_cells("/home/jovyan/work/test_files/runtime_files/output_mesh_iron.vtk",points,cells, point_data)

    return

def main(test_name):

    coord_sys = coordinate_setup()

    n_np_xyz, n_idx, n_n = node_setup(test_name)
    e_np_map, e_idx, e_n = element_setup(test_name)

    basis = basis_setup()
    region = region_setup(coord_sys)
    mesh_e, mesh = mesh_setup(region, n_n, e_n, basis, e_idx, e_np_map)
    decomp = decomposition_setup(mesh)

    geo_field = geometric_setup(region, decomp, n_n, n_idx, n_np_xyz)
    mat_field = material_setup(decomp, geo_field, region)

    eqs_set_field, eqs_set = equations_setup(region, geo_field, mat_field)
    dep_field, geo_field, eqs_set, eqs = dependent_setup(region, decomp, geo_field, eqs_set)
    problem, ctrl_loop = problem_setup()
    problem, nl_solver, li_solver, solver, solver_eqs = solver_setup(problem, eqs_set)
    bcs, solver_eqs = boundary_conditions_setup(solver_eqs, dep_field, n_n, n_np_xyz)

    try:
        problem.Solve()
    except:
        print("Failed")

    output(region, decomp, dep_field, mesh, geo_field, n_n, e_n, mesh_e)

    # Wrap it up
    problem.Destroy()
    coord_sys.Destroy()
    region.Destroy()
    basis.Destroy()
    iron.Finalise()

if __name__ == '__main__':
    print(os.getcwd())
    test_name = "cylinder_hexa"
    main(test_name)