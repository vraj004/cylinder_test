#!/usr/bin/env python

#> \file
#> \author Vijay Rajagopal, Liam Murray
#> \brief This is an example script to solve a finite elasticity equation using OpenCMISS calls in python.
#>
#> \section LICENSE
#>
#> Version: MPL 1.1/GPL 2.0/LGPL 2.1
#>
#> The contents of this file are subject to the Mozilla Public License
#> Version 1.1 (the "License"); you may not use this file except in
#> compliance with the License. You may obtain a copy of the License at
#> http://www.mozilla.org/MPL/
#>
#> Software distributed under the License is distributed on an "AS IS"
#> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
#> License for the specific language governing rights and limitations
#> under the License.
#>
#> The Original Code is OpenCMISS
#>
#> The Initial Developer of the Original Code is University of Auckland,
#> Auckland, New Zealand and University of Oxford, Oxford, United
#> Kingdom. Portions created by the University of Auckland and University
#> of Oxford are Copyright (C) 2007 by the University of Auckland and
#> the University of Oxford. All Rights Reserved.
#>
#> Contributor(s): 
#>
#> Alternatively, the contents of this file may be used under the terms of
#> either the GNU General Public License Version 2 or later (the "GPL"), or
#> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
#> in which case the provisions of the GPL or the LGPL are applicable instead
#> of those above. if you wish to allow use of your version of this file only
#> under the terms of either the GPL or the LGPL, and not to allow others to
#> use your version of this file under the terms of the MPL, indicate your
#> decision by deleting the provisions above and replace them with the notice
#> and other provisions required by the GPL or the LGPL. if you do not delete
#> the provisions above, a recipient may use your version of this file under
#> the terms of any one of the MPL, the GPL or the LGPL.
#>

#> \example FiniteElasticity/UniAxialExtension/src/UniAxialExtensionExample.py
## Example script to solve a finite elasticity equation using OpenCMISS calls in python.
## \par Latest Builds:
## \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/UniAxialExtension/build-intel'>Linux Intel Build</a>
## \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/UniAxialExtension/build-gnu'>Linux GNU Build</a>
#<

#> Main script
# Add Python bindings directory to PATH

#sys.path.append(os.sep.join((os.environ['OPENCMISS_ROOT'],'cm','bindings','python')))

# Intialise OpenCMISS
from opencmiss.iron import iron
import math
import numpy as np
import meshio
import generateMesh


#Setting up the vtk output function before use after solve.

def vtkoutput(totalNumberOfNodes, totalNumberOfElements, mesh,geo_field,dep_field):
    X = 1
    Y = 2
    Z = 3
    meshNodes = iron.MeshNodes()
    mesh.NodesGet(1,meshNodes)
    n_list = []
    n_bef = []
    n_aft = []
    for i in range(0, totalNumberOfNodes, 1):
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

    e_list = []
    e_list.append(
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 25, 26, 27, 28, 29, 30, 31, 32, 33, 49, 50, 51, 52, 53, 54, 55, 56, 57])
    e_list.append(
        [7, 8, 9, 10, 11, 12, 13, 14, 15, 31, 32, 33, 34, 35, 36, 37, 38, 39, 55, 56, 57, 58, 59, 60, 61, 62, 63])
    e_list.append(
        [13, 14, 15, 16, 17, 18, 19, 20, 21, 37, 38, 39, 40, 41, 42, 43, 44, 45, 61, 62, 63, 64, 65, 66, 67, 68, 69])
    e_list.append(
        [19, 20, 21, 22, 23, 24, 1, 2, 3, 43, 44, 45, 46, 47, 48, 25, 26, 27, 67, 68, 69, 70, 71, 72, 49, 50, 51])

    e_list_iron = np.array(e_list)[:,:]-1


    # Iron Numbering for 27?
    #   *  z = 0           z = 0.5         z = 1          
    #   *  6--7 --8     * 15--16--17    x 24--25--26
    #   *  |      |     *  |      |     x  |      |
    #   *  3  4   5     * 12  13  14    x 21  22  23
    #   *  |      |     *  |      |     x  |      |
    #   *  0--1 --2     *  9--10--11    x 18--19--20

    # Vijay Numbering for Hexa-27
    #   *  z = 0           z = 0.5         z = 1    
    #   *  0-- 9--18    *  1--10--19    *  2--11--20     
    #   *  |      |     *  |      |     *  |      |       
    #   *  3  12  21    *  4  13  22    *  5  14  23     
    #   *  |      |     *  |      |     *  |      |     
    #   *  6--15--24    *  7--16--25    *  8--17--26      

    # VTK Numbering for Hexa-27
    #   *  z = 0           z = 0.5         z = 1    
    #   *  3--10--2     * 19--23--18    x  7--14--6
    #   *  |      |     *  |      |     x  |      |
    #   * 11  24  9     * 20  26  21    x 15  25  13
    #   *  |      |     *  |      |     x  |      |
    #   *  0-- 8--1     * 16--22--17    x  4--12--5 

    # IRON_VTK = [
    #     0, 2, 8, 6, 18, 20, 26, 24, 
    #     1, 5, 7, 3, 19, 13, 25, 21, 
    #     9, 11, 17, 15, 12, 14, 10, 16, 
    #     4, 22, 13
    # ]

    VIJAY_VTK = [
        6, 24, 18, 0, 8, 26, 20, 2, 
        15, 21, 9, 3, 17, 23, 11, 5, 
        7, 25, 19, 1, 4, 22, 16, 10, 
        12, 14, 13
    ]

    e_list_vtk = e_list_iron[:, VIJAY_VTK]

    # +============+
    # Store data files
    # +============+
    node_file = open('output_mesh.node', 'w')
    node_file.writelines([str(line) + "\n" for line in n_list])
    node_file.close()
    elem_file=open('output_mesh.ele','w')
    elem_file.writelines([str(line) + "\n" for line in e_list])
    elem_file.close()
    bef_def = np.array(n_bef)
    aft_def = np.array(n_aft)
    np.save('output_mesh_before_def.npy',bef_def)
    np.save('output_mesh_after_def.npy',aft_def)
    # +============+
    # VTK export
    # +============+

    meshio.write_points_cells(
        "output.vtk",
        bef_def,
        [("hexahedron27", e_list_vtk)] + [("hexahedron27", e_list_vtk)],
        {"deformed": aft_def-bef_def}
    )

    return

# Set problem parameters
#cylinder geometry inputs
r_inner = 1.0
r_outer = 1.5
z_length = 1.0
#material parameters
c10 = 2.0
c01 = 6.0
#loading conditions
pressure_internal = 5.0
pressure_external = -1.0
stretch_ratio = 0.8
twist_angle = 0.0 #math.pi/6.0  #in radians
#mesh parameters
#set number of elements in each direction of cylinder.
numberOfRadialElements = 1
numberOfCircumferentialElements = 4
numberOfZElements = 2
#set geometric basis interpolation to one of:
# 1 - linear lagrange
# 2 - quadratic lagrange
# 3 - cubic lagrange
# 4 - hermite cubic
# 5 - hermite quintic
# 6 - hermite septum
# 7 - linear simplex
# 8 - quadratic simplex
# 9 - cubic simplex

InterpolationType = 2
#decide if you want to use a pressure basis for the pressure field
UsePressureBasis = False
#From here the code sets up the problem with the inputs from above.

# Set OpenCMISS parameters
coordinateSystemUserNumber = 1
regionUserNumber = 1
geometricbasisUserNumber = 1

pressureBasisUserNumber = 2
generatedMeshUserNumber = 1
meshUserNumber = 1
decompositionUserNumber = 1
geometricFieldUserNumber = 1
fibreFieldUserNumber = 2
materialFieldUserNumber = 3
dependentFieldUserNumber = 4
equationsSetUserNumber = 1
equationsSetFieldUserNumber = 5
problemUserNumber = 1

# Get the number of computational nodes and this computational node number
numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()
computationalNodeNumber = iron.ComputationalNodeNumberGet()

# Create a 3D rectangular cartesian coordinate system
coordinateSystem = iron.CoordinateSystem()
coordinateSystem.CreateStart(coordinateSystemUserNumber)
coordinateSystem.DimensionSet(3)
coordinateSystem.CreateFinish()

# Create a region and assign the coordinate system to the region
region = iron.Region()
region.CreateStart(regionUserNumber,iron.WorldRegion)
region.LabelSet("Region")
region.coordinateSystem = coordinateSystem
region.CreateFinish()

# Define bases

#setting up basis function related settings and accordingly number of nodes and number of elements.
if(InterpolationType==2):
    totalNumberOfElements = numberOfRadialElements*numberOfCircumferentialElements*numberOfZElements
    if InterpolationType==1:
        totalNumberOfNodes = 8*totalNumberOfElements-(4*totalNumberOfElements)
    elif InterpolationType==2:
        totalNumberOfNodes = 27*totalNumberOfElements-(9*totalNumberOfElements)
if(UsePressureBasis):
    numberOfMeshComponents = 2
    if InterpolationType==1:
        totalNumberOfPressureNodes = 8*totalNumberOfElements-(4*totalNumberOfElements)
    elif (InterpolationType==2):
        print("Using quadratic lagrange pressure basis is not supported, please change to linear interpolation or no pressure basis (only element varying)")
        exit()
else:
    numberOfMeshComponents = 1

#hard coding two parameters because the cylinder is 3D and lower number of gauss points per xi direction than 4 leads to significant model simulation errors.
numberOfXi = 3
NumberOfGaussXi = 4

geometricBasis = iron.Basis()
geometricBasis.CreateStart(geometricbasisUserNumber)
if InterpolationType in (1,2,3,4):
    geometricBasis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
elif InterpolationType in (7,8,9):
    geometricBasis.type = iron.BasisTypes.SIMPLEX
geometricBasis.numberOfXi = numberOfXi
if (InterpolationType==1):
    geometricBasis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*numberOfXi
elif (InterpolationType==2):
    geometricBasis.interpolationXi = [iron.BasisInterpolationSpecifications.QUADRATIC_LAGRANGE]*numberOfXi

geometricBasis.quadratureNumberOfGaussXi = [NumberOfGaussXi]*numberOfXi
geometricBasis.CreateFinish()

if(UsePressureBasis):
    # Define pressure basis
    pressureBasis = iron.Basis()
    pressureBasis.CreateStart(pressureBasisUserNumber)
    if InterpolationType in (1,2,3,4):
        pressureBasis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
    elif InterpolationType in (7,8,9):
        pressureBasis.type = iron.BasisTypes.SIMPLEX
    pressureBasis.numberOfXi = numberOfXi
    if InterpolationType==1:
        pressureBasis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*numberOfXi
    elif InterpolationType==2: #this should be unnecessary because i have caught this case above
        print("Using quadratic lagrange pressure basis is not supported, please change to linear interpolation or no pressure basis (only element varying)")
        exit()
    if(NumberOfGaussXi>0):
        pressureBasis.quadratureNumberOfGaussXi = [NumberOfGaussXi]*numberOfXi
    pressureBasis.CreateFinish()

# Start the creation of a manually generated mesh in the region
mesh = iron.Mesh()
mesh.CreateStart(meshUserNumber,region,numberOfXi)
mesh.NumberOfComponentsSet(numberOfMeshComponents)
mesh.NumberOfElementsSet(totalNumberOfElements)

#Define nodes for the mesh
nodes = iron.Nodes()
nodes.CreateStart(region,totalNumberOfNodes)
nodes.CreateFinish()
node_list, node_idx_list, top_node_list,bottom_node_list, yfix_node_list, xfix_node_list, internal_node_list,outer_node_list = generateMesh.annulus(r_inner, r_outer, z_length, numberOfRadialElements, numberOfCircumferentialElements, numberOfZElements, InterpolationType)

# Define elements for the mesh
elements = iron.MeshElements()
meshComponentNumber = 1
elements.CreateStart(mesh,meshComponentNumber,geometricBasis)


#els = [
#    list(range(1, 1+9, 1)) + list(range(25, 25+9, 1)) + list(range(49, 49+9, 1)),
#    list(range(7, 7+9, 1)) + list(range(31, 31+9, 1)) + list(range(55, 55+9, 1)),
#    list(range(13, 13+9, 1)) + list(range(37, 37+9, 1)) + list(range(61, 61+9, 1)),
#    [19, 20, 21, 22, 23, 24, 1, 2, 3] + [43, 44, 45, 46, 47, 48, 25, 26, 27] + [67, 68, 69, 70, 71, 72, 49, 50, 51]
#]

elements.CreateFinish()
mesh.CreateFinish() 

# Create a decomposition for the mesh
decomposition = iron.Decomposition()
decomposition.CreateStart(decompositionUserNumber,mesh)
decomposition.type = iron.DecompositionTypes.CALCULATED
decomposition.numberOfDomains = numberOfComputationalNodes
decomposition.CreateFinish()

# Create a field for the geometry
geometricField = iron.Field()
geometricField.CreateStart(geometricFieldUserNumber,region)
geometricField.MeshDecompositionSet(decomposition)
geometricField.TypeSet(iron.FieldTypes.GEOMETRIC)
geometricField.VariableLabelSet(iron.FieldVariableTypes.U,"Geometry")
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,1)
if InterpolationType == 4:
    geometricField.fieldScalingType = iron.FieldScalingTypes.ARITHMETIC_MEAN
geometricField.CreateFinish()

# Update the geometric field parameters manually for a cube with a cube-shaped hole in it.
#outer cube unit length = 2.0 mm
#inner cubic hole unit length = 1.0 mm
#Number of nodes = 16
geometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
for nodeidx in range(0,totalNumberOfNodes,1):
    nodeNum = nodeidx+1
    nodex = node_list[nodeidx][0]
    nodey = node_list[nodeidx][1]
    nodez = node_list[nodeidx][2]
    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,nodeNum,1,nodex)
    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,nodeNum,2,nodey)
    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,nodeNum,3,nodez)

geometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

# Create a fibre field and attach it to the geometric field
fibreField = iron.Field()
fibreField.CreateStart(fibreFieldUserNumber,region)
fibreField.TypeSet(iron.FieldTypes.FIBRE)
fibreField.MeshDecompositionSet(decomposition)
fibreField.GeometricFieldSet(geometricField)
fibreField.VariableLabelSet(iron.FieldVariableTypes.U,"Fibre")
if InterpolationType == 4:
    fibreField.fieldScalingType = iron.FieldScalingTypes.ARITHMETIC_MEAN
fibreField.CreateFinish()

# Create the material field
materialField = iron.Field()
materialField.CreateStart(materialFieldUserNumber,region)
materialField.TypeSet(iron.FieldTypes.MATERIAL)
materialField.MeshDecompositionSet(decomposition)
materialField.GeometricFieldSet(geometricField)
materialField.NumberOfVariablesSet(1)
materialField.NumberOfComponentsSet(iron.FieldVariableTypes.U,2)
materialField.VariableLabelSet(iron.FieldVariableTypes.U,"Material")
materialField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
materialField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
if InterpolationType == 4:
    materialField.fieldScalingType = iron.FieldScalingTypes.ARITHMETIC_MEAN
materialField.CreateFinish()

# Set Mooney-Rivlin constants c10 and c01 respectively.
iron.Field.ComponentValuesInitialiseDP(
    materialField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,c01)
iron.Field.ComponentValuesInitialiseDP(
    materialField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,c10)

# Create the dependent field
dependentField = iron.Field()
dependentField.CreateStart(dependentFieldUserNumber,region)
dependentField.VariableLabelSet(iron.FieldVariableTypes.U,"Dependent")
dependentField.TypeSet(iron.FieldTypes.GEOMETRIC_GENERAL)  
dependentField.MeshDecompositionSet(decomposition)
dependentField.GeometricFieldSet(geometricField) 
dependentField.DependentTypeSet(iron.FieldDependentTypes.DEPENDENT) 
dependentField.NumberOfVariablesSet(2)
dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.U,4)
dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.DELUDELN,4)
dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,1)  
dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,1,1)
dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,2,1)
dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,3,1)  
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U,4,iron.FieldInterpolationTypes.ELEMENT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.DELUDELN,4,iron.FieldInterpolationTypes.ELEMENT_BASED)
if(UsePressureBasis):
    # Set the pressure to be nodally based and use the second mesh component
    if InterpolationType == 4:
        dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U,4,iron.FieldInterpolationTypes.NODE_BASED)
        dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.DELUDELN,4,iron.FieldInterpolationTypes.NODE_BASED)
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,4,2)
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,4,2)
if InterpolationType == 4:
    dependentField.fieldScalingType = iron.FieldScalingTypes.ARITHMETIC_MEAN
dependentField.CreateFinish()

# Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
iron.Field.ParametersToFieldParametersComponentCopy(
    geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,
    dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1)
iron.Field.ParametersToFieldParametersComponentCopy(
    geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,
    dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2)
iron.Field.ParametersToFieldParametersComponentCopy(
    geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3,
    dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3)
iron.Field.ComponentValuesInitialiseDP(
    dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,4,1.2345678806304932)

# Create the equations_set
equationsSetField = iron.Field()
equationsSet = iron.EquationsSet()
equationsSetSpecification = [iron.EquationsSetClasses.ELASTICITY,
    iron.EquationsSetTypes.FINITE_ELASTICITY,
    iron.EquationsSetSubtypes.MOONEY_RIVLIN]
equationsSet.CreateStart(equationsSetUserNumber,region,fibreField,
    equationsSetSpecification, equationsSetFieldUserNumber, equationsSetField)
equationsSet.CreateFinish()

equationsSet.MaterialsCreateStart(materialFieldUserNumber,materialField)
equationsSet.MaterialsCreateFinish()

equationsSet.DependentCreateStart(dependentFieldUserNumber,dependentField)
equationsSet.DependentCreateFinish()

# Create equations
equations = iron.Equations()
equationsSet.EquationsCreateStart(equations)
equations.sparsityType = iron.EquationsSparsityTypes.SPARSE
equations.outputType = iron.EquationsOutputTypes.NONE
equationsSet.EquationsCreateFinish()

# Define the problem
problem = iron.Problem()
problemSpecification = [iron.ProblemClasses.ELASTICITY,
        iron.ProblemTypes.FINITE_ELASTICITY,
        iron.ProblemSubtypes.NONE]
problem.CreateStart(problemUserNumber, problemSpecification)
problem.CreateFinish()

# Create control loops
problem.ControlLoopCreateStart()
problem.ControlLoopCreateFinish()

# Create problem solver
nonLinearSolver = iron.Solver()
linearSolver = iron.Solver()
problem.SolversCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,nonLinearSolver)
nonLinearSolver.outputType = iron.SolverOutputTypes.PROGRESS
nonLinearSolver.NewtonJacobianCalculationTypeSet(iron.JacobianCalculationTypes.FD)
nonLinearSolver.NewtonLinearSolverGet(linearSolver)
linearSolver.linearType = iron.LinearSolverTypes.DIRECT
#linearSolver.libraryType = iron.SolverLibraries.LAPACK
problem.SolversCreateFinish()

# Create solver equations and add equations set to solver equations
solver = iron.Solver()
solverEquations = iron.SolverEquations()
problem.SolverEquationsCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,solver)
solver.SolverEquationsGet(solverEquations)
solverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
equationsSetIndex = solverEquations.EquationsSetAdd(equationsSet)
problem.SolverEquationsCreateFinish()

# Prescribe boundary conditions (absolute nodal parameters)
boundaryConditions = iron.BoundaryConditions()
solverEquations.BoundaryConditionsCreateStart(boundaryConditions)

for nodeidx in range(0,len(bottom_node_list),1):
    nodeNum = bottom_node_list[nodeidx][0]
    print(nodeNum)
    boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,nodeNum,3,iron.BoundaryConditionsTypes.FIXED_INCREMENTED,0.0)
for nodeidx in range(0,len(yfix_node_list),1):
    nodeNum = yfix_node_list[nodeidx][0]
    boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,nodeNum,2,iron.BoundaryConditionsTypes.FIXED_INCREMENTED,0.0)
for nodeidx in range(0,len(xfix_node_list),1):
    nodeNum = xfix_node_list[nodeidx][0]
    boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,nodeNum,1,iron.BoundaryConditionsTypes.FIXED_INCREMENTED,0.0)

#Pressure boundary conditions on the internal faces.
for nodeidx in range(0,len(internal_node_list),1):
    nodeNum = internal_node_list[nodeidx][0]
    boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.DELUDELN,1,1,nodeNum,3,iron.BoundaryConditionsTypes.PRESSURE_INCREMENTED,pressure_internal)

#Pressure boundary conditions on the external faces.
for nodeidx in range(0,len(outer_node_list),1):
    nodeNum = outer_node_list[nodeidx][0]
    boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.DELUDELN,1,1,nodeNum,3,iron.BoundaryConditionsTypes.PRESSURE_INCREMENTED,pressure_external)


#Apply stretch boundary condition on top faces.
for nodeidx in range(0,len(top_node_list),1):
    nodeNum = top_node_list[nodeidx][0]
    stretch_disp = stretch_ratio*z_length-z_length
    boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,nodeNum,3,iron.BoundaryConditionsTypes.FIXED_INCREMENTED,stretch_disp)

solverEquations.BoundaryConditionsCreateFinish()

# Solve the problem
problem.Solve()

vtkoutput(totalNumberOfNodes,totalNumberOfElements,mesh,geometricField,dependentField)
# Export results
fields = iron.Fields()
fields.CreateRegion(region)
fields.NodesExport("CylinderInflation","FORTRAN")
fields.ElementsExport("CylinderInflation","FORTRAN")
fields.Finalise()

