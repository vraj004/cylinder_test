#!/usr/bin/env python

#> \file
#> \author Chris Bradley
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
import sys, os

#sys.path.append(os.sep.join((os.environ['OPENCMISS_ROOT'],'cm','bindings','python')))

# Intialise OpenCMISS
from opencmiss.iron import iron
import math
import numpy as np
import meshio

# Set problem parameters
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
    e_list.append([1, 2, 3, 4, 5, 6, 7, 8])
    e_list.append([9, 10, 1, 2, 11, 12, 5, 6])
    e_list.append([13, 14, 9, 10, 15, 16, 11, 12])
    e_list.append([3, 4, 13, 14, 7, 8, 15, 16])

    e_list_iron = np.array(e_list)[:,:]-1
    e_list_vtk = e_list_iron[:, [0, 1, 3, 2,4, 5, 7, 6]]
    print(e_list_vtk)

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
        [("hexahedron", e_list_vtk)] + [("hexahedron", e_list_vtk)],
        {"deformed": aft_def-bef_def}
    )

    return
UsePressureBasis = False
NumberOfGaussXi = 2

coordinateSystemUserNumber = 1
regionUserNumber = 1
basisUserNumber = 1
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
numberGlobalZElements=1
# Set all diganostic levels on for testing
iron.DiagnosticsSetOn(iron.DiagnosticTypes.ALL,[1,2,3,4,5],"Diagnostics",["DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE"])

totalNumberOfNodes=72
totalNumberOfElements=4
InterpolationType = 1
if(UsePressureBasis):
  numberOfMeshComponents = 2
else:
  numberOfMeshComponents = 1
if(numberGlobalZElements==0):
    numberOfXi = 2
else:
    numberOfXi = 3

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

# Define basis
basis = iron.Basis()
basis.CreateStart(basisUserNumber)
if InterpolationType in (1,2,3,4):
    basis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
elif InterpolationType in (7,8,9):
    basis.type = iron.BasisTypes.SIMPLEX
basis.numberOfXi = numberOfXi
basis.interpolationXi = [iron.BasisInterpolationSpecifications.QUADRATIC_LAGRANGE]*numberOfXi
if(NumberOfGaussXi>0):
    basis.quadratureNumberOfGaussXi = [NumberOfGaussXi]*numberOfXi
basis.CreateFinish()

if(UsePressureBasis):
    # Define pressure basis
    pressureBasis = iron.Basis()
    pressureBasis.CreateStart(pressureBasisUserNumber)
    if InterpolationType in (1,2,3,4):
        pressureBasis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
    elif InterpolationType in (7,8,9):
        pressureBasis.type = iron.BasisTypes.SIMPLEX
    pressureBasis.numberOfXi = numberOfXi
    pressureBasis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*numberOfXi
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
#create list of nodes and their coordinates in the order you are going to define the elements. This will make element definition easier
node_list = []
elem_theta_delta = math.pi/2.0
theta_delta = math.pi/4.0
xorig = 0.0
yorig = 0.0
zorig = 0.0
theta_orig=2*math.pi
r_inner = 1.0
r_outer = 1.5
r_delta = (r_outer-r_inner)/2.0
r_mid = r_inner+r_delta
z_delta = 0.5
for ridx in range(0,3,1):
    r = r_inner+ridx*r_delta
    for thetaidx in range(0,8,1):
        theta = theta_orig-thetaidx*theta_delta
        for zidx in range(0,3,1):
            z = zorig+zidx*z_delta
            node_list.append(
                [
                    xorig+r*math.cos(theta),
                    yorig+r*math.sin(theta),
                    zorig+z
                ]
            )
nodesPerR = 24
# Define the element connectivity
elements = iron.MeshElements()
meshComponentNumber=1
elements.CreateStart(mesh,meshComponentNumber,basis)

elem_node_list.append(1,2,3,4,5,6,7,8,9,25,26,27,28,29,30,31,32,33,49,50,51,52,53,54,55,56,57)
elem_node_list.append(7,8,9,10,11,12,13,14,15,31,32,33,34,35,36,37,38,39,55,56,57,58,59,60,61,62,63)
elem_node_list.append(13,14,15,16,17,18,19,20,21,37,38,39,40,41,42,43,44,45,61,62,63,64,65,66,67,68,69)
elem_node_list.append(19,20,21,22,23,24,1,2,3,43,44,45,46,47,48,25,26,27,67,68,69,70,71,72,49,50,51)
elements.NodesSet(1, elem_node_list[0])
elements.NodesSet(2, elem_node_list[1])
elements.NodesSet(3, elem_node_list[2])
elements.NodesSet(4, elem_node_list[3])

mesh.ElementsSet(elemidx+1,iron.ElementsQuadrilateral,elem_node_list)


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
# node 1
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,1,1,0.0)
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,1,2,0.0)
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,1,3,0.0)
# node 2
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,2,1,0.0)
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,2,2,0.0)
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,2,3,1.0)
# node 3
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,3,1,-1.0)
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,3,2,0.0)
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,3,3,0.0)
# node 4
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,4,1, -1.0)
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,4,2,0.0)
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,4,3,1.0)
# node 5
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,5,1,0.5)
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,5,2,-0.5)
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,5,3,0.0)
# node 6
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,6,1,0.5)
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,6,2,-0.5)
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,6,3,1.0)
# node 7
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,7,1,-1.5)
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,7,2,-0.5)
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,7,3,0.0)
# node 8
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,8,1,-1.5)
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,8,2,-0.5)
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,8,3,1.0)

 # node 9
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,9,1,0.0)
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,9,2,1.0)
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,9,3,0.0)
 # node 10
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,10,1,0.0)
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,10,2,1.0)
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,10,3,1.0)
 # node 11
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,11,1,0.5)
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,11,2,1.5)
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,11,3,0.0)
 # node 12
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,12,1,0.5)
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,12,2,1.5)
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,12,3,1.0)
# # node 13
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,13,1,-1.0)
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,13,2,1.0)
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,13,3,0.0)
# node 14
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,14,1,-1.0)
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,14,2,1.0)
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,14,3,1.0)
# node 15
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,15,1,-1.5)
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,15,2,1.5)
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,15,3,0.0)
# node 16
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,16,1,-1.5)
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,16,2,1.5)
geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,16,3,1.0)
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
    materialField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,2.0)
iron.Field.ComponentValuesInitialiseDP(
    materialField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,6.0)

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

boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,7,1,iron.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,8,1,iron.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,15,1,iron.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,16,1,iron.BoundaryConditionsTypes.FIXED,0.0)

# Set y=0 nodes to no y displacement
boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,5,2,iron.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,6,2,iron.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,7,2,iron.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,8,2,iron.BoundaryConditionsTypes.FIXED,0.0)

# Set z=0 nodes to no y displacement
boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,1,3,iron.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,3,3,iron.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,5,3,iron.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,7,3,iron.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,9,3,iron.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,11,3,iron.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,13,3,iron.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,15,3,iron.BoundaryConditionsTypes.FIXED,0.0)

#Pressure boundary conditions on the internal faces.
boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.DELUDELN,1,1,1,3,iron.BoundaryConditionsTypes.PRESSURE_INCREMENTED,10.0)
boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.DELUDELN,1,1,2,3,iron.BoundaryConditionsTypes.PRESSURE_INCREMENTED,10.0)
#boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.DELUDELN,1,1,3,3,iron.BoundaryConditionsTypes.PRESSURE_INCREMENTED,0.5)
#boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.DELUDELN,1,1,4,3,iron.BoundaryConditionsTypes.PRESSURE_INCREMENTED,0.5)
boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.DELUDELN,1,1,9,3,iron.BoundaryConditionsTypes.PRESSURE_INCREMENTED,10.0)
boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.DELUDELN,1,1,10,3,iron.BoundaryConditionsTypes.PRESSURE_INCREMENTED,10.0)
#boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.DELUDELN,1,1,13,3,iron.BoundaryConditionsTypes.PRESSURE_INCREMENTED,0.0)
#boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.DELUDELN,1,1,14,3,iron.BoundaryConditionsTypes.PRESSURE_INCREMENTED,0.0)


solverEquations.BoundaryConditionsCreateFinish()

# Solve the problem
problem.Solve()

vtkoutput(totalNumberOfNodes,totalNumberOfElements,mesh,geometricField,dependentField)
# Export results
fields = iron.Fields()
fields.CreateRegion(region)
fields.NodesExport("UniAxialExtension","FORTRAN")
fields.ElementsExport("UniAxialExtension","FORTRAN")
fields.Finalise()

