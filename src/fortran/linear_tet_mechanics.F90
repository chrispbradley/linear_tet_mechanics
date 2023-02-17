!> Main program
PROGRAM LinearTetMechanicsExample

  USE OpenCMISS
  USE OpenCMISS_Iron
#ifndef NOMPIMOD
  USE MPI
#endif

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif

  !Test program parameters

  INTEGER(CMISSIntg), PARAMETER :: CONTEXT_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: COORDINATE_SYSTEM_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: NUMBER_OF_SPATIAL_COORDINATES=3
  INTEGER(CMISSIntg), PARAMETER :: REGION_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: BASIS_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: MESH_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: DECOMPOSITION_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: DECOMPOSER_USER_NUMBER=1

  INTEGER(CMISSIntg), PARAMETER :: NUMBER_OF_XI_COORDINATES=3
  INTEGER(CMISSIntg), PARAMETER :: TOTAL_NUMBER_OF_NODES=4
  INTEGER(CMISSIntg), PARAMETER :: NUMBER_OF_MESH_DIMENSIONS=3
  INTEGER(CMISSIntg), PARAMETER :: NUMBER_OF_MESH_COMPONENTS=1
  INTEGER(CMISSIntg), PARAMETER :: TOTAL_NUMBER_OF_ELEMENTS=1
  INTEGER(CMISSIntg), PARAMETER :: MESH_COMPONENT_NUMBER=1

  INTEGER(CMISSIntg), PARAMETER :: FIELD_GEOMETRY_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: FIELD_GEOMETRY_NUMBER_OF_VARIABLES=1
  INTEGER(CMISSIntg), PARAMETER :: FIELD_GEOMETRY_NUMBER_OF_COMPONENTS=3

  INTEGER(CMISSIntg), PARAMETER :: FIELD_FIBRE_USER_NUMBER=2
  INTEGER(CMISSIntg), PARAMETER :: FIELD_FIBRE_NUMBER_OF_VARIABLES=1
  INTEGER(CMISSIntg), PARAMETER :: FIELD_FIBRE_NUMBER_OF_COMPONENTS=3

  INTEGER(CMISSIntg), PARAMETER :: FIELD_MATERIAL_USER_NUMBER=3
  INTEGER(CMISSIntg), PARAMETER :: FIELD_MATERIAL_NUMBER_OF_VARIABLES=1
  INTEGER(CMISSIntg), PARAMETER :: FIELD_MATERIAL_NUMBER_OF_COMPONENTS=2

  INTEGER(CMISSIntg), PARAMETER :: FIELD_DEPENDENT_USER_NUMBER=4
  INTEGER(CMISSIntg), PARAMETER :: FIELD_DEPENDENT_NUMBER_OF_VARIABLES=2
  INTEGER(CMISSIntg), PARAMETER :: FIELD_DEPENDENT_NUMBER_OF_COMPONENTS=4

  INTEGER(CMISSIntg), PARAMETER :: EQUATIONS_SET_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: EQUATIONS_SET_FIELD_USER_NUMBER=13
  INTEGER(CMISSIntg), PARAMETER :: PROBLEM_USER_NUMBER=1

  !Program types

  !Program variables
  LOGICAL  :: directoryExists = .FALSE.

  INTEGER(CMISSIntg) :: numberGlobalXElements,numberGlobalYElements,numberGlobalZElements
  INTEGER(CMISSIntg) :: numberOfComputationalNodes,computationalNodeNumber
  INTEGER(CMISSIntg) :: decompositionIndex,equationsSetIndex  

  !OpenCMISS variables
  TYPE(cmfe_BasisType) :: basis
  TYPE(cmfe_BoundaryConditionsType) :: boundaryConditions
  TYPE(cmfe_ComputationEnvironmentType) :: computationEnvironment
  TYPE(cmfe_ContextType) :: context
  TYPE(cmfe_CoordinateSystemType) :: coordinateSystem
  TYPE(cmfe_MeshType) :: mesh
  TYPE(cmfe_DecompositionType) :: decomposition
  TYPE(cmfe_DecomposerType) :: decomposer
  TYPE(cmfe_EquationsType) :: equations
  TYPE(cmfe_EquationsSetType) :: equationsSet
  TYPE(cmfe_FieldType) :: geometricField,fibreField,materialField,dependentField,equationsSetField
  TYPE(cmfe_FieldsType) :: fields
  TYPE(cmfe_MeshElementsType) :: meshElements
  TYPE(cmfe_NodesType) :: nodes
  TYPE(cmfe_ProblemType) :: problem
  TYPE(cmfe_RegionType) :: region,worldRegion
  TYPE(cmfe_SolverType) :: solver,linearSolver
  TYPE(cmfe_SolverEquationsType) :: solverEquations
  TYPE(cmfe_WorkGroupType) :: worldWorkGroup

  !Generic CMISS variables
  INTEGER(CMISSIntg) :: err

  !Intialise OpenCMISS
  CALL cmfe_Initialise(err)
  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,err)
  !Set all diganostic levels on for testing
  CALL cmfe_DiagnosticsSetOn(CMFE_FROM_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["FiniteElasticity_FiniteElementResidualEvaluate"],err) !CMFE_ALL_DIAG_TYPE
  !Create a context
  CALL cmfe_Context_Initialise(context,err)
  CALL cmfe_Context_Create(CONTEXT_USER_NUMBER,context,err)
  CALL cmfe_Region_Initialise(worldRegion,err)
  CALL cmfe_Context_WorldRegionGet(context,worldRegion,err)

  WRITE(*,'(A)') "Program starting."

  !Get the number of computational nodes and this computational node number
  CALL cmfe_ComputationEnvironment_Initialise(computationEnvironment,err)
  CALL cmfe_Context_ComputationEnvironmentGet(context,computationEnvironment,err)
  
  CALL cmfe_WorkGroup_Initialise(worldWorkGroup,err)
  CALL cmfe_ComputationEnvironment_WorldWorkGroupGet(computationEnvironment,worldWorkGroup,err)
  CALL cmfe_WorkGroup_NumberOfGroupNodesGet(worldWorkGroup,numberOfComputationalNodes,err)
  CALL cmfe_WorkGroup_GroupNodeNumberGet(worldWorkGroup,computationalNodeNumber,err)

  numberGlobalXElements=1
  numberGlobalYElements=1
  numberGlobalZElements=1   

  !Create a CS - default is 3D rectangular cartesian CS with 0,0,0 as origin
  CALL cmfe_CoordinateSystem_Initialise(coordinateSystem,err)
  CALL cmfe_CoordinateSystem_CreateStart(COORDINATE_SYSTEM_USER_NUMBER,context,coordinateSystem,err)
  CALL cmfe_CoordinateSystem_TypeSet(coordinateSystem,CMFE_COORDINATE_RECTANGULAR_CARTESIAN_TYPE,err)
  CALL cmfe_CoordinateSystem_DimensionSet(coordinateSystem,NUMBER_OF_SPATIAL_COORDINATES,err)
  CALL cmfe_CoordinateSystem_OriginSet(coordinateSystem,[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP],err)
  CALL cmfe_CoordinateSystem_CreateFinish(coordinateSystem,err)

  !Create a region and assign the CS to the region
  CALL cmfe_Region_Initialise(region,err)
  CALL cmfe_Region_CreateStart(REGION_USER_NUMBER,worldRegion,region,err)
  CALL cmfe_Region_LabelSet(region,"Region",err)
  CALL cmfe_Region_CoordinateSystemSet(region,coordinateSystem,err)
  CALL cmfe_Region_CreateFinish(region,err)

  !Define basis function - Simplex tri-linear
  CALL cmfe_Basis_Initialise(basis,err)
  CALL cmfe_Basis_CreateStart(BASIS_USER_NUMBER,context,basis,err)
  CALL cmfe_Basis_TypeSet(basis,CMFE_BASIS_SIMPLEX_TYPE,err)
  CALL cmfe_Basis_NumberOfXiSet(basis,NUMBER_OF_XI_COORDINATES,err)
  CALL cmfe_Basis_InterpolationXiSet(basis,[CMFE_BASIS_LINEAR_SIMPLEX_INTERPOLATION, &
    & CMFE_BASIS_LINEAR_SIMPLEX_INTERPOLATION,CMFE_BASIS_LINEAR_SIMPLEX_INTERPOLATION],err)
  CALL cmfe_Basis_CreateFinish(basis,err)

  !Create a mesh
  CALL cmfe_Mesh_Initialise(mesh,err)
  CALL cmfe_Mesh_CreateStart(MESH_USER_NUMBER,region,NUMBER_OF_MESH_DIMENSIONS,mesh,err)

  CALL cmfe_Mesh_NumberOfComponentsSet(mesh,NUMBER_OF_MESH_COMPONENTS,err)
  CALL cmfe_Mesh_NumberOfElementsSet(mesh,TOTAL_NUMBER_OF_ELEMENTS,err)

  !Define nodes for the mesh
  CALL cmfe_Nodes_Initialise(nodes,err)
  CALL cmfe_Nodes_CreateStart(region,TOTAL_NUMBER_OF_NODES,nodes,err)
  CALL cmfe_Nodes_CreateFinish(nodes,err)

  CALL cmfe_MeshElements_Initialise(meshElements,err)
  CALL cmfe_MeshElements_CreateStart(mesh,MESH_COMPONENT_NUMBER,basis,meshElements,err)
  CALL cmfe_MeshElements_NodesSet(meshElements,1,[1,2,3,4],err)
  CALL cmfe_MeshElements_CreateFinish(meshElements,err)

  CALL cmfe_Mesh_CreateFinish(mesh,err)

  !Create a decomposition
  CALL cmfe_Decomposition_Initialise(decomposition,err)
  CALL cmfe_Decomposition_CreateStart(DECOMPOSITION_USER_NUMBER,mesh,decomposition,err)
  CALL cmfe_Decomposition_CreateFinish(decomposition,err)

  !Decompose
  CALL cmfe_Decomposer_Initialise(decomposer,err)
  CALL cmfe_Decomposer_CreateStart(DECOMPOSER_USER_NUMBER,region,worldWorkGroup,decomposer,err)
  !Add in the decomposition
  CALL cmfe_Decomposer_DecompositionAdd(decomposer,decomposition,decompositionIndex,err)
  !Finish the decomposer
  CALL cmfe_Decomposer_CreateFinish(decomposer,err)
  
  !Create a field to put the geometry (default is geometry)
  CALL cmfe_Field_Initialise(geometricField,err)
  CALL cmfe_Field_CreateStart(FIELD_GEOMETRY_USER_NUMBER,region,geometricField,err)
  CALL cmfe_Field_DecompositionSet(geometricField,decomposition,err)
  CALL cmfe_Field_TypeSet(geometricField,CMFE_FIELD_GEOMETRIC_TYPE,err)
  CALL cmfe_Field_NumberOfVariablesSet(geometricField,FIELD_GEOMETRY_NUMBER_OF_VARIABLES,err)
  CALL cmfe_Field_VariableLabelSet(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,"Geometry",err)
  CALL cmfe_Field_NumberOfComponentsSet(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,FIELD_GEOMETRY_NUMBER_OF_COMPONENTS,err)
  CALL cmfe_Field_ComponentMeshComponentSet(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,1,MESH_COMPONENT_NUMBER,err)
  CALL cmfe_Field_ComponentMeshComponentSet(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,2,MESH_COMPONENT_NUMBER,err)
  CALL cmfe_Field_ComponentMeshComponentSet(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,3,MESH_COMPONENT_NUMBER,err)
  CALL cmfe_Field_CreateFinish(geometricField,err)

  !node 1
  CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,1,1, &
    & 0.0_CMISSRP,err)
  CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,1,2, &
    & 0.0_CMISSRP,err)
  CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,1,3, &
    & 0.0_CMISSRP,err)
  !node 2
  CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,2,1, &
    & 1.0_CMISSRP,err)
  CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,2,2, &
    & 0.0_CMISSRP,err)
  CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,2,3, &
    & 0.0_CMISSRP,err)
  !node 3
  CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,3,1, &
    & 0.0_CMISSRP,err)
  CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,3,2, &
    & 1.0_CMISSRP,err)
  CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,3,3, &
    & 0.0_CMISSRP,err)
  !node 4
  CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,4,1, &
    & 0.0_CMISSRP,err)
  CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,4,2, &
    & 0.0_CMISSRP,err)
  CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,4,3, &
    & 1.0_CMISSRP,err)

  !Create a fibre field and attach it to the geometric field
  CALL cmfe_Field_Initialise(fibreField,err)
  CALL cmfe_Field_CreateStart(FIELD_FIBRE_USER_NUMBER,region,fibreField,err)
  CALL cmfe_Field_TypeSet(fibreField,CMFE_FIELD_FIBRE_TYPE,err)
  CALL cmfe_Field_DecompositionSet(fibreField,decomposition,err)
  CALL cmfe_Field_GeometricFieldSet(fibreField,geometricField,err)
  CALL cmfe_Field_NumberOfVariablesSet(fibreField,FIELD_FIBRE_NUMBER_OF_VARIABLES,err)
  CALL cmfe_Field_VariableLabelSet(fibreField,CMFE_FIELD_U_VARIABLE_TYPE,"Fibre",err)
  CALL cmfe_Field_NumberOfComponentsSet(fibreField,CMFE_FIELD_U_VARIABLE_TYPE,FIELD_FIBRE_NUMBER_OF_COMPONENTS,err)
  CALL cmfe_Field_ComponentMeshComponentSet(fibreField,CMFE_FIELD_U_VARIABLE_TYPE,1,MESH_COMPONENT_NUMBER,err)
  CALL cmfe_Field_ComponentMeshComponentSet(fibreField,CMFE_FIELD_U_VARIABLE_TYPE,2,MESH_COMPONENT_NUMBER,err)
  CALL cmfe_Field_ComponentMeshComponentSet(fibreField,CMFE_FIELD_U_VARIABLE_TYPE,3,MESH_COMPONENT_NUMBER,err)
  CALL cmfe_Field_CreateFinish(fibreField,err)

  !Create a material field and attach it to the geometric field
  CALL cmfe_Field_Initialise(materialField,err)
  CALL cmfe_Field_CreateStart(FIELD_MATERIAL_USER_NUMBER,region,materialField,err)
  CALL cmfe_Field_TypeSet(materialField,CMFE_FIELD_MATERIAL_TYPE,err)
  CALL cmfe_Field_DecompositionSet(materialField,decomposition,err)
  CALL cmfe_Field_GeometricFieldSet(materialField,geometricField,err)
  CALL cmfe_Field_NumberOfVariablesSet(materialField,FIELD_MATERIAL_NUMBER_OF_VARIABLES,err)
  CALL cmfe_Field_VariableLabelSet(materialField,CMFE_FIELD_U_VARIABLE_TYPE,"Material",err)
  CALL cmfe_Field_NumberOfComponentsSet(materialField,CMFE_FIELD_U_VARIABLE_TYPE,FIELD_MATERIAL_NUMBER_OF_COMPONENTS,err)
  CALL cmfe_Field_ComponentMeshComponentSet(materialField,CMFE_FIELD_U_VARIABLE_TYPE,1,MESH_COMPONENT_NUMBER,err)
  CALL cmfe_Field_ComponentMeshComponentSet(materialField,CMFE_FIELD_U_VARIABLE_TYPE,2,MESH_COMPONENT_NUMBER,err)
  CALL cmfe_Field_CreateFinish(materialField,err)

  !Set Mooney-Rivlin constants c10 and c01 to 2.0 and 6.0 respectively.
  CALL cmfe_Field_ComponentValuesInitialise(materialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,2.0_CMISSRP,err)
  CALL cmfe_Field_ComponentValuesInitialise(materialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,6.0_CMISSRP,err)

  !Create a dependent field with two variables and four components
  CALL cmfe_Field_Initialise(dependentField,err)
  CALL cmfe_Field_CreateStart(FIELD_DEPENDENT_USER_NUMBER,region,dependentField,err)
  CALL cmfe_Field_TypeSet(dependentField,CMFE_FIELD_GEOMETRIC_GENERAL_TYPE,err)
  CALL cmfe_Field_DecompositionSet(dependentField,decomposition,err)
  CALL cmfe_Field_GeometricFieldSet(dependentField,geometricField,err)
  CALL cmfe_Field_DependentTypeSet(dependentField,CMFE_FIELD_DEPENDENT_TYPE,err)
  CALL cmfe_Field_NumberOfVariablesSet(dependentField,FIELD_DEPENDENT_NUMBER_OF_VARIABLES,err)
  CALL cmfe_Field_NumberOfComponentsSet(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,FIELD_DEPENDENT_NUMBER_OF_COMPONENTS,err)
  CALL cmfe_Field_NumberOfComponentsSet(dependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DEPENDENT_NUMBER_OF_COMPONENTS,err)
  CALL cmfe_Field_VariableLabelSet(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,"Dependent",err)
  CALL cmfe_Field_ComponentMeshComponentSet(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,MESH_COMPONENT_NUMBER,err)
  CALL cmfe_Field_ComponentMeshComponentSet(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,2,MESH_COMPONENT_NUMBER,err)
  CALL cmfe_Field_ComponentMeshComponentSet(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,3,MESH_COMPONENT_NUMBER,err)  
  CALL cmfe_Field_ComponentInterpolationSet(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,4,CMFE_FIELD_ELEMENT_BASED_INTERPOLATION, &
    & err)
  CALL cmfe_Field_ComponentMeshComponentSet(dependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,MESH_COMPONENT_NUMBER,err)
  CALL cmfe_Field_ComponentMeshComponentSet(dependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,2,MESH_COMPONENT_NUMBER,err)
  CALL cmfe_Field_ComponentMeshComponentSet(dependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,3,MESH_COMPONENT_NUMBER,err)  
  CALL cmfe_Field_ComponentInterpolationSet(dependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4, &
    & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,err)
  CALL cmfe_Field_CreateFinish(dependentField,err)  

  !Create the equations_set
  CALL cmfe_Field_Initialise(equationsSetField,err)
  CALL cmfe_EquationsSet_Initialise(equationsSet,err)
  CALL cmfe_EquationsSet_CreateStart(EQUATIONS_SET_USER_NUMBER,region,fibreField,[CMFE_EQUATIONS_SET_ELASTICITY_CLASS, &
    & CMFE_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMFE_EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE],EQUATIONS_SET_FIELD_USER_NUMBER, &
    & equationsSetField,equationsSet,err)
  CALL cmfe_EquationsSet_CreateFinish(equationsSet,err)

  CALL cmfe_EquationsSet_DependentCreateStart(equationsSet,FIELD_DEPENDENT_USER_NUMBER,dependentField,err) 
  CALL cmfe_EquationsSet_DependentCreateFinish(equationsSet,err)

  CALL cmfe_EquationsSet_MaterialsCreateStart(equationsSet,FIELD_MATERIAL_USER_NUMBER,materialField,err)  
  CALL cmfe_EquationsSet_MaterialsCreateFinish(equationsSet,err)

  !Create the equations set equations
  CALL cmfe_Equations_Initialise(equations,err)
  CALL cmfe_EquationsSet_EquationsCreateStart(equationsSet,equations,err)
  CALL cmfe_Equations_SparsityTypeSet(equations,CMFE_EQUATIONS_SPARSE_MATRICES,err)
  CALL cmfe_Equations_OutputTypeSet(equations,CMFE_EQUATIONS_NO_OUTPUT,err)
  !CALL cmfe_Equations_OutputTypeSet(equations,CMFE_EQUATIONS_MATRIX_OUTPUT,err)
  CALL cmfe_EquationsSet_EquationsCreateFinish(equationsSet,err)   

  !Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,err)
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 2,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,err)
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 3,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,err)
  CALL cmfe_Field_ComponentValuesInitialise(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4,-8.0_CMISSRP, &
    & err)

  !Define the problem
  CALL cmfe_Problem_Initialise(problem,err)
  CALL cmfe_Problem_CreateStart(PROBLEM_USER_NUMBER,context,[CMFE_PROBLEM_ELASTICITY_CLASS,CMFE_PROBLEM_FINITE_ELASTICITY_TYPE, &
    & CMFE_PROBLEM_STATIC_FINITE_ELASTICITY_SUBTYPE],problem,err)
  CALL cmfe_Problem_CreateFinish(problem,err)

  !Create the problem control loop
  CALL cmfe_Problem_ControlLoopCreateStart(problem,err)
  CALL cmfe_Problem_ControlLoopCreateFinish(problem,err)

  !Create the problem solvers
  CALL cmfe_Solver_Initialise(solver,err)
  CALL cmfe_Solver_Initialise(linearSolver,err)
  CALL cmfe_Problem_SolversCreateStart(problem,err)
  CALL cmfe_Problem_SolverGet(problem,CMFE_CONTROL_LOOP_NODE,1,solver,err)
  !CALL cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_PROGRESS_OUTPUT,err)
  CALL cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_MATRIX_OUTPUT,err)
  CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(solver,CMFE_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,err)
  CALL cmfe_Solver_NewtonLinearSolverGet(solver,linearSolver,err)
  CALL cmfe_Solver_LinearTypeSet(linearSolver,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,err)
  CALL cmfe_Problem_SolversCreateFinish(problem,err)

  !Create the problem solver equations
  CALL cmfe_Solver_Initialise(solver,err)
  CALL cmfe_SolverEquations_Initialise(solverEquations,err)
  CALL cmfe_Problem_SolverEquationsCreateStart(problem,err)   
  CALL cmfe_Problem_SolverGet(problem,CMFE_CONTROL_LOOP_NODE,1,solver,err)
  CALL cmfe_Solver_SolverEquationsGet(solver,solverEquations,err)
  CALL cmfe_SolverEquations_SparsityTypeSet(solverEquations,CMFE_SOLVER_SPARSE_MATRICES,err)
  CALL cmfe_SolverEquations_EquationsSetAdd(solverEquations,equationsSet,equationsSetIndex,err)
  CALL cmfe_Problem_SolverEquationsCreateFinish(problem,err)

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL cmfe_BoundaryConditions_Initialise(boundaryConditions,err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(solverEquations,boundaryConditions,err)

  !Fix nodes 1,3,4 at x=0 and nodes 2 at x=1.1
  CALL cmfe_BoundaryConditions_SetNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,1,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,err)
  CALL cmfe_BoundaryConditions_SetNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,2,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED,1.1_CMISSRP,err)
  CALL cmfe_BoundaryConditions_SetNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,3,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,err)
  CALL cmfe_BoundaryConditions_SetNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,4,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,err)

  !Fix nodes 1,2,4 at y=0
  CALL cmfe_BoundaryConditions_SetNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,1,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,err)
  CALL cmfe_BoundaryConditions_SetNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,2,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,err)
  CALL cmfe_BoundaryConditions_SetNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,4,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,err)

  !Fix nodes 1,2,3 at z=0
  CALL cmfe_BoundaryConditions_SetNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,1,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,err)
  CALL cmfe_BoundaryConditions_SetNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,2,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,err)
  CALL cmfe_BoundaryConditions_SetNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,3,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,err)

  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(solverEquations,err)

  !Solve problem
  CALL cmfe_Problem_Solve(problem,err)

  INQUIRE(file="./results", exist=directoryExists)
  IF (.NOT.directoryExists) THEN
    CALL EXECUTE_COMMAND_LINE("mkdir ./results")
  ENDIF

  !Output solution  
  CALL cmfe_Fields_Initialise(fields,err)
  CALL cmfe_Fields_Create(region,fields,err)
  CALL cmfe_Fields_NodesExport(fields,"./results/LinearTet","FORTRAN",err)
  CALL cmfe_Fields_ElementsExport(fields,"./results/LinearTet","FORTRAN",err)
  CALL cmfe_Fields_Finalise(fields,err)

  !Destroy the context
  CALL cmfe_Context_Destroy(context,err)
  !Finalise OpenCMISS
  CALL cmfe_Finalise(err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

END PROGRAM LinearTetMechanicsExample

