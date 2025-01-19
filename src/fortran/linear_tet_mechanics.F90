!> Main program
PROGRAM LinearTetMechanicsExample

  USE OpenCMISS

  IMPLICIT NONE

  !Test program parameters

  INTEGER(OC_Intg), PARAMETER :: CONTEXT_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: COORDINATE_SYSTEM_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: NUMBER_OF_SPATIAL_COORDINATES=3
  INTEGER(OC_Intg), PARAMETER :: REGION_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: BASIS_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: MESH_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: DECOMPOSITION_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: DECOMPOSER_USER_NUMBER=1

  INTEGER(OC_Intg), PARAMETER :: NUMBER_OF_XI_COORDINATES=3
  INTEGER(OC_Intg), PARAMETER :: TOTAL_NUMBER_OF_NODES=4
  INTEGER(OC_Intg), PARAMETER :: NUMBER_OF_MESH_DIMENSIONS=3
  INTEGER(OC_Intg), PARAMETER :: NUMBER_OF_MESH_COMPONENTS=1
  INTEGER(OC_Intg), PARAMETER :: TOTAL_NUMBER_OF_ELEMENTS=1
  INTEGER(OC_Intg), PARAMETER :: MESH_COMPONENT_NUMBER=1

  INTEGER(OC_Intg), PARAMETER :: FIELD_GEOMETRY_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: FIELD_GEOMETRY_NUMBER_OF_VARIABLES=1
  INTEGER(OC_Intg), PARAMETER :: FIELD_GEOMETRY_NUMBER_OF_COMPONENTS=3

  INTEGER(OC_Intg), PARAMETER :: FIELD_FIBRE_USER_NUMBER=2
  INTEGER(OC_Intg), PARAMETER :: FIELD_FIBRE_NUMBER_OF_VARIABLES=1
  INTEGER(OC_Intg), PARAMETER :: FIELD_FIBRE_NUMBER_OF_COMPONENTS=3

  INTEGER(OC_Intg), PARAMETER :: FIELD_MATERIAL_USER_NUMBER=3
  INTEGER(OC_Intg), PARAMETER :: FIELD_MATERIAL_NUMBER_OF_VARIABLES=1
  INTEGER(OC_Intg), PARAMETER :: FIELD_MATERIAL_NUMBER_OF_COMPONENTS=2

  INTEGER(OC_Intg), PARAMETER :: FIELD_DEPENDENT_USER_NUMBER=4
  INTEGER(OC_Intg), PARAMETER :: FIELD_DEPENDENT_NUMBER_OF_VARIABLES=2
  INTEGER(OC_Intg), PARAMETER :: FIELD_DEPENDENT_NUMBER_OF_COMPONENTS=4

  INTEGER(OC_Intg), PARAMETER :: EQUATIONS_SET_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: EQUATIONS_SET_FIELD_USER_NUMBER=13
  INTEGER(OC_Intg), PARAMETER :: PROBLEM_USER_NUMBER=1

  !Program types

  !Program variables
  LOGICAL  :: directoryExists = .FALSE.

  INTEGER(OC_Intg) :: numberGlobalXElements,numberGlobalYElements,numberGlobalZElements
  INTEGER(OC_Intg) :: numberOfComputationalNodes,computationalNodeNumber
  INTEGER(OC_Intg) :: decompositionIndex,equationsSetIndex  

  !OpenCMISS variables
  TYPE(OC_BasisType) :: basis
  TYPE(OC_BoundaryConditionsType) :: boundaryConditions
  TYPE(OC_ComputationEnvironmentType) :: computationEnvironment
  TYPE(OC_ContextType) :: context
  TYPE(OC_CoordinateSystemType) :: coordinateSystem
  TYPE(OC_MeshType) :: mesh
  TYPE(OC_DecompositionType) :: decomposition
  TYPE(OC_DecomposerType) :: decomposer
  TYPE(OC_EquationsType) :: equations
  TYPE(OC_EquationsSetType) :: equationsSet
  TYPE(OC_FieldType) :: geometricField,fibreField,materialField,dependentField,equationsSetField
  TYPE(OC_FieldsType) :: fields
  TYPE(OC_MeshElementsType) :: meshElements
  TYPE(OC_NodesType) :: nodes
  TYPE(OC_ProblemType) :: problem
  TYPE(OC_RegionType) :: region,worldRegion
  TYPE(OC_SolverType) :: solver,linearSolver
  TYPE(OC_SolverEquationsType) :: solverEquations
  TYPE(OC_WorkGroupType) :: worldWorkGroup

  !Generic CMISS variables
  INTEGER(OC_Intg) :: err

  !Intialise OpenCMISS
  CALL OC_Initialise(err)
  CALL OC_ErrorHandlingModeSet(OC_ERRORS_TRAP_ERROR,err)
  !Set all diganostic levels on for testing
  CALL OC_DiagnosticsSetOn(OC_FROM_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["FiniteElasticity_FiniteElementResidualEvaluate"],err) !OC_ALL_DIAG_TYPE
  !Create a context
  CALL OC_Context_Initialise(context,err)
  CALL OC_Context_Create(CONTEXT_USER_NUMBER,context,err)
  CALL OC_Region_Initialise(worldRegion,err)
  CALL OC_Context_WorldRegionGet(context,worldRegion,err)

  WRITE(*,'(A)') "Program starting."

  !Get the number of computational nodes and this computational node number
  CALL OC_ComputationEnvironment_Initialise(computationEnvironment,err)
  CALL OC_Context_ComputationEnvironmentGet(context,computationEnvironment,err)
  
  CALL OC_WorkGroup_Initialise(worldWorkGroup,err)
  CALL OC_ComputationEnvironment_WorldWorkGroupGet(computationEnvironment,worldWorkGroup,err)
  CALL OC_WorkGroup_NumberOfGroupNodesGet(worldWorkGroup,numberOfComputationalNodes,err)
  CALL OC_WorkGroup_GroupNodeNumberGet(worldWorkGroup,computationalNodeNumber,err)

  numberGlobalXElements=1
  numberGlobalYElements=1
  numberGlobalZElements=1   

  !Create a CS - default is 3D rectangular cartesian CS with 0,0,0 as origin
  CALL OC_CoordinateSystem_Initialise(coordinateSystem,err)
  CALL OC_CoordinateSystem_CreateStart(COORDINATE_SYSTEM_USER_NUMBER,context,coordinateSystem,err)
  CALL OC_CoordinateSystem_TypeSet(coordinateSystem,OC_COORDINATE_RECTANGULAR_CARTESIAN_TYPE,err)
  CALL OC_CoordinateSystem_DimensionSet(coordinateSystem,NUMBER_OF_SPATIAL_COORDINATES,err)
  CALL OC_CoordinateSystem_OriginSet(coordinateSystem,[0.0_OC_RP,0.0_OC_RP,0.0_OC_RP],err)
  CALL OC_CoordinateSystem_CreateFinish(coordinateSystem,err)

  !Create a region and assign the CS to the region
  CALL OC_Region_Initialise(region,err)
  CALL OC_Region_CreateStart(REGION_USER_NUMBER,worldRegion,region,err)
  CALL OC_Region_LabelSet(region,"Region",err)
  CALL OC_Region_CoordinateSystemSet(region,coordinateSystem,err)
  CALL OC_Region_CreateFinish(region,err)

  !Define basis function - Simplex tri-linear
  CALL OC_Basis_Initialise(basis,err)
  CALL OC_Basis_CreateStart(BASIS_USER_NUMBER,context,basis,err)
  CALL OC_Basis_TypeSet(basis,OC_BASIS_SIMPLEX_TYPE,err)
  CALL OC_Basis_NumberOfXiSet(basis,NUMBER_OF_XI_COORDINATES,err)
  CALL OC_Basis_InterpolationXiSet(basis,[OC_BASIS_LINEAR_SIMPLEX_INTERPOLATION, &
    & OC_BASIS_LINEAR_SIMPLEX_INTERPOLATION,OC_BASIS_LINEAR_SIMPLEX_INTERPOLATION],err)
  CALL OC_Basis_CreateFinish(basis,err)

  !Create a mesh
  CALL OC_Mesh_Initialise(mesh,err)
  CALL OC_Mesh_CreateStart(MESH_USER_NUMBER,region,NUMBER_OF_MESH_DIMENSIONS,mesh,err)

  CALL OC_Mesh_NumberOfComponentsSet(mesh,NUMBER_OF_MESH_COMPONENTS,err)
  CALL OC_Mesh_NumberOfElementsSet(mesh,TOTAL_NUMBER_OF_ELEMENTS,err)

  !Define nodes for the mesh
  CALL OC_Nodes_Initialise(nodes,err)
  CALL OC_Nodes_CreateStart(region,TOTAL_NUMBER_OF_NODES,nodes,err)
  CALL OC_Nodes_CreateFinish(nodes,err)

  CALL OC_MeshElements_Initialise(meshElements,err)
  CALL OC_MeshElements_CreateStart(mesh,MESH_COMPONENT_NUMBER,basis,meshElements,err)
  CALL OC_MeshElements_NodesSet(meshElements,1,[1,2,3,4],err)
  CALL OC_MeshElements_CreateFinish(meshElements,err)

  CALL OC_Mesh_CreateFinish(mesh,err)

  !Create a decomposition
  CALL OC_Decomposition_Initialise(decomposition,err)
  CALL OC_Decomposition_CreateStart(DECOMPOSITION_USER_NUMBER,mesh,decomposition,err)
  CALL OC_Decomposition_CreateFinish(decomposition,err)

  !Decompose
  CALL OC_Decomposer_Initialise(decomposer,err)
  CALL OC_Decomposer_CreateStart(DECOMPOSER_USER_NUMBER,region,worldWorkGroup,decomposer,err)
  !Add in the decomposition
  CALL OC_Decomposer_DecompositionAdd(decomposer,decomposition,decompositionIndex,err)
  !Finish the decomposer
  CALL OC_Decomposer_CreateFinish(decomposer,err)
  
  !Create a field to put the geometry (default is geometry)
  CALL OC_Field_Initialise(geometricField,err)
  CALL OC_Field_CreateStart(FIELD_GEOMETRY_USER_NUMBER,region,geometricField,err)
  CALL OC_Field_DecompositionSet(geometricField,decomposition,err)
  CALL OC_Field_TypeSet(geometricField,OC_FIELD_GEOMETRIC_TYPE,err)
  CALL OC_Field_NumberOfVariablesSet(geometricField,FIELD_GEOMETRY_NUMBER_OF_VARIABLES,err)
  CALL OC_Field_VariableLabelSet(geometricField,OC_FIELD_U_VARIABLE_TYPE,"Geometry",err)
  CALL OC_Field_NumberOfComponentsSet(geometricField,OC_FIELD_U_VARIABLE_TYPE,FIELD_GEOMETRY_NUMBER_OF_COMPONENTS,err)
  CALL OC_Field_ComponentMeshComponentSet(geometricField,OC_FIELD_U_VARIABLE_TYPE,1,MESH_COMPONENT_NUMBER,err)
  CALL OC_Field_ComponentMeshComponentSet(geometricField,OC_FIELD_U_VARIABLE_TYPE,2,MESH_COMPONENT_NUMBER,err)
  CALL OC_Field_ComponentMeshComponentSet(geometricField,OC_FIELD_U_VARIABLE_TYPE,3,MESH_COMPONENT_NUMBER,err)
  CALL OC_Field_CreateFinish(geometricField,err)

  !node 1
  CALL OC_Field_ParameterSetUpdateNode(geometricField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,1,1,1, &
    & 0.0_OC_RP,err)
  CALL OC_Field_ParameterSetUpdateNode(geometricField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,1,1,2, &
    & 0.0_OC_RP,err)
  CALL OC_Field_ParameterSetUpdateNode(geometricField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,1,1,3, &
    & 0.0_OC_RP,err)
  !node 2
  CALL OC_Field_ParameterSetUpdateNode(geometricField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,1,2,1, &
    & 1.0_OC_RP,err)
  CALL OC_Field_ParameterSetUpdateNode(geometricField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,1,2,2, &
    & 0.0_OC_RP,err)
  CALL OC_Field_ParameterSetUpdateNode(geometricField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,1,2,3, &
    & 0.0_OC_RP,err)
  !node 3
  CALL OC_Field_ParameterSetUpdateNode(geometricField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,1,3,1, &
    & 0.0_OC_RP,err)
  CALL OC_Field_ParameterSetUpdateNode(geometricField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,1,3,2, &
    & 1.0_OC_RP,err)
  CALL OC_Field_ParameterSetUpdateNode(geometricField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,1,3,3, &
    & 0.0_OC_RP,err)
  !node 4
  CALL OC_Field_ParameterSetUpdateNode(geometricField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,1,4,1, &
    & 0.0_OC_RP,err)
  CALL OC_Field_ParameterSetUpdateNode(geometricField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,1,4,2, &
    & 0.0_OC_RP,err)
  CALL OC_Field_ParameterSetUpdateNode(geometricField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,1,4,3, &
    & 1.0_OC_RP,err)

  !Create a fibre field and attach it to the geometric field
  CALL OC_Field_Initialise(fibreField,err)
  CALL OC_Field_CreateStart(FIELD_FIBRE_USER_NUMBER,region,fibreField,err)
  CALL OC_Field_TypeSet(fibreField,OC_FIELD_FIBRE_TYPE,err)
  CALL OC_Field_DecompositionSet(fibreField,decomposition,err)
  CALL OC_Field_GeometricFieldSet(fibreField,geometricField,err)
  CALL OC_Field_NumberOfVariablesSet(fibreField,FIELD_FIBRE_NUMBER_OF_VARIABLES,err)
  CALL OC_Field_VariableLabelSet(fibreField,OC_FIELD_U_VARIABLE_TYPE,"Fibre",err)
  CALL OC_Field_NumberOfComponentsSet(fibreField,OC_FIELD_U_VARIABLE_TYPE,FIELD_FIBRE_NUMBER_OF_COMPONENTS,err)
  CALL OC_Field_ComponentMeshComponentSet(fibreField,OC_FIELD_U_VARIABLE_TYPE,1,MESH_COMPONENT_NUMBER,err)
  CALL OC_Field_ComponentMeshComponentSet(fibreField,OC_FIELD_U_VARIABLE_TYPE,2,MESH_COMPONENT_NUMBER,err)
  CALL OC_Field_ComponentMeshComponentSet(fibreField,OC_FIELD_U_VARIABLE_TYPE,3,MESH_COMPONENT_NUMBER,err)
  CALL OC_Field_CreateFinish(fibreField,err)

  !Create a material field and attach it to the geometric field
  CALL OC_Field_Initialise(materialField,err)
  CALL OC_Field_CreateStart(FIELD_MATERIAL_USER_NUMBER,region,materialField,err)
  CALL OC_Field_TypeSet(materialField,OC_FIELD_MATERIAL_TYPE,err)
  CALL OC_Field_DecompositionSet(materialField,decomposition,err)
  CALL OC_Field_GeometricFieldSet(materialField,geometricField,err)
  CALL OC_Field_NumberOfVariablesSet(materialField,FIELD_MATERIAL_NUMBER_OF_VARIABLES,err)
  CALL OC_Field_VariableLabelSet(materialField,OC_FIELD_U_VARIABLE_TYPE,"Material",err)
  CALL OC_Field_NumberOfComponentsSet(materialField,OC_FIELD_U_VARIABLE_TYPE,FIELD_MATERIAL_NUMBER_OF_COMPONENTS,err)
  CALL OC_Field_ComponentMeshComponentSet(materialField,OC_FIELD_U_VARIABLE_TYPE,1,MESH_COMPONENT_NUMBER,err)
  CALL OC_Field_ComponentMeshComponentSet(materialField,OC_FIELD_U_VARIABLE_TYPE,2,MESH_COMPONENT_NUMBER,err)
  CALL OC_Field_CreateFinish(materialField,err)

  !Set Mooney-Rivlin constants c10 and c01 to 2.0 and 6.0 respectively.
  CALL OC_Field_ComponentValuesInitialise(materialField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,2.0_OC_RP,err)
  CALL OC_Field_ComponentValuesInitialise(materialField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,2,6.0_OC_RP,err)

  !Create a dependent field with two variables and four components
  CALL OC_Field_Initialise(dependentField,err)
  CALL OC_Field_CreateStart(FIELD_DEPENDENT_USER_NUMBER,region,dependentField,err)
  CALL OC_Field_TypeSet(dependentField,OC_FIELD_GEOMETRIC_GENERAL_TYPE,err)
  CALL OC_Field_DecompositionSet(dependentField,decomposition,err)
  CALL OC_Field_GeometricFieldSet(dependentField,geometricField,err)
  CALL OC_Field_DependentTypeSet(dependentField,OC_FIELD_DEPENDENT_TYPE,err)
  CALL OC_Field_NumberOfVariablesSet(dependentField,FIELD_DEPENDENT_NUMBER_OF_VARIABLES,err)
  CALL OC_Field_VariableTypesSet(dependentField,[OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_T_VARIABLE_TYPE],err)
  CALL OC_Field_NumberOfComponentsSet(dependentField,OC_FIELD_U_VARIABLE_TYPE,FIELD_DEPENDENT_NUMBER_OF_COMPONENTS,err)
  CALL OC_Field_NumberOfComponentsSet(dependentField,OC_FIELD_T_VARIABLE_TYPE,FIELD_DEPENDENT_NUMBER_OF_COMPONENTS,err)
  CALL OC_Field_VariableLabelSet(dependentField,OC_FIELD_U_VARIABLE_TYPE,"Dependent",err)
  CALL OC_Field_ComponentMeshComponentSet(dependentField,OC_FIELD_U_VARIABLE_TYPE,1,MESH_COMPONENT_NUMBER,err)
  CALL OC_Field_ComponentMeshComponentSet(dependentField,OC_FIELD_U_VARIABLE_TYPE,2,MESH_COMPONENT_NUMBER,err)
  CALL OC_Field_ComponentMeshComponentSet(dependentField,OC_FIELD_U_VARIABLE_TYPE,3,MESH_COMPONENT_NUMBER,err)  
  CALL OC_Field_ComponentInterpolationSet(dependentField,OC_FIELD_U_VARIABLE_TYPE,4,OC_FIELD_ELEMENT_BASED_INTERPOLATION, &
    & err)
  CALL OC_Field_ComponentMeshComponentSet(dependentField,OC_FIELD_T_VARIABLE_TYPE,1,MESH_COMPONENT_NUMBER,err)
  CALL OC_Field_ComponentMeshComponentSet(dependentField,OC_FIELD_T_VARIABLE_TYPE,2,MESH_COMPONENT_NUMBER,err)
  CALL OC_Field_ComponentMeshComponentSet(dependentField,OC_FIELD_T_VARIABLE_TYPE,3,MESH_COMPONENT_NUMBER,err)  
  CALL OC_Field_ComponentInterpolationSet(dependentField,OC_FIELD_T_VARIABLE_TYPE,4, &
    & OC_FIELD_ELEMENT_BASED_INTERPOLATION,err)
  CALL OC_Field_CreateFinish(dependentField,err)  

  !Create the equations_set
  CALL OC_Field_Initialise(equationsSetField,err)
  CALL OC_EquationsSet_Initialise(equationsSet,err)
  CALL OC_EquationsSet_CreateStart(EQUATIONS_SET_USER_NUMBER,region,fibreField,[OC_EQUATIONS_SET_ELASTICITY_CLASS, &
    & OC_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,OC_EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE],EQUATIONS_SET_FIELD_USER_NUMBER, &
    & equationsSetField,equationsSet,err)
  CALL OC_EquationsSet_CreateFinish(equationsSet,err)

  CALL OC_EquationsSet_DependentCreateStart(equationsSet,FIELD_DEPENDENT_USER_NUMBER,dependentField,err) 
  CALL OC_EquationsSet_DependentCreateFinish(equationsSet,err)

  CALL OC_EquationsSet_MaterialsCreateStart(equationsSet,FIELD_MATERIAL_USER_NUMBER,materialField,err)  
  CALL OC_EquationsSet_MaterialsCreateFinish(equationsSet,err)

  !Create the equations set equations
  CALL OC_Equations_Initialise(equations,err)
  CALL OC_EquationsSet_EquationsCreateStart(equationsSet,equations,err)
  CALL OC_Equations_SparsityTypeSet(equations,OC_EQUATIONS_SPARSE_MATRICES,err)
  CALL OC_Equations_OutputTypeSet(equations,OC_EQUATIONS_NO_OUTPUT,err)
  !CALL OC_Equations_OutputTypeSet(equations,OC_EQUATIONS_MATRIX_OUTPUT,err)
  CALL OC_EquationsSet_EquationsCreateFinish(equationsSet,err)   

  !Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
  CALL OC_Field_ParametersToFieldParametersComponentCopy(geometricField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
    & 1,dependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,err)
  CALL OC_Field_ParametersToFieldParametersComponentCopy(geometricField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
    & 2,dependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,2,err)
  CALL OC_Field_ParametersToFieldParametersComponentCopy(geometricField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
    & 3,dependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,3,err)
  CALL OC_Field_ComponentValuesInitialise(dependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,4,-8.0_OC_RP, &
    & err)

  !Define the problem
  CALL OC_Problem_Initialise(problem,err)
  CALL OC_Problem_CreateStart(PROBLEM_USER_NUMBER,context,[OC_PROBLEM_ELASTICITY_CLASS,OC_PROBLEM_FINITE_ELASTICITY_TYPE, &
    & OC_PROBLEM_STATIC_FINITE_ELASTICITY_SUBTYPE],problem,err)
  CALL OC_Problem_CreateFinish(problem,err)

  !Create the problem control loop
  CALL OC_Problem_ControlLoopCreateStart(problem,err)
  CALL OC_Problem_ControlLoopCreateFinish(problem,err)

  !Create the problem solvers
  CALL OC_Solver_Initialise(solver,err)
  CALL OC_Solver_Initialise(linearSolver,err)
  CALL OC_Problem_SolversCreateStart(problem,err)
  CALL OC_Problem_SolverGet(problem,OC_CONTROL_LOOP_NODE,1,solver,err)
  !CALL OC_Solver_OutputTypeSet(solver,OC_SOLVER_PROGRESS_OUTPUT,err)
  CALL OC_Solver_OutputTypeSet(solver,OC_SOLVER_MATRIX_OUTPUT,err)
  CALL OC_Solver_NewtonJacobianCalculationTypeSet(solver,OC_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED,err)
  CALL OC_Solver_NewtonLinearSolverGet(solver,linearSolver,err)
  CALL OC_Solver_LinearTypeSet(linearSolver,OC_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,err)
  CALL OC_Problem_SolversCreateFinish(problem,err)

  !Create the problem solver equations
  CALL OC_Solver_Initialise(solver,err)
  CALL OC_SolverEquations_Initialise(solverEquations,err)
  CALL OC_Problem_SolverEquationsCreateStart(problem,err)   
  CALL OC_Problem_SolverGet(problem,OC_CONTROL_LOOP_NODE,1,solver,err)
  CALL OC_Solver_SolverEquationsGet(solver,solverEquations,err)
  CALL OC_SolverEquations_SparsityTypeSet(solverEquations,OC_SOLVER_SPARSE_MATRICES,err)
  CALL OC_SolverEquations_EquationsSetAdd(solverEquations,equationsSet,equationsSetIndex,err)
  CALL OC_Problem_SolverEquationsCreateFinish(problem,err)

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL OC_BoundaryConditions_Initialise(boundaryConditions,err)
  CALL OC_SolverEquations_BoundaryConditionsCreateStart(solverEquations,boundaryConditions,err)

  !Fix nodes 1,3,4 at x=0 and nodes 2 at x=1.1
  CALL OC_BoundaryConditions_SetNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE,1,1,1,1, &
    & OC_BOUNDARY_CONDITION_FIXED,0.0_OC_RP,err)
  CALL OC_BoundaryConditions_SetNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE,1,1,2,1, &
    & OC_BOUNDARY_CONDITION_FIXED,1.1_OC_RP,err)
  CALL OC_BoundaryConditions_SetNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE,1,1,3,1, &
    & OC_BOUNDARY_CONDITION_FIXED,0.0_OC_RP,err)
  CALL OC_BoundaryConditions_SetNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE,1,1,4,1, &
    & OC_BOUNDARY_CONDITION_FIXED,0.0_OC_RP,err)

  !Fix nodes 1,2,4 at y=0
  CALL OC_BoundaryConditions_SetNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE,1,1,1,2, &
    & OC_BOUNDARY_CONDITION_FIXED,0.0_OC_RP,err)
  CALL OC_BoundaryConditions_SetNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE,1,1,2,2, &
    & OC_BOUNDARY_CONDITION_FIXED,0.0_OC_RP,err)
  CALL OC_BoundaryConditions_SetNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE,1,1,4,2, &
    & OC_BOUNDARY_CONDITION_FIXED,0.0_OC_RP,err)

  !Fix nodes 1,2,3 at z=0
  CALL OC_BoundaryConditions_SetNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE,1,1,1,3, &
    & OC_BOUNDARY_CONDITION_FIXED,0.0_OC_RP,err)
  CALL OC_BoundaryConditions_SetNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE,1,1,2,3, &
    & OC_BOUNDARY_CONDITION_FIXED,0.0_OC_RP,err)
  CALL OC_BoundaryConditions_SetNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE,1,1,3,3, &
    & OC_BOUNDARY_CONDITION_FIXED,0.0_OC_RP,err)

  CALL OC_SolverEquations_BoundaryConditionsCreateFinish(solverEquations,err)

  !Solve problem
  CALL OC_Problem_Solve(problem,err)

  INQUIRE(file="./results", exist=directoryExists)
  IF (.NOT.directoryExists) THEN
    CALL EXECUTE_COMMAND_LINE("mkdir ./results")
  ENDIF

  !Output solution  
  CALL OC_Fields_Initialise(fields,err)
  CALL OC_Fields_Create(region,fields,err)
  CALL OC_Fields_NodesExport(fields,"./results/LinearTet","FORTRAN",err)
  CALL OC_Fields_ElementsExport(fields,"./results/LinearTet","FORTRAN",err)
  CALL OC_Fields_Finalise(fields,err)

  !Destroy the context
  CALL OC_Context_Destroy(context,err)
  !Finalise OpenCMISS
  CALL OC_Finalise(err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

END PROGRAM LinearTetMechanicsExample

