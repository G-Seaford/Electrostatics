PROGRAM MAIN

    USE CommonDataStructure
    USE SolveGaussSeidel

    IMPLICIT NONE

    TYPE(CommonData)    :: Simulation_Data
    INTEGER             :: status

    CALL InitialiseCommonData(Simulation_Data, status)

    CALL SolveSystem(Simulation_Data)

    CALL CleanUpData(Simulation_Data)

END PROGRAM MAIN