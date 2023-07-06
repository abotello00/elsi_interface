!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for computing estimates of the bounds of the spectrum of a matrix.
MODULE EigenBoundsModule
  USE DataTypesModule, ONLY : NTREAL, MPINTREAL
  USE LoggingModule, ONLY : EnterSubLog, ExitSubLog, WriteElement, &
       & WriteListElement, WriteHeader
  USE PMatrixMemoryPoolModule, ONLY : MatrixMemoryPool_p, &
       & DestructMatrixMemoryPool
  USE PSMatrixAlgebraModule, ONLY : MatrixMultiply, MatrixNorm, DotMatrix, &
       & IncrementMatrix, ScaleMatrix
  USE PSMatrixModule, ONLY : Matrix_ps, ConstructEmptyMatrix, CopyMatrix, &
       & DestructMatrix, GetMatrixTripletList, FillMatrixFromTripletList
  USE SolverParametersModule, ONLY : SolverParameters_t, PrintParameters, &
       & DestructSolverParameters, ConstructSolverParameters, &
       & CopySolverParameters
  USE TripletListModule, ONLY : TripletList_r, TripletList_c, &
       & AppendToTripletList, DestructTripletList, ConstructTripletList
  USE TripletModule, ONLY : Triplet_r
  USE NTMPIModule
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: GershgorinBounds
  PUBLIC :: PowerBounds
CONTAINS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute a bounds on the minimum and maximum eigenvalue of a matrix.
  !> Uses the Gershgorin theorem.
  SUBROUTINE GershgorinBounds(this, min_value, max_value)
    !> The matrix to compute the min/max of.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> A lower bound on the eigenspectrum.
    REAL(NTREAL), INTENT(OUT) :: min_value
    !> An uppder bound on the eigenspectrum.
    REAL(NTREAL), INTENT(OUT) :: max_value
    !! Local Data
    TYPE(TripletList_r) :: tlist_r
    TYPE(TripletList_c) :: tlist_c
    !! Local Data
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: per_column_min
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: per_column_max
    !! Counters/Temporary
    INTEGER :: II
    INTEGER :: local_column
    INTEGER :: ierr

    IF (this%is_complex) THEN


       !! Allocate Space For Result
       ALLOCATE(per_column_min(this%local_columns))
       ALLOCATE(per_column_max(this%local_columns))

       !! Compute The Local Contribution
       per_column_min = 0
       per_column_max = 0
       CALL GetMatrixTripletList(this, tlist_c)
       DO II = 1, tlist_c%CurrentSize
          local_column = tlist_c%DATA(II)%index_column - &
               & this%start_column + 1
          IF (tlist_c%DATA(II)%index_row .EQ. tlist_c%DATA(II)%index_column) THEN
             per_column_min(local_column) = per_column_min(local_column) + &
                  & REAL(tlist_c%DATA(II)%point_value, KIND = NTREAL)
             per_column_max(local_column) = per_column_max(local_column) + &
                  & REAL(tlist_c%DATA(II)%point_value, KIND = NTREAL)
          ELSE
             per_column_min(local_column) = per_column_min(local_column) - &
                  & ABS(tlist_c%DATA(II)%point_value)
             per_column_max(local_column) = per_column_max(local_column) + &
                  & ABS(tlist_c%DATA(II)%point_value)
          END IF
       END DO

       !! Sum Along Columns
       CALL MPI_Allreduce(MPI_IN_PLACE, per_column_min, SIZE(per_column_min), &
            & MPINTREAL, MPI_SUM, this%process_grid%column_comm, ierr)
       CALL MPI_Allreduce(MPI_IN_PLACE, per_column_max, SIZE(per_column_max), &
            & MPINTREAL, MPI_SUM, this%process_grid%column_comm, ierr)

       min_value = MINVAL(per_column_min)
       max_value = MAXVAL(per_column_max)

       CALL MPI_Allreduce(MPI_IN_PLACE, min_value, 1, MPINTREAL, MPI_MIN, &
            & this%process_grid%row_comm, ierr)
       CALL MPI_Allreduce(MPI_IN_PLACE, max_value, 1, MPINTREAL, MPI_MAX, &
            & this%process_grid%row_comm, ierr)

       CALL DestructTripletList(tlist_c)
       DEALLOCATE(per_column_min)
       DEALLOCATE(per_column_max)

    ELSE


       !! Allocate Space For Result
       ALLOCATE(per_column_min(this%local_columns))
       ALLOCATE(per_column_max(this%local_columns))

       !! Compute The Local Contribution
       per_column_min = 0
       per_column_max = 0
       CALL GetMatrixTripletList(this, tlist_r)
       DO II = 1, tlist_r%CurrentSize
          local_column = tlist_r%DATA(II)%index_column - &
               & this%start_column + 1
          IF (tlist_r%DATA(II)%index_row .EQ. tlist_r%DATA(II)%index_column) THEN
             per_column_min(local_column) = per_column_min(local_column) + &
                  & REAL(tlist_r%DATA(II)%point_value, KIND = NTREAL)
             per_column_max(local_column) = per_column_max(local_column) + &
                  & REAL(tlist_r%DATA(II)%point_value, KIND = NTREAL)
          ELSE
             per_column_min(local_column) = per_column_min(local_column) - &
                  & ABS(tlist_r%DATA(II)%point_value)
             per_column_max(local_column) = per_column_max(local_column) + &
                  & ABS(tlist_r%DATA(II)%point_value)
          END IF
       END DO

       !! Sum Along Columns
       CALL MPI_Allreduce(MPI_IN_PLACE, per_column_min, SIZE(per_column_min), &
            & MPINTREAL, MPI_SUM, this%process_grid%column_comm, ierr)
       CALL MPI_Allreduce(MPI_IN_PLACE, per_column_max, SIZE(per_column_max), &
            & MPINTREAL, MPI_SUM, this%process_grid%column_comm, ierr)

       min_value = MINVAL(per_column_min)
       max_value = MAXVAL(per_column_max)

       CALL MPI_Allreduce(MPI_IN_PLACE, min_value, 1, MPINTREAL, MPI_MIN, &
            & this%process_grid%row_comm, ierr)
       CALL MPI_Allreduce(MPI_IN_PLACE, max_value, 1, MPINTREAL, MPI_MAX, &
            & this%process_grid%row_comm, ierr)

       CALL DestructTripletList(tlist_r)
       DEALLOCATE(per_column_min)
       DEALLOCATE(per_column_max)

    END IF
  END SUBROUTINE GershgorinBounds
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute a bounds on the maximum eigenvalue of a matrix.
  !> Uses The Power Method.
  SUBROUTINE PowerBounds(this, max_value, solver_parameters_in)
    !> The matrix to compute the min/max of.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> An upper bound on the eigenspectrum.
    REAL(NTREAL), INTENT(OUT) :: max_value
    !> The parameters for this calculation.
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: param
    !! Local Data
    TYPE(Matrix_ps) :: vector, vector2
    REAL(NTREAL) :: scale_value
    REAL(NTREAL) :: norm_value
    TYPE(TripletList_r) :: temp_list
    TYPE(Triplet_r) :: temp_triplet
    INTEGER :: II
    TYPE(MatrixMemoryPool_p) :: pool

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       CALL CopySolverParameters(solver_parameters_in, param)
    ELSE
       CALL ConstructSolverParameters(param)
       param%max_iterations = 10
    END IF

    IF (param%be_verbose) THEN
       CALL WriteHeader("Power Bounds Solver")
       CALL EnterSubLog
       CALL PrintParameters(param)
    END IF

    !! Diagonal matrices serve as vectors.
    CALL ConstructEmptyMatrix(vector, this)
    CALL ConstructEmptyMatrix(vector2, this)

    !! Guess Vector
    CALL ConstructTripletList(temp_list)
    IF (this%process_grid%global_rank .EQ. 0) THEN
       temp_triplet%index_row = 1
       temp_triplet%index_column = 1
       temp_triplet%point_value = 1.0_NTREAL
       CALL AppendToTripletList(temp_list, temp_triplet)
    END IF
    CALL FillMatrixFromTripletList(vector, temp_list)

    !! Iterate
    IF (param%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    II = 1
    norm_value = param%converge_diff + 1.0_NTREAL
    DO II = 1, param%max_iterations
       IF (param%be_verbose .AND. II .GT. 1) THEN
          CALL WriteListElement(key = "Convergence", VALUE = norm_value)
       END IF

       !! x = Ax
       CALL MatrixMultiply(this, vector, vector2, &
            & threshold_in = param%threshold, memory_pool_in = pool)
       !! x = x/||x||
       scale_value = 1.0 / MatrixNorm(vector2)
       CALL ScaleMatrix(vector2, scale_value)

       !! Check if Converged
       CALL IncrementMatrix(vector2, vector, alpha_in = -1.0_NTREAL)
       norm_value = MatrixNorm(vector)

       CALL CopyMatrix(vector2, vector)

       IF (norm_value .LE. param%converge_diff) THEN
          EXIT
       END IF
    END DO
    IF (param%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key = "Total Iterations", VALUE = II - 1)
    END IF

    !! Compute The Largest Eigenvalue
    CALL DotMatrix(vector, vector, scale_value)
    CALL MatrixMultiply(this, vector, vector2, &
         & threshold_in = param%threshold, memory_pool_in = pool)
    CALL DotMatrix(vector, vector2, max_value)
    max_value = max_value / scale_value

    IF (param%be_verbose) THEN
       CALL WriteElement(key = "Max Eigen Value", VALUE = max_value)
       CALL ExitSubLog
    END IF

    !! Cleanup
    CALL DestructMatrix(vector)
    CALL DestructMatrix(vector2)
    CALL DestructMatrixMemoryPool(pool)
    CALL DestructSolverParameters(param)
  END SUBROUTINE PowerBounds
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE EigenBoundsModule
