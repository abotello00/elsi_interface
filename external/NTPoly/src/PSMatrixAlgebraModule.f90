!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Performing Distributed Sparse Matrix Algebra Operations.
MODULE PSMatrixAlgebraModule
  USE DataTypesModule, ONLY : NTREAL, MPINTREAL, NTCOMPLEX, MPINTCOMPLEX
  USE GemmTasksModule
  USE MatrixReduceModule, ONLY : ReduceHelper_t, ReduceAndComposeMatrixSizes, &
       & ReduceAndComposeMatrixData, ReduceAndComposeMatrixCleanup, &
       & ReduceAndSumMatrixSizes, ReduceAndSumMatrixData, &
       & ReduceAndSumMatrixCleanup, TestReduceSizeRequest, &
       & TestReduceInnerRequest, TestReduceDataRequest
  USE PMatrixMemoryPoolModule, ONLY : MatrixMemoryPool_p, &
       & CheckMemoryPoolValidity, DestructMatrixMemoryPool, &
       & ConstructMatrixMemoryPool
  USE PSMatrixModule, ONLY : Matrix_ps, ConstructEmptyMatrix, CopyMatrix, &
       & DestructMatrix, ConvertMatrixToComplex, ConjugateMatrix, &
       & MergeMatrixLocalBlocks, IsIdentity, TransposeMatrix, &
       & SplitMatrixToLocalBlocks
  USE SMatrixAlgebraModule, ONLY : MatrixMultiply, MatrixGrandSum, &
       & PairwiseMultiplyMatrix, IncrementMatrix, ScaleMatrix, &
       & MatrixColumnNorm, MatrixDiagonalScale
  USE SMatrixModule, ONLY : Matrix_lsr, Matrix_lsc, DestructMatrix, CopyMatrix,&
       & TransposeMatrix, ComposeMatrixColumns, MatrixToTripletList
  USE TripletListModule, ONLY : TripletList_r, TripletList_c, &
       & ConstructTripletList, AppendToTripletList, DestructTripletList, &
       & GetTripletAt
  USE TripletModule, ONLY : Triplet_r, Triplet_c
  USE NTMPIModule
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: MatrixSigma
  PUBLIC :: MatrixMultiply
  PUBLIC :: MatrixGrandSum
  PUBLIC :: PairwiseMultiplyMatrix
  PUBLIC :: MatrixNorm
  PUBLIC :: DotMatrix
  PUBLIC :: IncrementMatrix
  PUBLIC :: ScaleMatrix
  PUBLIC :: MatrixTrace
  PUBLIC :: SimilarityTransform
  PUBLIC :: MeasureAsymmetry
  PUBLIC :: SymmetrizeMatrix
  PUBLIC :: MatrixDiagonalScale
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE MatrixSigma
     MODULE PROCEDURE MatrixSigma_ps
  END INTERFACE MatrixSigma
  INTERFACE MatrixMultiply
     MODULE PROCEDURE MatrixMultiply_ps
  END INTERFACE MatrixMultiply
  INTERFACE MatrixGrandSum
     MODULE PROCEDURE MatrixGrandSum_psr
     MODULE PROCEDURE MatrixGrandSum_psc
  END INTERFACE MatrixGrandSum
  INTERFACE PairwiseMultiplyMatrix
     MODULE PROCEDURE PairwiseMultiplyMatrix_ps
  END INTERFACE PairwiseMultiplyMatrix
  INTERFACE MatrixNorm
     MODULE PROCEDURE MatrixNorm_ps
  END INTERFACE MatrixNorm
  INTERFACE DotMatrix
     MODULE PROCEDURE DotMatrix_psr
     MODULE PROCEDURE DotMatrix_psc
  END INTERFACE DotMatrix
  INTERFACE IncrementMatrix
     MODULE PROCEDURE IncrementMatrix_ps
  END INTERFACE IncrementMatrix
  INTERFACE ScaleMatrix
     MODULE PROCEDURE ScaleMatrix_psr
     MODULE PROCEDURE ScaleMatrix_psc
  END INTERFACE ScaleMatrix
  INTERFACE MatrixDiagonalScale
     MODULE PROCEDURE MatrixDiagonalScale_psr
     MODULE PROCEDURE MatrixDiagonalScale_psc
  END INTERFACE MatrixDiagonalScale
  INTERFACE MatrixTrace
     MODULE PROCEDURE MatrixTrace_psr
  END INTERFACE MatrixTrace
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute sigma for the inversion method.
  !> See \cite ozaki2001efficient for details.
  SUBROUTINE MatrixSigma_ps(this, sigma_value)
    !> The matrix to compute the sigma value of.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> Sigma
    REAL(NTREAL), INTENT(OUT) :: sigma_value
    !! Local Data
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: column_sigma_contribution
    !! Counters/Temporary
    INTEGER :: II, JJ
    TYPE(Matrix_lsr) :: merged_local_data_r
    TYPE(Matrix_lsc) :: merged_local_data_c
    INTEGER :: ierr

    IF (this%is_complex) THEN


       !! Merge all the local data
       CALL MergeMatrixLocalBlocks(this, merged_local_data_c)

       ALLOCATE(column_sigma_contribution(merged_local_data_c%columns))
       column_sigma_contribution = 0
       DO II = 1, merged_local_data_c%columns
          DO JJ = merged_local_data_c%outer_index(II), merged_local_data_c%outer_index(II + 1) - 1
             column_sigma_contribution(II) = column_sigma_contribution(II) + &
                  & ABS(merged_local_data_c%values(JJ + 1))
          END DO
       END DO
       CALL MPI_Allreduce(MPI_IN_PLACE, column_sigma_contribution, &
            & merged_local_data_c%columns, MPINTREAL, MPI_SUM, &
            & this%process_grid%column_comm, ierr)
       CALL MPI_Allreduce(MAXVAL(column_sigma_contribution), sigma_value, 1, &
            & MPINTREAL, MPI_MAX, this%process_grid%row_comm, ierr)
       sigma_value = 1.0_NTREAL / (sigma_value**2)

       DEALLOCATE(column_sigma_contribution)
       CALL DestructMatrix(merged_local_data_c)

    ELSE


       !! Merge all the local data
       CALL MergeMatrixLocalBlocks(this, merged_local_data_r)

       ALLOCATE(column_sigma_contribution(merged_local_data_r%columns))
       column_sigma_contribution = 0
       DO II = 1, merged_local_data_r%columns
          DO JJ = merged_local_data_r%outer_index(II), merged_local_data_r%outer_index(II + 1) - 1
             column_sigma_contribution(II) = column_sigma_contribution(II) + &
                  & ABS(merged_local_data_r%values(JJ + 1))
          END DO
       END DO
       CALL MPI_Allreduce(MPI_IN_PLACE, column_sigma_contribution, &
            & merged_local_data_r%columns, MPINTREAL, MPI_SUM, &
            & this%process_grid%column_comm, ierr)
       CALL MPI_Allreduce(MAXVAL(column_sigma_contribution), sigma_value, 1, &
            & MPINTREAL, MPI_MAX, this%process_grid%row_comm, ierr)
       sigma_value = 1.0_NTREAL / (sigma_value**2)

       DEALLOCATE(column_sigma_contribution)
       CALL DestructMatrix(merged_local_data_r)

    ENDIF
  END SUBROUTINE MatrixSigma_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Multiply two matrices together, and add to the third.
  !> C := alpha*matA*matB+ beta*matC
  SUBROUTINE MatrixMultiply_ps(matA, matB, matC, alpha_in, beta_in, &
       & threshold_in, memory_pool_in)
    !> Matrix A.
    TYPE(Matrix_ps), INTENT(IN)        :: matA
    !> Matrix B.
    TYPE(Matrix_ps), INTENT(IN)        :: matB
    !> matC = alpha*matA*matB + beta*matC
    TYPE(Matrix_ps), INTENT(INOUT)     :: matC
    !> Scales the multiplication
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: alpha_in
    !> Scales matrix we sum on to.
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: beta_in
    !> For flushing values to zero.
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: threshold_in
    !> A memory pool for the calculation.
    TYPE(MatrixMemoryPool_p), OPTIONAL, INTENT(INOUT) :: memory_pool_in
    !! Local Versions of Optional Parameter
    TYPE(Matrix_ps) :: matAConverted
    TYPE(Matrix_ps) :: matBConverted
    REAL(NTREAL) :: alpha
    REAL(NTREAL) :: beta
    REAL(NTREAL) :: threshold
    TYPE(MatrixMemoryPool_p) :: memory_pool

    !! Handle the optional parameters
    IF (.NOT. PRESENT(alpha_in)) THEN
       alpha = 1.0_NTREAL
    ELSE
       alpha = alpha_in
    END IF
    IF (.NOT. PRESENT(beta_in)) THEN
       beta = 0.0_NTREAL
    ELSE
       beta = beta_in
    END IF
    IF (.NOT. PRESENT(threshold_in)) THEN
       threshold = 0.0_NTREAL
    ELSE
       threshold = threshold_in
    END IF

    !! Setup Memory Pool
    IF (PRESENT(memory_pool_in)) THEN
       IF (matA%is_complex) THEN
          IF (.NOT. CheckMemoryPoolValidity(memory_pool_in, matA)) THEN
             CALL DestructMatrixMemoryPool(memory_pool_in)
             CALL ConstructMatrixMemoryPool(memory_pool_in, matA)
          END IF
       ELSE
          IF (.NOT. CheckMemoryPoolValidity(memory_pool_in, matB)) THEN
             CALL DestructMatrixMemoryPool(memory_pool_in)
             CALL ConstructMatrixMemoryPool(memory_pool_in, matB)
          END IF
       END IF
    ELSE
       IF (matA%is_complex) THEN
          CALL ConstructMatrixMemoryPool(memory_pool, matA)
       ELSE
          CALL ConstructMatrixMemoryPool(memory_pool, matB)
       END IF
    END IF

    !! Perform Upcasting
    IF (matB%is_complex .AND. .NOT. matA%is_complex) THEN
       CALL ConvertMatrixToComplex(matA, matAConverted)
       IF (PRESENT(memory_pool_in)) THEN
          CALL MatrixMultiply_psc(matAConverted, matB, matC, alpha, beta, &
               & threshold, memory_pool_in)
       ELSE
          CALL MatrixMultiply_psc(matAConverted, matB, matC, alpha, beta, &
               & threshold, memory_pool)
       END IF
    ELSE IF (matA%is_complex .AND. .NOT. matB%is_complex) THEN
       CALL ConvertMatrixToComplex(matB, matBConverted)
       IF (PRESENT(memory_pool_in)) THEN
          CALL MatrixMultiply_psc(matA, matBConverted, matC, alpha, beta, &
               & threshold, memory_pool_in)
       ELSE
          CALL MatrixMultiply_psc(matA, matBConverted, matC, alpha, beta, &
               & threshold, memory_pool)
       END IF
    ELSE IF (matA%is_complex .AND. matB%is_complex) THEN
       IF (PRESENT(memory_pool_in)) THEN
          CALL MatrixMultiply_psc(matA, matB, matC, alpha, beta, &
               & threshold, memory_pool_in)
       ELSE
          CALL MatrixMultiply_psc(matA, matB, matC, alpha, beta, &
               & threshold, memory_pool)
       END IF
    ELSE
       IF (PRESENT(memory_pool_in)) THEN
          CALL MatrixMultiply_psr(matA, matB, matC, alpha, beta, &
               & threshold, memory_pool_in)
       ELSE
          CALL MatrixMultiply_psr(matA, matB, matC, alpha, beta, &
               & threshold, memory_pool)
       END IF
    END IF

    CALL DestructMatrixMemoryPool(memory_pool)
    CALL DestructMatrix(matAConverted)
    CALL DestructMatrix(matBConverted)

  END SUBROUTINE MatrixMultiply_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> The actual implementation of matrix multiply is here. Takes the
  !> same parameters as the standard multiply, but nothing is optional.
  SUBROUTINE MatrixMultiply_psr(matA, matB, matC, alpha, beta, &
       & threshold, memory_pool)
    !! Parameters
    TYPE(Matrix_ps), INTENT(IN)    :: matA
    TYPE(Matrix_ps), INTENT(IN)    :: matB
    TYPE(Matrix_ps), INTENT(INOUT) :: matC
    REAL(NTREAL), INTENT(IN) :: alpha
    REAL(NTREAL), INTENT(IN) :: beta
    REAL(NTREAL), INTENT(IN) :: threshold
    TYPE(MatrixMemoryPool_p), INTENT(INOUT) :: memory_pool
    !! Temporary Matrices
    TYPE(Matrix_lsr), DIMENSION(:,:), ALLOCATABLE :: AdjacentABlocks
    TYPE(Matrix_lsr), DIMENSION(:), ALLOCATABLE :: LocalRowContribution
    TYPE(Matrix_lsr), DIMENSION(:), ALLOCATABLE :: GatheredRowContribution
    TYPE(Matrix_lsr), DIMENSION(:), ALLOCATABLE :: GatheredRowContributionT
    TYPE(Matrix_lsr), DIMENSION(:,:), ALLOCATABLE :: TransposedBBlocks
    TYPE(Matrix_lsr), DIMENSION(:), ALLOCATABLE :: LocalColumnContribution
    TYPE(Matrix_lsr), DIMENSION(:), ALLOCATABLE :: GatheredColumnContribution
    TYPE(Matrix_lsr), DIMENSION(:,:), ALLOCATABLE :: SliceContribution




    !! Communication Helpers
    TYPE(ReduceHelper_t), DIMENSION(:), ALLOCATABLE :: row_helper
    TYPE(ReduceHelper_t), DIMENSION(:), ALLOCATABLE :: column_helper
    TYPE(ReduceHelper_t), DIMENSION(:, :), ALLOCATABLE :: slice_helper
    !! For Iterating Over Local Blocks
    INTEGER :: II, II2, II2_range
    INTEGER :: JJ, JJ2, JJ2_range
    INTEGER :: duplicate_start_column, duplicate_offset_column
    INTEGER :: duplicate_start_row, duplicate_offset_row
    REAL(NTREAL) :: working_threshold
    !! Scheduling the A work
    INTEGER, DIMENSION(:), ALLOCATABLE :: ATasks
    INTEGER :: ATasks_completed
    !! Scheduling the B work
    INTEGER, DIMENSION(:), ALLOCATABLE :: BTasks
    INTEGER :: BTasks_completed
    !! Scheduling the AB work
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: ABTasks
    INTEGER :: ABTasks_completed
    !! Temporary AB matrix for scaling.
    TYPE(Matrix_ps) :: matAB

    !! The threshold needs to be smaller if we are doing a sliced version
    !! because you might flush a value that would be kept in the summed version.
    IF (matA%process_grid%num_process_slices .GT. 1) THEN
       working_threshold = threshold / (matA%process_grid%num_process_slices*1000)
    ELSE
       working_threshold = threshold
    END IF

    !! Construct The Temporary Matrices
    CALL ConstructEmptyMatrix(matAB, matA)

    ALLOCATE(AdjacentABlocks(matAB%process_grid%number_of_blocks_rows, &
         & matAB%process_grid%number_of_blocks_columns / &
         & matAB%process_grid%num_process_slices))
    ALLOCATE(LocalRowContribution(matAB%process_grid%number_of_blocks_rows))
    ALLOCATE(GatheredRowContribution(matAB%process_grid%number_of_blocks_rows))
    ALLOCATE(GatheredRowContributionT(matAB%process_grid%number_of_blocks_rows))

    ALLOCATE(TransposedBBlocks(matAB%process_grid%number_of_blocks_rows / &
         & matAB%process_grid%num_process_slices, &
         & matAB%process_grid%number_of_blocks_columns))
    ALLOCATE(LocalColumnContribution(&
         & matAB%process_grid%number_of_blocks_columns))
    ALLOCATE(GatheredColumnContribution(&
         & matAB%process_grid%number_of_blocks_columns))
    ALLOCATE(SliceContribution(matAB%process_grid%number_of_blocks_rows, &
         & matAB%process_grid%number_of_blocks_columns))

    !! Helpers
    ALLOCATE(row_helper(matAB%process_grid%number_of_blocks_rows))
    ALLOCATE(column_helper(matAB%process_grid%number_of_blocks_columns))
    ALLOCATE(slice_helper(matAB%process_grid%number_of_blocks_rows, &
         & matAB%process_grid%number_of_blocks_columns))

    !! Construct the task queues
    ALLOCATE(ATasks(matAB%process_grid%number_of_blocks_rows))
    DO II = 1, matAB%process_grid%number_of_blocks_rows
       ATasks(II) = LocalGatherA
    END DO
    ALLOCATE(BTasks(matAB%process_grid%number_of_blocks_columns))
    DO JJ = 1, matAB%process_grid%number_of_blocks_columns
       BTasks(JJ) = LocalGatherB
    END DO
    ALLOCATE(ABTasks(matAB%process_grid%number_of_blocks_rows, &
         & matAB%process_grid%number_of_blocks_columns))
    DO JJ = 1, matAB%process_grid%number_of_blocks_columns
       DO II = 1, matAB%process_grid%number_of_blocks_rows
          ABTasks(II,JJ) = AwaitingAB
       END DO
    END DO

    !! Setup A Tasks
    duplicate_start_column = matAB%process_grid%my_slice + 1
    duplicate_offset_column = matAB%process_grid%num_process_slices

    !! Setup B Tasks
    duplicate_start_row = matAB%process_grid%my_slice + 1
    duplicate_offset_row = matAB%process_grid%num_process_slices

    !! Run A Tasks
    ATasks_completed = 0
    BTasks_completed = 0
    ABTasks_completed = 0

    !$OMP PARALLEL
    !$OMP MASTER
    DO WHILE (ATasks_completed .LT. SIZE(ATasks) .OR. &
         & BTasks_completed .LT. SIZE(BTasks) .OR. &
         & ABTasks_completed .LT. SIZE(ABTasks))
       DO II = 1, matAB%process_grid%number_of_blocks_rows
          SELECT CASE (ATasks(II))
          CASE(LocalGatherA)
             ATasks(II) = TaskRunningA
             !$OMP TASK DEFAULT(SHARED), PRIVATE(JJ2, JJ2_range), FIRSTPRIVATE(II)
             !! First Align The Data We Are Working With
             JJ2_range = matAB%process_grid%number_of_blocks_columns / &
                  & matAB%process_grid%num_process_slices
             DO JJ2 = 1, JJ2_range
                CALL CopyMatrix(matA%local_data_r(II, &
                     & duplicate_start_column + &
                     & duplicate_offset_column * (JJ2 - 1)),&
                     & AdjacentABlocks(II, JJ2))
             END DO
             !! Then Do A Local Gather and Cleanup
             CALL ComposeMatrixColumns(AdjacentABlocks(II, :), &
                  & LocalRowContribution(II))
             DO JJ2 = 1, JJ2_range
                CALL DestructMatrix(AdjacentABlocks(II, JJ2))
             END DO
             ATasks(II) = SendSizeA
             !$OMP END TASK
          CASE(SendSizeA)
             !! Then Start A Global Gather
             CALL ReduceAndComposeMatrixSizes(LocalRowContribution(II), &
                  & matAB%process_grid%blocked_row_comm(II), &
                  & GatheredRowContribution(II), row_helper(II))
             ATasks(II) = ComposeA
          CASE(ComposeA)
             IF (TestReduceSizeRequest(row_helper(II))) THEN
                CALL ReduceAndComposeMatrixData(LocalRowContribution(II), &
                     & matAB%process_grid%blocked_row_comm(II), &
                     & GatheredRowContribution(II), row_helper(II))
                ATasks(II) = WaitInnerA
             END IF
          CASE(WaitInnerA)
             IF (TestReduceInnerRequest(row_helper(II))) THEN
                ATasks(II) = WaitDataA
             END IF
          CASE(WaitDataA)
             IF (TestReduceDataRequest(row_helper(II))) THEN
                ATasks(II) = AdjustIndicesA
             END IF
          CASE(AdjustIndicesA)
             ATasks(II) = TaskRunningA
             !$OMP TASK DEFAULT(SHARED), FIRSTPRIVATE(II)
             CALL ReduceAndComposeMatrixCleanup(LocalRowContribution(II), &
                  & GatheredRowContribution(II), row_helper(II))
             CALL DestructMatrix(LocalRowContribution(II))
             CALL TransposeMatrix(GatheredRowContribution(II), &
                  & GatheredRowContributionT(II))
             CALL DestructMatrix(GatheredRowContribution(II))
             ATasks(II) = CleanupA
             !$OMP END TASK
          CASE(CleanupA)
             ATasks(II) = FinishedA
             ATasks_completed = ATasks_completed + 1
          END SELECT
       END DO
       !! B Tasks
       DO JJ = 1 , matAB%process_grid%number_of_blocks_columns
          SELECT CASE (BTasks(JJ))
          CASE(LocalGatherB)
             BTasks(JJ) = TaskRunningB
             !$OMP TASK DEFAULT(SHARED), PRIVATE(II2, II2_range), FIRSTPRIVATE(JJ)
             !! First Transpose The Data We Are Working With
             II2_range = matAB%process_grid%number_of_blocks_rows / &
                  & matAB%process_grid%num_process_slices
             DO II2 = 1, II2_range
                CALL TransposeMatrix(matB%local_data_r(duplicate_start_row + &
                     & duplicate_offset_row * (II2 - 1), JJ), &
                     & TransposedBBlocks(II2, JJ))
             END DO
             !! Then Do A Local Gather and Cleanup
             CALL ComposeMatrixColumns(TransposedBBlocks(:, JJ), &
                  & LocalColumnContribution(JJ))
             DO II2 = 1, II2_range
                CALL DestructMatrix(TransposedBBlocks(II2, JJ))
             END DO
             BTasks(JJ) = SendSizeB
             !$OMP END TASK
          CASE(SendSizeB)
             !! Then A Global Gather
             CALL ReduceAndComposeMatrixSizes(LocalColumnContribution(JJ), &
                  & matAB%process_grid%blocked_column_comm(JJ), &
                  & GatheredColumnContribution(JJ), column_helper(JJ))
             BTasks(JJ) = LocalComposeB
          CASE(LocalComposeB)
             IF (TestReduceSizeRequest(column_helper(JJ))) THEN
                CALL ReduceAndComposeMatrixData(LocalColumnContribution(JJ),&
                     & matAB%process_grid%blocked_column_comm(JJ), &
                     & GatheredColumnContribution(JJ), column_helper(JJ))
                BTasks(JJ) = WaitInnerB
             END IF
          CASE(WaitInnerB)
             IF (TestReduceInnerRequest(column_helper(JJ))) THEN
                BTasks(JJ) = WaitDataB
             END IF
          CASE(WaitDataB)
             IF (TestReduceDataRequest(column_helper(JJ))) THEN
                BTasks(JJ) = AdjustIndicesB
             END IF
          CASE(AdjustIndicesB)
             BTasks(JJ) = TaskRunningB
             !$OMP TASK DEFAULT(SHARED), FIRSTPRIVATE(JJ)
             CALL ReduceAndComposeMatrixCleanup(LocalColumnContribution(JJ), &
                  & GatheredColumnContribution(JJ), column_helper(JJ))
             CALL DestructMatrix(LocalColumnContribution(JJ))
             BTasks(JJ) = CleanupB
             !$OMP END TASK
          CASE(CleanupB)
             BTasks(JJ) = FinishedB
             BTasks_completed = BTasks_completed + 1
          END SELECT
       END DO
       !! AB Tasks
       DO II = 1 , matAB%process_grid%number_of_blocks_rows
          DO JJ = 1, matAB%process_grid%number_of_blocks_columns
             SELECT CASE(ABTasks(II, JJ))
             CASE (AwaitingAB)
                IF (ATasks(II) .EQ. FinishedA .AND. &
                     & BTasks(JJ) .EQ. FinishedB) THEN
                   ABTasks(II, JJ) = GemmAB
                END IF
             CASE (GemmAB)
                ABTasks(II, JJ) = TaskRunningAB
                !$OMP TASK DEFAULT(shared), FIRSTPRIVATE(II, JJ)
                CALL MatrixMultiply(GatheredRowContributionT(II), &
                     & GatheredColumnContribution(JJ), &
                     & SliceContribution(II, JJ), &
                     & IsATransposed_in = .TRUE., IsBTransposed_in = .TRUE., &
                     & alpha_in = alpha, threshold_in = working_threshold, &
                     & blocked_memory_pool_in = memory_pool%grid_r(II, JJ))
                !! We can exit early if there is only one process slice
                IF (matAB%process_grid%num_process_slices .EQ. 1) THEN
                   ABTasks(II,JJ) = CleanupAB
                   CALL CopyMatrix(SliceContribution(II, JJ), matAB%local_data_r(II, JJ))
                   CALL DestructMatrix(SliceContribution(II, JJ))
                ELSE
                   ABTasks(II, JJ) = SendSizeAB
                END IF
                !$OMP END TASK
             CASE(SendSizeAB)
                CALL ReduceAndSumMatrixSizes(SliceContribution(II, JJ),&
                     & matAB%process_grid%blocked_between_slice_comm(II, JJ), &
                     & matAB%local_data_r(II, JJ), slice_helper(II, JJ))
                ABTasks(II, JJ) = GatherAndSumAB
             CASE (GatherAndSumAB)
                IF (TestReduceSizeRequest(slice_helper(II, JJ))) THEN
                   CALL ReduceAndSumMatrixData(SliceContribution(II, JJ), &
                        & matAB%process_grid%blocked_between_slice_comm(II, JJ), &
                        & matAB%local_data_r(II, JJ), slice_helper(II, JJ))
                   ABTasks(II, JJ) = WaitInnerAB
                END IF
             CASE (WaitInnerAB)
                IF (TestReduceInnerRequest(slice_helper(II, JJ))) THEN
                   ABTasks(II, JJ) = WaitDataAB
                END IF
             CASE (WaitDataAB)
                IF (TestReduceDataRequest(slice_helper(II, JJ))) THEN
                   ABTasks(II, JJ) = LocalSumAB
                END IF
             CASE(LocalSumAB)
                ABTasks(II, JJ) = TaskRunningAB
                !$OMP TASK DEFAULT(SHARED), FIRSTPRIVATE(II, JJ)
                CALL ReduceAndSumMatrixCleanup(SliceContribution(II, JJ), &
                     & matAB%local_data_r(II, JJ), threshold, slice_helper(II, JJ))
                CALL DestructMatrix(SliceContribution(II, JJ))
                ABTasks(II, JJ) = CleanupAB
                !$OMP END TASK
             CASE(CleanupAB)
                ABTasks(II, JJ) = FinishedAB
                ABTasks_completed = ABTasks_completed + 1
             END SELECT
          END DO
       END DO
       !! Prevent deadlock in the case where the number of tasks is capped.
       IF (matA%process_grid%omp_max_threads .EQ. 1) THEN
          !$OMP taskwait
       END IF
    END DO
    !$OMP END MASTER
    !$OMP END PARALLEL

    !! Cleanup
    DEALLOCATE(row_helper)
    DEALLOCATE(column_helper)
    DEALLOCATE(slice_helper)
    DEALLOCATE(ATasks)
    DEALLOCATE(BTasks)
    DEALLOCATE(ABTasks)

    !! Deallocate Buffers From A
    DO II = 1, matAB%process_grid%number_of_blocks_rows
       DO JJ2 = 1, matAB%process_grid%number_of_blocks_columns / &
            & matAB%process_grid%num_process_slices
          CALL DestructMatrix(AdjacentABlocks(II, JJ2))
       END DO
       CALL DestructMatrix(LocalRowContribution(II))
       CALL DestructMatrix(GatheredRowContribution(II))
    END DO
    DEALLOCATE(AdjacentABlocks)
    DEALLOCATE(LocalRowContribution)
    DEALLOCATE(GatheredRowContribution)
    !! Deallocate Buffers From B
    DO JJ = 1, matAB%process_grid%number_of_blocks_columns
       DO II2 = 1, matAB%process_grid%number_of_blocks_rows / &
            & matAB%process_grid%num_process_slices
          CALL DestructMatrix(TransposedBBlocks(II2, JJ))
       END DO
       CALL DestructMatrix(LocalColumnContribution(JJ))
    END DO
    DEALLOCATE(TransposedBBlocks)
    DEALLOCATE(LocalColumnContribution)
    !! Deallocate Buffers From Multiplying The Block
    DO II = 1, matAB%process_grid%number_of_blocks_rows
       CALL DestructMatrix(GatheredRowContributionT(II))
    END DO
    DO JJ = 1, matAB%process_grid%number_of_blocks_columns
       CALL DestructMatrix(GatheredColumnContribution(JJ))
    END DO
    DEALLOCATE(GatheredRowContributionT)
    DEALLOCATE(GatheredColumnContribution)
    !! Deallocate Buffers From Sum
    DO JJ = 1, matAB%process_grid%number_of_blocks_columns
       DO II = 1, matAB%process_grid%number_of_blocks_rows
          CALL DestructMatrix(SliceContribution(II, JJ))
       END DO
    END DO
    DEALLOCATE(SliceContribution)

    !! Copy to output matrix.
    IF (ABS(beta) .LT. TINY(beta)) THEN
       CALL CopyMatrix(matAB, matC)
    ELSE
       CALL ScaleMatrix(MatC, beta)
       CALL IncrementMatrix(MatAB, MatC)
    END IF
    CALL DestructMatrix(matAB)


  END SUBROUTINE MatrixMultiply_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> The actual implementation of matrix multiply is here. Takes the
  !> same parameters as the standard multiply, but nothing is optional.
  SUBROUTINE MatrixMultiply_psc(matA, matB, matC, alpha, beta, &
       & threshold, memory_pool)
    !! Parameters
    TYPE(Matrix_ps), INTENT(IN)    :: matA
    TYPE(Matrix_ps), INTENT(IN)    :: matB
    TYPE(Matrix_ps), INTENT(INOUT) :: matC
    REAL(NTREAL), INTENT(IN) :: alpha
    REAL(NTREAL), INTENT(IN) :: beta
    REAL(NTREAL), INTENT(IN) :: threshold
    TYPE(MatrixMemoryPool_p), INTENT(INOUT) :: memory_pool
    !! Temporary Matrices
    TYPE(Matrix_lsc), DIMENSION(:,:), ALLOCATABLE :: AdjacentABlocks
    TYPE(Matrix_lsc), DIMENSION(:), ALLOCATABLE :: LocalRowContribution
    TYPE(Matrix_lsc), DIMENSION(:), ALLOCATABLE :: GatheredRowContribution
    TYPE(Matrix_lsc), DIMENSION(:), ALLOCATABLE :: GatheredRowContributionT
    TYPE(Matrix_lsc), DIMENSION(:,:), ALLOCATABLE :: TransposedBBlocks
    TYPE(Matrix_lsc), DIMENSION(:), ALLOCATABLE :: LocalColumnContribution
    TYPE(Matrix_lsc), DIMENSION(:), ALLOCATABLE :: GatheredColumnContribution
    TYPE(Matrix_lsc), DIMENSION(:,:), ALLOCATABLE :: SliceContribution




    !! Communication Helpers
    TYPE(ReduceHelper_t), DIMENSION(:), ALLOCATABLE :: row_helper
    TYPE(ReduceHelper_t), DIMENSION(:), ALLOCATABLE :: column_helper
    TYPE(ReduceHelper_t), DIMENSION(:, :), ALLOCATABLE :: slice_helper
    !! For Iterating Over Local Blocks
    INTEGER :: II, II2, II2_range
    INTEGER :: JJ, JJ2, JJ2_range
    INTEGER :: duplicate_start_column, duplicate_offset_column
    INTEGER :: duplicate_start_row, duplicate_offset_row
    REAL(NTREAL) :: working_threshold
    !! Scheduling the A work
    INTEGER, DIMENSION(:), ALLOCATABLE :: ATasks
    INTEGER :: ATasks_completed
    !! Scheduling the B work
    INTEGER, DIMENSION(:), ALLOCATABLE :: BTasks
    INTEGER :: BTasks_completed
    !! Scheduling the AB work
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: ABTasks
    INTEGER :: ABTasks_completed
    !! Temporary AB matrix for scaling.
    TYPE(Matrix_ps) :: matAB

    !! The threshold needs to be smaller if we are doing a sliced version
    !! because you might flush a value that would be kept in the summed version.
    IF (matA%process_grid%num_process_slices .GT. 1) THEN
       working_threshold = threshold / (matA%process_grid%num_process_slices*1000)
    ELSE
       working_threshold = threshold
    END IF

    !! Construct The Temporary Matrices
    CALL ConstructEmptyMatrix(matAB, matA)

    ALLOCATE(AdjacentABlocks(matAB%process_grid%number_of_blocks_rows, &
         & matAB%process_grid%number_of_blocks_columns / &
         & matAB%process_grid%num_process_slices))
    ALLOCATE(LocalRowContribution(matAB%process_grid%number_of_blocks_rows))
    ALLOCATE(GatheredRowContribution(matAB%process_grid%number_of_blocks_rows))
    ALLOCATE(GatheredRowContributionT(matAB%process_grid%number_of_blocks_rows))

    ALLOCATE(TransposedBBlocks(matAB%process_grid%number_of_blocks_rows / &
         & matAB%process_grid%num_process_slices, &
         & matAB%process_grid%number_of_blocks_columns))
    ALLOCATE(LocalColumnContribution(&
         & matAB%process_grid%number_of_blocks_columns))
    ALLOCATE(GatheredColumnContribution(&
         & matAB%process_grid%number_of_blocks_columns))
    ALLOCATE(SliceContribution(matAB%process_grid%number_of_blocks_rows, &
         & matAB%process_grid%number_of_blocks_columns))

    !! Helpers
    ALLOCATE(row_helper(matAB%process_grid%number_of_blocks_rows))
    ALLOCATE(column_helper(matAB%process_grid%number_of_blocks_columns))
    ALLOCATE(slice_helper(matAB%process_grid%number_of_blocks_rows, &
         & matAB%process_grid%number_of_blocks_columns))

    !! Construct the task queues
    ALLOCATE(ATasks(matAB%process_grid%number_of_blocks_rows))
    DO II = 1, matAB%process_grid%number_of_blocks_rows
       ATasks(II) = LocalGatherA
    END DO
    ALLOCATE(BTasks(matAB%process_grid%number_of_blocks_columns))
    DO JJ = 1, matAB%process_grid%number_of_blocks_columns
       BTasks(JJ) = LocalGatherB
    END DO
    ALLOCATE(ABTasks(matAB%process_grid%number_of_blocks_rows, &
         & matAB%process_grid%number_of_blocks_columns))
    DO JJ = 1, matAB%process_grid%number_of_blocks_columns
       DO II = 1, matAB%process_grid%number_of_blocks_rows
          ABTasks(II,JJ) = AwaitingAB
       END DO
    END DO

    !! Setup A Tasks
    duplicate_start_column = matAB%process_grid%my_slice + 1
    duplicate_offset_column = matAB%process_grid%num_process_slices

    !! Setup B Tasks
    duplicate_start_row = matAB%process_grid%my_slice + 1
    duplicate_offset_row = matAB%process_grid%num_process_slices

    !! Run A Tasks
    ATasks_completed = 0
    BTasks_completed = 0
    ABTasks_completed = 0

    !$OMP PARALLEL
    !$OMP MASTER
    DO WHILE (ATasks_completed .LT. SIZE(ATasks) .OR. &
         & BTasks_completed .LT. SIZE(BTasks) .OR. &
         & ABTasks_completed .LT. SIZE(ABTasks))
       DO II = 1, matAB%process_grid%number_of_blocks_rows
          SELECT CASE (ATasks(II))
          CASE(LocalGatherA)
             ATasks(II) = TaskRunningA
             !$OMP TASK DEFAULT(SHARED), PRIVATE(JJ2, JJ2_range), FIRSTPRIVATE(II)
             !! First Align The Data We Are Working With
             JJ2_range = matAB%process_grid%number_of_blocks_columns / &
                  & matAB%process_grid%num_process_slices
             DO JJ2 = 1, JJ2_range
                CALL CopyMatrix(matA%local_data_c(II, &
                     & duplicate_start_column + &
                     & duplicate_offset_column * (JJ2 - 1)),&
                     & AdjacentABlocks(II, JJ2))
             END DO
             !! Then Do A Local Gather and Cleanup
             CALL ComposeMatrixColumns(AdjacentABlocks(II, :), &
                  & LocalRowContribution(II))
             DO JJ2 = 1, JJ2_range
                CALL DestructMatrix(AdjacentABlocks(II, JJ2))
             END DO
             ATasks(II) = SendSizeA
             !$OMP END TASK
          CASE(SendSizeA)
             !! Then Start A Global Gather
             CALL ReduceAndComposeMatrixSizes(LocalRowContribution(II), &
                  & matAB%process_grid%blocked_row_comm(II), &
                  & GatheredRowContribution(II), row_helper(II))
             ATasks(II) = ComposeA
          CASE(ComposeA)
             IF (TestReduceSizeRequest(row_helper(II))) THEN
                CALL ReduceAndComposeMatrixData(LocalRowContribution(II), &
                     & matAB%process_grid%blocked_row_comm(II), &
                     & GatheredRowContribution(II), row_helper(II))
                ATasks(II) = WaitInnerA
             END IF
          CASE(WaitInnerA)
             IF (TestReduceInnerRequest(row_helper(II))) THEN
                ATasks(II) = WaitDataA
             END IF
          CASE(WaitDataA)
             IF (TestReduceDataRequest(row_helper(II))) THEN
                ATasks(II) = AdjustIndicesA
             END IF
          CASE(AdjustIndicesA)
             ATasks(II) = TaskRunningA
             !$OMP TASK DEFAULT(SHARED), FIRSTPRIVATE(II)
             CALL ReduceAndComposeMatrixCleanup(LocalRowContribution(II), &
                  & GatheredRowContribution(II), row_helper(II))
             CALL DestructMatrix(LocalRowContribution(II))
             CALL TransposeMatrix(GatheredRowContribution(II), &
                  & GatheredRowContributionT(II))
             CALL DestructMatrix(GatheredRowContribution(II))
             ATasks(II) = CleanupA
             !$OMP END TASK
          CASE(CleanupA)
             ATasks(II) = FinishedA
             ATasks_completed = ATasks_completed + 1
          END SELECT
       END DO
       !! B Tasks
       DO JJ = 1 , matAB%process_grid%number_of_blocks_columns
          SELECT CASE (BTasks(JJ))
          CASE(LocalGatherB)
             BTasks(JJ) = TaskRunningB
             !$OMP TASK DEFAULT(SHARED), PRIVATE(II2, II2_range), FIRSTPRIVATE(JJ)
             !! First Transpose The Data We Are Working With
             II2_range = matAB%process_grid%number_of_blocks_rows / &
                  & matAB%process_grid%num_process_slices
             DO II2 = 1, II2_range
                CALL TransposeMatrix(matB%local_data_c(duplicate_start_row + &
                     & duplicate_offset_row * (II2 - 1), JJ), &
                     & TransposedBBlocks(II2, JJ))
             END DO
             !! Then Do A Local Gather and Cleanup
             CALL ComposeMatrixColumns(TransposedBBlocks(:, JJ), &
                  & LocalColumnContribution(JJ))
             DO II2 = 1, II2_range
                CALL DestructMatrix(TransposedBBlocks(II2, JJ))
             END DO
             BTasks(JJ) = SendSizeB
             !$OMP END TASK
          CASE(SendSizeB)
             !! Then A Global Gather
             CALL ReduceAndComposeMatrixSizes(LocalColumnContribution(JJ), &
                  & matAB%process_grid%blocked_column_comm(JJ), &
                  & GatheredColumnContribution(JJ), column_helper(JJ))
             BTasks(JJ) = LocalComposeB
          CASE(LocalComposeB)
             IF (TestReduceSizeRequest(column_helper(JJ))) THEN
                CALL ReduceAndComposeMatrixData(LocalColumnContribution(JJ),&
                     & matAB%process_grid%blocked_column_comm(JJ), &
                     & GatheredColumnContribution(JJ), column_helper(JJ))
                BTasks(JJ) = WaitInnerB
             END IF
          CASE(WaitInnerB)
             IF (TestReduceInnerRequest(column_helper(JJ))) THEN
                BTasks(JJ) = WaitDataB
             END IF
          CASE(WaitDataB)
             IF (TestReduceDataRequest(column_helper(JJ))) THEN
                BTasks(JJ) = AdjustIndicesB
             END IF
          CASE(AdjustIndicesB)
             BTasks(JJ) = TaskRunningB
             !$OMP TASK DEFAULT(SHARED), FIRSTPRIVATE(JJ)
             CALL ReduceAndComposeMatrixCleanup(LocalColumnContribution(JJ), &
                  & GatheredColumnContribution(JJ), column_helper(JJ))
             CALL DestructMatrix(LocalColumnContribution(JJ))
             BTasks(JJ) = CleanupB
             !$OMP END TASK
          CASE(CleanupB)
             BTasks(JJ) = FinishedB
             BTasks_completed = BTasks_completed + 1
          END SELECT
       END DO
       !! AB Tasks
       DO II = 1 , matAB%process_grid%number_of_blocks_rows
          DO JJ = 1, matAB%process_grid%number_of_blocks_columns
             SELECT CASE(ABTasks(II, JJ))
             CASE (AwaitingAB)
                IF (ATasks(II) .EQ. FinishedA .AND. &
                     & BTasks(JJ) .EQ. FinishedB) THEN
                   ABTasks(II, JJ) = GemmAB
                END IF
             CASE (GemmAB)
                ABTasks(II, JJ) = TaskRunningAB
                !$OMP TASK DEFAULT(shared), FIRSTPRIVATE(II, JJ)
                CALL MatrixMultiply(GatheredRowContributionT(II), &
                     & GatheredColumnContribution(JJ), &
                     & SliceContribution(II, JJ), &
                     & IsATransposed_in = .TRUE., IsBTransposed_in = .TRUE., &
                     & alpha_in = alpha, threshold_in = working_threshold, &
                     & blocked_memory_pool_in = memory_pool%grid_c(II, JJ))
                !! We can exit early if there is only one process slice
                IF (matAB%process_grid%num_process_slices .EQ. 1) THEN
                   ABTasks(II,JJ) = CleanupAB
                   CALL CopyMatrix(SliceContribution(II, JJ), matAB%local_data_c(II, JJ))
                   CALL DestructMatrix(SliceContribution(II, JJ))
                ELSE
                   ABTasks(II, JJ) = SendSizeAB
                END IF
                !$OMP END TASK
             CASE(SendSizeAB)
                CALL ReduceAndSumMatrixSizes(SliceContribution(II, JJ),&
                     & matAB%process_grid%blocked_between_slice_comm(II, JJ), &
                     & matAB%local_data_c(II, JJ), slice_helper(II, JJ))
                ABTasks(II, JJ) = GatherAndSumAB
             CASE (GatherAndSumAB)
                IF (TestReduceSizeRequest(slice_helper(II, JJ))) THEN
                   CALL ReduceAndSumMatrixData(SliceContribution(II, JJ), &
                        & matAB%process_grid%blocked_between_slice_comm(II, JJ), &
                        & matAB%local_data_c(II, JJ), slice_helper(II, JJ))
                   ABTasks(II, JJ) = WaitInnerAB
                END IF
             CASE (WaitInnerAB)
                IF (TestReduceInnerRequest(slice_helper(II, JJ))) THEN
                   ABTasks(II, JJ) = WaitDataAB
                END IF
             CASE (WaitDataAB)
                IF (TestReduceDataRequest(slice_helper(II, JJ))) THEN
                   ABTasks(II, JJ) = LocalSumAB
                END IF
             CASE(LocalSumAB)
                ABTasks(II, JJ) = TaskRunningAB
                !$OMP TASK DEFAULT(SHARED), FIRSTPRIVATE(II, JJ)
                CALL ReduceAndSumMatrixCleanup(SliceContribution(II, JJ), &
                     & matAB%local_data_c(II, JJ), threshold, slice_helper(II, JJ))
                CALL DestructMatrix(SliceContribution(II, JJ))
                ABTasks(II, JJ) = CleanupAB
                !$OMP END TASK
             CASE(CleanupAB)
                ABTasks(II, JJ) = FinishedAB
                ABTasks_completed = ABTasks_completed + 1
             END SELECT
          END DO
       END DO
       !! Prevent deadlock in the case where the number of tasks is capped.
       IF (matA%process_grid%omp_max_threads .EQ. 1) THEN
          !$OMP taskwait
       END IF
    END DO
    !$OMP END MASTER
    !$OMP END PARALLEL

    !! Cleanup
    DEALLOCATE(row_helper)
    DEALLOCATE(column_helper)
    DEALLOCATE(slice_helper)
    DEALLOCATE(ATasks)
    DEALLOCATE(BTasks)
    DEALLOCATE(ABTasks)

    !! Deallocate Buffers From A
    DO II = 1, matAB%process_grid%number_of_blocks_rows
       DO JJ2 = 1, matAB%process_grid%number_of_blocks_columns / &
            & matAB%process_grid%num_process_slices
          CALL DestructMatrix(AdjacentABlocks(II, JJ2))
       END DO
       CALL DestructMatrix(LocalRowContribution(II))
       CALL DestructMatrix(GatheredRowContribution(II))
    END DO
    DEALLOCATE(AdjacentABlocks)
    DEALLOCATE(LocalRowContribution)
    DEALLOCATE(GatheredRowContribution)
    !! Deallocate Buffers From B
    DO JJ = 1, matAB%process_grid%number_of_blocks_columns
       DO II2 = 1, matAB%process_grid%number_of_blocks_rows / &
            & matAB%process_grid%num_process_slices
          CALL DestructMatrix(TransposedBBlocks(II2, JJ))
       END DO
       CALL DestructMatrix(LocalColumnContribution(JJ))
    END DO
    DEALLOCATE(TransposedBBlocks)
    DEALLOCATE(LocalColumnContribution)
    !! Deallocate Buffers From Multiplying The Block
    DO II = 1, matAB%process_grid%number_of_blocks_rows
       CALL DestructMatrix(GatheredRowContributionT(II))
    END DO
    DO JJ = 1, matAB%process_grid%number_of_blocks_columns
       CALL DestructMatrix(GatheredColumnContribution(JJ))
    END DO
    DEALLOCATE(GatheredRowContributionT)
    DEALLOCATE(GatheredColumnContribution)
    !! Deallocate Buffers From Sum
    DO JJ = 1, matAB%process_grid%number_of_blocks_columns
       DO II = 1, matAB%process_grid%number_of_blocks_rows
          CALL DestructMatrix(SliceContribution(II, JJ))
       END DO
    END DO
    DEALLOCATE(SliceContribution)

    !! Copy to output matrix.
    IF (ABS(beta) .LT. TINY(beta)) THEN
       CALL CopyMatrix(matAB, matC)
    ELSE
       CALL ScaleMatrix(MatC, beta)
       CALL IncrementMatrix(MatAB, MatC)
    END IF
    CALL DestructMatrix(matAB)


  END SUBROUTINE MatrixMultiply_psc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sum up the elements in a matrix into a single value.
  SUBROUTINE MatrixGrandSum_psr(this, sum)
    !> The matrix to compute.
    TYPE(Matrix_ps), INTENT(IN)  :: this
    !> The sum of all elements.
    REAL(NTREAL), INTENT(OUT) :: sum
    !! Local Data
    INTEGER :: II, JJ
    REAL(NTREAL) :: temp_r
    COMPLEX(NTCOMPLEX) :: temp_c
    INTEGER :: ierr


    IF (this%is_complex) THEN



       sum = 0
       DO JJ = 1, this%process_grid%number_of_blocks_columns
          DO II = 1, this%process_grid%number_of_blocks_rows
             CALL MatrixGrandSum(this%local_data_c(II, JJ), temp_c)
             sum = sum + REAL(temp_c, KIND = NTREAL)
          END DO
       END DO

       !! Sum Among Process Slice
       CALL MPI_Allreduce(MPI_IN_PLACE, sum, 1, MPINTREAL, &
            & MPI_SUM, this%process_grid%within_slice_comm, ierr)


    ELSE



       sum = 0
       DO JJ = 1, this%process_grid%number_of_blocks_columns
          DO II = 1, this%process_grid%number_of_blocks_rows
             CALL MatrixGrandSum(this%local_data_r(II, JJ), temp_r)
             sum = sum + REAL(temp_r, KIND = NTREAL)
          END DO
       END DO

       !! Sum Among Process Slice
       CALL MPI_Allreduce(MPI_IN_PLACE, sum, 1, MPINTREAL, &
            & MPI_SUM, this%process_grid%within_slice_comm, ierr)


    END IF

  END SUBROUTINE MatrixGrandSum_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sum up the elements in a matrix into a single value.
  SUBROUTINE MatrixGrandSum_psc(this, sum)
    !> The matrix to compute.
    TYPE(Matrix_ps), INTENT(IN)  :: this
    !> The sum of all elements.
    COMPLEX(NTCOMPLEX), INTENT(OUT) :: sum
    !! Local Data
    INTEGER :: II, JJ
    REAL(NTREAL) :: temp_r
    COMPLEX(NTCOMPLEX) :: temp_c
    INTEGER :: ierr


    IF (this%is_complex) THEN



       sum = 0
       DO JJ = 1, this%process_grid%number_of_blocks_columns
          DO II = 1, this%process_grid%number_of_blocks_rows
             CALL MatrixGrandSum(this%local_data_c(II, JJ), temp_c)
             sum = sum + REAL(temp_c, KIND = NTREAL)
          END DO
       END DO

       !! Sum Among Process Slice
       CALL MPI_Allreduce(MPI_IN_PLACE, sum, 1, MPINTCOMPLEX, &
            & MPI_SUM, this%process_grid%within_slice_comm, ierr)


    ELSE



       sum = 0
       DO JJ = 1, this%process_grid%number_of_blocks_columns
          DO II = 1, this%process_grid%number_of_blocks_rows
             CALL MatrixGrandSum(this%local_data_r(II, JJ), temp_r)
             sum = sum + REAL(temp_r, KIND = NTREAL)
          END DO
       END DO

       !! Sum Among Process Slice
       CALL MPI_Allreduce(MPI_IN_PLACE, sum, 1, MPINTCOMPLEX, &
            & MPI_SUM, this%process_grid%within_slice_comm, ierr)


    END IF

  END SUBROUTINE MatrixGrandSum_psc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Elementwise multiplication. C_ij = A_ij * B_ij.
  !> Also known as a Hadamard product.
  RECURSIVE SUBROUTINE PairwiseMultiplyMatrix_ps(matA, matB, matC)
    !> Matrix A.
    TYPE(Matrix_ps), INTENT(IN)  :: matA
    !> Matrix B.
    TYPE(Matrix_ps), INTENT(IN)  :: matB
    !> matC = MatA mult MatB.
    TYPE(Matrix_ps), INTENT(INOUT)  :: matC
    !! Local Data
    TYPE(Matrix_ps) :: converted_matrix
    INTEGER :: II, JJ

    IF (matA%is_complex .AND. .NOT. matB%is_complex) THEN
       CALL ConvertMatrixToComplex(matB, converted_matrix)
       CALL PairwiseMultiplyMatrix(matA, converted_matrix, matC)
       CALL DestructMatrix(converted_matrix)
    ELSE IF (.NOT. matA%is_complex .AND. matB%is_complex) THEN
       CALL ConvertMatrixToComplex(matA, converted_matrix)
       CALL PairwiseMultiplyMatrix(converted_matrix, matB, matC)
       CALL DestructMatrix(converted_matrix)
    ELSE IF (matA%is_complex .AND. matB%is_complex) THEN


       CALL ConstructEmptyMatrix(matC, matA%actual_matrix_dimension, &
            & matA%process_grid, matA%is_complex)

       !$omp parallel
       !$omp do collapse(2)
       DO JJ = 1, matA%process_grid%number_of_blocks_columns
          DO II = 1, matA%process_grid%number_of_blocks_rows
             CALL PairwiseMultiplyMatrix(matA%local_data_c(II, JJ), matB%local_data_c(II, JJ), &
                  & matC%local_data_c(II, JJ))
          END DO
       END DO
       !$omp end do
       !$omp end parallel

    ELSE


       CALL ConstructEmptyMatrix(matC, matA%actual_matrix_dimension, &
            & matA%process_grid, matA%is_complex)

       !$omp parallel
       !$omp do collapse(2)
       DO JJ = 1, matA%process_grid%number_of_blocks_columns
          DO II = 1, matA%process_grid%number_of_blocks_rows
             CALL PairwiseMultiplyMatrix(matA%local_data_r(II, JJ), matB%local_data_r(II, JJ), &
                  & matC%local_data_r(II, JJ))
          END DO
       END DO
       !$omp end do
       !$omp end parallel

    END IF
  END SUBROUTINE PairwiseMultiplyMatrix_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the norm of a distributed sparse matrix along the rows.
  FUNCTION MatrixNorm_ps(this) RESULT(norm_value)
    !> The matrix to compute the norm of.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> The norm value of the full distributed sparse matrix.
    REAL(NTREAL) :: norm_value
    !! Local Data
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: local_norm
    TYPE(Matrix_lsr) :: merged_local_data_r
    TYPE(Matrix_lsc) :: merged_local_data_c
    INTEGER :: ierr

    IF (this%is_complex) THEN


       !! Merge all the local data
       CALL MergeMatrixLocalBlocks(this, merged_local_data_c)
       ALLOCATE(local_norm(merged_local_data_c%columns))

       !! Sum Along Columns
       CALL MatrixColumnNorm(merged_local_data_c, local_norm)
       CALL MPI_Allreduce(MPI_IN_PLACE, local_norm, SIZE(local_norm), &
            & MPINTREAL, MPI_SUM, this%process_grid%column_comm, ierr)

       !! Find Max Value Amonst Columns
       norm_value = MAXVAL(local_norm)
       CALL MPI_Allreduce(MPI_IN_PLACE, norm_value,1, MPINTREAL, MPI_MAX, &
            & this%process_grid%row_comm, ierr)

       CALL DestructMatrix(merged_local_data_c)
       DEALLOCATE(local_norm)

    ELSE


       !! Merge all the local data
       CALL MergeMatrixLocalBlocks(this, merged_local_data_r)
       ALLOCATE(local_norm(merged_local_data_r%columns))

       !! Sum Along Columns
       CALL MatrixColumnNorm(merged_local_data_r, local_norm)
       CALL MPI_Allreduce(MPI_IN_PLACE, local_norm, SIZE(local_norm), &
            & MPINTREAL, MPI_SUM, this%process_grid%column_comm, ierr)

       !! Find Max Value Amonst Columns
       norm_value = MAXVAL(local_norm)
       CALL MPI_Allreduce(MPI_IN_PLACE, norm_value,1, MPINTREAL, MPI_MAX, &
            & this%process_grid%row_comm, ierr)

       CALL DestructMatrix(merged_local_data_r)
       DEALLOCATE(local_norm)

    END IF
  END FUNCTION MatrixNorm_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> product = dot(Matrix A,Matrix B)
  !> Note that a dot product is the sum of elementwise multiplication, not
  !> traditional matrix multiplication.
  SUBROUTINE DotMatrix_psr(matA, matB, product)
    !> Matrix A.
    TYPE(Matrix_ps), INTENT(IN)  :: matA
    !> Matrix B.
    TYPE(Matrix_ps), INTENT(IN)  :: matB
    !> The dot product.
    REAL(NTREAL), INTENT(OUT) :: product


    !! Local Data
    TYPE(Matrix_ps) :: matAH
    TYPE(Matrix_ps) :: matC

    IF (matA%is_complex) THEN
       CALL CopyMatrix(matA, matAH)
       CALL ConjugateMatrix(matAH)
       CALL PairwiseMultiplyMatrix(matAH, matB, matC)
       CALL DestructMatrix(matAH)
    ELSE
       CALL PairwiseMultiplyMatrix(matA, matB, matC)
    END IF

    CALL MatrixGrandSum(matC, product)
    CALL DestructMatrix(matC)
  END SUBROUTINE DotMatrix_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> product = dot(Matrix A,Matrix B)
  !> Note that a dot product is the sum of elementwise multiplication, not
  !> traditional matrix multiplication.
  SUBROUTINE DotMatrix_psc(matA, matB, product)
    !> Matrix A.
    TYPE(Matrix_ps), INTENT(IN)  :: matA
    !> Matrix B.
    TYPE(Matrix_ps), INTENT(IN)  :: matB
    !> The dot product.
    COMPLEX(NTCOMPLEX), INTENT(OUT) :: product


    !! Local Data
    TYPE(Matrix_ps) :: matAH
    TYPE(Matrix_ps) :: matC

    IF (matA%is_complex) THEN
       CALL CopyMatrix(matA, matAH)
       CALL ConjugateMatrix(matAH)
       CALL PairwiseMultiplyMatrix(matAH, matB, matC)
       CALL DestructMatrix(matAH)
    ELSE
       CALL PairwiseMultiplyMatrix(matA, matB, matC)
    END IF

    CALL MatrixGrandSum(matC, product)
    CALL DestructMatrix(matC)
  END SUBROUTINE DotMatrix_psc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Matrix B = alpha*Matrix A + Matrix B (AXPY)
  !> This will utilize the sparse vector increment routine.
  RECURSIVE SUBROUTINE IncrementMatrix_ps(matA, matB, alpha_in, threshold_in)
    !> Matrix A.
    TYPE(Matrix_ps), INTENT(IN)  :: matA
    !> Matrix B.
    TYPE(Matrix_ps), INTENT(INOUT)  :: matB
    !> Multiplier (default = 1.0).
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: alpha_in
    !> For flushing values to zero (default = 0).
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: threshold_in
    !! Local Data
    TYPE(Matrix_ps) :: converted_matrix
    REAL(NTREAL) :: alpha
    REAL(NTREAL) :: threshold
    INTEGER :: II, JJ

    !! Optional Parameters
    IF (.NOT. PRESENT(alpha_in)) THEN
       alpha = 1.0_NTREAL
    ELSE
       alpha = alpha_in
    END IF
    IF (.NOT. PRESENT(threshold_in)) THEN
       threshold = 0.0_NTREAL
    ELSE
       threshold = threshold_in
    END IF

    IF (matA%is_complex .AND. .NOT. matB%is_complex) THEN
       CALL ConvertMatrixToComplex(matB, converted_matrix)
       CALL IncrementMatrix(matA, converted_matrix, alpha, threshold)
       CALL CopyMatrix(converted_matrix, matB)
    ELSE IF (.NOT. matA%is_complex .AND. matB%is_complex) THEN
       CALL ConvertMatrixToComplex(matA, converted_matrix)
       CALL IncrementMatrix(converted_matrix, matB, alpha, threshold)
    ELSE IF (matA%is_complex .AND. matB%is_complex) THEN


       !$omp parallel
       !$omp do collapse(2)
       DO JJ = 1, matA%process_grid%number_of_blocks_columns
          DO II = 1, matA%process_grid%number_of_blocks_rows
             CALL IncrementMatrix(matA%local_data_c(II, JJ), matB%local_data_c(II, JJ), alpha, &
                  & threshold)
          END DO
       END DO
       !$omp end do
       !$omp end parallel

    ELSE


       !$omp parallel
       !$omp do collapse(2)
       DO JJ = 1, matA%process_grid%number_of_blocks_columns
          DO II = 1, matA%process_grid%number_of_blocks_rows
             CALL IncrementMatrix(matA%local_data_r(II, JJ), matB%local_data_r(II, JJ), alpha, &
                  & threshold)
          END DO
       END DO
       !$omp end do
       !$omp end parallel

    END IF

    CALL DestructMatrix(converted_matrix)

  END SUBROUTINE IncrementMatrix_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Will scale a distributed sparse matrix by a constant.
  SUBROUTINE ScaleMatrix_psr(this, constant)
    !> Matrix to scale.
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    !> A constant scale factor.
    REAL(NTREAL), INTENT(IN) :: constant
    !! Local Data
    INTEGER :: II, JJ

    IF (this%is_complex) THEN


       !$omp parallel
       !$omp do collapse(2)
       DO JJ = 1, this%process_grid%number_of_blocks_columns
          DO II = 1, this%process_grid%number_of_blocks_rows
             CALL ScaleMatrix(this%local_data_c(II, JJ), constant)
          END DO
       END DO
       !$omp end do
       !$omp end parallel

    ELSE


       !$omp parallel
       !$omp do collapse(2)
       DO JJ = 1, this%process_grid%number_of_blocks_columns
          DO II = 1, this%process_grid%number_of_blocks_rows
             CALL ScaleMatrix(this%local_data_r(II, JJ), constant)
          END DO
       END DO
       !$omp end do
       !$omp end parallel

    END IF

  END SUBROUTINE ScaleMatrix_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Will scale a distributed sparse matrix by a constant.
  RECURSIVE SUBROUTINE ScaleMatrix_psc(this, constant)
    !> Matrix to scale.
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    !> A constant scale factor.
    COMPLEX(NTCOMPLEX), INTENT(IN) :: constant
    !! Local Data
    TYPE(Matrix_ps) :: this_c
    INTEGER :: II, JJ

    IF (this%is_complex) THEN


       !$omp parallel
       !$omp do collapse(2)
       DO JJ = 1, this%process_grid%number_of_blocks_columns
          DO II = 1, this%process_grid%number_of_blocks_rows
             CALL ScaleMatrix(this%local_data_c(II, JJ), constant)
          END DO
       END DO
       !$omp end do
       !$omp end parallel

    ELSE
       CALL ConvertMatrixToComplex(this, this_c)
       CALL ScaleMatrix_psc(this_c, constant)
       CALL CopyMatrix(this_c, this)
       CALL DestructMatrix(this_c)
    END IF

  END SUBROUTINE ScaleMatrix_psc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Will scale a distributed sparse matrix by a constant.
  SUBROUTINE MatrixDiagonalScale_psr(this, tlist)
    !> Matrix to scale.
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    !> A constant scale factor.
    TYPE(TripletList_r), INTENT(IN) :: tlist
    !! Local Data
    TYPE(Matrix_lsr) :: lmat
    TYPE(TripletList_r) :: filtered
    TYPE(Triplet_r) :: trip


    INTEGER :: II, col

    !! Merge to the local block
    CALL MergeMatrixLocalBlocks(this, lmat)

    !! Filter out the triplets that aren't stored locally
    CALL ConstructTripletList(filtered)
    DO II = 1, tlist%CurrentSize
       CALL GetTripletAt(tlist, II, trip)
       col = trip%index_column
       IF (col .GE. this%start_column .AND. col .LT. this%end_column) THEN
          trip%index_column = trip%index_column - this%start_column + 1
          trip%index_row = trip%index_column
          CALL AppendToTripletList(filtered, trip)
       END IF
    END DO

    !! Scale
    CALL MatrixDiagonalScale(lmat, filtered)

    !! Split
    CALL SplitMatrixToLocalBlocks(this, lmat)
    CALL DestructMatrix(lmat)
    CALL DestructTripletList(filtered)
  END SUBROUTINE MatrixDiagonalScale_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Will scale a distributed sparse matrix by a constant.
  RECURSIVE SUBROUTINE MatrixDiagonalScale_psc(this, tlist)
    !> Matrix to scale.
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    !> A constant scale factor.
    TYPE(TripletList_c), INTENT(IN) :: tlist
    !! Local Data
    TYPE(Matrix_lsc) :: lmat
    TYPE(TripletList_c) :: filtered
    TYPE(Triplet_c) :: trip


    INTEGER :: II, col

    !! Merge to the local block
    CALL MergeMatrixLocalBlocks(this, lmat)

    !! Filter out the triplets that aren't stored locally
    CALL ConstructTripletList(filtered)
    DO II = 1, tlist%CurrentSize
       CALL GetTripletAt(tlist, II, trip)
       col = trip%index_column
       IF (col .GE. this%start_column .AND. col .LT. this%end_column) THEN
          trip%index_column = trip%index_column - this%start_column + 1
          trip%index_row = trip%index_column
          CALL AppendToTripletList(filtered, trip)
       END IF
    END DO

    !! Scale
    CALL MatrixDiagonalScale(lmat, filtered)

    !! Split
    CALL SplitMatrixToLocalBlocks(this, lmat)
    CALL DestructMatrix(lmat)
    CALL DestructTripletList(filtered)
  END SUBROUTINE MatrixDiagonalScale_psc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the trace of the matrix.
  SUBROUTINE MatrixTrace_psr(this, trace_value)
    !> The matrix to compute the trace of.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> The trace value of the full distributed sparse matrix.
    REAL(NTREAL), INTENT(OUT) :: trace_value
    !! Local data
    TYPE(TripletList_r) :: triplet_list_r
    TYPE(TripletList_c) :: triplet_list_c
    !! Counters/Temporary
    INTEGER :: II
    TYPE(Matrix_lsr) :: merged_local_data_r
    TYPE(Matrix_lsc) :: merged_local_data_c
    INTEGER :: ierr

    IF (this%is_complex) THEN




       !! Merge all the local data
       CALL MergeMatrixLocalBlocks(this, merged_local_data_c)

       !! Compute The Local Contribution
       trace_value = 0
       CALL MatrixToTripletList(merged_local_data_c, triplet_list_c)
       DO II = 1, triplet_list_c%CurrentSize
          IF (this%start_row + triplet_list_c%DATA(II)%index_row .EQ. &
               & this%start_column + triplet_list_c%DATA(II)%index_column) THEN
             trace_value = trace_value + &
                  & REAL(triplet_list_c%DATA(II)%point_value, NTREAL)
          END IF
       END DO

       !! Sum Among Process Slice
       CALL MPI_Allreduce(MPI_IN_PLACE, trace_value, 1, MPINTREAL, &
            & MPI_SUM, this%process_grid%within_slice_comm, ierr)

       CALL DestructMatrix(merged_local_data_c)



    ELSE




       !! Merge all the local data
       CALL MergeMatrixLocalBlocks(this, merged_local_data_r)

       !! Compute The Local Contribution
       trace_value = 0
       CALL MatrixToTripletList(merged_local_data_r, triplet_list_r)
       DO II = 1, triplet_list_r%CurrentSize
          IF (this%start_row + triplet_list_r%DATA(II)%index_row .EQ. &
               & this%start_column + triplet_list_r%DATA(II)%index_column) THEN
             trace_value = trace_value + &
                  & REAL(triplet_list_r%DATA(II)%point_value, NTREAL)
          END IF
       END DO

       !! Sum Among Process Slice
       CALL MPI_Allreduce(MPI_IN_PLACE, trace_value, 1, MPINTREAL, &
            & MPI_SUM, this%process_grid%within_slice_comm, ierr)

       CALL DestructMatrix(merged_local_data_r)



    END IF
  END SUBROUTINE MatrixTrace_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Measure the asymmetry of a matrix NORM(A - A.T)
  FUNCTION MeasureAsymmetry(this) RESULT(norm_value)
    !> The matrix to measure
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> The norm value of the full distributed sparse matrix.
    REAL(NTREAL) :: norm_value
    !! Local variables
    TYPE(Matrix_ps) :: tmat

    CALL TransposeMatrix(this, tmat)
    CALL ConjugateMatrix(tmat)
    CALL IncrementMatrix(this, tmat, alpha_in=-1.0_NTREAL)
    norm_value = MatrixNorm(tmat)
    CALL DestructMatrix(tmat)

  END FUNCTION MeasureAsymmetry
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Make the matrix symmetric
  SUBROUTINE SymmetrizeMatrix(this)
    !> The matrix to symmetrize
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    !! Local variables
    TYPE(Matrix_ps) :: tmat

    CALL TransposeMatrix(this, tmat)
    CALL ConjugateMatrix(tmat)
    CALL IncrementMatrix(tmat, this, alpha_in=1.0_NTREAL)
    CALL ScaleMatrix(this, 0.5_NTREAL)
    CALL DestructMatrix(tmat)

  END SUBROUTINE SymmetrizeMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Transform a matrix B = P * A * P^-1
  !! This routine will check if P is the identity matrix, and if so
  !! just return A.
  SUBROUTINE SimilarityTransform(A, P, PInv, ResMat, pool_in, threshold_in)
    !> The matrix to transform
    TYPE(Matrix_ps), INTENT(IN) :: A
    !> The left matrix.
    TYPE(Matrix_ps), INTENT(IN) :: P
    !> The right matrix.
    TYPE(Matrix_ps), INTENT(IN) :: PInv
    !> The computed matrix P * A * P^-1
    TYPE(Matrix_ps), INTENT(INOUT) :: ResMat
    !> A matrix memory pool.
    TYPE(MatrixMemoryPool_p), INTENT(INOUT), OPTIONAL :: pool_in
    !> The threshold for removing small elements.
    REAL(NTREAL), INTENT(IN), OPTIONAL :: threshold_in
    !! Local variables
    TYPE(MatrixMemoryPool_p) :: pool
    TYPE(Matrix_ps) :: TempMat
    REAL(NTREAL) :: threshold

    !! Optional Parameters
    IF (.NOT. PRESENT(threshold_in)) THEN
       threshold = 0.0_NTREAL
    ELSE
       threshold = threshold_in
    END IF
    IF (.NOT. PRESENT(pool_in)) THEN
       CALL ConstructMatrixMemoryPool(pool, A)
    END IF

    !! Check if P is the identity matrix, if so we can exit early.
    IF (IsIdentity(P)) THEN
       CALL CopyMatrix(A, ResMat)
    ELSE
       !! Compute
       IF (PRESENT(pool_in)) THEN
          CALL MatrixMultiply(P, A, TempMat, &
               & threshold_in = threshold, memory_pool_in = pool_in)
          CALL MatrixMultiply(TempMat, PInv, ResMat, &
               & threshold_in = threshold, memory_pool_in = pool_in)
       ELSE
          CALL MatrixMultiply(P, A, TempMat, &
               & threshold_in = threshold, memory_pool_in = pool)
          CALL MatrixMultiply(TempMat, PInv, ResMat, &
               & threshold_in = threshold, memory_pool_in = pool)
       END IF
    END IF

    !! Cleanup
    IF (.NOT. PRESENT(pool_in)) THEN
       CALL DestructMatrixMemoryPool(pool)
    END IF
    CALL DestructMatrix(TempMat)
  END SUBROUTINE SimilarityTransform
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE PSMatrixAlgebraModule
