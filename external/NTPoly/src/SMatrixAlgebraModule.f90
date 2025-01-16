!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for performing linear algebra using sparse matrices.
MODULE SMatrixAlgebraModule
  USE DataTypesModule, ONLY : NTREAL, NTCOMPLEX
  USE DMatrixModule, ONLY : Matrix_ldr, Matrix_ldc, ConstructMatrixDFromS, &
       & ConstructMatrixSFromD, CopyMatrix, MultiplyMatrix, TransposeMatrix, &
       & DestructMatrix
  USE MatrixMemoryPoolModule, ONLY : MatrixMemoryPool_lr, MatrixMemoryPool_lc, &
       & DestructMatrixMemoryPool, CheckMemoryPoolValidity, SetPoolSparsity, &
       & ConstructMatrixMemoryPool
  USE SMatrixModule, ONLY: Matrix_lsr, Matrix_lsc, DestructMatrix, CopyMatrix, &
       & TransposeMatrix, ConjugateMatrix, ConstructMatrixFromTripletList, &
       & ConstructEmptyMatrix
  USE SVectorModule, ONLY : AddSparseVectors, PairwiseMultiplyVectors
  USE TripletListModule, ONLY: TripletList_r, TripletList_c, SortTripletList, &
       & DestructTripletList, ConstructTripletList
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ScaleMatrix
  PUBLIC :: IncrementMatrix
  PUBLIC :: DotMatrix
  PUBLIC :: PairwiseMultiplyMatrix
  PUBLIC :: MatrixMultiply
  PUBLIC :: MatrixColumnNorm
  PUBLIC :: MatrixNorm
  PUBLIC :: MatrixGrandSum
  PUBLIC :: MatrixDiagonalScale
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE ScaleMatrix
     MODULE PROCEDURE ScaleMatrix_lsr
     MODULE PROCEDURE ScaleMatrix_lsc
     MODULE PROCEDURE ScaleMatrix_lsc_c
  END INTERFACE ScaleMatrix
  INTERFACE IncrementMatrix
     MODULE PROCEDURE IncrementMatrix_lsr
     MODULE PROCEDURE IncrementMatrix_lsc
  END INTERFACE IncrementMatrix
  INTERFACE DotMatrix
     MODULE PROCEDURE DotMatrix_lsr
     MODULE PROCEDURE DotMatrix_lsc
  END INTERFACE DotMatrix
  INTERFACE PairwiseMultiplyMatrix
     MODULE PROCEDURE PairwiseMultiplyMatrix_lsr
     MODULE PROCEDURE PairwiseMultiplyMatrix_lsc
  END INTERFACE PairwiseMultiplyMatrix
  INTERFACE MatrixMultiply
     MODULE PROCEDURE GemmMatrix_lsr
     MODULE PROCEDURE GemmMatrix_lsc
  END INTERFACE MatrixMultiply
  INTERFACE MatrixColumnNorm
     MODULE PROCEDURE MatrixColumnNorm_lsr
     MODULE PROCEDURE MatrixColumnNorm_lsc
  END INTERFACE MatrixColumnNorm
  INTERFACE MatrixNorm
     MODULE PROCEDURE MatrixNorm_lsr
     MODULE PROCEDURE MatrixNorm_lsc
  END INTERFACE MatrixNorm
  INTERFACE MatrixGrandSum
     MODULE PROCEDURE MatrixGrandSum_lsr
     MODULE PROCEDURE MatrixGrandSum_lsc
  END INTERFACE MatrixGrandSum
  INTERFACE MultiplyBlock
     MODULE PROCEDURE MultiplyBlock_lsr
     MODULE PROCEDURE MultiplyBlock_lsc
  END INTERFACE MultiplyBlock
  INTERFACE PruneList
     MODULE PROCEDURE PruneList_lsr
     MODULE PROCEDURE PruneList_lsc
  END INTERFACE PruneList
  INTERFACE SparseBranch
     MODULE PROCEDURE SparseBranch_lsr
     MODULE PROCEDURE SparseBranch_lsc
  END INTERFACE SparseBranch
  INTERFACE DenseBranch
     MODULE PROCEDURE DenseBranch_lsr
     MODULE PROCEDURE DenseBranch_lsc
  END INTERFACE DenseBranch
  INTERFACE MatrixDiagonalScale
     MODULE PROCEDURE MatrixDiagonalScale_lsr
     MODULE PROCEDURE MatrixDiagonalScale_lsc
  END INTERFACE MatrixDiagonalScale
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Will scale a sparse matrix by a constant.
  PURE SUBROUTINE ScaleMatrix_lsr(matA,constant)
    !> Matrix A.
    TYPE(Matrix_lsr), INTENT(INOUT) :: matA
    !> Constant scale factor.
    REAL(NTREAL), INTENT(IN) :: constant


    matA%values(:) = constant * matA%values
  END SUBROUTINE ScaleMatrix_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Will scale a sparse matrix by a constant.
  PURE SUBROUTINE ScaleMatrix_lsc(matA,constant)
    !> Matrix A.
    TYPE(Matrix_lsc), INTENT(INOUT) :: matA
    !> Constant scale factor.
    REAL(NTREAL), INTENT(IN) :: constant


    matA%values(:) = constant * matA%values
  END SUBROUTINE ScaleMatrix_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Will scale a sparse matrix by a constant.
  PURE SUBROUTINE ScaleMatrix_lsc_c(matA,constant)
    !> Matrix A.
    TYPE(Matrix_lsc), INTENT(INOUT) :: matA
    !> Constant scale factor.
    COMPLEX(NTCOMPLEX), INTENT(IN) :: constant


    matA%values(:) = constant * matA%values
  END SUBROUTINE ScaleMatrix_lsc_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Matrix B = alpha*Matrix A + Matrix B (AXPY).
  !> This will utilize the sparse vector addition routine.
  PURE SUBROUTINE IncrementMatrix_lsr(matA, matB, alpha_in, threshold_in)
    !> Matrix A.
    TYPE(Matrix_lsr), INTENT(IN)  :: matA
    !> Matrix B.
    TYPE(Matrix_lsr), INTENT(INOUT) :: matB
    !> Multiplier (default = 1.0).
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: alpha_in
    !> For flushing values to zero (default = 0).
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: threshold_in
    !! Local Variables
    TYPE(Matrix_lsr) :: matC


    !! Counter Variables
    INTEGER :: II
    INTEGER :: inner_a, inner_b
    INTEGER :: total_a, total_b, total_c
    !! Temporary Variables
    INTEGER :: indices_added_into_c
    REAL(NTREAL) :: alpha
    REAL(NTREAL) :: threshold
    INTEGER :: size_of_a, size_of_b

    !! Process Optional Parameters
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

    size_of_a = matA%outer_index(matA%columns + 1)

    !! Allocate sufficient space for matC
    CALL ConstructEmptyMatrix(matC, matA%rows, matA%columns)
    IF (ALLOCATED(matB%values)) THEN
       size_of_b = matB%outer_index(matB%columns + 1)
       ALLOCATE(matC%inner_index(size_of_a + size_of_b))
       ALLOCATE(matC%values(size_of_a + size_of_b))
    ELSE
       ALLOCATE(matC%inner_index(size_of_a))
       ALLOCATE(matC%values(size_of_a))
    END IF

    !! Perform loops
    total_a = 1
    total_b = 1
    total_c = 1
    DO II = 1, matA%columns
       !! Inner counters
       inner_a = matA%outer_index(II + 1) - matA%outer_index(II)
       inner_b = matB%outer_index(II+1) - matB%outer_index(II)
       CALL AddSparseVectors(&
            matA%inner_index(total_a:total_a + inner_a - 1), &
            matA%values(total_a:total_a + inner_a - 1), &
            matB%inner_index(total_b:total_b + inner_b - 1), &
            matB%values(total_b:total_b + inner_b - 1), &
            matC%inner_index(total_c:), matC%values(total_c:), &
            indices_added_into_c, alpha, threshold)
       matC%outer_index(II + 1) = matC%outer_index(II) + indices_added_into_c
       total_a = total_a + inner_a
       total_b = total_b + inner_b
       total_c = total_c + indices_added_into_c
    END DO

    !! Cleanup
    CALL DestructMatrix(matB)
    CALL ConstructEmptyMatrix(matB, matC%rows, matC%columns)
    matB%outer_index(:) = matC%outer_index
    ALLOCATE(matB%inner_index(matC%outer_index(matC%columns + 1)))
    ALLOCATE(matB%values(matC%outer_index(matC%columns + 1)))
    matB%inner_index(:) = matC%inner_index(:matC%outer_index(matC%columns + 1))
    matB%values(:) = matC%values(:matC%outer_index(matC%columns + 1))
    CALL DestructMatrix(matC)
  END SUBROUTINE IncrementMatrix_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Matrix B = alpha*Matrix A + Matrix B (AXPY).
  !> This will utilize the sparse vector addition routine.
  PURE SUBROUTINE IncrementMatrix_lsc(matA, matB, alpha_in, threshold_in)
    !> Matrix A.
    TYPE(Matrix_lsc), INTENT(IN)  :: matA
    !> Matrix B.
    TYPE(Matrix_lsc), INTENT(INOUT) :: matB
    !> Multiplier (default = 1.0).
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: alpha_in
    !> For flushing values to zero (default = 0).
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: threshold_in
    !! Local Variables
    TYPE(Matrix_lsc) :: matC


    !! Counter Variables
    INTEGER :: II
    INTEGER :: inner_a, inner_b
    INTEGER :: total_a, total_b, total_c
    !! Temporary Variables
    INTEGER :: indices_added_into_c
    REAL(NTREAL) :: alpha
    REAL(NTREAL) :: threshold
    INTEGER :: size_of_a, size_of_b

    !! Process Optional Parameters
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

    size_of_a = matA%outer_index(matA%columns + 1)

    !! Allocate sufficient space for matC
    CALL ConstructEmptyMatrix(matC, matA%rows, matA%columns)
    IF (ALLOCATED(matB%values)) THEN
       size_of_b = matB%outer_index(matB%columns + 1)
       ALLOCATE(matC%inner_index(size_of_a + size_of_b))
       ALLOCATE(matC%values(size_of_a + size_of_b))
    ELSE
       ALLOCATE(matC%inner_index(size_of_a))
       ALLOCATE(matC%values(size_of_a))
    END IF

    !! Perform loops
    total_a = 1
    total_b = 1
    total_c = 1
    DO II = 1, matA%columns
       !! Inner counters
       inner_a = matA%outer_index(II + 1) - matA%outer_index(II)
       inner_b = matB%outer_index(II+1) - matB%outer_index(II)
       CALL AddSparseVectors(&
            matA%inner_index(total_a:total_a + inner_a - 1), &
            matA%values(total_a:total_a + inner_a - 1), &
            matB%inner_index(total_b:total_b + inner_b - 1), &
            matB%values(total_b:total_b + inner_b - 1), &
            matC%inner_index(total_c:), matC%values(total_c:), &
            indices_added_into_c, alpha, threshold)
       matC%outer_index(II + 1) = matC%outer_index(II) + indices_added_into_c
       total_a = total_a + inner_a
       total_b = total_b + inner_b
       total_c = total_c + indices_added_into_c
    END DO

    !! Cleanup
    CALL DestructMatrix(matB)
    CALL ConstructEmptyMatrix(matB, matC%rows, matC%columns)
    matB%outer_index(:) = matC%outer_index
    ALLOCATE(matB%inner_index(matC%outer_index(matC%columns + 1)))
    ALLOCATE(matB%values(matC%outer_index(matC%columns + 1)))
    matB%inner_index(:) = matC%inner_index(:matC%outer_index(matC%columns + 1))
    matB%values(:) = matC%values(:matC%outer_index(matC%columns + 1))
    CALL DestructMatrix(matC)
  END SUBROUTINE IncrementMatrix_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Pairwise Multiply two matrices.
  !> This will utilize the sparse vector pairwise multiply routine.
  PURE SUBROUTINE PairwiseMultiplyMatrix_lsr(matA, matB, matC)
    !> Matrix A.
    TYPE(Matrix_lsr), INTENT(IN)  :: matA
    !> Matrix B.
    TYPE(Matrix_lsr), INTENT(IN) :: matB
    !> matC = MatA mult MatB.
    TYPE(Matrix_lsr), INTENT(INOUT) :: matC
    !! Local Variables
    TYPE(Matrix_lsr) :: TempMat


    !! Counter Variables
    INTEGER :: II
    INTEGER :: inner_a, inner_b
    INTEGER :: total_a, total_b, total_c
    !! Temporary Variables
    INTEGER :: indices_added_into_c
    INTEGER :: size_of_a, size_of_b

    CALL ConstructEmptyMatrix(TempMat, matA%rows, matA%columns)
    size_of_a = matA%outer_index(matA%columns + 1)
    size_of_b = matB%outer_index(matB%columns + 1)
    ALLOCATE(TempMat%inner_index(MIN(size_of_a, size_of_b)))
    ALLOCATE(TempMat%values(MIN(size_of_a, size_of_b)))

    !! Perform loops
    total_a = 1
    total_b = 1
    total_c = 1
    DO II = 1, matA%columns
       !! Inner counters
       inner_a = matA%outer_index(II + 1) - matA%outer_index(II)
       inner_b = matB%outer_index(II + 1) - matB%outer_index(II)
       CALL PairwiseMultiplyVectors(&
            matA%inner_index(total_a:total_a + inner_a - 1), &
            matA%values(total_a:total_a + inner_a - 1), &
            matB%inner_index(total_b:total_b + inner_b - 1), &
            matB%values(total_b:total_b + inner_b - 1), &
            TempMat%inner_index(total_c:), &
            TempMat%values(total_c:), &
            indices_added_into_c)
       TempMat%outer_index(II + 1) = &
            & TempMat%outer_index(II) + indices_added_into_c
       total_a = total_a + inner_a
       total_b = total_b + inner_b
       total_c = total_c + indices_added_into_c
    END DO

    !! Cleanup
    CALL DestructMatrix(matC)
    CALL ConstructEmptyMatrix(matC, TempMat%rows, TempMat%columns)
    matC%outer_index(:) = TempMat%outer_index
    ALLOCATE(matC%inner_index(TempMat%outer_index(TempMat%columns+1)))
    ALLOCATE(matC%values(TempMat%outer_index(TempMat%columns+1)))
    matC%inner_index(:) = TempMat%inner_index(&
         & :TempMat%outer_index(TempMat%columns+1))
    matC%values(:) = TempMat%values(:TempMat%outer_index(TempMat%columns+1))
    CALL DestructMatrix(TempMat)
  END SUBROUTINE PairwiseMultiplyMatrix_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Pairwise Multiply two matrices.
  !> This will utilize the sparse vector pairwise routine.
  PURE SUBROUTINE PairwiseMultiplyMatrix_lsc(matA, matB, matC)
    !> Matrix A.
    TYPE(Matrix_lsc), INTENT(IN)  :: matA
    !> Matrix B.
    TYPE(Matrix_lsc), INTENT(IN) :: matB
    !> matC = MatA mult MatB.
    TYPE(Matrix_lsc), INTENT(INOUT) :: matC
    !! Local Variables
    TYPE(Matrix_lsc) :: TempMat


    !! Counter Variables
    INTEGER :: II
    INTEGER :: inner_a, inner_b
    INTEGER :: total_a, total_b, total_c
    !! Temporary Variables
    INTEGER :: indices_added_into_c
    INTEGER :: size_of_a, size_of_b

    CALL ConstructEmptyMatrix(TempMat, matA%rows, matA%columns)
    size_of_a = matA%outer_index(matA%columns + 1)
    size_of_b = matB%outer_index(matB%columns + 1)
    ALLOCATE(TempMat%inner_index(MIN(size_of_a, size_of_b)))
    ALLOCATE(TempMat%values(MIN(size_of_a, size_of_b)))

    !! Perform loops
    total_a = 1
    total_b = 1
    total_c = 1
    DO II = 1, matA%columns
       !! Inner counters
       inner_a = matA%outer_index(II + 1) - matA%outer_index(II)
       inner_b = matB%outer_index(II + 1) - matB%outer_index(II)
       CALL PairwiseMultiplyVectors(&
            matA%inner_index(total_a:total_a + inner_a - 1), &
            matA%values(total_a:total_a + inner_a - 1), &
            matB%inner_index(total_b:total_b + inner_b - 1), &
            matB%values(total_b:total_b + inner_b - 1), &
            TempMat%inner_index(total_c:), &
            TempMat%values(total_c:), &
            indices_added_into_c)
       TempMat%outer_index(II + 1) = &
            & TempMat%outer_index(II) + indices_added_into_c
       total_a = total_a + inner_a
       total_b = total_b + inner_b
       total_c = total_c + indices_added_into_c
    END DO

    !! Cleanup
    CALL DestructMatrix(matC)
    CALL ConstructEmptyMatrix(matC, TempMat%rows, TempMat%columns)
    matC%outer_index(:) = TempMat%outer_index
    ALLOCATE(matC%inner_index(TempMat%outer_index(TempMat%columns+1)))
    ALLOCATE(matC%values(TempMat%outer_index(TempMat%columns+1)))
    matC%inner_index(:) = TempMat%inner_index(&
         & :TempMat%outer_index(TempMat%columns+1))
    matC%values(:) = TempMat%values(:TempMat%outer_index(TempMat%columns+1))
    CALL DestructMatrix(TempMat)
  END SUBROUTINE PairwiseMultiplyMatrix_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Product = sum(MatA[ij]*MatB[ij])
  PURE SUBROUTINE DotMatrix_lsr(matA, matB, product)
    !> Matrix A.
    TYPE(Matrix_lsr), INTENT(IN) :: matA
    !> Matrix B.
    TYPE(Matrix_lsr), INTENT(IN) :: matB
    !> Dot product.
    REAL(NTREAL), INTENT(OUT) :: product
    !! Local Variables
    TYPE(Matrix_lsr) :: matC

    CALL PairwiseMultiplyMatrix(matA, matB, matC)

    CALL MatrixGrandSum(matC, product)
    CALL DestructMatrix(matC)

  END SUBROUTINE DotMatrix_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Product = sum(MatA^H[ij]*MatB[ij])
  PURE SUBROUTINE DotMatrix_lsc(matA, matB, product)
    !> Matrix A.
    TYPE(Matrix_lsc), INTENT(IN) :: matA
    !> Matrix B.
    TYPE(Matrix_lsc), INTENT(IN) :: matB
    !> Dot product.
    COMPLEX(NTCOMPLEX), INTENT(OUT) :: product
    !! Local Variables
    TYPE(Matrix_lsc) :: matC
    TYPE(Matrix_lsc) :: matAH

    CALL CopyMatrix(matA, matAH)
    CALL ConjugateMatrix(matAH)

    CALL PairwiseMultiplyMatrix(matAH, matB, matC)
    CALL MatrixGrandSum(matC, product)

    CALL DestructMatrix(matC)
    CALL DestructMatrix(matAH)

  END SUBROUTINE DotMatrix_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Multiply two matrices together, and add to the third.
  !> C := alpha*matA*op( matB ) + beta*matC
  SUBROUTINE GemmMatrix_lsr(matA, matB, matC, IsATransposed_in, &
       & IsBTransposed_in, alpha_in, beta_in, threshold_in, &
       & blocked_memory_pool_in)
    !> Matrix A.
    TYPE(Matrix_lsr), INTENT(IN)  :: matA
    !> Matrix B.
    TYPE(Matrix_lsr), INTENT(IN)  :: matB
    !> matC = alpha*matA*op( matB ) + beta*matC.
    TYPE(Matrix_lsr), INTENT(INOUT) :: matC
    !> True if A is already transposed.
    LOGICAL, OPTIONAL, INTENT(IN) :: IsATransposed_in
    !> True if B is already transposed.
    LOGICAL, OPTIONAL, INTENT(IN) :: IsBTransposed_in
    !> Scales the multiplication.
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: alpha_in
    !> Scales matrix we sum on to.
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: beta_in
    !> For flushing values to zero. Default value is 0.0.
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: threshold_in
    !> An optional memory pool for doing the calculation.
    TYPE(MatrixMemoryPool_lr), OPTIONAL, &
         & INTENT(INOUT), TARGET :: blocked_memory_pool_in
    !! Intermediate Data
    TYPE(Matrix_lsr) :: matAB
    LOGICAL :: IsATransposed, IsBTransposed
    REAL(NTREAL) :: alpha
    REAL(NTREAL) :: beta
    REAL(NTREAL) :: threshold
    TYPE(MatrixMemoryPool_lr) :: blocked_memory_pool


    REAL(NTREAL), PARAMETER :: sparsity_threshold = 0.1_NTREAL
    !! Counters and temporary data
    INTEGER :: mat_c_columns, mat_c_rows
    !! For Efficiency Purposes
    REAL(NTREAL) :: sparsity_a, sparsity_b
    REAL(NTREAL) :: sparsity_estimate
    LOGICAL :: pool_flag

    !! Process Optional Parameters
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
    IF (.NOT. PRESENT(IsATransposed_in)) THEN
       IsATransposed = .FALSE.
    ELSE
       IsATransposed = IsATransposed_in
    END IF
    IF (.NOT. PRESENT(IsBTransposed_in)) THEN
       IsBTransposed = .FALSE.
    ELSE
       IsBTransposed = IsBTransposed_in
    END IF
    IF (.NOT. PRESENT(threshold_in)) THEN
       threshold = 0.0
    ELSE
       threshold = threshold_in
    END IF

    !! Storage details for result matrix
    IF (IsATransposed) THEN
       mat_c_rows = matA%columns
    ELSE
       mat_c_rows = matA%rows
    END IF
    IF (IsBTransposed) THEN
       mat_c_columns = matB%rows
    ELSE
       mat_c_columns = matB%columns
    END IF

    !! Initialization of Memory
    sparsity_a = DBLE(SIZE(matA%values)) / (matA%rows * matA%columns)
    sparsity_b = DBLE(SIZE(matB%values)) / (matB%rows * matB%columns)
    sparsity_estimate = 4*MAX(sparsity_a, sparsity_b)
    IF (sparsity_estimate > 1.0) THEN
       sparsity_estimate = 1.0
    ELSE IF (sparsity_estimate < 1e-8) THEN
       sparsity_estimate = 1e-8
    END IF

    !! Decide whether to do dense or sparse version.
    IF (MIN(sparsity_a, sparsity_b) .GT. sparsity_threshold) THEN
       CALL DenseBranch(matA, matB, matAB, IsATransposed, IsBTransposed, &
            & alpha, threshold)
    ELSE
       !! Setup the memory pool
       IF (.NOT. PRESENT(blocked_memory_pool_in)) THEN
          CALL ConstructMatrixMemoryPool(blocked_memory_pool, mat_c_columns, &
               & mat_c_rows, sparsity_estimate)
          pool_flag = .FALSE.
       ELSEIF (.NOT. CheckMemoryPoolValidity(blocked_memory_pool_in, &
            & mat_c_columns, mat_c_rows)) THEN
          CALL DestructMatrixMemoryPool(blocked_memory_pool_in)
          CALL ConstructMatrixMemoryPool(blocked_memory_pool_in, mat_c_columns, &
               & mat_c_rows, sparsity_estimate)
          pool_flag = .TRUE.
       ELSE
          CALL SetPoolSparsity(blocked_memory_pool_in, sparsity_estimate)
          pool_flag = .TRUE.
       END IF
       !! Multiply
       IF (pool_flag) THEN
          CALL SparseBranch(matA, matB, matAB, IsATransposed, IsBTransposed, &
               & alpha, threshold, blocked_memory_pool_in)
       ELSE
          CALL SparseBranch(matA, matB, matAB, IsATransposed, IsBTransposed, &
               & alpha, threshold, blocked_memory_pool)
       END IF
    END IF

    !! Handle the add part of GEMM
    IF (PRESENT(beta_in)) THEN
       IF (ABS(beta_in) .GT. 0) THEN
          CALL ScaleMatrix(matC, beta)
          CALL IncrementMatrix(matAB, matC)
       ELSE
          CALL CopyMatrix(matAB, matC)
       END IF
    ELSE
       CALL CopyMatrix(matAB, matC)
    END IF

    CALL DestructMatrix(matAB)
    CALL DestructMatrixMemoryPool(blocked_memory_pool)
  END SUBROUTINE GemmMatrix_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Multiply two matrices together, and add to the third.
  !> C := alpha*matA*op( matB ) + beta*matC
  SUBROUTINE GemmMatrix_lsc(matA, matB, matC, IsATransposed_in, &
       & IsBTransposed_in, alpha_in, beta_in, threshold_in, &
       & blocked_memory_pool_in)
    !> Matrix A.
    TYPE(Matrix_lsc), INTENT(IN)  :: matA
    !> Matrix B.
    TYPE(Matrix_lsc), INTENT(IN)  :: matB
    !> matC = alpha*matA*op( matB ) + beta*matC.
    TYPE(Matrix_lsc), INTENT(INOUT) :: matC
    !> True if A is already transposed.
    LOGICAL, OPTIONAL, INTENT(IN) :: IsATransposed_in
    !> True if B is already transposed.
    LOGICAL, OPTIONAL, INTENT(IN) :: IsBTransposed_in
    !> Scales the multiplication.
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: alpha_in
    !> Scales matrix we sum on to.
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: beta_in
    !> For flushing values to zero. Default value is 0.0.
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: threshold_in
    !> An optional memory pool for doing the calculation.
    TYPE(MatrixMemoryPool_lc), OPTIONAL, &
         & INTENT(INOUT), TARGET :: blocked_memory_pool_in
    !! Intermediate Data
    TYPE(Matrix_lsc) :: matAB
    LOGICAL :: IsATransposed, IsBTransposed
    REAL(NTREAL) :: alpha
    REAL(NTREAL) :: beta
    REAL(NTREAL) :: threshold
    TYPE(MatrixMemoryPool_lc) :: blocked_memory_pool


    REAL(NTREAL), PARAMETER :: sparsity_threshold = 0.1_NTREAL
    !! Counters and temporary data
    INTEGER :: mat_c_columns, mat_c_rows
    !! For Efficiency Purposes
    REAL(NTREAL) :: sparsity_a, sparsity_b
    REAL(NTREAL) :: sparsity_estimate
    LOGICAL :: pool_flag

    !! Process Optional Parameters
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
    IF (.NOT. PRESENT(IsATransposed_in)) THEN
       IsATransposed = .FALSE.
    ELSE
       IsATransposed = IsATransposed_in
    END IF
    IF (.NOT. PRESENT(IsBTransposed_in)) THEN
       IsBTransposed = .FALSE.
    ELSE
       IsBTransposed = IsBTransposed_in
    END IF
    IF (.NOT. PRESENT(threshold_in)) THEN
       threshold = 0.0
    ELSE
       threshold = threshold_in
    END IF

    !! Storage details for result matrix
    IF (IsATransposed) THEN
       mat_c_rows = matA%columns
    ELSE
       mat_c_rows = matA%rows
    END IF
    IF (IsBTransposed) THEN
       mat_c_columns = matB%rows
    ELSE
       mat_c_columns = matB%columns
    END IF

    !! Initialization of Memory
    sparsity_a = DBLE(SIZE(matA%values)) / (matA%rows * matA%columns)
    sparsity_b = DBLE(SIZE(matB%values)) / (matB%rows * matB%columns)
    sparsity_estimate = 4*MAX(sparsity_a, sparsity_b)
    IF (sparsity_estimate > 1.0) THEN
       sparsity_estimate = 1.0
    ELSE IF (sparsity_estimate < 1e-8) THEN
       sparsity_estimate = 1e-8
    END IF

    !! Decide whether to do dense or sparse version.
    IF (MIN(sparsity_a, sparsity_b) .GT. sparsity_threshold) THEN
       CALL DenseBranch(matA, matB, matAB, IsATransposed, IsBTransposed, &
            & alpha, threshold)
    ELSE
       !! Setup the memory pool
       IF (.NOT. PRESENT(blocked_memory_pool_in)) THEN
          CALL ConstructMatrixMemoryPool(blocked_memory_pool, mat_c_columns, &
               & mat_c_rows, sparsity_estimate)
          pool_flag = .FALSE.
       ELSEIF (.NOT. CheckMemoryPoolValidity(blocked_memory_pool_in, &
            & mat_c_columns, mat_c_rows)) THEN
          CALL DestructMatrixMemoryPool(blocked_memory_pool_in)
          CALL ConstructMatrixMemoryPool(blocked_memory_pool_in, mat_c_columns, &
               & mat_c_rows, sparsity_estimate)
          pool_flag = .TRUE.
       ELSE
          CALL SetPoolSparsity(blocked_memory_pool_in, sparsity_estimate)
          pool_flag = .TRUE.
       END IF
       !! Multiply
       IF (pool_flag) THEN
          CALL SparseBranch(matA, matB, matAB, IsATransposed, IsBTransposed, &
               & alpha, threshold, blocked_memory_pool_in)
       ELSE
          CALL SparseBranch(matA, matB, matAB, IsATransposed, IsBTransposed, &
               & alpha, threshold, blocked_memory_pool)
       END IF
    END IF

    !! Handle the add part of GEMM
    IF (PRESENT(beta_in)) THEN
       IF (ABS(beta_in) .GT. 0) THEN
          CALL ScaleMatrix(matC, beta)
          CALL IncrementMatrix(matAB, matC)
       ELSE
          CALL CopyMatrix(matAB, matC)
       END IF
    ELSE
       CALL CopyMatrix(matAB, matC)
    END IF

    CALL DestructMatrix(matAB)
    CALL DestructMatrixMemoryPool(blocked_memory_pool)
  END SUBROUTINE GemmMatrix_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the norm of a sparse matrix along the columns.
  PURE SUBROUTINE MatrixColumnNorm_lsr(this, norm_per_column)
    !> The matrix to compute the norm of.
    TYPE(Matrix_lsr), INTENT(IN) :: this
    !> The norm value for each column in this matrix.
    REAL(NTREAL), DIMENSION(this%columns), INTENT(OUT) :: norm_per_column
    !! Local Data
    REAL(NTREAL) :: temp_value


    !! Local Data
    INTEGER :: II, JJ
    INTEGER :: elements_per_inner

    !! Allocate Space For Result
    norm_per_column = 0

    !! Iterate Over Local Data
    DO II = 1, this%columns
       elements_per_inner = this%outer_index(II + 1) - this%outer_index(II)
       DO JJ = 1, elements_per_inner
          temp_value = this%values(this%outer_index(II) + JJ)
          norm_per_column(II) = norm_per_column(II) + ABS(temp_value)
       END DO
    END DO
  END SUBROUTINE MatrixColumnNorm_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the norm of a sparse matrix along the columns.
  PURE SUBROUTINE MatrixColumnNorm_lsc(this, norm_per_column)
    !> The matrix to compute the norm of.
    TYPE(Matrix_lsc), INTENT(IN) :: this
    !> The norm value for each column in this matrix.
    REAL(NTREAL), DIMENSION(this%columns), INTENT(OUT) :: norm_per_column
    !! Local Data
    COMPLEX(NTCOMPLEX)  :: temp_value


    !! Local Data
    INTEGER :: II, JJ
    INTEGER :: elements_per_inner

    !! Allocate Space For Result
    norm_per_column = 0

    !! Iterate Over Local Data
    DO II = 1, this%columns
       elements_per_inner = this%outer_index(II + 1) - this%outer_index(II)
       DO JJ = 1, elements_per_inner
          temp_value = this%values(this%outer_index(II) + JJ)
          norm_per_column(II) = norm_per_column(II) + ABS(temp_value)
       END DO
    END DO
  END SUBROUTINE MatrixColumnNorm_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the 1 norm of a sparse matrix.
  PURE FUNCTION MatrixNorm_lsr(this) RESULT(norm)
    !> The matrix to compute the norm of.
    TYPE(Matrix_lsr), INTENT(IN) :: this
    !> The norm of the matrix.
    REAL(NTREAL) :: norm
    !! Local Variables
    REAL(NTREAL), DIMENSION(this%columns) :: column


    CALL MatrixColumnNorm(this, column)
    norm = MAXVAL(column)

  END FUNCTION MatrixNorm_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the 1 norm of a sparse matrix.
  PURE FUNCTION MatrixNorm_lsc(this) RESULT(norm)
    !> The matrix to compute the norm of.
    TYPE(Matrix_lsc), INTENT(IN) :: this
    !> The norm of the matrix.
    REAL(NTREAL) :: norm
    !! Local Variables
    REAL(NTREAL), DIMENSION(this%columns) :: column


    CALL MatrixColumnNorm(this, column)
    norm = MAXVAL(column)

  END FUNCTION MatrixNorm_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sum the elements of a matrix
  PURE SUBROUTINE MatrixGrandSum_lsr(this, sum_value)
    !> The matrix to sum
    TYPE(Matrix_lsr), INTENT(IN) :: this
    !> The sum of the matrix elements
    REAL(NTREAL), INTENT(OUT) :: sum_value


    sum_value = SUM(this%values)

  END SUBROUTINE MatrixGrandSum_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sum the elements of a matrix
  PURE SUBROUTINE MatrixGrandSum_lsc(this, sum_value)
    !> The matrix to sum
    TYPE(Matrix_lsc), INTENT(IN) :: this
    !> The sum of the matrix elements
    COMPLEX(NTCOMPLEX), INTENT(OUT) :: sum_value


    sum_value = SUM(this%values)

  END SUBROUTINE MatrixGrandSum_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Calculates the matrix product if using sparse-sparse algorithm.
  PURE SUBROUTINE SparseBranch_lsr(matA, matB, matC, IsATransposed, &
       & IsBTransposed, alpha, threshold, blocked_memory_pool)
    !> Matrix A.
    TYPE(Matrix_lsr), INTENT(IN)  :: matA
    !> Matrix B
    TYPE(Matrix_lsr), INTENT(IN)  :: matB
    !> matC = alpha*matA*op( matB ) + beta*matC.
    TYPE(Matrix_lsr), INTENT(INOUT) :: matC
    !> True if A is transposed.
    LOGICAL, INTENT(IN) :: IsATransposed
    !> True if B is transposed.
    LOGICAL, INTENT(IN) :: IsBTransposed
    !> Scaling value.
    REAL(NTREAL), INTENT(IN) :: alpha
    !> Threshold for flushing values.
    REAL(NTREAL), INTENT(IN) :: threshold
    !> Memory pool.
    TYPE(MatrixMemoryPool_lr), INTENT(INOUT) :: blocked_memory_pool
    !! Local Data
    TYPE(Matrix_lsr) :: matAT, matBT


    !! Block A and B
    IF (.NOT. IsATransposed) THEN
       CALL TransposeMatrix(matA, matAT)
    END IF
    IF (.NOT. IsBTransposed) THEN
       CALL TransposeMatrix(matB, matBT)
    END IF

    IF (IsATransposed .AND. IsBTransposed) THEN
       CALL MultiplyBlock(matA, matB, blocked_memory_pool)
    ELSEIF (IsATransposed) THEN
       CALL MultiplyBlock(matA, matBT, blocked_memory_pool)
    ELSEIF (IsBTransposed) THEN
       CALL MultiplyBlock(matAT, matB, blocked_memory_pool)
    ELSE
       CALL MultiplyBlock(matAT, matBT, blocked_memory_pool)
    END IF

    !! Go from triplets to return matrix
    CALL PruneList(blocked_memory_pool, alpha, threshold, &
         & blocked_memory_pool%columns, blocked_memory_pool%rows, matC)
  END SUBROUTINE SparseBranch_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Calculates the matrix product if using the sparse-sparse algorithm.
  PURE SUBROUTINE SparseBranch_lsc(matA, matB, matC, IsATransposed, &
       & IsBTransposed, alpha, threshold, blocked_memory_pool)
    !> Matrix A.
    TYPE(Matrix_lsc), INTENT(IN)  :: matA
    !> Matrix B
    TYPE(Matrix_lsc), INTENT(IN)  :: matB
    !> matC = alpha*matA*op( matB ) + beta*matC.
    TYPE(Matrix_lsc), INTENT(INOUT) :: matC
    !> True if A is transposed.
    LOGICAL, INTENT(IN) :: IsATransposed
    !> True if B is transposed.
    LOGICAL, INTENT(IN) :: IsBTransposed
    !> Scaling value.
    REAL(NTREAL), INTENT(IN) :: alpha
    !> Threshold for flushing values.
    REAL(NTREAL), INTENT(IN) :: threshold
    !> Memory pool.
    TYPE(MatrixMemoryPool_lc), INTENT(INOUT) :: blocked_memory_pool
    !! Local Data
    TYPE(Matrix_lsc) :: matAT, matBT


    !! Block A and B
    IF (.NOT. IsATransposed) THEN
       CALL TransposeMatrix(matA, matAT)
    END IF
    IF (.NOT. IsBTransposed) THEN
       CALL TransposeMatrix(matB, matBT)
    END IF

    IF (IsATransposed .AND. IsBTransposed) THEN
       CALL MultiplyBlock(matA, matB, blocked_memory_pool)
    ELSEIF (IsATransposed) THEN
       CALL MultiplyBlock(matA, matBT, blocked_memory_pool)
    ELSEIF (IsBTransposed) THEN
       CALL MultiplyBlock(matAT, matB, blocked_memory_pool)
    ELSE
       CALL MultiplyBlock(matAT, matBT, blocked_memory_pool)
    END IF

    !! Go from triplets to return matrix
    CALL PruneList(blocked_memory_pool, alpha, threshold, &
         & blocked_memory_pool%columns, blocked_memory_pool%rows, matC)
  END SUBROUTINE SparseBranch_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Calculate the matrix product using the dense-dense algorithm.
  SUBROUTINE DenseBranch_lsr(matA, matB, matC, IsATransposed, IsBTransposed, &
       & alpha, threshold)
    !> Matrix A.
    TYPE(Matrix_lsr), INTENT(IN)  :: matA
    !> Matrix B
    TYPE(Matrix_lsr), INTENT(IN)  :: matB
    !> matC = alpha*matA*op( matB ) + beta*matC.
    TYPE(Matrix_lsr), INTENT(INOUT) :: matC
    !> True if A is transposed.
    LOGICAL, INTENT(IN) :: IsATransposed
    !> True if B is transposed.
    LOGICAL, INTENT(IN) :: IsBTransposed
    !> Scaling value.
    REAL(NTREAL), INTENT(IN) :: alpha
    !> Threshold for flushing values.
    REAL(NTREAL), INTENT(IN) :: threshold
    !! Local Data
    TYPE(Matrix_ldr) :: DenseA
    TYPE(Matrix_ldr) :: DenseB
    TYPE(Matrix_ldr) :: DenseC


    !! Convert Forward
    CALL ConstructMatrixDFromS(matA, DenseA)
    CALL ConstructMatrixDFromS(matB, DenseB)

    !! Multiply
    CALL MultiplyMatrix(DenseA, DenseB, DenseC, &
         & IsATransposed_in = IsATransposed, IsBTransposed_in = IsBTransposed)

    !! Cleanup Intermediate
    CALL DestructMatrix(DenseA)
    CALL DestructMatrix(DenseB)

    !! Convert Back
    CALL ConstructMatrixSFromD(DenseC, matC, threshold)
    CALL ScaleMatrix(matC, alpha)

    !! Cleanup
    CALL DestructMatrix(DenseC)
  END SUBROUTINE DenseBranch_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Calculate the matrix product using the dense-dense algorithm.
  SUBROUTINE DenseBranch_lsc(matA, matB, matC, IsATransposed, IsBTransposed, &
       & alpha, threshold)
    !> Matrix A.
    TYPE(Matrix_lsc), INTENT(IN)  :: matA
    !> Matrix B
    TYPE(Matrix_lsc), INTENT(IN)  :: matB
    !> matC = alpha*matA*op( matB ) + beta*matC.
    TYPE(Matrix_lsc), INTENT(INOUT) :: matC
    !> True if A is transposed.
    LOGICAL, INTENT(IN) :: IsATransposed
    !> True if B is transposed.
    LOGICAL, INTENT(IN) :: IsBTransposed
    !> Scaling value.
    REAL(NTREAL), INTENT(IN) :: alpha
    !> Threshold for flushing values.
    REAL(NTREAL), INTENT(IN) :: threshold
    !! Local Data
    TYPE(Matrix_ldc) :: DenseA
    TYPE(Matrix_ldc) :: DenseB
    TYPE(Matrix_ldc) :: DenseC


    !! Convert Forward
    CALL ConstructMatrixDFromS(matA, DenseA)
    CALL ConstructMatrixDFromS(matB, DenseB)

    !! Multiply
    CALL MultiplyMatrix(DenseA, DenseB, DenseC, &
         & IsATransposed_in = IsATransposed, IsBTransposed_in = IsBTransposed)

    !! Cleanup Intermediate
    CALL DestructMatrix(DenseA)
    CALL DestructMatrix(DenseB)

    !! Convert Back
    CALL ConstructMatrixSFromD(DenseC, matC, threshold)
    CALL ScaleMatrix(matC, alpha)

    !! Cleanup
    CALL DestructMatrix(DenseC)
  END SUBROUTINE DenseBranch_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Multiplies a single block fo sparse-sparse.
  PURE SUBROUTINE MultiplyBlock_lsr(matAT,matBT,memorypool)
    !> Matrix A, already transposed.
    TYPE(Matrix_lsr), INTENT(IN)  :: matAT
    !> Matrix B, already transposed.
    TYPE(Matrix_lsr), INTENT(IN)  :: matBT
    !> Memory pool to multiply into.
    TYPE(MatrixMemoryPool_lr), INTENT(INOUT) :: memorypool
    !! Temp Variables
    REAL(NTREAL) :: val_a, val_b, val_c


    INTEGER :: inserted_vals
    INTEGER :: idx_a, idx_b, idx_hash
    INTEGER :: elements_per_inner_a, elements_per_inner_b
    LOGICAL :: is_set
    !! Counters
    INTEGER :: II, AA, BB

    !! Multiply
    DO II = 1, matAT%columns
       elements_per_inner_a = matAT%outer_index(II + 1) - &
            & matAT%outer_index(II)
       DO AA = 1, elements_per_inner_a
          val_a = matAT%values(matAT%outer_index(II) + AA)
          idx_a = matAT%inner_index(matAT%outer_index(II) + AA)
          elements_per_inner_b = matBT%outer_index(idx_a + 1) - &
               & matBT%outer_index(idx_a)
          DO BB = 1, elements_per_inner_b
             idx_b = matBT%inner_index(matBT%outer_index(idx_a) + BB)
             val_b = matBT%values(matBT%outer_index(idx_a)+ BB)
             val_c = memorypool%value_array(idx_b, II)
             is_set = memorypool%dirty_array(idx_b, II)
             IF (is_set .EQV. .FALSE.) THEN
                memorypool%dirty_array(idx_b, II) = .TRUE.
                idx_hash = (idx_b - 1) / memorypool%hash_size
                inserted_vals = & 
                     & memorypool%inserted_per_bucket(idx_hash + 1, II) + 1
                memorypool%inserted_per_bucket(idx_hash + 1, II) = &
                     & inserted_vals
                memorypool%hash_index(&
                     & inserted_vals + idx_hash * memorypool%hash_size, &
                     & II) = idx_b
             END IF
             memorypool%value_array(idx_b, II) = val_c + val_a * val_b
          END DO
       END DO
    END DO
  END SUBROUTINE MultiplyBlock_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Multiplies a single block fo sparse-sparse.
  PURE SUBROUTINE MultiplyBlock_lsc(matAT,matBT,memorypool)
    !> Matrix A, already transposed.
    TYPE(Matrix_lsc), INTENT(IN)  :: matAT
    !> Matrix B, already transposed.
    TYPE(Matrix_lsc), INTENT(IN)  :: matBT
    !> Memory pool to multiply into.
    TYPE(MatrixMemoryPool_lc), INTENT(INOUT) :: memorypool
    !! Temp Variables
    COMPLEX(NTCOMPLEX) :: val_a, val_b, val_c


    INTEGER :: inserted_vals
    INTEGER :: idx_a, idx_b, idx_hash
    INTEGER :: elements_per_inner_a, elements_per_inner_b
    LOGICAL :: is_set
    !! Counters
    INTEGER :: II, AA, BB

    !! Multiply
    DO II = 1, matAT%columns
       elements_per_inner_a = matAT%outer_index(II + 1) - &
            & matAT%outer_index(II)
       DO AA = 1, elements_per_inner_a
          val_a = matAT%values(matAT%outer_index(II) + AA)
          idx_a = matAT%inner_index(matAT%outer_index(II) + AA)
          elements_per_inner_b = matBT%outer_index(idx_a + 1) - &
               & matBT%outer_index(idx_a)
          DO BB = 1, elements_per_inner_b
             idx_b = matBT%inner_index(matBT%outer_index(idx_a) + BB)
             val_b = matBT%values(matBT%outer_index(idx_a)+ BB)
             val_c = memorypool%value_array(idx_b, II)
             is_set = memorypool%dirty_array(idx_b, II)
             IF (is_set .EQV. .FALSE.) THEN
                memorypool%dirty_array(idx_b, II) = .TRUE.
                idx_hash = (idx_b - 1) / memorypool%hash_size
                inserted_vals = & 
                     & memorypool%inserted_per_bucket(idx_hash + 1, II) + 1
                memorypool%inserted_per_bucket(idx_hash + 1, II) = &
                     & inserted_vals
                memorypool%hash_index(&
                     & inserted_vals + idx_hash * memorypool%hash_size, &
                     & II) = idx_b
             END IF
             memorypool%value_array(idx_b, II) = val_c + val_a * val_b
          END DO
       END DO
    END DO
  END SUBROUTINE MultiplyBlock_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Prunes out the values of the hash table into the matrix.
  PURE SUBROUTINE PruneList_lsr(memorypool,alpha,threshold, mat_c_columns, &
       & mat_c_rows, matAB)
    !> Memory pool to prune from.
    TYPE(MatrixMemoryPool_lr), INTENT(INOUT) :: memorypool
    !> Scaling value.
    REAL(NTREAL), INTENT(IN) :: alpha
    !> Threshold for flushing values to zero.
    REAL(NTREAL), INTENT(IN) :: threshold
    !> Size of the matrix we computed (columns).
    INTEGER, INTENT(IN) :: mat_c_columns
    !> Size of the matrix we computed (rows).
    INTEGER, INTENT(IN) :: mat_c_rows
    !> Sparse matrix to prune out into.
    TYPE(Matrix_lsr), INTENT(INOUT) :: matAB
    !! Local data
    REAL(NTREAL) :: working_value
    TYPE(TripletList_r) :: unsorted_pruned_list
    TYPE(TripletList_r) :: sorted_pruned_list


    !! Local data
    INTEGER :: working_col
    INTEGER :: values_per_hash
    INTEGER :: PII, HII, RII, CII

    !! Loop over the hash structure
    PII = 1
    DO RII = 1, mat_c_rows
       DO CII = 1, (mat_c_columns - 1) / memorypool%hash_size + 1
          !! Sort the elements in a hash
          values_per_hash = memorypool%inserted_per_bucket(CII, RII)
          memorypool%inserted_per_bucket(CII, RII) = 0
          !! Copy them
          DO HII = 1, values_per_hash
             working_col = memorypool%hash_index(HII + &
                  & (CII - 1) * memorypool%hash_size, RII)
             working_value = &
                  & memorypool%value_array(working_col, RII)
             memorypool%value_array(working_col, RII) = 0
             memorypool%dirty_array(working_col, RII) = .FALSE.
             !! If above threshold, insert
             IF (ABS(alpha*working_value) .GT. threshold) THEN
                memorypool%pruned_list(PII)%point_value = alpha*working_value
                memorypool%pruned_list(PII)%index_column = working_col
                memorypool%pruned_list(PII)%index_row = RII
                PII = PII + 1
             END IF
          END DO
       END DO
    END DO

    !! Convert to matrix
    CALL ConstructTripletList(unsorted_pruned_list, PII - 1)
    unsorted_pruned_list%DATA(:) = memorypool%pruned_list(1:PII - 1)
    CALL SortTripletList(unsorted_pruned_list, mat_c_columns, mat_c_rows, &
         & sorted_pruned_list, bubble_in = .TRUE.)
    CALL ConstructMatrixFromTripletList(matAB, sorted_pruned_list, mat_c_rows, &
         & mat_c_columns)
    CALL DestructTripletList(sorted_pruned_list)
    CALL DestructTripletList(unsorted_pruned_list)
  END SUBROUTINE PruneList_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Prunes out the values of the hash table into the matrix.
  PURE SUBROUTINE PruneList_lsc(memorypool,alpha,threshold, &
       & mat_c_columns, mat_c_rows, matAB)
    !> Memory pool to prune from.
    TYPE(MatrixMemoryPool_lc), INTENT(INOUT) :: memorypool
    !> Scaling value.
    REAL(NTREAL), INTENT(IN) :: alpha
    !> Threshold for flushing values to zero.
    REAL(NTREAL), INTENT(IN) :: threshold
    !> Size of the matrix we computed (columns).
    INTEGER, INTENT(IN) :: mat_c_columns
    !> Size of the matrix we computed (rows).
    INTEGER, INTENT(IN) :: mat_c_rows
    !> Sparse matrix to prune out into.
    TYPE(Matrix_lsc), INTENT(INOUT) :: matAB
    !! Local data
    COMPLEX(NTCOMPLEX) :: working_value
    TYPE(TripletList_c) :: unsorted_pruned_list
    TYPE(TripletList_c) :: sorted_pruned_list


    !! Local data
    INTEGER :: working_col
    INTEGER :: values_per_hash
    INTEGER :: PII, HII, RII, CII

    !! Loop over the hash structure
    PII = 1
    DO RII = 1, mat_c_rows
       DO CII = 1, (mat_c_columns - 1) / memorypool%hash_size + 1
          !! Sort the elements in a hash
          values_per_hash = memorypool%inserted_per_bucket(CII, RII)
          memorypool%inserted_per_bucket(CII, RII) = 0
          !! Copy them
          DO HII = 1, values_per_hash
             working_col = memorypool%hash_index(HII + &
                  & (CII - 1) * memorypool%hash_size, RII)
             working_value = &
                  & memorypool%value_array(working_col, RII)
             memorypool%value_array(working_col, RII) = 0
             memorypool%dirty_array(working_col, RII) = .FALSE.
             !! If above threshold, insert
             IF (ABS(alpha*working_value) .GT. threshold) THEN
                memorypool%pruned_list(PII)%point_value = alpha*working_value
                memorypool%pruned_list(PII)%index_column = working_col
                memorypool%pruned_list(PII)%index_row = RII
                PII = PII + 1
             END IF
          END DO
       END DO
    END DO

    !! Convert to matrix
    CALL ConstructTripletList(unsorted_pruned_list, PII - 1)
    unsorted_pruned_list%DATA(:) = memorypool%pruned_list(1:PII - 1)
    CALL SortTripletList(unsorted_pruned_list, mat_c_columns, mat_c_rows, &
         & sorted_pruned_list, bubble_in = .TRUE.)
    CALL ConstructMatrixFromTripletList(matAB, sorted_pruned_list, mat_c_rows, &
         & mat_c_columns)
    CALL DestructTripletList(sorted_pruned_list)
    CALL DestructTripletList(unsorted_pruned_list)
  END SUBROUTINE PruneList_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Scale a matrix using a diagonal matrix (triplet list form).
  SUBROUTINE MatrixDiagonalScale_lsr(mat, tlist)
    !> The matrix to scale.
    TYPE(Matrix_lsr), INTENT(INOUT)  :: mat
    !> The diagonal matrix.
    TYPE(TripletList_r), INTENT(IN)  :: tlist
    !! Intermediate Data
    REAL(NTREAL) :: val



    INTEGER :: col, II

    DO II = 1, tlist%CurrentSize
       col = tlist%DATA(II)%index_column
       val = tlist%DATA(II)%point_value
       mat%values(mat%outer_index(col) + 1:mat%outer_index(col + 1)) = &
            val * mat%values(mat%outer_index(col) + 1:mat%outer_index(col + 1))
    END DO
  END SUBROUTINE MatrixDiagonalScale_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Scale a matrix using a diagonal matrix (triplet list form).
  SUBROUTINE MatrixDiagonalScale_lsc(mat, tlist)
    !> The matrix to scale.
    TYPE(Matrix_lsc), INTENT(INOUT)  :: mat
    !> The diagonal matrix.
    TYPE(TripletList_c), INTENT(IN)  :: tlist
    !! Intermediate Data
    COMPLEX(NTCOMPLEX) :: val



    INTEGER :: col, II

    DO II = 1, tlist%CurrentSize
       col = tlist%DATA(II)%index_column
       val = tlist%DATA(II)%point_value
       mat%values(mat%outer_index(col) + 1:mat%outer_index(col + 1)) = &
            val * mat%values(mat%outer_index(col) + 1:mat%outer_index(col + 1))
    END DO
  END SUBROUTINE MatrixDiagonalScale_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE SMatrixAlgebraModule
