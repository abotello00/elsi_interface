!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for handling locally stored CSR matrices.
MODULE SMatrixModule
  USE DataTypesModule, ONLY: NTREAL, NTCOMPLEX, NTLONG
  USE MatrixMarketModule, ONLY : ParseMMHeader, WriteMMSize, WriteMMLine, &
       & MAX_LINE_LENGTH
  USE TripletListModule, ONLY: TripletList_r, TripletList_c, SortTripletList, &
       & DestructTripletList, ConstructTripletList, AppendToTripletList, &
       & SymmetrizeTripletList, ConvertTripletListType
  USE TripletModule, ONLY : Triplet_r, Triplet_c
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A datatype for storing a local, real CSR matrix.
  TYPE, PUBLIC :: Matrix_lsr
     INTEGER, DIMENSION(:), ALLOCATABLE :: outer_index !< Outer indices
     INTEGER, DIMENSION(:), ALLOCATABLE :: inner_index !< Inner indices
     REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: values !< Values
     INTEGER :: rows !< Matrix dimension: rows
     INTEGER :: columns !< Matrix dimension: columns
  END TYPE Matrix_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A datatype for storing a local, complex CSR matrix.
  TYPE, PUBLIC :: Matrix_lsc
     INTEGER, DIMENSION(:), ALLOCATABLE :: outer_index !< Outer indices
     INTEGER, DIMENSION(:), ALLOCATABLE :: inner_index !< Inner indices
     COMPLEX(NTCOMPLEX), DIMENSION(:), ALLOCATABLE :: values !< Values
     INTEGER :: rows !< Matrix dimension: rows
     INTEGER :: columns !< Matrix dimension: columns
  END TYPE Matrix_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Construct/Destruct
  PUBLIC :: ConstructEmptyMatrix
  PUBLIC :: ConstructMatrixFromFile
  PUBLIC :: ConstructMatrixFromTripletList
  PUBLIC :: DestructMatrix
  PUBLIC :: CopyMatrix
  !! Basic Accessors
  PUBLIC :: GetMatrixRows
  PUBLIC :: GetMatrixColumns
  PUBLIC :: ExtractMatrixRow
  PUBLIC :: ExtractMatrixColumn
  !! Routines for splitting and composing
  PUBLIC :: SplitMatrix
  PUBLIC :: SplitMatrixColumns
  PUBLIC :: ComposeMatrix
  PUBLIC :: ComposeMatrixColumns
  !! ETC
  PUBLIC :: ConvertMatrixType
  PUBLIC :: TransposeMatrix
  PUBLIC :: ConjugateMatrix
  PUBLIC :: PrintMatrix
  PUBLIC :: MatrixToTripletList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE ConstructEmptyMatrix
     MODULE PROCEDURE ConstructEmptyMatrixSub_lsr
     MODULE PROCEDURE ConstructEmptyMatrixSub_lsc
  END INTERFACE ConstructEmptyMatrix
  INTERFACE ConstructMatrixFromFile
     MODULE PROCEDURE ConstructMatrixFromFileSub_lsr
     MODULE PROCEDURE ConstructMatrixFromFileSub_lsc
  END INTERFACE ConstructMatrixFromFile
  INTERFACE ConstructMatrixFromTripletList
     MODULE PROCEDURE ConstructMatrixFromTripletListSub_lsr
     MODULE PROCEDURE ConstructMatrixFromTripletListSub_lsc
  END INTERFACE ConstructMatrixFromTripletList
  INTERFACE DestructMatrix
     MODULE PROCEDURE DestructMatrix_lsr
     MODULE PROCEDURE DestructMatrix_lsc
  END INTERFACE DestructMatrix
  INTERFACE CopyMatrix
     MODULE PROCEDURE CopyMatrix_lsr
     MODULE PROCEDURE CopyMatrix_lsc
  END INTERFACE CopyMatrix
  INTERFACE GetMatrixRows
     MODULE PROCEDURE GetMatrixRows_lsr
     MODULE PROCEDURE GetMatrixRows_lsc
  END INTERFACE GetMatrixRows
  INTERFACE GetMatrixColumns
     MODULE PROCEDURE GetMatrixColumns_lsr
     MODULE PROCEDURE GetMatrixColumns_lsc
  END INTERFACE GetMatrixColumns
  INTERFACE ExtractMatrixRow
     MODULE PROCEDURE ExtractMatrixRow_lsr
     MODULE PROCEDURE ExtractMatrixRow_lsc
  END INTERFACE ExtractMatrixRow
  INTERFACE ExtractMatrixColumn
     MODULE PROCEDURE ExtractMatrixColumn_lsr
     MODULE PROCEDURE ExtractMatrixColumn_lsc
  END INTERFACE ExtractMatrixColumn
  INTERFACE SplitMatrix
     MODULE PROCEDURE SplitMatrix_lsr
     MODULE PROCEDURE SplitMatrix_lsc
  END INTERFACE SplitMatrix
  INTERFACE SplitMatrixColumns
     MODULE PROCEDURE SplitMatrixColumns_lsr
     MODULE PROCEDURE SplitMatrixColumns_lsc
  END INTERFACE SplitMatrixColumns
  INTERFACE ComposeMatrix
     MODULE PROCEDURE ComposeMatrix_lsr
     MODULE PROCEDURE ComposeMatrix_lsc
  END INTERFACE ComposeMatrix
  INTERFACE ComposeMatrixColumns
     MODULE PROCEDURE ComposeMatrixColumns_lsr
     MODULE PROCEDURE ComposeMatrixColumns_lsc
  END INTERFACE ComposeMatrixColumns
  INTERFACE TransposeMatrix
     MODULE PROCEDURE TransposeMatrix_lsr
     MODULE PROCEDURE TransposeMatrix_lsc
  END INTERFACE TransposeMatrix
  INTERFACE ConjugateMatrix
     MODULE PROCEDURE ConjugateMatrix_lsc
  END INTERFACE ConjugateMatrix
  INTERFACE PrintMatrix
     MODULE PROCEDURE PrintMatrix_lsr
     MODULE PROCEDURE PrintMatrix_lsc
  END INTERFACE PrintMatrix
  INTERFACE MatrixToTripletList
     MODULE PROCEDURE MatrixToTripletList_lsr
     MODULE PROCEDURE MatrixToTripletList_lsc
  END INTERFACE MatrixToTripletList
  INTERFACE ConvertMatrixType
     MODULE PROCEDURE ConvertMatrixType_lsrtolsc
     MODULE PROCEDURE ConvertMatrixType_lsctolsr
  END INTERFACE ConvertMatrixType
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A subroutine type wrapper for the constructor.
  PURE SUBROUTINE ConstructEmptyMatrixSub_lsr(this, rows, columns, zero_in)
    !> The matrix to construct.
    TYPE(Matrix_lsr), INTENT(INOUT) :: this
    !> The number of matrix columns.
    INTEGER, INTENT(IN) :: columns
    !> The number of matrix rows.
    INTEGER, INTENT(IN) :: rows
    !> Whether to set the matrix to zero.
    LOGICAL, INTENT(IN), OPTIONAL :: zero_in


    CALL DestructMatrix(this)

    this%rows = rows
    this%columns = columns
    ALLOCATE(this%outer_index(this%columns + 1))
    this%outer_index = 0

    IF (PRESENT(zero_in)) THEN
       IF (zero_in) THEN
          ALLOCATE(this%inner_index(0))
          ALLOCATE(this%values(0))
       END IF
    END IF
  END SUBROUTINE ConstructEmptyMatrixSub_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A subroutine type wrapper for the constructor.
  PURE SUBROUTINE ConstructEmptyMatrixSub_lsc(this, rows, columns, zero_in)
    !> The matrix to construct.
    TYPE(Matrix_lsc), INTENT(INOUT) :: this
    !> The number of matrix columns.
    INTEGER, INTENT(IN) :: columns
    !> The number of matrix rows.
    INTEGER, INTENT(IN) :: rows
    !> Whether to set the matrix to zero.
    LOGICAL, INTENT(IN), OPTIONAL :: zero_in


    CALL DestructMatrix(this)

    this%rows = rows
    this%columns = columns
    ALLOCATE(this%outer_index(this%columns + 1))
    this%outer_index = 0

    IF (PRESENT(zero_in)) THEN
       IF (zero_in) THEN
          ALLOCATE(this%inner_index(0))
          ALLOCATE(this%values(0))
       END IF
    END IF
  END SUBROUTINE ConstructEmptyMatrixSub_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Subroutine wrapper for the construct from file function.
  SUBROUTINE ConstructMatrixFromFileSub_lsr(this, file_name)
    !> The matrix being constructed.
    TYPE(Matrix_lsr), INTENT(INOUT) :: this
    !> Name of the file.
    CHARACTER(len=*), INTENT(IN) :: file_name
    !! About the matrix market file.
    INTEGER :: sparsity_type, data_type, pattern_type
    !! Local Data
    TYPE(TripletList_r) :: tlist
    TYPE(TripletList_r) :: sorted_tlist
    TYPE(Triplet_r) :: temp


    !! Local Data
    INTEGER :: rows, columns, total_values
    CHARACTER(LEN = 81) :: input_buffer
    INTEGER, PARAMETER :: file_handler = 16
    LOGICAL :: found_comment_line
    LOGICAL :: error_occured
    INTEGER :: II

    OPEN(file_handler, file = file_name, status='old')

    !! Parse the header.
    READ(file_handler,fmt='(A)') input_buffer
    error_occured = ParseMMHeader(input_buffer, sparsity_type, data_type, &
         & pattern_type)

    !! Extra Comment Lines
    found_comment_line = .TRUE.
    DO WHILE(found_comment_line)
       READ(file_handler, fmt = '(A)') input_buffer
       IF (.NOT. input_buffer(1:1) .EQ. '%') THEN
          found_comment_line = .FALSE.
       END IF
    END DO

    !! Main data
    READ(input_buffer,*) rows, columns, total_values
    CALL ConstructTripletList(tlist)

    !! Read Values
    DO II = 1, total_values




       READ(file_handler,*) temp%index_row, temp%index_column, temp%point_value

       CALL AppendToTripletList(tlist, temp)
    END DO

    CLOSE(file_handler)
    CALL SymmetrizeTripletList(tlist, pattern_type)
    CALL SortTripletList(tlist, columns, rows, sorted_tlist)
    CALL ConstructMatrixFromTripletList(this, sorted_tlist, rows, columns)

    CALL DestructTripletList(tlist)
    CALL DestructTripletList(sorted_tlist)
  END SUBROUTINE ConstructMatrixFromFileSub_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE ConstructMatrixFromFileSub_lsc(this, file_name)
    !> The matrix being constructed.
    TYPE(Matrix_lsc), INTENT(INOUT) :: this
    !> Name of the file.
    CHARACTER(len=*), INTENT(IN) :: file_name
    !! About the matrix market file.
    INTEGER :: sparsity_type, data_type, pattern_type
    !! Local Data
    TYPE(TripletList_c) :: tlist
    TYPE(TripletList_c) :: sorted_tlist
    TYPE(Triplet_c) :: temp
    REAL(NTREAL) :: real_val, comp_val



    !! Local Data
    INTEGER :: rows, columns, total_values
    CHARACTER(LEN = 81) :: input_buffer
    INTEGER, PARAMETER :: file_handler = 16
    LOGICAL :: found_comment_line
    LOGICAL :: error_occured
    INTEGER :: II

    OPEN(file_handler, file = file_name, status='old')

    !! Parse the header.
    READ(file_handler,fmt='(A)') input_buffer
    error_occured = ParseMMHeader(input_buffer, sparsity_type, data_type, &
         & pattern_type)

    !! Extra Comment Lines
    found_comment_line = .TRUE.
    DO WHILE(found_comment_line)
       READ(file_handler, fmt = '(A)') input_buffer
       IF (.NOT. input_buffer(1:1) .EQ. '%') THEN
          found_comment_line = .FALSE.
       END IF
    END DO

    !! Main data
    READ(input_buffer,*) rows, columns, total_values
    CALL ConstructTripletList(tlist)

    !! Read Values
    DO II = 1, total_values

       READ(file_handler,*) temp%index_row, temp%index_column, real_val, comp_val
       temp%point_value = CMPLX(real_val, comp_val, KIND = NTCOMPLEX)



       CALL AppendToTripletList(tlist, temp)
    END DO

    CLOSE(file_handler)
    CALL SymmetrizeTripletList(tlist, pattern_type)
    CALL SortTripletList(tlist, columns, rows, sorted_tlist)
    CALL ConstructMatrixFromTripletList(this, sorted_tlist, rows, columns)

    CALL DestructTripletList(tlist)
    CALL DestructTripletList(sorted_tlist)

  END SUBROUTINE ConstructMatrixFromFileSub_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A subroutine wrapper for the triplet list based constructor.
  PURE SUBROUTINE ConstructMatrixFromTripletListSub_lsr(this, triplet_list, &
       & rows, columns)
    !> The matrix being constructed
    TYPE(Matrix_lsr), INTENT(INOUT) :: this
    !> A list of triplet values. They must be sorted.
    TYPE(TripletList_r), INTENT(IN) :: triplet_list
    !> Number of matrix rows
    INTEGER, INTENT(IN) :: rows
    !> Number of matrix columns
    INTEGER, INTENT(IN) :: columns



    !! Local Data
    INTEGER :: II, outer_idx

    CALL DestructMatrix(this)

    this%rows = rows
    this%columns = columns

    !! Allocate
    ALLOCATE(this%outer_index(this%columns + 1))
    this%outer_index = 0
    ALLOCATE(this%inner_index(triplet_list%CurrentSize))
    ALLOCATE(this%values(triplet_list%CurrentSize))

    !! Insert Values
    outer_idx = 1
    DO II = 1, triplet_list%CurrentSize
       !! Moving on to the next column?
       DO WHILE(.NOT. triplet_list%DATA(II)%index_column .EQ. outer_idx)
          outer_idx = outer_idx + 1
          this%outer_index(outer_idx + 1) = this%outer_index(outer_idx)
       END DO
       this%outer_index(outer_idx + 1)=this%outer_index(outer_idx + 1) + 1
       !! Insert inner index and value
       this%inner_index(II) = triplet_list%DATA(II)%index_row
       this%values(II) = triplet_list%DATA(II)%point_value
    END DO

    !! Fill In The Rest Of The Outer Values
    DO outer_idx = outer_idx + 2, this%columns + 1
       this%outer_index(outer_idx) = this%outer_index(outer_idx - 1)
    END DO

  END SUBROUTINE ConstructMatrixFromTripletListSub_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A subroutine wrapper for the triplet list based constructor.
  PURE SUBROUTINE ConstructMatrixFromTripletListSub_lsc(this, triplet_list, &
       & rows, columns)
    !> The matrix being constructed
    TYPE(Matrix_lsc), INTENT(INOUT) :: this
    !> A list of triplet values. They must be sorted.
    TYPE(TripletList_c), INTENT(IN) :: triplet_list
    !> Number of matrix rows
    INTEGER, INTENT(IN) :: rows
    !> Number of matrix columns
    INTEGER, INTENT(IN) :: columns


    !! Local Data
    INTEGER :: II, outer_idx

    CALL DestructMatrix(this)

    this%rows = rows
    this%columns = columns

    !! Allocate
    ALLOCATE(this%outer_index(this%columns + 1))
    this%outer_index = 0
    ALLOCATE(this%inner_index(triplet_list%CurrentSize))
    ALLOCATE(this%values(triplet_list%CurrentSize))

    !! Insert Values
    outer_idx = 1
    DO II = 1, triplet_list%CurrentSize
       !! Moving on to the next column?
       DO WHILE(.NOT. triplet_list%DATA(II)%index_column .EQ. outer_idx)
          outer_idx = outer_idx + 1
          this%outer_index(outer_idx + 1) = this%outer_index(outer_idx)
       END DO
       this%outer_index(outer_idx + 1)=this%outer_index(outer_idx + 1) + 1
       !! Insert inner index and value
       this%inner_index(II) = triplet_list%DATA(II)%index_row
       this%values(II) = triplet_list%DATA(II)%point_value
    END DO

    !! Fill In The Rest Of The Outer Values
    DO outer_idx = outer_idx + 2, this%columns + 1
       this%outer_index(outer_idx) = this%outer_index(outer_idx - 1)
    END DO
  END SUBROUTINE ConstructMatrixFromTripletListSub_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Explicitly destruct a sparse matrix.
  PURE SUBROUTINE DestructMatrix_lsr(this)
    !> The matrix to free up.
    TYPE(Matrix_lsr), INTENT(INOUT) :: this


    IF (ALLOCATED(this%outer_index)) THEN
       DEALLOCATE(this%outer_index)
    END IF
    IF (ALLOCATED(this%inner_index)) THEN
       DEALLOCATE(this%inner_index)
    END IF
    IF (ALLOCATED(this%values)) THEN
       DEALLOCATE(this%values)
    END IF
  END SUBROUTINE DestructMatrix_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Explicitly destruct a sparse matrix.
  PURE SUBROUTINE DestructMatrix_lsc(this)
    !> The matrix to free up.
    TYPE(Matrix_lsc), INTENT(INOUT) :: this


    IF (ALLOCATED(this%outer_index)) THEN
       DEALLOCATE(this%outer_index)
    END IF
    IF (ALLOCATED(this%inner_index)) THEN
       DEALLOCATE(this%inner_index)
    END IF
    IF (ALLOCATED(this%values)) THEN
       DEALLOCATE(this%values)
    END IF
  END SUBROUTINE DestructMatrix_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Copy a sparse matrix in a safe way.
  PURE SUBROUTINE CopyMatrix_lsr(matA, matB)
    !> Matrix to copy
    TYPE(Matrix_lsr), INTENT(IN) :: matA
    !> matB = matA
    TYPE(Matrix_lsr), INTENT(INOUT) :: matB


    CALL DestructMatrix(matB)
    matB = matA
  END SUBROUTINE CopyMatrix_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Copy a sparse matrix in a safe way.
  PURE SUBROUTINE CopyMatrix_lsc(matA, matB)
    !> Matrix to copy
    TYPE(Matrix_lsc), INTENT(IN) :: matA
    !> matB = matA
    TYPE(Matrix_lsc), INTENT(INOUT) :: matB


    CALL DestructMatrix(matB)
    matB = matA
  END SUBROUTINE CopyMatrix_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the number of rows of a matrix.
  PURE FUNCTION GetMatrixRows_lsr(this) RESULT(rows)
    !> The matrix.
    TYPE(Matrix_lsr), INTENT(IN) :: this
    !> The number of rows.
    INTEGER :: rows


    rows = this%rows
  END FUNCTION GetMatrixRows_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the number of rows of a matrix.
  PURE FUNCTION GetMatrixRows_lsc(this) RESULT(rows)
    !> The matrix.
    TYPE(Matrix_lsc), INTENT(IN) :: this
    !> The number of rows.
    INTEGER :: rows


    rows = this%rows
  END FUNCTION GetMatrixRows_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the number of columns of a matrix.
  PURE FUNCTION GetMatrixColumns_lsr(this) RESULT(columns)
    !! The matrix.
    TYPE(Matrix_lsr), INTENT(IN) :: this
    !> The number of columns.
    INTEGER :: columns


    columns = this%columns
  END FUNCTION GetMatrixColumns_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the number of columns of a matrix.
  PURE FUNCTION GetMatrixColumns_lsc(this) RESULT(columns)
    !! The matrix.
    TYPE(Matrix_lsc), INTENT(IN) :: this
    !> The number of columns.
    INTEGER :: columns


    columns = this%columns
  END FUNCTION GetMatrixColumns_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Extract a row from the matrix.
  PURE SUBROUTINE ExtractMatrixRow_lsr(this, row_number, row_out)
    !> The matrix to extract from.
    TYPE(Matrix_lsr), INTENT(IN) :: this
    !> The row to extract.
    INTEGER, INTENT(IN) :: row_number
    !> The matrix representing that row.
    TYPE(Matrix_lsr), INTENT(INOUT) :: row_out
    !! Temporary Variables
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: value_buffer


    !! Temporary Variables
    INTEGER :: values_found
    INTEGER :: elements_per_inner
    INTEGER :: II, JJ, KK

    !! Fill a value buffer
    CALL ConstructEmptyMatrix(row_out, 1, this%columns)
    ALLOCATE(value_buffer(this%columns))
    values_found = 0
    KK = 1
    row_out%outer_index(1) = 0
    DO II = 1, this%columns
       row_out%outer_index(II + 1) = row_out%outer_index(II + 1) + &
            & row_out%outer_index(II)
       elements_per_inner = this%outer_index(II + 1) - this%outer_index(II)
       DO JJ = 1, elements_per_inner
          IF (this%inner_index(KK) .EQ. row_number) THEN
             values_found = values_found + 1
             value_buffer(values_found) = this%values(KK)
             row_out%outer_index(II + 1) = row_out%outer_index(II + 1) + 1
          END IF
          KK = KK + 1
       END DO
    END DO

    !! Copy To Actual Matrix
    ALLOCATE(row_out%inner_index(values_found))
    row_out%inner_index = 1
    ALLOCATE(row_out%values(values_found))
    row_out%values(:) = value_buffer(:values_found)

    !! Cleanup
    DEALLOCATE(value_buffer)
  END SUBROUTINE ExtractMatrixRow_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Extract a row from the matrix.
  PURE SUBROUTINE ExtractMatrixRow_lsc(this, row_number, row_out)
    !> The matrix to extract from.
    TYPE(Matrix_lsc), INTENT(IN) :: this
    !> The row to extract.
    INTEGER, INTENT(IN) :: row_number
    !> The matrix representing that row.
    TYPE(Matrix_lsc), INTENT(INOUT) :: row_out
    !! Temporary Variables
    COMPLEX(NTCOMPLEX), DIMENSION(:), ALLOCATABLE :: value_buffer


    !! Temporary Variables
    INTEGER :: values_found
    INTEGER :: elements_per_inner
    INTEGER :: II, JJ, KK

    !! Fill a value buffer
    CALL ConstructEmptyMatrix(row_out, 1, this%columns)
    ALLOCATE(value_buffer(this%columns))
    values_found = 0
    KK = 1
    row_out%outer_index(1) = 0
    DO II = 1, this%columns
       row_out%outer_index(II + 1) = row_out%outer_index(II + 1) + &
            & row_out%outer_index(II)
       elements_per_inner = this%outer_index(II + 1) - this%outer_index(II)
       DO JJ = 1, elements_per_inner
          IF (this%inner_index(KK) .EQ. row_number) THEN
             values_found = values_found + 1
             value_buffer(values_found) = this%values(KK)
             row_out%outer_index(II + 1) = row_out%outer_index(II + 1) + 1
          END IF
          KK = KK + 1
       END DO
    END DO

    !! Copy To Actual Matrix
    ALLOCATE(row_out%inner_index(values_found))
    row_out%inner_index = 1
    ALLOCATE(row_out%values(values_found))
    row_out%values(:) = value_buffer(:values_found)

    !! Cleanup
    DEALLOCATE(value_buffer)
  END SUBROUTINE ExtractMatrixRow_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Extract a column from the matrix.
  PURE SUBROUTINE ExtractMatrixColumn_lsr(this, column_number, column_out)
    !> The matrix to extract from.
    TYPE(Matrix_lsr), INTENT(IN) :: this
    !> The column to extract.
    INTEGER, INTENT(IN) :: column_number
    !> The column representing that row.
    TYPE(Matrix_lsr), INTENT(INOUT) :: column_out


    !! Local variables
    INTEGER :: number_of_values
    INTEGER :: start_index
    INTEGER :: II

    !! Allocate Memory
    CALL ConstructEmptyMatrix(column_out, this%rows, 1)
    start_index = this%outer_index(column_number)
    number_of_values = this%outer_index(column_number + 1) - &
         & this%outer_index(column_number)
    ALLOCATE(column_out%inner_index(number_of_values))
    ALLOCATE(column_out%values(number_of_values))

    !! Copy Values
    column_out%outer_index(1) = 0
    column_out%outer_index(2) = number_of_values
    DO II = 1, number_of_values
       column_out%inner_index(II) = this%inner_index(start_index + II)
       column_out%values(II) = this%values(start_index + II)
    END DO
  END SUBROUTINE ExtractMatrixColumn_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Extract a column from the matrix.
  PURE SUBROUTINE ExtractMatrixColumn_lsc(this, column_number, column_out)
    !> The matrix to extract from.
    TYPE(Matrix_lsc), INTENT(IN) :: this
    !> The column to extract.
    INTEGER, INTENT(IN) :: column_number
    !> The column representing that row.
    TYPE(Matrix_lsc), INTENT(INOUT) :: column_out


    !! Local variables
    INTEGER :: number_of_values
    INTEGER :: start_index
    INTEGER :: II

    !! Allocate Memory
    CALL ConstructEmptyMatrix(column_out, this%rows, 1)
    start_index = this%outer_index(column_number)
    number_of_values = this%outer_index(column_number + 1) - &
         & this%outer_index(column_number)
    ALLOCATE(column_out%inner_index(number_of_values))
    ALLOCATE(column_out%values(number_of_values))

    !! Copy Values
    column_out%outer_index(1) = 0
    column_out%outer_index(2) = number_of_values
    DO II = 1, number_of_values
       column_out%inner_index(II) = this%inner_index(start_index + II)
       column_out%values(II) = this%values(start_index + II)
    END DO
  END SUBROUTINE ExtractMatrixColumn_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Transpose a sparse matrix and return it in a separate matrix.
  !> The current implementation has you go from matrix to triplet list,
  !> triplet list to transposed triplet list. The triplet list must then be
  !> sorted and then the return matrix is constructed.
  PURE SUBROUTINE TransposeMatrix_lsr(this, matT)
    !> The matrix to be transposed.
    TYPE(Matrix_lsr), INTENT(IN)  :: this
    !> The input matrix transposed.
    TYPE(Matrix_lsr), INTENT(INOUT) :: matT


    !! Local Data
    INTEGER, DIMENSION(:), ALLOCATABLE :: values_per_row
    INTEGER, DIMENSION(:), ALLOCATABLE :: offset_array
    !! Temporary Variables
    INTEGER :: II, JJ
    INTEGER :: inner_index, insert_pt, this_offset
    INTEGER :: num_values, elements_per_inner

    !! Allocate New Matrix
    num_values = this%outer_index(this%columns + 1)
    CALL ConstructEmptyMatrix(matT, this%columns, this%rows)
    ALLOCATE(matT%inner_index(num_values))
    ALLOCATE(matT%values(num_values))

    !! Temporary Arrays
    ALLOCATE(values_per_row(this%rows))
    ALLOCATE(offset_array(this%rows))

    !! Count the values per row
    values_per_row = 0
    DO II = 1, num_values
       inner_index = this%inner_index(II)
       values_per_row(inner_index) = values_per_row(inner_index) + 1
    END DO
    offset_array(1) = 0
    DO II = 2, this%rows
       offset_array(II) = offset_array(II - 1) + values_per_row(II - 1)
    END DO

    !! Insert
    matT%outer_index(:this%rows) = offset_array(:this%rows)
    matT%outer_index(this%rows + 1) = offset_array(this%rows) + &
         & values_per_row(this%rows)
    DO II = 1, this%columns
       elements_per_inner = this%outer_index(II + 1) - this%outer_index(II)
       this_offset = this%outer_index(II)
       DO JJ = 1, elements_per_inner
          inner_index = this%inner_index(this_offset + JJ)
          insert_pt = offset_array(inner_index)+1
          matT%inner_index(insert_pt) = II
          matT%values(insert_pt) = this%values(this_offset + JJ)
          offset_array(inner_index) = offset_array(inner_index) +1
       END DO
    END DO

    !! Cleanup
    DEALLOCATE(values_per_row)
    DEALLOCATE(offset_array)
  END SUBROUTINE TransposeMatrix_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Transpose a sparse matrix and return it in a separate matrix.
  !> The current implementation has you go from matrix to triplet list,
  !> triplet list to transposed triplet list. The triplet list must then be
  !> sorted and then the return matrix is constructed.
  PURE SUBROUTINE TransposeMatrix_lsc(this, matT)
    !> The matrix to be transposed.
    TYPE(Matrix_lsc), INTENT(IN)  :: this
    !> The input matrix transposed.
    TYPE(Matrix_lsc), INTENT(INOUT) :: matT


    !! Local Data
    INTEGER, DIMENSION(:), ALLOCATABLE :: values_per_row
    INTEGER, DIMENSION(:), ALLOCATABLE :: offset_array
    !! Temporary Variables
    INTEGER :: II, JJ
    INTEGER :: inner_index, insert_pt, this_offset
    INTEGER :: num_values, elements_per_inner

    !! Allocate New Matrix
    num_values = this%outer_index(this%columns + 1)
    CALL ConstructEmptyMatrix(matT, this%columns, this%rows)
    ALLOCATE(matT%inner_index(num_values))
    ALLOCATE(matT%values(num_values))

    !! Temporary Arrays
    ALLOCATE(values_per_row(this%rows))
    ALLOCATE(offset_array(this%rows))

    !! Count the values per row
    values_per_row = 0
    DO II = 1, num_values
       inner_index = this%inner_index(II)
       values_per_row(inner_index) = values_per_row(inner_index) + 1
    END DO
    offset_array(1) = 0
    DO II = 2, this%rows
       offset_array(II) = offset_array(II - 1) + values_per_row(II - 1)
    END DO

    !! Insert
    matT%outer_index(:this%rows) = offset_array(:this%rows)
    matT%outer_index(this%rows + 1) = offset_array(this%rows) + &
         & values_per_row(this%rows)
    DO II = 1, this%columns
       elements_per_inner = this%outer_index(II + 1) - this%outer_index(II)
       this_offset = this%outer_index(II)
       DO JJ = 1, elements_per_inner
          inner_index = this%inner_index(this_offset + JJ)
          insert_pt = offset_array(inner_index)+1
          matT%inner_index(insert_pt) = II
          matT%values(insert_pt) = this%values(this_offset + JJ)
          offset_array(inner_index) = offset_array(inner_index) +1
       END DO
    END DO

    !! Cleanup
    DEALLOCATE(values_per_row)
    DEALLOCATE(offset_array)
  END SUBROUTINE TransposeMatrix_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Create a big matrix from an array of matrices by putting them one next
  !> to another.
  PURE SUBROUTINE ComposeMatrix_lsr(mat_array, block_rows, block_columns, &
       & out_matrix)
    !> The number of rows of the array of blocks.
    INTEGER, INTENT(IN) :: block_rows
    !> The number of columns of the array of blocks.
    INTEGER, INTENT(IN) :: block_columns
    !> 2d array of matrices to compose.
    TYPE(Matrix_lsr), DIMENSION(:, :), INTENT(IN) :: mat_array
    !> The composed matrix.
    TYPE(Matrix_lsr), INTENT(INOUT) :: out_matrix
    !! Local Data
    TYPE(Matrix_lsr), DIMENSION(block_columns) :: merged_columns
    TYPE(Matrix_lsr) :: Temp
    TYPE(Matrix_lsr), DIMENSION(block_rows, block_columns) :: mat_t


    INTEGER :: II, JJ

    !! First transpose the matrices
    DO JJ = 1, block_columns
       DO II = 1, block_rows
          CALL TransposeMatrix(mat_array(II, JJ), mat_t(II, JJ))
       END DO
    END DO

    !! Next merge the columns
    DO JJ = 1, block_columns
       CALL ComposeMatrixColumns(mat_t(:, JJ), Temp)
       CALL TransposeMatrix(Temp, merged_columns(JJ))
    END DO

    !! Final Merge
    CALL ComposeMatrixColumns(merged_columns, out_matrix)

    !! Cleanup
    DO JJ = 1, block_columns
       DO II = 1, block_rows
          CALL DestructMatrix(mat_t(II, JJ))
       END DO
    END DO
    DO JJ = 1, block_columns
       CALL DestructMatrix(merged_columns(JJ))
    END DO
  END SUBROUTINE ComposeMatrix_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Create a big matrix from an array of matrices by putting them one next
  !> to another.
  PURE SUBROUTINE ComposeMatrix_lsc(mat_array, block_rows, block_columns, &
       & out_matrix)
    !> The number of rows of the array of blocks.
    INTEGER, INTENT(IN) :: block_rows
    !> The number of columns of the array of blocks.
    INTEGER, INTENT(IN) :: block_columns
    !> 2d array of matrices to compose.
    TYPE(Matrix_lsc), DIMENSION(:, :), INTENT(IN) :: mat_array
    !> The composed matrix.
    TYPE(Matrix_lsc), INTENT(INOUT) :: out_matrix
    !! Local Data
    TYPE(Matrix_lsc), DIMENSION(block_columns) :: merged_columns
    TYPE(Matrix_lsc) :: Temp
    TYPE(Matrix_lsc), DIMENSION(block_rows, block_columns) :: mat_t


    INTEGER :: II, JJ

    !! First transpose the matrices
    DO JJ = 1, block_columns
       DO II = 1, block_rows
          CALL TransposeMatrix(mat_array(II, JJ), mat_t(II, JJ))
       END DO
    END DO

    !! Next merge the columns
    DO JJ = 1, block_columns
       CALL ComposeMatrixColumns(mat_t(:, JJ), Temp)
       CALL TransposeMatrix(Temp, merged_columns(JJ))
    END DO

    !! Final Merge
    CALL ComposeMatrixColumns(merged_columns, out_matrix)

    !! Cleanup
    DO JJ = 1, block_columns
       DO II = 1, block_rows
          CALL DestructMatrix(mat_t(II, JJ))
       END DO
    END DO
    DO JJ = 1, block_columns
       CALL DestructMatrix(merged_columns(JJ))
    END DO
  END SUBROUTINE ComposeMatrix_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Create a big Matrix C = [Matrix 1 | Matrix 1, ...] where the columns of
  !> the first matrix are followed by the columns of the matrices in the list.
  PURE SUBROUTINE ComposeMatrixColumns_lsr(mat_list, out_matrix)
    !> A list of matrices to compose.
    TYPE(Matrix_lsr), DIMENSION(:), INTENT(IN) :: mat_list
    !> out_matrix = [Matrix 1 | Matrix 2, ...].
    TYPE(Matrix_lsr), INTENT(INOUT) :: out_matrix


    !! Local Variables
    INTEGER :: total_columns, total_values
    INTEGER :: inner_start, inner_length
    INTEGER :: outer_start, outer_length
    INTEGER :: outer_offset
    INTEGER :: II
    INTEGER :: size_of_mat

    CALL DestructMatrix(out_matrix)

    !! Figure Out The Sizes
    total_columns = 0
    total_values  = 0
    DO II = LBOUND(mat_list, dim = 1), UBOUND(mat_list, dim = 1)
       total_columns = total_columns + mat_list(II)%columns
       size_of_mat = mat_list(II)%outer_index(mat_list(II)%columns + 1)
       total_values  = total_values + size_of_mat
    END DO

    !! Allocate The Space
    CALL ConstructEmptyMatrix(out_matrix, mat_list(1)%rows, total_columns)
    ALLOCATE(out_matrix%inner_index(total_values))
    ALLOCATE(out_matrix%values(total_values))

    !! Fill In The Values
    inner_start = 1
    outer_start = 1
    outer_offset = 0
    DO II = LBOUND(mat_list, dim = 1), UBOUND(mat_list, dim = 1)
       !! Inner indices and values
       size_of_mat = mat_list(II)%outer_index(mat_list(II)%columns + 1)
       inner_length = size_of_mat
       out_matrix%inner_index(inner_start:inner_start + inner_length - 1) = &
            & mat_list(II)%inner_index
       out_matrix%values(inner_start:inner_start + inner_length - 1) = &
            & mat_list(II)%values
       inner_start = inner_start + inner_length
       !! Outer Indices
       outer_length = mat_list(II)%columns + 1
       out_matrix%outer_index(outer_start:outer_start+outer_length - 1) = &
            & mat_list(II)%outer_index + outer_offset
       outer_start = outer_start + outer_length - 1
       outer_offset = out_matrix%outer_index(outer_start)
    END DO
  END SUBROUTINE ComposeMatrixColumns_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Create a big Matrix C = [Matrix 1 | Matrix 1, ...] where the columns of
  !> the first matrix are followed by the columns of the matrices in the list.
  PURE SUBROUTINE ComposeMatrixColumns_lsc(mat_list, out_matrix)
    !> A list of matrices to compose.
    TYPE(Matrix_lsc), DIMENSION(:), INTENT(IN) :: mat_list
    !> out_matrix = [Matrix 1 | Matrix 2, ...].
    TYPE(Matrix_lsc), INTENT(INOUT) :: out_matrix


    !! Local Variables
    INTEGER :: total_columns, total_values
    INTEGER :: inner_start, inner_length
    INTEGER :: outer_start, outer_length
    INTEGER :: outer_offset
    INTEGER :: II
    INTEGER :: size_of_mat

    CALL DestructMatrix(out_matrix)

    !! Figure Out The Sizes
    total_columns = 0
    total_values  = 0
    DO II = LBOUND(mat_list, dim = 1), UBOUND(mat_list, dim = 1)
       total_columns = total_columns + mat_list(II)%columns
       size_of_mat = mat_list(II)%outer_index(mat_list(II)%columns + 1)
       total_values  = total_values + size_of_mat
    END DO

    !! Allocate The Space
    CALL ConstructEmptyMatrix(out_matrix, mat_list(1)%rows, total_columns)
    ALLOCATE(out_matrix%inner_index(total_values))
    ALLOCATE(out_matrix%values(total_values))

    !! Fill In The Values
    inner_start = 1
    outer_start = 1
    outer_offset = 0
    DO II = LBOUND(mat_list, dim = 1), UBOUND(mat_list, dim = 1)
       !! Inner indices and values
       size_of_mat = mat_list(II)%outer_index(mat_list(II)%columns + 1)
       inner_length = size_of_mat
       out_matrix%inner_index(inner_start:inner_start + inner_length - 1) = &
            & mat_list(II)%inner_index
       out_matrix%values(inner_start:inner_start + inner_length - 1) = &
            & mat_list(II)%values
       inner_start = inner_start + inner_length
       !! Outer Indices
       outer_length = mat_list(II)%columns + 1
       out_matrix%outer_index(outer_start:outer_start+outer_length - 1) = &
            & mat_list(II)%outer_index + outer_offset
       outer_start = outer_start + outer_length - 1
       outer_offset = out_matrix%outer_index(outer_start)
    END DO
  END SUBROUTINE ComposeMatrixColumns_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Split a sparse matrix into an array of sparse matrices.
  PURE SUBROUTINE SplitMatrix_lsr(this, block_rows, block_columns, &
       & split_array, block_size_row_in, block_size_column_in)
    !> The matrix to split.
    TYPE(Matrix_lsr), INTENT(IN) :: this
    !> Number of rows to split the matrix into.
    INTEGER, INTENT(IN) :: block_rows
    !> Number of columns to split the matrix into.
    INTEGER, INTENT(IN) :: block_columns
    !> A COLUMNxROW array for the output to go into.
    TYPE(Matrix_lsr), DIMENSION(:,:), INTENT(INOUT) :: split_array
    !> Specifies the size of the  rows.
    INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: block_size_row_in
    !> Specifies the size of the columns.
    INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: block_size_column_in
    !! Local Data
    TYPE(Matrix_lsr), DIMENSION(block_columns) :: column_split
    TYPE(Matrix_lsr), DIMENSION(block_rows) :: row_split
    TYPE(Matrix_lsr) :: Temp


    !! Local Data
    INTEGER, DIMENSION(block_rows) :: block_size_row
    INTEGER, DIMENSION(block_columns) :: block_size_column
    !! Temporary Variables
    INTEGER :: divisor_row, divisor_column
    INTEGER :: II, JJ

    !! Calculate the split sizes
    IF (PRESENT(block_size_row_in)) THEN
       block_size_row = block_size_row_in
    ELSE
       divisor_row = this%rows / block_rows
       block_size_row = divisor_row
       block_size_row(block_rows) = this%rows - divisor_row*(block_rows - 1)
    END IF
    IF (PRESENT(block_size_column_in)) THEN
       block_size_column = block_size_column_in
    ELSE
       divisor_column = this%columns / block_columns
       block_size_column = divisor_column
       block_size_column(block_columns) = this%columns - &
            & divisor_column*(block_columns - 1)
    END IF

    !! First split by columns which is easy with the CSR format
    CALL SplitMatrixColumns(this, block_columns, block_size_column, &
         & column_split)

    !! Now Split By Rows
    DO JJ = 1, block_columns
       CALL TransposeMatrix(column_split(JJ), Temp)
       CALL SplitMatrixColumns(Temp, block_rows, block_size_row, &
            & row_split)
       !! Copy into output array
       DO II = 1, block_rows
          CALL TransposeMatrix(row_split(II), split_array(II, JJ))
       END DO
    END DO

    !! Cleanup
    CALL DestructMatrix(Temp)
    DO II = 1, block_rows
       CALL DestructMatrix(row_split(II))
    END DO
    DO II = 1, block_columns
       CALL DestructMatrix(column_split(II))
    END DO
  END SUBROUTINE SplitMatrix_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Split a sparse matrix into an array of sparse matrices.
  PURE SUBROUTINE SplitMatrix_lsc(this, block_rows, block_columns, &
       & split_array, block_size_row_in, block_size_column_in)
    !> The matrix to split.
    TYPE(Matrix_lsc), INTENT(IN) :: this
    !> Number of rows to split the matrix into.
    INTEGER, INTENT(IN) :: block_rows
    !> Number of columns to split the matrix into.
    INTEGER, INTENT(IN) :: block_columns
    !> A COLUMNxROW array for the output to go into.
    TYPE(Matrix_lsc), DIMENSION(:,:), INTENT(INOUT) :: split_array
    !> Specifies the size of the  rows.
    INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: block_size_row_in
    !> Specifies the size of the columns.
    INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: block_size_column_in
    !! Local Data
    TYPE(Matrix_lsc), DIMENSION(block_columns) :: column_split
    TYPE(Matrix_lsc), DIMENSION(block_rows) :: row_split
    TYPE(Matrix_lsc) :: Temp


    !! Local Data
    INTEGER, DIMENSION(block_rows) :: block_size_row
    INTEGER, DIMENSION(block_columns) :: block_size_column
    !! Temporary Variables
    INTEGER :: divisor_row, divisor_column
    INTEGER :: II, JJ

    !! Calculate the split sizes
    IF (PRESENT(block_size_row_in)) THEN
       block_size_row = block_size_row_in
    ELSE
       divisor_row = this%rows / block_rows
       block_size_row = divisor_row
       block_size_row(block_rows) = this%rows - divisor_row*(block_rows - 1)
    END IF
    IF (PRESENT(block_size_column_in)) THEN
       block_size_column = block_size_column_in
    ELSE
       divisor_column = this%columns / block_columns
       block_size_column = divisor_column
       block_size_column(block_columns) = this%columns - &
            & divisor_column*(block_columns - 1)
    END IF

    !! First split by columns which is easy with the CSR format
    CALL SplitMatrixColumns(this, block_columns, block_size_column, &
         & column_split)

    !! Now Split By Rows
    DO JJ = 1, block_columns
       CALL TransposeMatrix(column_split(JJ), Temp)
       CALL SplitMatrixColumns(Temp, block_rows, block_size_row, &
            & row_split)
       !! Copy into output array
       DO II = 1, block_rows
          CALL TransposeMatrix(row_split(II), split_array(II, JJ))
       END DO
    END DO

    !! Cleanup
    CALL DestructMatrix(Temp)
    DO II = 1, block_rows
       CALL DestructMatrix(row_split(II))
    END DO
    DO II = 1, block_columns
       CALL DestructMatrix(column_split(II))
    END DO
  END SUBROUTINE SplitMatrix_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Split a matrix into into small blocks based on the specified offsets.
  PURE SUBROUTINE SplitMatrixColumns_lsr(this, num_blocks, block_sizes, &
       & split_list)
    !> This matrix to perform this operation on.
    TYPE(Matrix_lsr), INTENT(IN) :: this
    !> Number of blocks to split into.
    INTEGER, INTENT(IN) :: num_blocks
    !> The sizes used for splitting.
    INTEGER, DIMENSION(num_blocks), INTENT(IN) :: block_sizes
    !> 1D array of blocks.
    TYPE(Matrix_lsr), DIMENSION(num_blocks), INTENT(INOUT) :: split_list


    !! Local Data
    INTEGER, DIMENSION(num_blocks+1) :: block_offsets
    !! Temporary variables
    INTEGER :: II, loffset, lcolumns, linner_offset, total_values

    !! Compute Offsets
    block_offsets(1) = 1
    DO II = 2, num_blocks + 1
       block_offsets(II) = block_offsets(II - 1) + block_sizes(II - 1)
    END DO

    !! Split up the columns
    DO II = 1, num_blocks
       !! Temporary variables
       loffset = block_offsets(II)
       lcolumns = block_sizes(II)
       linner_offset = this%outer_index(loffset) + 1

       !! Construct
       CALL ConstructEmptyMatrix(split_list(II), this%rows, lcolumns)

       !! Copy Outer Index
       split_list(II)%outer_index(:) = &
            & this%outer_index(loffset:loffset + lcolumns)
       split_list(II)%outer_index(:) = split_list(II)%outer_index(:) - &
            & split_list(II)%outer_index(1)
       total_values = split_list(II)%outer_index(lcolumns + 1)

       !! Copy Inner Indices and Values
       IF (total_values .GT. 0) THEN
          ALLOCATE(split_list(II)%inner_index(total_values))
          split_list(II)%inner_index(:) = &
               & this%inner_index(linner_offset:linner_offset + total_values - 1)
          ALLOCATE(split_list(II)%values(total_values))
          split_list(II)%values(:) = &
               & this%values(linner_offset:linner_offset + total_values - 1)
       END IF
    END DO
  END SUBROUTINE SplitMatrixColumns_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Split a matrix into into small blocks based on the specified offsets.
  PURE SUBROUTINE SplitMatrixColumns_lsc(this, num_blocks, block_sizes, &
       & split_list)
    !> This matrix to perform this operation on.
    TYPE(Matrix_lsc), INTENT(IN) :: this
    !> Number of blocks to split into.
    INTEGER, INTENT(IN) :: num_blocks
    !> The sizes used for splitting.
    INTEGER, DIMENSION(num_blocks), INTENT(IN) :: block_sizes
    !> 1D array of blocks.
    TYPE(Matrix_lsc), DIMENSION(num_blocks), INTENT(INOUT) :: split_list


    !! Local Data
    INTEGER, DIMENSION(num_blocks+1) :: block_offsets
    !! Temporary variables
    INTEGER :: II, loffset, lcolumns, linner_offset, total_values

    !! Compute Offsets
    block_offsets(1) = 1
    DO II = 2, num_blocks + 1
       block_offsets(II) = block_offsets(II - 1) + block_sizes(II - 1)
    END DO

    !! Split up the columns
    DO II = 1, num_blocks
       !! Temporary variables
       loffset = block_offsets(II)
       lcolumns = block_sizes(II)
       linner_offset = this%outer_index(loffset) + 1

       !! Construct
       CALL ConstructEmptyMatrix(split_list(II), this%rows, lcolumns)

       !! Copy Outer Index
       split_list(II)%outer_index(:) = &
            & this%outer_index(loffset:loffset + lcolumns)
       split_list(II)%outer_index(:) = split_list(II)%outer_index(:) - &
            & split_list(II)%outer_index(1)
       total_values = split_list(II)%outer_index(lcolumns + 1)

       !! Copy Inner Indices and Values
       IF (total_values .GT. 0) THEN
          ALLOCATE(split_list(II)%inner_index(total_values))
          split_list(II)%inner_index(:) = &
               & this%inner_index(linner_offset:linner_offset + total_values - 1)
          ALLOCATE(split_list(II)%values(total_values))
          split_list(II)%values(:) = &
               & this%values(linner_offset:linner_offset + total_values - 1)
       END IF
    END DO
  END SUBROUTINE SplitMatrixColumns_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a triplet list from a matrix.
  PURE SUBROUTINE MatrixToTripletList_lsr(this, triplet_list)
    !> The matrix to construct the triplet list from.
    TYPE(Matrix_lsr), INTENT(IN) :: this
    !> The triplet list we created.
    TYPE(TripletList_r), INTENT(INOUT) :: triplet_list


    !! Helper variables
    INTEGER :: II, JJ, KK
    INTEGER :: elements_per_inner
    INTEGER :: size_of_this

    size_of_this = this%outer_index(this%columns + 1)
    CALL ConstructTripletList(triplet_list, size_of_this)

    KK = 1
    DO II = 1, this%columns
       elements_per_inner = this%outer_index(II + 1) - this%outer_index(II)
       DO JJ = 1, elements_per_inner
          triplet_list%DATA(KK)%index_column = II
          triplet_list%DATA(KK)%index_row = this%inner_index(KK)
          triplet_list%DATA(KK)%point_value = this%values(KK)
          KK = KK + 1
       END DO
    END DO
  END SUBROUTINE MatrixToTripletList_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a triplet list from a matrix.
  PURE SUBROUTINE MatrixToTripletList_lsc(this, triplet_list)
    !> The matrix to construct the triplet list from.
    TYPE(Matrix_lsc), INTENT(IN) :: this
    !> The triplet list we created.
    TYPE(TripletList_c), INTENT(INOUT) :: triplet_list


    !! Helper variables
    INTEGER :: II, JJ, KK
    INTEGER :: elements_per_inner
    INTEGER :: size_of_this

    size_of_this = this%outer_index(this%columns + 1)
    CALL ConstructTripletList(triplet_list, size_of_this)

    KK = 1
    DO II = 1, this%columns
       elements_per_inner = this%outer_index(II + 1) - this%outer_index(II)
       DO JJ = 1, elements_per_inner
          triplet_list%DATA(KK)%index_column = II
          triplet_list%DATA(KK)%index_row = this%inner_index(KK)
          triplet_list%DATA(KK)%point_value = this%values(KK)
          KK = KK + 1
       END DO
    END DO
  END SUBROUTINE MatrixToTripletList_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Print out a sparse matrix to the console.
  SUBROUTINE PrintMatrix_lsr(this, file_name_in)
    !> The matrix to be printed.
    TYPE(Matrix_lsr), INTENT(IN) :: this
    !> Optionally you can pass a file to print to.
    CHARACTER(len=*), OPTIONAL, INTENT(IN) :: file_name_in
    !! Local Data
    TYPE(TripletList_r) :: triplet_list


    !! Local Data
    INTEGER :: file_handler
    INTEGER :: size_of_this
    CHARACTER(LEN = MAX_LINE_LENGTH) :: tempstr
    INTEGER :: II

    !! Process Optional Parameters
    IF (PRESENT(file_name_in)) THEN
       file_handler = 16
       OPEN(unit = file_handler, file = file_name_in)
    ELSE
       file_handler = 6
    END IF

    !! Print
    CALL MatrixToTripletList(this, triplet_list)

    size_of_this = this%outer_index(this%columns + 1)




    WRITE(file_handler,'(A)') "%%MatrixMarket matrix coordinate real general"


    WRITE(file_handler,'(A)') "%"
    CALL WriteMMSize(tempstr, this%rows, this%columns, &
         & INT(size_of_this, KIND = NTLONG))
    WRITE(file_handler,'(A)') ADJUSTL(TRIM(tempstr))
    DO II = 1, size_of_this







       CALL WriteMMLine(tempstr, triplet_list%DATA(II)%index_row, &
            & triplet_list%DATA(II)%index_column, &
            & triplet_list%DATA(II)%point_value)
       WRITE(file_handler,'(A)') ADJUSTL(TRIM(tempstr))

    END DO

    IF (PRESENT(file_name_in)) THEN
       CLOSE(file_handler)
    END IF
    CALL DestructTripletList(triplet_list)

  END SUBROUTINE PrintMatrix_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Print out a sparse matrix to the console.
  SUBROUTINE PrintMatrix_lsc(this, file_name_in)
    !> The matrix to be printed.
    TYPE(Matrix_lsc), INTENT(IN) :: this
    !> Optionally you can pass a file to print to.
    CHARACTER(len=*), OPTIONAL, INTENT(IN) :: file_name_in
    !! Local Data
    TYPE(TripletList_c) :: triplet_list



    !! Local Data
    INTEGER :: file_handler
    INTEGER :: size_of_this
    CHARACTER(LEN = MAX_LINE_LENGTH) :: tempstr
    INTEGER :: II

    !! Process Optional Parameters
    IF (PRESENT(file_name_in)) THEN
       file_handler = 16
       OPEN(unit = file_handler, file = file_name_in)
    ELSE
       file_handler = 6
    END IF

    !! Print
    CALL MatrixToTripletList(this, triplet_list)

    size_of_this = this%outer_index(this%columns + 1)


    WRITE(file_handler,'(A)') "%%MatrixMarket matrix coordinate complex general"




    WRITE(file_handler,'(A)') "%"
    CALL WriteMMSize(tempstr, this%rows, this%columns, &
         & INT(size_of_this, KIND = NTLONG))
    WRITE(file_handler,'(A)') ADJUSTL(TRIM(tempstr))
    DO II = 1, size_of_this

       CALL WriteMMLine(tempstr, triplet_list%DATA(II)%index_row, &
            & triplet_list%DATA(II)%index_column, &
            & REAL(triplet_list%DATA(II)%point_value), &
            & AIMAG(triplet_list%DATA(II)%point_value))
       WRITE(file_handler,'(A)') ADJUSTL(TRIM(tempstr))






    END DO

    IF (PRESENT(file_name_in)) THEN
       CLOSE(file_handler)
    END IF
    CALL DestructTripletList(triplet_list)

  END SUBROUTINE PrintMatrix_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Every value in the matrix is changed into its complex conjugate.
  PURE SUBROUTINE ConjugateMatrix_lsc(this)
    !> The matrix to compute the complex conjugate of.
    TYPE(Matrix_lsc), INTENT(INOUT) :: this

    this%values(:) = CONJG(this%values)
  END SUBROUTINE ConjugateMatrix_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Convert a complex matrix to a real matrix.
  SUBROUTINE ConvertMatrixType_lsrtolsc(cin, rout)
    !> The starting matrix.
    TYPE(Matrix_lsc), INTENT(IN)    :: cin
    !> Real valued matrix.
    TYPE(Matrix_lsr), INTENT(INOUT) :: rout
    !! Local Variables
    TYPE(TripletList_c) :: in_list
    TYPE(TripletList_r) :: out_list

    CALL MatrixToTripletList(cin, in_list)
    CALL ConvertTripletListType(in_list, out_list)
    CALL ConstructMatrixFromTripletList(rout, out_list, cin%rows, cin%columns)

    CALL DestructTripletList(in_list)
    CALL DestructTripletList(out_list)

  END SUBROUTINE ConvertMatrixType_lsrtolsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Convert a real matrix to a complex matrix.
  SUBROUTINE ConvertMatrixType_lsctolsr(rin, cout)
    !> The starting matrix.
    TYPE(Matrix_lsr), INTENT(IN)    :: rin
    !> The complex valued matrix.
    TYPE(Matrix_lsc), INTENT(INOUT) :: cout
    !! Local Variables
    TYPE(TripletList_r) :: in_list
    TYPE(TripletList_c) :: out_list

    CALL MatrixToTripletList(rin, in_list)
    CALL ConvertTripletListType(in_list, out_list)
    CALL ConstructMatrixFromTripletList(cout, out_list, rin%rows, rin%columns)

    CALL DestructTripletList(in_list)
    CALL DestructTripletList(out_list)

  END SUBROUTINE ConvertMatrixType_lsctolsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE SMatrixModule
