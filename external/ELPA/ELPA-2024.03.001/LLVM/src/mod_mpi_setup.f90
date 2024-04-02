












module elpa_mpi_setup
  !use precision
  use iso_c_binding

  type :: elpa_mpi_setup_t
    logical                        :: useMPI

    integer(kind=c_int)            :: mpi_comm_parent, mpi_comm_cols, mpi_comm_rows
    integer(kind=c_int)            :: mpi_comm_parentExternal, mpi_comm_colsExternal, mpi_comm_rowsExternal

    integer(kind=c_int)            :: nRanks_comm_parent
    integer(kind=c_int)            :: nRanks_comm_rows
    integer(kind=c_int)            :: nRanks_comm_cols

    integer(kind=c_int)            :: myRank_comm_parent
    integer(kind=c_int)            :: myRank_comm_rows
    integer(kind=c_int)            :: myRank_comm_cols
    
    integer(kind=c_int)            :: nRanksExternal_comm_parent
    integer(kind=c_int)            :: nRanksExternal_comm_rows
    integer(kind=c_int)            :: nRanksExternal_comm_cols

    integer(kind=c_int)            :: myRankExternal_comm_parent
    integer(kind=c_int)            :: myRankExternal_comm_rows
    integer(kind=c_int)            :: myRankExternal_comm_cols


  end type

end module

