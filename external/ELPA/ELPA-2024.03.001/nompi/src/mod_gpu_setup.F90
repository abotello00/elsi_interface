












module elpa_gpu_setup
  !use precision
  use iso_c_binding

  type :: elpa_gpu_setup_t
    integer(kind=c_int)            :: use_gpu_vendor
    logical                        :: gpuIsAssigned

    ! per task information should be stored elsewhere
    integer(kind=C_intptr_T), allocatable :: gpublasHandleArray(:)
    integer(kind=C_intptr_T), allocatable :: gpusolverHandleArray(:)
    integer(kind=c_int), allocatable      :: gpuDeviceArray(:)
    integer(kind=c_intptr_t)              :: my_stream

    integer(kind=C_intptr_T), allocatable :: cublasHandleArray(:)
    integer(kind=C_intptr_T), allocatable :: cusolverHandleArray(:)
    integer(kind=c_int), allocatable      :: cudaDeviceArray(:)

    integer(kind=C_intptr_T), allocatable :: rocblasHandleArray(:)
    integer(kind=C_intptr_T), allocatable :: rocsolverHandleArray(:)
    integer(kind=c_int), allocatable      :: hipDeviceArray(:)

    integer(kind=C_intptr_T), allocatable :: syclHandleArray(:)
    integer(kind=C_intptr_T), allocatable :: syclsolverHandleArray(:)
    integer(kind=c_int), allocatable      :: syclDeviceArray(:)

    integer(kind=C_intptr_T), allocatable :: openmpOffloadHandleArray(:)
    integer(kind=C_intptr_T), allocatable :: openmpOffloadsolverHandleArray(:)
    integer(kind=c_int), allocatable      :: openmpOffloadDeviceArray(:)

    logical                               :: gpuAlreadySet
    integer(kind=c_int)        :: gpuDevAttrMaxThreadsPerBlock  = 0
    integer(kind=c_int)        :: gpuDevAttrMaxBlockDimX        = 1
    integer(kind=c_int)        :: gpuDevAttrMaxBlockDimY        = 2
    integer(kind=c_int)        :: gpuDevAttrMaxBlockDimZ        = 3
    integer(kind=c_int)        :: gpuDevAttrMaxGridDimX         = 4
    integer(kind=c_int)        :: gpuDevAttrMaxGridDimY         = 5
    integer(kind=c_int)        :: gpuDevAttrMaxGridDimZ         = 6
    integer(kind=c_int)        :: gpuDevAttrWarpSize            = 7
    integer(kind=c_int)        :: gpuDevAttrMultiProcessorCount = 8


    integer(kind=c_int)                   :: gpublasVersion
    integer(kind=c_int)                   :: rocblasVersion
    integer(kind=c_int)                   :: cublasVersion

    integer(kind=c_int)                   :: gpusPerNode

    integer(kind=c_int)                   :: nvidiaSMcount
    integer(kind=c_int)                   :: amdSMcount
    integer(kind=c_int)                   :: gpuSMcount

    integer(kind=c_int)                   :: nvidiaMaxThreadsPerBlock
    integer(kind=c_int)                   :: amdMaxThreadsPerBlock
    integer(kind=c_int)                   :: gpuMaxThreadsPerBlock

    integer(kind=c_int)                   :: nvidiaDevMaxBlockDimX
    integer(kind=c_int)                   :: amdDevMaxBlockDimX
    integer(kind=c_int)                   :: gpuDevMaxBlockDimX

    integer(kind=c_int)                   :: nvidiaDevMaxBlockDimY
    integer(kind=c_int)                   :: amdDevMaxBlockDimY
    integer(kind=c_int)                   :: gpuDevMaxBlockDimY

    integer(kind=c_int)                   :: nvidiaDevMaxBlockDimZ
    integer(kind=c_int)                   :: amdDevMaxBlockDimZ
    integer(kind=c_int)                   :: gpuDevMaxBlockDimZ

    integer(kind=c_int)                   :: nvidiaDevMaxGridDimX
    integer(kind=c_int)                   :: amdDevMaxGridDimX
    integer(kind=c_int)                   :: gpuDevMaxGridDimX

    integer(kind=c_int)                   :: nvidiaDevMaxGridDimY
    integer(kind=c_int)                   :: amdDevMaxGridDimY
    integer(kind=c_int)                   :: gpuDevMaxGridDimY

    integer(kind=c_int)                   :: nvidiaDevMaxGridDimZ
    integer(kind=c_int)                   :: amdDevMaxGridDimZ
    integer(kind=c_int)                   :: gpuDevMaxGridDimZ

    integer(kind=c_int)                   :: nvidiaDevWarpSize
    integer(kind=c_int)                   :: amdDevWarpSize
    integer(kind=c_int)                   :: gpuDevWarpSize

    integer(kind=c_intptr_t)              :: ccl_comm_rows, ccl_comm_cols, ccl_comm_all
  end type

end module

