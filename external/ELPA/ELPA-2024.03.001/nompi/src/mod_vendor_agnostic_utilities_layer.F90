












module elpa_gpu_util
  use precision
  use iso_c_binding

!cannot use "../src/GPU/../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)

  interface gpu_memcpy_async_and_stream_synchronize
    module procedure gpu_memcpy_async_and_stream_synchronize_double_scalar
    module procedure gpu_memcpy_async_and_stream_synchronize_single_scalar
    module procedure gpu_memcpy_async_and_stream_synchronize_complex_scalar
    module procedure gpu_memcpy_async_and_stream_synchronize_complex_single_scalar
    module procedure gpu_memcpy_async_and_stream_synchronize_double_1d
    module procedure gpu_memcpy_async_and_stream_synchronize_single_1d
    module procedure gpu_memcpy_async_and_stream_synchronize_complex_1d
    module procedure gpu_memcpy_async_and_stream_synchronize_complex_single_1d
    module procedure gpu_memcpy_async_and_stream_synchronize_double_2d
    module procedure gpu_memcpy_async_and_stream_synchronize_single_2d
    module procedure gpu_memcpy_async_and_stream_synchronize_complex_2d
    module procedure gpu_memcpy_async_and_stream_synchronize_complex_single_2d
    module procedure gpu_memcpy_async_and_stream_synchronize_double_3d
    module procedure gpu_memcpy_async_and_stream_synchronize_single_3d
    module procedure gpu_memcpy_async_and_stream_synchronize_complex_3d
    module procedure gpu_memcpy_async_and_stream_synchronize_complex_single_3d
  end interface

  contains


    subroutine gpu_memcpy_async_and_stream_synchronize_double_scalar &
                (errormessage, devPtr, devPtrOffset, &
                 hostarray, &
                num, direction, my_stream, doSyncBefore, doSyncAfter, doSyncDefault)

      use iso_c_binding
      use elpa_gpu
      use elpa_utilities

      implicit none
      integer(kind=c_intptr_t), intent(in)          :: my_stream
      logical                                       :: successGPU
      character(len=*)                              :: errormessage
      integer(kind=c_intptr_t)                      :: devPtr
      integer(kind=c_intptr_t)                      :: devPtrOffset

      logical                                       :: doSyncBefore
      logical                                       :: doSyncAfter
      logical                                       :: doSyncDefault

      real(kind=c_double)                           :: &
                                             hostArray
      integer(kind=c_intptr_t)                      :: num
      integer(kind=c_int), intent(in)               :: direction

      if (doSyncBefore) then
        successGPU = gpu_stream_synchronize(my_stream)
        call check_memcpy_GPU_f(trim(errormessage), 99,  successGPU)
      endif

      if (direction .eq. gpuMemcpyHostToDevice) then
        successGPU = gpu_memcpy_async(devPtr+devPtrOffset, &
                    int(loc(hostarray),kind=c_intptr_t), &
                                      num, direction, my_stream)
      else if (direction .eq. gpuMemcpyDeviceToHost) then
        successGPU = gpu_memcpy_async(&
                    int(loc(hostarray),kind=c_intptr_t), &
                                      devPtr+devPtrOffset, num, direction, my_stream)
      else
        print *,"gpu_memcpy_async_and_stream_synchronize: unknown error"
        stop 1
      endif

      if (doSyncAfter) then
        successGPU = gpu_stream_synchronize(my_stream)
        call check_memcpy_GPU_f(trim(errormessage), 117,  successGPU)
      endif
      if (doSyncDefault) then
        ! synchronize streamsPerThread; maybe not neccessary
        successGPU = gpu_stream_synchronize()
        call check_memcpy_GPU_f(trim(errormessage), 122,  successGPU)
      endif
    end


    subroutine gpu_memcpy_async_and_stream_synchronize_single_scalar &
                (errormessage, devPtr, devPtrOffset, &
                 hostarray, &
                num, direction, my_stream, doSyncBefore, doSyncAfter, doSyncDefault)

      use iso_c_binding
      use elpa_gpu
      use elpa_utilities

      implicit none
      integer(kind=c_intptr_t), intent(in)          :: my_stream
      logical                                       :: successGPU
      character(len=*)                              :: errormessage
      integer(kind=c_intptr_t)                      :: devPtr
      integer(kind=c_intptr_t)                      :: devPtrOffset

      logical                                       :: doSyncBefore
      logical                                       :: doSyncAfter
      logical                                       :: doSyncDefault

      real(kind=c_float)                            :: &
                                             hostArray
      integer(kind=c_intptr_t)                      :: num
      integer(kind=c_int), intent(in)               :: direction

      if (doSyncBefore) then
        successGPU = gpu_stream_synchronize(my_stream)
        call check_memcpy_GPU_f(trim(errormessage), 154,  successGPU)
      endif

      if (direction .eq. gpuMemcpyHostToDevice) then
        successGPU = gpu_memcpy_async(devPtr+devPtrOffset, &
                    int(loc(hostarray),kind=c_intptr_t), &
                                      num, direction, my_stream)
      else if (direction .eq. gpuMemcpyDeviceToHost) then
        successGPU = gpu_memcpy_async(&
                    int(loc(hostarray),kind=c_intptr_t), &
                                      devPtr+devPtrOffset, num, direction, my_stream)
      else
        print *,"gpu_memcpy_async_and_stream_synchronize: unknown error"
        stop 1
      endif

      if (doSyncAfter) then
        successGPU = gpu_stream_synchronize(my_stream)
        call check_memcpy_GPU_f(trim(errormessage), 172,  successGPU)
      endif
      if (doSyncDefault) then
        ! synchronize streamsPerThread; maybe not neccessary
        successGPU = gpu_stream_synchronize()
        call check_memcpy_GPU_f(trim(errormessage), 177,  successGPU)
      endif
    end


    subroutine gpu_memcpy_async_and_stream_synchronize_complex_scalar &
                (errormessage, devPtr, devPtrOffset, &
                 hostarray, &
                num, direction, my_stream, doSyncBefore, doSyncAfter, doSyncDefault)

      use iso_c_binding
      use elpa_gpu
      use elpa_utilities

      implicit none
      integer(kind=c_intptr_t), intent(in)          :: my_stream
      logical                                       :: successGPU
      character(len=*)                              :: errormessage
      integer(kind=c_intptr_t)                      :: devPtr
      integer(kind=c_intptr_t)                      :: devPtrOffset

      logical                                       :: doSyncBefore
      logical                                       :: doSyncAfter
      logical                                       :: doSyncDefault

      complex(kind=c_double_complex)                :: &
                                             hostArray
      integer(kind=c_intptr_t)                      :: num
      integer(kind=c_int), intent(in)               :: direction

      if (doSyncBefore) then
        successGPU = gpu_stream_synchronize(my_stream)
        call check_memcpy_GPU_f(trim(errormessage), 209,  successGPU)
      endif

      if (direction .eq. gpuMemcpyHostToDevice) then
        successGPU = gpu_memcpy_async(devPtr+devPtrOffset, &
                    int(loc(hostarray),kind=c_intptr_t), &
                                      num, direction, my_stream)
      else if (direction .eq. gpuMemcpyDeviceToHost) then
        successGPU = gpu_memcpy_async(&
                    int(loc(hostarray),kind=c_intptr_t), &
                                      devPtr+devPtrOffset, num, direction, my_stream)
      else
        print *,"gpu_memcpy_async_and_stream_synchronize: unknown error"
        stop 1
      endif

      if (doSyncAfter) then
        successGPU = gpu_stream_synchronize(my_stream)
        call check_memcpy_GPU_f(trim(errormessage), 227,  successGPU)
      endif
      if (doSyncDefault) then
        ! synchronize streamsPerThread; maybe not neccessary
        successGPU = gpu_stream_synchronize()
        call check_memcpy_GPU_f(trim(errormessage), 232,  successGPU)
      endif
    end


    subroutine gpu_memcpy_async_and_stream_synchronize_complex_single_scalar &
                (errormessage, devPtr, devPtrOffset, &
                 hostarray, &
                num, direction, my_stream, doSyncBefore, doSyncAfter, doSyncDefault)

      use iso_c_binding
      use elpa_gpu
      use elpa_utilities

      implicit none
      integer(kind=c_intptr_t), intent(in)          :: my_stream
      logical                                       :: successGPU
      character(len=*)                              :: errormessage
      integer(kind=c_intptr_t)                      :: devPtr
      integer(kind=c_intptr_t)                      :: devPtrOffset

      logical                                       :: doSyncBefore
      logical                                       :: doSyncAfter
      logical                                       :: doSyncDefault

      complex(kind=c_float_complex)                 :: &
                                             hostArray
      integer(kind=c_intptr_t)                      :: num
      integer(kind=c_int), intent(in)               :: direction

      if (doSyncBefore) then
        successGPU = gpu_stream_synchronize(my_stream)
        call check_memcpy_GPU_f(trim(errormessage), 264,  successGPU)
      endif

      if (direction .eq. gpuMemcpyHostToDevice) then
        successGPU = gpu_memcpy_async(devPtr+devPtrOffset, &
                    int(loc(hostarray),kind=c_intptr_t), &
                                      num, direction, my_stream)
      else if (direction .eq. gpuMemcpyDeviceToHost) then
        successGPU = gpu_memcpy_async(&
                    int(loc(hostarray),kind=c_intptr_t), &
                                      devPtr+devPtrOffset, num, direction, my_stream)
      else
        print *,"gpu_memcpy_async_and_stream_synchronize: unknown error"
        stop 1
      endif

      if (doSyncAfter) then
        successGPU = gpu_stream_synchronize(my_stream)
        call check_memcpy_GPU_f(trim(errormessage), 282,  successGPU)
      endif
      if (doSyncDefault) then
        ! synchronize streamsPerThread; maybe not neccessary
        successGPU = gpu_stream_synchronize()
        call check_memcpy_GPU_f(trim(errormessage), 287,  successGPU)
      endif
    end


    subroutine gpu_memcpy_async_and_stream_synchronize_double_1d &
                (errormessage, devPtr, devPtrOffset, &
                 hostarray, &
                 off1, &
                num, direction, my_stream, doSyncBefore, doSyncAfter, doSyncDefault)

      use iso_c_binding
      use elpa_gpu
      use elpa_utilities

      implicit none
      integer(kind=c_intptr_t), intent(in)          :: my_stream
      logical                                       :: successGPU
      character(len=*)                              :: errormessage
      integer(kind=c_intptr_t)                      :: devPtr
      integer(kind=c_intptr_t)                      :: devPtrOffset

      integer(kind=ik), intent(in)                  :: off1
      logical                                       :: doSyncBefore
      logical                                       :: doSyncAfter
      logical                                       :: doSyncDefault

      real(kind=c_double)                           :: &
                                         hostArray(:)
      integer(kind=c_intptr_t)                      :: num
      integer(kind=c_int), intent(in)               :: direction

      if (doSyncBefore) then
        successGPU = gpu_stream_synchronize(my_stream)
        call check_memcpy_GPU_f(trim(errormessage), 321,  successGPU)
      endif

      if (direction .eq. gpuMemcpyHostToDevice) then
        successGPU = gpu_memcpy_async(devPtr+devPtrOffset, &
                    int(loc(hostarray(off1)),kind=c_intptr_t), &
                                      num, direction, my_stream)
      else if (direction .eq. gpuMemcpyDeviceToHost) then
        successGPU = gpu_memcpy_async(&
                    int(loc(hostarray(off1)),kind=c_intptr_t), &
                                      devPtr+devPtrOffset, num, direction, my_stream)
      else
        print *,"gpu_memcpy_async_and_stream_synchronize: unknown error"
        stop 1
      endif

      if (doSyncAfter) then
        successGPU = gpu_stream_synchronize(my_stream)
        call check_memcpy_GPU_f(trim(errormessage), 339,  successGPU)
      endif
      if (doSyncDefault) then
        ! synchronize streamsPerThread; maybe not neccessary
        successGPU = gpu_stream_synchronize()
        call check_memcpy_GPU_f(trim(errormessage), 344,  successGPU)
      endif
    end


    subroutine gpu_memcpy_async_and_stream_synchronize_single_1d &
                (errormessage, devPtr, devPtrOffset, &
                 hostarray, &
                 off1, &
                num, direction, my_stream, doSyncBefore, doSyncAfter, doSyncDefault)

      use iso_c_binding
      use elpa_gpu
      use elpa_utilities

      implicit none
      integer(kind=c_intptr_t), intent(in)          :: my_stream
      logical                                       :: successGPU
      character(len=*)                              :: errormessage
      integer(kind=c_intptr_t)                      :: devPtr
      integer(kind=c_intptr_t)                      :: devPtrOffset

      integer(kind=ik), intent(in)                  :: off1
      logical                                       :: doSyncBefore
      logical                                       :: doSyncAfter
      logical                                       :: doSyncDefault

      real(kind=c_float)                            :: &
                                         hostArray(:)
      integer(kind=c_intptr_t)                      :: num
      integer(kind=c_int), intent(in)               :: direction

      if (doSyncBefore) then
        successGPU = gpu_stream_synchronize(my_stream)
        call check_memcpy_GPU_f(trim(errormessage), 378,  successGPU)
      endif

      if (direction .eq. gpuMemcpyHostToDevice) then
        successGPU = gpu_memcpy_async(devPtr+devPtrOffset, &
                    int(loc(hostarray(off1)),kind=c_intptr_t), &
                                      num, direction, my_stream)
      else if (direction .eq. gpuMemcpyDeviceToHost) then
        successGPU = gpu_memcpy_async(&
                    int(loc(hostarray(off1)),kind=c_intptr_t), &
                                      devPtr+devPtrOffset, num, direction, my_stream)
      else
        print *,"gpu_memcpy_async_and_stream_synchronize: unknown error"
        stop 1
      endif

      if (doSyncAfter) then
        successGPU = gpu_stream_synchronize(my_stream)
        call check_memcpy_GPU_f(trim(errormessage), 396,  successGPU)
      endif
      if (doSyncDefault) then
        ! synchronize streamsPerThread; maybe not neccessary
        successGPU = gpu_stream_synchronize()
        call check_memcpy_GPU_f(trim(errormessage), 401,  successGPU)
      endif
    end


    subroutine gpu_memcpy_async_and_stream_synchronize_complex_1d &
                (errormessage, devPtr, devPtrOffset, &
                 hostarray, &
                 off1, &
                num, direction, my_stream, doSyncBefore, doSyncAfter, doSyncDefault)

      use iso_c_binding
      use elpa_gpu
      use elpa_utilities

      implicit none
      integer(kind=c_intptr_t), intent(in)          :: my_stream
      logical                                       :: successGPU
      character(len=*)                              :: errormessage
      integer(kind=c_intptr_t)                      :: devPtr
      integer(kind=c_intptr_t)                      :: devPtrOffset

      integer(kind=ik), intent(in)                  :: off1
      logical                                       :: doSyncBefore
      logical                                       :: doSyncAfter
      logical                                       :: doSyncDefault

      complex(kind=c_double_complex)                :: &
                                         hostArray(:)
      integer(kind=c_intptr_t)                      :: num
      integer(kind=c_int), intent(in)               :: direction

      if (doSyncBefore) then
        successGPU = gpu_stream_synchronize(my_stream)
        call check_memcpy_GPU_f(trim(errormessage), 435,  successGPU)
      endif

      if (direction .eq. gpuMemcpyHostToDevice) then
        successGPU = gpu_memcpy_async(devPtr+devPtrOffset, &
                    int(loc(hostarray(off1)),kind=c_intptr_t), &
                                      num, direction, my_stream)
      else if (direction .eq. gpuMemcpyDeviceToHost) then
        successGPU = gpu_memcpy_async(&
                    int(loc(hostarray(off1)),kind=c_intptr_t), &
                                      devPtr+devPtrOffset, num, direction, my_stream)
      else
        print *,"gpu_memcpy_async_and_stream_synchronize: unknown error"
        stop 1
      endif

      if (doSyncAfter) then
        successGPU = gpu_stream_synchronize(my_stream)
        call check_memcpy_GPU_f(trim(errormessage), 453,  successGPU)
      endif
      if (doSyncDefault) then
        ! synchronize streamsPerThread; maybe not neccessary
        successGPU = gpu_stream_synchronize()
        call check_memcpy_GPU_f(trim(errormessage), 458,  successGPU)
      endif
    end


    subroutine gpu_memcpy_async_and_stream_synchronize_complex_single_1d &
                (errormessage, devPtr, devPtrOffset, &
                 hostarray, &
                 off1, &
                num, direction, my_stream, doSyncBefore, doSyncAfter, doSyncDefault)

      use iso_c_binding
      use elpa_gpu
      use elpa_utilities

      implicit none
      integer(kind=c_intptr_t), intent(in)          :: my_stream
      logical                                       :: successGPU
      character(len=*)                              :: errormessage
      integer(kind=c_intptr_t)                      :: devPtr
      integer(kind=c_intptr_t)                      :: devPtrOffset

      integer(kind=ik), intent(in)                  :: off1
      logical                                       :: doSyncBefore
      logical                                       :: doSyncAfter
      logical                                       :: doSyncDefault

      complex(kind=c_float_complex)                 :: &
                                         hostArray(:)
      integer(kind=c_intptr_t)                      :: num
      integer(kind=c_int), intent(in)               :: direction

      if (doSyncBefore) then
        successGPU = gpu_stream_synchronize(my_stream)
        call check_memcpy_GPU_f(trim(errormessage), 492,  successGPU)
      endif

      if (direction .eq. gpuMemcpyHostToDevice) then
        successGPU = gpu_memcpy_async(devPtr+devPtrOffset, &
                    int(loc(hostarray(off1)),kind=c_intptr_t), &
                                      num, direction, my_stream)
      else if (direction .eq. gpuMemcpyDeviceToHost) then
        successGPU = gpu_memcpy_async(&
                    int(loc(hostarray(off1)),kind=c_intptr_t), &
                                      devPtr+devPtrOffset, num, direction, my_stream)
      else
        print *,"gpu_memcpy_async_and_stream_synchronize: unknown error"
        stop 1
      endif

      if (doSyncAfter) then
        successGPU = gpu_stream_synchronize(my_stream)
        call check_memcpy_GPU_f(trim(errormessage), 510,  successGPU)
      endif
      if (doSyncDefault) then
        ! synchronize streamsPerThread; maybe not neccessary
        successGPU = gpu_stream_synchronize()
        call check_memcpy_GPU_f(trim(errormessage), 515,  successGPU)
      endif
    end


    subroutine gpu_memcpy_async_and_stream_synchronize_double_2d &
                (errormessage, devPtr, devPtrOffset, &
                 hostarray, &
                 off1, off2, &
                num, direction, my_stream, doSyncBefore, doSyncAfter, doSyncDefault)

      use iso_c_binding
      use elpa_gpu
      use elpa_utilities

      implicit none
      integer(kind=c_intptr_t), intent(in)          :: my_stream
      logical                                       :: successGPU
      character(len=*)                              :: errormessage
      integer(kind=c_intptr_t)                      :: devPtr
      integer(kind=c_intptr_t)                      :: devPtrOffset

      integer(kind=ik), intent(in)                  :: off1
      integer(kind=ik), intent(in)                  :: off2
      logical                                       :: doSyncBefore
      logical                                       :: doSyncAfter
      logical                                       :: doSyncDefault

      real(kind=c_double)                           :: &
                                         hostArray(:,:)
      integer(kind=c_intptr_t)                      :: num
      integer(kind=c_int), intent(in)               :: direction

      if (doSyncBefore) then
        successGPU = gpu_stream_synchronize(my_stream)
        call check_memcpy_GPU_f(trim(errormessage), 550,  successGPU)
      endif

      if (direction .eq. gpuMemcpyHostToDevice) then
        successGPU = gpu_memcpy_async(devPtr+devPtrOffset, &
                    int(loc(hostarray(off1,off2)),kind=c_intptr_t), &
                                      num, direction, my_stream)
      else if (direction .eq. gpuMemcpyDeviceToHost) then
        successGPU = gpu_memcpy_async(&
                    int(loc(hostarray(off1,off2)),kind=c_intptr_t), &
                                      devPtr+devPtrOffset, num, direction, my_stream)
      else
        print *,"gpu_memcpy_async_and_stream_synchronize: unknown error"
        stop 1
      endif

      if (doSyncAfter) then
        successGPU = gpu_stream_synchronize(my_stream)
        call check_memcpy_GPU_f(trim(errormessage), 568,  successGPU)
      endif
      if (doSyncDefault) then
        ! synchronize streamsPerThread; maybe not neccessary
        successGPU = gpu_stream_synchronize()
        call check_memcpy_GPU_f(trim(errormessage), 573,  successGPU)
      endif
    end


    subroutine gpu_memcpy_async_and_stream_synchronize_single_2d &
                (errormessage, devPtr, devPtrOffset, &
                 hostarray, &
                 off1, off2, &
                num, direction, my_stream, doSyncBefore, doSyncAfter, doSyncDefault)

      use iso_c_binding
      use elpa_gpu
      use elpa_utilities

      implicit none
      integer(kind=c_intptr_t), intent(in)          :: my_stream
      logical                                       :: successGPU
      character(len=*)                              :: errormessage
      integer(kind=c_intptr_t)                      :: devPtr
      integer(kind=c_intptr_t)                      :: devPtrOffset

      integer(kind=ik), intent(in)                  :: off1
      integer(kind=ik), intent(in)                  :: off2
      logical                                       :: doSyncBefore
      logical                                       :: doSyncAfter
      logical                                       :: doSyncDefault

      real(kind=c_float)                            :: &
                                         hostArray(:,:)
      integer(kind=c_intptr_t)                      :: num
      integer(kind=c_int), intent(in)               :: direction

      if (doSyncBefore) then
        successGPU = gpu_stream_synchronize(my_stream)
        call check_memcpy_GPU_f(trim(errormessage), 608,  successGPU)
      endif

      if (direction .eq. gpuMemcpyHostToDevice) then
        successGPU = gpu_memcpy_async(devPtr+devPtrOffset, &
                    int(loc(hostarray(off1,off2)),kind=c_intptr_t), &
                                      num, direction, my_stream)
      else if (direction .eq. gpuMemcpyDeviceToHost) then
        successGPU = gpu_memcpy_async(&
                    int(loc(hostarray(off1,off2)),kind=c_intptr_t), &
                                      devPtr+devPtrOffset, num, direction, my_stream)
      else
        print *,"gpu_memcpy_async_and_stream_synchronize: unknown error"
        stop 1
      endif

      if (doSyncAfter) then
        successGPU = gpu_stream_synchronize(my_stream)
        call check_memcpy_GPU_f(trim(errormessage), 626,  successGPU)
      endif
      if (doSyncDefault) then
        ! synchronize streamsPerThread; maybe not neccessary
        successGPU = gpu_stream_synchronize()
        call check_memcpy_GPU_f(trim(errormessage), 631,  successGPU)
      endif
    end


    subroutine gpu_memcpy_async_and_stream_synchronize_complex_2d &
                (errormessage, devPtr, devPtrOffset, &
                 hostarray, &
                 off1, off2, &
                num, direction, my_stream, doSyncBefore, doSyncAfter, doSyncDefault)

      use iso_c_binding
      use elpa_gpu
      use elpa_utilities

      implicit none
      integer(kind=c_intptr_t), intent(in)          :: my_stream
      logical                                       :: successGPU
      character(len=*)                              :: errormessage
      integer(kind=c_intptr_t)                      :: devPtr
      integer(kind=c_intptr_t)                      :: devPtrOffset

      integer(kind=ik), intent(in)                  :: off1
      integer(kind=ik), intent(in)                  :: off2
      logical                                       :: doSyncBefore
      logical                                       :: doSyncAfter
      logical                                       :: doSyncDefault

      complex(kind=c_double_complex)                :: &
                                         hostArray(:,:)
      integer(kind=c_intptr_t)                      :: num
      integer(kind=c_int), intent(in)               :: direction

      if (doSyncBefore) then
        successGPU = gpu_stream_synchronize(my_stream)
        call check_memcpy_GPU_f(trim(errormessage), 666,  successGPU)
      endif

      if (direction .eq. gpuMemcpyHostToDevice) then
        successGPU = gpu_memcpy_async(devPtr+devPtrOffset, &
                    int(loc(hostarray(off1,off2)),kind=c_intptr_t), &
                                      num, direction, my_stream)
      else if (direction .eq. gpuMemcpyDeviceToHost) then
        successGPU = gpu_memcpy_async(&
                    int(loc(hostarray(off1,off2)),kind=c_intptr_t), &
                                      devPtr+devPtrOffset, num, direction, my_stream)
      else
        print *,"gpu_memcpy_async_and_stream_synchronize: unknown error"
        stop 1
      endif

      if (doSyncAfter) then
        successGPU = gpu_stream_synchronize(my_stream)
        call check_memcpy_GPU_f(trim(errormessage), 684,  successGPU)
      endif
      if (doSyncDefault) then
        ! synchronize streamsPerThread; maybe not neccessary
        successGPU = gpu_stream_synchronize()
        call check_memcpy_GPU_f(trim(errormessage), 689,  successGPU)
      endif
    end


    subroutine gpu_memcpy_async_and_stream_synchronize_complex_single_2d &
                (errormessage, devPtr, devPtrOffset, &
                 hostarray, &
                 off1, off2, &
                num, direction, my_stream, doSyncBefore, doSyncAfter, doSyncDefault)

      use iso_c_binding
      use elpa_gpu
      use elpa_utilities

      implicit none
      integer(kind=c_intptr_t), intent(in)          :: my_stream
      logical                                       :: successGPU
      character(len=*)                              :: errormessage
      integer(kind=c_intptr_t)                      :: devPtr
      integer(kind=c_intptr_t)                      :: devPtrOffset

      integer(kind=ik), intent(in)                  :: off1
      integer(kind=ik), intent(in)                  :: off2
      logical                                       :: doSyncBefore
      logical                                       :: doSyncAfter
      logical                                       :: doSyncDefault

      complex(kind=c_float_complex)                 :: &
                                         hostArray(:,:)
      integer(kind=c_intptr_t)                      :: num
      integer(kind=c_int), intent(in)               :: direction

      if (doSyncBefore) then
        successGPU = gpu_stream_synchronize(my_stream)
        call check_memcpy_GPU_f(trim(errormessage), 724,  successGPU)
      endif

      if (direction .eq. gpuMemcpyHostToDevice) then
        successGPU = gpu_memcpy_async(devPtr+devPtrOffset, &
                    int(loc(hostarray(off1,off2)),kind=c_intptr_t), &
                                      num, direction, my_stream)
      else if (direction .eq. gpuMemcpyDeviceToHost) then
        successGPU = gpu_memcpy_async(&
                    int(loc(hostarray(off1,off2)),kind=c_intptr_t), &
                                      devPtr+devPtrOffset, num, direction, my_stream)
      else
        print *,"gpu_memcpy_async_and_stream_synchronize: unknown error"
        stop 1
      endif

      if (doSyncAfter) then
        successGPU = gpu_stream_synchronize(my_stream)
        call check_memcpy_GPU_f(trim(errormessage), 742,  successGPU)
      endif
      if (doSyncDefault) then
        ! synchronize streamsPerThread; maybe not neccessary
        successGPU = gpu_stream_synchronize()
        call check_memcpy_GPU_f(trim(errormessage), 747,  successGPU)
      endif
    end


    subroutine gpu_memcpy_async_and_stream_synchronize_double_3d &
                (errormessage, devPtr, devPtrOffset, &
                 hostarray, &
                 off1, off2, off3, &
                num, direction, my_stream, doSyncBefore, doSyncAfter, doSyncDefault)

      use iso_c_binding
      use elpa_gpu
      use elpa_utilities

      implicit none
      integer(kind=c_intptr_t), intent(in)          :: my_stream
      logical                                       :: successGPU
      character(len=*)                              :: errormessage
      integer(kind=c_intptr_t)                      :: devPtr
      integer(kind=c_intptr_t)                      :: devPtrOffset

      integer(kind=ik), intent(in)                  :: off1
      integer(kind=ik), intent(in)                  :: off2
      integer(kind=ik), intent(in)                  :: off3
      logical                                       :: doSyncBefore
      logical                                       :: doSyncAfter
      logical                                       :: doSyncDefault

      real(kind=c_double)                           :: &
                                         hostArray(:,:,:)
      integer(kind=c_intptr_t)                      :: num
      integer(kind=c_int), intent(in)               :: direction

      if (doSyncBefore) then
        successGPU = gpu_stream_synchronize(my_stream)
        call check_memcpy_GPU_f(trim(errormessage), 783,  successGPU)
      endif

      if (direction .eq. gpuMemcpyHostToDevice) then
        successGPU = gpu_memcpy_async(devPtr+devPtrOffset, &
                    int(loc(hostarray(off1,off2,off3)),kind=c_intptr_t), &
                                      num, direction, my_stream)
      else if (direction .eq. gpuMemcpyDeviceToHost) then
        successGPU = gpu_memcpy_async(&
                    int(loc(hostarray(off1,off2,off3)),kind=c_intptr_t), &
                                      devPtr+devPtrOffset, num, direction, my_stream)
      else
        print *,"gpu_memcpy_async_and_stream_synchronize: unknown error"
        stop 1
      endif

      if (doSyncAfter) then
        successGPU = gpu_stream_synchronize(my_stream)
        call check_memcpy_GPU_f(trim(errormessage), 801,  successGPU)
      endif
      if (doSyncDefault) then
        ! synchronize streamsPerThread; maybe not neccessary
        successGPU = gpu_stream_synchronize()
        call check_memcpy_GPU_f(trim(errormessage), 806,  successGPU)
      endif
    end


    subroutine gpu_memcpy_async_and_stream_synchronize_single_3d &
                (errormessage, devPtr, devPtrOffset, &
                 hostarray, &
                 off1, off2, off3, &
                num, direction, my_stream, doSyncBefore, doSyncAfter, doSyncDefault)

      use iso_c_binding
      use elpa_gpu
      use elpa_utilities

      implicit none
      integer(kind=c_intptr_t), intent(in)          :: my_stream
      logical                                       :: successGPU
      character(len=*)                              :: errormessage
      integer(kind=c_intptr_t)                      :: devPtr
      integer(kind=c_intptr_t)                      :: devPtrOffset

      integer(kind=ik), intent(in)                  :: off1
      integer(kind=ik), intent(in)                  :: off2
      integer(kind=ik), intent(in)                  :: off3
      logical                                       :: doSyncBefore
      logical                                       :: doSyncAfter
      logical                                       :: doSyncDefault

      real(kind=c_float)                            :: &
                                         hostArray(:,:,:)
      integer(kind=c_intptr_t)                      :: num
      integer(kind=c_int), intent(in)               :: direction

      if (doSyncBefore) then
        successGPU = gpu_stream_synchronize(my_stream)
        call check_memcpy_GPU_f(trim(errormessage), 842,  successGPU)
      endif

      if (direction .eq. gpuMemcpyHostToDevice) then
        successGPU = gpu_memcpy_async(devPtr+devPtrOffset, &
                    int(loc(hostarray(off1,off2,off3)),kind=c_intptr_t), &
                                      num, direction, my_stream)
      else if (direction .eq. gpuMemcpyDeviceToHost) then
        successGPU = gpu_memcpy_async(&
                    int(loc(hostarray(off1,off2,off3)),kind=c_intptr_t), &
                                      devPtr+devPtrOffset, num, direction, my_stream)
      else
        print *,"gpu_memcpy_async_and_stream_synchronize: unknown error"
        stop 1
      endif

      if (doSyncAfter) then
        successGPU = gpu_stream_synchronize(my_stream)
        call check_memcpy_GPU_f(trim(errormessage), 860,  successGPU)
      endif
      if (doSyncDefault) then
        ! synchronize streamsPerThread; maybe not neccessary
        successGPU = gpu_stream_synchronize()
        call check_memcpy_GPU_f(trim(errormessage), 865,  successGPU)
      endif
    end


    subroutine gpu_memcpy_async_and_stream_synchronize_complex_3d &
                (errormessage, devPtr, devPtrOffset, &
                 hostarray, &
                 off1, off2, off3, &
                num, direction, my_stream, doSyncBefore, doSyncAfter, doSyncDefault)

      use iso_c_binding
      use elpa_gpu
      use elpa_utilities

      implicit none
      integer(kind=c_intptr_t), intent(in)          :: my_stream
      logical                                       :: successGPU
      character(len=*)                              :: errormessage
      integer(kind=c_intptr_t)                      :: devPtr
      integer(kind=c_intptr_t)                      :: devPtrOffset

      integer(kind=ik), intent(in)                  :: off1
      integer(kind=ik), intent(in)                  :: off2
      integer(kind=ik), intent(in)                  :: off3
      logical                                       :: doSyncBefore
      logical                                       :: doSyncAfter
      logical                                       :: doSyncDefault

      complex(kind=c_double_complex)                :: &
                                         hostArray(:,:,:)
      integer(kind=c_intptr_t)                      :: num
      integer(kind=c_int), intent(in)               :: direction

      if (doSyncBefore) then
        successGPU = gpu_stream_synchronize(my_stream)
        call check_memcpy_GPU_f(trim(errormessage), 901,  successGPU)
      endif

      if (direction .eq. gpuMemcpyHostToDevice) then
        successGPU = gpu_memcpy_async(devPtr+devPtrOffset, &
                    int(loc(hostarray(off1,off2,off3)),kind=c_intptr_t), &
                                      num, direction, my_stream)
      else if (direction .eq. gpuMemcpyDeviceToHost) then
        successGPU = gpu_memcpy_async(&
                    int(loc(hostarray(off1,off2,off3)),kind=c_intptr_t), &
                                      devPtr+devPtrOffset, num, direction, my_stream)
      else
        print *,"gpu_memcpy_async_and_stream_synchronize: unknown error"
        stop 1
      endif

      if (doSyncAfter) then
        successGPU = gpu_stream_synchronize(my_stream)
        call check_memcpy_GPU_f(trim(errormessage), 919,  successGPU)
      endif
      if (doSyncDefault) then
        ! synchronize streamsPerThread; maybe not neccessary
        successGPU = gpu_stream_synchronize()
        call check_memcpy_GPU_f(trim(errormessage), 924,  successGPU)
      endif
    end


    subroutine gpu_memcpy_async_and_stream_synchronize_complex_single_3d &
                (errormessage, devPtr, devPtrOffset, &
                 hostarray, &
                 off1, off2, off3, &
                num, direction, my_stream, doSyncBefore, doSyncAfter, doSyncDefault)

      use iso_c_binding
      use elpa_gpu
      use elpa_utilities

      implicit none
      integer(kind=c_intptr_t), intent(in)          :: my_stream
      logical                                       :: successGPU
      character(len=*)                              :: errormessage
      integer(kind=c_intptr_t)                      :: devPtr
      integer(kind=c_intptr_t)                      :: devPtrOffset

      integer(kind=ik), intent(in)                  :: off1
      integer(kind=ik), intent(in)                  :: off2
      integer(kind=ik), intent(in)                  :: off3
      logical                                       :: doSyncBefore
      logical                                       :: doSyncAfter
      logical                                       :: doSyncDefault

      complex(kind=c_float_complex)                 :: &
                                         hostArray(:,:,:)
      integer(kind=c_intptr_t)                      :: num
      integer(kind=c_int), intent(in)               :: direction

      if (doSyncBefore) then
        successGPU = gpu_stream_synchronize(my_stream)
        call check_memcpy_GPU_f(trim(errormessage), 960,  successGPU)
      endif

      if (direction .eq. gpuMemcpyHostToDevice) then
        successGPU = gpu_memcpy_async(devPtr+devPtrOffset, &
                    int(loc(hostarray(off1,off2,off3)),kind=c_intptr_t), &
                                      num, direction, my_stream)
      else if (direction .eq. gpuMemcpyDeviceToHost) then
        successGPU = gpu_memcpy_async(&
                    int(loc(hostarray(off1,off2,off3)),kind=c_intptr_t), &
                                      devPtr+devPtrOffset, num, direction, my_stream)
      else
        print *,"gpu_memcpy_async_and_stream_synchronize: unknown error"
        stop 1
      endif

      if (doSyncAfter) then
        successGPU = gpu_stream_synchronize(my_stream)
        call check_memcpy_GPU_f(trim(errormessage), 978,  successGPU)
      endif
      if (doSyncDefault) then
        ! synchronize streamsPerThread; maybe not neccessary
        successGPU = gpu_stream_synchronize()
        call check_memcpy_GPU_f(trim(errormessage), 983,  successGPU)
      endif
    end


end module

