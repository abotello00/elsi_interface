










!    Copyright 2021, A. Marek
!
!    This file is part of ELPA.
!
!    The ELPA library was originally created by the ELPA consortium,
!    consisting of the following organizations:
!
!    - Rechenzentrum Garching der Max-Planck-Gesellschaft (RZG),
!    - Bergische Universität Wuppertal, Lehrstuhl für angewandte
!      Informatik,
!    - Technische Universität München, Lehrstuhl für Informatik mit
!      Schwerpunkt Wissenschaftliches Rechnen ,
!    - Fritz-Haber-Institut, Berlin, Abt. Theorie,
!    - Max-Plack-Institut für Mathematik in den Naturwissenschaften,
!      Leipzig, Abt. Komplexe Strukutren in Biologie und Kognition,
!      and
!    - IBM Deutschland GmbH
!
!
!    More information can be found here:
!    http://elpa.mpcdf.mpg.de/
!
!    ELPA is free software: you can redistribute it and/or modify
!    it under the terms of the version 3 of the license of the
!    GNU Lesser General Public License as published by the Free
!    Software Foundation.
!
!    ELPA is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with ELPA.  If not, see <http://www.gnu.org/licenses/>
!
!    ELPA reflects a substantial effort on the part of the original
!    ELPA consortium, and we ask you to respect the spirit of the
!    license that we chose: i.e., please contribute any changes you
!    may have back to the original ELPA library distribution, and keep
!    any derivatives of ELPA under the same license that we chose for
!    the original distribution, the GNU Lesser General Public License.
!
! This file was written by A. Marek, MPCDF

module elpa_cholesky
  use, intrinsic :: iso_c_binding
  use precision
  implicit none

  public

  public :: elpa_cholesky_a_h_a_real_double_impl       !< Cholesky factorization of a double-precision real matrix
  public :: elpa_cholesky_d_ptr_real_double_impl       !< Cholesky factorization of a double-precision real matrix

  public :: elpa_cholesky_a_h_a_complex_double_impl    !< Cholesky factorization of a double-precision complex matrix
  public :: elpa_cholesky_d_ptr_complex_double_impl    !< Cholesky factorization of a double-precision complex matrix

  public :: elpa_cholesky_a_h_a_real_single_impl       !< Cholesky factorization of a single-precision real matrix
  public :: elpa_cholesky_d_ptr_real_single_impl       !< Cholesky factorization of a single-precision real matrix


  public :: elpa_cholesky_a_h_a_complex_single_impl    !< Cholesky factorization of a single-precision complex matrix
  public :: elpa_cholesky_d_ptr_complex_single_impl    !< Cholesky factorization of a single-precision complex matrix

  contains

















!> \brief  elpa_cholesky_a_h_a_real_double_impl: Cholesky factorization of a double-precision real symmetric matrix
!> \details
!> \param  obj                    elpa_t object contains:
!> \param     - obj%na            Order of matrix
!> \param     - obj%local_nrows   Leading dimension of a
!> \param     - obj%local_ncols   local columns of matrix a
!> \param     - obj%nblk          blocksize of cyclic distribution, must be the same in both directions!
!> \param     - obj%mpi_comm_rows MPI communicator for rows
!> \param     - obj%mpi_comm_cols MPI communicator for columns
!> \param     - obj%wantDebug     logical, more debug information on failure
!> \param  a(lda,matrixCols)      Distributed matrix which should be inverted
!>                                Distribution is like in Scalapack.
!>                                Only upper triangle needs to be set.
!>                                The lower triangle is not referenced.
!> \result succes                 logical, reports success or failure
   function elpa_cholesky_a_h_a_real_double_impl (obj, a) result(success)
!    This file is part of ELPA.
!
!    The ELPA library was originally created by the ELPA consortium,
!    consisting of the following organizations:
!
!    - Max Planck Computing and Data Facility (MPCDF), formerly known as
!      Rechenzentrum Garching der Max-Planck-Gesellschaft (RZG),
!    - Bergische Universität Wuppertal, Lehrstuhl für angewandte
!      Informatik,
!    - Technische Universität München, Lehrstuhl für Informatik mit
!      Schwerpunkt Wissenschaftliches Rechnen ,
!    - Fritz-Haber-Institut, Berlin, Abt. Theorie,
!    - Max-Plack-Institut für Mathematik in den Naturwissenschaften,
!      Leipzig, Abt. Komplexe Strukutren in Biologie und Kognition,
!      and
!    - IBM Deutschland GmbH
!
!    This particular source code file contains additions, changes and
!    enhancements authored by Intel Corporation which is not part of
!    the ELPA consortium.
!
!    More information can be found here:
!    http://elpa.mpcdf.mpg.de/
!
!    ELPA is free software: you can redistribute it and/or modify
!    it under the terms of the version 3 of the license of the
!    GNU Lesser General Public License as published by the Free
!    Software Foundation.
!
!    ELPA is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with ELPA.  If not, see <http://www.gnu.org/licenses/>
!
!    ELPA reflects a substantial effort on the part of the original
!    ELPA consortium, and we ask you to respect the spirit of the
!    license that we chose: i.e., please contribute any changes you
!    may have back to the original ELPA library distribution, and keep
!    any derivatives of ELPA under the same license that we chose for
!    the original distribution, the GNU Lesser General Public License.




!cannot use "../src/cholesky/./../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)



  use elpa1_compute
  use elpa_utilities
  use elpa_mpi
  use precision
  use elpa_abstract_impl
  use elpa_omp
  use elpa_blas_interfaces
  use elpa_gpu
  use mod_check_for_gpu
  use invert_trm_gpu !, only : gpu_copy_&
                     !        &double&
                     !        &_tmp1_tmp2, &
                     !        gpu_copy_&
                     !        &double&
                     !        &_a_tmp1
  use cholesky_gpu
  use mod_query_gpu_usage
  implicit none
!    Copyright 2011, A. Marek
!
!    This file is part of ELPA.
!
!    The ELPA library was originally created by the ELPA consortium,
!    consisting of the following organizations:
!
!    - Max Planck Computing and Data Facility (MPCDF), formerly known as
!      Rechenzentrum Garching der Max-Planck-Gesellschaft (RZG),
!    - Bergische Universität Wuppertal, Lehrstuhl für angewandte
!      Informatik,
!    - Technische Universität München, Lehrstuhl für Informatik mit
!      Schwerpunkt Wissenschaftliches Rechnen ,
!    - Fritz-Haber-Institut, Berlin, Abt. Theorie,
!    - Max-Plack-Institut für Mathematik in den Naturwissenschaften,
!      Leipzig, Abt. Komplexe Strukutren in Biologie und Kognition,
!      and
!    - IBM Deutschland GmbH
!
!    This particular source code file contains additions, changes and
!    enhancements authored by Intel Corporation which is not part of
!    the ELPA consortium.
!
!    More information can be found here:
!    http://elpa.mpcdf.mpg.de/
!
!    ELPA is free software: you can redistribute it and/or modify
!    it under the terms of the version 3 of the license of the
!    GNU Lesser General Public License as published by the Free
!    Software Foundation.
!
!    ELPA is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with ELPA.  If not, see <http://www.gnu.org/licenses/>
!
!    ELPA reflects a substantial effort on the part of the original
!    ELPA consortium, and we ask you to respect the spirit of the
!    license that we chose: i.e., please contribute any changes you
!    may have back to the original ELPA library distribution, and keep
!    any derivatives of ELPA under the same license that we chose for
!    the original distribution, the GNU Lesser General Public License.
!
!
  integer, parameter :: rk = C_DOUBLE
  integer, parameter :: rck = C_DOUBLE
  real(kind=rck), parameter      :: ZERO=0.0_rk, ONE = 1.0_rk

  class(elpa_abstract_impl_t), intent(inout) :: obj
  integer(kind=ik)                           :: na, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols, mpi_comm_all
  real(kind=rck)                    :: a(obj%local_nrows,*)
  integer(kind=ik)                           :: my_prow, my_pcol, np_rows, np_cols, myid
  integer(kind=MPI_KIND)                     :: mpierr, my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI, myidMPI
  integer(kind=ik)                           :: l_cols, l_rows, l_col1, l_row1, l_colx, l_rowx
  integer(kind=ik)                           :: n, nc, i, info
  integer(kind=BLAS_KIND)                    :: infoBLAS
  integer(kind=ik)                           :: lcs, lce, lrs, lre
  integer(kind=ik)                           :: tile_size, l_rows_tile, l_cols_tile

  real(kind=rck), allocatable       :: tmp1(:), tmp2(:,:), tmatr(:,:), tmatc(:,:)
  logical                                    :: wantDebug
  logical                                    :: success
  integer(kind=ik)                           :: istat, debug, error
  character(200)                             :: errorMessage
  integer(kind=ik)                           :: nrThreads, limitThreads
  character(20)                              :: gpuString
  logical                                    :: successGPU
  logical                                    :: useGPU
  logical                                    :: useCCL
  integer(kind=c_int)                        :: numGPU
  integer(kind=c_intptr_t)                   :: num
  integer(kind=c_intptr_t)                   :: tmp1_dev, tmatc_dev, tmatr_dev, a_dev, tmp2_dev
  integer(kind=c_intptr_t)                   :: info_dev, info_abs_accumulated_dev
  type(c_ptr)                                :: tmp1_mpi_dev
  real(kind=rck), pointer           :: tmp1_mpi_fortran_ptr(:,:)
  integer(kind=c_intptr_t)                   :: a_off, tmatc_off, tmatr_off
  type(c_ptr)                                :: tmatc_mpi_dev
  real(kind=rck), pointer           :: tmatc_mpi_fortran_ptr(:,:)

  integer(kind=c_intptr_t), parameter        :: size_of_datatype = size_of_&
                                                            &double&
                                                            &_&
                                                            &real

  integer(kind=c_intptr_t)                   :: gpublasHandle, gpusolverHandle, my_stream, offset
  integer(kind=c_int)                        :: gpu_cholesky
  integer(kind=ik)                           :: blocking

  success = .true.
  useGPU = .false.
  useCCL = .false.



  ! check whether the above setting should be overriden
  if (obj%is_set("gpu_cholesky") == 1) then
    call obj%get("gpu_cholesky", gpu_cholesky, error)
    if (error .ne. ELPA_OK) then
      print *,"Problem getting option for gpu_cholesky. Aborting..."
      stop 1
    endif
    if (useGPU .and. gpu_cholesky .eq. 0) then
      useGPU = .false.
    else if (.not.(useGPU) .and. gpu_cholesky .eq. 1) then
      useGPU = .true.
    else 
    endif
  else 
    ! no override by user 
  endif

  if (.not.(useGPU)) then
  endif

  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  
  call obj%timer%start("elpa_cholesky_&
  &real&
  &_&
  &double&
  &"//gpuString)

  nrThreads=1

  na         = obj%na
  matrixRows = obj%local_nrows
  nblk       = obj%nblk
  matrixCols = obj%local_ncols

  call obj%get("debug",debug,error)
  if (error .ne. ELPA_OK) then
    write(error_unit,*) "ELPA_CHOLESKY: Problem getting option for debug settings. Aborting..."
    success = .false.
    return
  endif
  if (debug == 1) then
    wantDebug = .true.
  else
    wantDebug = .false.
  endif

  mpi_comm_all    = obj%mpi_setup%mpi_comm_parent
  mpi_comm_cols   = obj%mpi_setup%mpi_comm_cols
  mpi_comm_rows   = obj%mpi_setup%mpi_comm_rows
    
  myid    = obj%mpi_setup%myRank_comm_parent
  my_prow = obj%mpi_setup%myRank_comm_rows
  my_pcol = obj%mpi_setup%myRank_comm_cols
        
  np_rows = obj%mpi_setup%nRanks_comm_rows
  np_cols = obj%mpi_setup%nRanks_comm_cols

  success = .true.

  ! Matrix is split into tiles; work is done only for tiles on the diagonal or above

  call obj%timer%start("prepare")

  if (obj%is_set("blocking_in_cholesky") == 1) then
    call obj%get("blocking_in_cholesky", blocking, error)
    if (error .ne. ELPA_OK) then
      write(error_unit,*) "ELPA_CHOLESKY: Problem in getting keyword 'blocking_in_cholesky'. Aborting..."
      stop 1
    endif
  else
    if (useGPU) then
      blocking = 1024
    else
      blocking = 128
    endif
    call obj%set("blocking_in_cholesky", blocking, error)
    if (error .ne. ELPA_OK) then
      write(error_unit,*) "ELPA_CHOLESKY: Problem in setting keyword 'blocking_in_cholesky'. Aborting..."
      stop 1
    endif
  endif 


  tile_size = nblk*least_common_multiple(np_rows,np_cols) ! minimum global tile size
  tile_size = ((blocking*max(np_rows,np_cols)-1)/tile_size+1)*tile_size ! make local tiles at least 128 wide

  l_rows_tile = tile_size/np_rows ! local rows of a tile
  l_cols_tile = tile_size/np_cols ! local cols of a tile

  l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a
  l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a

  if (useGPU) then
    call obj%timer%start("check_for_gpu")
    if (check_for_gpu(obj, myid, numGPU)) then
      ! set the neccessary parameters
      call set_gpu_parameters()
    else
      print *,"ELPA_CHOLESKY: GPUs are requested but not detected! Aborting..."
      success = .false.
      return
    endif
    call obj%timer%stop("check_for_gpu")
  else ! useGPU
  endif ! useGPU


  if (useGPU) then

    successGPU = gpu_malloc(tmp1_dev, nblk*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmp1_dev", 335,  successGPU)

    successGPU = gpu_memset(tmp1_dev, 0, nblk*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmp1_dev", 346,  successGPU)

    successGPU = gpu_malloc(tmp2_dev, nblk*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmp2_dev", 350,  successGPU)

    successGPU = gpu_memset(tmp2_dev, 0, nblk*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmp2_dev", 361,  successGPU)

    successGPU = gpu_malloc(tmatc_dev, l_cols*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmatc_dev", 365,  successGPU)

    successGPU = gpu_memset(tmatc_dev, 0, l_cols*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmatc_dev", 376,  successGPU)

    successGPU = gpu_malloc(tmatr_dev, l_rows*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmatr_dev", 380,  successGPU)

    successGPU = gpu_memset(tmatr_dev, 0, l_rows*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmatr_dev", 391,  successGPU)

    successGPU = gpu_malloc(a_dev, matrixRows*matrixCols*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: a_dev", 396,  successGPU)
  endif ! useGPU


  allocate(tmp1(nblk*nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmp1", 441,  istat,  errorMessage)

  allocate(tmp2(nblk,nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmp2", 444,  istat,  errorMessage)



  tmp1 = 0
  tmp2 = 0

  allocate(tmatr(l_rows,nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmatr", 460,  istat,  errorMessage)

  allocate(tmatc(l_cols,nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmatc", 463,  istat,  errorMessage)


  tmatr = 0
  tmatc = 0

  if (useGPU) then
    num = matrixRows*matrixCols* size_of_datatype
    successGPU = gpu_memcpy(a_dev, int(loc(a(1,1)),kind=c_intptr_t), &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
    call check_memcpy_GPU_f("elpa_cholesky 1: memcpy a-> a_dev", 495,  successGPU)
  endif


  call obj%timer%stop("prepare")
  call obj%timer%start("loop1")
  do n = 1, na, nblk


    ! Calculate first local row and column of the still remaining matrix
    ! on the local processor

    l_row1 = local_index(n, my_prow, np_rows, nblk, +1)
    l_col1 = local_index(n, my_pcol, np_cols, nblk, +1)

    l_rowx = local_index(n+nblk, my_prow, np_rows, nblk, +1)
    l_colx = local_index(n+nblk, my_pcol, np_cols, nblk, +1)

    if (n+nblk > na) then

      ! This is the last step, just do a Cholesky-Factorization
      ! of the remaining block

      if (useGPU) then
        if (my_prow==prow(n, nblk, np_rows) .and. my_pcol==pcol(n, nblk, np_cols)) then

          call obj%timer%start("lapack")
          num = matrixRows*matrixCols* size_of_datatype
          successGPU = gpu_memcpy(int(loc(a(1,1)),kind=c_intptr_t), a_dev,  &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a", 582,  successGPU)

          call DPOTRF('U', int(na-n+1,kind=BLAS_KIND), a(l_row1,l_col1), &
                             int(matrixRows,kind=BLAS_KIND), infoBLAS )
          info = int(infoBLAS,kind=ik)

          num = matrixRows*matrixCols* size_of_datatype
          successGPU = gpu_memcpy(a_dev, int(loc(a(1,1)),kind=c_intptr_t), &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a", 600,  successGPU)

          call obj%timer%stop("lapack")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &real&

            &: Error in dpotrf: ",info
            success = .false. ! "
            return
          endif ! info
        endif ! (my_prow==prow(n, nblk, np_rows) .and. my_pcol==pcol(n, nblk, np_cols))

      else ! useGPU
        if (my_prow==prow(n, nblk, np_rows) .and. my_pcol==pcol(n, nblk, np_cols)) then
          call obj%timer%start("blas")

          call DPOTRF('U', int(na-n+1,kind=BLAS_KIND), a(l_row1,l_col1), &
                             int(matrixRows,kind=BLAS_KIND), infoBLAS )
          info = int(infoBLAS,kind=ik)
          call obj%timer%stop("blas")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &real&

            &: Error in dpotrf: ",info
            success = .false. ! "
            return
          endif

        endif
      endif ! useGPU
      
      exit ! Loop
    endif ! (n+nblk > na) 

    ! This is not the last step

    if (my_prow==prow(n, nblk, np_rows)) then

      if (my_pcol==pcol(n, nblk, np_cols)) then

        if (useGPU) then
          call obj%timer%start("lapack")

          num = matrixRows*matrixCols* size_of_datatype

          successGPU = gpu_memcpy(int(loc(a(1,1)),kind=c_intptr_t), a_dev,  &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a", 742,  successGPU)


          call DPOTRF('U', int(nblk,kind=BLAS_KIND), a(l_row1,l_col1), &
                               int(matrixRows,kind=BLAS_KIND) , infoBLAS )
          info = int(infoBLAS,kind=ik)

          num = matrixRows*matrixCols* size_of_datatype
          successGPU = gpu_memcpy(a_dev, int(loc(a(1,1)),kind=c_intptr_t), &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a", 769,  successGPU)

          call obj%timer%stop("lapack")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &real&

            &: Error in dpotrf: ",info
            success = .false. ! "
            return
          endif ! info
        else ! useGPU
          ! The process owning the upper left remaining block does the
          ! Cholesky-Factorization of this block
          call obj%timer%start("blas")

          call DPOTRF('U', int(nblk,kind=BLAS_KIND), a(l_row1,l_col1), &
                               int(matrixRows,kind=BLAS_KIND) , infoBLAS )
          info = int(infoBLAS,kind=ik)
          call obj%timer%stop("blas")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &real&

            &: Error in dpotrf 2: ",info
            success = .false. ! "
            return
          endif ! info/=0
        endif ! useGPU
 
        if (useGPU) then
          my_stream = obj%gpu_setup%my_stream
          call gpu_copy_double_a_tmp1 (a_dev, tmp1_dev, l_row1, l_col1, matrixRows, nblk, my_stream)
        else ! useGPU
          nc = 0
          do i=1,nblk
            tmp1(nc+1:nc+i) = a(l_row1:l_row1+i-1,l_col1+i-1)
            nc = nc+i
          enddo
        endif ! useGPU
      endif ! (my_pcol==pcol(n, nblk, np_cols))


      if (useGPU) then
        my_stream = obj%gpu_setup%my_stream
        call gpu_copy_double_tmp1_tmp2 (tmp1_dev, tmp2_dev, nblk, nblk, my_stream)
      else ! useGPU
        nc = 0
        do i=1,nblk
          tmp2(1:i,i) = tmp1(nc+1:nc+i)
          nc = nc+i
        enddo
      endif ! useGPU

      if (useGPU) then
        call obj%timer%start("gpublas")
        gpublasHandle = obj%gpu_setup%gpublasHandleArray(0)
        if (l_cols-l_colx+1 > 0) then
          a_off = (l_row1-1 + (l_colx-1)*matrixRows) * size_of_datatype
          call gpublas_DTRSM('L', 'U', 'T', 'N', nblk, l_cols-l_colx+1, ONE, &
                            tmp2_dev, nblk, a_dev+a_off, matrixRows, gpublasHandle)
          if (wantDebug .and. gpu_vendor() /= SYCL_GPU) successGPU = gpu_DeviceSynchronize()
        endif
        call obj%timer%stop("gpublas")

      else ! useGPU

        call obj%timer%start("lapack")
        if (l_cols-l_colx+1>0) then
          call DTRSM('L', 'U', 'T', 'N', int(nblk,kind=BLAS_KIND),  &
                            int(l_cols-l_colx+1,kind=BLAS_KIND), ONE, tmp2, &
                            int(ubound(tmp2,dim=1),kind=BLAS_KIND), a(l_row1,l_colx), int(matrixRows,kind=BLAS_KIND) )
        endif
        call obj%timer%stop("lapack")
      endif ! useGPU
    endif ! (my_prow==prow(n, nblk, np_rows))


    if (useGPU) then
      if (my_prow==prow(n, nblk, np_rows)) then
        ! if l_cols-l_colx+1 == 0 kernel launch with 0 blocks => raises error
        call obj%timer%start("a_tmatc_kernel")
        if (l_cols-l_colx+1>0) then
           my_stream = obj%gpu_setup%my_stream
           call gpu_copy_double_a_tmatc(a_dev, tmatc_dev, nblk, matrixRows, l_cols, l_colx, l_row1, my_stream)
        endif
        call obj%timer%stop("a_tmatc_kernel")
      endif
    else ! useGPU
      do i=1,nblk
        if (my_prow==prow(n, nblk, np_rows)) tmatc(l_colx:l_cols,i) = a(l_row1+i-1,l_colx:l_cols)
      enddo
    endif ! useGPU

    if (useGPU  .and. .not. useCCL) then
      num = l_cols*nblk*size_of_datatype
      successGPU = gpu_memcpy(int(loc(tmatc),kind=c_intptr_t), tmatc_dev, num, &
                              gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("elpa_cholesky: tmatc_dev to tmatc", 1201,  successGPU)
    endif


    if (useGPU .and. .not. useCCL) then

      ! not needed if transpose vec CCL
      num = l_rows*nblk*size_of_datatype
      successGPU = gpu_memcpy(int(loc(tmatr),kind=c_intptr_t), tmatr_dev, num, &
                              gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("elpa_cholesky: tmatr_dev to tmatr", 1220,  successGPU)
    endif ! (useGPU .and. .not. useCCL)


    if (useCCL) then

    else ! useCCL
      call elpa_transpose_vectors_&
      &real&
      &_&
      &double &
      (obj, tmatc, ubound(tmatc,dim=1), mpi_comm_cols, &
      tmatr, ubound(tmatr,dim=1), mpi_comm_rows, &
      n, na, nblk, nblk, nrThreads, .false., success)
    endif ! useCCL

    if (useGPU .and. .not. useCCL) then

      num = l_rows*nblk*size_of_datatype
      successGPU = gpu_memcpy(tmatr_dev, int(loc(tmatr),kind=c_intptr_t), num, &
                              gpuMemcpyHostToDevice)
      call check_memcpy_GPU_f("elpa_cholesky: tmat to tmatr_dev", 1271,  successGPU)
    endif

    if (useGPU) then
      do i=0,(na-1)/tile_size
        lcs = max(l_colx,i*l_cols_tile+1)
        lce = min(l_cols,(i+1)*l_cols_tile)
        lrs = l_rowx
        lre = min(l_rows,(i+1)*l_rows_tile)
        if (lce  < lcs .or. lre < lrs) cycle
        call obj%timer%start("gpublas")
        gpublasHandle = obj%gpu_setup%gpublasHandleArray(0)
        tmatr_off = (lrs-1 + (1-1)*l_rows) * size_of_datatype
        tmatc_off = (lcs-1 + (1-1)*l_cols) * size_of_datatype
        a_off = (lrs-1 + (lcs-1)*matrixRows) * size_of_datatype
        call gpublas_DGEMM('N', 'T', lre-lrs+1, lce-lcs+1, nblk, &
                            -ONE, tmatr_dev+tmatr_off, l_rows, tmatc_dev+tmatc_off, l_cols, ONE, &
                            a_dev+a_off, matrixRows, gpublasHandle)
        if (wantDebug .and. gpu_vendor() /= SYCL_GPU) successGPU = gpu_DeviceSynchronize()
        call obj%timer%stop("gpublas")
      enddo
    else !useGPU
      do i=0,(na-1)/tile_size
        lcs = max(l_colx,i*l_cols_tile+1)
        lce = min(l_cols,(i+1)*l_cols_tile)
        lrs = l_rowx
        lre = min(l_rows,(i+1)*l_rows_tile)
        if (lce<lcs .or. lre<lrs) cycle
        call obj%timer%start("blas")
        call DGEMM('N', 'T', int(lre-lrs+1,kind=BLAS_KIND), int(lce-lcs+1,kind=BLAS_KIND), &
                            int(nblk,kind=BLAS_KIND), -ONE,  &
                            tmatr(lrs,1), int(ubound(tmatr,dim=1),kind=BLAS_KIND), tmatc(lcs,1), &
                            int(ubound(tmatc,dim=1),kind=BLAS_KIND), &
                            ONE, a(lrs,lcs), int(matrixRows,kind=BLAS_KIND))
        call obj%timer%stop("blas")
      enddo
    endif ! useGPU


  enddo ! n = 1, na, nblk

  call obj%timer%stop("loop1")
  
  call obj%timer%start("deallocate")
  if (useGPU) then



    successGPU = gpu_free(tmp1_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmp1_dev", 1351,  successGPU)

    successGPU = gpu_free(tmp2_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmp1_dev", 1354,  successGPU)

    successGPU = gpu_free(tmatc_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmatc_dev", 1357,  successGPU)

    successGPU = gpu_free(tmatr_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmatr_dev", 1360,  successGPU)

  endif


  deallocate(tmp1, tmp2, tmatr, tmatc, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("elpa_cholesky: tmp1, tmp2, tmatr, tmatc", 1384,  istat,  errorMessage)

  call obj%timer%stop("deallocate")


  call obj%timer%start("loop2")

  call obj%timer%start("copy1")

  if (useGPU) then
    num = matrixRows*matrixCols* size_of_datatype
    successGPU = gpu_memcpy(int(loc(a(1,1)),kind=c_intptr_t), a_dev,  &
                     matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
    call check_memcpy_GPU_f("elpa_cholesky: memcpy 2 a-> a_dev", 1405,  successGPU)
  endif
  ! Set the lower triangle to 0, it contains garbage (form the above matrix multiplications)

  call obj%timer%stop("copy1")
  do i=1,na
    if (my_pcol==pcol(i, nblk, np_cols)) then
      ! column i is on local processor
      l_col1 = local_index(i  , my_pcol, np_cols, nblk, +1) ! local column number
      l_row1 = local_index(i+1, my_prow, np_rows, nblk, +1) ! first row below diagonal
      if (useGPU) then
        call obj%timer%start("copy2")
!#ifndef DEVICE_POINTER
        ! set to zero maybe a kerne would be better here
          offset = (l_row1-1 + matrixRows * (l_col1-1)) * size_of_datatype
          num = (l_rows-l_row1+1)*size_of_datatype
          successGPU = gpu_memset(a_dev+offset, 0, num)
          call check_memcpy_GPU_f("elpa_cholesky: memset a_dev", 1450,  successGPU)
!#else 
!          !offset = (l_row1-1 + matrixRows * (l_col1-1)) * size_of_datatype
!          !num = (l_rows-l_row1+1)*size_of_datatype
!#ifdef WITH_GPU_STREAMS
!          !my_stream = obj%gpu_setup%my_stream
!          !successGPU = gpu_memset_async(a_dev+offset, 0, num, my_stream)
!          !call check_memcpy_GPU_f("elpa_cholesky: memset a_dev", 1458,  successGPU)
!
!          !successGPU = gpu_stream_synchronize(my_stream)
!          !call check_memcpy_GPU_f("elpa_cholesky: memset", 1461,  successGPU)
!#else
!          !successGPU = gpu_memset(a_dev+offset, 0, num)
!          !call check_memcpy_GPU_f("elpa_cholesky: memset tmp1_dev", 1464,  successGPU)
!#endif
!
!
!!#if defined(WITH_NVIDIA_CUSOLVER) || defined(WITH_AMD_ROCSOLVER)
!!        successGPU = gpu_memcpy(int(loc(a_tmp(1,1)),kind=c_intptr_t), a_dev, &
!!                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
!!        call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a_tmp", 1471,  successGPU)
!!#endif
!!        a_tmp(l_row1:l_rows,l_col1) = 0
!!#if defined(WITH_NVIDIA_CUSOLVER) || defined(WITH_AMD_ROCSOLVER)
!!        successGPU = gpu_memcpy(a_dev, int(loc(a_tmp(1,1)),kind=c_intptr_t), &
!!                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
!!        call check_memcpy_GPU_f("elpa_cholesky: memcpy a_tmp-> a_dev", 1477,  successGPU)
!#endif

!#endif 

  call obj%timer%stop("copy2")

      else ! useGPU
        a(l_row1:l_rows,l_col1) = 0
      endif ! useGPU
    endif ! (my_pcol==pcol(i, nblk, np_cols))
  enddo
  !endif ! useGPU

  call obj%timer%stop("loop2")


  call obj%timer%start("cleanup")
  if (useGPU) then
    ! copy back
    successGPU = gpu_memcpy(int(loc(a(1,1)),kind=c_intptr_t), a_dev,  &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
    call check_memcpy_GPU_f("elpa_cholesky: memcpy a-> d_dev", 1502,  successGPU)

    successGPU = gpu_free(a_dev)
    call check_dealloc_GPU_f("elpa_cholesky: a_dev", 1505,  successGPU)


  endif

  ! restore original OpenMP settings
  call obj%timer%stop("cleanup")
  call obj%timer%stop("elpa_cholesky_&
  &real&
  &_&
  &double&
  &"//gpuString)

    end function elpa_cholesky_a_h_a_real_double_impl



















!> \brief  elpa_cholesky_d_ptr_real_double_impl: Cholesky factorization of a double-precision real symmetric matrix
!> \details
!> \param  obj                    elpa_t object contains:
!> \param     - obj%na            Order of matrix
!> \param     - obj%local_nrows   Leading dimension of a
!> \param     - obj%local_ncols   local columns of matrix a
!> \param     - obj%nblk          blocksize of cyclic distribution, must be the same in both directions!
!> \param     - obj%mpi_comm_rows MPI communicator for rows
!> \param     - obj%mpi_comm_cols MPI communicator for columns
!> \param     - obj%wantDebug     logical, more debug information on failure
!> \param  aDev(lda,matrixCols)   Distributed matrix which should be inverted as type(c_ptr) living on device
!>                                Distribution is like in Scalapack.
!>                                Only upper triangle needs to be set.
!>                                The lower triangle is not referenced.
!> \result succes                 logical, reports success or failure
   function elpa_cholesky_d_ptr_real_double_impl (obj, aDev) result(success)
!    This file is part of ELPA.
!
!    The ELPA library was originally created by the ELPA consortium,
!    consisting of the following organizations:
!
!    - Max Planck Computing and Data Facility (MPCDF), formerly known as
!      Rechenzentrum Garching der Max-Planck-Gesellschaft (RZG),
!    - Bergische Universität Wuppertal, Lehrstuhl für angewandte
!      Informatik,
!    - Technische Universität München, Lehrstuhl für Informatik mit
!      Schwerpunkt Wissenschaftliches Rechnen ,
!    - Fritz-Haber-Institut, Berlin, Abt. Theorie,
!    - Max-Plack-Institut für Mathematik in den Naturwissenschaften,
!      Leipzig, Abt. Komplexe Strukutren in Biologie und Kognition,
!      and
!    - IBM Deutschland GmbH
!
!    This particular source code file contains additions, changes and
!    enhancements authored by Intel Corporation which is not part of
!    the ELPA consortium.
!
!    More information can be found here:
!    http://elpa.mpcdf.mpg.de/
!
!    ELPA is free software: you can redistribute it and/or modify
!    it under the terms of the version 3 of the license of the
!    GNU Lesser General Public License as published by the Free
!    Software Foundation.
!
!    ELPA is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with ELPA.  If not, see <http://www.gnu.org/licenses/>
!
!    ELPA reflects a substantial effort on the part of the original
!    ELPA consortium, and we ask you to respect the spirit of the
!    license that we chose: i.e., please contribute any changes you
!    may have back to the original ELPA library distribution, and keep
!    any derivatives of ELPA under the same license that we chose for
!    the original distribution, the GNU Lesser General Public License.




!cannot use "../src/cholesky/./../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)



  use elpa1_compute
  use elpa_utilities
  use elpa_mpi
  use precision
  use elpa_abstract_impl
  use elpa_omp
  use elpa_blas_interfaces
  use elpa_gpu
  use mod_check_for_gpu
  use invert_trm_gpu !, only : gpu_copy_&
                     !        &double&
                     !        &_tmp1_tmp2, &
                     !        gpu_copy_&
                     !        &double&
                     !        &_a_tmp1
  use cholesky_gpu
  use mod_query_gpu_usage
  implicit none
!    Copyright 2011, A. Marek
!
!    This file is part of ELPA.
!
!    The ELPA library was originally created by the ELPA consortium,
!    consisting of the following organizations:
!
!    - Max Planck Computing and Data Facility (MPCDF), formerly known as
!      Rechenzentrum Garching der Max-Planck-Gesellschaft (RZG),
!    - Bergische Universität Wuppertal, Lehrstuhl für angewandte
!      Informatik,
!    - Technische Universität München, Lehrstuhl für Informatik mit
!      Schwerpunkt Wissenschaftliches Rechnen ,
!    - Fritz-Haber-Institut, Berlin, Abt. Theorie,
!    - Max-Plack-Institut für Mathematik in den Naturwissenschaften,
!      Leipzig, Abt. Komplexe Strukutren in Biologie und Kognition,
!      and
!    - IBM Deutschland GmbH
!
!    This particular source code file contains additions, changes and
!    enhancements authored by Intel Corporation which is not part of
!    the ELPA consortium.
!
!    More information can be found here:
!    http://elpa.mpcdf.mpg.de/
!
!    ELPA is free software: you can redistribute it and/or modify
!    it under the terms of the version 3 of the license of the
!    GNU Lesser General Public License as published by the Free
!    Software Foundation.
!
!    ELPA is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with ELPA.  If not, see <http://www.gnu.org/licenses/>
!
!    ELPA reflects a substantial effort on the part of the original
!    ELPA consortium, and we ask you to respect the spirit of the
!    license that we chose: i.e., please contribute any changes you
!    may have back to the original ELPA library distribution, and keep
!    any derivatives of ELPA under the same license that we chose for
!    the original distribution, the GNU Lesser General Public License.
!
!
  integer, parameter :: rk = C_DOUBLE
  integer, parameter :: rck = C_DOUBLE
  real(kind=rck), parameter      :: ZERO=0.0_rk, ONE = 1.0_rk

  class(elpa_abstract_impl_t), intent(inout) :: obj
  integer(kind=ik)                           :: na, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols, mpi_comm_all
  type(c_ptr)                                :: aDev
!#if !defined(WITH_NVIDIA_CUSOLVER) && !defined(WITH_AMD_ROCSOLVER)
  real(kind=rck), allocatable       :: a_tmp(:,:)
!#endif
  integer(kind=ik)                           :: my_prow, my_pcol, np_rows, np_cols, myid
  integer(kind=MPI_KIND)                     :: mpierr, my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI, myidMPI
  integer(kind=ik)                           :: l_cols, l_rows, l_col1, l_row1, l_colx, l_rowx
  integer(kind=ik)                           :: n, nc, i, info
  integer(kind=BLAS_KIND)                    :: infoBLAS
  integer(kind=ik)                           :: lcs, lce, lrs, lre
  integer(kind=ik)                           :: tile_size, l_rows_tile, l_cols_tile

  real(kind=rck), allocatable       :: tmp1(:), tmp2(:,:), tmatr(:,:), tmatc(:,:)
  logical                                    :: wantDebug
  logical                                    :: success
  integer(kind=ik)                           :: istat, debug, error
  character(200)                             :: errorMessage
  integer(kind=ik)                           :: nrThreads, limitThreads
  character(20)                              :: gpuString
  logical                                    :: successGPU
  logical                                    :: useGPU
  logical                                    :: useCCL
  integer(kind=c_int)                        :: numGPU
  integer(kind=c_intptr_t)                   :: num
  integer(kind=c_intptr_t)                   :: tmp1_dev, tmatc_dev, tmatr_dev, a_dev, tmp2_dev
  integer(kind=c_intptr_t)                   :: info_dev, info_abs_accumulated_dev
  type(c_ptr)                                :: tmp1_mpi_dev
  real(kind=rck), pointer           :: tmp1_mpi_fortran_ptr(:,:)
  integer(kind=c_intptr_t)                   :: a_off, tmatc_off, tmatr_off
  type(c_ptr)                                :: tmatc_mpi_dev
  real(kind=rck), pointer           :: tmatc_mpi_fortran_ptr(:,:)

  integer(kind=c_intptr_t), parameter        :: size_of_datatype = size_of_&
                                                            &double&
                                                            &_&
                                                            &real

  integer(kind=c_intptr_t)                   :: gpublasHandle, gpusolverHandle, my_stream, offset
  integer(kind=c_int)                        :: gpu_cholesky
  integer(kind=ik)                           :: blocking

  success = .true.
  useGPU = .false.
  useCCL = .false.

  useGPU = .true.

  if (.not.(useGPU)) then
    print *,"You used the interface for device pointers for elpa_cholesky but did not specify GPU usage!. Aborting..."
    stop 1
  endif

  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  
  call obj%timer%start("elpa_cholesky_&
  &real&
  &_&
  &double&
  &"//gpuString)

  nrThreads=1

  na         = obj%na
  matrixRows = obj%local_nrows
  nblk       = obj%nblk
  matrixCols = obj%local_ncols

  call obj%get("debug",debug,error)
  if (error .ne. ELPA_OK) then
    write(error_unit,*) "ELPA_CHOLESKY: Problem getting option for debug settings. Aborting..."
    success = .false.
    return
  endif
  if (debug == 1) then
    wantDebug = .true.
  else
    wantDebug = .false.
  endif

  mpi_comm_all    = obj%mpi_setup%mpi_comm_parent
  mpi_comm_cols   = obj%mpi_setup%mpi_comm_cols
  mpi_comm_rows   = obj%mpi_setup%mpi_comm_rows
    
  myid    = obj%mpi_setup%myRank_comm_parent
  my_prow = obj%mpi_setup%myRank_comm_rows
  my_pcol = obj%mpi_setup%myRank_comm_cols
        
  np_rows = obj%mpi_setup%nRanks_comm_rows
  np_cols = obj%mpi_setup%nRanks_comm_cols

  success = .true.

  ! Matrix is split into tiles; work is done only for tiles on the diagonal or above

  call obj%timer%start("prepare")

  if (obj%is_set("blocking_in_cholesky") == 1) then
    call obj%get("blocking_in_cholesky", blocking, error)
    if (error .ne. ELPA_OK) then
      write(error_unit,*) "ELPA_CHOLESKY: Problem in getting keyword 'blocking_in_cholesky'. Aborting..."
      stop 1
    endif
  else
    if (useGPU) then
      blocking = 1024
    else
      blocking = 128
    endif
    call obj%set("blocking_in_cholesky", blocking, error)
    if (error .ne. ELPA_OK) then
      write(error_unit,*) "ELPA_CHOLESKY: Problem in setting keyword 'blocking_in_cholesky'. Aborting..."
      stop 1
    endif
  endif 


  tile_size = nblk*least_common_multiple(np_rows,np_cols) ! minimum global tile size
  tile_size = ((blocking*max(np_rows,np_cols)-1)/tile_size+1)*tile_size ! make local tiles at least 128 wide

  l_rows_tile = tile_size/np_rows ! local rows of a tile
  l_cols_tile = tile_size/np_cols ! local cols of a tile

  l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a
  l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a

  if (useGPU) then
    call obj%timer%start("check_for_gpu")
    if (check_for_gpu(obj, myid, numGPU)) then
      ! set the neccessary parameters
      call set_gpu_parameters()
    else
      print *,"ELPA_CHOLESKY: GPUs are requested but not detected! Aborting..."
      success = .false.
      return
    endif
    call obj%timer%stop("check_for_gpu")
  else ! useGPU
  endif ! useGPU


  if (useGPU) then

    successGPU = gpu_malloc(tmp1_dev, nblk*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmp1_dev", 335,  successGPU)

    successGPU = gpu_memset(tmp1_dev, 0, nblk*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmp1_dev", 346,  successGPU)

    successGPU = gpu_malloc(tmp2_dev, nblk*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmp2_dev", 350,  successGPU)

    successGPU = gpu_memset(tmp2_dev, 0, nblk*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmp2_dev", 361,  successGPU)

    successGPU = gpu_malloc(tmatc_dev, l_cols*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmatc_dev", 365,  successGPU)

    successGPU = gpu_memset(tmatc_dev, 0, l_cols*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmatc_dev", 376,  successGPU)

    successGPU = gpu_malloc(tmatr_dev, l_rows*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmatr_dev", 380,  successGPU)

    successGPU = gpu_memset(tmatr_dev, 0, l_rows*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmatr_dev", 391,  successGPU)

    a_dev = transfer(aDev, a_dev)
!#if !defined(WITH_NVIDIA_CUSOLVER) && !defined(WITH_AMD_ROCSOLVER)
    allocate(a_tmp(obj%local_nrows,obj%local_ncols), stat=istat, errmsg=errorMessage)
    call check_allocate_f("elpa_cholesky: a_tmp", 407,  istat,  errorMessage)
!#endif 
  endif ! useGPU


  allocate(tmp1(nblk*nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmp1", 441,  istat,  errorMessage)

  allocate(tmp2(nblk,nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmp2", 444,  istat,  errorMessage)



  tmp1 = 0
  tmp2 = 0

  allocate(tmatr(l_rows,nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmatr", 460,  istat,  errorMessage)

  allocate(tmatc(l_cols,nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmatc", 463,  istat,  errorMessage)


  tmatr = 0
  tmatc = 0



  call obj%timer%stop("prepare")
  call obj%timer%start("loop1")
  do n = 1, na, nblk


    ! Calculate first local row and column of the still remaining matrix
    ! on the local processor

    l_row1 = local_index(n, my_prow, np_rows, nblk, +1)
    l_col1 = local_index(n, my_pcol, np_cols, nblk, +1)

    l_rowx = local_index(n+nblk, my_prow, np_rows, nblk, +1)
    l_colx = local_index(n+nblk, my_pcol, np_cols, nblk, +1)

    if (n+nblk > na) then

      ! This is the last step, just do a Cholesky-Factorization
      ! of the remaining block

      if (useGPU) then
        if (my_prow==prow(n, nblk, np_rows) .and. my_pcol==pcol(n, nblk, np_cols)) then

          call obj%timer%start("lapack")
          num = matrixRows*matrixCols* size_of_datatype

          successGPU = gpu_memcpy(int(loc(a_tmp(1,1)),kind=c_intptr_t), a_dev,  &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a_tmp", 615,  successGPU)

          call DPOTRF('U', int(na-n+1,kind=BLAS_KIND), a_tmp(l_row1,l_col1), &
                             int(matrixRows,kind=BLAS_KIND), infoBLAS )
          info = int(infoBLAS,kind=ik)

          num = matrixRows*matrixCols* size_of_datatype
          successGPU = gpu_memcpy(a_dev, int(loc(a_tmp(1,1)),kind=c_intptr_t), &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a_tmp", 633,  successGPU)

          call obj%timer%stop("lapack")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &real&

            &: Error in dpotrf: ",info
            success = .false. ! "
            return
          endif ! info
        endif ! (my_prow==prow(n, nblk, np_rows) .and. my_pcol==pcol(n, nblk, np_cols))

      else ! useGPU
        if (my_prow==prow(n, nblk, np_rows) .and. my_pcol==pcol(n, nblk, np_cols)) then
          call obj%timer%start("blas")

          info = int(infoBLAS,kind=ik)
          call obj%timer%stop("blas")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &real&

            &: Error in dpotrf: ",info
            success = .false. ! "
            return
          endif

        endif
      endif ! useGPU
      
      exit ! Loop
    endif ! (n+nblk > na) 

    ! This is not the last step

    if (my_prow==prow(n, nblk, np_rows)) then

      if (my_pcol==pcol(n, nblk, np_cols)) then

        if (useGPU) then
          call obj%timer%start("lapack")

          num = matrixRows*matrixCols* size_of_datatype

          successGPU = gpu_memcpy(int(loc(a_tmp(1,1)),kind=c_intptr_t), a_dev,  &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a_tmp", 783,  successGPU)

          call DPOTRF('U', int(nblk,kind=BLAS_KIND), a_tmp(l_row1,l_col1), &
                               int(matrixRows,kind=BLAS_KIND) , infoBLAS )
          info = int(infoBLAS,kind=ik)
          num = matrixRows*matrixCols* size_of_datatype
          successGPU = gpu_memcpy(a_dev, int(loc(a_tmp(1,1)),kind=c_intptr_t), &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a_tmp", 799,  successGPU)

          call obj%timer%stop("lapack")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &real&

            &: Error in dpotrf: ",info
            success = .false. ! "
            return
          endif ! info
        else ! useGPU
          ! The process owning the upper left remaining block does the
          ! Cholesky-Factorization of this block
          call obj%timer%start("blas")

          info = int(infoBLAS,kind=ik)
          call obj%timer%stop("blas")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &real&

            &: Error in dpotrf 2: ",info
            success = .false. ! "
            return
          endif ! info/=0
        endif ! useGPU
 
        if (useGPU) then
          my_stream = obj%gpu_setup%my_stream
          call gpu_copy_double_a_tmp1 (a_dev, tmp1_dev, l_row1, l_col1, matrixRows, nblk, my_stream)
        else ! useGPU
          nc = 0
          do i=1,nblk
            nc = nc+i
          enddo
        endif ! useGPU
      endif ! (my_pcol==pcol(n, nblk, np_cols))


      if (useGPU) then
        my_stream = obj%gpu_setup%my_stream
        call gpu_copy_double_tmp1_tmp2 (tmp1_dev, tmp2_dev, nblk, nblk, my_stream)
      else ! useGPU
        nc = 0
        do i=1,nblk
          tmp2(1:i,i) = tmp1(nc+1:nc+i)
          nc = nc+i
        enddo
      endif ! useGPU

      if (useGPU) then
        call obj%timer%start("gpublas")
        gpublasHandle = obj%gpu_setup%gpublasHandleArray(0)
        if (l_cols-l_colx+1 > 0) then
          a_off = (l_row1-1 + (l_colx-1)*matrixRows) * size_of_datatype
          call gpublas_DTRSM('L', 'U', 'T', 'N', nblk, l_cols-l_colx+1, ONE, &
                            tmp2_dev, nblk, a_dev+a_off, matrixRows, gpublasHandle)
          if (wantDebug .and. gpu_vendor() /= SYCL_GPU) successGPU = gpu_DeviceSynchronize()
        endif
        call obj%timer%stop("gpublas")

      else ! useGPU

        call obj%timer%start("lapack")
        call obj%timer%stop("lapack")
      endif ! useGPU
    endif ! (my_prow==prow(n, nblk, np_rows))


    if (useGPU) then
      if (my_prow==prow(n, nblk, np_rows)) then
        ! if l_cols-l_colx+1 == 0 kernel launch with 0 blocks => raises error
        call obj%timer%start("a_tmatc_kernel")
        if (l_cols-l_colx+1>0) then
           my_stream = obj%gpu_setup%my_stream
           call gpu_copy_double_a_tmatc(a_dev, tmatc_dev, nblk, matrixRows, l_cols, l_colx, l_row1, my_stream)
        endif
        call obj%timer%stop("a_tmatc_kernel")
      endif
    else ! useGPU
      do i=1,nblk
      enddo
    endif ! useGPU

    if (useGPU  .and. .not. useCCL) then
      num = l_cols*nblk*size_of_datatype
      successGPU = gpu_memcpy(int(loc(tmatc),kind=c_intptr_t), tmatc_dev, num, &
                              gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("elpa_cholesky: tmatc_dev to tmatc", 1201,  successGPU)
    endif


    if (useGPU .and. .not. useCCL) then

      ! not needed if transpose vec CCL
      num = l_rows*nblk*size_of_datatype
      successGPU = gpu_memcpy(int(loc(tmatr),kind=c_intptr_t), tmatr_dev, num, &
                              gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("elpa_cholesky: tmatr_dev to tmatr", 1220,  successGPU)
    endif ! (useGPU .and. .not. useCCL)


    if (useCCL) then

    else ! useCCL
      call elpa_transpose_vectors_&
      &real&
      &_&
      &double &
      (obj, tmatc, ubound(tmatc,dim=1), mpi_comm_cols, &
      tmatr, ubound(tmatr,dim=1), mpi_comm_rows, &
      n, na, nblk, nblk, nrThreads, .false., success)
    endif ! useCCL

    if (useGPU .and. .not. useCCL) then

      num = l_rows*nblk*size_of_datatype
      successGPU = gpu_memcpy(tmatr_dev, int(loc(tmatr),kind=c_intptr_t), num, &
                              gpuMemcpyHostToDevice)
      call check_memcpy_GPU_f("elpa_cholesky: tmat to tmatr_dev", 1271,  successGPU)
    endif

    if (useGPU) then
      do i=0,(na-1)/tile_size
        lcs = max(l_colx,i*l_cols_tile+1)
        lce = min(l_cols,(i+1)*l_cols_tile)
        lrs = l_rowx
        lre = min(l_rows,(i+1)*l_rows_tile)
        if (lce  < lcs .or. lre < lrs) cycle
        call obj%timer%start("gpublas")
        gpublasHandle = obj%gpu_setup%gpublasHandleArray(0)
        tmatr_off = (lrs-1 + (1-1)*l_rows) * size_of_datatype
        tmatc_off = (lcs-1 + (1-1)*l_cols) * size_of_datatype
        a_off = (lrs-1 + (lcs-1)*matrixRows) * size_of_datatype
        call gpublas_DGEMM('N', 'T', lre-lrs+1, lce-lcs+1, nblk, &
                            -ONE, tmatr_dev+tmatr_off, l_rows, tmatc_dev+tmatc_off, l_cols, ONE, &
                            a_dev+a_off, matrixRows, gpublasHandle)
        if (wantDebug .and. gpu_vendor() /= SYCL_GPU) successGPU = gpu_DeviceSynchronize()
        call obj%timer%stop("gpublas")
      enddo
    else !useGPU
      do i=0,(na-1)/tile_size
        lcs = max(l_colx,i*l_cols_tile+1)
        lce = min(l_cols,(i+1)*l_cols_tile)
        lrs = l_rowx
        lre = min(l_rows,(i+1)*l_rows_tile)
        if (lce<lcs .or. lre<lrs) cycle
        call obj%timer%start("blas")
        call obj%timer%stop("blas")
      enddo
    endif ! useGPU


  enddo ! n = 1, na, nblk

  call obj%timer%stop("loop1")
  
  call obj%timer%start("deallocate")
  if (useGPU) then



    successGPU = gpu_free(tmp1_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmp1_dev", 1351,  successGPU)

    successGPU = gpu_free(tmp2_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmp1_dev", 1354,  successGPU)

    successGPU = gpu_free(tmatc_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmatc_dev", 1357,  successGPU)

    successGPU = gpu_free(tmatr_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmatr_dev", 1360,  successGPU)

  endif


  deallocate(tmp1, tmp2, tmatr, tmatc, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("elpa_cholesky: tmp1, tmp2, tmatr, tmatc", 1384,  istat,  errorMessage)

  call obj%timer%stop("deallocate")


  call obj%timer%start("loop2")

  call obj%timer%start("copy1")

!  if (useGPU) then
!          ! is this needed: guess not
!#if !defined(WITH_NVIDIA_CUSOLVER) && !defined(WITH_AMD_ROCSOLVER)
!    num = matrixRows*matrixCols* size_of_datatype
!#ifdef WITH_GPU_STREAMS 
!    my_stream = obj%gpu_setup%my_stream
!    call gpu_memcpy_async_and_stream_synchronize &
!    ("elpa_choleksky: a_dev -> a_tmp", a_dev, 0_c_intptr_t, &
!                                    a_tmp(1:obj%local_nrows,1:obj%local_ncols), &
!                                    1, 1, num, gpuMemcpyDeviceToHost, my_stream, .false., .true., .false.)
!#else
!    successGPU = gpu_memcpy(int(loc(a_tmp(1,1)),kind=c_intptr_t), a_dev,  &
!                     matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
!    call check_memcpy_GPU_f("elpa_cholesky: memcpy 2 a-> a_dev", 1422,  successGPU)
!#endif
!#endif 
!  endif
  ! Set the lower triangle to 0, it contains garbage (form the above matrix multiplications)

  call obj%timer%stop("copy1")
  do i=1,na
    if (my_pcol==pcol(i, nblk, np_cols)) then
      ! column i is on local processor
      l_col1 = local_index(i  , my_pcol, np_cols, nblk, +1) ! local column number
      l_row1 = local_index(i+1, my_prow, np_rows, nblk, +1) ! first row below diagonal
      if (useGPU) then
        call obj%timer%start("copy2")
!#ifndef 
        ! set to zero maybe a kerne would be better here
          offset = (l_row1-1 + matrixRows * (l_col1-1)) * size_of_datatype
          num = (l_rows-l_row1+1)*size_of_datatype
          successGPU = gpu_memset(a_dev+offset, 0, num)
          call check_memcpy_GPU_f("elpa_cholesky: memset a_dev", 1450,  successGPU)
!#else 
!          !offset = (l_row1-1 + matrixRows * (l_col1-1)) * size_of_datatype
!          !num = (l_rows-l_row1+1)*size_of_datatype
!#ifdef WITH_GPU_STREAMS
!          !my_stream = obj%gpu_setup%my_stream
!          !successGPU = gpu_memset_async(a_dev+offset, 0, num, my_stream)
!          !call check_memcpy_GPU_f("elpa_cholesky: memset a_dev", 1458,  successGPU)
!
!          !successGPU = gpu_stream_synchronize(my_stream)
!          !call check_memcpy_GPU_f("elpa_cholesky: memset", 1461,  successGPU)
!#else
!          !successGPU = gpu_memset(a_dev+offset, 0, num)
!          !call check_memcpy_GPU_f("elpa_cholesky: memset tmp1_dev", 1464,  successGPU)
!#endif
!
!
!!#if defined(WITH_NVIDIA_CUSOLVER) || defined(WITH_AMD_ROCSOLVER)
!!        successGPU = gpu_memcpy(int(loc(a_tmp(1,1)),kind=c_intptr_t), a_dev, &
!!                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
!!        call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a_tmp", 1471,  successGPU)
!!#endif
!!        a_tmp(l_row1:l_rows,l_col1) = 0
!!#if defined(WITH_NVIDIA_CUSOLVER) || defined(WITH_AMD_ROCSOLVER)
!!        successGPU = gpu_memcpy(a_dev, int(loc(a_tmp(1,1)),kind=c_intptr_t), &
!!                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
!!        call check_memcpy_GPU_f("elpa_cholesky: memcpy a_tmp-> a_dev", 1477,  successGPU)
!#endif

!#endif 

  call obj%timer%stop("copy2")

      else ! useGPU
      endif ! useGPU
    endif ! (my_pcol==pcol(i, nblk, np_cols))
  enddo
  !endif ! useGPU

  call obj%timer%stop("loop2")


  call obj%timer%start("cleanup")
!#if !defined(WITH_NVIDIA_CUSOLVER) && !defined(WITH_AMD_ROCSOLVER)
!    num = matrixRows*matrixCols* size_of_datatype
!#ifdef WITH_GPU_STREAMS
!    my_stream = obj%gpu_setup%my_stream
!    call gpu_memcpy_async_and_stream_synchronize &
!    ("elpa_choleksky: a_tmp -> a_dev", a_dev, 0_c_intptr_t, &
!                                    a_tmp(1:obj%local_nrows,1:obj%local_ncols), &
!                                    1, 1, num, gpuMemcpyHostToDevice, my_stream, .false., .true., .false.)
!#else
!    successGPU = gpu_memcpy(a_dev, int(loc(a_tmp(1,1)),kind=c_intptr_t), &
!                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
!    call check_memcpy_GPU_f("elpa_cholesky: memcpy a_tmp-> a_dev", 1525,  successGPU)
!#endif
!#endif 

!#if !defined(WITH_NVIDIA_CUSOLVER) && !defined(WITH_AMD_ROCSOLVER)

    deallocate(a_tmp, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("elpa_cholesky: a_tmp", 1536,  istat,  errorMessage)
!#endif   

  ! restore original OpenMP settings
  call obj%timer%stop("cleanup")
  call obj%timer%stop("elpa_cholesky_&
  &real&
  &_&
  &double&
  &"//gpuString)

    end function elpa_cholesky_d_ptr_real_double_impl






















!> \brief  elpa_cholesky_a_h_a_real_single_impl: Cholesky factorization of a double-precision real symmetric matrix
!> \details
!> \param  obj                    elpa_t object contains:
!> \param     - obj%na            Order of matrix
!> \param     - obj%local_nrows   Leading dimension of a
!> \param     - obj%local_ncols   local columns of matrix a
!> \param     - obj%nblk          blocksize of cyclic distribution, must be the same in both directions!
!> \param     - obj%mpi_comm_rows MPI communicator for rows
!> \param     - obj%mpi_comm_cols MPI communicator for columns
!> \param     - obj%wantDebug     logical, more debug information on failure
!> \param  a(lda,matrixCols)      Distributed matrix which should be inverted
!>                                Distribution is like in Scalapack.
!>                                Only upper triangle needs to be set.
!>                                The lower triangle is not referenced.
!> \result succes                 logical, reports success or failure
   function elpa_cholesky_a_h_a_real_single_impl(obj, a) result(success)
!    This file is part of ELPA.
!
!    The ELPA library was originally created by the ELPA consortium,
!    consisting of the following organizations:
!
!    - Max Planck Computing and Data Facility (MPCDF), formerly known as
!      Rechenzentrum Garching der Max-Planck-Gesellschaft (RZG),
!    - Bergische Universität Wuppertal, Lehrstuhl für angewandte
!      Informatik,
!    - Technische Universität München, Lehrstuhl für Informatik mit
!      Schwerpunkt Wissenschaftliches Rechnen ,
!    - Fritz-Haber-Institut, Berlin, Abt. Theorie,
!    - Max-Plack-Institut für Mathematik in den Naturwissenschaften,
!      Leipzig, Abt. Komplexe Strukutren in Biologie und Kognition,
!      and
!    - IBM Deutschland GmbH
!
!    This particular source code file contains additions, changes and
!    enhancements authored by Intel Corporation which is not part of
!    the ELPA consortium.
!
!    More information can be found here:
!    http://elpa.mpcdf.mpg.de/
!
!    ELPA is free software: you can redistribute it and/or modify
!    it under the terms of the version 3 of the license of the
!    GNU Lesser General Public License as published by the Free
!    Software Foundation.
!
!    ELPA is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with ELPA.  If not, see <http://www.gnu.org/licenses/>
!
!    ELPA reflects a substantial effort on the part of the original
!    ELPA consortium, and we ask you to respect the spirit of the
!    license that we chose: i.e., please contribute any changes you
!    may have back to the original ELPA library distribution, and keep
!    any derivatives of ELPA under the same license that we chose for
!    the original distribution, the GNU Lesser General Public License.




!cannot use "../src/cholesky/./../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)



  use elpa1_compute
  use elpa_utilities
  use elpa_mpi
  use precision
  use elpa_abstract_impl
  use elpa_omp
  use elpa_blas_interfaces
  use elpa_gpu
  use mod_check_for_gpu
  use invert_trm_gpu !, only : gpu_copy_&
                     !        &single&
                     !        &_tmp1_tmp2, &
                     !        gpu_copy_&
                     !        &single&
                     !        &_a_tmp1
  use cholesky_gpu
  use mod_query_gpu_usage
  implicit none
!    Copyright 2011, A. Marek
!
!    This file is part of ELPA.
!
!    The ELPA library was originally created by the ELPA consortium,
!    consisting of the following organizations:
!
!    - Max Planck Computing and Data Facility (MPCDF), formerly known as
!      Rechenzentrum Garching der Max-Planck-Gesellschaft (RZG),
!    - Bergische Universität Wuppertal, Lehrstuhl für angewandte
!      Informatik,
!    - Technische Universität München, Lehrstuhl für Informatik mit
!      Schwerpunkt Wissenschaftliches Rechnen ,
!    - Fritz-Haber-Institut, Berlin, Abt. Theorie,
!    - Max-Plack-Institut für Mathematik in den Naturwissenschaften,
!      Leipzig, Abt. Komplexe Strukutren in Biologie und Kognition,
!      and
!    - IBM Deutschland GmbH
!
!    This particular source code file contains additions, changes and
!    enhancements authored by Intel Corporation which is not part of
!    the ELPA consortium.
!
!    More information can be found here:
!    http://elpa.mpcdf.mpg.de/
!
!    ELPA is free software: you can redistribute it and/or modify
!    it under the terms of the version 3 of the license of the
!    GNU Lesser General Public License as published by the Free
!    Software Foundation.
!
!    ELPA is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with ELPA.  If not, see <http://www.gnu.org/licenses/>
!
!    ELPA reflects a substantial effort on the part of the original
!    ELPA consortium, and we ask you to respect the spirit of the
!    license that we chose: i.e., please contribute any changes you
!    may have back to the original ELPA library distribution, and keep
!    any derivatives of ELPA under the same license that we chose for
!    the original distribution, the GNU Lesser General Public License.
!
!
  integer, parameter :: rk = C_FLOAT
  integer, parameter :: rck = C_FLOAT
  real(kind=rck), parameter      :: ZERO=0.0_rk, ONE = 1.0_rk

  class(elpa_abstract_impl_t), intent(inout) :: obj
  integer(kind=ik)                           :: na, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols, mpi_comm_all
  real(kind=rck)                    :: a(obj%local_nrows,*)
  integer(kind=ik)                           :: my_prow, my_pcol, np_rows, np_cols, myid
  integer(kind=MPI_KIND)                     :: mpierr, my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI, myidMPI
  integer(kind=ik)                           :: l_cols, l_rows, l_col1, l_row1, l_colx, l_rowx
  integer(kind=ik)                           :: n, nc, i, info
  integer(kind=BLAS_KIND)                    :: infoBLAS
  integer(kind=ik)                           :: lcs, lce, lrs, lre
  integer(kind=ik)                           :: tile_size, l_rows_tile, l_cols_tile

  real(kind=rck), allocatable       :: tmp1(:), tmp2(:,:), tmatr(:,:), tmatc(:,:)
  logical                                    :: wantDebug
  logical                                    :: success
  integer(kind=ik)                           :: istat, debug, error
  character(200)                             :: errorMessage
  integer(kind=ik)                           :: nrThreads, limitThreads
  character(20)                              :: gpuString
  logical                                    :: successGPU
  logical                                    :: useGPU
  logical                                    :: useCCL
  integer(kind=c_int)                        :: numGPU
  integer(kind=c_intptr_t)                   :: num
  integer(kind=c_intptr_t)                   :: tmp1_dev, tmatc_dev, tmatr_dev, a_dev, tmp2_dev
  integer(kind=c_intptr_t)                   :: info_dev, info_abs_accumulated_dev
  type(c_ptr)                                :: tmp1_mpi_dev
  real(kind=rck), pointer           :: tmp1_mpi_fortran_ptr(:,:)
  integer(kind=c_intptr_t)                   :: a_off, tmatc_off, tmatr_off
  type(c_ptr)                                :: tmatc_mpi_dev
  real(kind=rck), pointer           :: tmatc_mpi_fortran_ptr(:,:)

  integer(kind=c_intptr_t), parameter        :: size_of_datatype = size_of_&
                                                            &single&
                                                            &_&
                                                            &real

  integer(kind=c_intptr_t)                   :: gpublasHandle, gpusolverHandle, my_stream, offset
  integer(kind=c_int)                        :: gpu_cholesky
  integer(kind=ik)                           :: blocking

  success = .true.
  useGPU = .false.
  useCCL = .false.



  ! check whether the above setting should be overriden
  if (obj%is_set("gpu_cholesky") == 1) then
    call obj%get("gpu_cholesky", gpu_cholesky, error)
    if (error .ne. ELPA_OK) then
      print *,"Problem getting option for gpu_cholesky. Aborting..."
      stop 1
    endif
    if (useGPU .and. gpu_cholesky .eq. 0) then
      useGPU = .false.
    else if (.not.(useGPU) .and. gpu_cholesky .eq. 1) then
      useGPU = .true.
    else 
    endif
  else 
    ! no override by user 
  endif

  if (.not.(useGPU)) then
  endif

  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  
  call obj%timer%start("elpa_cholesky_&
  &real&
  &_&
  &single&
  &"//gpuString)

  nrThreads=1

  na         = obj%na
  matrixRows = obj%local_nrows
  nblk       = obj%nblk
  matrixCols = obj%local_ncols

  call obj%get("debug",debug,error)
  if (error .ne. ELPA_OK) then
    write(error_unit,*) "ELPA_CHOLESKY: Problem getting option for debug settings. Aborting..."
    success = .false.
    return
  endif
  if (debug == 1) then
    wantDebug = .true.
  else
    wantDebug = .false.
  endif

  mpi_comm_all    = obj%mpi_setup%mpi_comm_parent
  mpi_comm_cols   = obj%mpi_setup%mpi_comm_cols
  mpi_comm_rows   = obj%mpi_setup%mpi_comm_rows
    
  myid    = obj%mpi_setup%myRank_comm_parent
  my_prow = obj%mpi_setup%myRank_comm_rows
  my_pcol = obj%mpi_setup%myRank_comm_cols
        
  np_rows = obj%mpi_setup%nRanks_comm_rows
  np_cols = obj%mpi_setup%nRanks_comm_cols

  success = .true.

  ! Matrix is split into tiles; work is done only for tiles on the diagonal or above

  call obj%timer%start("prepare")

  if (obj%is_set("blocking_in_cholesky") == 1) then
    call obj%get("blocking_in_cholesky", blocking, error)
    if (error .ne. ELPA_OK) then
      write(error_unit,*) "ELPA_CHOLESKY: Problem in getting keyword 'blocking_in_cholesky'. Aborting..."
      stop 1
    endif
  else
    if (useGPU) then
      blocking = 1024
    else
      blocking = 128
    endif
    call obj%set("blocking_in_cholesky", blocking, error)
    if (error .ne. ELPA_OK) then
      write(error_unit,*) "ELPA_CHOLESKY: Problem in setting keyword 'blocking_in_cholesky'. Aborting..."
      stop 1
    endif
  endif 


  tile_size = nblk*least_common_multiple(np_rows,np_cols) ! minimum global tile size
  tile_size = ((blocking*max(np_rows,np_cols)-1)/tile_size+1)*tile_size ! make local tiles at least 128 wide

  l_rows_tile = tile_size/np_rows ! local rows of a tile
  l_cols_tile = tile_size/np_cols ! local cols of a tile

  l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a
  l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a

  if (useGPU) then
    call obj%timer%start("check_for_gpu")
    if (check_for_gpu(obj, myid, numGPU)) then
      ! set the neccessary parameters
      call set_gpu_parameters()
    else
      print *,"ELPA_CHOLESKY: GPUs are requested but not detected! Aborting..."
      success = .false.
      return
    endif
    call obj%timer%stop("check_for_gpu")
  else ! useGPU
  endif ! useGPU


  if (useGPU) then

    successGPU = gpu_malloc(tmp1_dev, nblk*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmp1_dev", 335,  successGPU)

    successGPU = gpu_memset(tmp1_dev, 0, nblk*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmp1_dev", 346,  successGPU)

    successGPU = gpu_malloc(tmp2_dev, nblk*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmp2_dev", 350,  successGPU)

    successGPU = gpu_memset(tmp2_dev, 0, nblk*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmp2_dev", 361,  successGPU)

    successGPU = gpu_malloc(tmatc_dev, l_cols*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmatc_dev", 365,  successGPU)

    successGPU = gpu_memset(tmatc_dev, 0, l_cols*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmatc_dev", 376,  successGPU)

    successGPU = gpu_malloc(tmatr_dev, l_rows*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmatr_dev", 380,  successGPU)

    successGPU = gpu_memset(tmatr_dev, 0, l_rows*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmatr_dev", 391,  successGPU)

    successGPU = gpu_malloc(a_dev, matrixRows*matrixCols*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: a_dev", 396,  successGPU)
  endif ! useGPU


  allocate(tmp1(nblk*nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmp1", 441,  istat,  errorMessage)

  allocate(tmp2(nblk,nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmp2", 444,  istat,  errorMessage)



  tmp1 = 0
  tmp2 = 0

  allocate(tmatr(l_rows,nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmatr", 460,  istat,  errorMessage)

  allocate(tmatc(l_cols,nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmatc", 463,  istat,  errorMessage)


  tmatr = 0
  tmatc = 0

  if (useGPU) then
    num = matrixRows*matrixCols* size_of_datatype
    successGPU = gpu_memcpy(a_dev, int(loc(a(1,1)),kind=c_intptr_t), &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
    call check_memcpy_GPU_f("elpa_cholesky 1: memcpy a-> a_dev", 495,  successGPU)
  endif


  call obj%timer%stop("prepare")
  call obj%timer%start("loop1")
  do n = 1, na, nblk


    ! Calculate first local row and column of the still remaining matrix
    ! on the local processor

    l_row1 = local_index(n, my_prow, np_rows, nblk, +1)
    l_col1 = local_index(n, my_pcol, np_cols, nblk, +1)

    l_rowx = local_index(n+nblk, my_prow, np_rows, nblk, +1)
    l_colx = local_index(n+nblk, my_pcol, np_cols, nblk, +1)

    if (n+nblk > na) then

      ! This is the last step, just do a Cholesky-Factorization
      ! of the remaining block

      if (useGPU) then
        if (my_prow==prow(n, nblk, np_rows) .and. my_pcol==pcol(n, nblk, np_cols)) then

          call obj%timer%start("lapack")
          num = matrixRows*matrixCols* size_of_datatype
          successGPU = gpu_memcpy(int(loc(a(1,1)),kind=c_intptr_t), a_dev,  &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a", 582,  successGPU)

          call SPOTRF('U', int(na-n+1,kind=BLAS_KIND), a(l_row1,l_col1), &
                             int(matrixRows,kind=BLAS_KIND), infoBLAS )
          info = int(infoBLAS,kind=ik)

          num = matrixRows*matrixCols* size_of_datatype
          successGPU = gpu_memcpy(a_dev, int(loc(a(1,1)),kind=c_intptr_t), &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a", 600,  successGPU)

          call obj%timer%stop("lapack")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &real&

            &: Error in dpotrf: ",info
            success = .false. ! "
            return
          endif ! info
        endif ! (my_prow==prow(n, nblk, np_rows) .and. my_pcol==pcol(n, nblk, np_cols))

      else ! useGPU
        if (my_prow==prow(n, nblk, np_rows) .and. my_pcol==pcol(n, nblk, np_cols)) then
          call obj%timer%start("blas")

          call SPOTRF('U', int(na-n+1,kind=BLAS_KIND), a(l_row1,l_col1), &
                             int(matrixRows,kind=BLAS_KIND), infoBLAS )
          info = int(infoBLAS,kind=ik)
          call obj%timer%stop("blas")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &real&

            &: Error in dpotrf: ",info
            success = .false. ! "
            return
          endif

        endif
      endif ! useGPU
      
      exit ! Loop
    endif ! (n+nblk > na) 

    ! This is not the last step

    if (my_prow==prow(n, nblk, np_rows)) then

      if (my_pcol==pcol(n, nblk, np_cols)) then

        if (useGPU) then
          call obj%timer%start("lapack")

          num = matrixRows*matrixCols* size_of_datatype

          successGPU = gpu_memcpy(int(loc(a(1,1)),kind=c_intptr_t), a_dev,  &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a", 742,  successGPU)


          call SPOTRF('U', int(nblk,kind=BLAS_KIND), a(l_row1,l_col1), &
                               int(matrixRows,kind=BLAS_KIND) , infoBLAS )
          info = int(infoBLAS,kind=ik)

          num = matrixRows*matrixCols* size_of_datatype
          successGPU = gpu_memcpy(a_dev, int(loc(a(1,1)),kind=c_intptr_t), &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a", 769,  successGPU)

          call obj%timer%stop("lapack")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &real&

            &: Error in dpotrf: ",info
            success = .false. ! "
            return
          endif ! info
        else ! useGPU
          ! The process owning the upper left remaining block does the
          ! Cholesky-Factorization of this block
          call obj%timer%start("blas")

          call SPOTRF('U', int(nblk,kind=BLAS_KIND), a(l_row1,l_col1), &
                               int(matrixRows,kind=BLAS_KIND) , infoBLAS )
          info = int(infoBLAS,kind=ik)
          call obj%timer%stop("blas")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &real&

            &: Error in dpotrf 2: ",info
            success = .false. ! "
            return
          endif ! info/=0
        endif ! useGPU
 
        if (useGPU) then
          my_stream = obj%gpu_setup%my_stream
          call gpu_copy_float_a_tmp1 (a_dev, tmp1_dev, l_row1, l_col1, matrixRows, nblk, my_stream)
        else ! useGPU
          nc = 0
          do i=1,nblk
            tmp1(nc+1:nc+i) = a(l_row1:l_row1+i-1,l_col1+i-1)
            nc = nc+i
          enddo
        endif ! useGPU
      endif ! (my_pcol==pcol(n, nblk, np_cols))


      if (useGPU) then
        my_stream = obj%gpu_setup%my_stream
        call gpu_copy_float_tmp1_tmp2 (tmp1_dev, tmp2_dev, nblk, nblk, my_stream)
      else ! useGPU
        nc = 0
        do i=1,nblk
          tmp2(1:i,i) = tmp1(nc+1:nc+i)
          nc = nc+i
        enddo
      endif ! useGPU

      if (useGPU) then
        call obj%timer%start("gpublas")
        gpublasHandle = obj%gpu_setup%gpublasHandleArray(0)
        if (l_cols-l_colx+1 > 0) then
          a_off = (l_row1-1 + (l_colx-1)*matrixRows) * size_of_datatype
          call gpublas_STRSM('L', 'U', 'T', 'N', nblk, l_cols-l_colx+1, ONE, &
                            tmp2_dev, nblk, a_dev+a_off, matrixRows, gpublasHandle)
          if (wantDebug .and. gpu_vendor() /= SYCL_GPU) successGPU = gpu_DeviceSynchronize()
        endif
        call obj%timer%stop("gpublas")

      else ! useGPU

        call obj%timer%start("lapack")
        if (l_cols-l_colx+1>0) then
          call STRSM('L', 'U', 'T', 'N', int(nblk,kind=BLAS_KIND),  &
                            int(l_cols-l_colx+1,kind=BLAS_KIND), ONE, tmp2, &
                            int(ubound(tmp2,dim=1),kind=BLAS_KIND), a(l_row1,l_colx), int(matrixRows,kind=BLAS_KIND) )
        endif
        call obj%timer%stop("lapack")
      endif ! useGPU
    endif ! (my_prow==prow(n, nblk, np_rows))


    if (useGPU) then
      if (my_prow==prow(n, nblk, np_rows)) then
        ! if l_cols-l_colx+1 == 0 kernel launch with 0 blocks => raises error
        call obj%timer%start("a_tmatc_kernel")
        if (l_cols-l_colx+1>0) then
           my_stream = obj%gpu_setup%my_stream
           call gpu_copy_float_a_tmatc(a_dev, tmatc_dev, nblk, matrixRows, l_cols, l_colx, l_row1, my_stream)
        endif
        call obj%timer%stop("a_tmatc_kernel")
      endif
    else ! useGPU
      do i=1,nblk
        if (my_prow==prow(n, nblk, np_rows)) tmatc(l_colx:l_cols,i) = a(l_row1+i-1,l_colx:l_cols)
      enddo
    endif ! useGPU

    if (useGPU  .and. .not. useCCL) then
      num = l_cols*nblk*size_of_datatype
      successGPU = gpu_memcpy(int(loc(tmatc),kind=c_intptr_t), tmatc_dev, num, &
                              gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("elpa_cholesky: tmatc_dev to tmatc", 1201,  successGPU)
    endif


    if (useGPU .and. .not. useCCL) then

      ! not needed if transpose vec CCL
      num = l_rows*nblk*size_of_datatype
      successGPU = gpu_memcpy(int(loc(tmatr),kind=c_intptr_t), tmatr_dev, num, &
                              gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("elpa_cholesky: tmatr_dev to tmatr", 1220,  successGPU)
    endif ! (useGPU .and. .not. useCCL)


    if (useCCL) then

    else ! useCCL
      call elpa_transpose_vectors_&
      &real&
      &_&
      &single &
      (obj, tmatc, ubound(tmatc,dim=1), mpi_comm_cols, &
      tmatr, ubound(tmatr,dim=1), mpi_comm_rows, &
      n, na, nblk, nblk, nrThreads, .false., success)
    endif ! useCCL

    if (useGPU .and. .not. useCCL) then

      num = l_rows*nblk*size_of_datatype
      successGPU = gpu_memcpy(tmatr_dev, int(loc(tmatr),kind=c_intptr_t), num, &
                              gpuMemcpyHostToDevice)
      call check_memcpy_GPU_f("elpa_cholesky: tmat to tmatr_dev", 1271,  successGPU)
    endif

    if (useGPU) then
      do i=0,(na-1)/tile_size
        lcs = max(l_colx,i*l_cols_tile+1)
        lce = min(l_cols,(i+1)*l_cols_tile)
        lrs = l_rowx
        lre = min(l_rows,(i+1)*l_rows_tile)
        if (lce  < lcs .or. lre < lrs) cycle
        call obj%timer%start("gpublas")
        gpublasHandle = obj%gpu_setup%gpublasHandleArray(0)
        tmatr_off = (lrs-1 + (1-1)*l_rows) * size_of_datatype
        tmatc_off = (lcs-1 + (1-1)*l_cols) * size_of_datatype
        a_off = (lrs-1 + (lcs-1)*matrixRows) * size_of_datatype
        call gpublas_SGEMM('N', 'T', lre-lrs+1, lce-lcs+1, nblk, &
                            -ONE, tmatr_dev+tmatr_off, l_rows, tmatc_dev+tmatc_off, l_cols, ONE, &
                            a_dev+a_off, matrixRows, gpublasHandle)
        if (wantDebug .and. gpu_vendor() /= SYCL_GPU) successGPU = gpu_DeviceSynchronize()
        call obj%timer%stop("gpublas")
      enddo
    else !useGPU
      do i=0,(na-1)/tile_size
        lcs = max(l_colx,i*l_cols_tile+1)
        lce = min(l_cols,(i+1)*l_cols_tile)
        lrs = l_rowx
        lre = min(l_rows,(i+1)*l_rows_tile)
        if (lce<lcs .or. lre<lrs) cycle
        call obj%timer%start("blas")
        call SGEMM('N', 'T', int(lre-lrs+1,kind=BLAS_KIND), int(lce-lcs+1,kind=BLAS_KIND), &
                            int(nblk,kind=BLAS_KIND), -ONE,  &
                            tmatr(lrs,1), int(ubound(tmatr,dim=1),kind=BLAS_KIND), tmatc(lcs,1), &
                            int(ubound(tmatc,dim=1),kind=BLAS_KIND), &
                            ONE, a(lrs,lcs), int(matrixRows,kind=BLAS_KIND))
        call obj%timer%stop("blas")
      enddo
    endif ! useGPU


  enddo ! n = 1, na, nblk

  call obj%timer%stop("loop1")
  
  call obj%timer%start("deallocate")
  if (useGPU) then



    successGPU = gpu_free(tmp1_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmp1_dev", 1351,  successGPU)

    successGPU = gpu_free(tmp2_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmp1_dev", 1354,  successGPU)

    successGPU = gpu_free(tmatc_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmatc_dev", 1357,  successGPU)

    successGPU = gpu_free(tmatr_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmatr_dev", 1360,  successGPU)

  endif


  deallocate(tmp1, tmp2, tmatr, tmatc, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("elpa_cholesky: tmp1, tmp2, tmatr, tmatc", 1384,  istat,  errorMessage)

  call obj%timer%stop("deallocate")


  call obj%timer%start("loop2")

  call obj%timer%start("copy1")

  if (useGPU) then
    num = matrixRows*matrixCols* size_of_datatype
    successGPU = gpu_memcpy(int(loc(a(1,1)),kind=c_intptr_t), a_dev,  &
                     matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
    call check_memcpy_GPU_f("elpa_cholesky: memcpy 2 a-> a_dev", 1405,  successGPU)
  endif
  ! Set the lower triangle to 0, it contains garbage (form the above matrix multiplications)

  call obj%timer%stop("copy1")
  do i=1,na
    if (my_pcol==pcol(i, nblk, np_cols)) then
      ! column i is on local processor
      l_col1 = local_index(i  , my_pcol, np_cols, nblk, +1) ! local column number
      l_row1 = local_index(i+1, my_prow, np_rows, nblk, +1) ! first row below diagonal
      if (useGPU) then
        call obj%timer%start("copy2")
!#ifndef DEVICE_POINTER
        ! set to zero maybe a kerne would be better here
          offset = (l_row1-1 + matrixRows * (l_col1-1)) * size_of_datatype
          num = (l_rows-l_row1+1)*size_of_datatype
          successGPU = gpu_memset(a_dev+offset, 0, num)
          call check_memcpy_GPU_f("elpa_cholesky: memset a_dev", 1450,  successGPU)
!#else 
!          !offset = (l_row1-1 + matrixRows * (l_col1-1)) * size_of_datatype
!          !num = (l_rows-l_row1+1)*size_of_datatype
!#ifdef WITH_GPU_STREAMS
!          !my_stream = obj%gpu_setup%my_stream
!          !successGPU = gpu_memset_async(a_dev+offset, 0, num, my_stream)
!          !call check_memcpy_GPU_f("elpa_cholesky: memset a_dev", 1458,  successGPU)
!
!          !successGPU = gpu_stream_synchronize(my_stream)
!          !call check_memcpy_GPU_f("elpa_cholesky: memset", 1461,  successGPU)
!#else
!          !successGPU = gpu_memset(a_dev+offset, 0, num)
!          !call check_memcpy_GPU_f("elpa_cholesky: memset tmp1_dev", 1464,  successGPU)
!#endif
!
!
!!#if defined(WITH_NVIDIA_CUSOLVER) || defined(WITH_AMD_ROCSOLVER)
!!        successGPU = gpu_memcpy(int(loc(a_tmp(1,1)),kind=c_intptr_t), a_dev, &
!!                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
!!        call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a_tmp", 1471,  successGPU)
!!#endif
!!        a_tmp(l_row1:l_rows,l_col1) = 0
!!#if defined(WITH_NVIDIA_CUSOLVER) || defined(WITH_AMD_ROCSOLVER)
!!        successGPU = gpu_memcpy(a_dev, int(loc(a_tmp(1,1)),kind=c_intptr_t), &
!!                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
!!        call check_memcpy_GPU_f("elpa_cholesky: memcpy a_tmp-> a_dev", 1477,  successGPU)
!#endif

!#endif 

  call obj%timer%stop("copy2")

      else ! useGPU
        a(l_row1:l_rows,l_col1) = 0
      endif ! useGPU
    endif ! (my_pcol==pcol(i, nblk, np_cols))
  enddo
  !endif ! useGPU

  call obj%timer%stop("loop2")


  call obj%timer%start("cleanup")
  if (useGPU) then
    ! copy back
    successGPU = gpu_memcpy(int(loc(a(1,1)),kind=c_intptr_t), a_dev,  &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
    call check_memcpy_GPU_f("elpa_cholesky: memcpy a-> d_dev", 1502,  successGPU)

    successGPU = gpu_free(a_dev)
    call check_dealloc_GPU_f("elpa_cholesky: a_dev", 1505,  successGPU)


  endif

  ! restore original OpenMP settings
  call obj%timer%stop("cleanup")
  call obj%timer%stop("elpa_cholesky_&
  &real&
  &_&
  &single&
  &"//gpuString)

    end function elpa_cholesky_a_h_a_real_single_impl






















!> \brief  elpa_cholesky_d_ptr_real_single_impl: Cholesky factorization of a double-precision real symmetric matrix
!> \details
!> \param  obj                    elpa_t object contains:
!> \param     - obj%na            Order of matrix
!> \param     - obj%local_nrows   Leading dimension of a
!> \param     - obj%local_ncols   local columns of matrix a
!> \param     - obj%nblk          blocksize of cyclic distribution, must be the same in both directions!
!> \param     - obj%mpi_comm_rows MPI communicator for rows
!> \param     - obj%mpi_comm_cols MPI communicator for columns
!> \param     - obj%wantDebug     logical, more debug information on failure
!> \param  aDev(lda,matrixCols)   Distributed matrix which should be inverted as type(c_ptr) living on a device
!>                                Distribution is like in Scalapack.
!>                                Only upper triangle needs to be set.
!>                                The lower triangle is not referenced.
!> \result succes                 logical, reports success or failure
   function elpa_cholesky_d_ptr_real_single_impl(obj, aDev) result(success)
!    This file is part of ELPA.
!
!    The ELPA library was originally created by the ELPA consortium,
!    consisting of the following organizations:
!
!    - Max Planck Computing and Data Facility (MPCDF), formerly known as
!      Rechenzentrum Garching der Max-Planck-Gesellschaft (RZG),
!    - Bergische Universität Wuppertal, Lehrstuhl für angewandte
!      Informatik,
!    - Technische Universität München, Lehrstuhl für Informatik mit
!      Schwerpunkt Wissenschaftliches Rechnen ,
!    - Fritz-Haber-Institut, Berlin, Abt. Theorie,
!    - Max-Plack-Institut für Mathematik in den Naturwissenschaften,
!      Leipzig, Abt. Komplexe Strukutren in Biologie und Kognition,
!      and
!    - IBM Deutschland GmbH
!
!    This particular source code file contains additions, changes and
!    enhancements authored by Intel Corporation which is not part of
!    the ELPA consortium.
!
!    More information can be found here:
!    http://elpa.mpcdf.mpg.de/
!
!    ELPA is free software: you can redistribute it and/or modify
!    it under the terms of the version 3 of the license of the
!    GNU Lesser General Public License as published by the Free
!    Software Foundation.
!
!    ELPA is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with ELPA.  If not, see <http://www.gnu.org/licenses/>
!
!    ELPA reflects a substantial effort on the part of the original
!    ELPA consortium, and we ask you to respect the spirit of the
!    license that we chose: i.e., please contribute any changes you
!    may have back to the original ELPA library distribution, and keep
!    any derivatives of ELPA under the same license that we chose for
!    the original distribution, the GNU Lesser General Public License.




!cannot use "../src/cholesky/./../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)



  use elpa1_compute
  use elpa_utilities
  use elpa_mpi
  use precision
  use elpa_abstract_impl
  use elpa_omp
  use elpa_blas_interfaces
  use elpa_gpu
  use mod_check_for_gpu
  use invert_trm_gpu !, only : gpu_copy_&
                     !        &single&
                     !        &_tmp1_tmp2, &
                     !        gpu_copy_&
                     !        &single&
                     !        &_a_tmp1
  use cholesky_gpu
  use mod_query_gpu_usage
  implicit none
!    Copyright 2011, A. Marek
!
!    This file is part of ELPA.
!
!    The ELPA library was originally created by the ELPA consortium,
!    consisting of the following organizations:
!
!    - Max Planck Computing and Data Facility (MPCDF), formerly known as
!      Rechenzentrum Garching der Max-Planck-Gesellschaft (RZG),
!    - Bergische Universität Wuppertal, Lehrstuhl für angewandte
!      Informatik,
!    - Technische Universität München, Lehrstuhl für Informatik mit
!      Schwerpunkt Wissenschaftliches Rechnen ,
!    - Fritz-Haber-Institut, Berlin, Abt. Theorie,
!    - Max-Plack-Institut für Mathematik in den Naturwissenschaften,
!      Leipzig, Abt. Komplexe Strukutren in Biologie und Kognition,
!      and
!    - IBM Deutschland GmbH
!
!    This particular source code file contains additions, changes and
!    enhancements authored by Intel Corporation which is not part of
!    the ELPA consortium.
!
!    More information can be found here:
!    http://elpa.mpcdf.mpg.de/
!
!    ELPA is free software: you can redistribute it and/or modify
!    it under the terms of the version 3 of the license of the
!    GNU Lesser General Public License as published by the Free
!    Software Foundation.
!
!    ELPA is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with ELPA.  If not, see <http://www.gnu.org/licenses/>
!
!    ELPA reflects a substantial effort on the part of the original
!    ELPA consortium, and we ask you to respect the spirit of the
!    license that we chose: i.e., please contribute any changes you
!    may have back to the original ELPA library distribution, and keep
!    any derivatives of ELPA under the same license that we chose for
!    the original distribution, the GNU Lesser General Public License.
!
!
  integer, parameter :: rk = C_FLOAT
  integer, parameter :: rck = C_FLOAT
  real(kind=rck), parameter      :: ZERO=0.0_rk, ONE = 1.0_rk

  class(elpa_abstract_impl_t), intent(inout) :: obj
  integer(kind=ik)                           :: na, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols, mpi_comm_all
  type(c_ptr)                                :: aDev
!#if !defined(WITH_NVIDIA_CUSOLVER) && !defined(WITH_AMD_ROCSOLVER)
  real(kind=rck), allocatable       :: a_tmp(:,:)
!#endif
  integer(kind=ik)                           :: my_prow, my_pcol, np_rows, np_cols, myid
  integer(kind=MPI_KIND)                     :: mpierr, my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI, myidMPI
  integer(kind=ik)                           :: l_cols, l_rows, l_col1, l_row1, l_colx, l_rowx
  integer(kind=ik)                           :: n, nc, i, info
  integer(kind=BLAS_KIND)                    :: infoBLAS
  integer(kind=ik)                           :: lcs, lce, lrs, lre
  integer(kind=ik)                           :: tile_size, l_rows_tile, l_cols_tile

  real(kind=rck), allocatable       :: tmp1(:), tmp2(:,:), tmatr(:,:), tmatc(:,:)
  logical                                    :: wantDebug
  logical                                    :: success
  integer(kind=ik)                           :: istat, debug, error
  character(200)                             :: errorMessage
  integer(kind=ik)                           :: nrThreads, limitThreads
  character(20)                              :: gpuString
  logical                                    :: successGPU
  logical                                    :: useGPU
  logical                                    :: useCCL
  integer(kind=c_int)                        :: numGPU
  integer(kind=c_intptr_t)                   :: num
  integer(kind=c_intptr_t)                   :: tmp1_dev, tmatc_dev, tmatr_dev, a_dev, tmp2_dev
  integer(kind=c_intptr_t)                   :: info_dev, info_abs_accumulated_dev
  type(c_ptr)                                :: tmp1_mpi_dev
  real(kind=rck), pointer           :: tmp1_mpi_fortran_ptr(:,:)
  integer(kind=c_intptr_t)                   :: a_off, tmatc_off, tmatr_off
  type(c_ptr)                                :: tmatc_mpi_dev
  real(kind=rck), pointer           :: tmatc_mpi_fortran_ptr(:,:)

  integer(kind=c_intptr_t), parameter        :: size_of_datatype = size_of_&
                                                            &single&
                                                            &_&
                                                            &real

  integer(kind=c_intptr_t)                   :: gpublasHandle, gpusolverHandle, my_stream, offset
  integer(kind=c_int)                        :: gpu_cholesky
  integer(kind=ik)                           :: blocking

  success = .true.
  useGPU = .false.
  useCCL = .false.

  useGPU = .true.

  if (.not.(useGPU)) then
    print *,"You used the interface for device pointers for elpa_cholesky but did not specify GPU usage!. Aborting..."
    stop 1
  endif

  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  
  call obj%timer%start("elpa_cholesky_&
  &real&
  &_&
  &single&
  &"//gpuString)

  nrThreads=1

  na         = obj%na
  matrixRows = obj%local_nrows
  nblk       = obj%nblk
  matrixCols = obj%local_ncols

  call obj%get("debug",debug,error)
  if (error .ne. ELPA_OK) then
    write(error_unit,*) "ELPA_CHOLESKY: Problem getting option for debug settings. Aborting..."
    success = .false.
    return
  endif
  if (debug == 1) then
    wantDebug = .true.
  else
    wantDebug = .false.
  endif

  mpi_comm_all    = obj%mpi_setup%mpi_comm_parent
  mpi_comm_cols   = obj%mpi_setup%mpi_comm_cols
  mpi_comm_rows   = obj%mpi_setup%mpi_comm_rows
    
  myid    = obj%mpi_setup%myRank_comm_parent
  my_prow = obj%mpi_setup%myRank_comm_rows
  my_pcol = obj%mpi_setup%myRank_comm_cols
        
  np_rows = obj%mpi_setup%nRanks_comm_rows
  np_cols = obj%mpi_setup%nRanks_comm_cols

  success = .true.

  ! Matrix is split into tiles; work is done only for tiles on the diagonal or above

  call obj%timer%start("prepare")

  if (obj%is_set("blocking_in_cholesky") == 1) then
    call obj%get("blocking_in_cholesky", blocking, error)
    if (error .ne. ELPA_OK) then
      write(error_unit,*) "ELPA_CHOLESKY: Problem in getting keyword 'blocking_in_cholesky'. Aborting..."
      stop 1
    endif
  else
    if (useGPU) then
      blocking = 1024
    else
      blocking = 128
    endif
    call obj%set("blocking_in_cholesky", blocking, error)
    if (error .ne. ELPA_OK) then
      write(error_unit,*) "ELPA_CHOLESKY: Problem in setting keyword 'blocking_in_cholesky'. Aborting..."
      stop 1
    endif
  endif 


  tile_size = nblk*least_common_multiple(np_rows,np_cols) ! minimum global tile size
  tile_size = ((blocking*max(np_rows,np_cols)-1)/tile_size+1)*tile_size ! make local tiles at least 128 wide

  l_rows_tile = tile_size/np_rows ! local rows of a tile
  l_cols_tile = tile_size/np_cols ! local cols of a tile

  l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a
  l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a

  if (useGPU) then
    call obj%timer%start("check_for_gpu")
    if (check_for_gpu(obj, myid, numGPU)) then
      ! set the neccessary parameters
      call set_gpu_parameters()
    else
      print *,"ELPA_CHOLESKY: GPUs are requested but not detected! Aborting..."
      success = .false.
      return
    endif
    call obj%timer%stop("check_for_gpu")
  else ! useGPU
  endif ! useGPU


  if (useGPU) then

    successGPU = gpu_malloc(tmp1_dev, nblk*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmp1_dev", 335,  successGPU)

    successGPU = gpu_memset(tmp1_dev, 0, nblk*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmp1_dev", 346,  successGPU)

    successGPU = gpu_malloc(tmp2_dev, nblk*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmp2_dev", 350,  successGPU)

    successGPU = gpu_memset(tmp2_dev, 0, nblk*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmp2_dev", 361,  successGPU)

    successGPU = gpu_malloc(tmatc_dev, l_cols*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmatc_dev", 365,  successGPU)

    successGPU = gpu_memset(tmatc_dev, 0, l_cols*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmatc_dev", 376,  successGPU)

    successGPU = gpu_malloc(tmatr_dev, l_rows*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmatr_dev", 380,  successGPU)

    successGPU = gpu_memset(tmatr_dev, 0, l_rows*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmatr_dev", 391,  successGPU)

    a_dev = transfer(aDev, a_dev)
!#if !defined(WITH_NVIDIA_CUSOLVER) && !defined(WITH_AMD_ROCSOLVER)
    allocate(a_tmp(obj%local_nrows,obj%local_ncols), stat=istat, errmsg=errorMessage)
    call check_allocate_f("elpa_cholesky: a_tmp", 407,  istat,  errorMessage)
!#endif 
  endif ! useGPU


  allocate(tmp1(nblk*nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmp1", 441,  istat,  errorMessage)

  allocate(tmp2(nblk,nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmp2", 444,  istat,  errorMessage)



  tmp1 = 0
  tmp2 = 0

  allocate(tmatr(l_rows,nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmatr", 460,  istat,  errorMessage)

  allocate(tmatc(l_cols,nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmatc", 463,  istat,  errorMessage)


  tmatr = 0
  tmatc = 0



  call obj%timer%stop("prepare")
  call obj%timer%start("loop1")
  do n = 1, na, nblk


    ! Calculate first local row and column of the still remaining matrix
    ! on the local processor

    l_row1 = local_index(n, my_prow, np_rows, nblk, +1)
    l_col1 = local_index(n, my_pcol, np_cols, nblk, +1)

    l_rowx = local_index(n+nblk, my_prow, np_rows, nblk, +1)
    l_colx = local_index(n+nblk, my_pcol, np_cols, nblk, +1)

    if (n+nblk > na) then

      ! This is the last step, just do a Cholesky-Factorization
      ! of the remaining block

      if (useGPU) then
        if (my_prow==prow(n, nblk, np_rows) .and. my_pcol==pcol(n, nblk, np_cols)) then

          call obj%timer%start("lapack")
          num = matrixRows*matrixCols* size_of_datatype

          successGPU = gpu_memcpy(int(loc(a_tmp(1,1)),kind=c_intptr_t), a_dev,  &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a_tmp", 615,  successGPU)

          call SPOTRF('U', int(na-n+1,kind=BLAS_KIND), a_tmp(l_row1,l_col1), &
                             int(matrixRows,kind=BLAS_KIND), infoBLAS )
          info = int(infoBLAS,kind=ik)

          num = matrixRows*matrixCols* size_of_datatype
          successGPU = gpu_memcpy(a_dev, int(loc(a_tmp(1,1)),kind=c_intptr_t), &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a_tmp", 633,  successGPU)

          call obj%timer%stop("lapack")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &real&

            &: Error in dpotrf: ",info
            success = .false. ! "
            return
          endif ! info
        endif ! (my_prow==prow(n, nblk, np_rows) .and. my_pcol==pcol(n, nblk, np_cols))

      else ! useGPU
        if (my_prow==prow(n, nblk, np_rows) .and. my_pcol==pcol(n, nblk, np_cols)) then
          call obj%timer%start("blas")

          info = int(infoBLAS,kind=ik)
          call obj%timer%stop("blas")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &real&

            &: Error in dpotrf: ",info
            success = .false. ! "
            return
          endif

        endif
      endif ! useGPU
      
      exit ! Loop
    endif ! (n+nblk > na) 

    ! This is not the last step

    if (my_prow==prow(n, nblk, np_rows)) then

      if (my_pcol==pcol(n, nblk, np_cols)) then

        if (useGPU) then
          call obj%timer%start("lapack")

          num = matrixRows*matrixCols* size_of_datatype

          successGPU = gpu_memcpy(int(loc(a_tmp(1,1)),kind=c_intptr_t), a_dev,  &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a_tmp", 783,  successGPU)

          call SPOTRF('U', int(nblk,kind=BLAS_KIND), a_tmp(l_row1,l_col1), &
                               int(matrixRows,kind=BLAS_KIND) , infoBLAS )
          info = int(infoBLAS,kind=ik)
          num = matrixRows*matrixCols* size_of_datatype
          successGPU = gpu_memcpy(a_dev, int(loc(a_tmp(1,1)),kind=c_intptr_t), &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a_tmp", 799,  successGPU)

          call obj%timer%stop("lapack")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &real&

            &: Error in dpotrf: ",info
            success = .false. ! "
            return
          endif ! info
        else ! useGPU
          ! The process owning the upper left remaining block does the
          ! Cholesky-Factorization of this block
          call obj%timer%start("blas")

          info = int(infoBLAS,kind=ik)
          call obj%timer%stop("blas")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &real&

            &: Error in dpotrf 2: ",info
            success = .false. ! "
            return
          endif ! info/=0
        endif ! useGPU
 
        if (useGPU) then
          my_stream = obj%gpu_setup%my_stream
          call gpu_copy_float_a_tmp1 (a_dev, tmp1_dev, l_row1, l_col1, matrixRows, nblk, my_stream)
        else ! useGPU
          nc = 0
          do i=1,nblk
            nc = nc+i
          enddo
        endif ! useGPU
      endif ! (my_pcol==pcol(n, nblk, np_cols))


      if (useGPU) then
        my_stream = obj%gpu_setup%my_stream
        call gpu_copy_float_tmp1_tmp2 (tmp1_dev, tmp2_dev, nblk, nblk, my_stream)
      else ! useGPU
        nc = 0
        do i=1,nblk
          tmp2(1:i,i) = tmp1(nc+1:nc+i)
          nc = nc+i
        enddo
      endif ! useGPU

      if (useGPU) then
        call obj%timer%start("gpublas")
        gpublasHandle = obj%gpu_setup%gpublasHandleArray(0)
        if (l_cols-l_colx+1 > 0) then
          a_off = (l_row1-1 + (l_colx-1)*matrixRows) * size_of_datatype
          call gpublas_STRSM('L', 'U', 'T', 'N', nblk, l_cols-l_colx+1, ONE, &
                            tmp2_dev, nblk, a_dev+a_off, matrixRows, gpublasHandle)
          if (wantDebug .and. gpu_vendor() /= SYCL_GPU) successGPU = gpu_DeviceSynchronize()
        endif
        call obj%timer%stop("gpublas")

      else ! useGPU

        call obj%timer%start("lapack")
        call obj%timer%stop("lapack")
      endif ! useGPU
    endif ! (my_prow==prow(n, nblk, np_rows))


    if (useGPU) then
      if (my_prow==prow(n, nblk, np_rows)) then
        ! if l_cols-l_colx+1 == 0 kernel launch with 0 blocks => raises error
        call obj%timer%start("a_tmatc_kernel")
        if (l_cols-l_colx+1>0) then
           my_stream = obj%gpu_setup%my_stream
           call gpu_copy_float_a_tmatc(a_dev, tmatc_dev, nblk, matrixRows, l_cols, l_colx, l_row1, my_stream)
        endif
        call obj%timer%stop("a_tmatc_kernel")
      endif
    else ! useGPU
      do i=1,nblk
      enddo
    endif ! useGPU

    if (useGPU  .and. .not. useCCL) then
      num = l_cols*nblk*size_of_datatype
      successGPU = gpu_memcpy(int(loc(tmatc),kind=c_intptr_t), tmatc_dev, num, &
                              gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("elpa_cholesky: tmatc_dev to tmatc", 1201,  successGPU)
    endif


    if (useGPU .and. .not. useCCL) then

      ! not needed if transpose vec CCL
      num = l_rows*nblk*size_of_datatype
      successGPU = gpu_memcpy(int(loc(tmatr),kind=c_intptr_t), tmatr_dev, num, &
                              gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("elpa_cholesky: tmatr_dev to tmatr", 1220,  successGPU)
    endif ! (useGPU .and. .not. useCCL)


    if (useCCL) then

    else ! useCCL
      call elpa_transpose_vectors_&
      &real&
      &_&
      &single &
      (obj, tmatc, ubound(tmatc,dim=1), mpi_comm_cols, &
      tmatr, ubound(tmatr,dim=1), mpi_comm_rows, &
      n, na, nblk, nblk, nrThreads, .false., success)
    endif ! useCCL

    if (useGPU .and. .not. useCCL) then

      num = l_rows*nblk*size_of_datatype
      successGPU = gpu_memcpy(tmatr_dev, int(loc(tmatr),kind=c_intptr_t), num, &
                              gpuMemcpyHostToDevice)
      call check_memcpy_GPU_f("elpa_cholesky: tmat to tmatr_dev", 1271,  successGPU)
    endif

    if (useGPU) then
      do i=0,(na-1)/tile_size
        lcs = max(l_colx,i*l_cols_tile+1)
        lce = min(l_cols,(i+1)*l_cols_tile)
        lrs = l_rowx
        lre = min(l_rows,(i+1)*l_rows_tile)
        if (lce  < lcs .or. lre < lrs) cycle
        call obj%timer%start("gpublas")
        gpublasHandle = obj%gpu_setup%gpublasHandleArray(0)
        tmatr_off = (lrs-1 + (1-1)*l_rows) * size_of_datatype
        tmatc_off = (lcs-1 + (1-1)*l_cols) * size_of_datatype
        a_off = (lrs-1 + (lcs-1)*matrixRows) * size_of_datatype
        call gpublas_SGEMM('N', 'T', lre-lrs+1, lce-lcs+1, nblk, &
                            -ONE, tmatr_dev+tmatr_off, l_rows, tmatc_dev+tmatc_off, l_cols, ONE, &
                            a_dev+a_off, matrixRows, gpublasHandle)
        if (wantDebug .and. gpu_vendor() /= SYCL_GPU) successGPU = gpu_DeviceSynchronize()
        call obj%timer%stop("gpublas")
      enddo
    else !useGPU
      do i=0,(na-1)/tile_size
        lcs = max(l_colx,i*l_cols_tile+1)
        lce = min(l_cols,(i+1)*l_cols_tile)
        lrs = l_rowx
        lre = min(l_rows,(i+1)*l_rows_tile)
        if (lce<lcs .or. lre<lrs) cycle
        call obj%timer%start("blas")
        call obj%timer%stop("blas")
      enddo
    endif ! useGPU


  enddo ! n = 1, na, nblk

  call obj%timer%stop("loop1")
  
  call obj%timer%start("deallocate")
  if (useGPU) then



    successGPU = gpu_free(tmp1_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmp1_dev", 1351,  successGPU)

    successGPU = gpu_free(tmp2_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmp1_dev", 1354,  successGPU)

    successGPU = gpu_free(tmatc_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmatc_dev", 1357,  successGPU)

    successGPU = gpu_free(tmatr_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmatr_dev", 1360,  successGPU)

  endif


  deallocate(tmp1, tmp2, tmatr, tmatc, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("elpa_cholesky: tmp1, tmp2, tmatr, tmatc", 1384,  istat,  errorMessage)

  call obj%timer%stop("deallocate")


  call obj%timer%start("loop2")

  call obj%timer%start("copy1")

!  if (useGPU) then
!          ! is this needed: guess not
!#if !defined(WITH_NVIDIA_CUSOLVER) && !defined(WITH_AMD_ROCSOLVER)
!    num = matrixRows*matrixCols* size_of_datatype
!#ifdef WITH_GPU_STREAMS 
!    my_stream = obj%gpu_setup%my_stream
!    call gpu_memcpy_async_and_stream_synchronize &
!    ("elpa_choleksky: a_dev -> a_tmp", a_dev, 0_c_intptr_t, &
!                                    a_tmp(1:obj%local_nrows,1:obj%local_ncols), &
!                                    1, 1, num, gpuMemcpyDeviceToHost, my_stream, .false., .true., .false.)
!#else
!    successGPU = gpu_memcpy(int(loc(a_tmp(1,1)),kind=c_intptr_t), a_dev,  &
!                     matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
!    call check_memcpy_GPU_f("elpa_cholesky: memcpy 2 a-> a_dev", 1422,  successGPU)
!#endif
!#endif 
!  endif
  ! Set the lower triangle to 0, it contains garbage (form the above matrix multiplications)

  call obj%timer%stop("copy1")
  do i=1,na
    if (my_pcol==pcol(i, nblk, np_cols)) then
      ! column i is on local processor
      l_col1 = local_index(i  , my_pcol, np_cols, nblk, +1) ! local column number
      l_row1 = local_index(i+1, my_prow, np_rows, nblk, +1) ! first row below diagonal
      if (useGPU) then
        call obj%timer%start("copy2")
!#ifndef 
        ! set to zero maybe a kerne would be better here
          offset = (l_row1-1 + matrixRows * (l_col1-1)) * size_of_datatype
          num = (l_rows-l_row1+1)*size_of_datatype
          successGPU = gpu_memset(a_dev+offset, 0, num)
          call check_memcpy_GPU_f("elpa_cholesky: memset a_dev", 1450,  successGPU)
!#else 
!          !offset = (l_row1-1 + matrixRows * (l_col1-1)) * size_of_datatype
!          !num = (l_rows-l_row1+1)*size_of_datatype
!#ifdef WITH_GPU_STREAMS
!          !my_stream = obj%gpu_setup%my_stream
!          !successGPU = gpu_memset_async(a_dev+offset, 0, num, my_stream)
!          !call check_memcpy_GPU_f("elpa_cholesky: memset a_dev", 1458,  successGPU)
!
!          !successGPU = gpu_stream_synchronize(my_stream)
!          !call check_memcpy_GPU_f("elpa_cholesky: memset", 1461,  successGPU)
!#else
!          !successGPU = gpu_memset(a_dev+offset, 0, num)
!          !call check_memcpy_GPU_f("elpa_cholesky: memset tmp1_dev", 1464,  successGPU)
!#endif
!
!
!!#if defined(WITH_NVIDIA_CUSOLVER) || defined(WITH_AMD_ROCSOLVER)
!!        successGPU = gpu_memcpy(int(loc(a_tmp(1,1)),kind=c_intptr_t), a_dev, &
!!                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
!!        call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a_tmp", 1471,  successGPU)
!!#endif
!!        a_tmp(l_row1:l_rows,l_col1) = 0
!!#if defined(WITH_NVIDIA_CUSOLVER) || defined(WITH_AMD_ROCSOLVER)
!!        successGPU = gpu_memcpy(a_dev, int(loc(a_tmp(1,1)),kind=c_intptr_t), &
!!                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
!!        call check_memcpy_GPU_f("elpa_cholesky: memcpy a_tmp-> a_dev", 1477,  successGPU)
!#endif

!#endif 

  call obj%timer%stop("copy2")

      else ! useGPU
      endif ! useGPU
    endif ! (my_pcol==pcol(i, nblk, np_cols))
  enddo
  !endif ! useGPU

  call obj%timer%stop("loop2")


  call obj%timer%start("cleanup")
!#if !defined(WITH_NVIDIA_CUSOLVER) && !defined(WITH_AMD_ROCSOLVER)
!    num = matrixRows*matrixCols* size_of_datatype
!#ifdef WITH_GPU_STREAMS
!    my_stream = obj%gpu_setup%my_stream
!    call gpu_memcpy_async_and_stream_synchronize &
!    ("elpa_choleksky: a_tmp -> a_dev", a_dev, 0_c_intptr_t, &
!                                    a_tmp(1:obj%local_nrows,1:obj%local_ncols), &
!                                    1, 1, num, gpuMemcpyHostToDevice, my_stream, .false., .true., .false.)
!#else
!    successGPU = gpu_memcpy(a_dev, int(loc(a_tmp(1,1)),kind=c_intptr_t), &
!                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
!    call check_memcpy_GPU_f("elpa_cholesky: memcpy a_tmp-> a_dev", 1525,  successGPU)
!#endif
!#endif 

!#if !defined(WITH_NVIDIA_CUSOLVER) && !defined(WITH_AMD_ROCSOLVER)

    deallocate(a_tmp, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("elpa_cholesky: a_tmp", 1536,  istat,  errorMessage)
!#endif   

  ! restore original OpenMP settings
  call obj%timer%stop("cleanup")
  call obj%timer%stop("elpa_cholesky_&
  &real&
  &_&
  &single&
  &"//gpuString)

    end function elpa_cholesky_d_ptr_real_single_impl























!> \brief  elpa_cholesky_a_h_a_complex_double_impl: Cholesky factorization of a double-precision complex hermitian matrix
!> \details
!> \param  obj                    elpa_t object contains:
!> \param     - obj%na            Order of matrix
!> \param     - obj%local_nrows   Leading dimension of a
!> \param     - obj%local_ncols   local columns of matrix a
!> \param     - obj%nblk          blocksize of cyclic distribution, must be the same in both directions!
!> \param     - obj%mpi_comm_rows MPI communicator for rows
!> \param     - obj%mpi_comm_cols MPI communicator for columns
!> \param     - obj%wantDebug     logical, more debug information on failure
!> \param  a(lda,matrixCols)      Distributed matrix which should be inverted
!>                                Distribution is like in Scalapack.
!>                                Only upper triangle needs to be set.
!>                                The lower triangle is not referenced.
!> \result succes                 logical, reports success or failure
    function elpa_cholesky_a_h_a_complex_double_impl(obj, a) result(success)

!    This file is part of ELPA.
!
!    The ELPA library was originally created by the ELPA consortium,
!    consisting of the following organizations:
!
!    - Max Planck Computing and Data Facility (MPCDF), formerly known as
!      Rechenzentrum Garching der Max-Planck-Gesellschaft (RZG),
!    - Bergische Universität Wuppertal, Lehrstuhl für angewandte
!      Informatik,
!    - Technische Universität München, Lehrstuhl für Informatik mit
!      Schwerpunkt Wissenschaftliches Rechnen ,
!    - Fritz-Haber-Institut, Berlin, Abt. Theorie,
!    - Max-Plack-Institut für Mathematik in den Naturwissenschaften,
!      Leipzig, Abt. Komplexe Strukutren in Biologie und Kognition,
!      and
!    - IBM Deutschland GmbH
!
!    This particular source code file contains additions, changes and
!    enhancements authored by Intel Corporation which is not part of
!    the ELPA consortium.
!
!    More information can be found here:
!    http://elpa.mpcdf.mpg.de/
!
!    ELPA is free software: you can redistribute it and/or modify
!    it under the terms of the version 3 of the license of the
!    GNU Lesser General Public License as published by the Free
!    Software Foundation.
!
!    ELPA is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with ELPA.  If not, see <http://www.gnu.org/licenses/>
!
!    ELPA reflects a substantial effort on the part of the original
!    ELPA consortium, and we ask you to respect the spirit of the
!    license that we chose: i.e., please contribute any changes you
!    may have back to the original ELPA library distribution, and keep
!    any derivatives of ELPA under the same license that we chose for
!    the original distribution, the GNU Lesser General Public License.




!cannot use "../src/cholesky/./../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)



  use elpa1_compute
  use elpa_utilities
  use elpa_mpi
  use precision
  use elpa_abstract_impl
  use elpa_omp
  use elpa_blas_interfaces
  use elpa_gpu
  use mod_check_for_gpu
  use invert_trm_gpu !, only : gpu_copy_&
                     !        &double&
                     !        &_tmp1_tmp2, &
                     !        gpu_copy_&
                     !        &double&
                     !        &_a_tmp1
  use cholesky_gpu
  use mod_query_gpu_usage
  implicit none
!    Copyright 2011, A. Marek
!
!    This file is part of ELPA.
!
!    The ELPA library was originally created by the ELPA consortium,
!    consisting of the following organizations:
!
!    - Max Planck Computing and Data Facility (MPCDF), formerly known as
!      Rechenzentrum Garching der Max-Planck-Gesellschaft (RZG),
!    - Bergische Universität Wuppertal, Lehrstuhl für angewandte
!      Informatik,
!    - Technische Universität München, Lehrstuhl für Informatik mit
!      Schwerpunkt Wissenschaftliches Rechnen ,
!    - Fritz-Haber-Institut, Berlin, Abt. Theorie,
!    - Max-Plack-Institut für Mathematik in den Naturwissenschaften,
!      Leipzig, Abt. Komplexe Strukutren in Biologie und Kognition,
!      and
!    - IBM Deutschland GmbH
!
!    This particular source code file contains additions, changes and
!    enhancements authored by Intel Corporation which is not part of
!    the ELPA consortium.
!
!    More information can be found here:
!    http://elpa.mpcdf.mpg.de/
!
!    ELPA is free software: you can redistribute it and/or modify
!    it under the terms of the version 3 of the license of the
!    GNU Lesser General Public License as published by the Free
!    Software Foundation.
!
!    ELPA is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with ELPA.  If not, see <http://www.gnu.org/licenses/>
!
!    ELPA reflects a substantial effort on the part of the original
!    ELPA consortium, and we ask you to respect the spirit of the
!    license that we chose: i.e., please contribute any changes you
!    may have back to the original ELPA library distribution, and keep
!    any derivatives of ELPA under the same license that we chose for
!    the original distribution, the GNU Lesser General Public License.
!
!

  integer, parameter :: rk = C_DOUBLE
  integer, parameter :: ck = C_DOUBLE_COMPLEX
  integer, parameter :: rck = C_DOUBLE_COMPLEX
  complex(kind=rck), parameter     :: ZERO = (0.0_rk,0.0_rk), ONE = (1.0_rk,0.0_rk)
  class(elpa_abstract_impl_t), intent(inout) :: obj
  integer(kind=ik)                           :: na, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols, mpi_comm_all
  complex(kind=rck)                    :: a(obj%local_nrows,*)
  integer(kind=ik)                           :: my_prow, my_pcol, np_rows, np_cols, myid
  integer(kind=MPI_KIND)                     :: mpierr, my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI, myidMPI
  integer(kind=ik)                           :: l_cols, l_rows, l_col1, l_row1, l_colx, l_rowx
  integer(kind=ik)                           :: n, nc, i, info
  integer(kind=BLAS_KIND)                    :: infoBLAS
  integer(kind=ik)                           :: lcs, lce, lrs, lre
  integer(kind=ik)                           :: tile_size, l_rows_tile, l_cols_tile

  complex(kind=rck), allocatable       :: tmp1(:), tmp2(:,:), tmatr(:,:), tmatc(:,:)
  logical                                    :: wantDebug
  logical                                    :: success
  integer(kind=ik)                           :: istat, debug, error
  character(200)                             :: errorMessage
  integer(kind=ik)                           :: nrThreads, limitThreads
  character(20)                              :: gpuString
  logical                                    :: successGPU
  logical                                    :: useGPU
  logical                                    :: useCCL
  integer(kind=c_int)                        :: numGPU
  integer(kind=c_intptr_t)                   :: num
  integer(kind=c_intptr_t)                   :: tmp1_dev, tmatc_dev, tmatr_dev, a_dev, tmp2_dev
  integer(kind=c_intptr_t)                   :: info_dev, info_abs_accumulated_dev
  type(c_ptr)                                :: tmp1_mpi_dev
  complex(kind=rck), pointer           :: tmp1_mpi_fortran_ptr(:,:)
  integer(kind=c_intptr_t)                   :: a_off, tmatc_off, tmatr_off
  type(c_ptr)                                :: tmatc_mpi_dev
  complex(kind=rck), pointer           :: tmatc_mpi_fortran_ptr(:,:)

  integer(kind=c_intptr_t), parameter        :: size_of_datatype = size_of_&
                                                            &double&
                                                            &_&
                                                            &complex

  integer(kind=c_intptr_t)                   :: gpublasHandle, gpusolverHandle, my_stream, offset
  integer(kind=c_int)                        :: gpu_cholesky
  integer(kind=ik)                           :: blocking

  success = .true.
  useGPU = .false.
  useCCL = .false.



  ! check whether the above setting should be overriden
  if (obj%is_set("gpu_cholesky") == 1) then
    call obj%get("gpu_cholesky", gpu_cholesky, error)
    if (error .ne. ELPA_OK) then
      print *,"Problem getting option for gpu_cholesky. Aborting..."
      stop 1
    endif
    if (useGPU .and. gpu_cholesky .eq. 0) then
      useGPU = .false.
    else if (.not.(useGPU) .and. gpu_cholesky .eq. 1) then
      useGPU = .true.
    else 
    endif
  else 
    ! no override by user 
  endif

  if (.not.(useGPU)) then
  endif

  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  
  call obj%timer%start("elpa_cholesky_&
  &complex&
  &_&
  &double&
  &"//gpuString)

  nrThreads=1

  na         = obj%na
  matrixRows = obj%local_nrows
  nblk       = obj%nblk
  matrixCols = obj%local_ncols

  call obj%get("debug",debug,error)
  if (error .ne. ELPA_OK) then
    write(error_unit,*) "ELPA_CHOLESKY: Problem getting option for debug settings. Aborting..."
    success = .false.
    return
  endif
  if (debug == 1) then
    wantDebug = .true.
  else
    wantDebug = .false.
  endif

  mpi_comm_all    = obj%mpi_setup%mpi_comm_parent
  mpi_comm_cols   = obj%mpi_setup%mpi_comm_cols
  mpi_comm_rows   = obj%mpi_setup%mpi_comm_rows
    
  myid    = obj%mpi_setup%myRank_comm_parent
  my_prow = obj%mpi_setup%myRank_comm_rows
  my_pcol = obj%mpi_setup%myRank_comm_cols
        
  np_rows = obj%mpi_setup%nRanks_comm_rows
  np_cols = obj%mpi_setup%nRanks_comm_cols

  success = .true.

  ! Matrix is split into tiles; work is done only for tiles on the diagonal or above

  call obj%timer%start("prepare")

  if (obj%is_set("blocking_in_cholesky") == 1) then
    call obj%get("blocking_in_cholesky", blocking, error)
    if (error .ne. ELPA_OK) then
      write(error_unit,*) "ELPA_CHOLESKY: Problem in getting keyword 'blocking_in_cholesky'. Aborting..."
      stop 1
    endif
  else
    if (useGPU) then
      blocking = 1024
    else
      blocking = 128
    endif
    call obj%set("blocking_in_cholesky", blocking, error)
    if (error .ne. ELPA_OK) then
      write(error_unit,*) "ELPA_CHOLESKY: Problem in setting keyword 'blocking_in_cholesky'. Aborting..."
      stop 1
    endif
  endif 


  tile_size = nblk*least_common_multiple(np_rows,np_cols) ! minimum global tile size
  tile_size = ((blocking*max(np_rows,np_cols)-1)/tile_size+1)*tile_size ! make local tiles at least 128 wide

  l_rows_tile = tile_size/np_rows ! local rows of a tile
  l_cols_tile = tile_size/np_cols ! local cols of a tile

  l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a
  l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a

  if (useGPU) then
    call obj%timer%start("check_for_gpu")
    if (check_for_gpu(obj, myid, numGPU)) then
      ! set the neccessary parameters
      call set_gpu_parameters()
    else
      print *,"ELPA_CHOLESKY: GPUs are requested but not detected! Aborting..."
      success = .false.
      return
    endif
    call obj%timer%stop("check_for_gpu")
  else ! useGPU
  endif ! useGPU


  if (useGPU) then

    successGPU = gpu_malloc(tmp1_dev, nblk*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmp1_dev", 335,  successGPU)

    successGPU = gpu_memset(tmp1_dev, 0, nblk*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmp1_dev", 346,  successGPU)

    successGPU = gpu_malloc(tmp2_dev, nblk*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmp2_dev", 350,  successGPU)

    successGPU = gpu_memset(tmp2_dev, 0, nblk*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmp2_dev", 361,  successGPU)

    successGPU = gpu_malloc(tmatc_dev, l_cols*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmatc_dev", 365,  successGPU)

    successGPU = gpu_memset(tmatc_dev, 0, l_cols*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmatc_dev", 376,  successGPU)

    successGPU = gpu_malloc(tmatr_dev, l_rows*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmatr_dev", 380,  successGPU)

    successGPU = gpu_memset(tmatr_dev, 0, l_rows*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmatr_dev", 391,  successGPU)

    successGPU = gpu_malloc(a_dev, matrixRows*matrixCols*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: a_dev", 396,  successGPU)
  endif ! useGPU


  allocate(tmp1(nblk*nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmp1", 441,  istat,  errorMessage)

  allocate(tmp2(nblk,nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmp2", 444,  istat,  errorMessage)



  tmp1 = 0
  tmp2 = 0

  allocate(tmatr(l_rows,nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmatr", 460,  istat,  errorMessage)

  allocate(tmatc(l_cols,nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmatc", 463,  istat,  errorMessage)


  tmatr = 0
  tmatc = 0

  if (useGPU) then
    num = matrixRows*matrixCols* size_of_datatype
    successGPU = gpu_memcpy(a_dev, int(loc(a(1,1)),kind=c_intptr_t), &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
    call check_memcpy_GPU_f("elpa_cholesky 1: memcpy a-> a_dev", 495,  successGPU)
  endif


  call obj%timer%stop("prepare")
  call obj%timer%start("loop1")
  do n = 1, na, nblk


    ! Calculate first local row and column of the still remaining matrix
    ! on the local processor

    l_row1 = local_index(n, my_prow, np_rows, nblk, +1)
    l_col1 = local_index(n, my_pcol, np_cols, nblk, +1)

    l_rowx = local_index(n+nblk, my_prow, np_rows, nblk, +1)
    l_colx = local_index(n+nblk, my_pcol, np_cols, nblk, +1)

    if (n+nblk > na) then

      ! This is the last step, just do a Cholesky-Factorization
      ! of the remaining block

      if (useGPU) then
        if (my_prow==prow(n, nblk, np_rows) .and. my_pcol==pcol(n, nblk, np_cols)) then

          call obj%timer%start("lapack")
          num = matrixRows*matrixCols* size_of_datatype
          successGPU = gpu_memcpy(int(loc(a(1,1)),kind=c_intptr_t), a_dev,  &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a", 582,  successGPU)

          call ZPOTRF('U', int(na-n+1,kind=BLAS_KIND), a(l_row1,l_col1), &
                             int(matrixRows,kind=BLAS_KIND), infoBLAS )
          info = int(infoBLAS,kind=ik)

          num = matrixRows*matrixCols* size_of_datatype
          successGPU = gpu_memcpy(a_dev, int(loc(a(1,1)),kind=c_intptr_t), &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a", 600,  successGPU)

          call obj%timer%stop("lapack")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &complex&

            &: Error in zpotrf: ",info
            success = .false. ! "
            return
          endif ! info
        endif ! (my_prow==prow(n, nblk, np_rows) .and. my_pcol==pcol(n, nblk, np_cols))

      else ! useGPU
        if (my_prow==prow(n, nblk, np_rows) .and. my_pcol==pcol(n, nblk, np_cols)) then
          call obj%timer%start("blas")

          call ZPOTRF('U', int(na-n+1,kind=BLAS_KIND), a(l_row1,l_col1), &
                             int(matrixRows,kind=BLAS_KIND), infoBLAS )
          info = int(infoBLAS,kind=ik)
          call obj%timer%stop("blas")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &complex&

            &: Error in zpotrf: ",info
            success = .false. ! "
            return
          endif

        endif
      endif ! useGPU
      
      exit ! Loop
    endif ! (n+nblk > na) 

    ! This is not the last step

    if (my_prow==prow(n, nblk, np_rows)) then

      if (my_pcol==pcol(n, nblk, np_cols)) then

        if (useGPU) then
          call obj%timer%start("lapack")

          num = matrixRows*matrixCols* size_of_datatype

          successGPU = gpu_memcpy(int(loc(a(1,1)),kind=c_intptr_t), a_dev,  &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a", 742,  successGPU)


          call ZPOTRF('U', int(nblk,kind=BLAS_KIND), a(l_row1,l_col1), &
                               int(matrixRows,kind=BLAS_KIND) , infoBLAS )
          info = int(infoBLAS,kind=ik)

          num = matrixRows*matrixCols* size_of_datatype
          successGPU = gpu_memcpy(a_dev, int(loc(a(1,1)),kind=c_intptr_t), &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a", 769,  successGPU)

          call obj%timer%stop("lapack")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &complex&

            &: Error in zpotrf: ",info
            success = .false. ! "
            return
          endif ! info
        else ! useGPU
          ! The process owning the upper left remaining block does the
          ! Cholesky-Factorization of this block
          call obj%timer%start("blas")

          call ZPOTRF('U', int(nblk,kind=BLAS_KIND), a(l_row1,l_col1), &
                               int(matrixRows,kind=BLAS_KIND) , infoBLAS )
          info = int(infoBLAS,kind=ik)
          call obj%timer%stop("blas")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &complex&

            &: Error in zpotrf 2: ",info

            success = .false. ! "
            return
          endif ! info/=0
        endif ! useGPU
 
        if (useGPU) then
          my_stream = obj%gpu_setup%my_stream
          call gpu_copy_double_complex_a_tmp1 (a_dev, tmp1_dev, l_row1, l_col1, matrixRows, nblk, my_stream)
        else ! useGPU
          nc = 0
          do i=1,nblk
            tmp1(nc+1:nc+i) = a(l_row1:l_row1+i-1,l_col1+i-1)
            nc = nc+i
          enddo
        endif ! useGPU
      endif ! (my_pcol==pcol(n, nblk, np_cols))


      if (useGPU) then
        my_stream = obj%gpu_setup%my_stream
        call gpu_copy_double_complex_tmp1_tmp2 (tmp1_dev, tmp2_dev, nblk, nblk, my_stream)
      else ! useGPU
        nc = 0
        do i=1,nblk
          tmp2(1:i,i) = tmp1(nc+1:nc+i)
          nc = nc+i
        enddo
      endif ! useGPU

      if (useGPU) then
        call obj%timer%start("gpublas")
        gpublasHandle = obj%gpu_setup%gpublasHandleArray(0)
        if (l_cols-l_colx+1 > 0) then
          a_off = (l_row1-1 + (l_colx-1)*matrixRows) * size_of_datatype
          call gpublas_ZTRSM('L', 'U', 'C', 'N', nblk, l_cols-l_colx+1, ONE, &
                            tmp2_dev, nblk, a_dev+a_off, matrixRows, gpublasHandle)
          if (wantDebug .and. gpu_vendor() /= SYCL_GPU) successGPU = gpu_DeviceSynchronize()
        endif
        call obj%timer%stop("gpublas")

      else ! useGPU

        call obj%timer%start("lapack")
        if (l_cols-l_colx+1>0) then
          call ZTRSM('L', 'U', 'C', 'N', int(nblk,kind=BLAS_KIND),  &
                            int(l_cols-l_colx+1,kind=BLAS_KIND), ONE, tmp2, &
                            int(ubound(tmp2,dim=1),kind=BLAS_KIND), a(l_row1,l_colx), int(matrixRows,kind=BLAS_KIND) )
        endif
        call obj%timer%stop("lapack")
      endif ! useGPU
    endif ! (my_prow==prow(n, nblk, np_rows))


    if (useGPU) then
      if (my_prow==prow(n, nblk, np_rows)) then
        ! if l_cols-l_colx+1 == 0 kernel launch with 0 blocks => raises error
        call obj%timer%start("a_tmatc_kernel")
        if (l_cols-l_colx+1>0) then
           my_stream = obj%gpu_setup%my_stream
           call gpu_copy_double_complex_a_tmatc(a_dev, tmatc_dev, nblk, matrixRows, l_cols, l_colx, l_row1, my_stream)
        endif
        call obj%timer%stop("a_tmatc_kernel")
      endif
    else ! useGPU
      do i=1,nblk
        if (my_prow==prow(n, nblk, np_rows)) tmatc(l_colx:l_cols,i) = conjg(a(l_row1+i-1,l_colx:l_cols))
      enddo
    endif ! useGPU

    if (useGPU  .and. .not. useCCL) then
      num = l_cols*nblk*size_of_datatype
      successGPU = gpu_memcpy(int(loc(tmatc),kind=c_intptr_t), tmatc_dev, num, &
                              gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("elpa_cholesky: tmatc_dev to tmatc", 1201,  successGPU)
    endif


    if (useGPU .and. .not. useCCL) then

      ! not needed if transpose vec CCL
      num = l_rows*nblk*size_of_datatype
      successGPU = gpu_memcpy(int(loc(tmatr),kind=c_intptr_t), tmatr_dev, num, &
                              gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("elpa_cholesky: tmatr_dev to tmatr", 1220,  successGPU)
    endif ! (useGPU .and. .not. useCCL)


    if (useCCL) then

    else ! useCCL
      call elpa_transpose_vectors_&
      &complex&
      &_&
      &double &
      (obj, tmatc, ubound(tmatc,dim=1), mpi_comm_cols, &
      tmatr, ubound(tmatr,dim=1), mpi_comm_rows, &
      n, na, nblk, nblk, nrThreads, .false., success)
    endif ! useCCL

    if (useGPU .and. .not. useCCL) then

      num = l_rows*nblk*size_of_datatype
      successGPU = gpu_memcpy(tmatr_dev, int(loc(tmatr),kind=c_intptr_t), num, &
                              gpuMemcpyHostToDevice)
      call check_memcpy_GPU_f("elpa_cholesky: tmat to tmatr_dev", 1271,  successGPU)
    endif

    if (useGPU) then
      do i=0,(na-1)/tile_size
        lcs = max(l_colx,i*l_cols_tile+1)
        lce = min(l_cols,(i+1)*l_cols_tile)
        lrs = l_rowx
        lre = min(l_rows,(i+1)*l_rows_tile)
        if (lce  < lcs .or. lre < lrs) cycle
        call obj%timer%start("gpublas")
        gpublasHandle = obj%gpu_setup%gpublasHandleArray(0)
        tmatr_off = (lrs-1 + (1-1)*l_rows) * size_of_datatype
        tmatc_off = (lcs-1 + (1-1)*l_cols) * size_of_datatype
        a_off = (lrs-1 + (lcs-1)*matrixRows) * size_of_datatype
        call gpublas_ZGEMM('N', 'C', lre-lrs+1, lce-lcs+1, nblk, &
                            -ONE, tmatr_dev+tmatr_off, l_rows, tmatc_dev+tmatc_off, l_cols, ONE, &
                            a_dev+a_off, matrixRows, gpublasHandle)
        if (wantDebug .and. gpu_vendor() /= SYCL_GPU) successGPU = gpu_DeviceSynchronize()
        call obj%timer%stop("gpublas")
      enddo
    else !useGPU
      do i=0,(na-1)/tile_size
        lcs = max(l_colx,i*l_cols_tile+1)
        lce = min(l_cols,(i+1)*l_cols_tile)
        lrs = l_rowx
        lre = min(l_rows,(i+1)*l_rows_tile)
        if (lce<lcs .or. lre<lrs) cycle
        call obj%timer%start("blas")
        call ZGEMM('N', 'C', int(lre-lrs+1,kind=BLAS_KIND), int(lce-lcs+1,kind=BLAS_KIND), &
                            int(nblk,kind=BLAS_KIND), -ONE,  &
                            tmatr(lrs,1), int(ubound(tmatr,dim=1),kind=BLAS_KIND), tmatc(lcs,1), &
                            int(ubound(tmatc,dim=1),kind=BLAS_KIND), &
                            ONE, a(lrs,lcs), int(matrixRows,kind=BLAS_KIND))
        call obj%timer%stop("blas")
      enddo
    endif ! useGPU


  enddo ! n = 1, na, nblk

  call obj%timer%stop("loop1")
  
  call obj%timer%start("deallocate")
  if (useGPU) then



    successGPU = gpu_free(tmp1_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmp1_dev", 1351,  successGPU)

    successGPU = gpu_free(tmp2_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmp1_dev", 1354,  successGPU)

    successGPU = gpu_free(tmatc_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmatc_dev", 1357,  successGPU)

    successGPU = gpu_free(tmatr_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmatr_dev", 1360,  successGPU)

  endif


  deallocate(tmp1, tmp2, tmatr, tmatc, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("elpa_cholesky: tmp1, tmp2, tmatr, tmatc", 1384,  istat,  errorMessage)

  call obj%timer%stop("deallocate")


  call obj%timer%start("loop2")

  call obj%timer%start("copy1")

  if (useGPU) then
    num = matrixRows*matrixCols* size_of_datatype
    successGPU = gpu_memcpy(int(loc(a(1,1)),kind=c_intptr_t), a_dev,  &
                     matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
    call check_memcpy_GPU_f("elpa_cholesky: memcpy 2 a-> a_dev", 1405,  successGPU)
  endif
  ! Set the lower triangle to 0, it contains garbage (form the above matrix multiplications)

  call obj%timer%stop("copy1")
  do i=1,na
    if (my_pcol==pcol(i, nblk, np_cols)) then
      ! column i is on local processor
      l_col1 = local_index(i  , my_pcol, np_cols, nblk, +1) ! local column number
      l_row1 = local_index(i+1, my_prow, np_rows, nblk, +1) ! first row below diagonal
      if (useGPU) then
        call obj%timer%start("copy2")
!#ifndef DEVICE_POINTER
        ! set to zero maybe a kerne would be better here
          offset = (l_row1-1 + matrixRows * (l_col1-1)) * size_of_datatype
          num = (l_rows-l_row1+1)*size_of_datatype
          successGPU = gpu_memset(a_dev+offset, 0, num)
          call check_memcpy_GPU_f("elpa_cholesky: memset a_dev", 1450,  successGPU)
!#else 
!          !offset = (l_row1-1 + matrixRows * (l_col1-1)) * size_of_datatype
!          !num = (l_rows-l_row1+1)*size_of_datatype
!#ifdef WITH_GPU_STREAMS
!          !my_stream = obj%gpu_setup%my_stream
!          !successGPU = gpu_memset_async(a_dev+offset, 0, num, my_stream)
!          !call check_memcpy_GPU_f("elpa_cholesky: memset a_dev", 1458,  successGPU)
!
!          !successGPU = gpu_stream_synchronize(my_stream)
!          !call check_memcpy_GPU_f("elpa_cholesky: memset", 1461,  successGPU)
!#else
!          !successGPU = gpu_memset(a_dev+offset, 0, num)
!          !call check_memcpy_GPU_f("elpa_cholesky: memset tmp1_dev", 1464,  successGPU)
!#endif
!
!
!!#if defined(WITH_NVIDIA_CUSOLVER) || defined(WITH_AMD_ROCSOLVER)
!!        successGPU = gpu_memcpy(int(loc(a_tmp(1,1)),kind=c_intptr_t), a_dev, &
!!                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
!!        call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a_tmp", 1471,  successGPU)
!!#endif
!!        a_tmp(l_row1:l_rows,l_col1) = 0
!!#if defined(WITH_NVIDIA_CUSOLVER) || defined(WITH_AMD_ROCSOLVER)
!!        successGPU = gpu_memcpy(a_dev, int(loc(a_tmp(1,1)),kind=c_intptr_t), &
!!                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
!!        call check_memcpy_GPU_f("elpa_cholesky: memcpy a_tmp-> a_dev", 1477,  successGPU)
!#endif

!#endif 

  call obj%timer%stop("copy2")

      else ! useGPU
        a(l_row1:l_rows,l_col1) = 0
      endif ! useGPU
    endif ! (my_pcol==pcol(i, nblk, np_cols))
  enddo
  !endif ! useGPU

  call obj%timer%stop("loop2")


  call obj%timer%start("cleanup")
  if (useGPU) then
    ! copy back
    successGPU = gpu_memcpy(int(loc(a(1,1)),kind=c_intptr_t), a_dev,  &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
    call check_memcpy_GPU_f("elpa_cholesky: memcpy a-> d_dev", 1502,  successGPU)

    successGPU = gpu_free(a_dev)
    call check_dealloc_GPU_f("elpa_cholesky: a_dev", 1505,  successGPU)


  endif

  ! restore original OpenMP settings
  call obj%timer%stop("cleanup")
  call obj%timer%stop("elpa_cholesky_&
  &complex&
  &_&
  &double&
  &"//gpuString)

    end function elpa_cholesky_a_h_a_complex_double_impl




















!> \brief  elpa_cholesky_d_ptr_complex_double_impl: Cholesky factorization of a double-precision complex hermitian matrix
!> \details
!> \param  obj                    elpa_t object contains:
!> \param     - obj%na            Order of matrix
!> \param     - obj%local_nrows   Leading dimension of a
!> \param     - obj%local_ncols   local columns of matrix a
!> \param     - obj%nblk          blocksize of cyclic distribution, must be the same in both directions!
!> \param     - obj%mpi_comm_rows MPI communicator for rows
!> \param     - obj%mpi_comm_cols MPI communicator for columns
!> \param     - obj%wantDebug     logical, more debug information on failure
!> \param  aDev(lda,matrixCols)   Distributed matrix which should be inverted as type(c_ptr) living on a device
!>                                Distribution is like in Scalapack.
!>                                Only upper triangle needs to be set.
!>                                The lower triangle is not referenced.
!> \result succes                 logical, reports success or failure
    function elpa_cholesky_d_ptr_complex_double_impl(obj, aDev) result(success)

!    This file is part of ELPA.
!
!    The ELPA library was originally created by the ELPA consortium,
!    consisting of the following organizations:
!
!    - Max Planck Computing and Data Facility (MPCDF), formerly known as
!      Rechenzentrum Garching der Max-Planck-Gesellschaft (RZG),
!    - Bergische Universität Wuppertal, Lehrstuhl für angewandte
!      Informatik,
!    - Technische Universität München, Lehrstuhl für Informatik mit
!      Schwerpunkt Wissenschaftliches Rechnen ,
!    - Fritz-Haber-Institut, Berlin, Abt. Theorie,
!    - Max-Plack-Institut für Mathematik in den Naturwissenschaften,
!      Leipzig, Abt. Komplexe Strukutren in Biologie und Kognition,
!      and
!    - IBM Deutschland GmbH
!
!    This particular source code file contains additions, changes and
!    enhancements authored by Intel Corporation which is not part of
!    the ELPA consortium.
!
!    More information can be found here:
!    http://elpa.mpcdf.mpg.de/
!
!    ELPA is free software: you can redistribute it and/or modify
!    it under the terms of the version 3 of the license of the
!    GNU Lesser General Public License as published by the Free
!    Software Foundation.
!
!    ELPA is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with ELPA.  If not, see <http://www.gnu.org/licenses/>
!
!    ELPA reflects a substantial effort on the part of the original
!    ELPA consortium, and we ask you to respect the spirit of the
!    license that we chose: i.e., please contribute any changes you
!    may have back to the original ELPA library distribution, and keep
!    any derivatives of ELPA under the same license that we chose for
!    the original distribution, the GNU Lesser General Public License.




!cannot use "../src/cholesky/./../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)



  use elpa1_compute
  use elpa_utilities
  use elpa_mpi
  use precision
  use elpa_abstract_impl
  use elpa_omp
  use elpa_blas_interfaces
  use elpa_gpu
  use mod_check_for_gpu
  use invert_trm_gpu !, only : gpu_copy_&
                     !        &double&
                     !        &_tmp1_tmp2, &
                     !        gpu_copy_&
                     !        &double&
                     !        &_a_tmp1
  use cholesky_gpu
  use mod_query_gpu_usage
  implicit none
!    Copyright 2011, A. Marek
!
!    This file is part of ELPA.
!
!    The ELPA library was originally created by the ELPA consortium,
!    consisting of the following organizations:
!
!    - Max Planck Computing and Data Facility (MPCDF), formerly known as
!      Rechenzentrum Garching der Max-Planck-Gesellschaft (RZG),
!    - Bergische Universität Wuppertal, Lehrstuhl für angewandte
!      Informatik,
!    - Technische Universität München, Lehrstuhl für Informatik mit
!      Schwerpunkt Wissenschaftliches Rechnen ,
!    - Fritz-Haber-Institut, Berlin, Abt. Theorie,
!    - Max-Plack-Institut für Mathematik in den Naturwissenschaften,
!      Leipzig, Abt. Komplexe Strukutren in Biologie und Kognition,
!      and
!    - IBM Deutschland GmbH
!
!    This particular source code file contains additions, changes and
!    enhancements authored by Intel Corporation which is not part of
!    the ELPA consortium.
!
!    More information can be found here:
!    http://elpa.mpcdf.mpg.de/
!
!    ELPA is free software: you can redistribute it and/or modify
!    it under the terms of the version 3 of the license of the
!    GNU Lesser General Public License as published by the Free
!    Software Foundation.
!
!    ELPA is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with ELPA.  If not, see <http://www.gnu.org/licenses/>
!
!    ELPA reflects a substantial effort on the part of the original
!    ELPA consortium, and we ask you to respect the spirit of the
!    license that we chose: i.e., please contribute any changes you
!    may have back to the original ELPA library distribution, and keep
!    any derivatives of ELPA under the same license that we chose for
!    the original distribution, the GNU Lesser General Public License.
!
!

  integer, parameter :: rk = C_DOUBLE
  integer, parameter :: ck = C_DOUBLE_COMPLEX
  integer, parameter :: rck = C_DOUBLE_COMPLEX
  complex(kind=rck), parameter     :: ZERO = (0.0_rk,0.0_rk), ONE = (1.0_rk,0.0_rk)
  class(elpa_abstract_impl_t), intent(inout) :: obj
  integer(kind=ik)                           :: na, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols, mpi_comm_all
  type(c_ptr)                                :: aDev
!#if !defined(WITH_NVIDIA_CUSOLVER) && !defined(WITH_AMD_ROCSOLVER)
  complex(kind=rck), allocatable       :: a_tmp(:,:)
!#endif
  integer(kind=ik)                           :: my_prow, my_pcol, np_rows, np_cols, myid
  integer(kind=MPI_KIND)                     :: mpierr, my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI, myidMPI
  integer(kind=ik)                           :: l_cols, l_rows, l_col1, l_row1, l_colx, l_rowx
  integer(kind=ik)                           :: n, nc, i, info
  integer(kind=BLAS_KIND)                    :: infoBLAS
  integer(kind=ik)                           :: lcs, lce, lrs, lre
  integer(kind=ik)                           :: tile_size, l_rows_tile, l_cols_tile

  complex(kind=rck), allocatable       :: tmp1(:), tmp2(:,:), tmatr(:,:), tmatc(:,:)
  logical                                    :: wantDebug
  logical                                    :: success
  integer(kind=ik)                           :: istat, debug, error
  character(200)                             :: errorMessage
  integer(kind=ik)                           :: nrThreads, limitThreads
  character(20)                              :: gpuString
  logical                                    :: successGPU
  logical                                    :: useGPU
  logical                                    :: useCCL
  integer(kind=c_int)                        :: numGPU
  integer(kind=c_intptr_t)                   :: num
  integer(kind=c_intptr_t)                   :: tmp1_dev, tmatc_dev, tmatr_dev, a_dev, tmp2_dev
  integer(kind=c_intptr_t)                   :: info_dev, info_abs_accumulated_dev
  type(c_ptr)                                :: tmp1_mpi_dev
  complex(kind=rck), pointer           :: tmp1_mpi_fortran_ptr(:,:)
  integer(kind=c_intptr_t)                   :: a_off, tmatc_off, tmatr_off
  type(c_ptr)                                :: tmatc_mpi_dev
  complex(kind=rck), pointer           :: tmatc_mpi_fortran_ptr(:,:)

  integer(kind=c_intptr_t), parameter        :: size_of_datatype = size_of_&
                                                            &double&
                                                            &_&
                                                            &complex

  integer(kind=c_intptr_t)                   :: gpublasHandle, gpusolverHandle, my_stream, offset
  integer(kind=c_int)                        :: gpu_cholesky
  integer(kind=ik)                           :: blocking

  success = .true.
  useGPU = .false.
  useCCL = .false.

  useGPU = .true.

  if (.not.(useGPU)) then
    print *,"You used the interface for device pointers for elpa_cholesky but did not specify GPU usage!. Aborting..."
    stop 1
  endif

  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  
  call obj%timer%start("elpa_cholesky_&
  &complex&
  &_&
  &double&
  &"//gpuString)

  nrThreads=1

  na         = obj%na
  matrixRows = obj%local_nrows
  nblk       = obj%nblk
  matrixCols = obj%local_ncols

  call obj%get("debug",debug,error)
  if (error .ne. ELPA_OK) then
    write(error_unit,*) "ELPA_CHOLESKY: Problem getting option for debug settings. Aborting..."
    success = .false.
    return
  endif
  if (debug == 1) then
    wantDebug = .true.
  else
    wantDebug = .false.
  endif

  mpi_comm_all    = obj%mpi_setup%mpi_comm_parent
  mpi_comm_cols   = obj%mpi_setup%mpi_comm_cols
  mpi_comm_rows   = obj%mpi_setup%mpi_comm_rows
    
  myid    = obj%mpi_setup%myRank_comm_parent
  my_prow = obj%mpi_setup%myRank_comm_rows
  my_pcol = obj%mpi_setup%myRank_comm_cols
        
  np_rows = obj%mpi_setup%nRanks_comm_rows
  np_cols = obj%mpi_setup%nRanks_comm_cols

  success = .true.

  ! Matrix is split into tiles; work is done only for tiles on the diagonal or above

  call obj%timer%start("prepare")

  if (obj%is_set("blocking_in_cholesky") == 1) then
    call obj%get("blocking_in_cholesky", blocking, error)
    if (error .ne. ELPA_OK) then
      write(error_unit,*) "ELPA_CHOLESKY: Problem in getting keyword 'blocking_in_cholesky'. Aborting..."
      stop 1
    endif
  else
    if (useGPU) then
      blocking = 1024
    else
      blocking = 128
    endif
    call obj%set("blocking_in_cholesky", blocking, error)
    if (error .ne. ELPA_OK) then
      write(error_unit,*) "ELPA_CHOLESKY: Problem in setting keyword 'blocking_in_cholesky'. Aborting..."
      stop 1
    endif
  endif 


  tile_size = nblk*least_common_multiple(np_rows,np_cols) ! minimum global tile size
  tile_size = ((blocking*max(np_rows,np_cols)-1)/tile_size+1)*tile_size ! make local tiles at least 128 wide

  l_rows_tile = tile_size/np_rows ! local rows of a tile
  l_cols_tile = tile_size/np_cols ! local cols of a tile

  l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a
  l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a

  if (useGPU) then
    call obj%timer%start("check_for_gpu")
    if (check_for_gpu(obj, myid, numGPU)) then
      ! set the neccessary parameters
      call set_gpu_parameters()
    else
      print *,"ELPA_CHOLESKY: GPUs are requested but not detected! Aborting..."
      success = .false.
      return
    endif
    call obj%timer%stop("check_for_gpu")
  else ! useGPU
  endif ! useGPU


  if (useGPU) then

    successGPU = gpu_malloc(tmp1_dev, nblk*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmp1_dev", 335,  successGPU)

    successGPU = gpu_memset(tmp1_dev, 0, nblk*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmp1_dev", 346,  successGPU)

    successGPU = gpu_malloc(tmp2_dev, nblk*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmp2_dev", 350,  successGPU)

    successGPU = gpu_memset(tmp2_dev, 0, nblk*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmp2_dev", 361,  successGPU)

    successGPU = gpu_malloc(tmatc_dev, l_cols*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmatc_dev", 365,  successGPU)

    successGPU = gpu_memset(tmatc_dev, 0, l_cols*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmatc_dev", 376,  successGPU)

    successGPU = gpu_malloc(tmatr_dev, l_rows*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmatr_dev", 380,  successGPU)

    successGPU = gpu_memset(tmatr_dev, 0, l_rows*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmatr_dev", 391,  successGPU)

    a_dev = transfer(aDev, a_dev)
!#if !defined(WITH_NVIDIA_CUSOLVER) && !defined(WITH_AMD_ROCSOLVER)
    allocate(a_tmp(obj%local_nrows,obj%local_ncols), stat=istat, errmsg=errorMessage)
    call check_allocate_f("elpa_cholesky: a_tmp", 407,  istat,  errorMessage)
!#endif 
  endif ! useGPU


  allocate(tmp1(nblk*nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmp1", 441,  istat,  errorMessage)

  allocate(tmp2(nblk,nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmp2", 444,  istat,  errorMessage)



  tmp1 = 0
  tmp2 = 0

  allocate(tmatr(l_rows,nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmatr", 460,  istat,  errorMessage)

  allocate(tmatc(l_cols,nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmatc", 463,  istat,  errorMessage)


  tmatr = 0
  tmatc = 0



  call obj%timer%stop("prepare")
  call obj%timer%start("loop1")
  do n = 1, na, nblk


    ! Calculate first local row and column of the still remaining matrix
    ! on the local processor

    l_row1 = local_index(n, my_prow, np_rows, nblk, +1)
    l_col1 = local_index(n, my_pcol, np_cols, nblk, +1)

    l_rowx = local_index(n+nblk, my_prow, np_rows, nblk, +1)
    l_colx = local_index(n+nblk, my_pcol, np_cols, nblk, +1)

    if (n+nblk > na) then

      ! This is the last step, just do a Cholesky-Factorization
      ! of the remaining block

      if (useGPU) then
        if (my_prow==prow(n, nblk, np_rows) .and. my_pcol==pcol(n, nblk, np_cols)) then

          call obj%timer%start("lapack")
          num = matrixRows*matrixCols* size_of_datatype

          successGPU = gpu_memcpy(int(loc(a_tmp(1,1)),kind=c_intptr_t), a_dev,  &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a_tmp", 615,  successGPU)

          call ZPOTRF('U', int(na-n+1,kind=BLAS_KIND), a_tmp(l_row1,l_col1), &
                             int(matrixRows,kind=BLAS_KIND), infoBLAS )
          info = int(infoBLAS,kind=ik)

          num = matrixRows*matrixCols* size_of_datatype
          successGPU = gpu_memcpy(a_dev, int(loc(a_tmp(1,1)),kind=c_intptr_t), &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a_tmp", 633,  successGPU)

          call obj%timer%stop("lapack")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &complex&

            &: Error in zpotrf: ",info
            success = .false. ! "
            return
          endif ! info
        endif ! (my_prow==prow(n, nblk, np_rows) .and. my_pcol==pcol(n, nblk, np_cols))

      else ! useGPU
        if (my_prow==prow(n, nblk, np_rows) .and. my_pcol==pcol(n, nblk, np_cols)) then
          call obj%timer%start("blas")

          info = int(infoBLAS,kind=ik)
          call obj%timer%stop("blas")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &complex&

            &: Error in zpotrf: ",info
            success = .false. ! "
            return
          endif

        endif
      endif ! useGPU
      
      exit ! Loop
    endif ! (n+nblk > na) 

    ! This is not the last step

    if (my_prow==prow(n, nblk, np_rows)) then

      if (my_pcol==pcol(n, nblk, np_cols)) then

        if (useGPU) then
          call obj%timer%start("lapack")

          num = matrixRows*matrixCols* size_of_datatype

          successGPU = gpu_memcpy(int(loc(a_tmp(1,1)),kind=c_intptr_t), a_dev,  &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a_tmp", 783,  successGPU)

          call ZPOTRF('U', int(nblk,kind=BLAS_KIND), a_tmp(l_row1,l_col1), &
                               int(matrixRows,kind=BLAS_KIND) , infoBLAS )
          info = int(infoBLAS,kind=ik)
          num = matrixRows*matrixCols* size_of_datatype
          successGPU = gpu_memcpy(a_dev, int(loc(a_tmp(1,1)),kind=c_intptr_t), &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a_tmp", 799,  successGPU)

          call obj%timer%stop("lapack")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &complex&

            &: Error in zpotrf: ",info
            success = .false. ! "
            return
          endif ! info
        else ! useGPU
          ! The process owning the upper left remaining block does the
          ! Cholesky-Factorization of this block
          call obj%timer%start("blas")

          info = int(infoBLAS,kind=ik)
          call obj%timer%stop("blas")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &complex&

            &: Error in zpotrf 2: ",info

            success = .false. ! "
            return
          endif ! info/=0
        endif ! useGPU
 
        if (useGPU) then
          my_stream = obj%gpu_setup%my_stream
          call gpu_copy_double_complex_a_tmp1 (a_dev, tmp1_dev, l_row1, l_col1, matrixRows, nblk, my_stream)
        else ! useGPU
          nc = 0
          do i=1,nblk
            nc = nc+i
          enddo
        endif ! useGPU
      endif ! (my_pcol==pcol(n, nblk, np_cols))


      if (useGPU) then
        my_stream = obj%gpu_setup%my_stream
        call gpu_copy_double_complex_tmp1_tmp2 (tmp1_dev, tmp2_dev, nblk, nblk, my_stream)
      else ! useGPU
        nc = 0
        do i=1,nblk
          tmp2(1:i,i) = tmp1(nc+1:nc+i)
          nc = nc+i
        enddo
      endif ! useGPU

      if (useGPU) then
        call obj%timer%start("gpublas")
        gpublasHandle = obj%gpu_setup%gpublasHandleArray(0)
        if (l_cols-l_colx+1 > 0) then
          a_off = (l_row1-1 + (l_colx-1)*matrixRows) * size_of_datatype
          call gpublas_ZTRSM('L', 'U', 'C', 'N', nblk, l_cols-l_colx+1, ONE, &
                            tmp2_dev, nblk, a_dev+a_off, matrixRows, gpublasHandle)
          if (wantDebug .and. gpu_vendor() /= SYCL_GPU) successGPU = gpu_DeviceSynchronize()
        endif
        call obj%timer%stop("gpublas")

      else ! useGPU

        call obj%timer%start("lapack")
        call obj%timer%stop("lapack")
      endif ! useGPU
    endif ! (my_prow==prow(n, nblk, np_rows))


    if (useGPU) then
      if (my_prow==prow(n, nblk, np_rows)) then
        ! if l_cols-l_colx+1 == 0 kernel launch with 0 blocks => raises error
        call obj%timer%start("a_tmatc_kernel")
        if (l_cols-l_colx+1>0) then
           my_stream = obj%gpu_setup%my_stream
           call gpu_copy_double_complex_a_tmatc(a_dev, tmatc_dev, nblk, matrixRows, l_cols, l_colx, l_row1, my_stream)
        endif
        call obj%timer%stop("a_tmatc_kernel")
      endif
    else ! useGPU
      do i=1,nblk
      enddo
    endif ! useGPU

    if (useGPU  .and. .not. useCCL) then
      num = l_cols*nblk*size_of_datatype
      successGPU = gpu_memcpy(int(loc(tmatc),kind=c_intptr_t), tmatc_dev, num, &
                              gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("elpa_cholesky: tmatc_dev to tmatc", 1201,  successGPU)
    endif


    if (useGPU .and. .not. useCCL) then

      ! not needed if transpose vec CCL
      num = l_rows*nblk*size_of_datatype
      successGPU = gpu_memcpy(int(loc(tmatr),kind=c_intptr_t), tmatr_dev, num, &
                              gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("elpa_cholesky: tmatr_dev to tmatr", 1220,  successGPU)
    endif ! (useGPU .and. .not. useCCL)


    if (useCCL) then

    else ! useCCL
      call elpa_transpose_vectors_&
      &complex&
      &_&
      &double &
      (obj, tmatc, ubound(tmatc,dim=1), mpi_comm_cols, &
      tmatr, ubound(tmatr,dim=1), mpi_comm_rows, &
      n, na, nblk, nblk, nrThreads, .false., success)
    endif ! useCCL

    if (useGPU .and. .not. useCCL) then

      num = l_rows*nblk*size_of_datatype
      successGPU = gpu_memcpy(tmatr_dev, int(loc(tmatr),kind=c_intptr_t), num, &
                              gpuMemcpyHostToDevice)
      call check_memcpy_GPU_f("elpa_cholesky: tmat to tmatr_dev", 1271,  successGPU)
    endif

    if (useGPU) then
      do i=0,(na-1)/tile_size
        lcs = max(l_colx,i*l_cols_tile+1)
        lce = min(l_cols,(i+1)*l_cols_tile)
        lrs = l_rowx
        lre = min(l_rows,(i+1)*l_rows_tile)
        if (lce  < lcs .or. lre < lrs) cycle
        call obj%timer%start("gpublas")
        gpublasHandle = obj%gpu_setup%gpublasHandleArray(0)
        tmatr_off = (lrs-1 + (1-1)*l_rows) * size_of_datatype
        tmatc_off = (lcs-1 + (1-1)*l_cols) * size_of_datatype
        a_off = (lrs-1 + (lcs-1)*matrixRows) * size_of_datatype
        call gpublas_ZGEMM('N', 'C', lre-lrs+1, lce-lcs+1, nblk, &
                            -ONE, tmatr_dev+tmatr_off, l_rows, tmatc_dev+tmatc_off, l_cols, ONE, &
                            a_dev+a_off, matrixRows, gpublasHandle)
        if (wantDebug .and. gpu_vendor() /= SYCL_GPU) successGPU = gpu_DeviceSynchronize()
        call obj%timer%stop("gpublas")
      enddo
    else !useGPU
      do i=0,(na-1)/tile_size
        lcs = max(l_colx,i*l_cols_tile+1)
        lce = min(l_cols,(i+1)*l_cols_tile)
        lrs = l_rowx
        lre = min(l_rows,(i+1)*l_rows_tile)
        if (lce<lcs .or. lre<lrs) cycle
        call obj%timer%start("blas")
        call obj%timer%stop("blas")
      enddo
    endif ! useGPU


  enddo ! n = 1, na, nblk

  call obj%timer%stop("loop1")
  
  call obj%timer%start("deallocate")
  if (useGPU) then



    successGPU = gpu_free(tmp1_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmp1_dev", 1351,  successGPU)

    successGPU = gpu_free(tmp2_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmp1_dev", 1354,  successGPU)

    successGPU = gpu_free(tmatc_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmatc_dev", 1357,  successGPU)

    successGPU = gpu_free(tmatr_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmatr_dev", 1360,  successGPU)

  endif


  deallocate(tmp1, tmp2, tmatr, tmatc, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("elpa_cholesky: tmp1, tmp2, tmatr, tmatc", 1384,  istat,  errorMessage)

  call obj%timer%stop("deallocate")


  call obj%timer%start("loop2")

  call obj%timer%start("copy1")

!  if (useGPU) then
!          ! is this needed: guess not
!#if !defined(WITH_NVIDIA_CUSOLVER) && !defined(WITH_AMD_ROCSOLVER)
!    num = matrixRows*matrixCols* size_of_datatype
!#ifdef WITH_GPU_STREAMS 
!    my_stream = obj%gpu_setup%my_stream
!    call gpu_memcpy_async_and_stream_synchronize &
!    ("elpa_choleksky: a_dev -> a_tmp", a_dev, 0_c_intptr_t, &
!                                    a_tmp(1:obj%local_nrows,1:obj%local_ncols), &
!                                    1, 1, num, gpuMemcpyDeviceToHost, my_stream, .false., .true., .false.)
!#else
!    successGPU = gpu_memcpy(int(loc(a_tmp(1,1)),kind=c_intptr_t), a_dev,  &
!                     matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
!    call check_memcpy_GPU_f("elpa_cholesky: memcpy 2 a-> a_dev", 1422,  successGPU)
!#endif
!#endif 
!  endif
  ! Set the lower triangle to 0, it contains garbage (form the above matrix multiplications)

  call obj%timer%stop("copy1")
  do i=1,na
    if (my_pcol==pcol(i, nblk, np_cols)) then
      ! column i is on local processor
      l_col1 = local_index(i  , my_pcol, np_cols, nblk, +1) ! local column number
      l_row1 = local_index(i+1, my_prow, np_rows, nblk, +1) ! first row below diagonal
      if (useGPU) then
        call obj%timer%start("copy2")
!#ifndef 
        ! set to zero maybe a kerne would be better here
          offset = (l_row1-1 + matrixRows * (l_col1-1)) * size_of_datatype
          num = (l_rows-l_row1+1)*size_of_datatype
          successGPU = gpu_memset(a_dev+offset, 0, num)
          call check_memcpy_GPU_f("elpa_cholesky: memset a_dev", 1450,  successGPU)
!#else 
!          !offset = (l_row1-1 + matrixRows * (l_col1-1)) * size_of_datatype
!          !num = (l_rows-l_row1+1)*size_of_datatype
!#ifdef WITH_GPU_STREAMS
!          !my_stream = obj%gpu_setup%my_stream
!          !successGPU = gpu_memset_async(a_dev+offset, 0, num, my_stream)
!          !call check_memcpy_GPU_f("elpa_cholesky: memset a_dev", 1458,  successGPU)
!
!          !successGPU = gpu_stream_synchronize(my_stream)
!          !call check_memcpy_GPU_f("elpa_cholesky: memset", 1461,  successGPU)
!#else
!          !successGPU = gpu_memset(a_dev+offset, 0, num)
!          !call check_memcpy_GPU_f("elpa_cholesky: memset tmp1_dev", 1464,  successGPU)
!#endif
!
!
!!#if defined(WITH_NVIDIA_CUSOLVER) || defined(WITH_AMD_ROCSOLVER)
!!        successGPU = gpu_memcpy(int(loc(a_tmp(1,1)),kind=c_intptr_t), a_dev, &
!!                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
!!        call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a_tmp", 1471,  successGPU)
!!#endif
!!        a_tmp(l_row1:l_rows,l_col1) = 0
!!#if defined(WITH_NVIDIA_CUSOLVER) || defined(WITH_AMD_ROCSOLVER)
!!        successGPU = gpu_memcpy(a_dev, int(loc(a_tmp(1,1)),kind=c_intptr_t), &
!!                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
!!        call check_memcpy_GPU_f("elpa_cholesky: memcpy a_tmp-> a_dev", 1477,  successGPU)
!#endif

!#endif 

  call obj%timer%stop("copy2")

      else ! useGPU
      endif ! useGPU
    endif ! (my_pcol==pcol(i, nblk, np_cols))
  enddo
  !endif ! useGPU

  call obj%timer%stop("loop2")


  call obj%timer%start("cleanup")
!#if !defined(WITH_NVIDIA_CUSOLVER) && !defined(WITH_AMD_ROCSOLVER)
!    num = matrixRows*matrixCols* size_of_datatype
!#ifdef WITH_GPU_STREAMS
!    my_stream = obj%gpu_setup%my_stream
!    call gpu_memcpy_async_and_stream_synchronize &
!    ("elpa_choleksky: a_tmp -> a_dev", a_dev, 0_c_intptr_t, &
!                                    a_tmp(1:obj%local_nrows,1:obj%local_ncols), &
!                                    1, 1, num, gpuMemcpyHostToDevice, my_stream, .false., .true., .false.)
!#else
!    successGPU = gpu_memcpy(a_dev, int(loc(a_tmp(1,1)),kind=c_intptr_t), &
!                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
!    call check_memcpy_GPU_f("elpa_cholesky: memcpy a_tmp-> a_dev", 1525,  successGPU)
!#endif
!#endif 

!#if !defined(WITH_NVIDIA_CUSOLVER) && !defined(WITH_AMD_ROCSOLVER)

    deallocate(a_tmp, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("elpa_cholesky: a_tmp", 1536,  istat,  errorMessage)
!#endif   

  ! restore original OpenMP settings
  call obj%timer%stop("cleanup")
  call obj%timer%stop("elpa_cholesky_&
  &complex&
  &_&
  &double&
  &"//gpuString)

    end function elpa_cholesky_d_ptr_complex_double_impl



















!> \brief  elpa_cholesky_a_h_a_complex_single_impl: Cholesky factorization of a single-precision complex hermitian matrix
!> \details
!> \param  obj                    elpa_t object contains:
!> \param     - obj%na            Order of matrix
!> \param     - obj%local_nrows   Leading dimension of a
!> \param     - obj%local_ncols   local columns of matrix a
!> \param     - obj%nblk          blocksize of cyclic distribution, must be the same in both directions!
!> \param     - obj%mpi_comm_rows MPI communicator for rows
!> \param     - obj%mpi_comm_cols MPI communicator for columns
!> \param     - obj%wantDebug     logical, more debug information on failure
!> \param  a(lda,matrixCols)      Distributed matrix which should be inverted
!>                                Distribution is like in Scalapack.
!>                                Only upper triangle needs to be set.
!>                                The lower triangle is not referenced.
!> \result succes                 logical, reports success or failure
    function elpa_cholesky_a_h_a_complex_single_impl(obj, a) result(success)

!    This file is part of ELPA.
!
!    The ELPA library was originally created by the ELPA consortium,
!    consisting of the following organizations:
!
!    - Max Planck Computing and Data Facility (MPCDF), formerly known as
!      Rechenzentrum Garching der Max-Planck-Gesellschaft (RZG),
!    - Bergische Universität Wuppertal, Lehrstuhl für angewandte
!      Informatik,
!    - Technische Universität München, Lehrstuhl für Informatik mit
!      Schwerpunkt Wissenschaftliches Rechnen ,
!    - Fritz-Haber-Institut, Berlin, Abt. Theorie,
!    - Max-Plack-Institut für Mathematik in den Naturwissenschaften,
!      Leipzig, Abt. Komplexe Strukutren in Biologie und Kognition,
!      and
!    - IBM Deutschland GmbH
!
!    This particular source code file contains additions, changes and
!    enhancements authored by Intel Corporation which is not part of
!    the ELPA consortium.
!
!    More information can be found here:
!    http://elpa.mpcdf.mpg.de/
!
!    ELPA is free software: you can redistribute it and/or modify
!    it under the terms of the version 3 of the license of the
!    GNU Lesser General Public License as published by the Free
!    Software Foundation.
!
!    ELPA is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with ELPA.  If not, see <http://www.gnu.org/licenses/>
!
!    ELPA reflects a substantial effort on the part of the original
!    ELPA consortium, and we ask you to respect the spirit of the
!    license that we chose: i.e., please contribute any changes you
!    may have back to the original ELPA library distribution, and keep
!    any derivatives of ELPA under the same license that we chose for
!    the original distribution, the GNU Lesser General Public License.




!cannot use "../src/cholesky/./../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)



  use elpa1_compute
  use elpa_utilities
  use elpa_mpi
  use precision
  use elpa_abstract_impl
  use elpa_omp
  use elpa_blas_interfaces
  use elpa_gpu
  use mod_check_for_gpu
  use invert_trm_gpu !, only : gpu_copy_&
                     !        &single&
                     !        &_tmp1_tmp2, &
                     !        gpu_copy_&
                     !        &single&
                     !        &_a_tmp1
  use cholesky_gpu
  use mod_query_gpu_usage
  implicit none
!    Copyright 2011, A. Marek
!
!    This file is part of ELPA.
!
!    The ELPA library was originally created by the ELPA consortium,
!    consisting of the following organizations:
!
!    - Max Planck Computing and Data Facility (MPCDF), formerly known as
!      Rechenzentrum Garching der Max-Planck-Gesellschaft (RZG),
!    - Bergische Universität Wuppertal, Lehrstuhl für angewandte
!      Informatik,
!    - Technische Universität München, Lehrstuhl für Informatik mit
!      Schwerpunkt Wissenschaftliches Rechnen ,
!    - Fritz-Haber-Institut, Berlin, Abt. Theorie,
!    - Max-Plack-Institut für Mathematik in den Naturwissenschaften,
!      Leipzig, Abt. Komplexe Strukutren in Biologie und Kognition,
!      and
!    - IBM Deutschland GmbH
!
!    This particular source code file contains additions, changes and
!    enhancements authored by Intel Corporation which is not part of
!    the ELPA consortium.
!
!    More information can be found here:
!    http://elpa.mpcdf.mpg.de/
!
!    ELPA is free software: you can redistribute it and/or modify
!    it under the terms of the version 3 of the license of the
!    GNU Lesser General Public License as published by the Free
!    Software Foundation.
!
!    ELPA is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with ELPA.  If not, see <http://www.gnu.org/licenses/>
!
!    ELPA reflects a substantial effort on the part of the original
!    ELPA consortium, and we ask you to respect the spirit of the
!    license that we chose: i.e., please contribute any changes you
!    may have back to the original ELPA library distribution, and keep
!    any derivatives of ELPA under the same license that we chose for
!    the original distribution, the GNU Lesser General Public License.
!
!

  integer, parameter :: rk = C_FLOAT
  integer, parameter :: ck = C_FLOAT_COMPLEX
  integer, parameter :: rck = C_FLOAT_COMPLEX
  complex(kind=rck), parameter     :: ZERO = (0.0_rk,0.0_rk), ONE = (1.0_rk,0.0_rk)
  class(elpa_abstract_impl_t), intent(inout) :: obj
  integer(kind=ik)                           :: na, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols, mpi_comm_all
  complex(kind=rck)                    :: a(obj%local_nrows,*)
  integer(kind=ik)                           :: my_prow, my_pcol, np_rows, np_cols, myid
  integer(kind=MPI_KIND)                     :: mpierr, my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI, myidMPI
  integer(kind=ik)                           :: l_cols, l_rows, l_col1, l_row1, l_colx, l_rowx
  integer(kind=ik)                           :: n, nc, i, info
  integer(kind=BLAS_KIND)                    :: infoBLAS
  integer(kind=ik)                           :: lcs, lce, lrs, lre
  integer(kind=ik)                           :: tile_size, l_rows_tile, l_cols_tile

  complex(kind=rck), allocatable       :: tmp1(:), tmp2(:,:), tmatr(:,:), tmatc(:,:)
  logical                                    :: wantDebug
  logical                                    :: success
  integer(kind=ik)                           :: istat, debug, error
  character(200)                             :: errorMessage
  integer(kind=ik)                           :: nrThreads, limitThreads
  character(20)                              :: gpuString
  logical                                    :: successGPU
  logical                                    :: useGPU
  logical                                    :: useCCL
  integer(kind=c_int)                        :: numGPU
  integer(kind=c_intptr_t)                   :: num
  integer(kind=c_intptr_t)                   :: tmp1_dev, tmatc_dev, tmatr_dev, a_dev, tmp2_dev
  integer(kind=c_intptr_t)                   :: info_dev, info_abs_accumulated_dev
  type(c_ptr)                                :: tmp1_mpi_dev
  complex(kind=rck), pointer           :: tmp1_mpi_fortran_ptr(:,:)
  integer(kind=c_intptr_t)                   :: a_off, tmatc_off, tmatr_off
  type(c_ptr)                                :: tmatc_mpi_dev
  complex(kind=rck), pointer           :: tmatc_mpi_fortran_ptr(:,:)

  integer(kind=c_intptr_t), parameter        :: size_of_datatype = size_of_&
                                                            &single&
                                                            &_&
                                                            &complex

  integer(kind=c_intptr_t)                   :: gpublasHandle, gpusolverHandle, my_stream, offset
  integer(kind=c_int)                        :: gpu_cholesky
  integer(kind=ik)                           :: blocking

  success = .true.
  useGPU = .false.
  useCCL = .false.



  ! check whether the above setting should be overriden
  if (obj%is_set("gpu_cholesky") == 1) then
    call obj%get("gpu_cholesky", gpu_cholesky, error)
    if (error .ne. ELPA_OK) then
      print *,"Problem getting option for gpu_cholesky. Aborting..."
      stop 1
    endif
    if (useGPU .and. gpu_cholesky .eq. 0) then
      useGPU = .false.
    else if (.not.(useGPU) .and. gpu_cholesky .eq. 1) then
      useGPU = .true.
    else 
    endif
  else 
    ! no override by user 
  endif

  if (.not.(useGPU)) then
  endif

  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  
  call obj%timer%start("elpa_cholesky_&
  &complex&
  &_&
  &single&
  &"//gpuString)

  nrThreads=1

  na         = obj%na
  matrixRows = obj%local_nrows
  nblk       = obj%nblk
  matrixCols = obj%local_ncols

  call obj%get("debug",debug,error)
  if (error .ne. ELPA_OK) then
    write(error_unit,*) "ELPA_CHOLESKY: Problem getting option for debug settings. Aborting..."
    success = .false.
    return
  endif
  if (debug == 1) then
    wantDebug = .true.
  else
    wantDebug = .false.
  endif

  mpi_comm_all    = obj%mpi_setup%mpi_comm_parent
  mpi_comm_cols   = obj%mpi_setup%mpi_comm_cols
  mpi_comm_rows   = obj%mpi_setup%mpi_comm_rows
    
  myid    = obj%mpi_setup%myRank_comm_parent
  my_prow = obj%mpi_setup%myRank_comm_rows
  my_pcol = obj%mpi_setup%myRank_comm_cols
        
  np_rows = obj%mpi_setup%nRanks_comm_rows
  np_cols = obj%mpi_setup%nRanks_comm_cols

  success = .true.

  ! Matrix is split into tiles; work is done only for tiles on the diagonal or above

  call obj%timer%start("prepare")

  if (obj%is_set("blocking_in_cholesky") == 1) then
    call obj%get("blocking_in_cholesky", blocking, error)
    if (error .ne. ELPA_OK) then
      write(error_unit,*) "ELPA_CHOLESKY: Problem in getting keyword 'blocking_in_cholesky'. Aborting..."
      stop 1
    endif
  else
    if (useGPU) then
      blocking = 1024
    else
      blocking = 128
    endif
    call obj%set("blocking_in_cholesky", blocking, error)
    if (error .ne. ELPA_OK) then
      write(error_unit,*) "ELPA_CHOLESKY: Problem in setting keyword 'blocking_in_cholesky'. Aborting..."
      stop 1
    endif
  endif 


  tile_size = nblk*least_common_multiple(np_rows,np_cols) ! minimum global tile size
  tile_size = ((blocking*max(np_rows,np_cols)-1)/tile_size+1)*tile_size ! make local tiles at least 128 wide

  l_rows_tile = tile_size/np_rows ! local rows of a tile
  l_cols_tile = tile_size/np_cols ! local cols of a tile

  l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a
  l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a

  if (useGPU) then
    call obj%timer%start("check_for_gpu")
    if (check_for_gpu(obj, myid, numGPU)) then
      ! set the neccessary parameters
      call set_gpu_parameters()
    else
      print *,"ELPA_CHOLESKY: GPUs are requested but not detected! Aborting..."
      success = .false.
      return
    endif
    call obj%timer%stop("check_for_gpu")
  else ! useGPU
  endif ! useGPU


  if (useGPU) then

    successGPU = gpu_malloc(tmp1_dev, nblk*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmp1_dev", 335,  successGPU)

    successGPU = gpu_memset(tmp1_dev, 0, nblk*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmp1_dev", 346,  successGPU)

    successGPU = gpu_malloc(tmp2_dev, nblk*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmp2_dev", 350,  successGPU)

    successGPU = gpu_memset(tmp2_dev, 0, nblk*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmp2_dev", 361,  successGPU)

    successGPU = gpu_malloc(tmatc_dev, l_cols*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmatc_dev", 365,  successGPU)

    successGPU = gpu_memset(tmatc_dev, 0, l_cols*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmatc_dev", 376,  successGPU)

    successGPU = gpu_malloc(tmatr_dev, l_rows*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmatr_dev", 380,  successGPU)

    successGPU = gpu_memset(tmatr_dev, 0, l_rows*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmatr_dev", 391,  successGPU)

    successGPU = gpu_malloc(a_dev, matrixRows*matrixCols*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: a_dev", 396,  successGPU)
  endif ! useGPU


  allocate(tmp1(nblk*nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmp1", 441,  istat,  errorMessage)

  allocate(tmp2(nblk,nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmp2", 444,  istat,  errorMessage)



  tmp1 = 0
  tmp2 = 0

  allocate(tmatr(l_rows,nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmatr", 460,  istat,  errorMessage)

  allocate(tmatc(l_cols,nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmatc", 463,  istat,  errorMessage)


  tmatr = 0
  tmatc = 0

  if (useGPU) then
    num = matrixRows*matrixCols* size_of_datatype
    successGPU = gpu_memcpy(a_dev, int(loc(a(1,1)),kind=c_intptr_t), &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
    call check_memcpy_GPU_f("elpa_cholesky 1: memcpy a-> a_dev", 495,  successGPU)
  endif


  call obj%timer%stop("prepare")
  call obj%timer%start("loop1")
  do n = 1, na, nblk


    ! Calculate first local row and column of the still remaining matrix
    ! on the local processor

    l_row1 = local_index(n, my_prow, np_rows, nblk, +1)
    l_col1 = local_index(n, my_pcol, np_cols, nblk, +1)

    l_rowx = local_index(n+nblk, my_prow, np_rows, nblk, +1)
    l_colx = local_index(n+nblk, my_pcol, np_cols, nblk, +1)

    if (n+nblk > na) then

      ! This is the last step, just do a Cholesky-Factorization
      ! of the remaining block

      if (useGPU) then
        if (my_prow==prow(n, nblk, np_rows) .and. my_pcol==pcol(n, nblk, np_cols)) then

          call obj%timer%start("lapack")
          num = matrixRows*matrixCols* size_of_datatype
          successGPU = gpu_memcpy(int(loc(a(1,1)),kind=c_intptr_t), a_dev,  &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a", 582,  successGPU)

          call CPOTRF('U', int(na-n+1,kind=BLAS_KIND), a(l_row1,l_col1), &
                             int(matrixRows,kind=BLAS_KIND), infoBLAS )
          info = int(infoBLAS,kind=ik)

          num = matrixRows*matrixCols* size_of_datatype
          successGPU = gpu_memcpy(a_dev, int(loc(a(1,1)),kind=c_intptr_t), &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a", 600,  successGPU)

          call obj%timer%stop("lapack")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &complex&

            &: Error in zpotrf: ",info
            success = .false. ! "
            return
          endif ! info
        endif ! (my_prow==prow(n, nblk, np_rows) .and. my_pcol==pcol(n, nblk, np_cols))

      else ! useGPU
        if (my_prow==prow(n, nblk, np_rows) .and. my_pcol==pcol(n, nblk, np_cols)) then
          call obj%timer%start("blas")

          call CPOTRF('U', int(na-n+1,kind=BLAS_KIND), a(l_row1,l_col1), &
                             int(matrixRows,kind=BLAS_KIND), infoBLAS )
          info = int(infoBLAS,kind=ik)
          call obj%timer%stop("blas")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &complex&

            &: Error in zpotrf: ",info
            success = .false. ! "
            return
          endif

        endif
      endif ! useGPU
      
      exit ! Loop
    endif ! (n+nblk > na) 

    ! This is not the last step

    if (my_prow==prow(n, nblk, np_rows)) then

      if (my_pcol==pcol(n, nblk, np_cols)) then

        if (useGPU) then
          call obj%timer%start("lapack")

          num = matrixRows*matrixCols* size_of_datatype

          successGPU = gpu_memcpy(int(loc(a(1,1)),kind=c_intptr_t), a_dev,  &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a", 742,  successGPU)


          call CPOTRF('U', int(nblk,kind=BLAS_KIND), a(l_row1,l_col1), &
                               int(matrixRows,kind=BLAS_KIND) , infoBLAS )
          info = int(infoBLAS,kind=ik)

          num = matrixRows*matrixCols* size_of_datatype
          successGPU = gpu_memcpy(a_dev, int(loc(a(1,1)),kind=c_intptr_t), &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a", 769,  successGPU)

          call obj%timer%stop("lapack")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &complex&

            &: Error in zpotrf: ",info
            success = .false. ! "
            return
          endif ! info
        else ! useGPU
          ! The process owning the upper left remaining block does the
          ! Cholesky-Factorization of this block
          call obj%timer%start("blas")

          call CPOTRF('U', int(nblk,kind=BLAS_KIND), a(l_row1,l_col1), &
                               int(matrixRows,kind=BLAS_KIND) , infoBLAS )
          info = int(infoBLAS,kind=ik)
          call obj%timer%stop("blas")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &complex&

            &: Error in zpotrf 2: ",info

            success = .false. ! "
            return
          endif ! info/=0
        endif ! useGPU
 
        if (useGPU) then
          my_stream = obj%gpu_setup%my_stream
          call gpu_copy_float_complex_a_tmp1 (a_dev, tmp1_dev, l_row1, l_col1, matrixRows, nblk, my_stream)
        else ! useGPU
          nc = 0
          do i=1,nblk
            tmp1(nc+1:nc+i) = a(l_row1:l_row1+i-1,l_col1+i-1)
            nc = nc+i
          enddo
        endif ! useGPU
      endif ! (my_pcol==pcol(n, nblk, np_cols))


      if (useGPU) then
        my_stream = obj%gpu_setup%my_stream
        call gpu_copy_float_complex_tmp1_tmp2 (tmp1_dev, tmp2_dev, nblk, nblk, my_stream)
      else ! useGPU
        nc = 0
        do i=1,nblk
          tmp2(1:i,i) = tmp1(nc+1:nc+i)
          nc = nc+i
        enddo
      endif ! useGPU

      if (useGPU) then
        call obj%timer%start("gpublas")
        gpublasHandle = obj%gpu_setup%gpublasHandleArray(0)
        if (l_cols-l_colx+1 > 0) then
          a_off = (l_row1-1 + (l_colx-1)*matrixRows) * size_of_datatype
          call gpublas_CTRSM('L', 'U', 'C', 'N', nblk, l_cols-l_colx+1, ONE, &
                            tmp2_dev, nblk, a_dev+a_off, matrixRows, gpublasHandle)
          if (wantDebug .and. gpu_vendor() /= SYCL_GPU) successGPU = gpu_DeviceSynchronize()
        endif
        call obj%timer%stop("gpublas")

      else ! useGPU

        call obj%timer%start("lapack")
        if (l_cols-l_colx+1>0) then
          call CTRSM('L', 'U', 'C', 'N', int(nblk,kind=BLAS_KIND),  &
                            int(l_cols-l_colx+1,kind=BLAS_KIND), ONE, tmp2, &
                            int(ubound(tmp2,dim=1),kind=BLAS_KIND), a(l_row1,l_colx), int(matrixRows,kind=BLAS_KIND) )
        endif
        call obj%timer%stop("lapack")
      endif ! useGPU
    endif ! (my_prow==prow(n, nblk, np_rows))


    if (useGPU) then
      if (my_prow==prow(n, nblk, np_rows)) then
        ! if l_cols-l_colx+1 == 0 kernel launch with 0 blocks => raises error
        call obj%timer%start("a_tmatc_kernel")
        if (l_cols-l_colx+1>0) then
           my_stream = obj%gpu_setup%my_stream
           call gpu_copy_float_complex_a_tmatc(a_dev, tmatc_dev, nblk, matrixRows, l_cols, l_colx, l_row1, my_stream)
        endif
        call obj%timer%stop("a_tmatc_kernel")
      endif
    else ! useGPU
      do i=1,nblk
        if (my_prow==prow(n, nblk, np_rows)) tmatc(l_colx:l_cols,i) = conjg(a(l_row1+i-1,l_colx:l_cols))
      enddo
    endif ! useGPU

    if (useGPU  .and. .not. useCCL) then
      num = l_cols*nblk*size_of_datatype
      successGPU = gpu_memcpy(int(loc(tmatc),kind=c_intptr_t), tmatc_dev, num, &
                              gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("elpa_cholesky: tmatc_dev to tmatc", 1201,  successGPU)
    endif


    if (useGPU .and. .not. useCCL) then

      ! not needed if transpose vec CCL
      num = l_rows*nblk*size_of_datatype
      successGPU = gpu_memcpy(int(loc(tmatr),kind=c_intptr_t), tmatr_dev, num, &
                              gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("elpa_cholesky: tmatr_dev to tmatr", 1220,  successGPU)
    endif ! (useGPU .and. .not. useCCL)


    if (useCCL) then

    else ! useCCL
      call elpa_transpose_vectors_&
      &complex&
      &_&
      &single &
      (obj, tmatc, ubound(tmatc,dim=1), mpi_comm_cols, &
      tmatr, ubound(tmatr,dim=1), mpi_comm_rows, &
      n, na, nblk, nblk, nrThreads, .false., success)
    endif ! useCCL

    if (useGPU .and. .not. useCCL) then

      num = l_rows*nblk*size_of_datatype
      successGPU = gpu_memcpy(tmatr_dev, int(loc(tmatr),kind=c_intptr_t), num, &
                              gpuMemcpyHostToDevice)
      call check_memcpy_GPU_f("elpa_cholesky: tmat to tmatr_dev", 1271,  successGPU)
    endif

    if (useGPU) then
      do i=0,(na-1)/tile_size
        lcs = max(l_colx,i*l_cols_tile+1)
        lce = min(l_cols,(i+1)*l_cols_tile)
        lrs = l_rowx
        lre = min(l_rows,(i+1)*l_rows_tile)
        if (lce  < lcs .or. lre < lrs) cycle
        call obj%timer%start("gpublas")
        gpublasHandle = obj%gpu_setup%gpublasHandleArray(0)
        tmatr_off = (lrs-1 + (1-1)*l_rows) * size_of_datatype
        tmatc_off = (lcs-1 + (1-1)*l_cols) * size_of_datatype
        a_off = (lrs-1 + (lcs-1)*matrixRows) * size_of_datatype
        call gpublas_CGEMM('N', 'C', lre-lrs+1, lce-lcs+1, nblk, &
                            -ONE, tmatr_dev+tmatr_off, l_rows, tmatc_dev+tmatc_off, l_cols, ONE, &
                            a_dev+a_off, matrixRows, gpublasHandle)
        if (wantDebug .and. gpu_vendor() /= SYCL_GPU) successGPU = gpu_DeviceSynchronize()
        call obj%timer%stop("gpublas")
      enddo
    else !useGPU
      do i=0,(na-1)/tile_size
        lcs = max(l_colx,i*l_cols_tile+1)
        lce = min(l_cols,(i+1)*l_cols_tile)
        lrs = l_rowx
        lre = min(l_rows,(i+1)*l_rows_tile)
        if (lce<lcs .or. lre<lrs) cycle
        call obj%timer%start("blas")
        call CGEMM('N', 'C', int(lre-lrs+1,kind=BLAS_KIND), int(lce-lcs+1,kind=BLAS_KIND), &
                            int(nblk,kind=BLAS_KIND), -ONE,  &
                            tmatr(lrs,1), int(ubound(tmatr,dim=1),kind=BLAS_KIND), tmatc(lcs,1), &
                            int(ubound(tmatc,dim=1),kind=BLAS_KIND), &
                            ONE, a(lrs,lcs), int(matrixRows,kind=BLAS_KIND))
        call obj%timer%stop("blas")
      enddo
    endif ! useGPU


  enddo ! n = 1, na, nblk

  call obj%timer%stop("loop1")
  
  call obj%timer%start("deallocate")
  if (useGPU) then



    successGPU = gpu_free(tmp1_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmp1_dev", 1351,  successGPU)

    successGPU = gpu_free(tmp2_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmp1_dev", 1354,  successGPU)

    successGPU = gpu_free(tmatc_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmatc_dev", 1357,  successGPU)

    successGPU = gpu_free(tmatr_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmatr_dev", 1360,  successGPU)

  endif


  deallocate(tmp1, tmp2, tmatr, tmatc, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("elpa_cholesky: tmp1, tmp2, tmatr, tmatc", 1384,  istat,  errorMessage)

  call obj%timer%stop("deallocate")


  call obj%timer%start("loop2")

  call obj%timer%start("copy1")

  if (useGPU) then
    num = matrixRows*matrixCols* size_of_datatype
    successGPU = gpu_memcpy(int(loc(a(1,1)),kind=c_intptr_t), a_dev,  &
                     matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
    call check_memcpy_GPU_f("elpa_cholesky: memcpy 2 a-> a_dev", 1405,  successGPU)
  endif
  ! Set the lower triangle to 0, it contains garbage (form the above matrix multiplications)

  call obj%timer%stop("copy1")
  do i=1,na
    if (my_pcol==pcol(i, nblk, np_cols)) then
      ! column i is on local processor
      l_col1 = local_index(i  , my_pcol, np_cols, nblk, +1) ! local column number
      l_row1 = local_index(i+1, my_prow, np_rows, nblk, +1) ! first row below diagonal
      if (useGPU) then
        call obj%timer%start("copy2")
!#ifndef DEVICE_POINTER
        ! set to zero maybe a kerne would be better here
          offset = (l_row1-1 + matrixRows * (l_col1-1)) * size_of_datatype
          num = (l_rows-l_row1+1)*size_of_datatype
          successGPU = gpu_memset(a_dev+offset, 0, num)
          call check_memcpy_GPU_f("elpa_cholesky: memset a_dev", 1450,  successGPU)
!#else 
!          !offset = (l_row1-1 + matrixRows * (l_col1-1)) * size_of_datatype
!          !num = (l_rows-l_row1+1)*size_of_datatype
!#ifdef WITH_GPU_STREAMS
!          !my_stream = obj%gpu_setup%my_stream
!          !successGPU = gpu_memset_async(a_dev+offset, 0, num, my_stream)
!          !call check_memcpy_GPU_f("elpa_cholesky: memset a_dev", 1458,  successGPU)
!
!          !successGPU = gpu_stream_synchronize(my_stream)
!          !call check_memcpy_GPU_f("elpa_cholesky: memset", 1461,  successGPU)
!#else
!          !successGPU = gpu_memset(a_dev+offset, 0, num)
!          !call check_memcpy_GPU_f("elpa_cholesky: memset tmp1_dev", 1464,  successGPU)
!#endif
!
!
!!#if defined(WITH_NVIDIA_CUSOLVER) || defined(WITH_AMD_ROCSOLVER)
!!        successGPU = gpu_memcpy(int(loc(a_tmp(1,1)),kind=c_intptr_t), a_dev, &
!!                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
!!        call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a_tmp", 1471,  successGPU)
!!#endif
!!        a_tmp(l_row1:l_rows,l_col1) = 0
!!#if defined(WITH_NVIDIA_CUSOLVER) || defined(WITH_AMD_ROCSOLVER)
!!        successGPU = gpu_memcpy(a_dev, int(loc(a_tmp(1,1)),kind=c_intptr_t), &
!!                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
!!        call check_memcpy_GPU_f("elpa_cholesky: memcpy a_tmp-> a_dev", 1477,  successGPU)
!#endif

!#endif 

  call obj%timer%stop("copy2")

      else ! useGPU
        a(l_row1:l_rows,l_col1) = 0
      endif ! useGPU
    endif ! (my_pcol==pcol(i, nblk, np_cols))
  enddo
  !endif ! useGPU

  call obj%timer%stop("loop2")


  call obj%timer%start("cleanup")
  if (useGPU) then
    ! copy back
    successGPU = gpu_memcpy(int(loc(a(1,1)),kind=c_intptr_t), a_dev,  &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
    call check_memcpy_GPU_f("elpa_cholesky: memcpy a-> d_dev", 1502,  successGPU)

    successGPU = gpu_free(a_dev)
    call check_dealloc_GPU_f("elpa_cholesky: a_dev", 1505,  successGPU)


  endif

  ! restore original OpenMP settings
  call obj%timer%stop("cleanup")
  call obj%timer%stop("elpa_cholesky_&
  &complex&
  &_&
  &single&
  &"//gpuString)

    end function elpa_cholesky_a_h_a_complex_single_impl




















!> \brief  elpa_cholesky_d_ptr_complex_single_impl: Cholesky factorization of a single-precision complex hermitian matrix
!> \details
!> \param  obj                    elpa_t object contains:
!> \param     - obj%na            Order of matrix
!> \param     - obj%local_nrows   Leading dimension of a
!> \param     - obj%local_ncols   local columns of matrix a
!> \param     - obj%nblk          blocksize of cyclic distribution, must be the same in both directions!
!> \param     - obj%mpi_comm_rows MPI communicator for rows
!> \param     - obj%mpi_comm_cols MPI communicator for columns
!> \param     - obj%wantDebug     logical, more debug information on failure
!> \param  aDev(lda,matrixCols)   Distributed matrix which should be inverted, as type(c_ptr) living on device
!>                                Distribution is like in Scalapack.
!>                                Only upper triangle needs to be set.
!>                                The lower triangle is not referenced.
!> \result succes                 logical, reports success or failure
    function elpa_cholesky_d_ptr_complex_single_impl(obj, aDev) result(success)

!    This file is part of ELPA.
!
!    The ELPA library was originally created by the ELPA consortium,
!    consisting of the following organizations:
!
!    - Max Planck Computing and Data Facility (MPCDF), formerly known as
!      Rechenzentrum Garching der Max-Planck-Gesellschaft (RZG),
!    - Bergische Universität Wuppertal, Lehrstuhl für angewandte
!      Informatik,
!    - Technische Universität München, Lehrstuhl für Informatik mit
!      Schwerpunkt Wissenschaftliches Rechnen ,
!    - Fritz-Haber-Institut, Berlin, Abt. Theorie,
!    - Max-Plack-Institut für Mathematik in den Naturwissenschaften,
!      Leipzig, Abt. Komplexe Strukutren in Biologie und Kognition,
!      and
!    - IBM Deutschland GmbH
!
!    This particular source code file contains additions, changes and
!    enhancements authored by Intel Corporation which is not part of
!    the ELPA consortium.
!
!    More information can be found here:
!    http://elpa.mpcdf.mpg.de/
!
!    ELPA is free software: you can redistribute it and/or modify
!    it under the terms of the version 3 of the license of the
!    GNU Lesser General Public License as published by the Free
!    Software Foundation.
!
!    ELPA is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with ELPA.  If not, see <http://www.gnu.org/licenses/>
!
!    ELPA reflects a substantial effort on the part of the original
!    ELPA consortium, and we ask you to respect the spirit of the
!    license that we chose: i.e., please contribute any changes you
!    may have back to the original ELPA library distribution, and keep
!    any derivatives of ELPA under the same license that we chose for
!    the original distribution, the GNU Lesser General Public License.




!cannot use "../src/cholesky/./../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)



  use elpa1_compute
  use elpa_utilities
  use elpa_mpi
  use precision
  use elpa_abstract_impl
  use elpa_omp
  use elpa_blas_interfaces
  use elpa_gpu
  use mod_check_for_gpu
  use invert_trm_gpu !, only : gpu_copy_&
                     !        &single&
                     !        &_tmp1_tmp2, &
                     !        gpu_copy_&
                     !        &single&
                     !        &_a_tmp1
  use cholesky_gpu
  use mod_query_gpu_usage
  implicit none
!    Copyright 2011, A. Marek
!
!    This file is part of ELPA.
!
!    The ELPA library was originally created by the ELPA consortium,
!    consisting of the following organizations:
!
!    - Max Planck Computing and Data Facility (MPCDF), formerly known as
!      Rechenzentrum Garching der Max-Planck-Gesellschaft (RZG),
!    - Bergische Universität Wuppertal, Lehrstuhl für angewandte
!      Informatik,
!    - Technische Universität München, Lehrstuhl für Informatik mit
!      Schwerpunkt Wissenschaftliches Rechnen ,
!    - Fritz-Haber-Institut, Berlin, Abt. Theorie,
!    - Max-Plack-Institut für Mathematik in den Naturwissenschaften,
!      Leipzig, Abt. Komplexe Strukutren in Biologie und Kognition,
!      and
!    - IBM Deutschland GmbH
!
!    This particular source code file contains additions, changes and
!    enhancements authored by Intel Corporation which is not part of
!    the ELPA consortium.
!
!    More information can be found here:
!    http://elpa.mpcdf.mpg.de/
!
!    ELPA is free software: you can redistribute it and/or modify
!    it under the terms of the version 3 of the license of the
!    GNU Lesser General Public License as published by the Free
!    Software Foundation.
!
!    ELPA is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with ELPA.  If not, see <http://www.gnu.org/licenses/>
!
!    ELPA reflects a substantial effort on the part of the original
!    ELPA consortium, and we ask you to respect the spirit of the
!    license that we chose: i.e., please contribute any changes you
!    may have back to the original ELPA library distribution, and keep
!    any derivatives of ELPA under the same license that we chose for
!    the original distribution, the GNU Lesser General Public License.
!
!

  integer, parameter :: rk = C_FLOAT
  integer, parameter :: ck = C_FLOAT_COMPLEX
  integer, parameter :: rck = C_FLOAT_COMPLEX
  complex(kind=rck), parameter     :: ZERO = (0.0_rk,0.0_rk), ONE = (1.0_rk,0.0_rk)
  class(elpa_abstract_impl_t), intent(inout) :: obj
  integer(kind=ik)                           :: na, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols, mpi_comm_all
  type(c_ptr)                                :: aDev
!#if !defined(WITH_NVIDIA_CUSOLVER) && !defined(WITH_AMD_ROCSOLVER)
  complex(kind=rck), allocatable       :: a_tmp(:,:)
!#endif
  integer(kind=ik)                           :: my_prow, my_pcol, np_rows, np_cols, myid
  integer(kind=MPI_KIND)                     :: mpierr, my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI, myidMPI
  integer(kind=ik)                           :: l_cols, l_rows, l_col1, l_row1, l_colx, l_rowx
  integer(kind=ik)                           :: n, nc, i, info
  integer(kind=BLAS_KIND)                    :: infoBLAS
  integer(kind=ik)                           :: lcs, lce, lrs, lre
  integer(kind=ik)                           :: tile_size, l_rows_tile, l_cols_tile

  complex(kind=rck), allocatable       :: tmp1(:), tmp2(:,:), tmatr(:,:), tmatc(:,:)
  logical                                    :: wantDebug
  logical                                    :: success
  integer(kind=ik)                           :: istat, debug, error
  character(200)                             :: errorMessage
  integer(kind=ik)                           :: nrThreads, limitThreads
  character(20)                              :: gpuString
  logical                                    :: successGPU
  logical                                    :: useGPU
  logical                                    :: useCCL
  integer(kind=c_int)                        :: numGPU
  integer(kind=c_intptr_t)                   :: num
  integer(kind=c_intptr_t)                   :: tmp1_dev, tmatc_dev, tmatr_dev, a_dev, tmp2_dev
  integer(kind=c_intptr_t)                   :: info_dev, info_abs_accumulated_dev
  type(c_ptr)                                :: tmp1_mpi_dev
  complex(kind=rck), pointer           :: tmp1_mpi_fortran_ptr(:,:)
  integer(kind=c_intptr_t)                   :: a_off, tmatc_off, tmatr_off
  type(c_ptr)                                :: tmatc_mpi_dev
  complex(kind=rck), pointer           :: tmatc_mpi_fortran_ptr(:,:)

  integer(kind=c_intptr_t), parameter        :: size_of_datatype = size_of_&
                                                            &single&
                                                            &_&
                                                            &complex

  integer(kind=c_intptr_t)                   :: gpublasHandle, gpusolverHandle, my_stream, offset
  integer(kind=c_int)                        :: gpu_cholesky
  integer(kind=ik)                           :: blocking

  success = .true.
  useGPU = .false.
  useCCL = .false.

  useGPU = .true.

  if (.not.(useGPU)) then
    print *,"You used the interface for device pointers for elpa_cholesky but did not specify GPU usage!. Aborting..."
    stop 1
  endif

  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  
  call obj%timer%start("elpa_cholesky_&
  &complex&
  &_&
  &single&
  &"//gpuString)

  nrThreads=1

  na         = obj%na
  matrixRows = obj%local_nrows
  nblk       = obj%nblk
  matrixCols = obj%local_ncols

  call obj%get("debug",debug,error)
  if (error .ne. ELPA_OK) then
    write(error_unit,*) "ELPA_CHOLESKY: Problem getting option for debug settings. Aborting..."
    success = .false.
    return
  endif
  if (debug == 1) then
    wantDebug = .true.
  else
    wantDebug = .false.
  endif

  mpi_comm_all    = obj%mpi_setup%mpi_comm_parent
  mpi_comm_cols   = obj%mpi_setup%mpi_comm_cols
  mpi_comm_rows   = obj%mpi_setup%mpi_comm_rows
    
  myid    = obj%mpi_setup%myRank_comm_parent
  my_prow = obj%mpi_setup%myRank_comm_rows
  my_pcol = obj%mpi_setup%myRank_comm_cols
        
  np_rows = obj%mpi_setup%nRanks_comm_rows
  np_cols = obj%mpi_setup%nRanks_comm_cols

  success = .true.

  ! Matrix is split into tiles; work is done only for tiles on the diagonal or above

  call obj%timer%start("prepare")

  if (obj%is_set("blocking_in_cholesky") == 1) then
    call obj%get("blocking_in_cholesky", blocking, error)
    if (error .ne. ELPA_OK) then
      write(error_unit,*) "ELPA_CHOLESKY: Problem in getting keyword 'blocking_in_cholesky'. Aborting..."
      stop 1
    endif
  else
    if (useGPU) then
      blocking = 1024
    else
      blocking = 128
    endif
    call obj%set("blocking_in_cholesky", blocking, error)
    if (error .ne. ELPA_OK) then
      write(error_unit,*) "ELPA_CHOLESKY: Problem in setting keyword 'blocking_in_cholesky'. Aborting..."
      stop 1
    endif
  endif 


  tile_size = nblk*least_common_multiple(np_rows,np_cols) ! minimum global tile size
  tile_size = ((blocking*max(np_rows,np_cols)-1)/tile_size+1)*tile_size ! make local tiles at least 128 wide

  l_rows_tile = tile_size/np_rows ! local rows of a tile
  l_cols_tile = tile_size/np_cols ! local cols of a tile

  l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a
  l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a

  if (useGPU) then
    call obj%timer%start("check_for_gpu")
    if (check_for_gpu(obj, myid, numGPU)) then
      ! set the neccessary parameters
      call set_gpu_parameters()
    else
      print *,"ELPA_CHOLESKY: GPUs are requested but not detected! Aborting..."
      success = .false.
      return
    endif
    call obj%timer%stop("check_for_gpu")
  else ! useGPU
  endif ! useGPU


  if (useGPU) then

    successGPU = gpu_malloc(tmp1_dev, nblk*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmp1_dev", 335,  successGPU)

    successGPU = gpu_memset(tmp1_dev, 0, nblk*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmp1_dev", 346,  successGPU)

    successGPU = gpu_malloc(tmp2_dev, nblk*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmp2_dev", 350,  successGPU)

    successGPU = gpu_memset(tmp2_dev, 0, nblk*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmp2_dev", 361,  successGPU)

    successGPU = gpu_malloc(tmatc_dev, l_cols*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmatc_dev", 365,  successGPU)

    successGPU = gpu_memset(tmatc_dev, 0, l_cols*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmatc_dev", 376,  successGPU)

    successGPU = gpu_malloc(tmatr_dev, l_rows*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmatr_dev", 380,  successGPU)

    successGPU = gpu_memset(tmatr_dev, 0, l_rows*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmatr_dev", 391,  successGPU)

    a_dev = transfer(aDev, a_dev)
!#if !defined(WITH_NVIDIA_CUSOLVER) && !defined(WITH_AMD_ROCSOLVER)
    allocate(a_tmp(obj%local_nrows,obj%local_ncols), stat=istat, errmsg=errorMessage)
    call check_allocate_f("elpa_cholesky: a_tmp", 407,  istat,  errorMessage)
!#endif 
  endif ! useGPU


  allocate(tmp1(nblk*nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmp1", 441,  istat,  errorMessage)

  allocate(tmp2(nblk,nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmp2", 444,  istat,  errorMessage)



  tmp1 = 0
  tmp2 = 0

  allocate(tmatr(l_rows,nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmatr", 460,  istat,  errorMessage)

  allocate(tmatc(l_cols,nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmatc", 463,  istat,  errorMessage)


  tmatr = 0
  tmatc = 0



  call obj%timer%stop("prepare")
  call obj%timer%start("loop1")
  do n = 1, na, nblk


    ! Calculate first local row and column of the still remaining matrix
    ! on the local processor

    l_row1 = local_index(n, my_prow, np_rows, nblk, +1)
    l_col1 = local_index(n, my_pcol, np_cols, nblk, +1)

    l_rowx = local_index(n+nblk, my_prow, np_rows, nblk, +1)
    l_colx = local_index(n+nblk, my_pcol, np_cols, nblk, +1)

    if (n+nblk > na) then

      ! This is the last step, just do a Cholesky-Factorization
      ! of the remaining block

      if (useGPU) then
        if (my_prow==prow(n, nblk, np_rows) .and. my_pcol==pcol(n, nblk, np_cols)) then

          call obj%timer%start("lapack")
          num = matrixRows*matrixCols* size_of_datatype

          successGPU = gpu_memcpy(int(loc(a_tmp(1,1)),kind=c_intptr_t), a_dev,  &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a_tmp", 615,  successGPU)

          call CPOTRF('U', int(na-n+1,kind=BLAS_KIND), a_tmp(l_row1,l_col1), &
                             int(matrixRows,kind=BLAS_KIND), infoBLAS )
          info = int(infoBLAS,kind=ik)

          num = matrixRows*matrixCols* size_of_datatype
          successGPU = gpu_memcpy(a_dev, int(loc(a_tmp(1,1)),kind=c_intptr_t), &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a_tmp", 633,  successGPU)

          call obj%timer%stop("lapack")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &complex&

            &: Error in zpotrf: ",info
            success = .false. ! "
            return
          endif ! info
        endif ! (my_prow==prow(n, nblk, np_rows) .and. my_pcol==pcol(n, nblk, np_cols))

      else ! useGPU
        if (my_prow==prow(n, nblk, np_rows) .and. my_pcol==pcol(n, nblk, np_cols)) then
          call obj%timer%start("blas")

          info = int(infoBLAS,kind=ik)
          call obj%timer%stop("blas")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &complex&

            &: Error in zpotrf: ",info
            success = .false. ! "
            return
          endif

        endif
      endif ! useGPU
      
      exit ! Loop
    endif ! (n+nblk > na) 

    ! This is not the last step

    if (my_prow==prow(n, nblk, np_rows)) then

      if (my_pcol==pcol(n, nblk, np_cols)) then

        if (useGPU) then
          call obj%timer%start("lapack")

          num = matrixRows*matrixCols* size_of_datatype

          successGPU = gpu_memcpy(int(loc(a_tmp(1,1)),kind=c_intptr_t), a_dev,  &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a_tmp", 783,  successGPU)

          call CPOTRF('U', int(nblk,kind=BLAS_KIND), a_tmp(l_row1,l_col1), &
                               int(matrixRows,kind=BLAS_KIND) , infoBLAS )
          info = int(infoBLAS,kind=ik)
          num = matrixRows*matrixCols* size_of_datatype
          successGPU = gpu_memcpy(a_dev, int(loc(a_tmp(1,1)),kind=c_intptr_t), &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a_tmp", 799,  successGPU)

          call obj%timer%stop("lapack")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &complex&

            &: Error in zpotrf: ",info
            success = .false. ! "
            return
          endif ! info
        else ! useGPU
          ! The process owning the upper left remaining block does the
          ! Cholesky-Factorization of this block
          call obj%timer%start("blas")

          info = int(infoBLAS,kind=ik)
          call obj%timer%stop("blas")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &complex&

            &: Error in zpotrf 2: ",info

            success = .false. ! "
            return
          endif ! info/=0
        endif ! useGPU
 
        if (useGPU) then
          my_stream = obj%gpu_setup%my_stream
          call gpu_copy_float_complex_a_tmp1 (a_dev, tmp1_dev, l_row1, l_col1, matrixRows, nblk, my_stream)
        else ! useGPU
          nc = 0
          do i=1,nblk
            nc = nc+i
          enddo
        endif ! useGPU
      endif ! (my_pcol==pcol(n, nblk, np_cols))


      if (useGPU) then
        my_stream = obj%gpu_setup%my_stream
        call gpu_copy_float_complex_tmp1_tmp2 (tmp1_dev, tmp2_dev, nblk, nblk, my_stream)
      else ! useGPU
        nc = 0
        do i=1,nblk
          tmp2(1:i,i) = tmp1(nc+1:nc+i)
          nc = nc+i
        enddo
      endif ! useGPU

      if (useGPU) then
        call obj%timer%start("gpublas")
        gpublasHandle = obj%gpu_setup%gpublasHandleArray(0)
        if (l_cols-l_colx+1 > 0) then
          a_off = (l_row1-1 + (l_colx-1)*matrixRows) * size_of_datatype
          call gpublas_CTRSM('L', 'U', 'C', 'N', nblk, l_cols-l_colx+1, ONE, &
                            tmp2_dev, nblk, a_dev+a_off, matrixRows, gpublasHandle)
          if (wantDebug .and. gpu_vendor() /= SYCL_GPU) successGPU = gpu_DeviceSynchronize()
        endif
        call obj%timer%stop("gpublas")

      else ! useGPU

        call obj%timer%start("lapack")
        call obj%timer%stop("lapack")
      endif ! useGPU
    endif ! (my_prow==prow(n, nblk, np_rows))


    if (useGPU) then
      if (my_prow==prow(n, nblk, np_rows)) then
        ! if l_cols-l_colx+1 == 0 kernel launch with 0 blocks => raises error
        call obj%timer%start("a_tmatc_kernel")
        if (l_cols-l_colx+1>0) then
           my_stream = obj%gpu_setup%my_stream
           call gpu_copy_float_complex_a_tmatc(a_dev, tmatc_dev, nblk, matrixRows, l_cols, l_colx, l_row1, my_stream)
        endif
        call obj%timer%stop("a_tmatc_kernel")
      endif
    else ! useGPU
      do i=1,nblk
      enddo
    endif ! useGPU

    if (useGPU  .and. .not. useCCL) then
      num = l_cols*nblk*size_of_datatype
      successGPU = gpu_memcpy(int(loc(tmatc),kind=c_intptr_t), tmatc_dev, num, &
                              gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("elpa_cholesky: tmatc_dev to tmatc", 1201,  successGPU)
    endif


    if (useGPU .and. .not. useCCL) then

      ! not needed if transpose vec CCL
      num = l_rows*nblk*size_of_datatype
      successGPU = gpu_memcpy(int(loc(tmatr),kind=c_intptr_t), tmatr_dev, num, &
                              gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("elpa_cholesky: tmatr_dev to tmatr", 1220,  successGPU)
    endif ! (useGPU .and. .not. useCCL)


    if (useCCL) then

    else ! useCCL
      call elpa_transpose_vectors_&
      &complex&
      &_&
      &single &
      (obj, tmatc, ubound(tmatc,dim=1), mpi_comm_cols, &
      tmatr, ubound(tmatr,dim=1), mpi_comm_rows, &
      n, na, nblk, nblk, nrThreads, .false., success)
    endif ! useCCL

    if (useGPU .and. .not. useCCL) then

      num = l_rows*nblk*size_of_datatype
      successGPU = gpu_memcpy(tmatr_dev, int(loc(tmatr),kind=c_intptr_t), num, &
                              gpuMemcpyHostToDevice)
      call check_memcpy_GPU_f("elpa_cholesky: tmat to tmatr_dev", 1271,  successGPU)
    endif

    if (useGPU) then
      do i=0,(na-1)/tile_size
        lcs = max(l_colx,i*l_cols_tile+1)
        lce = min(l_cols,(i+1)*l_cols_tile)
        lrs = l_rowx
        lre = min(l_rows,(i+1)*l_rows_tile)
        if (lce  < lcs .or. lre < lrs) cycle
        call obj%timer%start("gpublas")
        gpublasHandle = obj%gpu_setup%gpublasHandleArray(0)
        tmatr_off = (lrs-1 + (1-1)*l_rows) * size_of_datatype
        tmatc_off = (lcs-1 + (1-1)*l_cols) * size_of_datatype
        a_off = (lrs-1 + (lcs-1)*matrixRows) * size_of_datatype
        call gpublas_CGEMM('N', 'C', lre-lrs+1, lce-lcs+1, nblk, &
                            -ONE, tmatr_dev+tmatr_off, l_rows, tmatc_dev+tmatc_off, l_cols, ONE, &
                            a_dev+a_off, matrixRows, gpublasHandle)
        if (wantDebug .and. gpu_vendor() /= SYCL_GPU) successGPU = gpu_DeviceSynchronize()
        call obj%timer%stop("gpublas")
      enddo
    else !useGPU
      do i=0,(na-1)/tile_size
        lcs = max(l_colx,i*l_cols_tile+1)
        lce = min(l_cols,(i+1)*l_cols_tile)
        lrs = l_rowx
        lre = min(l_rows,(i+1)*l_rows_tile)
        if (lce<lcs .or. lre<lrs) cycle
        call obj%timer%start("blas")
        call obj%timer%stop("blas")
      enddo
    endif ! useGPU


  enddo ! n = 1, na, nblk

  call obj%timer%stop("loop1")
  
  call obj%timer%start("deallocate")
  if (useGPU) then



    successGPU = gpu_free(tmp1_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmp1_dev", 1351,  successGPU)

    successGPU = gpu_free(tmp2_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmp1_dev", 1354,  successGPU)

    successGPU = gpu_free(tmatc_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmatc_dev", 1357,  successGPU)

    successGPU = gpu_free(tmatr_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmatr_dev", 1360,  successGPU)

  endif


  deallocate(tmp1, tmp2, tmatr, tmatc, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("elpa_cholesky: tmp1, tmp2, tmatr, tmatc", 1384,  istat,  errorMessage)

  call obj%timer%stop("deallocate")


  call obj%timer%start("loop2")

  call obj%timer%start("copy1")

!  if (useGPU) then
!          ! is this needed: guess not
!#if !defined(WITH_NVIDIA_CUSOLVER) && !defined(WITH_AMD_ROCSOLVER)
!    num = matrixRows*matrixCols* size_of_datatype
!#ifdef WITH_GPU_STREAMS 
!    my_stream = obj%gpu_setup%my_stream
!    call gpu_memcpy_async_and_stream_synchronize &
!    ("elpa_choleksky: a_dev -> a_tmp", a_dev, 0_c_intptr_t, &
!                                    a_tmp(1:obj%local_nrows,1:obj%local_ncols), &
!                                    1, 1, num, gpuMemcpyDeviceToHost, my_stream, .false., .true., .false.)
!#else
!    successGPU = gpu_memcpy(int(loc(a_tmp(1,1)),kind=c_intptr_t), a_dev,  &
!                     matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
!    call check_memcpy_GPU_f("elpa_cholesky: memcpy 2 a-> a_dev", 1422,  successGPU)
!#endif
!#endif 
!  endif
  ! Set the lower triangle to 0, it contains garbage (form the above matrix multiplications)

  call obj%timer%stop("copy1")
  do i=1,na
    if (my_pcol==pcol(i, nblk, np_cols)) then
      ! column i is on local processor
      l_col1 = local_index(i  , my_pcol, np_cols, nblk, +1) ! local column number
      l_row1 = local_index(i+1, my_prow, np_rows, nblk, +1) ! first row below diagonal
      if (useGPU) then
        call obj%timer%start("copy2")
!#ifndef 
        ! set to zero maybe a kerne would be better here
          offset = (l_row1-1 + matrixRows * (l_col1-1)) * size_of_datatype
          num = (l_rows-l_row1+1)*size_of_datatype
          successGPU = gpu_memset(a_dev+offset, 0, num)
          call check_memcpy_GPU_f("elpa_cholesky: memset a_dev", 1450,  successGPU)
!#else 
!          !offset = (l_row1-1 + matrixRows * (l_col1-1)) * size_of_datatype
!          !num = (l_rows-l_row1+1)*size_of_datatype
!#ifdef WITH_GPU_STREAMS
!          !my_stream = obj%gpu_setup%my_stream
!          !successGPU = gpu_memset_async(a_dev+offset, 0, num, my_stream)
!          !call check_memcpy_GPU_f("elpa_cholesky: memset a_dev", 1458,  successGPU)
!
!          !successGPU = gpu_stream_synchronize(my_stream)
!          !call check_memcpy_GPU_f("elpa_cholesky: memset", 1461,  successGPU)
!#else
!          !successGPU = gpu_memset(a_dev+offset, 0, num)
!          !call check_memcpy_GPU_f("elpa_cholesky: memset tmp1_dev", 1464,  successGPU)
!#endif
!
!
!!#if defined(WITH_NVIDIA_CUSOLVER) || defined(WITH_AMD_ROCSOLVER)
!!        successGPU = gpu_memcpy(int(loc(a_tmp(1,1)),kind=c_intptr_t), a_dev, &
!!                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
!!        call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a_tmp", 1471,  successGPU)
!!#endif
!!        a_tmp(l_row1:l_rows,l_col1) = 0
!!#if defined(WITH_NVIDIA_CUSOLVER) || defined(WITH_AMD_ROCSOLVER)
!!        successGPU = gpu_memcpy(a_dev, int(loc(a_tmp(1,1)),kind=c_intptr_t), &
!!                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
!!        call check_memcpy_GPU_f("elpa_cholesky: memcpy a_tmp-> a_dev", 1477,  successGPU)
!#endif

!#endif 

  call obj%timer%stop("copy2")

      else ! useGPU
      endif ! useGPU
    endif ! (my_pcol==pcol(i, nblk, np_cols))
  enddo
  !endif ! useGPU

  call obj%timer%stop("loop2")


  call obj%timer%start("cleanup")
!#if !defined(WITH_NVIDIA_CUSOLVER) && !defined(WITH_AMD_ROCSOLVER)
!    num = matrixRows*matrixCols* size_of_datatype
!#ifdef WITH_GPU_STREAMS
!    my_stream = obj%gpu_setup%my_stream
!    call gpu_memcpy_async_and_stream_synchronize &
!    ("elpa_choleksky: a_tmp -> a_dev", a_dev, 0_c_intptr_t, &
!                                    a_tmp(1:obj%local_nrows,1:obj%local_ncols), &
!                                    1, 1, num, gpuMemcpyHostToDevice, my_stream, .false., .true., .false.)
!#else
!    successGPU = gpu_memcpy(a_dev, int(loc(a_tmp(1,1)),kind=c_intptr_t), &
!                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
!    call check_memcpy_GPU_f("elpa_cholesky: memcpy a_tmp-> a_dev", 1525,  successGPU)
!#endif
!#endif 

!#if !defined(WITH_NVIDIA_CUSOLVER) && !defined(WITH_AMD_ROCSOLVER)

    deallocate(a_tmp, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("elpa_cholesky: a_tmp", 1536,  istat,  errorMessage)
!#endif   

  ! restore original OpenMP settings
  call obj%timer%stop("cleanup")
  call obj%timer%stop("elpa_cholesky_&
  &complex&
  &_&
  &single&
  &"//gpuString)

    end function elpa_cholesky_d_ptr_complex_single_impl



end module
