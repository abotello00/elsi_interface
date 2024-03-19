










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
! ELPA1 -- Faster replacements for ScaLAPACK symmetric eigenvalue routines
!
! Copyright of the original code rests with the authors inside the ELPA
! consortium. The copyright of any additional modifications shall rest
! with their original authors, but shall adhere to the licensing terms
! distributed along with the original code in the file "COPYING".

!> \brief Fortran module which contains the source of ELPA 1stage
module elpa1_compute
  use elpa_utilities
  use elpa_mpi
  implicit none

  PRIVATE ! set default to private

  public :: tridiag_real_double               ! Transform real symmetric matrix to tridiagonal form
  public :: tridiag_real
  public :: trans_ev_real_double              ! Transform real eigenvectors of a tridiagonal matrix back
  public :: trans_ev_real

  !public :: solve_tridi_double
  public :: solve_tridi_double_impl

  interface tridiag_real
    module procedure tridiag_real_double
  end interface

  interface trans_ev_real
    module procedure trans_ev_real_double
  end interface

  public :: tridiag_real_single        ! Transform real single-precision symmetric matrix to tridiagonal form
  public :: trans_ev_real_single       ! Transform real  single-precision eigenvectors of a tridiagonal matrix back
  !public :: solve_tridi_single
  public :: solve_tridi_single_impl

  public :: tridiag_complex_double            ! Transform complex hermitian matrix to tridiagonal form
  public :: tridiag_complex
  public :: trans_ev_complex_double           ! Transform eigenvectors of a tridiagonal matrix back
  public :: trans_ev_complex

  interface tridiag_complex
    module procedure tridiag_complex_double
  end interface

  interface trans_ev_complex
    module procedure trans_ev_complex_double
  end interface

  public :: tridiag_complex_single     ! Transform complex single-precision hermitian matrix to tridiagonal form
  public :: trans_ev_complex_single    ! Transform complex single-precision eigenvectors of a tridiagonal matrix back

  public :: hh_transform_real_double
  public :: hh_transform_real
  public :: elpa_reduce_add_vectors_real_double
  public :: elpa_reduce_add_vectors_real
  public :: elpa_transpose_vectors_real_double
  public :: elpa_transpose_vectors_ss_real_double
  public :: elpa_transpose_vectors_real
  public :: elpa_transpose_vectors_ss_real

  interface hh_transform_real
    module procedure hh_transform_real_double
  end interface

  interface elpa_reduce_add_vectors_real
    module procedure elpa_reduce_add_vectors_real_double
  end interface

  interface elpa_transpose_vectors_real
    module procedure elpa_transpose_vectors_real_double
  end interface
  
  interface elpa_transpose_vectors_ss_real
    module procedure elpa_transpose_vectors_ss_real_double
  end interface

  public :: hh_transform_real_single
  public :: elpa_reduce_add_vectors_real_single
  public :: elpa_transpose_vectors_real_single
  public :: elpa_transpose_vectors_ss_real_single

  public :: hh_transform_complex_double
  public :: hh_transform_complex
  public :: elpa_reduce_add_vectors_complex_double
  public :: elpa_reduce_add_vectors_complex
  public :: elpa_transpose_vectors_complex_double
  public :: elpa_transpose_vectors_ss_complex_double
  public :: elpa_transpose_vectors_complex
  public :: elpa_transpose_vectors_ss_complex

  interface hh_transform_complex
    module procedure hh_transform_complex_double
  end interface

  interface elpa_reduce_add_vectors_complex
    module procedure elpa_reduce_add_vectors_complex_double
  end interface

  interface elpa_transpose_vectors_complex
    module procedure elpa_transpose_vectors_complex_double
  end interface
  
  interface elpa_transpose_vectors_ss_complex
    module procedure elpa_transpose_vectors_ss_complex_double
  end interface

  public :: hh_transform_complex_single
  public :: elpa_reduce_add_vectors_complex_single
  public :: elpa_transpose_vectors_complex_single
  public :: elpa_transpose_vectors_ss_complex_single


  contains

!________________________________________________________________
! elpa_transpose_vectors_template.F90
! elpa_reduce_add_vectors_template.F90

! real double precision
!#define DOUBLE_PRECISION_REAL 1




















!cannot use "../src/elpa1/../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)




subroutine elpa_transpose_vectors_&
&real&
&_&
&double &
(obj, vmat_s, ld_s, comm_s, vmat_t, ld_t, comm_t, & 
 nvs, nvr, nvc, nblk, nrThreads, comm_s_isRows, success)

!-------------------------------------------------------------------------------
! This routine transposes an array of vectors which are distributed in
! communicator comm_s into its transposed form distributed in communicator comm_t.
! There must be an identical copy of vmat_s in every communicator comm_s.
! After this routine, there is an identical copy of vmat_t in every communicator comm_t.
!
! vmat_s    original array of vectors
! ld_s      leading dimension of vmat_s
! comm_s    communicator over which vmat_s is distributed
! vmat_t    array of vectors in transposed form
! ld_t      leading dimension of vmat_t
! comm_t    communicator over which vmat_t is distributed
! nvs       global index where to start in vmat_s/vmat_t
!           Please note: this is kind of a hint, some values before nvs will be
!           accessed in vmat_s/put into vmat_t
! nvr       global length of vmat_s/vmat_t
! nvc       number of columns in vmat_s/vmat_t
! nblk      block size of block cyclic distribution
!
!-------------------------------------------------------------------------------
  use precision
  use elpa_abstract_impl
  use elpa_mpi

  implicit none
  class(elpa_abstract_impl_t), intent(inout) :: obj
  integer(kind=ik), intent(in)                      :: ld_s, comm_s, ld_t, comm_t, nvs, nvr, nvc, nblk
  real(kind=c_double), intent(in)   :: vmat_s(ld_s,nvc)
  real(kind=c_double), intent(inout):: vmat_t(ld_t,nvc)

  real(kind=c_double), allocatable  :: aux(:)
  integer(kind=ik)                                  :: myps, mypt, nps, npt
  integer(kind=MPI_KIND)                            :: mypsMPI, myptMPI, npsMPI, nptMPI
  integer(kind=ik)                                  :: n, lc, k, i, ips, ipt, ns, nl
  integer(kind=MPI_KIND)                            :: mpierr
  integer(kind=ik)                                  :: lcm_s_t, nblks_tot, nblks_comm, nblks_skip
  integer(kind=ik)                                  :: aux_stride
  integer(kind=ik), intent(in)                      :: nrThreads
  integer(kind=ik)                                  :: istat
  character(200)                                    :: errorMessage

  integer(kind=MPI_KIND)                            :: bcast_request1
  logical                                           :: useNonBlockingCollectives
  logical                                           :: useNonBlockingCollectivesRows
  logical                                           :: useNonBlockingCollectivesCols
  logical, intent(in)                               :: comm_s_isRows
  integer(kind=c_int)                               :: non_blocking_collectives_rows, error, &
                                                       non_blocking_collectives_cols
  logical                                           :: success

  integer(kind=ik)                                  :: mpi_comm_all, my_mpi_rank, transposed_mpi_rank, message_size, matrix_order
  integer(kind=ik)                                  :: ld_st, solver

  success = .true.

  call obj%timer%start("&
          &elpa_transpose_vectors_&
  &real&
  &" // &
  &"_double" &
  )

  call obj%get("nbc_row_transpose_vectors", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    write(error_unit,*) "Problem setting option for non blocking collectives for_rows in transpose_vectors. Aborting..."
    success = .false.
     call obj%timer%stop("&
          &elpa_transpose_vectors_&
     &real&
     &" // &
     &"_double" &
    )
    return
  endif

  call obj%get("nbc_col_transpose_vectors", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    write(error_unit,*) "Problem setting option for non blocking collectives for_cols in transpose_vectors. Aborting..."
    success = .false.
     call obj%timer%stop("&
          &elpa_transpose_vectors_&
     &real&
     &" // &
     &"_double" &
    )
    return
  endif

  if (non_blocking_collectives_rows .eq. 1) then
    useNonBlockingCollectivesRows = .true.
  else
    useNonBlockingCollectivesRows = .false.
  endif

  if (non_blocking_collectives_cols .eq. 1) then
    useNonBlockingCollectivesCols = .true.
  else
    useNonBlockingCollectivesCols = .false.
  endif

  if (comm_s_isRows) then
    useNonBlockingCollectives = useNonBlockingCollectivesRows
  else
    useNonBlockingCollectives = useNonBlockingCollectivesCols
  endif

  call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(comm_s,kind=MPI_KIND),mypsMPI, mpierr)
  call mpi_comm_size(int(comm_s,kind=MPI_KIND),npsMPI ,mpierr)
  call mpi_comm_rank(int(comm_t,kind=MPI_KIND),myptMPI, mpierr)
  call mpi_comm_size(int(comm_t,kind=MPI_KIND),nptMPI ,mpierr)
  myps = int(mypsMPI,kind=c_int)
  nps = int(npsMPI,kind=c_int)
  mypt = int(myptMPI,kind=c_int)
  npt = int(nptMPI,kind=c_int)
  call obj%timer%stop("mpi_communication")

  ! TODO_23_11
  ! this codepath doesn't work for ELPA2, Cholesky, maybe smth else. Fix it, along with optimization of Cholesky-GPU
  ! because there nvc>1 and ld_s != ld_t (so, we can't make a contigous-memory MPI_Send call)

  ! The basic idea of this routine is that for every block (in the block cyclic
  ! distribution), the processor within comm_t which owns the diagonal
  ! broadcasts its values of vmat_s to all processors within comm_t.
  ! Of course this has not to be done for every block separately, since
  ! the communictation pattern repeats in the global matrix after
  ! the least common multiple of (nps,npt) blocks

  lcm_s_t   = least_common_multiple(nps,npt) ! least common multiple of nps, npt

  nblks_tot = (nvr+nblk-1)/nblk ! number of blocks corresponding to nvr

  ! Get the number of blocks to be skipped at the begin.
  ! This must be a multiple of lcm_s_t (else it is getting complicated),
  ! thus some elements before nvs will be accessed/set.

  nblks_skip = ((nvs-1)/(nblk*lcm_s_t))*lcm_s_t

  allocate(aux( ((nblks_tot-nblks_skip+lcm_s_t-1)/lcm_s_t) * nblk * nvc ), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_transpose_vectors: aux", 275,  istat,  errorMessage)
  do n = 0, lcm_s_t-1

    ips = mod(n,nps)
    ipt = mod(n,npt)

    if (mypt == ipt) then ! (mypt == ipt)

      nblks_comm = (nblks_tot-nblks_skip-n+lcm_s_t-1)/lcm_s_t
      aux_stride = nblk * nblks_comm
!      if(nblks_comm==0) cycle
      if (nblks_comm .ne. 0) then
        if (myps == ips) then
!          k = 0
          do lc=1,nvc
            do i = nblks_skip+n, nblks_tot-1, lcm_s_t
              k = (i - nblks_skip - n)/lcm_s_t * nblk + (lc - 1) * aux_stride
              ns = (i/nps)*nblk ! local start of block i
              nl = min(nvr-i*nblk,nblk) ! length
              aux(k+1:k+nl) = vmat_s(ns+1:ns+nl,lc)
!              k = k+nblk
            enddo
          enddo
        endif ! (myps == ips)



!        k = 0
        do lc=1,nvc
          do i = nblks_skip+n, nblks_tot-1, lcm_s_t
            k = (i - nblks_skip - n)/lcm_s_t * nblk + (lc - 1) * aux_stride
            ns = (i/npt)*nblk ! local start of block i
            nl = min(nvr-i*nblk,nblk) ! length
            vmat_t(ns+1:ns+nl,lc) = aux(k+1:k+nl)
!            k = k+nblk
          enddo
        enddo
      endif
    endif ! (mypt == ipt)

  enddo
  deallocate(aux, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("elpa_transpose_vectors: aux", 374,  istat,  errorMessage)

  call obj%timer%stop("&
  &elpa_transpose_vectors_&
  &real&
  &" // &
  &"_double" &
  )

end subroutine




!cannot use "../src/elpa1/../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)




subroutine elpa_transpose_vectors_ss_&
&real&
&_&
&double &
(obj, vmat_s, ld_s, comm_s, vmat_t, ld_t, comm_t, & 
 nvs, nvr, nvc, nblk, nrThreads, comm_s_isRows, success)

!-------------------------------------------------------------------------------
! This routine transposes an array of vectors which are distributed in
! communicator comm_s into its transposed form distributed in communicator comm_t.
! There must be an identical copy of vmat_s in every communicator comm_s.
! After this routine, there is an identical copy of vmat_t in every communicator comm_t.
!
! vmat_s    original array of vectors
! ld_s      leading dimension of vmat_s
! comm_s    communicator over which vmat_s is distributed
! vmat_t    array of vectors in transposed form
! ld_t      leading dimension of vmat_t
! comm_t    communicator over which vmat_t is distributed
! nvs       global index where to start in vmat_s/vmat_t
!           Please note: this is kind of a hint, some values before nvs will be
!           accessed in vmat_s/put into vmat_t
! nvr       global length of vmat_s/vmat_t
! nvc       number of columns in vmat_s/vmat_t
! nblk      block size of block cyclic distribution
!
!-------------------------------------------------------------------------------
  use precision
  use elpa_abstract_impl
  use elpa_mpi

  implicit none
  class(elpa_abstract_impl_t), intent(inout) :: obj
  integer(kind=ik), intent(in)                      :: ld_s, comm_s, ld_t, comm_t, nvs, nvr, nvc, nblk
  real(kind=c_double), intent(in)   :: vmat_s(ld_s,nvc)
  real(kind=c_double), intent(inout):: vmat_t(ld_t,nvc)

  real(kind=c_double), allocatable  :: aux(:)
  integer(kind=ik)                                  :: myps, mypt, nps, npt
  integer(kind=MPI_KIND)                            :: mypsMPI, myptMPI, npsMPI, nptMPI
  integer(kind=ik)                                  :: n, lc, k, i, ips, ipt, ns, nl
  integer(kind=MPI_KIND)                            :: mpierr
  integer(kind=ik)                                  :: lcm_s_t, nblks_tot, nblks_comm, nblks_skip
  integer(kind=ik)                                  :: aux_stride
  integer(kind=ik), intent(in)                      :: nrThreads
  integer(kind=ik)                                  :: istat
  character(200)                                    :: errorMessage

  integer(kind=MPI_KIND)                            :: bcast_request1
  logical                                           :: useNonBlockingCollectives
  logical                                           :: useNonBlockingCollectivesRows
  logical                                           :: useNonBlockingCollectivesCols
  logical, intent(in)                               :: comm_s_isRows
  integer(kind=c_int)                               :: non_blocking_collectives_rows, error, &
                                                       non_blocking_collectives_cols
  logical                                           :: success

  integer(kind=ik)                                  :: mpi_comm_all, my_mpi_rank, transposed_mpi_rank, message_size, matrix_order
  integer(kind=ik)                                  :: ld_st, solver

  success = .true.

  call obj%timer%start("&
          &elpa_transpose_vectors_ss_&
  &real&
  &" // &
  &"_double" &
  )

  call obj%get("nbc_row_transpose_vectors", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    write(error_unit,*) "Problem setting option for non blocking collectives for_rows in transpose_vectors. Aborting..."
    success = .false.
     call obj%timer%stop("&
          &elpa_transpose_vectors_ss_&
     &real&
     &" // &
     &"_double" &
    )
    return
  endif

  call obj%get("nbc_col_transpose_vectors", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    write(error_unit,*) "Problem setting option for non blocking collectives for_cols in transpose_vectors. Aborting..."
    success = .false.
     call obj%timer%stop("&
          &elpa_transpose_vectors_ss_&
     &real&
     &" // &
     &"_double" &
    )
    return
  endif

  if (non_blocking_collectives_rows .eq. 1) then
    useNonBlockingCollectivesRows = .true.
  else
    useNonBlockingCollectivesRows = .false.
  endif

  if (non_blocking_collectives_cols .eq. 1) then
    useNonBlockingCollectivesCols = .true.
  else
    useNonBlockingCollectivesCols = .false.
  endif

  if (comm_s_isRows) then
    useNonBlockingCollectives = useNonBlockingCollectivesRows
  else
    useNonBlockingCollectives = useNonBlockingCollectivesCols
  endif

  call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(comm_s,kind=MPI_KIND),mypsMPI, mpierr)
  call mpi_comm_size(int(comm_s,kind=MPI_KIND),npsMPI ,mpierr)
  call mpi_comm_rank(int(comm_t,kind=MPI_KIND),myptMPI, mpierr)
  call mpi_comm_size(int(comm_t,kind=MPI_KIND),nptMPI ,mpierr)
  myps = int(mypsMPI,kind=c_int)
  nps = int(npsMPI,kind=c_int)
  mypt = int(myptMPI,kind=c_int)
  npt = int(nptMPI,kind=c_int)
  call obj%timer%stop("mpi_communication")

  ! TODO_23_11
  ! this codepath doesn't work for ELPA2, Cholesky, maybe smth else. Fix it, along with optimization of Cholesky-GPU
  ! because there nvc>1 and ld_s != ld_t (so, we can't make a contigous-memory MPI_Send call)

  ! The basic idea of this routine is that for every block (in the block cyclic
  ! distribution), the processor within comm_t which owns the diagonal
  ! broadcasts its values of vmat_s to all processors within comm_t.
  ! Of course this has not to be done for every block separately, since
  ! the communictation pattern repeats in the global matrix after
  ! the least common multiple of (nps,npt) blocks

  lcm_s_t   = least_common_multiple(nps,npt) ! least common multiple of nps, npt

  nblks_tot = (nvr+nblk-1)/nblk ! number of blocks corresponding to nvr

  ! Get the number of blocks to be skipped at the begin.
  ! This must be a multiple of lcm_s_t (else it is getting complicated),
  ! thus some elements before nvs will be accessed/set.

  nblks_skip = ((nvs-1)/(nblk*lcm_s_t))*lcm_s_t

  allocate(aux( ((nblks_tot-nblks_skip+lcm_s_t-1)/lcm_s_t) * nblk * nvc ), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_transpose_vectors: aux", 275,  istat,  errorMessage)
  do n = 0, lcm_s_t-1

    ips = mod(n,nps)
    ipt = mod(n,npt)

    if (mypt == ipt) then ! (mypt == ipt)

      nblks_comm = (nblks_tot-nblks_skip-n+lcm_s_t-1)/lcm_s_t
      aux_stride = nblk * nblks_comm
!      if(nblks_comm==0) cycle
      if (nblks_comm .ne. 0) then
        if (myps == ips) then
!          k = 0
          do lc=1,nvc
            do i = nblks_skip+n, nblks_tot-1, lcm_s_t
              k = (i - nblks_skip - n)/lcm_s_t * nblk + (lc - 1) * aux_stride
              ns = (i/nps)*nblk ! local start of block i
              nl = min(nvr-i*nblk,nblk) ! length
              aux(k+1:k+nl) = vmat_s(ns+1:ns+nl,lc)
!              k = k+nblk
            enddo
          enddo
        endif ! (myps == ips)



!        k = 0
        do lc=1,nvc
          do i = nblks_skip+n, nblks_tot-1, lcm_s_t
            k = (i - nblks_skip - n)/lcm_s_t * nblk + (lc - 1) * aux_stride
            ns = (i/npt)*nblk ! local start of block i
            nl = min(nvr-i*nblk,nblk) ! length
            vmat_t(ns+1:ns+nl,lc) = - aux(k+1:k+nl)
!            k = k+nblk
          enddo
        enddo
      endif
    endif ! (mypt == ipt)

  enddo
  deallocate(aux, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("elpa_transpose_vectors: aux", 374,  istat,  errorMessage)

  call obj%timer%stop("&
  &elpa_transpose_vectors_ss_&
  &real&
  &" // &
  &"_double" &
  )

end subroutine




!cannot use "../src/elpa1/../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)


subroutine elpa_reduce_add_vectors_&
&real&
&_&
&double &
(obj, vmat_s, ld_s, comm_s, vmat_t, ld_t, comm_t, nvr, nvc, nblk, nrThreads)

!-------------------------------------------------------------------------------
! This routine does a reduce of all vectors in vmat_s over the communicator comm_t.
! The result of the reduce is gathered on the processors owning the diagonal
! and added to the array of vectors vmat_t (which is distributed over comm_t).
!
! Opposed to elpa_transpose_vectors, there is NO identical copy of vmat_s
! in the different members within vmat_t (else a reduce wouldn't be necessary).
! After this routine, an allreduce of vmat_t has to be done.
!
! vmat_s    array of vectors to be reduced and added
! ld_s      leading dimension of vmat_s
! comm_s    communicator over which vmat_s is distributed
! vmat_t    array of vectors to which vmat_s is added
! ld_t      leading dimension of vmat_t
! comm_t    communicator over which vmat_t is distributed
! nvr       global length of vmat_s/vmat_t
! nvc       number of columns in vmat_s/vmat_t
! nblk      block size of block cyclic distribution
!
!-------------------------------------------------------------------------------

  use precision
  use elpa_mpi
  use elpa_abstract_impl
  implicit none

  class(elpa_abstract_impl_t), intent(inout)         :: obj
  integer(kind=ik), intent(in)                       :: ld_s, comm_s, ld_t, comm_t, nvr, nvc, nblk
  real(kind=c_double), intent(in)    :: vmat_s(ld_s,nvc)
  real(kind=c_double), intent(inout) :: vmat_t(ld_t,nvc)

  real(kind=c_double), allocatable   :: aux1(:), aux2(:)
  integer(kind=ik)                                   :: myps, mypt, nps, npt
  integer(kind=MPI_KIND)                             :: mypsMPI, npsMPI, myptMPI, nptMPI
  integer(kind=ik)                                   :: n, lc, k, i, ips, ipt, ns, nl
  integer(kind=MPI_KIND)                             :: mpierr
  integer(kind=ik)                                   :: lcm_s_t, nblks_tot
  integer(kind=ik)                                   :: aux_stride
  integer(kind=ik), intent(in)                       :: nrThreads
  integer(kind=ik)                                   :: istat
  character(200)                                     :: errorMessage

  call obj%timer%start("elpa_reduce_add_vectors_&
  &real&
  &" // &
  &"_double" &
  )

  call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(comm_s,kind=MPI_KIND), mypsMPI, mpierr)
  call mpi_comm_size(int(comm_s,kind=MPI_KIND), npsMPI,  mpierr)
  call mpi_comm_rank(int(comm_t,kind=MPI_KIND), myptMPI, mpierr)
  call mpi_comm_size(int(comm_t,kind=MPI_KIND), nptMPI ,mpierr)
  myps = int(mypsMPI,kind=c_int)
  nps = int(npsMPI,kind=c_int)
  mypt = int(myptMPI,kind=c_int)
  npt = int(nptMPI,kind=c_int)

  call obj%timer%stop("mpi_communication")

  ! Look to elpa_transpose_vectors for the basic idea!

  ! The communictation pattern repeats in the global matrix after
  ! the least common multiple of (nps,npt) blocks

  lcm_s_t   = least_common_multiple(nps,npt) ! least common multiple of nps, npt

  nblks_tot = (nvr+nblk-1)/nblk ! number of blocks corresponding to nvr

  allocate(aux1( ((nblks_tot+lcm_s_t-1)/lcm_s_t) * nblk * nvc ), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_reduce_add: aux1", 129,  istat,  errorMessage)

  allocate(aux2( ((nblks_tot+lcm_s_t-1)/lcm_s_t) * nblk * nvc ), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_reduce_add: aux2", 132,  istat,  errorMessage)
  aux1(:) = 0
  aux2(:) = 0
  do n = 0, lcm_s_t-1

    ips = mod(n,nps)
    ipt = mod(n,npt)

    aux_stride = nblk * ((nblks_tot - n + lcm_s_t - 1)/lcm_s_t)

    if (myps == ips) then

!      k = 0
      do lc=1,nvc
        do i = n, nblks_tot-1, lcm_s_t
          k = (i - n)/lcm_s_t * nblk + (lc - 1) * aux_stride
          ns = (i/nps)*nblk ! local start of block i
          nl = min(nvr-i*nblk,nblk) ! length
          aux1(k+1:k+nl) = vmat_s(ns+1:ns+nl,lc)
!          k = k+nblk
        enddo
      enddo

      k = nvc * aux_stride

      if (k>0) aux2 = aux1

      if (mypt == ipt) then
!        k = 0
        do lc=1,nvc
          do  i = n, nblks_tot-1, lcm_s_t
            k = (i - n)/lcm_s_t * nblk + (lc - 1) * aux_stride
            ns = (i/npt)*nblk ! local start of block i
            nl = min(nvr-i*nblk,nblk) ! length
            vmat_t(ns+1:ns+nl,lc) = vmat_t(ns+1:ns+nl,lc) + aux2(k+1:k+nl)
!            k = k+nblk
          enddo
        enddo
      endif ! (mypt == ipt)

    endif ! (myps == ips)

  enddo

  deallocate(aux1, aux2, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("elpa_reduce_add: aux1, aux2", 218,  istat,  errorMessage)

  call obj%timer%stop("elpa_reduce_add_vectors_&
  &real&
  &" // &
  &"_double" &
  )
end subroutine



! real single precision






















!cannot use "../src/elpa1/../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)




subroutine elpa_transpose_vectors_&
&real&
&_&
&single &
(obj, vmat_s, ld_s, comm_s, vmat_t, ld_t, comm_t, & 
 nvs, nvr, nvc, nblk, nrThreads, comm_s_isRows, success)

!-------------------------------------------------------------------------------
! This routine transposes an array of vectors which are distributed in
! communicator comm_s into its transposed form distributed in communicator comm_t.
! There must be an identical copy of vmat_s in every communicator comm_s.
! After this routine, there is an identical copy of vmat_t in every communicator comm_t.
!
! vmat_s    original array of vectors
! ld_s      leading dimension of vmat_s
! comm_s    communicator over which vmat_s is distributed
! vmat_t    array of vectors in transposed form
! ld_t      leading dimension of vmat_t
! comm_t    communicator over which vmat_t is distributed
! nvs       global index where to start in vmat_s/vmat_t
!           Please note: this is kind of a hint, some values before nvs will be
!           accessed in vmat_s/put into vmat_t
! nvr       global length of vmat_s/vmat_t
! nvc       number of columns in vmat_s/vmat_t
! nblk      block size of block cyclic distribution
!
!-------------------------------------------------------------------------------
  use precision
  use elpa_abstract_impl
  use elpa_mpi

  implicit none
  class(elpa_abstract_impl_t), intent(inout) :: obj
  integer(kind=ik), intent(in)                      :: ld_s, comm_s, ld_t, comm_t, nvs, nvr, nvc, nblk
  real(kind=c_float), intent(in)   :: vmat_s(ld_s,nvc)
  real(kind=c_float), intent(inout):: vmat_t(ld_t,nvc)

  real(kind=c_float), allocatable  :: aux(:)
  integer(kind=ik)                                  :: myps, mypt, nps, npt
  integer(kind=MPI_KIND)                            :: mypsMPI, myptMPI, npsMPI, nptMPI
  integer(kind=ik)                                  :: n, lc, k, i, ips, ipt, ns, nl
  integer(kind=MPI_KIND)                            :: mpierr
  integer(kind=ik)                                  :: lcm_s_t, nblks_tot, nblks_comm, nblks_skip
  integer(kind=ik)                                  :: aux_stride
  integer(kind=ik), intent(in)                      :: nrThreads
  integer(kind=ik)                                  :: istat
  character(200)                                    :: errorMessage

  integer(kind=MPI_KIND)                            :: bcast_request1
  logical                                           :: useNonBlockingCollectives
  logical                                           :: useNonBlockingCollectivesRows
  logical                                           :: useNonBlockingCollectivesCols
  logical, intent(in)                               :: comm_s_isRows
  integer(kind=c_int)                               :: non_blocking_collectives_rows, error, &
                                                       non_blocking_collectives_cols
  logical                                           :: success

  integer(kind=ik)                                  :: mpi_comm_all, my_mpi_rank, transposed_mpi_rank, message_size, matrix_order
  integer(kind=ik)                                  :: ld_st, solver

  success = .true.

  call obj%timer%start("&
          &elpa_transpose_vectors_&
  &real&
  &" // &
  &"_single" &
  )

  call obj%get("nbc_row_transpose_vectors", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    write(error_unit,*) "Problem setting option for non blocking collectives for_rows in transpose_vectors. Aborting..."
    success = .false.
     call obj%timer%stop("&
          &elpa_transpose_vectors_&
     &real&
     &" // &
     &"_single" &
    )
    return
  endif

  call obj%get("nbc_col_transpose_vectors", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    write(error_unit,*) "Problem setting option for non blocking collectives for_cols in transpose_vectors. Aborting..."
    success = .false.
     call obj%timer%stop("&
          &elpa_transpose_vectors_&
     &real&
     &" // &
     &"_single" &
    )
    return
  endif

  if (non_blocking_collectives_rows .eq. 1) then
    useNonBlockingCollectivesRows = .true.
  else
    useNonBlockingCollectivesRows = .false.
  endif

  if (non_blocking_collectives_cols .eq. 1) then
    useNonBlockingCollectivesCols = .true.
  else
    useNonBlockingCollectivesCols = .false.
  endif

  if (comm_s_isRows) then
    useNonBlockingCollectives = useNonBlockingCollectivesRows
  else
    useNonBlockingCollectives = useNonBlockingCollectivesCols
  endif

  call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(comm_s,kind=MPI_KIND),mypsMPI, mpierr)
  call mpi_comm_size(int(comm_s,kind=MPI_KIND),npsMPI ,mpierr)
  call mpi_comm_rank(int(comm_t,kind=MPI_KIND),myptMPI, mpierr)
  call mpi_comm_size(int(comm_t,kind=MPI_KIND),nptMPI ,mpierr)
  myps = int(mypsMPI,kind=c_int)
  nps = int(npsMPI,kind=c_int)
  mypt = int(myptMPI,kind=c_int)
  npt = int(nptMPI,kind=c_int)
  call obj%timer%stop("mpi_communication")

  ! TODO_23_11
  ! this codepath doesn't work for ELPA2, Cholesky, maybe smth else. Fix it, along with optimization of Cholesky-GPU
  ! because there nvc>1 and ld_s != ld_t (so, we can't make a contigous-memory MPI_Send call)

  ! The basic idea of this routine is that for every block (in the block cyclic
  ! distribution), the processor within comm_t which owns the diagonal
  ! broadcasts its values of vmat_s to all processors within comm_t.
  ! Of course this has not to be done for every block separately, since
  ! the communictation pattern repeats in the global matrix after
  ! the least common multiple of (nps,npt) blocks

  lcm_s_t   = least_common_multiple(nps,npt) ! least common multiple of nps, npt

  nblks_tot = (nvr+nblk-1)/nblk ! number of blocks corresponding to nvr

  ! Get the number of blocks to be skipped at the begin.
  ! This must be a multiple of lcm_s_t (else it is getting complicated),
  ! thus some elements before nvs will be accessed/set.

  nblks_skip = ((nvs-1)/(nblk*lcm_s_t))*lcm_s_t

  allocate(aux( ((nblks_tot-nblks_skip+lcm_s_t-1)/lcm_s_t) * nblk * nvc ), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_transpose_vectors: aux", 275,  istat,  errorMessage)
  do n = 0, lcm_s_t-1

    ips = mod(n,nps)
    ipt = mod(n,npt)

    if (mypt == ipt) then ! (mypt == ipt)

      nblks_comm = (nblks_tot-nblks_skip-n+lcm_s_t-1)/lcm_s_t
      aux_stride = nblk * nblks_comm
!      if(nblks_comm==0) cycle
      if (nblks_comm .ne. 0) then
        if (myps == ips) then
!          k = 0
          do lc=1,nvc
            do i = nblks_skip+n, nblks_tot-1, lcm_s_t
              k = (i - nblks_skip - n)/lcm_s_t * nblk + (lc - 1) * aux_stride
              ns = (i/nps)*nblk ! local start of block i
              nl = min(nvr-i*nblk,nblk) ! length
              aux(k+1:k+nl) = vmat_s(ns+1:ns+nl,lc)
!              k = k+nblk
            enddo
          enddo
        endif ! (myps == ips)



!        k = 0
        do lc=1,nvc
          do i = nblks_skip+n, nblks_tot-1, lcm_s_t
            k = (i - nblks_skip - n)/lcm_s_t * nblk + (lc - 1) * aux_stride
            ns = (i/npt)*nblk ! local start of block i
            nl = min(nvr-i*nblk,nblk) ! length
            vmat_t(ns+1:ns+nl,lc) = aux(k+1:k+nl)
!            k = k+nblk
          enddo
        enddo
      endif
    endif ! (mypt == ipt)

  enddo
  deallocate(aux, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("elpa_transpose_vectors: aux", 374,  istat,  errorMessage)

  call obj%timer%stop("&
  &elpa_transpose_vectors_&
  &real&
  &" // &
  &"_single" &
  )

end subroutine




!cannot use "../src/elpa1/../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)




subroutine elpa_transpose_vectors_ss_&
&real&
&_&
&single &
(obj, vmat_s, ld_s, comm_s, vmat_t, ld_t, comm_t, & 
 nvs, nvr, nvc, nblk, nrThreads, comm_s_isRows, success)

!-------------------------------------------------------------------------------
! This routine transposes an array of vectors which are distributed in
! communicator comm_s into its transposed form distributed in communicator comm_t.
! There must be an identical copy of vmat_s in every communicator comm_s.
! After this routine, there is an identical copy of vmat_t in every communicator comm_t.
!
! vmat_s    original array of vectors
! ld_s      leading dimension of vmat_s
! comm_s    communicator over which vmat_s is distributed
! vmat_t    array of vectors in transposed form
! ld_t      leading dimension of vmat_t
! comm_t    communicator over which vmat_t is distributed
! nvs       global index where to start in vmat_s/vmat_t
!           Please note: this is kind of a hint, some values before nvs will be
!           accessed in vmat_s/put into vmat_t
! nvr       global length of vmat_s/vmat_t
! nvc       number of columns in vmat_s/vmat_t
! nblk      block size of block cyclic distribution
!
!-------------------------------------------------------------------------------
  use precision
  use elpa_abstract_impl
  use elpa_mpi

  implicit none
  class(elpa_abstract_impl_t), intent(inout) :: obj
  integer(kind=ik), intent(in)                      :: ld_s, comm_s, ld_t, comm_t, nvs, nvr, nvc, nblk
  real(kind=c_float), intent(in)   :: vmat_s(ld_s,nvc)
  real(kind=c_float), intent(inout):: vmat_t(ld_t,nvc)

  real(kind=c_float), allocatable  :: aux(:)
  integer(kind=ik)                                  :: myps, mypt, nps, npt
  integer(kind=MPI_KIND)                            :: mypsMPI, myptMPI, npsMPI, nptMPI
  integer(kind=ik)                                  :: n, lc, k, i, ips, ipt, ns, nl
  integer(kind=MPI_KIND)                            :: mpierr
  integer(kind=ik)                                  :: lcm_s_t, nblks_tot, nblks_comm, nblks_skip
  integer(kind=ik)                                  :: aux_stride
  integer(kind=ik), intent(in)                      :: nrThreads
  integer(kind=ik)                                  :: istat
  character(200)                                    :: errorMessage

  integer(kind=MPI_KIND)                            :: bcast_request1
  logical                                           :: useNonBlockingCollectives
  logical                                           :: useNonBlockingCollectivesRows
  logical                                           :: useNonBlockingCollectivesCols
  logical, intent(in)                               :: comm_s_isRows
  integer(kind=c_int)                               :: non_blocking_collectives_rows, error, &
                                                       non_blocking_collectives_cols
  logical                                           :: success

  integer(kind=ik)                                  :: mpi_comm_all, my_mpi_rank, transposed_mpi_rank, message_size, matrix_order
  integer(kind=ik)                                  :: ld_st, solver

  success = .true.

  call obj%timer%start("&
          &elpa_transpose_vectors_ss_&
  &real&
  &" // &
  &"_single" &
  )

  call obj%get("nbc_row_transpose_vectors", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    write(error_unit,*) "Problem setting option for non blocking collectives for_rows in transpose_vectors. Aborting..."
    success = .false.
     call obj%timer%stop("&
          &elpa_transpose_vectors_ss_&
     &real&
     &" // &
     &"_single" &
    )
    return
  endif

  call obj%get("nbc_col_transpose_vectors", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    write(error_unit,*) "Problem setting option for non blocking collectives for_cols in transpose_vectors. Aborting..."
    success = .false.
     call obj%timer%stop("&
          &elpa_transpose_vectors_ss_&
     &real&
     &" // &
     &"_single" &
    )
    return
  endif

  if (non_blocking_collectives_rows .eq. 1) then
    useNonBlockingCollectivesRows = .true.
  else
    useNonBlockingCollectivesRows = .false.
  endif

  if (non_blocking_collectives_cols .eq. 1) then
    useNonBlockingCollectivesCols = .true.
  else
    useNonBlockingCollectivesCols = .false.
  endif

  if (comm_s_isRows) then
    useNonBlockingCollectives = useNonBlockingCollectivesRows
  else
    useNonBlockingCollectives = useNonBlockingCollectivesCols
  endif

  call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(comm_s,kind=MPI_KIND),mypsMPI, mpierr)
  call mpi_comm_size(int(comm_s,kind=MPI_KIND),npsMPI ,mpierr)
  call mpi_comm_rank(int(comm_t,kind=MPI_KIND),myptMPI, mpierr)
  call mpi_comm_size(int(comm_t,kind=MPI_KIND),nptMPI ,mpierr)
  myps = int(mypsMPI,kind=c_int)
  nps = int(npsMPI,kind=c_int)
  mypt = int(myptMPI,kind=c_int)
  npt = int(nptMPI,kind=c_int)
  call obj%timer%stop("mpi_communication")

  ! TODO_23_11
  ! this codepath doesn't work for ELPA2, Cholesky, maybe smth else. Fix it, along with optimization of Cholesky-GPU
  ! because there nvc>1 and ld_s != ld_t (so, we can't make a contigous-memory MPI_Send call)

  ! The basic idea of this routine is that for every block (in the block cyclic
  ! distribution), the processor within comm_t which owns the diagonal
  ! broadcasts its values of vmat_s to all processors within comm_t.
  ! Of course this has not to be done for every block separately, since
  ! the communictation pattern repeats in the global matrix after
  ! the least common multiple of (nps,npt) blocks

  lcm_s_t   = least_common_multiple(nps,npt) ! least common multiple of nps, npt

  nblks_tot = (nvr+nblk-1)/nblk ! number of blocks corresponding to nvr

  ! Get the number of blocks to be skipped at the begin.
  ! This must be a multiple of lcm_s_t (else it is getting complicated),
  ! thus some elements before nvs will be accessed/set.

  nblks_skip = ((nvs-1)/(nblk*lcm_s_t))*lcm_s_t

  allocate(aux( ((nblks_tot-nblks_skip+lcm_s_t-1)/lcm_s_t) * nblk * nvc ), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_transpose_vectors: aux", 275,  istat,  errorMessage)
  do n = 0, lcm_s_t-1

    ips = mod(n,nps)
    ipt = mod(n,npt)

    if (mypt == ipt) then ! (mypt == ipt)

      nblks_comm = (nblks_tot-nblks_skip-n+lcm_s_t-1)/lcm_s_t
      aux_stride = nblk * nblks_comm
!      if(nblks_comm==0) cycle
      if (nblks_comm .ne. 0) then
        if (myps == ips) then
!          k = 0
          do lc=1,nvc
            do i = nblks_skip+n, nblks_tot-1, lcm_s_t
              k = (i - nblks_skip - n)/lcm_s_t * nblk + (lc - 1) * aux_stride
              ns = (i/nps)*nblk ! local start of block i
              nl = min(nvr-i*nblk,nblk) ! length
              aux(k+1:k+nl) = vmat_s(ns+1:ns+nl,lc)
!              k = k+nblk
            enddo
          enddo
        endif ! (myps == ips)



!        k = 0
        do lc=1,nvc
          do i = nblks_skip+n, nblks_tot-1, lcm_s_t
            k = (i - nblks_skip - n)/lcm_s_t * nblk + (lc - 1) * aux_stride
            ns = (i/npt)*nblk ! local start of block i
            nl = min(nvr-i*nblk,nblk) ! length
            vmat_t(ns+1:ns+nl,lc) = - aux(k+1:k+nl)
!            k = k+nblk
          enddo
        enddo
      endif
    endif ! (mypt == ipt)

  enddo
  deallocate(aux, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("elpa_transpose_vectors: aux", 374,  istat,  errorMessage)

  call obj%timer%stop("&
  &elpa_transpose_vectors_ss_&
  &real&
  &" // &
  &"_single" &
  )

end subroutine




!cannot use "../src/elpa1/../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)


subroutine elpa_reduce_add_vectors_&
&real&
&_&
&single &
(obj, vmat_s, ld_s, comm_s, vmat_t, ld_t, comm_t, nvr, nvc, nblk, nrThreads)

!-------------------------------------------------------------------------------
! This routine does a reduce of all vectors in vmat_s over the communicator comm_t.
! The result of the reduce is gathered on the processors owning the diagonal
! and added to the array of vectors vmat_t (which is distributed over comm_t).
!
! Opposed to elpa_transpose_vectors, there is NO identical copy of vmat_s
! in the different members within vmat_t (else a reduce wouldn't be necessary).
! After this routine, an allreduce of vmat_t has to be done.
!
! vmat_s    array of vectors to be reduced and added
! ld_s      leading dimension of vmat_s
! comm_s    communicator over which vmat_s is distributed
! vmat_t    array of vectors to which vmat_s is added
! ld_t      leading dimension of vmat_t
! comm_t    communicator over which vmat_t is distributed
! nvr       global length of vmat_s/vmat_t
! nvc       number of columns in vmat_s/vmat_t
! nblk      block size of block cyclic distribution
!
!-------------------------------------------------------------------------------

  use precision
  use elpa_mpi
  use elpa_abstract_impl
  implicit none

  class(elpa_abstract_impl_t), intent(inout)         :: obj
  integer(kind=ik), intent(in)                       :: ld_s, comm_s, ld_t, comm_t, nvr, nvc, nblk
  real(kind=c_float), intent(in)    :: vmat_s(ld_s,nvc)
  real(kind=c_float), intent(inout) :: vmat_t(ld_t,nvc)

  real(kind=c_float), allocatable   :: aux1(:), aux2(:)
  integer(kind=ik)                                   :: myps, mypt, nps, npt
  integer(kind=MPI_KIND)                             :: mypsMPI, npsMPI, myptMPI, nptMPI
  integer(kind=ik)                                   :: n, lc, k, i, ips, ipt, ns, nl
  integer(kind=MPI_KIND)                             :: mpierr
  integer(kind=ik)                                   :: lcm_s_t, nblks_tot
  integer(kind=ik)                                   :: aux_stride
  integer(kind=ik), intent(in)                       :: nrThreads
  integer(kind=ik)                                   :: istat
  character(200)                                     :: errorMessage

  call obj%timer%start("elpa_reduce_add_vectors_&
  &real&
  &" // &
  &"_single" &
  )

  call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(comm_s,kind=MPI_KIND), mypsMPI, mpierr)
  call mpi_comm_size(int(comm_s,kind=MPI_KIND), npsMPI,  mpierr)
  call mpi_comm_rank(int(comm_t,kind=MPI_KIND), myptMPI, mpierr)
  call mpi_comm_size(int(comm_t,kind=MPI_KIND), nptMPI ,mpierr)
  myps = int(mypsMPI,kind=c_int)
  nps = int(npsMPI,kind=c_int)
  mypt = int(myptMPI,kind=c_int)
  npt = int(nptMPI,kind=c_int)

  call obj%timer%stop("mpi_communication")

  ! Look to elpa_transpose_vectors for the basic idea!

  ! The communictation pattern repeats in the global matrix after
  ! the least common multiple of (nps,npt) blocks

  lcm_s_t   = least_common_multiple(nps,npt) ! least common multiple of nps, npt

  nblks_tot = (nvr+nblk-1)/nblk ! number of blocks corresponding to nvr

  allocate(aux1( ((nblks_tot+lcm_s_t-1)/lcm_s_t) * nblk * nvc ), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_reduce_add: aux1", 129,  istat,  errorMessage)

  allocate(aux2( ((nblks_tot+lcm_s_t-1)/lcm_s_t) * nblk * nvc ), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_reduce_add: aux2", 132,  istat,  errorMessage)
  aux1(:) = 0
  aux2(:) = 0
  do n = 0, lcm_s_t-1

    ips = mod(n,nps)
    ipt = mod(n,npt)

    aux_stride = nblk * ((nblks_tot - n + lcm_s_t - 1)/lcm_s_t)

    if (myps == ips) then

!      k = 0
      do lc=1,nvc
        do i = n, nblks_tot-1, lcm_s_t
          k = (i - n)/lcm_s_t * nblk + (lc - 1) * aux_stride
          ns = (i/nps)*nblk ! local start of block i
          nl = min(nvr-i*nblk,nblk) ! length
          aux1(k+1:k+nl) = vmat_s(ns+1:ns+nl,lc)
!          k = k+nblk
        enddo
      enddo

      k = nvc * aux_stride

      if (k>0) aux2 = aux1

      if (mypt == ipt) then
!        k = 0
        do lc=1,nvc
          do  i = n, nblks_tot-1, lcm_s_t
            k = (i - n)/lcm_s_t * nblk + (lc - 1) * aux_stride
            ns = (i/npt)*nblk ! local start of block i
            nl = min(nvr-i*nblk,nblk) ! length
            vmat_t(ns+1:ns+nl,lc) = vmat_t(ns+1:ns+nl,lc) + aux2(k+1:k+nl)
!            k = k+nblk
          enddo
        enddo
      endif ! (mypt == ipt)

    endif ! (myps == ips)

  enddo

  deallocate(aux1, aux2, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("elpa_reduce_add: aux1, aux2", 218,  istat,  errorMessage)

  call obj%timer%stop("elpa_reduce_add_vectors_&
  &real&
  &" // &
  &"_single" &
  )
end subroutine



! complex double precision






















!cannot use "../src/elpa1/../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)




subroutine elpa_transpose_vectors_&
&complex&
&_&
&double &
(obj, vmat_s, ld_s, comm_s, vmat_t, ld_t, comm_t, & 
 nvs, nvr, nvc, nblk, nrThreads, comm_s_isRows, success)

!-------------------------------------------------------------------------------
! This routine transposes an array of vectors which are distributed in
! communicator comm_s into its transposed form distributed in communicator comm_t.
! There must be an identical copy of vmat_s in every communicator comm_s.
! After this routine, there is an identical copy of vmat_t in every communicator comm_t.
!
! vmat_s    original array of vectors
! ld_s      leading dimension of vmat_s
! comm_s    communicator over which vmat_s is distributed
! vmat_t    array of vectors in transposed form
! ld_t      leading dimension of vmat_t
! comm_t    communicator over which vmat_t is distributed
! nvs       global index where to start in vmat_s/vmat_t
!           Please note: this is kind of a hint, some values before nvs will be
!           accessed in vmat_s/put into vmat_t
! nvr       global length of vmat_s/vmat_t
! nvc       number of columns in vmat_s/vmat_t
! nblk      block size of block cyclic distribution
!
!-------------------------------------------------------------------------------
  use precision
  use elpa_abstract_impl
  use elpa_mpi

  implicit none
  class(elpa_abstract_impl_t), intent(inout) :: obj
  integer(kind=ik), intent(in)                      :: ld_s, comm_s, ld_t, comm_t, nvs, nvr, nvc, nblk
  complex(kind=c_double), intent(in)   :: vmat_s(ld_s,nvc)
  complex(kind=c_double), intent(inout):: vmat_t(ld_t,nvc)

  complex(kind=c_double), allocatable  :: aux(:)
  integer(kind=ik)                                  :: myps, mypt, nps, npt
  integer(kind=MPI_KIND)                            :: mypsMPI, myptMPI, npsMPI, nptMPI
  integer(kind=ik)                                  :: n, lc, k, i, ips, ipt, ns, nl
  integer(kind=MPI_KIND)                            :: mpierr
  integer(kind=ik)                                  :: lcm_s_t, nblks_tot, nblks_comm, nblks_skip
  integer(kind=ik)                                  :: aux_stride
  integer(kind=ik), intent(in)                      :: nrThreads
  integer(kind=ik)                                  :: istat
  character(200)                                    :: errorMessage

  integer(kind=MPI_KIND)                            :: bcast_request1
  logical                                           :: useNonBlockingCollectives
  logical                                           :: useNonBlockingCollectivesRows
  logical                                           :: useNonBlockingCollectivesCols
  logical, intent(in)                               :: comm_s_isRows
  integer(kind=c_int)                               :: non_blocking_collectives_rows, error, &
                                                       non_blocking_collectives_cols
  logical                                           :: success

  integer(kind=ik)                                  :: mpi_comm_all, my_mpi_rank, transposed_mpi_rank, message_size, matrix_order
  integer(kind=ik)                                  :: ld_st, solver

  success = .true.

  call obj%timer%start("&
          &elpa_transpose_vectors_&
  &complex&
  &" // &
  &"_double" &
  )

  call obj%get("nbc_row_transpose_vectors", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    write(error_unit,*) "Problem setting option for non blocking collectives for_rows in transpose_vectors. Aborting..."
    success = .false.
     call obj%timer%stop("&
          &elpa_transpose_vectors_&
     &complex&
     &" // &
     &"_double" &
    )
    return
  endif

  call obj%get("nbc_col_transpose_vectors", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    write(error_unit,*) "Problem setting option for non blocking collectives for_cols in transpose_vectors. Aborting..."
    success = .false.
     call obj%timer%stop("&
          &elpa_transpose_vectors_&
     &complex&
     &" // &
     &"_double" &
    )
    return
  endif

  if (non_blocking_collectives_rows .eq. 1) then
    useNonBlockingCollectivesRows = .true.
  else
    useNonBlockingCollectivesRows = .false.
  endif

  if (non_blocking_collectives_cols .eq. 1) then
    useNonBlockingCollectivesCols = .true.
  else
    useNonBlockingCollectivesCols = .false.
  endif

  if (comm_s_isRows) then
    useNonBlockingCollectives = useNonBlockingCollectivesRows
  else
    useNonBlockingCollectives = useNonBlockingCollectivesCols
  endif

  call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(comm_s,kind=MPI_KIND),mypsMPI, mpierr)
  call mpi_comm_size(int(comm_s,kind=MPI_KIND),npsMPI ,mpierr)
  call mpi_comm_rank(int(comm_t,kind=MPI_KIND),myptMPI, mpierr)
  call mpi_comm_size(int(comm_t,kind=MPI_KIND),nptMPI ,mpierr)
  myps = int(mypsMPI,kind=c_int)
  nps = int(npsMPI,kind=c_int)
  mypt = int(myptMPI,kind=c_int)
  npt = int(nptMPI,kind=c_int)
  call obj%timer%stop("mpi_communication")

  ! TODO_23_11
  ! this codepath doesn't work for ELPA2, Cholesky, maybe smth else. Fix it, along with optimization of Cholesky-GPU
  ! because there nvc>1 and ld_s != ld_t (so, we can't make a contigous-memory MPI_Send call)

  ! The basic idea of this routine is that for every block (in the block cyclic
  ! distribution), the processor within comm_t which owns the diagonal
  ! broadcasts its values of vmat_s to all processors within comm_t.
  ! Of course this has not to be done for every block separately, since
  ! the communictation pattern repeats in the global matrix after
  ! the least common multiple of (nps,npt) blocks

  lcm_s_t   = least_common_multiple(nps,npt) ! least common multiple of nps, npt

  nblks_tot = (nvr+nblk-1)/nblk ! number of blocks corresponding to nvr

  ! Get the number of blocks to be skipped at the begin.
  ! This must be a multiple of lcm_s_t (else it is getting complicated),
  ! thus some elements before nvs will be accessed/set.

  nblks_skip = ((nvs-1)/(nblk*lcm_s_t))*lcm_s_t

  allocate(aux( ((nblks_tot-nblks_skip+lcm_s_t-1)/lcm_s_t) * nblk * nvc ), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_transpose_vectors: aux", 275,  istat,  errorMessage)
  do n = 0, lcm_s_t-1

    ips = mod(n,nps)
    ipt = mod(n,npt)

    if (mypt == ipt) then ! (mypt == ipt)

      nblks_comm = (nblks_tot-nblks_skip-n+lcm_s_t-1)/lcm_s_t
      aux_stride = nblk * nblks_comm
!      if(nblks_comm==0) cycle
      if (nblks_comm .ne. 0) then
        if (myps == ips) then
!          k = 0
          do lc=1,nvc
            do i = nblks_skip+n, nblks_tot-1, lcm_s_t
              k = (i - nblks_skip - n)/lcm_s_t * nblk + (lc - 1) * aux_stride
              ns = (i/nps)*nblk ! local start of block i
              nl = min(nvr-i*nblk,nblk) ! length
              aux(k+1:k+nl) = vmat_s(ns+1:ns+nl,lc)
!              k = k+nblk
            enddo
          enddo
        endif ! (myps == ips)



!        k = 0
        do lc=1,nvc
          do i = nblks_skip+n, nblks_tot-1, lcm_s_t
            k = (i - nblks_skip - n)/lcm_s_t * nblk + (lc - 1) * aux_stride
            ns = (i/npt)*nblk ! local start of block i
            nl = min(nvr-i*nblk,nblk) ! length
            vmat_t(ns+1:ns+nl,lc) = aux(k+1:k+nl)
!            k = k+nblk
          enddo
        enddo
      endif
    endif ! (mypt == ipt)

  enddo
  deallocate(aux, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("elpa_transpose_vectors: aux", 374,  istat,  errorMessage)

  call obj%timer%stop("&
  &elpa_transpose_vectors_&
  &complex&
  &" // &
  &"_double" &
  )

end subroutine




!cannot use "../src/elpa1/../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)




subroutine elpa_transpose_vectors_ss_&
&complex&
&_&
&double &
(obj, vmat_s, ld_s, comm_s, vmat_t, ld_t, comm_t, & 
 nvs, nvr, nvc, nblk, nrThreads, comm_s_isRows, success)

!-------------------------------------------------------------------------------
! This routine transposes an array of vectors which are distributed in
! communicator comm_s into its transposed form distributed in communicator comm_t.
! There must be an identical copy of vmat_s in every communicator comm_s.
! After this routine, there is an identical copy of vmat_t in every communicator comm_t.
!
! vmat_s    original array of vectors
! ld_s      leading dimension of vmat_s
! comm_s    communicator over which vmat_s is distributed
! vmat_t    array of vectors in transposed form
! ld_t      leading dimension of vmat_t
! comm_t    communicator over which vmat_t is distributed
! nvs       global index where to start in vmat_s/vmat_t
!           Please note: this is kind of a hint, some values before nvs will be
!           accessed in vmat_s/put into vmat_t
! nvr       global length of vmat_s/vmat_t
! nvc       number of columns in vmat_s/vmat_t
! nblk      block size of block cyclic distribution
!
!-------------------------------------------------------------------------------
  use precision
  use elpa_abstract_impl
  use elpa_mpi

  implicit none
  class(elpa_abstract_impl_t), intent(inout) :: obj
  integer(kind=ik), intent(in)                      :: ld_s, comm_s, ld_t, comm_t, nvs, nvr, nvc, nblk
  complex(kind=c_double), intent(in)   :: vmat_s(ld_s,nvc)
  complex(kind=c_double), intent(inout):: vmat_t(ld_t,nvc)

  complex(kind=c_double), allocatable  :: aux(:)
  integer(kind=ik)                                  :: myps, mypt, nps, npt
  integer(kind=MPI_KIND)                            :: mypsMPI, myptMPI, npsMPI, nptMPI
  integer(kind=ik)                                  :: n, lc, k, i, ips, ipt, ns, nl
  integer(kind=MPI_KIND)                            :: mpierr
  integer(kind=ik)                                  :: lcm_s_t, nblks_tot, nblks_comm, nblks_skip
  integer(kind=ik)                                  :: aux_stride
  integer(kind=ik), intent(in)                      :: nrThreads
  integer(kind=ik)                                  :: istat
  character(200)                                    :: errorMessage

  integer(kind=MPI_KIND)                            :: bcast_request1
  logical                                           :: useNonBlockingCollectives
  logical                                           :: useNonBlockingCollectivesRows
  logical                                           :: useNonBlockingCollectivesCols
  logical, intent(in)                               :: comm_s_isRows
  integer(kind=c_int)                               :: non_blocking_collectives_rows, error, &
                                                       non_blocking_collectives_cols
  logical                                           :: success

  integer(kind=ik)                                  :: mpi_comm_all, my_mpi_rank, transposed_mpi_rank, message_size, matrix_order
  integer(kind=ik)                                  :: ld_st, solver

  success = .true.

  call obj%timer%start("&
          &elpa_transpose_vectors_ss_&
  &complex&
  &" // &
  &"_double" &
  )

  call obj%get("nbc_row_transpose_vectors", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    write(error_unit,*) "Problem setting option for non blocking collectives for_rows in transpose_vectors. Aborting..."
    success = .false.
     call obj%timer%stop("&
          &elpa_transpose_vectors_ss_&
     &complex&
     &" // &
     &"_double" &
    )
    return
  endif

  call obj%get("nbc_col_transpose_vectors", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    write(error_unit,*) "Problem setting option for non blocking collectives for_cols in transpose_vectors. Aborting..."
    success = .false.
     call obj%timer%stop("&
          &elpa_transpose_vectors_ss_&
     &complex&
     &" // &
     &"_double" &
    )
    return
  endif

  if (non_blocking_collectives_rows .eq. 1) then
    useNonBlockingCollectivesRows = .true.
  else
    useNonBlockingCollectivesRows = .false.
  endif

  if (non_blocking_collectives_cols .eq. 1) then
    useNonBlockingCollectivesCols = .true.
  else
    useNonBlockingCollectivesCols = .false.
  endif

  if (comm_s_isRows) then
    useNonBlockingCollectives = useNonBlockingCollectivesRows
  else
    useNonBlockingCollectives = useNonBlockingCollectivesCols
  endif

  call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(comm_s,kind=MPI_KIND),mypsMPI, mpierr)
  call mpi_comm_size(int(comm_s,kind=MPI_KIND),npsMPI ,mpierr)
  call mpi_comm_rank(int(comm_t,kind=MPI_KIND),myptMPI, mpierr)
  call mpi_comm_size(int(comm_t,kind=MPI_KIND),nptMPI ,mpierr)
  myps = int(mypsMPI,kind=c_int)
  nps = int(npsMPI,kind=c_int)
  mypt = int(myptMPI,kind=c_int)
  npt = int(nptMPI,kind=c_int)
  call obj%timer%stop("mpi_communication")

  ! TODO_23_11
  ! this codepath doesn't work for ELPA2, Cholesky, maybe smth else. Fix it, along with optimization of Cholesky-GPU
  ! because there nvc>1 and ld_s != ld_t (so, we can't make a contigous-memory MPI_Send call)

  ! The basic idea of this routine is that for every block (in the block cyclic
  ! distribution), the processor within comm_t which owns the diagonal
  ! broadcasts its values of vmat_s to all processors within comm_t.
  ! Of course this has not to be done for every block separately, since
  ! the communictation pattern repeats in the global matrix after
  ! the least common multiple of (nps,npt) blocks

  lcm_s_t   = least_common_multiple(nps,npt) ! least common multiple of nps, npt

  nblks_tot = (nvr+nblk-1)/nblk ! number of blocks corresponding to nvr

  ! Get the number of blocks to be skipped at the begin.
  ! This must be a multiple of lcm_s_t (else it is getting complicated),
  ! thus some elements before nvs will be accessed/set.

  nblks_skip = ((nvs-1)/(nblk*lcm_s_t))*lcm_s_t

  allocate(aux( ((nblks_tot-nblks_skip+lcm_s_t-1)/lcm_s_t) * nblk * nvc ), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_transpose_vectors: aux", 275,  istat,  errorMessage)
  do n = 0, lcm_s_t-1

    ips = mod(n,nps)
    ipt = mod(n,npt)

    if (mypt == ipt) then ! (mypt == ipt)

      nblks_comm = (nblks_tot-nblks_skip-n+lcm_s_t-1)/lcm_s_t
      aux_stride = nblk * nblks_comm
!      if(nblks_comm==0) cycle
      if (nblks_comm .ne. 0) then
        if (myps == ips) then
!          k = 0
          do lc=1,nvc
            do i = nblks_skip+n, nblks_tot-1, lcm_s_t
              k = (i - nblks_skip - n)/lcm_s_t * nblk + (lc - 1) * aux_stride
              ns = (i/nps)*nblk ! local start of block i
              nl = min(nvr-i*nblk,nblk) ! length
              aux(k+1:k+nl) = vmat_s(ns+1:ns+nl,lc)
!              k = k+nblk
            enddo
          enddo
        endif ! (myps == ips)



!        k = 0
        do lc=1,nvc
          do i = nblks_skip+n, nblks_tot-1, lcm_s_t
            k = (i - nblks_skip - n)/lcm_s_t * nblk + (lc - 1) * aux_stride
            ns = (i/npt)*nblk ! local start of block i
            nl = min(nvr-i*nblk,nblk) ! length
            vmat_t(ns+1:ns+nl,lc) = - aux(k+1:k+nl)
!            k = k+nblk
          enddo
        enddo
      endif
    endif ! (mypt == ipt)

  enddo
  deallocate(aux, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("elpa_transpose_vectors: aux", 374,  istat,  errorMessage)

  call obj%timer%stop("&
  &elpa_transpose_vectors_ss_&
  &complex&
  &" // &
  &"_double" &
  )

end subroutine




!cannot use "../src/elpa1/../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)


subroutine elpa_reduce_add_vectors_&
&complex&
&_&
&double &
(obj, vmat_s, ld_s, comm_s, vmat_t, ld_t, comm_t, nvr, nvc, nblk, nrThreads)

!-------------------------------------------------------------------------------
! This routine does a reduce of all vectors in vmat_s over the communicator comm_t.
! The result of the reduce is gathered on the processors owning the diagonal
! and added to the array of vectors vmat_t (which is distributed over comm_t).
!
! Opposed to elpa_transpose_vectors, there is NO identical copy of vmat_s
! in the different members within vmat_t (else a reduce wouldn't be necessary).
! After this routine, an allreduce of vmat_t has to be done.
!
! vmat_s    array of vectors to be reduced and added
! ld_s      leading dimension of vmat_s
! comm_s    communicator over which vmat_s is distributed
! vmat_t    array of vectors to which vmat_s is added
! ld_t      leading dimension of vmat_t
! comm_t    communicator over which vmat_t is distributed
! nvr       global length of vmat_s/vmat_t
! nvc       number of columns in vmat_s/vmat_t
! nblk      block size of block cyclic distribution
!
!-------------------------------------------------------------------------------

  use precision
  use elpa_mpi
  use elpa_abstract_impl
  implicit none

  class(elpa_abstract_impl_t), intent(inout)         :: obj
  integer(kind=ik), intent(in)                       :: ld_s, comm_s, ld_t, comm_t, nvr, nvc, nblk
  complex(kind=c_double), intent(in)    :: vmat_s(ld_s,nvc)
  complex(kind=c_double), intent(inout) :: vmat_t(ld_t,nvc)

  complex(kind=c_double), allocatable   :: aux1(:), aux2(:)
  integer(kind=ik)                                   :: myps, mypt, nps, npt
  integer(kind=MPI_KIND)                             :: mypsMPI, npsMPI, myptMPI, nptMPI
  integer(kind=ik)                                   :: n, lc, k, i, ips, ipt, ns, nl
  integer(kind=MPI_KIND)                             :: mpierr
  integer(kind=ik)                                   :: lcm_s_t, nblks_tot
  integer(kind=ik)                                   :: aux_stride
  integer(kind=ik), intent(in)                       :: nrThreads
  integer(kind=ik)                                   :: istat
  character(200)                                     :: errorMessage

  call obj%timer%start("elpa_reduce_add_vectors_&
  &complex&
  &" // &
  &"_double" &
  )

  call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(comm_s,kind=MPI_KIND), mypsMPI, mpierr)
  call mpi_comm_size(int(comm_s,kind=MPI_KIND), npsMPI,  mpierr)
  call mpi_comm_rank(int(comm_t,kind=MPI_KIND), myptMPI, mpierr)
  call mpi_comm_size(int(comm_t,kind=MPI_KIND), nptMPI ,mpierr)
  myps = int(mypsMPI,kind=c_int)
  nps = int(npsMPI,kind=c_int)
  mypt = int(myptMPI,kind=c_int)
  npt = int(nptMPI,kind=c_int)

  call obj%timer%stop("mpi_communication")

  ! Look to elpa_transpose_vectors for the basic idea!

  ! The communictation pattern repeats in the global matrix after
  ! the least common multiple of (nps,npt) blocks

  lcm_s_t   = least_common_multiple(nps,npt) ! least common multiple of nps, npt

  nblks_tot = (nvr+nblk-1)/nblk ! number of blocks corresponding to nvr

  allocate(aux1( ((nblks_tot+lcm_s_t-1)/lcm_s_t) * nblk * nvc ), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_reduce_add: aux1", 129,  istat,  errorMessage)

  allocate(aux2( ((nblks_tot+lcm_s_t-1)/lcm_s_t) * nblk * nvc ), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_reduce_add: aux2", 132,  istat,  errorMessage)
  aux1(:) = 0
  aux2(:) = 0
  do n = 0, lcm_s_t-1

    ips = mod(n,nps)
    ipt = mod(n,npt)

    aux_stride = nblk * ((nblks_tot - n + lcm_s_t - 1)/lcm_s_t)

    if (myps == ips) then

!      k = 0
      do lc=1,nvc
        do i = n, nblks_tot-1, lcm_s_t
          k = (i - n)/lcm_s_t * nblk + (lc - 1) * aux_stride
          ns = (i/nps)*nblk ! local start of block i
          nl = min(nvr-i*nblk,nblk) ! length
          aux1(k+1:k+nl) = vmat_s(ns+1:ns+nl,lc)
!          k = k+nblk
        enddo
      enddo

      k = nvc * aux_stride

      if (k>0) aux2 = aux1

      if (mypt == ipt) then
!        k = 0
        do lc=1,nvc
          do  i = n, nblks_tot-1, lcm_s_t
            k = (i - n)/lcm_s_t * nblk + (lc - 1) * aux_stride
            ns = (i/npt)*nblk ! local start of block i
            nl = min(nvr-i*nblk,nblk) ! length
            vmat_t(ns+1:ns+nl,lc) = vmat_t(ns+1:ns+nl,lc) + aux2(k+1:k+nl)
!            k = k+nblk
          enddo
        enddo
      endif ! (mypt == ipt)

    endif ! (myps == ips)

  enddo

  deallocate(aux1, aux2, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("elpa_reduce_add: aux1, aux2", 218,  istat,  errorMessage)

  call obj%timer%stop("elpa_reduce_add_vectors_&
  &complex&
  &" // &
  &"_double" &
  )
end subroutine



! complex single precision





















!cannot use "../src/elpa1/../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)




subroutine elpa_transpose_vectors_&
&complex&
&_&
&single &
(obj, vmat_s, ld_s, comm_s, vmat_t, ld_t, comm_t, & 
 nvs, nvr, nvc, nblk, nrThreads, comm_s_isRows, success)

!-------------------------------------------------------------------------------
! This routine transposes an array of vectors which are distributed in
! communicator comm_s into its transposed form distributed in communicator comm_t.
! There must be an identical copy of vmat_s in every communicator comm_s.
! After this routine, there is an identical copy of vmat_t in every communicator comm_t.
!
! vmat_s    original array of vectors
! ld_s      leading dimension of vmat_s
! comm_s    communicator over which vmat_s is distributed
! vmat_t    array of vectors in transposed form
! ld_t      leading dimension of vmat_t
! comm_t    communicator over which vmat_t is distributed
! nvs       global index where to start in vmat_s/vmat_t
!           Please note: this is kind of a hint, some values before nvs will be
!           accessed in vmat_s/put into vmat_t
! nvr       global length of vmat_s/vmat_t
! nvc       number of columns in vmat_s/vmat_t
! nblk      block size of block cyclic distribution
!
!-------------------------------------------------------------------------------
  use precision
  use elpa_abstract_impl
  use elpa_mpi

  implicit none
  class(elpa_abstract_impl_t), intent(inout) :: obj
  integer(kind=ik), intent(in)                      :: ld_s, comm_s, ld_t, comm_t, nvs, nvr, nvc, nblk
  complex(kind=c_float), intent(in)   :: vmat_s(ld_s,nvc)
  complex(kind=c_float), intent(inout):: vmat_t(ld_t,nvc)

  complex(kind=c_float), allocatable  :: aux(:)
  integer(kind=ik)                                  :: myps, mypt, nps, npt
  integer(kind=MPI_KIND)                            :: mypsMPI, myptMPI, npsMPI, nptMPI
  integer(kind=ik)                                  :: n, lc, k, i, ips, ipt, ns, nl
  integer(kind=MPI_KIND)                            :: mpierr
  integer(kind=ik)                                  :: lcm_s_t, nblks_tot, nblks_comm, nblks_skip
  integer(kind=ik)                                  :: aux_stride
  integer(kind=ik), intent(in)                      :: nrThreads
  integer(kind=ik)                                  :: istat
  character(200)                                    :: errorMessage

  integer(kind=MPI_KIND)                            :: bcast_request1
  logical                                           :: useNonBlockingCollectives
  logical                                           :: useNonBlockingCollectivesRows
  logical                                           :: useNonBlockingCollectivesCols
  logical, intent(in)                               :: comm_s_isRows
  integer(kind=c_int)                               :: non_blocking_collectives_rows, error, &
                                                       non_blocking_collectives_cols
  logical                                           :: success

  integer(kind=ik)                                  :: mpi_comm_all, my_mpi_rank, transposed_mpi_rank, message_size, matrix_order
  integer(kind=ik)                                  :: ld_st, solver

  success = .true.

  call obj%timer%start("&
          &elpa_transpose_vectors_&
  &complex&
  &" // &
  &"_single" &
  )

  call obj%get("nbc_row_transpose_vectors", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    write(error_unit,*) "Problem setting option for non blocking collectives for_rows in transpose_vectors. Aborting..."
    success = .false.
     call obj%timer%stop("&
          &elpa_transpose_vectors_&
     &complex&
     &" // &
     &"_single" &
    )
    return
  endif

  call obj%get("nbc_col_transpose_vectors", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    write(error_unit,*) "Problem setting option for non blocking collectives for_cols in transpose_vectors. Aborting..."
    success = .false.
     call obj%timer%stop("&
          &elpa_transpose_vectors_&
     &complex&
     &" // &
     &"_single" &
    )
    return
  endif

  if (non_blocking_collectives_rows .eq. 1) then
    useNonBlockingCollectivesRows = .true.
  else
    useNonBlockingCollectivesRows = .false.
  endif

  if (non_blocking_collectives_cols .eq. 1) then
    useNonBlockingCollectivesCols = .true.
  else
    useNonBlockingCollectivesCols = .false.
  endif

  if (comm_s_isRows) then
    useNonBlockingCollectives = useNonBlockingCollectivesRows
  else
    useNonBlockingCollectives = useNonBlockingCollectivesCols
  endif

  call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(comm_s,kind=MPI_KIND),mypsMPI, mpierr)
  call mpi_comm_size(int(comm_s,kind=MPI_KIND),npsMPI ,mpierr)
  call mpi_comm_rank(int(comm_t,kind=MPI_KIND),myptMPI, mpierr)
  call mpi_comm_size(int(comm_t,kind=MPI_KIND),nptMPI ,mpierr)
  myps = int(mypsMPI,kind=c_int)
  nps = int(npsMPI,kind=c_int)
  mypt = int(myptMPI,kind=c_int)
  npt = int(nptMPI,kind=c_int)
  call obj%timer%stop("mpi_communication")

  ! TODO_23_11
  ! this codepath doesn't work for ELPA2, Cholesky, maybe smth else. Fix it, along with optimization of Cholesky-GPU
  ! because there nvc>1 and ld_s != ld_t (so, we can't make a contigous-memory MPI_Send call)

  ! The basic idea of this routine is that for every block (in the block cyclic
  ! distribution), the processor within comm_t which owns the diagonal
  ! broadcasts its values of vmat_s to all processors within comm_t.
  ! Of course this has not to be done for every block separately, since
  ! the communictation pattern repeats in the global matrix after
  ! the least common multiple of (nps,npt) blocks

  lcm_s_t   = least_common_multiple(nps,npt) ! least common multiple of nps, npt

  nblks_tot = (nvr+nblk-1)/nblk ! number of blocks corresponding to nvr

  ! Get the number of blocks to be skipped at the begin.
  ! This must be a multiple of lcm_s_t (else it is getting complicated),
  ! thus some elements before nvs will be accessed/set.

  nblks_skip = ((nvs-1)/(nblk*lcm_s_t))*lcm_s_t

  allocate(aux( ((nblks_tot-nblks_skip+lcm_s_t-1)/lcm_s_t) * nblk * nvc ), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_transpose_vectors: aux", 275,  istat,  errorMessage)
  do n = 0, lcm_s_t-1

    ips = mod(n,nps)
    ipt = mod(n,npt)

    if (mypt == ipt) then ! (mypt == ipt)

      nblks_comm = (nblks_tot-nblks_skip-n+lcm_s_t-1)/lcm_s_t
      aux_stride = nblk * nblks_comm
!      if(nblks_comm==0) cycle
      if (nblks_comm .ne. 0) then
        if (myps == ips) then
!          k = 0
          do lc=1,nvc
            do i = nblks_skip+n, nblks_tot-1, lcm_s_t
              k = (i - nblks_skip - n)/lcm_s_t * nblk + (lc - 1) * aux_stride
              ns = (i/nps)*nblk ! local start of block i
              nl = min(nvr-i*nblk,nblk) ! length
              aux(k+1:k+nl) = vmat_s(ns+1:ns+nl,lc)
!              k = k+nblk
            enddo
          enddo
        endif ! (myps == ips)



!        k = 0
        do lc=1,nvc
          do i = nblks_skip+n, nblks_tot-1, lcm_s_t
            k = (i - nblks_skip - n)/lcm_s_t * nblk + (lc - 1) * aux_stride
            ns = (i/npt)*nblk ! local start of block i
            nl = min(nvr-i*nblk,nblk) ! length
            vmat_t(ns+1:ns+nl,lc) = aux(k+1:k+nl)
!            k = k+nblk
          enddo
        enddo
      endif
    endif ! (mypt == ipt)

  enddo
  deallocate(aux, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("elpa_transpose_vectors: aux", 374,  istat,  errorMessage)

  call obj%timer%stop("&
  &elpa_transpose_vectors_&
  &complex&
  &" // &
  &"_single" &
  )

end subroutine




!cannot use "../src/elpa1/../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)




subroutine elpa_transpose_vectors_ss_&
&complex&
&_&
&single &
(obj, vmat_s, ld_s, comm_s, vmat_t, ld_t, comm_t, & 
 nvs, nvr, nvc, nblk, nrThreads, comm_s_isRows, success)

!-------------------------------------------------------------------------------
! This routine transposes an array of vectors which are distributed in
! communicator comm_s into its transposed form distributed in communicator comm_t.
! There must be an identical copy of vmat_s in every communicator comm_s.
! After this routine, there is an identical copy of vmat_t in every communicator comm_t.
!
! vmat_s    original array of vectors
! ld_s      leading dimension of vmat_s
! comm_s    communicator over which vmat_s is distributed
! vmat_t    array of vectors in transposed form
! ld_t      leading dimension of vmat_t
! comm_t    communicator over which vmat_t is distributed
! nvs       global index where to start in vmat_s/vmat_t
!           Please note: this is kind of a hint, some values before nvs will be
!           accessed in vmat_s/put into vmat_t
! nvr       global length of vmat_s/vmat_t
! nvc       number of columns in vmat_s/vmat_t
! nblk      block size of block cyclic distribution
!
!-------------------------------------------------------------------------------
  use precision
  use elpa_abstract_impl
  use elpa_mpi

  implicit none
  class(elpa_abstract_impl_t), intent(inout) :: obj
  integer(kind=ik), intent(in)                      :: ld_s, comm_s, ld_t, comm_t, nvs, nvr, nvc, nblk
  complex(kind=c_float), intent(in)   :: vmat_s(ld_s,nvc)
  complex(kind=c_float), intent(inout):: vmat_t(ld_t,nvc)

  complex(kind=c_float), allocatable  :: aux(:)
  integer(kind=ik)                                  :: myps, mypt, nps, npt
  integer(kind=MPI_KIND)                            :: mypsMPI, myptMPI, npsMPI, nptMPI
  integer(kind=ik)                                  :: n, lc, k, i, ips, ipt, ns, nl
  integer(kind=MPI_KIND)                            :: mpierr
  integer(kind=ik)                                  :: lcm_s_t, nblks_tot, nblks_comm, nblks_skip
  integer(kind=ik)                                  :: aux_stride
  integer(kind=ik), intent(in)                      :: nrThreads
  integer(kind=ik)                                  :: istat
  character(200)                                    :: errorMessage

  integer(kind=MPI_KIND)                            :: bcast_request1
  logical                                           :: useNonBlockingCollectives
  logical                                           :: useNonBlockingCollectivesRows
  logical                                           :: useNonBlockingCollectivesCols
  logical, intent(in)                               :: comm_s_isRows
  integer(kind=c_int)                               :: non_blocking_collectives_rows, error, &
                                                       non_blocking_collectives_cols
  logical                                           :: success

  integer(kind=ik)                                  :: mpi_comm_all, my_mpi_rank, transposed_mpi_rank, message_size, matrix_order
  integer(kind=ik)                                  :: ld_st, solver

  success = .true.

  call obj%timer%start("&
          &elpa_transpose_vectors_ss_&
  &complex&
  &" // &
  &"_single" &
  )

  call obj%get("nbc_row_transpose_vectors", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    write(error_unit,*) "Problem setting option for non blocking collectives for_rows in transpose_vectors. Aborting..."
    success = .false.
     call obj%timer%stop("&
          &elpa_transpose_vectors_ss_&
     &complex&
     &" // &
     &"_single" &
    )
    return
  endif

  call obj%get("nbc_col_transpose_vectors", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    write(error_unit,*) "Problem setting option for non blocking collectives for_cols in transpose_vectors. Aborting..."
    success = .false.
     call obj%timer%stop("&
          &elpa_transpose_vectors_ss_&
     &complex&
     &" // &
     &"_single" &
    )
    return
  endif

  if (non_blocking_collectives_rows .eq. 1) then
    useNonBlockingCollectivesRows = .true.
  else
    useNonBlockingCollectivesRows = .false.
  endif

  if (non_blocking_collectives_cols .eq. 1) then
    useNonBlockingCollectivesCols = .true.
  else
    useNonBlockingCollectivesCols = .false.
  endif

  if (comm_s_isRows) then
    useNonBlockingCollectives = useNonBlockingCollectivesRows
  else
    useNonBlockingCollectives = useNonBlockingCollectivesCols
  endif

  call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(comm_s,kind=MPI_KIND),mypsMPI, mpierr)
  call mpi_comm_size(int(comm_s,kind=MPI_KIND),npsMPI ,mpierr)
  call mpi_comm_rank(int(comm_t,kind=MPI_KIND),myptMPI, mpierr)
  call mpi_comm_size(int(comm_t,kind=MPI_KIND),nptMPI ,mpierr)
  myps = int(mypsMPI,kind=c_int)
  nps = int(npsMPI,kind=c_int)
  mypt = int(myptMPI,kind=c_int)
  npt = int(nptMPI,kind=c_int)
  call obj%timer%stop("mpi_communication")

  ! TODO_23_11
  ! this codepath doesn't work for ELPA2, Cholesky, maybe smth else. Fix it, along with optimization of Cholesky-GPU
  ! because there nvc>1 and ld_s != ld_t (so, we can't make a contigous-memory MPI_Send call)

  ! The basic idea of this routine is that for every block (in the block cyclic
  ! distribution), the processor within comm_t which owns the diagonal
  ! broadcasts its values of vmat_s to all processors within comm_t.
  ! Of course this has not to be done for every block separately, since
  ! the communictation pattern repeats in the global matrix after
  ! the least common multiple of (nps,npt) blocks

  lcm_s_t   = least_common_multiple(nps,npt) ! least common multiple of nps, npt

  nblks_tot = (nvr+nblk-1)/nblk ! number of blocks corresponding to nvr

  ! Get the number of blocks to be skipped at the begin.
  ! This must be a multiple of lcm_s_t (else it is getting complicated),
  ! thus some elements before nvs will be accessed/set.

  nblks_skip = ((nvs-1)/(nblk*lcm_s_t))*lcm_s_t

  allocate(aux( ((nblks_tot-nblks_skip+lcm_s_t-1)/lcm_s_t) * nblk * nvc ), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_transpose_vectors: aux", 275,  istat,  errorMessage)
  do n = 0, lcm_s_t-1

    ips = mod(n,nps)
    ipt = mod(n,npt)

    if (mypt == ipt) then ! (mypt == ipt)

      nblks_comm = (nblks_tot-nblks_skip-n+lcm_s_t-1)/lcm_s_t
      aux_stride = nblk * nblks_comm
!      if(nblks_comm==0) cycle
      if (nblks_comm .ne. 0) then
        if (myps == ips) then
!          k = 0
          do lc=1,nvc
            do i = nblks_skip+n, nblks_tot-1, lcm_s_t
              k = (i - nblks_skip - n)/lcm_s_t * nblk + (lc - 1) * aux_stride
              ns = (i/nps)*nblk ! local start of block i
              nl = min(nvr-i*nblk,nblk) ! length
              aux(k+1:k+nl) = vmat_s(ns+1:ns+nl,lc)
!              k = k+nblk
            enddo
          enddo
        endif ! (myps == ips)



!        k = 0
        do lc=1,nvc
          do i = nblks_skip+n, nblks_tot-1, lcm_s_t
            k = (i - nblks_skip - n)/lcm_s_t * nblk + (lc - 1) * aux_stride
            ns = (i/npt)*nblk ! local start of block i
            nl = min(nvr-i*nblk,nblk) ! length
            vmat_t(ns+1:ns+nl,lc) = - aux(k+1:k+nl)
!            k = k+nblk
          enddo
        enddo
      endif
    endif ! (mypt == ipt)

  enddo
  deallocate(aux, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("elpa_transpose_vectors: aux", 374,  istat,  errorMessage)

  call obj%timer%stop("&
  &elpa_transpose_vectors_ss_&
  &complex&
  &" // &
  &"_single" &
  )

end subroutine




!cannot use "../src/elpa1/../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)


subroutine elpa_reduce_add_vectors_&
&complex&
&_&
&single &
(obj, vmat_s, ld_s, comm_s, vmat_t, ld_t, comm_t, nvr, nvc, nblk, nrThreads)

!-------------------------------------------------------------------------------
! This routine does a reduce of all vectors in vmat_s over the communicator comm_t.
! The result of the reduce is gathered on the processors owning the diagonal
! and added to the array of vectors vmat_t (which is distributed over comm_t).
!
! Opposed to elpa_transpose_vectors, there is NO identical copy of vmat_s
! in the different members within vmat_t (else a reduce wouldn't be necessary).
! After this routine, an allreduce of vmat_t has to be done.
!
! vmat_s    array of vectors to be reduced and added
! ld_s      leading dimension of vmat_s
! comm_s    communicator over which vmat_s is distributed
! vmat_t    array of vectors to which vmat_s is added
! ld_t      leading dimension of vmat_t
! comm_t    communicator over which vmat_t is distributed
! nvr       global length of vmat_s/vmat_t
! nvc       number of columns in vmat_s/vmat_t
! nblk      block size of block cyclic distribution
!
!-------------------------------------------------------------------------------

  use precision
  use elpa_mpi
  use elpa_abstract_impl
  implicit none

  class(elpa_abstract_impl_t), intent(inout)         :: obj
  integer(kind=ik), intent(in)                       :: ld_s, comm_s, ld_t, comm_t, nvr, nvc, nblk
  complex(kind=c_float), intent(in)    :: vmat_s(ld_s,nvc)
  complex(kind=c_float), intent(inout) :: vmat_t(ld_t,nvc)

  complex(kind=c_float), allocatable   :: aux1(:), aux2(:)
  integer(kind=ik)                                   :: myps, mypt, nps, npt
  integer(kind=MPI_KIND)                             :: mypsMPI, npsMPI, myptMPI, nptMPI
  integer(kind=ik)                                   :: n, lc, k, i, ips, ipt, ns, nl
  integer(kind=MPI_KIND)                             :: mpierr
  integer(kind=ik)                                   :: lcm_s_t, nblks_tot
  integer(kind=ik)                                   :: aux_stride
  integer(kind=ik), intent(in)                       :: nrThreads
  integer(kind=ik)                                   :: istat
  character(200)                                     :: errorMessage

  call obj%timer%start("elpa_reduce_add_vectors_&
  &complex&
  &" // &
  &"_single" &
  )

  call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(comm_s,kind=MPI_KIND), mypsMPI, mpierr)
  call mpi_comm_size(int(comm_s,kind=MPI_KIND), npsMPI,  mpierr)
  call mpi_comm_rank(int(comm_t,kind=MPI_KIND), myptMPI, mpierr)
  call mpi_comm_size(int(comm_t,kind=MPI_KIND), nptMPI ,mpierr)
  myps = int(mypsMPI,kind=c_int)
  nps = int(npsMPI,kind=c_int)
  mypt = int(myptMPI,kind=c_int)
  npt = int(nptMPI,kind=c_int)

  call obj%timer%stop("mpi_communication")

  ! Look to elpa_transpose_vectors for the basic idea!

  ! The communictation pattern repeats in the global matrix after
  ! the least common multiple of (nps,npt) blocks

  lcm_s_t   = least_common_multiple(nps,npt) ! least common multiple of nps, npt

  nblks_tot = (nvr+nblk-1)/nblk ! number of blocks corresponding to nvr

  allocate(aux1( ((nblks_tot+lcm_s_t-1)/lcm_s_t) * nblk * nvc ), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_reduce_add: aux1", 129,  istat,  errorMessage)

  allocate(aux2( ((nblks_tot+lcm_s_t-1)/lcm_s_t) * nblk * nvc ), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_reduce_add: aux2", 132,  istat,  errorMessage)
  aux1(:) = 0
  aux2(:) = 0
  do n = 0, lcm_s_t-1

    ips = mod(n,nps)
    ipt = mod(n,npt)

    aux_stride = nblk * ((nblks_tot - n + lcm_s_t - 1)/lcm_s_t)

    if (myps == ips) then

!      k = 0
      do lc=1,nvc
        do i = n, nblks_tot-1, lcm_s_t
          k = (i - n)/lcm_s_t * nblk + (lc - 1) * aux_stride
          ns = (i/nps)*nblk ! local start of block i
          nl = min(nvr-i*nblk,nblk) ! length
          aux1(k+1:k+nl) = vmat_s(ns+1:ns+nl,lc)
!          k = k+nblk
        enddo
      enddo

      k = nvc * aux_stride

      if (k>0) aux2 = aux1

      if (mypt == ipt) then
!        k = 0
        do lc=1,nvc
          do  i = n, nblks_tot-1, lcm_s_t
            k = (i - n)/lcm_s_t * nblk + (lc - 1) * aux_stride
            ns = (i/npt)*nblk ! local start of block i
            nl = min(nvr-i*nblk,nblk) ! length
            vmat_t(ns+1:ns+nl,lc) = vmat_t(ns+1:ns+nl,lc) + aux2(k+1:k+nl)
!            k = k+nblk
          enddo
        enddo
      endif ! (mypt == ipt)

    endif ! (myps == ips)

  enddo

  deallocate(aux1, aux2, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("elpa_reduce_add: aux1, aux2", 218,  istat,  errorMessage)

  call obj%timer%stop("elpa_reduce_add_vectors_&
  &complex&
  &" // &
  &"_single" &
  )
end subroutine



!________________________________________________________________
! ./GPU/elpa_gpu_ccl_transpose_vectors_template.F90
! ./GPU/elpa_gpu_ccl_reduce_add_vectors_template.F90


!________________________________________________________________
! elpa1_compute_template.F90

! real double precision





















!cannot use "../src/elpa1/../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)










!> \brief Reduces a distributed symmetric matrix to tridiagonal form (like Scalapack Routine PDSYTRD)
!>
!  Parameters
!
!> \param obj	      object of elpa_type
!> \param na          Order of matrix
!>
!> \param a_mat(matrixRows,matrixCols)    Distributed matrix which should be reduced.
!>              Distribution is like in Scalapack.
!>              Opposed to PDSYTRD, a(:,:) must be set completely (upper and lower half)
!>              a(:,:) is overwritten on exit with the Householder vectors
!>
!> \param matrixRows         Leading dimension of a
!>
!> \param nblk        blocksize of cyclic distribution, must be the same in both directions!
!>
!> \param matrixCols  local columns of matrix
!>
!> \param mpi_comm_rows        MPI-Communicator for rows
!> \param mpi_comm_cols        MPI-Communicator for columns
!>
!> \param d_vec(na)       Diagonal elements (returned), identical on all processors
!>
!> \param e_vec(na)       Off-Diagonal elements (returned), identical on all processors
!>
!> \param tau(na)     Factors for the Householder vectors (returned), needed for back transformation
!>
!> \param useGPU      If true,  GPU version of the subroutine will be used
!> \param wantDebug   if true more debug information
!>
subroutine tridiag_&
  &real&
  &_&
  &double &
  (obj, na, a_mat, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols, d_vec, e_vec, tau, useGPU, wantDebug, &
   max_threads_in, isSkewsymmetric, success)
  use, intrinsic :: iso_c_binding
  use precision
  use elpa_abstract_impl
  use matrix_plot
  use elpa_omp
  use elpa_blas_interfaces
  use elpa_gpu
  use elpa_gpu_util
  use tridiag_gpu

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

  class(elpa_abstract_impl_t), intent(inout)    :: obj
  integer(kind=ik), intent(in)                  :: na, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols
  logical, intent(in)                           :: useGPU, wantDebug
  logical, intent(in)                           :: isSkewsymmetric

  logical                                       :: useCCL=.false.

  real(kind=rck), intent(out)          :: tau(na)
  real(kind=rck), intent(inout)        :: a_mat(matrixRows,*)
  real(kind=rk), intent(out)                    :: d_vec(na)
  real(kind=rk), intent(out)                    :: e_vec(na)
  integer(kind=ik)                              :: max_stored_uv = 32 ! TODO_23_11 - make it tunable instead of hard-coded
  logical,          parameter                   :: mat_vec_as_one_block = .true.

  ! id in processor row and column and total numbers of processor rows and columns
  integer(kind=ik)                              :: my_prow, my_pcol, np_rows, np_cols
  integer(kind=MPI_KIND)                        :: my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI
  integer(kind=MPI_KIND)                        :: mpierr
  integer(kind=ik)                              :: totalblocks, max_loc_block_rows, max_loc_block_cols, max_local_rows, &
                                                   max_local_cols
  ! updated after each istep (in the main cycle) to contain number of
  ! local columns and rows of the remaining part of the matrix
  integer(kind=ik)                              :: l_cols, l_rows
  integer(kind=ik)                              :: n_stored_vecs
  integer(kind=ik)                              :: isOurProcessRowInt ! TODO_23_11 - get rid of it
  logical                                       :: isOurProcessRow, isOurProcessCol, isOurProcessCol_prev


  integer(kind=C_intptr_T)                      :: a_dev, v_row_dev, v_col_dev, u_row_dev, u_col_dev, vu_stored_rows_dev, &
                                                   uv_stored_cols_dev, d_vec_dev, e_vec_dev, tau_dev
  logical                                       :: successGPU

  integer(kind=ik)                              :: istep, i, j, l_col_beg, l_col_end, l_row_beg, l_row_end
  integer(kind=ik)                              :: tile_size, l_rows_per_tile, l_cols_per_tile
  integer(kind=c_intptr_t)                      :: offset_dev

  integer(kind=ik), intent(in)                  :: max_threads_in
  integer(kind=ik)                              :: max_threads

  real(kind=rk)                                 :: vnorm2
  real(kind=rck)                       :: vav, x, aux1(2), vrl, xf, conjg_tau, dot_prod
  real(kind=rck), allocatable          :: aux(:)

  integer(kind=c_intptr_t)                      :: aux_dev, aux1_dev, aux_complex_dev, vav_dev, vav_host_or_dev, dot_prod_dev, & 
                                                   xf_dev, a_updated_element_dev, xf_host_or_dev, tau_istep_host_or_dev
  integer(kind=c_intptr_t)                      :: vnorm2_dev, vrl_dev ! alias pointers for aux1_dev(1), aux1_dev(2)


  integer(kind=c_intptr_t)                      :: num
  real(kind=rck), pointer              :: v_row(:), & ! used to store calculated Householder Vector
                                                   v_col(:)   ! the same Vector, but transposed 
  real(kind=rck), pointer              :: u_row_debug(:), & ! used to store calculated Householder Vector
                                                   u_col_debug(:)   ! the same Vector, but transposed 
  real(kind=rck), pointer              :: u_col(:), u_row(:)

  ! the following two matrices store pairs of vectors v and u calculated in each step
  ! at most max_stored_uv Vector pairs are stored, than the matrix A_i is explicitli updated
  ! u and v are stored both in row and Vector forms
  ! pattern: v1,u1,v2,u2,v3,u3,....
  ! todo: It is little bit confusing, I think, that variables _row actually store columns and vice versa
  !real(kind=rck), pointer             :: vu_stored_rows(:,:)
  real(kind=rck), allocatable         :: vu_stored_rows(:,:)
  ! pattern: u1,v1,u2,v2,u3,v3,....
  real(kind=rck), allocatable         :: uv_stored_cols(:,:)


  type(c_ptr)                                   :: v_row_host, v_col_host
  type(c_ptr)                                   :: u_row_host, u_col_host
  !type(c_ptr)                                   :: vu_stored_rows_host, uv_stored_cols_host
  real(kind=rk), allocatable                    :: tmp_real(:)
  integer(kind=ik)                              :: min_tile_size, error
  integer(kind=ik)                              :: istat
  character(200)                                :: errorMessage
  character(20)                                 :: gpuString
  integer(kind=ik)                              :: nblockEnd
  integer(kind=c_intptr_t), parameter           :: size_of_datatype = size_of_&
                                                                      &double&
                                                                      &_&
                                                                      &real
  integer(kind=c_intptr_t), parameter           :: size_of_datatype_real = size_of_&
                                                                      &double&
                                                                      &_real
  integer(kind=MPI_KIND)                        :: bcast_request1, bcast_request2, bcast_request3
  integer(kind=MPI_KIND)                        :: allreduce_request1, allreduce_request2, allreduce_request3
  integer(kind=MPI_KIND)                        :: allreduce_request4, allreduce_request5, allreduce_request6, &
                                                   allreduce_request7
  logical                                       :: useNonBlockingCollectivesCols
  logical                                       :: useNonBlockingCollectivesRows
  integer(kind=c_int)                           :: non_blocking_collectives_rows, non_blocking_collectives_cols
  logical                                       :: success

  integer(kind=c_intptr_t)                      :: gpuHandle, my_stream
  integer(kind=c_int) :: pointerMode

  integer(kind=ik)                              :: string_length, sm_count


  allocate(aux(2*max_stored_uv), stat=istat, errmsg=errorMessage)

  success = .true.

  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  if (useGPU) then
    max_threads=1
  else
    max_threads=max_threads_in
  endif

  call obj%timer%start("tridiag_&
  &real&
  &" // &
  "_double" // &
  gpuString )

  call obj%get("nbc_row_elpa1_full_to_tridi", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    write(error_unit,*) "Problem setting option for non blocking collectives for rows in elpa1_tridiag. Aborting..."
    success = .false.
    call obj%timer%stop("tridiag_&
    &real&
    &" // &
    "_double" // &
    gpuString )
    return
  endif

  call obj%get("nbc_col_elpa1_full_to_tridi", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    write(error_unit,*) "Problem setting option for non blocking collectives for cols in elpa1_tridiag. Aborting..."
    success = .false.
    call obj%timer%stop("tridiag_&
    &real&
    &" // &
    "_double" // &
    gpuString )
    return
  endif

  if (non_blocking_collectives_rows .eq. 1) then
    useNonBlockingCollectivesRows = .true.
  else
    useNonBlockingCollectivesRows = .false.
  endif

  if (non_blocking_collectives_cols .eq. 1) then
    useNonBlockingCollectivesCols = .true.
  else
    useNonBlockingCollectivesCols = .false.
  endif


  !mpi_comm_all    = obj%mpi_setup%mpi_comm_parent
  !mpi_comm_cols   = obj%mpi_setup%mpi_comm_cols
  !mpi_comm_rows   = obj%mpi_setup%mpi_comm_rows

  my_prow = obj%mpi_setup%myRank_comm_rows
  my_pcol = obj%mpi_setup%myRank_comm_cols

  np_rows = obj%mpi_setup%nRanks_comm_rows
  np_cols = obj%mpi_setup%nRanks_comm_cols


  !if (wantDebug) call obj%timer%start("mpi_communication")
  !call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND), my_prowMPI, mpierr)
  !call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND), np_rowsMPI, mpierr)
  !call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND), my_pcolMPI, mpierr)
  !call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND), np_colsMPI, mpierr)


  !my_prow = int(my_prowMPI, kind=c_int)
  !np_rows = int(np_rowsMPI, kind=c_int)
  !my_pcol = int(my_pcolMPI, kind=c_int)
  !np_cols = int(np_colsMPI, kind=c_int)
  !if (wantDebug) call obj%timer%stop("mpi_communication")

  ! Matrix is split into tiles; work is done only for tiles on the diagonal or above
  ! seems that tile is a square submatrix, consisting by several blocks
  ! it is a smallest possible square submatrix, where blocks being distributed among
  ! processors are "aligned" in both rows and columns
  !  -----------------
  ! | 1 4 | 1 4 | 1 4 | ...
  ! | 2 5 | 2 5 | 2 5 | ...
  ! | 3 6 | 3 6 | 3 6 | ...
  !  ----------------- ...
  ! | 1 4 | 1 4 | 1 4 | ...
  ! | 2 5 | 2 5 | 2 5 | ...
  ! | 3 6 | 3 6 | 3 6 | ...
  !  ----------------- .
  !   : :   : :   : :    .
  !   : :   : :   : :      .
  !
  ! this is a tile, where each number represents block, assigned to a processor with the shown number
  ! size of this small block is nblk
  ! Image is for situation with 6 processors, 3 processor rows and 2 columns
  ! tile_size is thus nblk * 6
  !
  tile_size = nblk*least_common_multiple(np_rows,np_cols) ! minimum global tile size

  ! make tile_size a smallest possible multiple of previously defined tile size, such that it is
  ! larger or equal to min_tile_size
  ! min_tile_size has been originally hardcoded as 128 * max(np_rows, np_cols), so it is now the implicit value
  ! it can, however, be set by the user
  call obj%get("min_tile_size", min_tile_size ,error)
  if (error .ne. ELPA_OK) then
    write(error_unit,*) "Problem setting option for min_tile_size. Aborting..."
    success = .false.
    call obj%timer%stop("tridiag_&
    &real&
    &" // &
    "_double" // &
    gpuString )
    return
  endif
  if(min_tile_size == 0) then
    ! not set by the user, use the default value
    min_tile_size = 128*max(np_rows, np_cols)
  endif
  tile_size = ((min_tile_size-1)/tile_size+1)*tile_size

  nblockEnd = 3

  l_rows_per_tile = tile_size/np_rows ! local rows of a tile
  l_cols_per_tile = tile_size/np_cols ! local cols of a tile

  totalblocks = (na-1)/nblk + 1
  max_loc_block_rows = (totalblocks-1)/np_rows + 1
  max_loc_block_cols = (totalblocks-1)/np_cols + 1

  ! localy owned submatrix has size at most max_local_rows x max_local_cols at each processor
  max_local_rows = max_loc_block_rows*nblk
  max_local_cols = max_loc_block_cols*nblk

  ! allocate memmory for vectors
  ! todo: It is little bit confusing, I think, that variables _row actually store columns and vice versa
  ! todo: if something has length max_local_rows, it is actually a column, no?
  ! todo: probably one should read it as v_row = Vector v distributed among rows

  allocate(uv_stored_cols(max_local_cols,2*max_stored_uv), stat=istat, errmsg=errorMessage)
  call check_alloc("tridiag_&
       &real ", "uv_stored_cols", istat, errorMessage)

  allocate(vu_stored_rows(max_local_rows,2*max_stored_uv), stat=istat, errmsg=errorMessage)
  call check_alloc("tridiag_&
       &real ", "vu_stored_rows", istat, errorMessage)

  if (useGPU) then

    ! allocate v_row 1 element longer to allow store and broadcast tau together with it
    if (gpu_vendor() /= OPENMP_OFFLOAD_GPU .and. gpu_vendor() /= SYCL_GPU) then
      num = (max_local_rows+1) * size_of_datatype
      successGPU = gpu_malloc_host(v_row_host, num)
      call check_host_alloc_GPU_f("tridiag: v_row_host", 444,  successGPU)
      call c_f_pointer(v_row_host,v_row,(/(max_local_rows+1)/))
    else
      allocate(v_row(max_local_rows+1))
    endif

    if (gpu_vendor() /= OPENMP_OFFLOAD_GPU .and. gpu_vendor() /= SYCL_GPU) then
      num = (max_local_cols) * size_of_datatype
      successGPU = gpu_malloc_host(v_col_host,num)
      call check_host_alloc_GPU_f("tridiag: v_col_host", 453,  successGPU)
      call c_f_pointer(v_col_host,v_col,(/(max_local_cols)/))
    else
      allocate(v_col(max_local_cols))
    endif

    if (gpu_vendor() /= OPENMP_OFFLOAD_GPU .and. gpu_vendor() /= SYCL_GPU) then
      num = (max_local_cols) * size_of_datatype
      successGPU = gpu_malloc_host(u_col_host,num)
      call check_host_alloc_GPU_f("tridiag: u_col_host", 462,  successGPU)
      call c_f_pointer(u_col_host,u_col,(/(max_local_cols)/))
    else
      allocate(u_col(max_local_cols))
    endif

    if (gpu_vendor() /= OPENMP_OFFLOAD_GPU .and. gpu_vendor() /= SYCL_GPU) then
      num = (max_local_rows) * size_of_datatype
      successGPU = gpu_malloc_host(u_row_host,num)
      call check_host_alloc_GPU_f("tridiag: u_row_host", 471,  successGPU)
      call c_f_pointer(u_row_host,u_row,(/(max_local_rows)/))
    else
      allocate(u_row(max_local_rows))
    endif

    
    if (gpu_vendor() /= OPENMP_OFFLOAD_GPU .and. gpu_vendor() /= SYCL_GPU) then
      num = (max_local_rows * 2*max_stored_uv) * size_of_datatype
      successGPU = gpu_host_register(int(loc(vu_stored_rows),kind=c_intptr_t), num, gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: vu_stored_rows", 481,  successGPU)

      num = (max_local_cols * 2*max_stored_uv) * size_of_datatype
      successGPU = gpu_host_register(int(loc(uv_stored_cols),kind=c_intptr_t), num, gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: uv_stored_cols", 485,  successGPU)

      num = (1 * 2*max_stored_uv) * size_of_datatype
      successGPU = gpu_host_register(int(loc(aux),kind=c_intptr_t), num, gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: aux", 489,  successGPU)

      num = na * size_of_datatype_real
      successGPU = gpu_host_register(int(loc(d_vec),kind=c_intptr_t), num, gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: d_vec", 493,  successGPU)

      num = na * size_of_datatype_real
      successGPU = gpu_host_register(int(loc(e_vec),kind=c_intptr_t), num, gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: e_vec", 497,  successGPU)

      num = na * size_of_datatype
      successGPU = gpu_host_register(int(loc(tau),kind=c_intptr_t), num, gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: tau", 501,  successGPU)

      num = 2 * size_of_datatype
      successGPU = gpu_host_register(int(loc(aux1),kind=c_intptr_t), num, gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: aux1", 505,  successGPU)

      num = 1 * size_of_datatype
      successGPU = gpu_host_register(int(loc(vav),kind=c_intptr_t), num, gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: vav", 509,  successGPU)

      num = 1 * size_of_datatype
      successGPU = gpu_host_register(int(loc(dot_prod),kind=c_intptr_t), num, gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: dot_prod", 513,  successGPU)

      num = 1 * size_of_datatype
      successGPU = gpu_host_register(int(loc(xf),kind=c_intptr_t), num, gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: xf", 517,  successGPU)
    endif ! gpu_vendor() /= OPENMP_OFFLOAD_GPU .and. gpu_vendor() /= SYCL_GPU
  else ! useGPU

    allocate(v_row(max_local_rows+1), stat=istat, errmsg=errorMessage)
    call check_alloc("tridiag_&
    &real ", "v_row", istat, errorMessage)

    allocate(v_col(max_local_cols), stat=istat, errmsg=errorMessage)
    call check_alloc("tridiag_&
     &real ", "v_col", istat, errorMessage)

    allocate(u_col(max_local_cols), stat=istat, errmsg=errorMessage)
    call check_alloc("tridiag_&
    &real ", "u_col", istat, errorMessage)

    allocate(u_row(max_local_rows), stat=istat, errmsg=errorMessage)
    call check_alloc("tridiag_&
    &real ", "u_row", istat, errorMessage)
      
  endif ! useGPU


  v_row = 0
  u_row = 0
  v_col = 0
  u_col = 0

  if (useGPU) then
    successGPU = gpu_malloc(v_row_dev, (max_local_rows+1) * size_of_datatype)
    call check_alloc_GPU_f("tridiag: v_row_dev", 556,  successGPU)

    successGPU = gpu_malloc(u_row_dev, max_local_rows * size_of_datatype)

    call check_alloc_GPU_f("tridiag: u_row_dev", 560,  successGPU)

    successGPU = gpu_malloc(v_col_dev, max_local_cols * size_of_datatype)
    call check_alloc_GPU_f("tridiag: v_col_dev", 563,  successGPU)

    successGPU = gpu_malloc(u_col_dev, max_local_cols * size_of_datatype)
    call check_alloc_GPU_f("tridiag: u_col_dev", 566,  successGPU)

    successGPU = gpu_malloc(vu_stored_rows_dev, max_local_rows * 2 * max_stored_uv * size_of_datatype)
    call check_alloc_GPU_f("tridiag: vu_stored_rows_dev", 569,  successGPU)

    successGPU = gpu_malloc(uv_stored_cols_dev, max_local_cols * 2 * max_stored_uv * size_of_datatype)
    call check_alloc_GPU_f("tridiag: uv_stored_cols_dev", 572,  successGPU)

    successGPU = gpu_malloc(d_vec_dev, na * size_of_datatype_real)
    call check_alloc_GPU_f("tridiag: d_vec_dev", 575,  successGPU)

    successGPU = gpu_malloc(e_vec_dev, na * size_of_datatype_real)
    call check_alloc_GPU_f("tridiag: e_vec_dev", 578,  successGPU)

    successGPU = gpu_malloc(tau_dev, na * size_of_datatype)
    call check_alloc_GPU_f("tridiag: tau_dev", 581,  successGPU)

    successGPU = gpu_malloc(aux_dev, 2*max_stored_uv * size_of_datatype)
    call check_alloc_GPU_f("tridiag: aux_dev", 584,  successGPU)

    successGPU = gpu_malloc(aux1_dev, 2 * size_of_datatype)
    call check_alloc_GPU_f("tridiag: aux1_dev", 587,  successGPU)

    successGPU = gpu_malloc(aux_complex_dev, 2 *max_stored_uv* size_of_datatype)
    call check_alloc_GPU_f("tridiag: aux_complex_dev", 590,  successGPU)

    successGPU = gpu_malloc(vav_dev, 1 * size_of_datatype)
    call check_alloc_GPU_f("tridiag: vav_dev", 593,  successGPU)

    successGPU = gpu_malloc(dot_prod_dev, 1 * size_of_datatype)
    call check_alloc_GPU_f("tridiag: dot_prod_dev", 596,  successGPU)

    successGPU = gpu_malloc(xf_dev, 1 * size_of_datatype)
    call check_alloc_GPU_f("tridiag: xf_dev", 599,  successGPU)

  endif !useGPU

  d_vec(:) = 0
  e_vec(:) = 0
  tau(:) = 0

  if (useGPU) then
    successGPU = gpu_memset(d_vec_dev, 0, na * size_of_datatype_real)
    call check_memcpy_GPU_f("tridiag: d_vec_dev", 630,  successGPU)

    successGPU = gpu_memset(e_vec_dev, 0, na * size_of_datatype_real)
    call check_memcpy_GPU_f("tridiag: e_vec_dev", 633,  successGPU)

    successGPU = gpu_memset(tau_dev, 0, na * size_of_datatype)
    call check_memcpy_GPU_f("tridiag: tau_dev", 636,  successGPU)
  endif

  n_stored_vecs = 0

  l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a_mat
  l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a_mat

  if (my_prow == prow(na, nblk, np_rows) .and. my_pcol == pcol(na, nblk, np_cols)) then
      d_vec(na) = a_mat(l_rows,l_cols)
  endif

  if (useGPU) then
    ! allocate memory for matrix A on the device and than copy the matrix

    num = matrixRows * matrixCols * size_of_datatype

    successGPU = gpu_malloc(a_dev, num)
    call check_alloc_GPU_f("tridiag: a_dev", 659,  successGPU)


    successGPU = gpu_memcpy(a_dev, int(loc(a_mat(1,1)),kind=c_intptr_t), &
                              num, gpuMemcpyHostToDevice)
    call check_memcpy_GPU_f("tridiag: a_dev", 678,  successGPU)

  endif ! useGPU

  ! main cycle of tridiagonalization
  ! in each step, 1 Householder Vector is calculated
  do istep = na, nblockEnd ,-1


    ! Calculate number of local rows and columns of the still remaining matrix
    ! on the local processor
    l_rows = local_index(istep-1, my_prow, np_rows, nblk, -1)
    l_cols = local_index(istep-1, my_pcol, np_cols, nblk, -1)

    ! Calculate Vector for Householder transformation on all procs
    ! owning column istep

    if (useGPU) then 
      ! copy l_cols + 1 column of a_dev to v_row_dev
      isOurProcessRow      = (my_prow == prow(istep-1, nblk, np_rows))
      isOurProcessCol      = (my_pcol == pcol(istep-1, nblk, np_cols))
      isOurProcessCol_prev = (my_pcol == pcol(istep  , nblk, np_cols)) ! isOurProcessCol from the previous step
      call gpu_copy_and_set_zeros_double(v_row_dev, a_dev, l_rows, l_cols, matrixRows, istep, &
                                            aux1_dev, vav_dev, d_vec_dev, &
                                            isOurProcessRow, isOurProcessCol, isOurProcessCol_prev, &
                                            isSkewsymmetric, useCCL, wantDebug, my_stream)
      if (gpu_vendor() /= SYCL_GPU) successGPU = gpu_DeviceSynchronize()
    endif ! useGPU

    if (my_pcol == pcol(istep, nblk, np_cols)) then

      ! Get Vector to be transformed; distribute last element and norm of
      ! remaining elements to all procs in current column

!             ! copy l_cols + 1 column of A to v_row
      ! if (useGPU) then

! #ifdef WITH_NVTX
!         call nvtxRangePush("memcpy new D-D a_dev(:,l_cols+1)->v_row_dev")
! #endif
!         ! TODO_23_11:  create a dev-dev copy kernel or merge it to another kernel
!         offset_dev = l_cols * matrixRows * size_of_datatype
!         successGPU = gpu_memcpy(v_row_dev, a_dev + offset_dev, (l_rows)* size_of_datatype, gpuMemcpyDeviceToDevice)
!         call check_memcpy_GPU_f("tridiag a_dev 1", 731,  successGPU)

! #ifdef WITH_NVTX
!         call nvtxRangePop()
! #endif

!       else ! useGPU
!         v_row(1:l_rows) = a_mat(1:l_rows,l_cols+1)
!       endif ! useGPU

      ! copy l_cols + 1 column of A to v_row
      if (.not. useGPU) then
        v_row(1:l_rows) = a_mat(1:l_rows,l_cols+1)
      endif ! useGPU

      if (n_stored_vecs > 0 .and. l_rows > 0) then
        if (useGPU) then
          if (wantDebug) call obj%timer%start("gpublas gemv skinny with copying")
          ! v_row_dev = vu_stored_rows_dev * uv_stored_cols_dev(l_cols+1,1:2*n_stored_vecs) + v_row_dev
          gpuHandle = obj%gpu_setup%gpublasHandleArray(0)
          call gpublas_DGEMV('N', l_rows, 2*n_stored_vecs,  &
                                    ONE, vu_stored_rows_dev, max_local_rows,  &
                                    uv_stored_cols_dev+(l_cols+1-1 +max_local_cols*(1-1))*size_of_datatype , max_local_cols, &
                                    ONE, v_row_dev, 1, gpuHandle)

          if (wantDebug .and. gpu_vendor() /= SYCL_GPU) successGPU = gpu_DeviceSynchronize()
          if (wantDebug) call obj%timer%stop("gpublas gemv skinny with copying")

        else ! useGPU

          if (wantDebug) call obj%timer%start("blas")
          ! v_row = vu_stored_rows * uv_stored_cols(l_cols+1,1:2*n_stored_vecs) + v_row
          call DGEMV('N',   &
                            int(l_rows,kind=BLAS_KIND), int(2*n_stored_vecs,kind=BLAS_KIND), &
                            ONE, vu_stored_rows, int(ubound(vu_stored_rows,dim=1),kind=BLAS_KIND), &
                            uv_stored_cols(l_cols+1,1), &
                            int(ubound(uv_stored_cols,dim=1),kind=BLAS_KIND), &
                            ONE, v_row, 1_BLAS_KIND)
          if (wantDebug) call obj%timer%stop("blas")
        endif ! useGPU
      endif ! (n_stored_vecs > 0 .and. l_rows > 0)

      if (useGPU) then
        if (my_prow == prow(istep-1, nblk, np_rows)) then
          isOurProcessRowInt = 1
        else
          isOurProcessRowInt = 0
        end if

        my_stream = obj%gpu_setup%my_stream
        call gpu_dot_product_and_assign_double(v_row_dev, l_rows, isOurProcessRowInt, aux1_dev, wantDebug, my_stream)
        if (.not. useCCL) then


          successGPU = gpu_memcpy(int(loc(aux1),kind=c_intptr_t), aux1_dev, 2*size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("tridiag: aux1_dev -> aux1", 819,  successGPU)


        endif ! .not. useCCL 
      endif ! useGPU  


      if (.not. useGPU) then
        if (my_prow == prow(istep-1, nblk, np_rows)) then
          aux1(1) = dot_product(v_row(1:l_rows-1),v_row(1:l_rows-1)) ! = "q"
          aux1(2) = v_row(l_rows) ! = "a_11" (or rather a_nn)
        else
          aux1(1) = dot_product(v_row(1:l_rows),v_row(1:l_rows))
          aux1(2) = 0.
        endif
      endif ! .not. useGPU 


      if (useCCL) then
      else ! useCCL
        vnorm2 = aux1(1)
        vrl    = aux1(2)
      endif ! useCCL


      ! Householder transformation
      if (useCCL) then
      else ! useCCL
        call hh_transform_real_&
                &double &
                (obj, vrl, vnorm2, xf, tau(istep), wantDebug)
      endif ! useCCL
      
      ! vrl is newly computed off-diagonal element of the final tridiagonal matrix
      if (my_prow == prow(istep-1, nblk, np_rows)) then
        if (.not. useCCL) then
          e_vec(istep-1) = vrl
        endif ! .not. useCCL
      endif

      if (.not. useGPU) then
        ! Scale v_row and store Householder Vector for back transformation
        v_row(1:l_rows) = v_row(1:l_rows) * xf

        if (my_prow == prow(istep-1, nblk, np_rows)) then
          v_row(l_rows) = 1.
        endif

        ! store Householder Vector for back transformation
        ! update a_mat
        a_mat(1:l_rows,l_cols+1) = v_row(1:l_rows)
      endif ! .not. useGPU 

      if (.not. useCCL) then
        ! add tau after the end of actuall v_row, to be broadcasted with it
        v_row(l_rows+1) = tau(istep)
      endif

      
      if (useGPU) then
        if (useCCL) then
          xf_host_or_dev = xf_dev
        else 
          xf_host_or_dev = int(loc(xf),kind=c_intptr_t)
        endif       
        
        isOurProcessRow = (my_prow == prow(istep-1, nblk, np_rows))
        call gpu_set_e_vec_scale_set_one_store_v_row_double(e_vec_dev, vrl_dev, a_dev, v_row_dev, tau_dev, xf_host_or_dev, & 
                                                  l_rows, l_cols, matrixRows, istep, isOurProcessRow, useCCL, wantDebug, my_stream)
      endif ! useGPU  


      if (useGPU .and. .not. useCCL) then      
        !v_row_dev -> v_row

        successGPU = gpu_memcpy(int(loc(v_row),kind=c_intptr_t), v_row_dev, l_rows*size_of_datatype, gpuMemcpyDeviceToHost)
        call check_memcpy_GPU_f("tridiag: v_row_dev -> v_row", 995,  successGPU)

      endif ! useGPU .and. .not. useCCL

    endif !(my_pcol == pcol(istep, nblk, np_cols))



    !recover tau, which has been broadcasted together with v_row
    if (.not. useCCL) then
      tau(istep) =  v_row(l_rows+1)
    endif

    ! Transpose Householder Vector v_row -> v_col
    if (useCCL) then
    else ! useCCL
      call elpa_transpose_vectors_&
          &real&
          &_&
          &double &
                (obj, v_row, ubound(v_row,dim=1), mpi_comm_rows, v_col, ubound(v_col,dim=1), mpi_comm_cols, &
                1, istep-1, 1, nblk, max_threads, .true., success)
      if (wantDebug .and. gpu_vendor() /= SYCL_GPU) successGPU = gpu_DeviceSynchronize()
      if (.not.(success)) then
        write(error_unit,*) "Error in elpa_transpose_vectors. Aborting!"
        return
      endif
    endif ! useCCL

    ! Calculate u = (A + VU**T + UV**T)*v // Dongarra 1987: "y = (A - UV**T - VU**T)*u"

    ! For cache efficiency, we use only the upper half of the matrix tiles for this,
    ! thus the result is partly in u_col(:) and partly in u_row(:)

    if (.not. useGPU) then
      u_col(1:l_cols) = 0
      u_row(1:l_rows) = 0
    endif

    if (useGPU .and. useCCL) then
      successGPU = gpu_memset(u_col_dev, 0, l_cols * size_of_datatype) ! TODO_23_11: omit this, but change gpublas_gemm to u_col_dev=a_dev^T*v_row_dev+0*u_col_dev?
      call check_memcpy_GPU_f("tridiag: u_col_dev", 1101,  successGPU)
    endif

    if (l_rows>0 .and. l_cols>0) then
      if (useGPU) then

        if (.not. mat_vec_as_one_block) then

          successGPU = gpu_memcpy(v_col_dev, int(loc(v_col(1)),kind=c_intptr_t), &
                        l_cols * size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("tridiag: v_col_dev", 1146,  successGPU)


        endif ! .not. mat_vec_as_one_block

        if (.not. useCCL) then

          successGPU = gpu_memcpy(v_row_dev, int(loc(v_row(1)),kind=c_intptr_t), &
                                    l_rows * size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("tridiag: v_row_dev", 1170,  successGPU)


        endif ! .not. useCCL
      endif ! useGPU


      if (.not. useGPU) then
        do i=0, (istep-2)/tile_size ! iteration over tiles
          l_col_beg = i*l_cols_per_tile+1
          l_col_end = min(l_cols,(i+1)*l_cols_per_tile)
          if (l_col_end < l_col_beg) cycle
          do j = 0, i
            l_row_beg = j*l_rows_per_tile+1
            l_row_end = min(l_rows,(j+1)*l_rows_per_tile)
            if (l_row_end < l_row_beg) cycle

            ! multiplication by blocks is efficient only for CPU
            ! for GPU we introduced 2 other ways, either by stripes (more simmilar to the original
            ! CPU implementation) or by one large matrix Vector multiply
            if (wantDebug) call obj%timer%start("blas")
            ! u_col = a_mat*v_row + u_col(=0)
            call DGEMV('T',  &
                        int(l_row_end-l_row_beg+1,kind=BLAS_KIND), int(l_col_end-l_col_beg+1,kind=BLAS_KIND), &
                        ONE, a_mat(l_row_beg, l_col_beg), int(matrixRows,kind=BLAS_KIND),         &
                        v_row(l_row_beg:max_local_rows+1), 1_BLAS_KIND,                           &
                        ONE, u_col(l_col_beg:max_local_cols), 1_BLAS_KIND)

            if (i/=j) then
              if (isSkewsymmetric) then
                call DGEMV('N',int(l_row_end-l_row_beg+1,kind=BLAS_KIND), int(l_col_end-l_col_beg+1,kind=BLAS_KIND), &
                                    -ONE, a_mat(l_row_beg,l_col_beg), int(matrixRows,kind=BLAS_KIND),               &
                                    v_col(l_col_beg:max_local_cols), 1_BLAS_KIND, ONE, u_row(l_row_beg:max_local_rows), &
                                    1_BLAS_KIND)

              else
                ! u_row = a_mat*v_col + u_row(=0)
                call DGEMV('N',int(l_row_end-l_row_beg+1,kind=BLAS_KIND), int(l_col_end-l_col_beg+1,kind=BLAS_KIND),  &
                                    ONE, a_mat(l_row_beg,l_col_beg), int(matrixRows,kind=BLAS_KIND),               &
                                    v_col(l_col_beg:max_local_cols), 1_BLAS_KIND, ONE, u_row(l_row_beg:max_local_rows), &
                                    1_BLAS_KIND)
              endif
            endif
            if (wantDebug) call obj%timer%stop("blas")

          enddo  ! j=0,i
        enddo  ! i=0,(istep-2)/tile_size
      endif ! .not. useGPU

      if (useGPU) then
        if (mat_vec_as_one_block) then
          ! Unlike for CPU, we (for each MPI thread) do just one large mat-vec multiplication
          ! this requires altering of the algorithm when later explicitly updating the matrix
          ! after max_stored_uv is reached : we need to update all tiles, not only those above diagonal
          if (wantDebug) call obj%timer%start("gpublas")
          
          gpuHandle = obj%gpu_setup%gpublasHandleArray(0) ! 0.3-0.5 us

              
          ! u_col_dev = a_dev^T*v_row_dev
          call gpublas_DGEMV('T', l_rows,l_cols,  &
                                    ONE, a_dev, matrixRows,                   &
                                    v_row_dev , 1,                          &
                                    ZERO, u_col_dev, 1, gpuHandle)
              
          if (wantDebug .and. gpu_vendor() /= SYCL_GPU) successGPU = gpu_DeviceSynchronize()
              
       ! todo: try with non transposed!!!
!                 if(i/=j) then
!                   call gpublas_DGEMV('N', l_row_end-l_row_beg+1,l_col_end-l_col_beg+1,  &
!                                             ONE, a_dev + offset_dev, matrixRows,                        &
!                                             v_col_dev + (l_col_beg - 1) *                      &
!                                             size_of_datatype, 1,                          &
!                                             ONE, u_row_dev + (l_row_beg - 1) *                 &
!                                             size_of_datatype, 1)
!                 endif
          if (wantDebug) call obj%timer%stop("gpublas")
        else  ! mat_vec_as_one_block
          !perform multiplication by stripes - it is faster than by blocks, since we call cublas with
          !larger matrices. In general, however, this algorithm is very simmilar to the one with CPU
          do i=0,(istep-2)/tile_size
            l_col_beg = i*l_cols_per_tile+1
            l_col_end = min(l_cols,(i+1)*l_cols_per_tile)
            if (l_col_end<l_col_beg) cycle

            l_row_beg = 1
            l_row_end = min(l_rows,(i+1)*l_rows_per_tile)
                  
            offset_dev = ((l_row_beg-1) + (l_col_beg - 1) * matrixRows) * size_of_datatype
  
            gpuHandle = obj%gpu_setup%gpublasHandleArray(0)
            call gpublas_DGEMV('T', &
                          l_row_end-l_row_beg+1, l_col_end-l_col_beg+1, &
                          ONE, a_dev + offset_dev, matrixRows,  &
                          v_row_dev + (l_row_beg - 1) * size_of_datatype, 1,  &
                          ONE, u_col_dev + (l_col_beg - 1) * size_of_datatype, 1, gpuHandle)
          enddo !i=0,(istep-2)/tile_size

          do i=0,(istep-2)/tile_size
            l_col_beg = i*l_cols_per_tile+1
            l_col_end = min(l_cols,(i+1)*l_cols_per_tile)
            if (l_col_end<l_col_beg) cycle

            l_row_beg = 1
            l_row_end = min(l_rows,i*l_rows_per_tile)
              
            offset_dev = ((l_row_beg-1) + (l_col_beg - 1) * matrixRows) * size_of_datatype

            if (isSkewsymmetric) then
                gpuHandle = obj%gpu_setup%gpublasHandleArray(0)
                call gpublas_DGEMV('N', l_row_end-l_row_beg+1, l_col_end-l_col_beg+1, &
                            -ONE, a_dev + offset_dev, matrixRows, &
                            v_col_dev + (l_col_beg - 1) * size_of_datatype,1, &
                            ONE, u_row_dev + (l_row_beg - 1) * size_of_datatype, 1, gpuHandle)
            else
                gpuHandle = obj%gpu_setup%gpublasHandleArray(0)
                call gpublas_DGEMV('N', l_row_end-l_row_beg+1, l_col_end-l_col_beg+1, &
                            ONE, a_dev + offset_dev, matrixRows, &
                            v_col_dev + (l_col_beg - 1) * size_of_datatype,1, &
                            ONE, u_row_dev + (l_row_beg - 1) * size_of_datatype, 1, gpuHandle)
            endif
          enddo ! i=0,(istep-2)/tile_size
        end if ! mat_vec_as_one_block / per stripes


        if (.not. mat_vec_as_one_block) then

          successGPU = gpu_memcpy(int(loc(u_row(1)),kind=c_intptr_t), u_row_dev, l_rows*size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("tridiag: u_row_dev 1", 1378,  successGPU)


        endif ! .not. mat_vec_as_one_block
      endif ! useGPU


      ! second calculate (VU**T + UV**T)*v part of (A + VU**T + UV**T)*v
      if (n_stored_vecs > 0) then
        if (.not. useGPU) then
          if (wantDebug) call obj%timer%start("blas")

      
          ! aux = vu_stored_rows^T*v_row
          call DGEMV('T',     &
                              int(l_rows,kind=BLAS_KIND), int(2*n_stored_vecs,kind=BLAS_KIND),   &
                              ONE, vu_stored_rows, int(ubound(vu_stored_rows,dim=1),kind=BLAS_KIND),   &
                              v_row,  1_BLAS_KIND, ZERO, aux, 1_BLAS_KIND)
              
          ! u_col = uv_stored_cols*aux + u_col
          call DGEMV('N', int(l_cols,kind=BLAS_KIND), int(2*n_stored_vecs,kind=BLAS_KIND),   &
                              ONE, uv_stored_cols, int(ubound(uv_stored_cols,dim=1),kind=BLAS_KIND),   &
                              aux, 1_BLAS_KIND, ONE, u_col,  1_BLAS_KIND)
          if (wantDebug) call obj%timer%stop("blas")
        else ! .not. useGPU
          if (wantDebug) call obj%timer%start("gpublas gemv x2 skinny")


          ! aux_dev = vu_stored_rows_dev^T*v_row_dev
          gpuHandle = obj%gpu_setup%gpublasHandleArray(0)
          call gpublas_DGEMV('T', l_rows, 2*n_stored_vecs, &
                                      ONE, vu_stored_rows_dev, max_local_rows,   &
                                      v_row_dev,  1, ZERO, aux_dev, 1, gpuHandle)
              
          ! u_col_dev = uv_stored_cols_dev*aux_dev + u_col_dev
          call gpublas_DGEMV('N', l_cols, 2*n_stored_vecs, &
                                      ONE, uv_stored_cols_dev, max_local_cols,   &
                                      aux_dev, 1, ONE, u_col_dev, 1, gpuHandle)

          if (wantDebug .and. gpu_vendor() /= SYCL_GPU) successGPU = gpu_DeviceSynchronize()

          if (wantDebug) call obj%timer%stop("gpublas gemv x2 skinny")  
        endif ! .not. useGPU
      endif ! n_stored_vecs > 0

    endif  ! (l_rows>0 .and. l_cols>0)

    if (useGPU .and. l_cols>0 .and. (.not. useCCL)) then

      successGPU = gpu_memcpy(int(loc(u_col(1)),kind=c_intptr_t), u_col_dev, l_cols*size_of_datatype, gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("tridiag: u_col_dev 1", 1462,  successGPU)

    endif ! useGPU


    if (useGPU .and. .not. useCCL) then 
        if (l_rows==0) then
        u_col(1:l_cols) = 0
      endif
    endif ! useGPU

    ! Sum up all u_row(:) parts along rows and add them to the u_col(:) parts
    ! on the processors containing the diagonal
    ! This is only necessary if u_row has been calculated, i.e. if the
    ! global tile size is smaller than the global remaining matrix
    
    ! in GPU case, u_row_dev=0 up to now, so elpa_reduce_add_vectors is skipped

    if (tile_size < istep-1 .and. .not. useGPU) then
      call elpa_reduce_add_vectors_&
            &real&
            &_&
            &double &
            (obj, u_row, ubound(u_row,dim=1), mpi_comm_rows, u_col, ubound(u_col,dim=1), &
            mpi_comm_cols, istep-1, 1, nblk, max_threads)
    endif ! (tile_size < istep-1 .and. .not. useGPU)

    ! Sum up all the u_col(:) parts, transpose u_col -> u_row

    if (l_cols > 0) then
    endif ! (l_cols > 0)
    
    ! Transpose Householder Vector u_col -> u_row
    if (useCCL) then
    else ! useCCL
      if (isSkewsymmetric) then
        call elpa_transpose_vectors_ss_&
        &real&
        &_&
        &double &
        (obj, u_col, ubound(u_col,dim=1), mpi_comm_cols, u_row, ubound(u_row,dim=1), &
          mpi_comm_rows, 1, istep-1, 1, nblk, max_threads, .false., success)
        if (.not.(success)) then
          write(error_unit,*) "Error in elpa_transpose_vectors_ss. Aborting!"
          return
        endif
      else ! isSkewsymmetric
        call elpa_transpose_vectors_&
        &real&
        &_&
        &double &
        (obj, u_col, ubound(u_col,dim=1), mpi_comm_cols, u_row, ubound(u_row,dim=1), &
        mpi_comm_rows, 1, istep-1, 1, nblk, max_threads, .false., success)
        if (.not.(success)) then
          write(error_unit,*) "Error in elpa_transpose_vectors. Aborting!"
          return
        endif
      endif ! isSkewsymmetric
    endif ! useCCL


    if (useGPU .and. .not. useCCL) then

      successGPU = gpu_memcpy(v_col_dev, int(loc(v_col(1)),kind=c_intptr_t), l_cols*size_of_datatype, gpuMemcpyHostToDevice)
      call check_memcpy_GPU_f("tridiag: v_col_dev", 1607,  successGPU)


      successGPU = gpu_memcpy(u_col_dev, int(loc(u_col(1)),kind=c_intptr_t), l_cols*size_of_datatype, gpuMemcpyHostToDevice)
      call check_memcpy_GPU_f("tridiag: u_col_dev", 1622,  successGPU)

    endif ! (useGPU .and. .not. useCCL)



    if (useGPU .and. .not. useCCL) then

      successGPU = gpu_memcpy(u_row_dev, int(loc(u_row(1)),kind=c_intptr_t), l_rows*size_of_datatype, gpuMemcpyHostToDevice)
      call check_memcpy_GPU_f("tridiag: u_row_dev", 1657,  successGPU)


      successGPU = gpu_memcpy(v_row_dev, int(loc(v_row(1)),kind=c_intptr_t), l_rows*size_of_datatype, gpuMemcpyHostToDevice)
      call check_memcpy_GPU_f("tridiag: v_row_dev", 1671,  successGPU)

    endif ! (useGPU .and. .not. useCCL)


    ! calculate u**T * v (same as v**T * (A + VU**T + UV**T) * v )
    if (.not. useGPU .or. (useGPU .and. .not. useCCL)) then
      vav = 0 ! x=0
      if (l_cols>0) vav = dot_product(v_col(1:l_cols), u_col(1:l_cols))
    endif



    ! store u and v in the matrices U and V
    ! these matrices are stored combined in one here

   if (.not. useCCL) then
      conjg_tau = tau(istep)
    endif ! .not. useCCL
    
    if (.not. useGPU) then
      if (l_rows > 0) then
        ! update vu_stored_rows
        vu_stored_rows(1:l_rows,2*n_stored_vecs+1) = conjg_tau*v_row(1:l_rows)
        vu_stored_rows(1:l_rows,2*n_stored_vecs+2) = 0.5*conjg_tau*vav*v_row(1:l_rows) - u_row(1:l_rows)
      endif
      if (l_cols > 0) then
        ! update uv_stored_cols
        uv_stored_cols(1:l_cols,2*n_stored_vecs+1) = 0.5*conjg_tau*vav*v_col(1:l_cols) - u_col(1:l_cols)
        uv_stored_cols(1:l_cols,2*n_stored_vecs+2) = conjg_tau*v_col(1:l_cols)
      endif
    endif ! .not. useGPU


    if (useGPU) then

      ! kernel: update vu_stored_rows_dev, uv_stored_cols
      ! then cpu's "store u,v in U,V" can be deleted. But we should take care of dot_prod below, where vu_stored_rows and uv_stored_cols are used
      if (useCCL) then
        vav_host_or_dev = vav_dev
        !tau_istep_host_or_dev = tau_dev + (istep-1)*size_of_datatype
        tau_istep_host_or_dev = v_row_dev + (l_rows+1-1)*size_of_datatype
      else ! useCCL
        vav_host_or_dev = int(loc(vav),kind=c_intptr_t)
        tau_istep_host_or_dev = int(loc(tau(istep)), kind=c_intptr_t)
      endif ! useCCL
      call gpu_store_u_v_in_uv_vu_double(vu_stored_rows_dev, uv_stored_cols_dev, v_row_dev, u_row_dev, &
                                          v_col_dev, u_col_dev, tau_dev, aux_complex_dev, &
                                          vav_host_or_dev, tau_istep_host_or_dev, &
                                          l_rows, l_cols, n_stored_vecs,  max_local_rows, max_local_cols, istep, &
                                          useCCL, wantDebug, my_stream)
    endif ! useGPU


    ! We have calculated another Householder Vector, number of implicitly stored increased
    n_stored_vecs = n_stored_vecs+1

    ! If the limit of max_stored_uv is reached, calculate A + VU**T + UV**T
    if (n_stored_vecs == max_stored_uv .or. istep == 3) then        
      
      if (.not. useGPU .OR. .not. mat_vec_as_one_block) then
        do i = 0, (istep-2)/tile_size
          ! go over tiles above (or on) the diagonal
          l_col_beg = i*l_cols_per_tile+1
          l_col_end = min(l_cols,(i+1)*l_cols_per_tile)
          l_row_beg = 1
          l_row_end = min(l_rows,(i+1)*l_rows_per_tile)
          if (l_col_end<l_col_beg .or. l_row_end<l_row_beg) then
            cycle
          endif

          if (useGPU) then
            if (.not. mat_vec_as_one_block) then
              ! if using mat-vec multiply by stripes, it is enough to update tiles above (or on) the diagonal only
              ! we than use the same calls as for CPU version
              if (wantDebug) call obj%timer%start("gpublas_gemm")
              ! update a_dev
              ! a_dev = vu_stored_rows_dev*uv_stored_cols_dev + a_dev
              gpuHandle = obj%gpu_setup%gpublasHandleArray(0)
              call gpublas_DGEMM('N', 'T',     &
                                      l_row_end-l_row_beg+1, l_col_end-l_col_beg+1, 2*n_stored_vecs,                      &
                                      ONE, vu_stored_rows_dev + (l_row_beg - 1) *                                         &
                                      size_of_datatype,  &
                                      max_local_rows, uv_stored_cols_dev + (l_col_beg - 1) *                              &
                                      size_of_datatype,  &
                                      max_local_cols, ONE, a_dev + ((l_row_beg - 1) + (l_col_beg - 1) * matrixRows) *     &
                                      size_of_datatype , matrixRows, gpuHandle)
              if (wantDebug .and. gpu_vendor() /= SYCL_GPU) successGPU = gpu_DeviceSynchronize()
              if (wantDebug) call obj%timer%stop("gpublas_gemm")
            endif ! .not. mat_vec_as_one_block
          else ! useGPU
            if (wantDebug) call obj%timer%start("blas_gemm")
            ! update a_mat
            ! a_mat = vu_stored_rows*uv_stored_cols + a_mat
            call DGEMM('N', 'T',                &
                                int(l_row_end-l_row_beg+1,kind=BLAS_KIND), int(l_col_end-l_col_beg+1,kind=BLAS_KIND), &
                                int(2*n_stored_vecs,kind=BLAS_KIND),    &
                                ONE, vu_stored_rows(l_row_beg:max_local_rows,1:2*max_stored_uv), &
                                int(ubound(vu_stored_rows,dim=1),kind=BLAS_KIND),   &
                                uv_stored_cols(l_col_beg,1), &
                                int(ubound(uv_stored_cols,dim=1),kind=BLAS_KIND),        &
                                ONE, a_mat(l_row_beg,l_col_beg), int(matrixRows,kind=BLAS_KIND))
            if (wantDebug) call obj%timer%stop("blas_gemm")
          endif !useGPU
        enddo ! i = 0, (istep-2)/tile_size

      else !.not. useGPU or .not. mat_vec_as_one_block (i.e. useGPU and mat_vec_as_one_block)

        !update whole (remaining) part of matrix, including tiles below diagonal
        !we can do that in one large cublas call
        if (wantDebug) call obj%timer%start("gpublas_gemm")


        gpuHandle = obj%gpu_setup%gpublasHandleArray(0)
        ! update a_dev
        ! a_dev = vu_stored_rows_dev*uv_stored_cols_dev + a_dev
        call gpublas_DGEMM('N', 'T', l_rows, l_cols, 2*n_stored_vecs,   &
                                  ONE, vu_stored_rows_dev, max_local_rows, &
                                  uv_stored_cols_dev, max_local_cols,  &
                                  ONE, a_dev, matrixRows, gpuHandle)

        if (wantDebug .and. gpu_vendor() /= SYCL_GPU) successGPU = gpu_DeviceSynchronize()
        if (wantDebug) call obj%timer%stop("gpublas_gemm")
        
        ! copy real(a_dev(l_rows,l_cols)) -> d_vec_dev(istep-1) for correct initial value of d_vec_dev before atomicAdd
        if ((.not. isSkewsymmetric) .and. &
            (my_prow == prow(istep-1, nblk, np_rows) .and. my_pcol == pcol(istep-1, nblk, np_cols))) then
          offset_dev = ((l_rows-1) + (l_cols-1)*matrixRows) * size_of_datatype
          successGPU = gpu_memcpy(d_vec_dev + (istep-2)*size_of_datatype_real, &
                                  a_dev + offset_dev, 1*size_of_datatype_real, gpuMemcpyDeviceToDevice)
          call check_memcpy_GPU_f("tridiag a_dev->d_vec_dev", 1863,  successGPU)
        endif

      endif !.not. useGPU or .not. mat_vec_as_one_block

      n_stored_vecs = 0
    endif ! (n_stored_vecs == max_stored_uv .or. istep == 3)

    if (my_prow == prow(istep-1, nblk, np_rows) .and. my_pcol == pcol(istep-1, nblk, np_cols)) then

      if (useGPU) then
        my_stream = obj%gpu_setup%my_stream
        call gpu_update_matrix_element_add_double(vu_stored_rows_dev, uv_stored_cols_dev, a_dev, d_vec_dev, &
                                            l_rows, l_cols, matrixRows, max_local_rows, max_local_cols, istep, n_stored_vecs, &
                                            isSkewsymmetric, wantDebug, my_stream)

      else ! useGPU
        if (n_stored_vecs > 0) then
          ! update a_mat (only one elememt!)
          dot_prod = dot_product(vu_stored_rows(l_rows,1:2*n_stored_vecs), uv_stored_cols(l_cols,1:2*n_stored_vecs))
          a_mat(l_rows,l_cols) = a_mat(l_rows,l_cols) + dot_prod
        endif
        if (isSkewsymmetric) then
          d_vec(istep-1) = 0.0_rk
        else
          d_vec(istep-1) = a_mat(l_rows,l_cols)
        endif
      endif ! useGPU

    endif ! (my_prow == prow(istep-1, nblk, np_rows) .and. my_pcol == pcol(istep-1, nblk, np_cols))


  enddo ! main cycle over istep=na,3,-1


  if (useGPU) then
    ! copy a_dev -> a_mat for backtransformation
    num = matrixRows * matrixCols * size_of_datatype
    successGPU = gpu_memcpy(int(loc(a_mat(1,1)),kind=c_intptr_t), a_dev, num, gpuMemcpyDeviceToHost)
    call check_memcpy_GPU_f("tridiag: a_dev", 1921,  successGPU)

  endif ! useGPU



  ! Store e_vec(1)

  if (my_prow==prow(1, nblk, np_rows) .and. my_pcol==pcol(2, nblk, np_cols)) then
    if (useGPU) then
      successGPU = gpu_memcpy(int(loc(e_vec(1)),kind=c_intptr_t), a_dev + (matrixRows * (l_cols - 1)) * size_of_datatype, &
                              1 * size_of_datatype, gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("tridiag: a_dev 7", 2026,  successGPU)
    else !useGPU
      e_vec(1) = a_mat(1,l_cols) ! use last l_cols value of loop above
    endif !useGPU
  endif ! if (my_prow==prow(1, nblk, np_rows) .and. my_pcol==pcol(2, nblk, np_cols))

  ! Store d_vec(1)
  if (my_prow==prow(1, nblk, np_rows) .and. my_pcol==pcol(1, nblk, np_cols)) then
    if(useGPU) then
      successGPU = gpu_memcpy(int(loc(d_vec(1)),kind=c_intptr_t), a_dev, 1 * size_of_datatype, gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("tridiag: a_dev 8", 2045,  successGPU)
    else !useGPU
      if (isSkewsymmetric) then
        d_vec(1) = 0.0_rk
      else
        d_vec(1) = a_mat(1,1)
      endif
    endif !useGPU
  endif ! (my_prow==prow(1, nblk, np_rows) .and. my_pcol==pcol(1, nblk, np_cols))

  if (useGPU) then
    offset_dev = 1 * size_of_datatype_real
    ! first and last elements of d_vec are treated separately
    successGPU = gpu_memcpy(int(loc(d_vec(2)),kind=c_intptr_t), &
                            d_vec_dev + offset_dev, (na-2) * size_of_datatype_real, gpuMemcpyDeviceToHost)
    call check_memcpy_GPU_f("tridiag: d_vec", 2062,  successGPU)

    if (useCCL) then
      ! e_vec(1) is treated separately
      offset_dev = 1 * size_of_datatype_real
      successGPU = gpu_memcpy(int(loc(e_vec(2)),kind=c_intptr_t), &
                              e_vec_dev + offset_dev, (na-1) * size_of_datatype_real, gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("tridiag: e_vec", 2069,  successGPU)

      ! tau(2) is treated separately, tau(1) is not used
      offset_dev = 2 * size_of_datatype
      successGPU = gpu_memcpy(int(loc(tau(3)),kind=c_intptr_t), &
                              tau_dev + offset_dev, (na-2) * size_of_datatype, gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("tridiag: tau", 2075,  successGPU)
    endif

    ! todo: should we leave a_mat on the device for further use?
    successGPU = gpu_free(a_dev)
    call check_dealloc_GPU_f("tridiag: a_dev 9", 2080,  successGPU)

    successGPU = gpu_free(v_row_dev)
    call check_dealloc_GPU_f("tridiag: v_row_dev", 2083,  successGPU)

    successGPU = gpu_free(u_row_dev)
    call check_dealloc_GPU_f("tridiag: (u_row_dev", 2086,  successGPU)

    successGPU = gpu_free(v_col_dev)
    call check_dealloc_GPU_f("tridiag: v_col_dev", 2089,  successGPU)

    successGPU = gpu_free(u_col_dev)
    call check_dealloc_GPU_f("tridiag: u_col_dev ", 2092,  successGPU)

    successGPU = gpu_free(vu_stored_rows_dev)
    call check_dealloc_GPU_f("tridiag: vu_stored_rows_dev ", 2095,  successGPU)

    successGPU = gpu_free(uv_stored_cols_dev)
    call check_dealloc_GPU_f("tridiag:uv_stored_cols_dev ", 2098,  successGPU)

    successGPU = gpu_free(d_vec_dev)
    call check_dealloc_GPU_f("tridiag: d_vec_dev", 2101,  successGPU)

    successGPU = gpu_free(e_vec_dev)
    call check_dealloc_GPU_f("tridiag: e_vec_dev", 2104,  successGPU)

    successGPU = gpu_free(tau_dev)
    call check_dealloc_GPU_f("tridiag: tau_dev", 2107,  successGPU)

    successGPU = gpu_free(aux_dev)
    call check_dealloc_GPU_f("tridiag: aux_dev", 2110,  successGPU)

    successGPU = gpu_free(aux1_dev)
    call check_dealloc_GPU_f("tridiag: aux1_dev", 2113,  successGPU)

    successGPU = gpu_free(aux_complex_dev)
    call check_dealloc_GPU_f("tridiag: aux_complex_dev", 2116,  successGPU)

    successGPU = gpu_free(vav_dev)
    call check_dealloc_GPU_f("tridiag: vav_dev", 2119,  successGPU)

    successGPU = gpu_free(dot_prod_dev)
    call check_dealloc_GPU_f("tridiag: dot_prod_dev", 2122,  successGPU)

    successGPU = gpu_free(xf_dev)
    call check_dealloc_GPU_f("tridiag: xf_dev", 2125,  successGPU)

  endif ! useGPU

  ! distribute the arrays d_vec and e_vec to all processors

  allocate(tmp_real(na), stat=istat, errmsg=errorMessage)
  call check_allocate_f("tridiag: tmp_real", 2138,  istat,  errorMessage)


  deallocate(tmp_real, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag: tmp_real", 2187,  istat,  errorMessage)

  if (useGPU) then

  else ! useGPU
    deallocate(v_row, v_col, u_row, u_col, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("tridiag: v_row, v_col, u_row, u_col", 2260,  istat,  errorMessage)
  endif ! useGPU

  deallocate(vu_stored_rows, uv_stored_cols, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag: vu_stored_rows, uv_stored_cols", 2264,  istat,  errorMessage)

  deallocate(aux, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag: aux", 2267,  istat,  errorMessage)



  call obj%timer%stop("tridiag_&
  &real&
  &" // &
  "_double" // &
  gpuString )

end subroutine tridiag_&
&real&
&_&
&double





!> \brief Transforms the eigenvectors of a tridiagonal matrix back
!>                     to the eigenvectors of the original matrix
!>                     (like Scalapack Routine PDORMTR)
!>
!  Parameters
!
!> \param na          Order of matrix a_mat, number of rows of matrix q_mat
!>
!> \param nqc         Number of columns of matrix q_mat
!>
!> \param a_mat(lda,matrixCols)  Matrix containing the Householder vectors (i.e. matrix a after tridiag_real)
!>                           Distribution is like in Scalapack.
!>
!> \param lda         Leading dimension of a_mat
!>
!> \param tau(na)     Factors of the Householder vectors
!>
!> \param q_mat           On input: Eigenvectors of tridiagonal matrix
!>                    On output: Transformed eigenvectors
!>                    Distribution is like in Scalapack.
!>
!> \param ldq         Leading dimension of q_mat
!>
!> \param nblk        blocksize of cyclic distribution, must be the same in both directions!
!>
!> \param matrixCols  local columns of matrix a_mat and q_mat
!>
!> \param mpi_comm_rows        MPI-Communicator for rows
!>
!> \param mpi_comm_cols        MPI-Communicator for columns
!>
!> \param useGPU      If true,  GPU version of the subroutine will be used
!>


subroutine trans_ev_&
&real&
&_&
&double &
(obj, na, nqc, a_mat, lda, tau, q_mat, ldq, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols, useGPU, success)
  use, intrinsic :: iso_c_binding
  use precision
  use elpa_abstract_impl
  use elpa_blas_interfaces
  use elpa_gpu
  use elpa_gpu_util

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

  class(elpa_abstract_impl_t), intent(inout)    :: obj
  integer(kind=ik), intent(in)                  :: na, nqc, lda, ldq, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols
  real(kind=rck), intent(in)           :: tau(na)

  real(kind=rck), intent(inout)        :: a_mat(lda,*)
  real(kind=rck), intent(inout)        :: q_mat(ldq,*)
  logical, intent(in)                           :: useGPU
  integer(kind=ik)                              :: max_stored_rows, max_stored_rows_fac

  integer(kind=ik)                              :: my_prow, my_pcol, np_rows, np_cols
  integer(kind=MPI_KIND)                        :: mpierr, my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI
  integer(kind=ik)                              :: totalblocks, max_blocks_row, max_blocks_col, max_local_rows, max_local_cols
  integer(kind=ik)                              :: l_cols, l_rows, l_colh, nstor
  integer(kind=ik)                              :: istep, n, nc, ic, ics, ice, nb, cur_pcol
  integer(kind=ik)                              :: hvn_ubnd, hvm_ubnd

  real(kind=rck), allocatable          :: hvb(:), hvm(:,:)
  real(kind=rck), pointer              :: tmp1(:), tmp2(:)
  real(kind=rck), allocatable          :: h1(:), h2(:), tmp_debug(:)
  real(kind=rck), pointer              :: tmat(:,:)
  real(kind=rck), pointer              :: hvm1(:)
  type(c_ptr)                                   :: tmp1_host, tmp2_host
  type(c_ptr)                                   :: hvm1_host, tmat_host

  integer(kind=ik)                              :: istat
  character(200)                                :: errorMessage
  character(20)                                 :: gpuString

  integer(kind=c_intptr_t)                      :: num
  integer(kind=C_intptr_T)                      :: q_dev, tmp_dev, hvm_dev, tmat_dev

  integer(kind=ik)                              :: blockStep
  logical                                       :: successGPU
  integer(kind=c_intptr_t), parameter           :: size_of_datatype = size_of_&
                                                                      &double&
                                                                      &_&
                                                                      &real
  integer(kind=ik)                              :: error
  integer(kind=MPI_KIND)                        :: bcast_request1, allreduce_request1, allreduce_request2
  logical                                       :: useNonBlockingCollectivesCols
  logical                                       :: useNonBlockingCollectivesRows
  integer(kind=c_int)                           :: non_blocking_collectives_rows, non_blocking_collectives_cols
  logical                                       :: success
  integer(kind=c_intptr_t)                      :: gpuHandle, my_stream


! NCCL requires streams
!#if defined (USE_CCL_TRANS_EV) && !defined(WITH_GPU_STREAMS)
!  successGPU = cuda_stream_create(obj%gpu_setup%my_stream)
!  if (.not.(successGPU)) then
!    print *,"Cannot create gpu stream handle"
!  endif
!  my_stream = obj%gpu_setup%my_stream
!#endif

  success = .true.

  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  call obj%timer%start("trans_ev_&
  &real&
  &" // &
  &"_double" //&
  gpuString)

  call obj%get("nbc_row_elpa1_tridi_to_full", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    write(error_unit,*) "Problem setting option for non blocking collectives for rows in elpa1_tridi_to_full. Aborting..."
    call obj%timer%stop("trans_ev_&
    &real&
    &" // &
    &"_double" //&
    gpuString)
    success = .false.
    return
  endif

  call obj%get("nbc_col_elpa1_tridi_to_full", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    write(error_unit,*) "Problem setting option for non blocking collectives for cols in elpa1_tridi_to_full. Aborting..."
    call obj%timer%stop("trans_ev_&
    &real&
    &" // &
    &"_double" //&
    gpuString)
    success = .false.
    return
  endif

  if (non_blocking_collectives_rows .eq. 1) then
    useNonBlockingCollectivesRows = .true.
  else
    useNonBlockingCollectivesRows = .false.
  endif

  if (non_blocking_collectives_cols .eq. 1) then
    useNonBlockingCollectivesCols = .true.
  else
    useNonBlockingCollectivesCols = .false.
  endif

  my_prow = obj%mpi_setup%myRank_comm_rows
  my_pcol = obj%mpi_setup%myRank_comm_cols

  np_rows = obj%mpi_setup%nRanks_comm_rows
  np_cols = obj%mpi_setup%nRanks_comm_cols


  !call obj%timer%start("mpi_communication")
  !call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND) ,my_prowMPI, mpierr)
  !call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND) ,np_rowsMPI, mpierr)
  !call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND) ,my_pcolMPI, mpierr)
  !call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND) ,np_colsMPI, mpierr)

  !my_prow = int(my_prowMPI, kind=c_int)
  !np_rows = int(np_rowsMPI, kind=c_int)
  !my_pcol = int(my_pcolMPI, kind=c_int)
  !np_cols = int(np_colsMPI, kind=c_int)
  !call obj%timer%stop("mpi_communication")

  call obj%get("max_stored_rows",max_stored_rows_fac, error)

  totalblocks = (na-1)/nblk + 1
  max_blocks_row = (totalblocks-1)/np_rows + 1
  max_blocks_col = ((nqc-1)/nblk)/np_cols + 1  ! Columns of q_mat!

  max_local_rows = max_blocks_row*nblk
  max_local_cols = max_blocks_col*nblk

  max_stored_rows = (max_stored_rows_fac/nblk+1)*nblk
 
  if (.not.(useGPU)) then
    allocate(tmat(max_stored_rows,max_stored_rows), stat=istat, errmsg=errorMessage)
    call check_alloc("trans_ev_&
    &real&
    &", "tmat", istat, errorMessage)

    allocate(tmp1(max_local_cols*max_stored_rows), stat=istat, errmsg=errorMessage)
    call check_alloc("trans_ev_&
    &real&
    &", "tmp1", istat, errorMessage)

    allocate(tmp2(max_local_cols*max_stored_rows), stat=istat, errmsg=errorMessage)
    call check_alloc("trans_ev_&
    &real&
    &", "tmp2", istat, errorMessage)
  endif

  allocate(h1(max_stored_rows*max_stored_rows), stat=istat, errmsg=errorMessage)
  call check_alloc("trans_ev_&
  &real&
  &", "h1", istat, errorMessage)

  allocate(h2(max_stored_rows*max_stored_rows), stat=istat, errmsg=errorMessage)
  call check_alloc("trans_ev_&
  &real&
  &", "h2", istat, errorMessage)

  allocate(hvb(max_local_rows*nblk), stat=istat, errmsg=errorMessage)
  call check_alloc("trans_ev_&
  &real&
  &", "hvn", istat, errorMessage)

  allocate(hvm(max_local_rows,max_stored_rows), stat=istat, errmsg=errorMessage)
  call check_alloc("trans_ev_&
  &real&
  &", "hvm", istat, errorMessage)

  hvm = 0   ! Must be set to 0 !!!
  hvb = 0   ! Safety only
  blockStep = nblk

  l_cols = local_index(nqc, my_pcol, np_cols, nblk, -1) ! Local columns of q_mat

  nstor = 0
  if (useGPU) then
    hvn_ubnd = 0
  endif

 
  if (useGPU) then
    ! todo: this is used only for copying hmv to device.. it should be possible to go without it
    !allocate(hvm1(max_local_rows*max_stored_rows), stat=istat, errmsg=errorMessage)
    !call check_alloc("trans_ev_&
    !&real&
    !&", "hvm1", istat, errorMessage)
    successGPU = gpu_malloc(tmat_dev, max_stored_rows * max_stored_rows * size_of_datatype)
    call check_alloc_GPU_f("trans_ev", 362,  successGPU)

    successGPU = gpu_malloc(hvm_dev, max_local_rows * max_stored_rows * size_of_datatype)
    call check_alloc_GPU_f("trans_ev", 365,  successGPU)

    successGPU = gpu_malloc(tmp_dev, max_local_cols * max_stored_rows * size_of_datatype)
    call check_alloc_GPU_f("trans_ev", 368,  successGPU)

    num = ldq * matrixCols * size_of_datatype
    successGPU = gpu_malloc(q_dev, num)
    call check_alloc_GPU_f("trans_ev", 372,  successGPU)
  

    successGPU = gpu_memcpy(q_dev, int(loc(q_mat(1,1)),kind=c_intptr_t), &
                  num, gpuMemcpyHostToDevice)
    call check_memcpy_GPU_f("trans_ev", 394,  successGPU)
  endif  ! useGPU

  do istep = 1, na, blockStep


    ics = MAX(istep,3)
    ice = MIN(istep+nblk-1,na)
    if (ice<ics) cycle

    cur_pcol = pcol(istep, nblk, np_cols)

    nb = 0
    do ic = ics, ice

      l_colh = local_index(ic  , my_pcol, np_cols, nblk, -1) ! Column of Householder Vector
      l_rows = local_index(ic-1, my_prow, np_rows, nblk, -1) ! # rows of Householder Vector


      if (my_pcol == cur_pcol) then
        hvb(nb+1:nb+l_rows) = a_mat(1:l_rows,l_colh)
        if (my_prow == prow(ic-1, nblk, np_rows)) then
          hvb(nb+l_rows) = 1.
        endif
      endif

      nb = nb+l_rows
    enddo


    nb = 0
    do ic = ics, ice
      l_rows = local_index(ic-1, my_prow, np_rows, nblk, -1) ! # rows of Householder Vector
      hvm(1:l_rows,nstor+1) = hvb(nb+1:nb+l_rows)
      if (useGPU) then
        hvm_ubnd = l_rows
      endif
      nstor = nstor+1
      nb = nb+l_rows
    enddo

    ! Please note: for smaller matix sizes (na/np_rows<=256), a value of 32 for nstor is enough!
    if (nstor+nblk > max_stored_rows .or. istep+nblk > na .or. (na/np_rows <= 256 .and. nstor >= 32)) then

      ! Calculate scalar products of stored vectors.
      ! This can be done in different ways, we use dsyrk or zherk

      tmat = 0
      call obj%timer%start("blas")
      if (l_rows>0) then
        call DSYRK('U', 'T',   &
                         int(nstor,kind=BLAS_KIND), int(l_rows,kind=BLAS_KIND), ONE, &
                         hvm, int(ubound(hvm,dim=1),kind=BLAS_KIND), ZERO, tmat, int(max_stored_rows,kind=BLAS_KIND))
      endif
      call obj%timer%stop("blas")
      nc = 0
      do n = 1, nstor-1
        h1(nc+1:nc+n) = tmat(1:n,n+1)
        nc = nc+n
      enddo

      if (nc > 0) h2 = h1

      ! Calculate triangular matrix T

      nc = 0
      tmat(1,1) = tau(ice-nstor+1)
      do n = 1, nstor-1
        call obj%timer%start("blas")
        call DTRMV('L', 'T' , 'N', int(n,kind=BLAS_KIND), tmat, &
                            int(max_stored_rows,kind=BLAS_KIND), h2(nc+1), 1_BLAS_KIND)
        call obj%timer%stop("blas")

        tmat(n+1,1:n) = &
        -h2(nc+1:nc+n)  &
        *tau(ice-nstor+n+1)

        tmat(n+1,n+1) = tau(ice-nstor+n+1)
        nc = nc+n
      enddo

      if (useGPU) then
        ! todo: is this reshape really neccessary?
        hvm1(1:hvm_ubnd*nstor) = reshape(hvm(1:hvm_ubnd,1:nstor), (/ hvm_ubnd*nstor /))

        !hvm_dev(1:hvm_ubnd*nstor) = hvm1(1:hvm_ubnd*nstor)
        successGPU = gpu_memcpy(hvm_dev, int(loc(hvm1(1)),kind=c_intptr_t),   &
                      hvm_ubnd * nstor * size_of_datatype, gpuMemcpyHostToDevice)

        call check_memcpy_GPU_f("trans_ev", 546,  successGPU)

        !tmat_dev = tmat
        successGPU = gpu_memcpy(tmat_dev, int(loc(tmat(1,1)),kind=c_intptr_t),   &
                      max_stored_rows * max_stored_rows * size_of_datatype, gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("trans_ev", 551,  successGPU)
      endif


      ! Q = Q - V * T * V**T * Q

      if (l_rows > 0) then
        if (useGPU) then
          call obj%timer%start("gpublas")
          gpuHandle = obj%gpu_setup%gpublasHandleArray(0)
          call gpublas_DGEMM('T', 'N',   &
                                   nstor, l_cols, l_rows, ONE, hvm_dev, hvm_ubnd,  &
                                   q_dev, ldq, ZERO, tmp_dev, nstor, gpuHandle)
          call obj%timer%stop("gpublas")
        else ! useGPU

          call obj%timer%start("blas")
          call DGEMM('T', 'N',  &
                              int(nstor,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), &
                              int(l_rows,kind=BLAS_KIND), ONE, hvm, int(ubound(hvm,dim=1),kind=BLAS_KIND), &
                              q_mat, int(ldq,kind=BLAS_KIND), ZERO, tmp1, int(nstor,kind=BLAS_KIND))
          call obj%timer%stop("blas")
        endif ! useGPU

      else !l_rows>0

        if (useGPU) then
          if (gpu_vendor() /= OPENMP_OFFLOAD_GPU) then
            successGPU = gpu_memset(tmp_dev, 0, l_cols * nstor * size_of_datatype)
            call check_memcpy_GPU_f("trans_ev", 590,  successGPU)
          else
            allocate(tmp_debug(l_cols * nstor))
            tmp_debug(:) = 0.
            successGPU = gpu_memcpy(tmp_dev, int(loc(tmp_debug),kind=c_intptr_t), &
                                    l_cols*nstor*size_of_datatype, gpuMemcpyHostToDevice)
            call check_memcpy_GPU_f("trans_ev", 596,  successGPU)
            deallocate(tmp_debug)
          endif
        else
          tmp1(1:l_cols*nstor) = 0
        endif
      endif  !l_rows>0

!     tmp2 = tmp1

      if (l_rows > 0) then
        if (useGPU) then
          call obj%timer%start("gpublas")
          gpuHandle = obj%gpu_setup%gpublasHandleArray(0)
          call gpublas_DTRMM('L', 'L', 'N', 'N',     &
                                   nstor, l_cols, ONE, tmat_dev, max_stored_rows,  &
                                   tmp_dev, nstor, gpuHandle)

          call gpublas_DGEMM('N', 'N' ,l_rows ,l_cols ,nstor,  &
                                   -ONE, hvm_dev, hvm_ubnd, tmp_dev, nstor,   &
                                   ONE, q_dev, ldq, gpuHandle)
          call obj%timer%stop("gpublas")
        else !useGPU
          call obj%timer%start("blas")

          call DTRMM('L', 'L', 'N', 'N', int(nstor,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND),   &
                              ONE, tmat, int(max_stored_rows,kind=BLAS_KIND), tmp1, int(nstor,kind=BLAS_KIND))
          call DGEMM('N', 'N', int(l_rows,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), &
                              int(nstor,kind=BLAS_KIND), -ONE, hvm, int(ubound(hvm,dim=1),kind=BLAS_KIND), &
                              tmp1, int(nstor,kind=BLAS_KIND), ONE, q_mat, int(ldq,kind=BLAS_KIND))
          call obj%timer%stop("blas")
        endif ! useGPU
      endif  ! l_rows>0
      nstor = 0
    endif  ! (nstor+nblk>max_stored_rows .or. istep+nblk>na .or. (na/np_rows<=256 .and. nstor>=32))

    
  enddo ! istep = 1, na, blockStep

  deallocate(h1, h2, hvb, hvm, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_&     &MATH_DATATYPE&     &: h1, h2, hvb, hvm", 848,  istat,  errorMessage)

  if (useGPU) then

    !q_mat = q_dev
    successGPU = gpu_memcpy(int(loc(q_mat(1,1)),kind=c_intptr_t), &
                  q_dev, ldq * matrixCols * size_of_datatype, gpuMemcpyDeviceToHost)
    call check_memcpy_GPU_f("trans_ev", 864,  successGPU)

    !deallocate(hvm1, stat=istat, errmsg=errorMessage)
    !if (istat .ne. 0) then
    !  print *,"trans_ev_&
    !  &real&
    !  &: error when deallocating hvm1 "//errorMessage
    !  stop 1
    !endif

    !deallocate(q_dev, tmp_dev, hvm_dev, tmat_dev)
    successGPU = gpu_free(q_dev)
    call check_dealloc_GPU_f("trans_ev", 906,  successGPU)

    successGPU = gpu_free(tmp_dev)
    call check_dealloc_GPU_f("trans_ev", 909,  successGPU)

    successGPU = gpu_free(hvm_dev)
    call check_dealloc_GPU_f("trans_ev", 912,  successGPU)

    successGPU = gpu_free(tmat_dev)
    call check_dealloc_GPU_f("trans_ev", 915,  successGPU)
  else ! useGPU
    deallocate(tmat, tmp1, tmp2, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("trans_ev_&     &MATH_DATATYPE&     &: tmat, tmp1, tmp2", 920,  istat,  errorMessage)
  endif ! useGPU

!#if defined (USE_CCL_TRANS_EV) && !defined(WITH_GPU_STREAMS)
!  successGPU = cuda_stream_destroy(obj%gpu_setup%my_stream)
!  if (.not.(successGPU)) then
!    print *,"Cannot destroy gpu stream handle"
!  endif
!#endif

  call obj%timer%stop("trans_ev_&
  &real&
  &" // &
  &"_double" // &
  gpuString )

end subroutine trans_ev_&
&real&
&_&
&double


! now comes a dirty hack:
! the file elpa1_solve_tridi_real_template.F90 must be included twice
! for the legacy and for the new API. In the new API, however, some routines
! must be named "..._impl"

!#include "elpa1_solve_tridi_real_template.F90"




!cannot use "../src/elpa1/../solve_tridi/../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)


subroutine solve_tridi_&
&double_impl &
    ( obj, na, nev, d, e, q, ldq, nblk, matrixCols, mpi_comm_all, mpi_comm_rows, &
                                           mpi_comm_cols, useGPU, wantDebug, success, max_threads )

      use precision
      use elpa_abstract_impl
      use merge_recursive
      use merge_systems
      use elpa_mpi
      use ELPA_utilities
      use distribute_global_column
      use elpa_mpi
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
      integer(kind=ik), intent(in)               :: na, nev, ldq, nblk, matrixCols, &
                                                    mpi_comm_all, mpi_comm_rows, mpi_comm_cols
      real(kind=rk8), intent(inout)    :: d(na), e(na)
      real(kind=rk8), intent(inout)    :: q(ldq,*)
      logical, intent(in)                        :: useGPU, wantDebug
      logical, intent(out)                       :: success

      integer(kind=ik)                           :: i, j, n, np, nc, nev1, l_cols, l_rows
      integer(kind=ik)                           :: my_prow, my_pcol, np_rows, np_cols
      integer(kind=MPI_KIND)                     :: mpierr, my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI
      integer(kind=ik), allocatable              :: limits(:), l_col(:), p_col(:), l_col_bc(:), p_col_bc(:)

      integer(kind=ik)                           :: istat
      character(200)                             :: errorMessage
      character(20)                              :: gpuString
      integer(kind=ik), intent(in)               :: max_threads

      if(useGPU) then
        gpuString = "_gpu"
      else
        gpuString = ""
      endif

      call obj%timer%start("solve_tridi" // "_double" // gpuString)

      call obj%timer%start("mpi_communication")
      call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND) ,my_prowMPI, mpierr)
      call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND) ,np_rowsMPI, mpierr)
      call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND) ,my_pcolMPI, mpierr)
      call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND) ,np_colsMPI, mpierr)

      my_prow = int(my_prowMPI,kind=c_int)
      np_rows = int(np_rowsMPI,kind=c_int)
      my_pcol = int(my_pcolMPI,kind=c_int)
      np_cols = int(np_colsMPI,kind=c_int)

      call obj%timer%stop("mpi_communication")

      success = .true.

      l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a and q
      l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local columns of q

      ! Set Q to 0
      q(1:l_rows, 1:l_cols) = 0.0_rk

      ! Get the limits of the subdivisons, each subdivison has as many cols
      ! as fit on the respective processor column

      allocate(limits(0:np_cols), stat=istat, errmsg=errorMessage)
      call check_allocate_f("solve_tridi: limits", 128,  istat,  errorMessage)

      limits(0) = 0
      do np=0,np_cols-1
        nc = local_index(na, np, np_cols, nblk, -1) ! number of columns on proc column np

        ! Check for the case that a column has have zero width.
        ! This is not supported!
        ! Scalapack supports it but delivers no results for these columns,
        ! which is rather annoying
        if (nc==0) then
          call obj%timer%stop("solve_tridi" // "_double")
          if (wantDebug) write(error_unit,*) 'ELPA1_solve_tridi: ERROR: Problem contains processor column with zero width'
          success = .false.
          return
        endif
        limits(np+1) = limits(np) + nc
      enddo

      ! Subdivide matrix by subtracting rank 1 modifications

      do i=1,np_cols-1
        n = limits(i)
        d(n) = d(n)-abs(e(n))
        d(n+1) = d(n+1)-abs(e(n))
      enddo

      ! Solve sub problems on processsor columns

      nc = limits(my_pcol) ! column after which my problem starts

      if (np_cols>1) then
        nev1 = l_cols ! all eigenvectors are needed
      else
        nev1 = MIN(nev,l_cols)
      endif
      call solve_tridi_col_&
           &double_impl &
             (obj, l_cols, nev1, nc, d(nc+1), e(nc+1), q, ldq, nblk,  &
                        matrixCols, mpi_comm_rows, useGPU, wantDebug, success, max_threads)
      if (.not.(success)) then
        call obj%timer%stop("solve_tridi" // "_double" // gpuString)
        return
      endif
      ! If there is only 1 processor column, we are done

      if (np_cols==1) then
        deallocate(limits, stat=istat, errmsg=errorMessage)
        call check_deallocate_f("solve_tridi: limits", 176,  istat,  errorMessage)

        call obj%timer%stop("solve_tridi" // "_double" // gpuString)
        return
      endif

      ! Set index arrays for Q columns

      ! Dense distribution scheme:

      allocate(l_col(na), stat=istat, errmsg=errorMessage)
      call check_allocate_f("solve_tridi: l_col", 187,  istat,  errorMessage)

      allocate(p_col(na), stat=istat, errmsg=errorMessage)
      call check_allocate_f("solve_tridi: p_col", 190,  istat,  errorMessage)

      n = 0
      do np=0,np_cols-1
        nc = local_index(na, np, np_cols, nblk, -1)
        do i=1,nc
          n = n+1
          l_col(n) = i
          p_col(n) = np
        enddo
      enddo

      ! Block cyclic distribution scheme, only nev columns are set:

      allocate(l_col_bc(na), stat=istat, errmsg=errorMessage)
      call check_allocate_f("solve_tridi: l_col_bc", 205,  istat,  errorMessage)

      allocate(p_col_bc(na), stat=istat, errmsg=errorMessage)
      call check_allocate_f("solve_tridi: p_col_bc", 208,  istat,  errorMessage)

      p_col_bc(:) = -1
      l_col_bc(:) = -1

      do i = 0, na-1, nblk*np_cols
        do j = 0, np_cols-1
          do n = 1, nblk
            if (i+j*nblk+n <= MIN(nev,na)) then
              p_col_bc(i+j*nblk+n) = j
              l_col_bc(i+j*nblk+n) = i/np_cols + n
             endif
           enddo
         enddo
      enddo

      ! Recursively merge sub problems
      call merge_recursive_&
           &double &
           (obj, 0, np_cols, ldq, matrixCols, nblk, &
           l_col, p_col, l_col_bc, p_col_bc, limits, &
           np_cols, na, q, d, e, &
           mpi_comm_all, mpi_comm_rows, mpi_comm_cols,&
           useGPU, wantDebug, success, max_threads)

      if (.not.(success)) then
        call obj%timer%stop("solve_tridi" // "_double" // gpuString)
        return
      endif

      deallocate(limits,l_col,p_col,l_col_bc,p_col_bc, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("solve_tridi: limits, l_col, p_col, l_col_bc, p_col_bc", 239,  istat,  errorMessage)

      call obj%timer%stop("solve_tridi" // "_double" // gpuString)
      return

    end subroutine solve_tridi_&
        &double_impl

    subroutine solve_tridi_col_&
    &double_impl &
      ( obj, na, nev, nqoff, d, e, q, ldq, nblk, matrixCols, mpi_comm_rows, useGPU, wantDebug, success, max_threads )

   ! Solves the symmetric, tridiagonal eigenvalue problem on one processor column
   ! with the divide and conquer method.
   ! Works best if the number of processor rows is a power of 2!
      use precision
      use elpa_abstract_impl
      use elpa_mpi
      use merge_systems
      use ELPA_utilities
      use distribute_global_column
      implicit none
      class(elpa_abstract_impl_t), intent(inout) :: obj

      integer(kind=ik)              :: na, nev, nqoff, ldq, nblk, matrixCols, mpi_comm_rows
      real(kind=rk8)      :: d(na), e(na)
      real(kind=rk8)      :: q(ldq,*)

      integer(kind=ik), parameter   :: min_submatrix_size = 16 ! Minimum size of the submatrices to be used

      real(kind=rk8), allocatable    :: qmat1(:,:), qmat2(:,:)
      integer(kind=ik)              :: i, n, np
      integer(kind=ik)              :: ndiv, noff, nmid, nlen, max_size
      integer(kind=ik)              :: my_prow, np_rows
      integer(kind=MPI_KIND)        :: mpierr, my_prowMPI, np_rowsMPI

      integer(kind=ik), allocatable :: limits(:), l_col(:), p_col_i(:), p_col_o(:)
      logical, intent(in)           :: useGPU, wantDebug
      logical, intent(out)          :: success
      integer(kind=ik)              :: istat
      character(200)                :: errorMessage

      integer(kind=ik), intent(in)  :: max_threads

      integer(kind=MPI_KIND)        :: bcast_request1, bcast_request2
      logical                       :: useNonBlockingCollectivesRows
      integer(kind=c_int)           :: non_blocking_collectives, error

      success = .true.

      call obj%timer%start("solve_tridi_col" // "_double")

      call obj%get("nbc_row_solve_tridi", non_blocking_collectives, error)
      if (error .ne. ELPA_OK) then
        write(error_unit,*) "Problem setting option for non blocking collectives for rows in solve_tridi. Aborting..."
        success = .false.
        return
      endif

      if (non_blocking_collectives .eq. 1) then
        useNonBlockingCollectivesRows = .true.
      else
        useNonBlockingCollectivesRows = .false.
      endif

      call obj%timer%start("mpi_communication")
      call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND), my_prowMPI, mpierr)
      call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND), np_rowsMPI, mpierr)

      my_prow = int(my_prowMPI,kind=c_int)
      np_rows = int(np_rowsMPI,kind=c_int)
      call obj%timer%stop("mpi_communication")
      success = .true.
      ! Calculate the number of subdivisions needed.

      n = na
      ndiv = 1
      do while(2*ndiv<=np_rows .and. n>2*min_submatrix_size)
        n = ((n+3)/4)*2 ! the bigger one of the two halves, we want EVEN boundaries
        ndiv = ndiv*2
      enddo

      ! If there is only 1 processor row and not all eigenvectors are needed
      ! and the matrix size is big enough, then use 2 subdivisions
      ! so that merge_systems is called once and only the needed
      ! eigenvectors are calculated for the final problem.

      if (np_rows==1 .and. nev<na .and. na>2*min_submatrix_size) ndiv = 2

      allocate(limits(0:ndiv), stat=istat, errmsg=errorMessage)
      call check_deallocate_f("solve_tridi_col: limits", 445,  istat,  errorMessage)

      limits(0) = 0
      limits(ndiv) = na

      n = ndiv
      do while(n>1)
        n = n/2 ! n is always a power of 2
        do i=0,ndiv-1,2*n
          ! We want to have even boundaries (for cache line alignments)
          limits(i+n) = limits(i) + ((limits(i+2*n)-limits(i)+3)/4)*2
        enddo
      enddo

      ! Calculate the maximum size of a subproblem

      max_size = 0
      do i=1,ndiv
        max_size = MAX(max_size,limits(i)-limits(i-1))
      enddo

      ! Subdivide matrix by subtracting rank 1 modifications

      do i=1,ndiv-1
        n = limits(i)
        d(n) = d(n)-abs(e(n))
        d(n+1) = d(n+1)-abs(e(n))
      enddo

      if (np_rows==1)    then

        ! For 1 processor row there may be 1 or 2 subdivisions
        do n=0,ndiv-1
          noff = limits(n)        ! Start of subproblem
          nlen = limits(n+1)-noff ! Size of subproblem

          call solve_tridi_single_problem_&
          &double_impl &
                                  (obj, nlen,d(noff+1),e(noff+1), &
                                    q(nqoff+noff+1,noff+1),ubound(q,dim=1), wantDebug, success)

          if (.not.(success)) return
        enddo

      else

        ! Solve sub problems in parallel with solve_tridi_single
        ! There is at maximum 1 subproblem per processor

        allocate(qmat1(max_size,max_size), stat=istat, errmsg=errorMessage)
        call check_deallocate_f("solve_tridi_col: qmat1", 495,  istat,  errorMessage)

        allocate(qmat2(max_size,max_size), stat=istat, errmsg=errorMessage)
        call check_deallocate_f("solve_tridi_col: qmat2", 498,  istat,  errorMessage)

        qmat1 = 0 ! Make sure that all elements are defined

        if (my_prow < ndiv) then

          noff = limits(my_prow)        ! Start of subproblem
          nlen = limits(my_prow+1)-noff ! Size of subproblem
          call solve_tridi_single_problem_&
          &double_impl &
                                    (obj, nlen,d(noff+1),e(noff+1),qmat1, &
                                    ubound(qmat1,dim=1), wantDebug, success)

          if (.not.(success)) return
        endif

        ! Fill eigenvectors in qmat1 into global matrix q

        do np = 0, ndiv-1

          noff = limits(np)
          nlen = limits(np+1)-noff
!          qmat2 = qmat1 ! is this correct
          do i=1,nlen

            call distribute_global_column_&
            &double &
                     (obj, qmat1(1,i), q(1,noff+i), nqoff+noff, nlen, my_prow, np_rows, nblk)
          enddo

        enddo

        deallocate(qmat1, qmat2, stat=istat, errmsg=errorMessage)
        call check_deallocate_f("solve_tridi_col: qmat1, qmat2", 561,  istat,  errorMessage)

      endif

      ! Allocate and set index arrays l_col and p_col

      allocate(l_col(na), p_col_i(na),  p_col_o(na), stat=istat, errmsg=errorMessage)
      call check_deallocate_f("solve_tridi_col: l_col, p_col_i, p_col_o", 568,  istat,  errorMessage)

      do i=1,na
        l_col(i) = i
        p_col_i(i) = 0
        p_col_o(i) = 0
      enddo

      ! Merge subproblems

      n = 1
      do while(n<ndiv) ! if ndiv==1, the problem was solved by single call to solve_tridi_single

        do i=0,ndiv-1,2*n

          noff = limits(i)
          nmid = limits(i+n) - noff
          nlen = limits(i+2*n) - noff

          if (nlen == na) then
            ! Last merge, set p_col_o=-1 for unneeded (output) eigenvectors
            p_col_o(nev+1:na) = -1
          endif
          call merge_systems_&
          &double &
                              (obj, nlen, nmid, d(noff+1), e(noff+nmid), q, ldq, nqoff+noff, nblk, &
                               matrixCols, int(mpi_comm_rows,kind=ik), int(mpi_comm_self,kind=ik), &
                               l_col(noff+1), p_col_i(noff+1), &
                               l_col(noff+1), p_col_o(noff+1), 0, 1, useGPU, wantDebug, success, max_threads)
          if (.not.(success)) return

        enddo

        n = 2*n

      enddo

      deallocate(limits, l_col, p_col_i, p_col_o, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("solve_tridi_col: limits, l_col, p_col_i, p_col_o", 606,  istat,  errorMessage)

      call obj%timer%stop("solve_tridi_col" // "_double")

    end subroutine solve_tridi_col_&
    &double_impl

    subroutine solve_tridi_single_problem_&
    &double_impl &
    (obj, nlen, d, e, q, ldq, wantDebug, success)

   ! Solves the symmetric, tridiagonal eigenvalue problem on a single processor.
   ! Takes precautions if DSTEDC fails or if the eigenvalues are not ordered correctly.
     use precision
     use elpa_abstract_impl
     use elpa_blas_interfaces
     use ELPA_utilities
     implicit none
     class(elpa_abstract_impl_t), intent(inout) :: obj
     integer(kind=ik)                         :: nlen, ldq
     real(kind=rk8)                 :: d(nlen), e(nlen), q(ldq,nlen)

     real(kind=rk8), allocatable    :: work(:), qtmp(:), ds(:), es(:)
     real(kind=rk8)                 :: dtmp

     integer(kind=ik)              :: i, j, lwork, liwork, info
     integer(kind=BLAS_KIND)       :: infoBLAS
     integer(kind=ik), allocatable :: iwork(:)

     logical, intent(in)           :: wantDebug
     logical, intent(out)          :: success
      integer(kind=ik)             :: istat
      character(200)               :: errorMessage

     call obj%timer%start("solve_tridi_single" // "_double")

     success = .true.
     allocate(ds(nlen), es(nlen), stat=istat, errmsg=errorMessage)
     call check_allocate_f("solve_tridi_single: ds, es", 644,  istat,  errorMessage)

     ! Save d and e for the case that dstedc fails

     ds(:) = d(:)
     es(:) = e(:)

     ! First try dstedc, this is normally faster but it may fail sometimes (why???)

     lwork = 1 + 4*nlen + nlen**2
     liwork =  3 + 5*nlen
     allocate(work(lwork), iwork(liwork), stat=istat, errmsg=errorMessage)
     call check_allocate_f("solve_tridi_single: work, iwork", 656,  istat,  errorMessage)
     call obj%timer%start("lapack")
     call DSTEDC('I', int(nlen,kind=BLAS_KIND), d, e, q, int(ldq,kind=BLAS_KIND),    &
                          work, int(lwork,kind=BLAS_KIND), int(iwork,kind=BLAS_KIND), int(liwork,kind=BLAS_KIND), &
                          infoBLAS)
     info = int(infoBLAS,kind=ik)
     call obj%timer%stop("lapack")

     if (info /= 0) then

       ! DSTEDC failed, try DSTEQR. The workspace is enough for DSTEQR.

       write(error_unit,'(a,i8,a)') 'Warning: Lapack routine DSTEDC failed, info= ',info,', Trying DSTEQR!'

       d(:) = ds(:)
       e(:) = es(:)
       call obj%timer%start("lapack")
       call DSTEQR('I', int(nlen,kind=BLAS_KIND), d, e, q, int(ldq,kind=BLAS_KIND), work, infoBLAS )
       info = int(infoBLAS,kind=ik)
       call obj%timer%stop("lapack")

       ! If DSTEQR fails also, we don't know what to do further ...

       if (info /= 0) then
         if (wantDebug) then
           write(error_unit,'(a,i8,a)') 'ELPA1_solve_tridi_single: ERROR: Lapack routine DSTEQR failed, info= ',info,', Aborting!'
         endif
         success = .false.
         return
       endif
     end if

       deallocate(work,iwork,ds,es, stat=istat, errmsg=errorMessage)
       call check_deallocate_f("solve_tridi_single: work, iwork, ds, es", 689,  istat,  errorMessage)

      ! Check if eigenvalues are monotonically increasing
      ! This seems to be not always the case  (in the IBM implementation of dstedc ???)

      do i=1,nlen-1
        if (d(i+1)<d(i)) then
          if (abs(d(i+1) - d(i)) / abs(d(i+1) + d(i)) > 1e-14_rk8) then
            write(error_unit,'(a,i8,2g25.16)') '***WARNING: Monotony error dste**:',i+1,d(i),d(i+1)
          else
            write(error_unit,'(a,i8,2g25.16)') 'Info: Monotony error dste{dc,qr}:',i+1,d(i),d(i+1)
            write(error_unit,'(a)') 'The eigenvalues from a lapack call are not sorted to machine precision.'
            write(error_unit,'(a)') 'In this extent, this is completely harmless.'
            write(error_unit,'(a)') 'Still, we keep this info message just in case.'
          end if
          allocate(qtmp(nlen), stat=istat, errmsg=errorMessage)
          call check_allocate_f("solve_tridi_single: qtmp", 709,  istat,  errorMessage)

          dtmp = d(i+1)
          qtmp(1:nlen) = q(1:nlen,i+1)
          do j=i,1,-1
            if (dtmp<d(j)) then
              d(j+1)        = d(j)
              q(1:nlen,j+1) = q(1:nlen,j)
            else
              exit ! Loop
            endif
          enddo
          d(j+1)        = dtmp
          q(1:nlen,j+1) = qtmp(1:nlen)
          deallocate(qtmp, stat=istat, errmsg=errorMessage)
          call check_deallocate_f("solve_tridi_single: qtmp", 724,  istat,  errorMessage)

       endif
     enddo
     call obj%timer%stop("solve_tridi_single" // "_double")

    end subroutine solve_tridi_single_problem_&
    &double_impl

!#include "elpa1_merge_systems_real_template.F90"






subroutine hh_transform_real_&
&double &
(obj, alpha, xnorm_sq, xf, tau, wantDebug)
  ! Similar to LAPACK routine DLARFP, but uses ||x||**2 instead of x(:)
  ! and returns the factor xf by which x has to be scaled.
  ! It also hasn't the special handling for numbers < 1.d-300 or > 1.d150
  ! since this would be expensive for the parallel implementation.
  use precision
  use elpa_abstract_impl
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

  class(elpa_abstract_impl_t), intent(inout)    :: obj
  logical, intent(in)                           :: wantDebug
  real(kind=rk), intent(inout)       :: alpha
  real(kind=rk), intent(in)          :: xnorm_sq
  real(kind=rk), intent(out)         :: xf, tau

  real(kind=rk)                      :: BETA

  if (wantDebug) call obj%timer%start("hh_transform_&
                   &real&
     	      &" // &
                   &"_double" )


  if ( XNORM_SQ==0.0_rk ) then

    if ( ALPHA>=0.0_rk ) then
      TAU = 0.0_rk
    else
      TAU = 2.0_rk
      ALPHA = -ALPHA
    endif
    XF = 0.0_rk

  else

    BETA = SIGN( SQRT( ALPHA**2 + XNORM_SQ ), ALPHA )
    ALPHA = ALPHA + BETA
    IF ( BETA<0 ) THEN
      BETA = -BETA
      TAU  = -ALPHA / BETA
    ELSE
      ALPHA = XNORM_SQ / ALPHA

      TAU = ALPHA / BETA
      ALPHA = -ALPHA
    END IF
    XF = 1.0_rk/ALPHA
    ALPHA = BETA
  endif

  if (wantDebug) call obj%timer%stop("hh_transform_&
  &real&
  &" // &
  &"_double" )

end subroutine hh_transform_real_&
    &double



! real single precision























!cannot use "../src/elpa1/../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)










!> \brief Reduces a distributed symmetric matrix to tridiagonal form (like Scalapack Routine PDSYTRD)
!>
!  Parameters
!
!> \param obj	      object of elpa_type
!> \param na          Order of matrix
!>
!> \param a_mat(matrixRows,matrixCols)    Distributed matrix which should be reduced.
!>              Distribution is like in Scalapack.
!>              Opposed to PDSYTRD, a(:,:) must be set completely (upper and lower half)
!>              a(:,:) is overwritten on exit with the Householder vectors
!>
!> \param matrixRows         Leading dimension of a
!>
!> \param nblk        blocksize of cyclic distribution, must be the same in both directions!
!>
!> \param matrixCols  local columns of matrix
!>
!> \param mpi_comm_rows        MPI-Communicator for rows
!> \param mpi_comm_cols        MPI-Communicator for columns
!>
!> \param d_vec(na)       Diagonal elements (returned), identical on all processors
!>
!> \param e_vec(na)       Off-Diagonal elements (returned), identical on all processors
!>
!> \param tau(na)     Factors for the Householder vectors (returned), needed for back transformation
!>
!> \param useGPU      If true,  GPU version of the subroutine will be used
!> \param wantDebug   if true more debug information
!>
subroutine tridiag_&
  &real&
  &_&
  &single &
  (obj, na, a_mat, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols, d_vec, e_vec, tau, useGPU, wantDebug, &
   max_threads_in, isSkewsymmetric, success)
  use, intrinsic :: iso_c_binding
  use precision
  use elpa_abstract_impl
  use matrix_plot
  use elpa_omp
  use elpa_blas_interfaces
  use elpa_gpu
  use elpa_gpu_util
  use tridiag_gpu

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

  class(elpa_abstract_impl_t), intent(inout)    :: obj
  integer(kind=ik), intent(in)                  :: na, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols
  logical, intent(in)                           :: useGPU, wantDebug
  logical, intent(in)                           :: isSkewsymmetric

  logical                                       :: useCCL=.false.

  real(kind=rck), intent(out)          :: tau(na)
  real(kind=rck), intent(inout)        :: a_mat(matrixRows,*)
  real(kind=rk), intent(out)                    :: d_vec(na)
  real(kind=rk), intent(out)                    :: e_vec(na)
  integer(kind=ik)                              :: max_stored_uv = 32 ! TODO_23_11 - make it tunable instead of hard-coded
  logical,          parameter                   :: mat_vec_as_one_block = .true.

  ! id in processor row and column and total numbers of processor rows and columns
  integer(kind=ik)                              :: my_prow, my_pcol, np_rows, np_cols
  integer(kind=MPI_KIND)                        :: my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI
  integer(kind=MPI_KIND)                        :: mpierr
  integer(kind=ik)                              :: totalblocks, max_loc_block_rows, max_loc_block_cols, max_local_rows, &
                                                   max_local_cols
  ! updated after each istep (in the main cycle) to contain number of
  ! local columns and rows of the remaining part of the matrix
  integer(kind=ik)                              :: l_cols, l_rows
  integer(kind=ik)                              :: n_stored_vecs
  integer(kind=ik)                              :: isOurProcessRowInt ! TODO_23_11 - get rid of it
  logical                                       :: isOurProcessRow, isOurProcessCol, isOurProcessCol_prev


  integer(kind=C_intptr_T)                      :: a_dev, v_row_dev, v_col_dev, u_row_dev, u_col_dev, vu_stored_rows_dev, &
                                                   uv_stored_cols_dev, d_vec_dev, e_vec_dev, tau_dev
  logical                                       :: successGPU

  integer(kind=ik)                              :: istep, i, j, l_col_beg, l_col_end, l_row_beg, l_row_end
  integer(kind=ik)                              :: tile_size, l_rows_per_tile, l_cols_per_tile
  integer(kind=c_intptr_t)                      :: offset_dev

  integer(kind=ik), intent(in)                  :: max_threads_in
  integer(kind=ik)                              :: max_threads

  real(kind=rk)                                 :: vnorm2
  real(kind=rck)                       :: vav, x, aux1(2), vrl, xf, conjg_tau, dot_prod
  real(kind=rck), allocatable          :: aux(:)

  integer(kind=c_intptr_t)                      :: aux_dev, aux1_dev, aux_complex_dev, vav_dev, vav_host_or_dev, dot_prod_dev, & 
                                                   xf_dev, a_updated_element_dev, xf_host_or_dev, tau_istep_host_or_dev
  integer(kind=c_intptr_t)                      :: vnorm2_dev, vrl_dev ! alias pointers for aux1_dev(1), aux1_dev(2)


  integer(kind=c_intptr_t)                      :: num
  real(kind=rck), pointer              :: v_row(:), & ! used to store calculated Householder Vector
                                                   v_col(:)   ! the same Vector, but transposed 
  real(kind=rck), pointer              :: u_row_debug(:), & ! used to store calculated Householder Vector
                                                   u_col_debug(:)   ! the same Vector, but transposed 
  real(kind=rck), pointer              :: u_col(:), u_row(:)

  ! the following two matrices store pairs of vectors v and u calculated in each step
  ! at most max_stored_uv Vector pairs are stored, than the matrix A_i is explicitli updated
  ! u and v are stored both in row and Vector forms
  ! pattern: v1,u1,v2,u2,v3,u3,....
  ! todo: It is little bit confusing, I think, that variables _row actually store columns and vice versa
  !real(kind=rck), pointer             :: vu_stored_rows(:,:)
  real(kind=rck), allocatable         :: vu_stored_rows(:,:)
  ! pattern: u1,v1,u2,v2,u3,v3,....
  real(kind=rck), allocatable         :: uv_stored_cols(:,:)


  type(c_ptr)                                   :: v_row_host, v_col_host
  type(c_ptr)                                   :: u_row_host, u_col_host
  !type(c_ptr)                                   :: vu_stored_rows_host, uv_stored_cols_host
  real(kind=rk), allocatable                    :: tmp_real(:)
  integer(kind=ik)                              :: min_tile_size, error
  integer(kind=ik)                              :: istat
  character(200)                                :: errorMessage
  character(20)                                 :: gpuString
  integer(kind=ik)                              :: nblockEnd
  integer(kind=c_intptr_t), parameter           :: size_of_datatype = size_of_&
                                                                      &single&
                                                                      &_&
                                                                      &real
  integer(kind=c_intptr_t), parameter           :: size_of_datatype_real = size_of_&
                                                                      &single&
                                                                      &_real
  integer(kind=MPI_KIND)                        :: bcast_request1, bcast_request2, bcast_request3
  integer(kind=MPI_KIND)                        :: allreduce_request1, allreduce_request2, allreduce_request3
  integer(kind=MPI_KIND)                        :: allreduce_request4, allreduce_request5, allreduce_request6, &
                                                   allreduce_request7
  logical                                       :: useNonBlockingCollectivesCols
  logical                                       :: useNonBlockingCollectivesRows
  integer(kind=c_int)                           :: non_blocking_collectives_rows, non_blocking_collectives_cols
  logical                                       :: success

  integer(kind=c_intptr_t)                      :: gpuHandle, my_stream
  integer(kind=c_int) :: pointerMode

  integer(kind=ik)                              :: string_length, sm_count


  allocate(aux(2*max_stored_uv), stat=istat, errmsg=errorMessage)

  success = .true.

  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  if (useGPU) then
    max_threads=1
  else
    max_threads=max_threads_in
  endif

  call obj%timer%start("tridiag_&
  &real&
  &" // &
  "_single" // &
  gpuString )

  call obj%get("nbc_row_elpa1_full_to_tridi", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    write(error_unit,*) "Problem setting option for non blocking collectives for rows in elpa1_tridiag. Aborting..."
    success = .false.
    call obj%timer%stop("tridiag_&
    &real&
    &" // &
    "_single" // &
    gpuString )
    return
  endif

  call obj%get("nbc_col_elpa1_full_to_tridi", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    write(error_unit,*) "Problem setting option for non blocking collectives for cols in elpa1_tridiag. Aborting..."
    success = .false.
    call obj%timer%stop("tridiag_&
    &real&
    &" // &
    "_single" // &
    gpuString )
    return
  endif

  if (non_blocking_collectives_rows .eq. 1) then
    useNonBlockingCollectivesRows = .true.
  else
    useNonBlockingCollectivesRows = .false.
  endif

  if (non_blocking_collectives_cols .eq. 1) then
    useNonBlockingCollectivesCols = .true.
  else
    useNonBlockingCollectivesCols = .false.
  endif


  !mpi_comm_all    = obj%mpi_setup%mpi_comm_parent
  !mpi_comm_cols   = obj%mpi_setup%mpi_comm_cols
  !mpi_comm_rows   = obj%mpi_setup%mpi_comm_rows

  my_prow = obj%mpi_setup%myRank_comm_rows
  my_pcol = obj%mpi_setup%myRank_comm_cols

  np_rows = obj%mpi_setup%nRanks_comm_rows
  np_cols = obj%mpi_setup%nRanks_comm_cols


  !if (wantDebug) call obj%timer%start("mpi_communication")
  !call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND), my_prowMPI, mpierr)
  !call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND), np_rowsMPI, mpierr)
  !call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND), my_pcolMPI, mpierr)
  !call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND), np_colsMPI, mpierr)


  !my_prow = int(my_prowMPI, kind=c_int)
  !np_rows = int(np_rowsMPI, kind=c_int)
  !my_pcol = int(my_pcolMPI, kind=c_int)
  !np_cols = int(np_colsMPI, kind=c_int)
  !if (wantDebug) call obj%timer%stop("mpi_communication")

  ! Matrix is split into tiles; work is done only for tiles on the diagonal or above
  ! seems that tile is a square submatrix, consisting by several blocks
  ! it is a smallest possible square submatrix, where blocks being distributed among
  ! processors are "aligned" in both rows and columns
  !  -----------------
  ! | 1 4 | 1 4 | 1 4 | ...
  ! | 2 5 | 2 5 | 2 5 | ...
  ! | 3 6 | 3 6 | 3 6 | ...
  !  ----------------- ...
  ! | 1 4 | 1 4 | 1 4 | ...
  ! | 2 5 | 2 5 | 2 5 | ...
  ! | 3 6 | 3 6 | 3 6 | ...
  !  ----------------- .
  !   : :   : :   : :    .
  !   : :   : :   : :      .
  !
  ! this is a tile, where each number represents block, assigned to a processor with the shown number
  ! size of this small block is nblk
  ! Image is for situation with 6 processors, 3 processor rows and 2 columns
  ! tile_size is thus nblk * 6
  !
  tile_size = nblk*least_common_multiple(np_rows,np_cols) ! minimum global tile size

  ! make tile_size a smallest possible multiple of previously defined tile size, such that it is
  ! larger or equal to min_tile_size
  ! min_tile_size has been originally hardcoded as 128 * max(np_rows, np_cols), so it is now the implicit value
  ! it can, however, be set by the user
  call obj%get("min_tile_size", min_tile_size ,error)
  if (error .ne. ELPA_OK) then
    write(error_unit,*) "Problem setting option for min_tile_size. Aborting..."
    success = .false.
    call obj%timer%stop("tridiag_&
    &real&
    &" // &
    "_single" // &
    gpuString )
    return
  endif
  if(min_tile_size == 0) then
    ! not set by the user, use the default value
    min_tile_size = 128*max(np_rows, np_cols)
  endif
  tile_size = ((min_tile_size-1)/tile_size+1)*tile_size

  nblockEnd = 3

  l_rows_per_tile = tile_size/np_rows ! local rows of a tile
  l_cols_per_tile = tile_size/np_cols ! local cols of a tile

  totalblocks = (na-1)/nblk + 1
  max_loc_block_rows = (totalblocks-1)/np_rows + 1
  max_loc_block_cols = (totalblocks-1)/np_cols + 1

  ! localy owned submatrix has size at most max_local_rows x max_local_cols at each processor
  max_local_rows = max_loc_block_rows*nblk
  max_local_cols = max_loc_block_cols*nblk

  ! allocate memmory for vectors
  ! todo: It is little bit confusing, I think, that variables _row actually store columns and vice versa
  ! todo: if something has length max_local_rows, it is actually a column, no?
  ! todo: probably one should read it as v_row = Vector v distributed among rows

  allocate(uv_stored_cols(max_local_cols,2*max_stored_uv), stat=istat, errmsg=errorMessage)
  call check_alloc("tridiag_&
       &real ", "uv_stored_cols", istat, errorMessage)

  allocate(vu_stored_rows(max_local_rows,2*max_stored_uv), stat=istat, errmsg=errorMessage)
  call check_alloc("tridiag_&
       &real ", "vu_stored_rows", istat, errorMessage)

  if (useGPU) then

    ! allocate v_row 1 element longer to allow store and broadcast tau together with it
    if (gpu_vendor() /= OPENMP_OFFLOAD_GPU .and. gpu_vendor() /= SYCL_GPU) then
      num = (max_local_rows+1) * size_of_datatype
      successGPU = gpu_malloc_host(v_row_host, num)
      call check_host_alloc_GPU_f("tridiag: v_row_host", 444,  successGPU)
      call c_f_pointer(v_row_host,v_row,(/(max_local_rows+1)/))
    else
      allocate(v_row(max_local_rows+1))
    endif

    if (gpu_vendor() /= OPENMP_OFFLOAD_GPU .and. gpu_vendor() /= SYCL_GPU) then
      num = (max_local_cols) * size_of_datatype
      successGPU = gpu_malloc_host(v_col_host,num)
      call check_host_alloc_GPU_f("tridiag: v_col_host", 453,  successGPU)
      call c_f_pointer(v_col_host,v_col,(/(max_local_cols)/))
    else
      allocate(v_col(max_local_cols))
    endif

    if (gpu_vendor() /= OPENMP_OFFLOAD_GPU .and. gpu_vendor() /= SYCL_GPU) then
      num = (max_local_cols) * size_of_datatype
      successGPU = gpu_malloc_host(u_col_host,num)
      call check_host_alloc_GPU_f("tridiag: u_col_host", 462,  successGPU)
      call c_f_pointer(u_col_host,u_col,(/(max_local_cols)/))
    else
      allocate(u_col(max_local_cols))
    endif

    if (gpu_vendor() /= OPENMP_OFFLOAD_GPU .and. gpu_vendor() /= SYCL_GPU) then
      num = (max_local_rows) * size_of_datatype
      successGPU = gpu_malloc_host(u_row_host,num)
      call check_host_alloc_GPU_f("tridiag: u_row_host", 471,  successGPU)
      call c_f_pointer(u_row_host,u_row,(/(max_local_rows)/))
    else
      allocate(u_row(max_local_rows))
    endif

    
    if (gpu_vendor() /= OPENMP_OFFLOAD_GPU .and. gpu_vendor() /= SYCL_GPU) then
      num = (max_local_rows * 2*max_stored_uv) * size_of_datatype
      successGPU = gpu_host_register(int(loc(vu_stored_rows),kind=c_intptr_t), num, gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: vu_stored_rows", 481,  successGPU)

      num = (max_local_cols * 2*max_stored_uv) * size_of_datatype
      successGPU = gpu_host_register(int(loc(uv_stored_cols),kind=c_intptr_t), num, gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: uv_stored_cols", 485,  successGPU)

      num = (1 * 2*max_stored_uv) * size_of_datatype
      successGPU = gpu_host_register(int(loc(aux),kind=c_intptr_t), num, gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: aux", 489,  successGPU)

      num = na * size_of_datatype_real
      successGPU = gpu_host_register(int(loc(d_vec),kind=c_intptr_t), num, gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: d_vec", 493,  successGPU)

      num = na * size_of_datatype_real
      successGPU = gpu_host_register(int(loc(e_vec),kind=c_intptr_t), num, gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: e_vec", 497,  successGPU)

      num = na * size_of_datatype
      successGPU = gpu_host_register(int(loc(tau),kind=c_intptr_t), num, gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: tau", 501,  successGPU)

      num = 2 * size_of_datatype
      successGPU = gpu_host_register(int(loc(aux1),kind=c_intptr_t), num, gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: aux1", 505,  successGPU)

      num = 1 * size_of_datatype
      successGPU = gpu_host_register(int(loc(vav),kind=c_intptr_t), num, gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: vav", 509,  successGPU)

      num = 1 * size_of_datatype
      successGPU = gpu_host_register(int(loc(dot_prod),kind=c_intptr_t), num, gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: dot_prod", 513,  successGPU)

      num = 1 * size_of_datatype
      successGPU = gpu_host_register(int(loc(xf),kind=c_intptr_t), num, gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: xf", 517,  successGPU)
    endif ! gpu_vendor() /= OPENMP_OFFLOAD_GPU .and. gpu_vendor() /= SYCL_GPU
  else ! useGPU

    allocate(v_row(max_local_rows+1), stat=istat, errmsg=errorMessage)
    call check_alloc("tridiag_&
    &real ", "v_row", istat, errorMessage)

    allocate(v_col(max_local_cols), stat=istat, errmsg=errorMessage)
    call check_alloc("tridiag_&
     &real ", "v_col", istat, errorMessage)

    allocate(u_col(max_local_cols), stat=istat, errmsg=errorMessage)
    call check_alloc("tridiag_&
    &real ", "u_col", istat, errorMessage)

    allocate(u_row(max_local_rows), stat=istat, errmsg=errorMessage)
    call check_alloc("tridiag_&
    &real ", "u_row", istat, errorMessage)
      
  endif ! useGPU


  v_row = 0
  u_row = 0
  v_col = 0
  u_col = 0

  if (useGPU) then
    successGPU = gpu_malloc(v_row_dev, (max_local_rows+1) * size_of_datatype)
    call check_alloc_GPU_f("tridiag: v_row_dev", 556,  successGPU)

    successGPU = gpu_malloc(u_row_dev, max_local_rows * size_of_datatype)

    call check_alloc_GPU_f("tridiag: u_row_dev", 560,  successGPU)

    successGPU = gpu_malloc(v_col_dev, max_local_cols * size_of_datatype)
    call check_alloc_GPU_f("tridiag: v_col_dev", 563,  successGPU)

    successGPU = gpu_malloc(u_col_dev, max_local_cols * size_of_datatype)
    call check_alloc_GPU_f("tridiag: u_col_dev", 566,  successGPU)

    successGPU = gpu_malloc(vu_stored_rows_dev, max_local_rows * 2 * max_stored_uv * size_of_datatype)
    call check_alloc_GPU_f("tridiag: vu_stored_rows_dev", 569,  successGPU)

    successGPU = gpu_malloc(uv_stored_cols_dev, max_local_cols * 2 * max_stored_uv * size_of_datatype)
    call check_alloc_GPU_f("tridiag: uv_stored_cols_dev", 572,  successGPU)

    successGPU = gpu_malloc(d_vec_dev, na * size_of_datatype_real)
    call check_alloc_GPU_f("tridiag: d_vec_dev", 575,  successGPU)

    successGPU = gpu_malloc(e_vec_dev, na * size_of_datatype_real)
    call check_alloc_GPU_f("tridiag: e_vec_dev", 578,  successGPU)

    successGPU = gpu_malloc(tau_dev, na * size_of_datatype)
    call check_alloc_GPU_f("tridiag: tau_dev", 581,  successGPU)

    successGPU = gpu_malloc(aux_dev, 2*max_stored_uv * size_of_datatype)
    call check_alloc_GPU_f("tridiag: aux_dev", 584,  successGPU)

    successGPU = gpu_malloc(aux1_dev, 2 * size_of_datatype)
    call check_alloc_GPU_f("tridiag: aux1_dev", 587,  successGPU)

    successGPU = gpu_malloc(aux_complex_dev, 2 *max_stored_uv* size_of_datatype)
    call check_alloc_GPU_f("tridiag: aux_complex_dev", 590,  successGPU)

    successGPU = gpu_malloc(vav_dev, 1 * size_of_datatype)
    call check_alloc_GPU_f("tridiag: vav_dev", 593,  successGPU)

    successGPU = gpu_malloc(dot_prod_dev, 1 * size_of_datatype)
    call check_alloc_GPU_f("tridiag: dot_prod_dev", 596,  successGPU)

    successGPU = gpu_malloc(xf_dev, 1 * size_of_datatype)
    call check_alloc_GPU_f("tridiag: xf_dev", 599,  successGPU)

  endif !useGPU

  d_vec(:) = 0
  e_vec(:) = 0
  tau(:) = 0

  if (useGPU) then
    successGPU = gpu_memset(d_vec_dev, 0, na * size_of_datatype_real)
    call check_memcpy_GPU_f("tridiag: d_vec_dev", 630,  successGPU)

    successGPU = gpu_memset(e_vec_dev, 0, na * size_of_datatype_real)
    call check_memcpy_GPU_f("tridiag: e_vec_dev", 633,  successGPU)

    successGPU = gpu_memset(tau_dev, 0, na * size_of_datatype)
    call check_memcpy_GPU_f("tridiag: tau_dev", 636,  successGPU)
  endif

  n_stored_vecs = 0

  l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a_mat
  l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a_mat

  if (my_prow == prow(na, nblk, np_rows) .and. my_pcol == pcol(na, nblk, np_cols)) then
      d_vec(na) = a_mat(l_rows,l_cols)
  endif

  if (useGPU) then
    ! allocate memory for matrix A on the device and than copy the matrix

    num = matrixRows * matrixCols * size_of_datatype

    successGPU = gpu_malloc(a_dev, num)
    call check_alloc_GPU_f("tridiag: a_dev", 659,  successGPU)


    successGPU = gpu_memcpy(a_dev, int(loc(a_mat(1,1)),kind=c_intptr_t), &
                              num, gpuMemcpyHostToDevice)
    call check_memcpy_GPU_f("tridiag: a_dev", 678,  successGPU)

  endif ! useGPU

  ! main cycle of tridiagonalization
  ! in each step, 1 Householder Vector is calculated
  do istep = na, nblockEnd ,-1


    ! Calculate number of local rows and columns of the still remaining matrix
    ! on the local processor
    l_rows = local_index(istep-1, my_prow, np_rows, nblk, -1)
    l_cols = local_index(istep-1, my_pcol, np_cols, nblk, -1)

    ! Calculate Vector for Householder transformation on all procs
    ! owning column istep

    if (useGPU) then 
      ! copy l_cols + 1 column of a_dev to v_row_dev
      isOurProcessRow      = (my_prow == prow(istep-1, nblk, np_rows))
      isOurProcessCol      = (my_pcol == pcol(istep-1, nblk, np_cols))
      isOurProcessCol_prev = (my_pcol == pcol(istep  , nblk, np_cols)) ! isOurProcessCol from the previous step
      call gpu_copy_and_set_zeros_float(v_row_dev, a_dev, l_rows, l_cols, matrixRows, istep, &
                                            aux1_dev, vav_dev, d_vec_dev, &
                                            isOurProcessRow, isOurProcessCol, isOurProcessCol_prev, &
                                            isSkewsymmetric, useCCL, wantDebug, my_stream)
      if (gpu_vendor() /= SYCL_GPU) successGPU = gpu_DeviceSynchronize()
    endif ! useGPU

    if (my_pcol == pcol(istep, nblk, np_cols)) then

      ! Get Vector to be transformed; distribute last element and norm of
      ! remaining elements to all procs in current column

!             ! copy l_cols + 1 column of A to v_row
      ! if (useGPU) then

! #ifdef WITH_NVTX
!         call nvtxRangePush("memcpy new D-D a_dev(:,l_cols+1)->v_row_dev")
! #endif
!         ! TODO_23_11:  create a dev-dev copy kernel or merge it to another kernel
!         offset_dev = l_cols * matrixRows * size_of_datatype
!         successGPU = gpu_memcpy(v_row_dev, a_dev + offset_dev, (l_rows)* size_of_datatype, gpuMemcpyDeviceToDevice)
!         call check_memcpy_GPU_f("tridiag a_dev 1", 731,  successGPU)

! #ifdef WITH_NVTX
!         call nvtxRangePop()
! #endif

!       else ! useGPU
!         v_row(1:l_rows) = a_mat(1:l_rows,l_cols+1)
!       endif ! useGPU

      ! copy l_cols + 1 column of A to v_row
      if (.not. useGPU) then
        v_row(1:l_rows) = a_mat(1:l_rows,l_cols+1)
      endif ! useGPU

      if (n_stored_vecs > 0 .and. l_rows > 0) then
        if (useGPU) then
          if (wantDebug) call obj%timer%start("gpublas gemv skinny with copying")
          ! v_row_dev = vu_stored_rows_dev * uv_stored_cols_dev(l_cols+1,1:2*n_stored_vecs) + v_row_dev
          gpuHandle = obj%gpu_setup%gpublasHandleArray(0)
          call gpublas_SGEMV('N', l_rows, 2*n_stored_vecs,  &
                                    ONE, vu_stored_rows_dev, max_local_rows,  &
                                    uv_stored_cols_dev+(l_cols+1-1 +max_local_cols*(1-1))*size_of_datatype , max_local_cols, &
                                    ONE, v_row_dev, 1, gpuHandle)

          if (wantDebug .and. gpu_vendor() /= SYCL_GPU) successGPU = gpu_DeviceSynchronize()
          if (wantDebug) call obj%timer%stop("gpublas gemv skinny with copying")

        else ! useGPU

          if (wantDebug) call obj%timer%start("blas")
          ! v_row = vu_stored_rows * uv_stored_cols(l_cols+1,1:2*n_stored_vecs) + v_row
          call SGEMV('N',   &
                            int(l_rows,kind=BLAS_KIND), int(2*n_stored_vecs,kind=BLAS_KIND), &
                            ONE, vu_stored_rows, int(ubound(vu_stored_rows,dim=1),kind=BLAS_KIND), &
                            uv_stored_cols(l_cols+1,1), &
                            int(ubound(uv_stored_cols,dim=1),kind=BLAS_KIND), &
                            ONE, v_row, 1_BLAS_KIND)
          if (wantDebug) call obj%timer%stop("blas")
        endif ! useGPU
      endif ! (n_stored_vecs > 0 .and. l_rows > 0)

      if (useGPU) then
        if (my_prow == prow(istep-1, nblk, np_rows)) then
          isOurProcessRowInt = 1
        else
          isOurProcessRowInt = 0
        end if

        my_stream = obj%gpu_setup%my_stream
        call gpu_dot_product_and_assign_float(v_row_dev, l_rows, isOurProcessRowInt, aux1_dev, wantDebug, my_stream)
        if (.not. useCCL) then


          successGPU = gpu_memcpy(int(loc(aux1),kind=c_intptr_t), aux1_dev, 2*size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("tridiag: aux1_dev -> aux1", 819,  successGPU)


        endif ! .not. useCCL 
      endif ! useGPU  


      if (.not. useGPU) then
        if (my_prow == prow(istep-1, nblk, np_rows)) then
          aux1(1) = dot_product(v_row(1:l_rows-1),v_row(1:l_rows-1)) ! = "q"
          aux1(2) = v_row(l_rows) ! = "a_11" (or rather a_nn)
        else
          aux1(1) = dot_product(v_row(1:l_rows),v_row(1:l_rows))
          aux1(2) = 0.
        endif
      endif ! .not. useGPU 


      if (useCCL) then
      else ! useCCL
        vnorm2 = aux1(1)
        vrl    = aux1(2)
      endif ! useCCL


      ! Householder transformation
      if (useCCL) then
      else ! useCCL
        call hh_transform_real_&
                &single &
                (obj, vrl, vnorm2, xf, tau(istep), wantDebug)
      endif ! useCCL
      
      ! vrl is newly computed off-diagonal element of the final tridiagonal matrix
      if (my_prow == prow(istep-1, nblk, np_rows)) then
        if (.not. useCCL) then
          e_vec(istep-1) = vrl
        endif ! .not. useCCL
      endif

      if (.not. useGPU) then
        ! Scale v_row and store Householder Vector for back transformation
        v_row(1:l_rows) = v_row(1:l_rows) * xf

        if (my_prow == prow(istep-1, nblk, np_rows)) then
          v_row(l_rows) = 1.
        endif

        ! store Householder Vector for back transformation
        ! update a_mat
        a_mat(1:l_rows,l_cols+1) = v_row(1:l_rows)
      endif ! .not. useGPU 

      if (.not. useCCL) then
        ! add tau after the end of actuall v_row, to be broadcasted with it
        v_row(l_rows+1) = tau(istep)
      endif

      
      if (useGPU) then
        if (useCCL) then
          xf_host_or_dev = xf_dev
        else 
          xf_host_or_dev = int(loc(xf),kind=c_intptr_t)
        endif       
        
        isOurProcessRow = (my_prow == prow(istep-1, nblk, np_rows))
        call gpu_set_e_vec_scale_set_one_store_v_row_float(e_vec_dev, vrl_dev, a_dev, v_row_dev, tau_dev, xf_host_or_dev, & 
                                                  l_rows, l_cols, matrixRows, istep, isOurProcessRow, useCCL, wantDebug, my_stream)
      endif ! useGPU  


      if (useGPU .and. .not. useCCL) then      
        !v_row_dev -> v_row

        successGPU = gpu_memcpy(int(loc(v_row),kind=c_intptr_t), v_row_dev, l_rows*size_of_datatype, gpuMemcpyDeviceToHost)
        call check_memcpy_GPU_f("tridiag: v_row_dev -> v_row", 995,  successGPU)

      endif ! useGPU .and. .not. useCCL

    endif !(my_pcol == pcol(istep, nblk, np_cols))



    !recover tau, which has been broadcasted together with v_row
    if (.not. useCCL) then
      tau(istep) =  v_row(l_rows+1)
    endif

    ! Transpose Householder Vector v_row -> v_col
    if (useCCL) then
    else ! useCCL
      call elpa_transpose_vectors_&
          &real&
          &_&
          &single &
                (obj, v_row, ubound(v_row,dim=1), mpi_comm_rows, v_col, ubound(v_col,dim=1), mpi_comm_cols, &
                1, istep-1, 1, nblk, max_threads, .true., success)
      if (wantDebug .and. gpu_vendor() /= SYCL_GPU) successGPU = gpu_DeviceSynchronize()
      if (.not.(success)) then
        write(error_unit,*) "Error in elpa_transpose_vectors. Aborting!"
        return
      endif
    endif ! useCCL

    ! Calculate u = (A + VU**T + UV**T)*v // Dongarra 1987: "y = (A - UV**T - VU**T)*u"

    ! For cache efficiency, we use only the upper half of the matrix tiles for this,
    ! thus the result is partly in u_col(:) and partly in u_row(:)

    if (.not. useGPU) then
      u_col(1:l_cols) = 0
      u_row(1:l_rows) = 0
    endif

    if (useGPU .and. useCCL) then
      successGPU = gpu_memset(u_col_dev, 0, l_cols * size_of_datatype) ! TODO_23_11: omit this, but change gpublas_gemm to u_col_dev=a_dev^T*v_row_dev+0*u_col_dev?
      call check_memcpy_GPU_f("tridiag: u_col_dev", 1101,  successGPU)
    endif

    if (l_rows>0 .and. l_cols>0) then
      if (useGPU) then

        if (.not. mat_vec_as_one_block) then

          successGPU = gpu_memcpy(v_col_dev, int(loc(v_col(1)),kind=c_intptr_t), &
                        l_cols * size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("tridiag: v_col_dev", 1146,  successGPU)


        endif ! .not. mat_vec_as_one_block

        if (.not. useCCL) then

          successGPU = gpu_memcpy(v_row_dev, int(loc(v_row(1)),kind=c_intptr_t), &
                                    l_rows * size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("tridiag: v_row_dev", 1170,  successGPU)


        endif ! .not. useCCL
      endif ! useGPU


      if (.not. useGPU) then
        do i=0, (istep-2)/tile_size ! iteration over tiles
          l_col_beg = i*l_cols_per_tile+1
          l_col_end = min(l_cols,(i+1)*l_cols_per_tile)
          if (l_col_end < l_col_beg) cycle
          do j = 0, i
            l_row_beg = j*l_rows_per_tile+1
            l_row_end = min(l_rows,(j+1)*l_rows_per_tile)
            if (l_row_end < l_row_beg) cycle

            ! multiplication by blocks is efficient only for CPU
            ! for GPU we introduced 2 other ways, either by stripes (more simmilar to the original
            ! CPU implementation) or by one large matrix Vector multiply
            if (wantDebug) call obj%timer%start("blas")
            ! u_col = a_mat*v_row + u_col(=0)
            call SGEMV('T',  &
                        int(l_row_end-l_row_beg+1,kind=BLAS_KIND), int(l_col_end-l_col_beg+1,kind=BLAS_KIND), &
                        ONE, a_mat(l_row_beg, l_col_beg), int(matrixRows,kind=BLAS_KIND),         &
                        v_row(l_row_beg:max_local_rows+1), 1_BLAS_KIND,                           &
                        ONE, u_col(l_col_beg:max_local_cols), 1_BLAS_KIND)

            if (i/=j) then
              if (isSkewsymmetric) then
                call SGEMV('N',int(l_row_end-l_row_beg+1,kind=BLAS_KIND), int(l_col_end-l_col_beg+1,kind=BLAS_KIND), &
                                    -ONE, a_mat(l_row_beg,l_col_beg), int(matrixRows,kind=BLAS_KIND),               &
                                    v_col(l_col_beg:max_local_cols), 1_BLAS_KIND, ONE, u_row(l_row_beg:max_local_rows), &
                                    1_BLAS_KIND)

              else
                ! u_row = a_mat*v_col + u_row(=0)
                call SGEMV('N',int(l_row_end-l_row_beg+1,kind=BLAS_KIND), int(l_col_end-l_col_beg+1,kind=BLAS_KIND),  &
                                    ONE, a_mat(l_row_beg,l_col_beg), int(matrixRows,kind=BLAS_KIND),               &
                                    v_col(l_col_beg:max_local_cols), 1_BLAS_KIND, ONE, u_row(l_row_beg:max_local_rows), &
                                    1_BLAS_KIND)
              endif
            endif
            if (wantDebug) call obj%timer%stop("blas")

          enddo  ! j=0,i
        enddo  ! i=0,(istep-2)/tile_size
      endif ! .not. useGPU

      if (useGPU) then
        if (mat_vec_as_one_block) then
          ! Unlike for CPU, we (for each MPI thread) do just one large mat-vec multiplication
          ! this requires altering of the algorithm when later explicitly updating the matrix
          ! after max_stored_uv is reached : we need to update all tiles, not only those above diagonal
          if (wantDebug) call obj%timer%start("gpublas")
          
          gpuHandle = obj%gpu_setup%gpublasHandleArray(0) ! 0.3-0.5 us

              
          ! u_col_dev = a_dev^T*v_row_dev
          call gpublas_SGEMV('T', l_rows,l_cols,  &
                                    ONE, a_dev, matrixRows,                   &
                                    v_row_dev , 1,                          &
                                    ZERO, u_col_dev, 1, gpuHandle)
              
          if (wantDebug .and. gpu_vendor() /= SYCL_GPU) successGPU = gpu_DeviceSynchronize()
              
       ! todo: try with non transposed!!!
!                 if(i/=j) then
!                   call gpublas_SGEMV('N', l_row_end-l_row_beg+1,l_col_end-l_col_beg+1,  &
!                                             ONE, a_dev + offset_dev, matrixRows,                        &
!                                             v_col_dev + (l_col_beg - 1) *                      &
!                                             size_of_datatype, 1,                          &
!                                             ONE, u_row_dev + (l_row_beg - 1) *                 &
!                                             size_of_datatype, 1)
!                 endif
          if (wantDebug) call obj%timer%stop("gpublas")
        else  ! mat_vec_as_one_block
          !perform multiplication by stripes - it is faster than by blocks, since we call cublas with
          !larger matrices. In general, however, this algorithm is very simmilar to the one with CPU
          do i=0,(istep-2)/tile_size
            l_col_beg = i*l_cols_per_tile+1
            l_col_end = min(l_cols,(i+1)*l_cols_per_tile)
            if (l_col_end<l_col_beg) cycle

            l_row_beg = 1
            l_row_end = min(l_rows,(i+1)*l_rows_per_tile)
                  
            offset_dev = ((l_row_beg-1) + (l_col_beg - 1) * matrixRows) * size_of_datatype
  
            gpuHandle = obj%gpu_setup%gpublasHandleArray(0)
            call gpublas_SGEMV('T', &
                          l_row_end-l_row_beg+1, l_col_end-l_col_beg+1, &
                          ONE, a_dev + offset_dev, matrixRows,  &
                          v_row_dev + (l_row_beg - 1) * size_of_datatype, 1,  &
                          ONE, u_col_dev + (l_col_beg - 1) * size_of_datatype, 1, gpuHandle)
          enddo !i=0,(istep-2)/tile_size

          do i=0,(istep-2)/tile_size
            l_col_beg = i*l_cols_per_tile+1
            l_col_end = min(l_cols,(i+1)*l_cols_per_tile)
            if (l_col_end<l_col_beg) cycle

            l_row_beg = 1
            l_row_end = min(l_rows,i*l_rows_per_tile)
              
            offset_dev = ((l_row_beg-1) + (l_col_beg - 1) * matrixRows) * size_of_datatype

            if (isSkewsymmetric) then
                gpuHandle = obj%gpu_setup%gpublasHandleArray(0)
                call gpublas_SGEMV('N', l_row_end-l_row_beg+1, l_col_end-l_col_beg+1, &
                            -ONE, a_dev + offset_dev, matrixRows, &
                            v_col_dev + (l_col_beg - 1) * size_of_datatype,1, &
                            ONE, u_row_dev + (l_row_beg - 1) * size_of_datatype, 1, gpuHandle)
            else
                gpuHandle = obj%gpu_setup%gpublasHandleArray(0)
                call gpublas_SGEMV('N', l_row_end-l_row_beg+1, l_col_end-l_col_beg+1, &
                            ONE, a_dev + offset_dev, matrixRows, &
                            v_col_dev + (l_col_beg - 1) * size_of_datatype,1, &
                            ONE, u_row_dev + (l_row_beg - 1) * size_of_datatype, 1, gpuHandle)
            endif
          enddo ! i=0,(istep-2)/tile_size
        end if ! mat_vec_as_one_block / per stripes


        if (.not. mat_vec_as_one_block) then

          successGPU = gpu_memcpy(int(loc(u_row(1)),kind=c_intptr_t), u_row_dev, l_rows*size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("tridiag: u_row_dev 1", 1378,  successGPU)


        endif ! .not. mat_vec_as_one_block
      endif ! useGPU


      ! second calculate (VU**T + UV**T)*v part of (A + VU**T + UV**T)*v
      if (n_stored_vecs > 0) then
        if (.not. useGPU) then
          if (wantDebug) call obj%timer%start("blas")

      
          ! aux = vu_stored_rows^T*v_row
          call SGEMV('T',     &
                              int(l_rows,kind=BLAS_KIND), int(2*n_stored_vecs,kind=BLAS_KIND),   &
                              ONE, vu_stored_rows, int(ubound(vu_stored_rows,dim=1),kind=BLAS_KIND),   &
                              v_row,  1_BLAS_KIND, ZERO, aux, 1_BLAS_KIND)
              
          ! u_col = uv_stored_cols*aux + u_col
          call SGEMV('N', int(l_cols,kind=BLAS_KIND), int(2*n_stored_vecs,kind=BLAS_KIND),   &
                              ONE, uv_stored_cols, int(ubound(uv_stored_cols,dim=1),kind=BLAS_KIND),   &
                              aux, 1_BLAS_KIND, ONE, u_col,  1_BLAS_KIND)
          if (wantDebug) call obj%timer%stop("blas")
        else ! .not. useGPU
          if (wantDebug) call obj%timer%start("gpublas gemv x2 skinny")


          ! aux_dev = vu_stored_rows_dev^T*v_row_dev
          gpuHandle = obj%gpu_setup%gpublasHandleArray(0)
          call gpublas_SGEMV('T', l_rows, 2*n_stored_vecs, &
                                      ONE, vu_stored_rows_dev, max_local_rows,   &
                                      v_row_dev,  1, ZERO, aux_dev, 1, gpuHandle)
              
          ! u_col_dev = uv_stored_cols_dev*aux_dev + u_col_dev
          call gpublas_SGEMV('N', l_cols, 2*n_stored_vecs, &
                                      ONE, uv_stored_cols_dev, max_local_cols,   &
                                      aux_dev, 1, ONE, u_col_dev, 1, gpuHandle)

          if (wantDebug .and. gpu_vendor() /= SYCL_GPU) successGPU = gpu_DeviceSynchronize()

          if (wantDebug) call obj%timer%stop("gpublas gemv x2 skinny")  
        endif ! .not. useGPU
      endif ! n_stored_vecs > 0

    endif  ! (l_rows>0 .and. l_cols>0)

    if (useGPU .and. l_cols>0 .and. (.not. useCCL)) then

      successGPU = gpu_memcpy(int(loc(u_col(1)),kind=c_intptr_t), u_col_dev, l_cols*size_of_datatype, gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("tridiag: u_col_dev 1", 1462,  successGPU)

    endif ! useGPU


    if (useGPU .and. .not. useCCL) then 
        if (l_rows==0) then
        u_col(1:l_cols) = 0
      endif
    endif ! useGPU

    ! Sum up all u_row(:) parts along rows and add them to the u_col(:) parts
    ! on the processors containing the diagonal
    ! This is only necessary if u_row has been calculated, i.e. if the
    ! global tile size is smaller than the global remaining matrix
    
    ! in GPU case, u_row_dev=0 up to now, so elpa_reduce_add_vectors is skipped

    if (tile_size < istep-1 .and. .not. useGPU) then
      call elpa_reduce_add_vectors_&
            &real&
            &_&
            &single &
            (obj, u_row, ubound(u_row,dim=1), mpi_comm_rows, u_col, ubound(u_col,dim=1), &
            mpi_comm_cols, istep-1, 1, nblk, max_threads)
    endif ! (tile_size < istep-1 .and. .not. useGPU)

    ! Sum up all the u_col(:) parts, transpose u_col -> u_row

    if (l_cols > 0) then
    endif ! (l_cols > 0)
    
    ! Transpose Householder Vector u_col -> u_row
    if (useCCL) then
    else ! useCCL
      if (isSkewsymmetric) then
        call elpa_transpose_vectors_ss_&
        &real&
        &_&
        &single &
        (obj, u_col, ubound(u_col,dim=1), mpi_comm_cols, u_row, ubound(u_row,dim=1), &
          mpi_comm_rows, 1, istep-1, 1, nblk, max_threads, .false., success)
        if (.not.(success)) then
          write(error_unit,*) "Error in elpa_transpose_vectors_ss. Aborting!"
          return
        endif
      else ! isSkewsymmetric
        call elpa_transpose_vectors_&
        &real&
        &_&
        &single &
        (obj, u_col, ubound(u_col,dim=1), mpi_comm_cols, u_row, ubound(u_row,dim=1), &
        mpi_comm_rows, 1, istep-1, 1, nblk, max_threads, .false., success)
        if (.not.(success)) then
          write(error_unit,*) "Error in elpa_transpose_vectors. Aborting!"
          return
        endif
      endif ! isSkewsymmetric
    endif ! useCCL


    if (useGPU .and. .not. useCCL) then

      successGPU = gpu_memcpy(v_col_dev, int(loc(v_col(1)),kind=c_intptr_t), l_cols*size_of_datatype, gpuMemcpyHostToDevice)
      call check_memcpy_GPU_f("tridiag: v_col_dev", 1607,  successGPU)


      successGPU = gpu_memcpy(u_col_dev, int(loc(u_col(1)),kind=c_intptr_t), l_cols*size_of_datatype, gpuMemcpyHostToDevice)
      call check_memcpy_GPU_f("tridiag: u_col_dev", 1622,  successGPU)

    endif ! (useGPU .and. .not. useCCL)



    if (useGPU .and. .not. useCCL) then

      successGPU = gpu_memcpy(u_row_dev, int(loc(u_row(1)),kind=c_intptr_t), l_rows*size_of_datatype, gpuMemcpyHostToDevice)
      call check_memcpy_GPU_f("tridiag: u_row_dev", 1657,  successGPU)


      successGPU = gpu_memcpy(v_row_dev, int(loc(v_row(1)),kind=c_intptr_t), l_rows*size_of_datatype, gpuMemcpyHostToDevice)
      call check_memcpy_GPU_f("tridiag: v_row_dev", 1671,  successGPU)

    endif ! (useGPU .and. .not. useCCL)


    ! calculate u**T * v (same as v**T * (A + VU**T + UV**T) * v )
    if (.not. useGPU .or. (useGPU .and. .not. useCCL)) then
      vav = 0 ! x=0
      if (l_cols>0) vav = dot_product(v_col(1:l_cols), u_col(1:l_cols))
    endif



    ! store u and v in the matrices U and V
    ! these matrices are stored combined in one here

   if (.not. useCCL) then
      conjg_tau = tau(istep)
    endif ! .not. useCCL
    
    if (.not. useGPU) then
      if (l_rows > 0) then
        ! update vu_stored_rows
        vu_stored_rows(1:l_rows,2*n_stored_vecs+1) = conjg_tau*v_row(1:l_rows)
        vu_stored_rows(1:l_rows,2*n_stored_vecs+2) = 0.5*conjg_tau*vav*v_row(1:l_rows) - u_row(1:l_rows)
      endif
      if (l_cols > 0) then
        ! update uv_stored_cols
        uv_stored_cols(1:l_cols,2*n_stored_vecs+1) = 0.5*conjg_tau*vav*v_col(1:l_cols) - u_col(1:l_cols)
        uv_stored_cols(1:l_cols,2*n_stored_vecs+2) = conjg_tau*v_col(1:l_cols)
      endif
    endif ! .not. useGPU


    if (useGPU) then

      ! kernel: update vu_stored_rows_dev, uv_stored_cols
      ! then cpu's "store u,v in U,V" can be deleted. But we should take care of dot_prod below, where vu_stored_rows and uv_stored_cols are used
      if (useCCL) then
        vav_host_or_dev = vav_dev
        !tau_istep_host_or_dev = tau_dev + (istep-1)*size_of_datatype
        tau_istep_host_or_dev = v_row_dev + (l_rows+1-1)*size_of_datatype
      else ! useCCL
        vav_host_or_dev = int(loc(vav),kind=c_intptr_t)
        tau_istep_host_or_dev = int(loc(tau(istep)), kind=c_intptr_t)
      endif ! useCCL
      call gpu_store_u_v_in_uv_vu_float(vu_stored_rows_dev, uv_stored_cols_dev, v_row_dev, u_row_dev, &
                                          v_col_dev, u_col_dev, tau_dev, aux_complex_dev, &
                                          vav_host_or_dev, tau_istep_host_or_dev, &
                                          l_rows, l_cols, n_stored_vecs,  max_local_rows, max_local_cols, istep, &
                                          useCCL, wantDebug, my_stream)
    endif ! useGPU


    ! We have calculated another Householder Vector, number of implicitly stored increased
    n_stored_vecs = n_stored_vecs+1

    ! If the limit of max_stored_uv is reached, calculate A + VU**T + UV**T
    if (n_stored_vecs == max_stored_uv .or. istep == 3) then        
      
      if (.not. useGPU .OR. .not. mat_vec_as_one_block) then
        do i = 0, (istep-2)/tile_size
          ! go over tiles above (or on) the diagonal
          l_col_beg = i*l_cols_per_tile+1
          l_col_end = min(l_cols,(i+1)*l_cols_per_tile)
          l_row_beg = 1
          l_row_end = min(l_rows,(i+1)*l_rows_per_tile)
          if (l_col_end<l_col_beg .or. l_row_end<l_row_beg) then
            cycle
          endif

          if (useGPU) then
            if (.not. mat_vec_as_one_block) then
              ! if using mat-vec multiply by stripes, it is enough to update tiles above (or on) the diagonal only
              ! we than use the same calls as for CPU version
              if (wantDebug) call obj%timer%start("gpublas_gemm")
              ! update a_dev
              ! a_dev = vu_stored_rows_dev*uv_stored_cols_dev + a_dev
              gpuHandle = obj%gpu_setup%gpublasHandleArray(0)
              call gpublas_SGEMM('N', 'T',     &
                                      l_row_end-l_row_beg+1, l_col_end-l_col_beg+1, 2*n_stored_vecs,                      &
                                      ONE, vu_stored_rows_dev + (l_row_beg - 1) *                                         &
                                      size_of_datatype,  &
                                      max_local_rows, uv_stored_cols_dev + (l_col_beg - 1) *                              &
                                      size_of_datatype,  &
                                      max_local_cols, ONE, a_dev + ((l_row_beg - 1) + (l_col_beg - 1) * matrixRows) *     &
                                      size_of_datatype , matrixRows, gpuHandle)
              if (wantDebug .and. gpu_vendor() /= SYCL_GPU) successGPU = gpu_DeviceSynchronize()
              if (wantDebug) call obj%timer%stop("gpublas_gemm")
            endif ! .not. mat_vec_as_one_block
          else ! useGPU
            if (wantDebug) call obj%timer%start("blas_gemm")
            ! update a_mat
            ! a_mat = vu_stored_rows*uv_stored_cols + a_mat
            call SGEMM('N', 'T',                &
                                int(l_row_end-l_row_beg+1,kind=BLAS_KIND), int(l_col_end-l_col_beg+1,kind=BLAS_KIND), &
                                int(2*n_stored_vecs,kind=BLAS_KIND),    &
                                ONE, vu_stored_rows(l_row_beg:max_local_rows,1:2*max_stored_uv), &
                                int(ubound(vu_stored_rows,dim=1),kind=BLAS_KIND),   &
                                uv_stored_cols(l_col_beg,1), &
                                int(ubound(uv_stored_cols,dim=1),kind=BLAS_KIND),        &
                                ONE, a_mat(l_row_beg,l_col_beg), int(matrixRows,kind=BLAS_KIND))
            if (wantDebug) call obj%timer%stop("blas_gemm")
          endif !useGPU
        enddo ! i = 0, (istep-2)/tile_size

      else !.not. useGPU or .not. mat_vec_as_one_block (i.e. useGPU and mat_vec_as_one_block)

        !update whole (remaining) part of matrix, including tiles below diagonal
        !we can do that in one large cublas call
        if (wantDebug) call obj%timer%start("gpublas_gemm")


        gpuHandle = obj%gpu_setup%gpublasHandleArray(0)
        ! update a_dev
        ! a_dev = vu_stored_rows_dev*uv_stored_cols_dev + a_dev
        call gpublas_SGEMM('N', 'T', l_rows, l_cols, 2*n_stored_vecs,   &
                                  ONE, vu_stored_rows_dev, max_local_rows, &
                                  uv_stored_cols_dev, max_local_cols,  &
                                  ONE, a_dev, matrixRows, gpuHandle)

        if (wantDebug .and. gpu_vendor() /= SYCL_GPU) successGPU = gpu_DeviceSynchronize()
        if (wantDebug) call obj%timer%stop("gpublas_gemm")
        
        ! copy real(a_dev(l_rows,l_cols)) -> d_vec_dev(istep-1) for correct initial value of d_vec_dev before atomicAdd
        if ((.not. isSkewsymmetric) .and. &
            (my_prow == prow(istep-1, nblk, np_rows) .and. my_pcol == pcol(istep-1, nblk, np_cols))) then
          offset_dev = ((l_rows-1) + (l_cols-1)*matrixRows) * size_of_datatype
          successGPU = gpu_memcpy(d_vec_dev + (istep-2)*size_of_datatype_real, &
                                  a_dev + offset_dev, 1*size_of_datatype_real, gpuMemcpyDeviceToDevice)
          call check_memcpy_GPU_f("tridiag a_dev->d_vec_dev", 1863,  successGPU)
        endif

      endif !.not. useGPU or .not. mat_vec_as_one_block

      n_stored_vecs = 0
    endif ! (n_stored_vecs == max_stored_uv .or. istep == 3)

    if (my_prow == prow(istep-1, nblk, np_rows) .and. my_pcol == pcol(istep-1, nblk, np_cols)) then

      if (useGPU) then
        my_stream = obj%gpu_setup%my_stream
        call gpu_update_matrix_element_add_float(vu_stored_rows_dev, uv_stored_cols_dev, a_dev, d_vec_dev, &
                                            l_rows, l_cols, matrixRows, max_local_rows, max_local_cols, istep, n_stored_vecs, &
                                            isSkewsymmetric, wantDebug, my_stream)

      else ! useGPU
        if (n_stored_vecs > 0) then
          ! update a_mat (only one elememt!)
          dot_prod = dot_product(vu_stored_rows(l_rows,1:2*n_stored_vecs), uv_stored_cols(l_cols,1:2*n_stored_vecs))
          a_mat(l_rows,l_cols) = a_mat(l_rows,l_cols) + dot_prod
        endif
        if (isSkewsymmetric) then
          d_vec(istep-1) = 0.0_rk
        else
          d_vec(istep-1) = a_mat(l_rows,l_cols)
        endif
      endif ! useGPU

    endif ! (my_prow == prow(istep-1, nblk, np_rows) .and. my_pcol == pcol(istep-1, nblk, np_cols))


  enddo ! main cycle over istep=na,3,-1


  if (useGPU) then
    ! copy a_dev -> a_mat for backtransformation
    num = matrixRows * matrixCols * size_of_datatype
    successGPU = gpu_memcpy(int(loc(a_mat(1,1)),kind=c_intptr_t), a_dev, num, gpuMemcpyDeviceToHost)
    call check_memcpy_GPU_f("tridiag: a_dev", 1921,  successGPU)

  endif ! useGPU



  ! Store e_vec(1)

  if (my_prow==prow(1, nblk, np_rows) .and. my_pcol==pcol(2, nblk, np_cols)) then
    if (useGPU) then
      successGPU = gpu_memcpy(int(loc(e_vec(1)),kind=c_intptr_t), a_dev + (matrixRows * (l_cols - 1)) * size_of_datatype, &
                              1 * size_of_datatype, gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("tridiag: a_dev 7", 2026,  successGPU)
    else !useGPU
      e_vec(1) = a_mat(1,l_cols) ! use last l_cols value of loop above
    endif !useGPU
  endif ! if (my_prow==prow(1, nblk, np_rows) .and. my_pcol==pcol(2, nblk, np_cols))

  ! Store d_vec(1)
  if (my_prow==prow(1, nblk, np_rows) .and. my_pcol==pcol(1, nblk, np_cols)) then
    if(useGPU) then
      successGPU = gpu_memcpy(int(loc(d_vec(1)),kind=c_intptr_t), a_dev, 1 * size_of_datatype, gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("tridiag: a_dev 8", 2045,  successGPU)
    else !useGPU
      if (isSkewsymmetric) then
        d_vec(1) = 0.0_rk
      else
        d_vec(1) = a_mat(1,1)
      endif
    endif !useGPU
  endif ! (my_prow==prow(1, nblk, np_rows) .and. my_pcol==pcol(1, nblk, np_cols))

  if (useGPU) then
    offset_dev = 1 * size_of_datatype_real
    ! first and last elements of d_vec are treated separately
    successGPU = gpu_memcpy(int(loc(d_vec(2)),kind=c_intptr_t), &
                            d_vec_dev + offset_dev, (na-2) * size_of_datatype_real, gpuMemcpyDeviceToHost)
    call check_memcpy_GPU_f("tridiag: d_vec", 2062,  successGPU)

    if (useCCL) then
      ! e_vec(1) is treated separately
      offset_dev = 1 * size_of_datatype_real
      successGPU = gpu_memcpy(int(loc(e_vec(2)),kind=c_intptr_t), &
                              e_vec_dev + offset_dev, (na-1) * size_of_datatype_real, gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("tridiag: e_vec", 2069,  successGPU)

      ! tau(2) is treated separately, tau(1) is not used
      offset_dev = 2 * size_of_datatype
      successGPU = gpu_memcpy(int(loc(tau(3)),kind=c_intptr_t), &
                              tau_dev + offset_dev, (na-2) * size_of_datatype, gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("tridiag: tau", 2075,  successGPU)
    endif

    ! todo: should we leave a_mat on the device for further use?
    successGPU = gpu_free(a_dev)
    call check_dealloc_GPU_f("tridiag: a_dev 9", 2080,  successGPU)

    successGPU = gpu_free(v_row_dev)
    call check_dealloc_GPU_f("tridiag: v_row_dev", 2083,  successGPU)

    successGPU = gpu_free(u_row_dev)
    call check_dealloc_GPU_f("tridiag: (u_row_dev", 2086,  successGPU)

    successGPU = gpu_free(v_col_dev)
    call check_dealloc_GPU_f("tridiag: v_col_dev", 2089,  successGPU)

    successGPU = gpu_free(u_col_dev)
    call check_dealloc_GPU_f("tridiag: u_col_dev ", 2092,  successGPU)

    successGPU = gpu_free(vu_stored_rows_dev)
    call check_dealloc_GPU_f("tridiag: vu_stored_rows_dev ", 2095,  successGPU)

    successGPU = gpu_free(uv_stored_cols_dev)
    call check_dealloc_GPU_f("tridiag:uv_stored_cols_dev ", 2098,  successGPU)

    successGPU = gpu_free(d_vec_dev)
    call check_dealloc_GPU_f("tridiag: d_vec_dev", 2101,  successGPU)

    successGPU = gpu_free(e_vec_dev)
    call check_dealloc_GPU_f("tridiag: e_vec_dev", 2104,  successGPU)

    successGPU = gpu_free(tau_dev)
    call check_dealloc_GPU_f("tridiag: tau_dev", 2107,  successGPU)

    successGPU = gpu_free(aux_dev)
    call check_dealloc_GPU_f("tridiag: aux_dev", 2110,  successGPU)

    successGPU = gpu_free(aux1_dev)
    call check_dealloc_GPU_f("tridiag: aux1_dev", 2113,  successGPU)

    successGPU = gpu_free(aux_complex_dev)
    call check_dealloc_GPU_f("tridiag: aux_complex_dev", 2116,  successGPU)

    successGPU = gpu_free(vav_dev)
    call check_dealloc_GPU_f("tridiag: vav_dev", 2119,  successGPU)

    successGPU = gpu_free(dot_prod_dev)
    call check_dealloc_GPU_f("tridiag: dot_prod_dev", 2122,  successGPU)

    successGPU = gpu_free(xf_dev)
    call check_dealloc_GPU_f("tridiag: xf_dev", 2125,  successGPU)

  endif ! useGPU

  ! distribute the arrays d_vec and e_vec to all processors

  allocate(tmp_real(na), stat=istat, errmsg=errorMessage)
  call check_allocate_f("tridiag: tmp_real", 2138,  istat,  errorMessage)


  deallocate(tmp_real, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag: tmp_real", 2187,  istat,  errorMessage)

  if (useGPU) then

  else ! useGPU
    deallocate(v_row, v_col, u_row, u_col, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("tridiag: v_row, v_col, u_row, u_col", 2260,  istat,  errorMessage)
  endif ! useGPU

  deallocate(vu_stored_rows, uv_stored_cols, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag: vu_stored_rows, uv_stored_cols", 2264,  istat,  errorMessage)

  deallocate(aux, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag: aux", 2267,  istat,  errorMessage)



  call obj%timer%stop("tridiag_&
  &real&
  &" // &
  "_single" // &
  gpuString )

end subroutine tridiag_&
&real&
&_&
&single





!> \brief Transforms the eigenvectors of a tridiagonal matrix back
!>                     to the eigenvectors of the original matrix
!>                     (like Scalapack Routine PDORMTR)
!>
!  Parameters
!
!> \param na          Order of matrix a_mat, number of rows of matrix q_mat
!>
!> \param nqc         Number of columns of matrix q_mat
!>
!> \param a_mat(lda,matrixCols)  Matrix containing the Householder vectors (i.e. matrix a after tridiag_real)
!>                           Distribution is like in Scalapack.
!>
!> \param lda         Leading dimension of a_mat
!>
!> \param tau(na)     Factors of the Householder vectors
!>
!> \param q_mat           On input: Eigenvectors of tridiagonal matrix
!>                    On output: Transformed eigenvectors
!>                    Distribution is like in Scalapack.
!>
!> \param ldq         Leading dimension of q_mat
!>
!> \param nblk        blocksize of cyclic distribution, must be the same in both directions!
!>
!> \param matrixCols  local columns of matrix a_mat and q_mat
!>
!> \param mpi_comm_rows        MPI-Communicator for rows
!>
!> \param mpi_comm_cols        MPI-Communicator for columns
!>
!> \param useGPU      If true,  GPU version of the subroutine will be used
!>


subroutine trans_ev_&
&real&
&_&
&single &
(obj, na, nqc, a_mat, lda, tau, q_mat, ldq, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols, useGPU, success)
  use, intrinsic :: iso_c_binding
  use precision
  use elpa_abstract_impl
  use elpa_blas_interfaces
  use elpa_gpu
  use elpa_gpu_util

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

  class(elpa_abstract_impl_t), intent(inout)    :: obj
  integer(kind=ik), intent(in)                  :: na, nqc, lda, ldq, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols
  real(kind=rck), intent(in)           :: tau(na)

  real(kind=rck), intent(inout)        :: a_mat(lda,*)
  real(kind=rck), intent(inout)        :: q_mat(ldq,*)
  logical, intent(in)                           :: useGPU
  integer(kind=ik)                              :: max_stored_rows, max_stored_rows_fac

  integer(kind=ik)                              :: my_prow, my_pcol, np_rows, np_cols
  integer(kind=MPI_KIND)                        :: mpierr, my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI
  integer(kind=ik)                              :: totalblocks, max_blocks_row, max_blocks_col, max_local_rows, max_local_cols
  integer(kind=ik)                              :: l_cols, l_rows, l_colh, nstor
  integer(kind=ik)                              :: istep, n, nc, ic, ics, ice, nb, cur_pcol
  integer(kind=ik)                              :: hvn_ubnd, hvm_ubnd

  real(kind=rck), allocatable          :: hvb(:), hvm(:,:)
  real(kind=rck), pointer              :: tmp1(:), tmp2(:)
  real(kind=rck), allocatable          :: h1(:), h2(:), tmp_debug(:)
  real(kind=rck), pointer              :: tmat(:,:)
  real(kind=rck), pointer              :: hvm1(:)
  type(c_ptr)                                   :: tmp1_host, tmp2_host
  type(c_ptr)                                   :: hvm1_host, tmat_host

  integer(kind=ik)                              :: istat
  character(200)                                :: errorMessage
  character(20)                                 :: gpuString

  integer(kind=c_intptr_t)                      :: num
  integer(kind=C_intptr_T)                      :: q_dev, tmp_dev, hvm_dev, tmat_dev

  integer(kind=ik)                              :: blockStep
  logical                                       :: successGPU
  integer(kind=c_intptr_t), parameter           :: size_of_datatype = size_of_&
                                                                      &single&
                                                                      &_&
                                                                      &real
  integer(kind=ik)                              :: error
  integer(kind=MPI_KIND)                        :: bcast_request1, allreduce_request1, allreduce_request2
  logical                                       :: useNonBlockingCollectivesCols
  logical                                       :: useNonBlockingCollectivesRows
  integer(kind=c_int)                           :: non_blocking_collectives_rows, non_blocking_collectives_cols
  logical                                       :: success
  integer(kind=c_intptr_t)                      :: gpuHandle, my_stream


! NCCL requires streams
!#if defined (USE_CCL_TRANS_EV) && !defined(WITH_GPU_STREAMS)
!  successGPU = cuda_stream_create(obj%gpu_setup%my_stream)
!  if (.not.(successGPU)) then
!    print *,"Cannot create gpu stream handle"
!  endif
!  my_stream = obj%gpu_setup%my_stream
!#endif

  success = .true.

  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  call obj%timer%start("trans_ev_&
  &real&
  &" // &
  &"_single" //&
  gpuString)

  call obj%get("nbc_row_elpa1_tridi_to_full", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    write(error_unit,*) "Problem setting option for non blocking collectives for rows in elpa1_tridi_to_full. Aborting..."
    call obj%timer%stop("trans_ev_&
    &real&
    &" // &
    &"_single" //&
    gpuString)
    success = .false.
    return
  endif

  call obj%get("nbc_col_elpa1_tridi_to_full", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    write(error_unit,*) "Problem setting option for non blocking collectives for cols in elpa1_tridi_to_full. Aborting..."
    call obj%timer%stop("trans_ev_&
    &real&
    &" // &
    &"_single" //&
    gpuString)
    success = .false.
    return
  endif

  if (non_blocking_collectives_rows .eq. 1) then
    useNonBlockingCollectivesRows = .true.
  else
    useNonBlockingCollectivesRows = .false.
  endif

  if (non_blocking_collectives_cols .eq. 1) then
    useNonBlockingCollectivesCols = .true.
  else
    useNonBlockingCollectivesCols = .false.
  endif

  my_prow = obj%mpi_setup%myRank_comm_rows
  my_pcol = obj%mpi_setup%myRank_comm_cols

  np_rows = obj%mpi_setup%nRanks_comm_rows
  np_cols = obj%mpi_setup%nRanks_comm_cols


  !call obj%timer%start("mpi_communication")
  !call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND) ,my_prowMPI, mpierr)
  !call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND) ,np_rowsMPI, mpierr)
  !call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND) ,my_pcolMPI, mpierr)
  !call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND) ,np_colsMPI, mpierr)

  !my_prow = int(my_prowMPI, kind=c_int)
  !np_rows = int(np_rowsMPI, kind=c_int)
  !my_pcol = int(my_pcolMPI, kind=c_int)
  !np_cols = int(np_colsMPI, kind=c_int)
  !call obj%timer%stop("mpi_communication")

  call obj%get("max_stored_rows",max_stored_rows_fac, error)

  totalblocks = (na-1)/nblk + 1
  max_blocks_row = (totalblocks-1)/np_rows + 1
  max_blocks_col = ((nqc-1)/nblk)/np_cols + 1  ! Columns of q_mat!

  max_local_rows = max_blocks_row*nblk
  max_local_cols = max_blocks_col*nblk

  max_stored_rows = (max_stored_rows_fac/nblk+1)*nblk
 
  if (.not.(useGPU)) then
    allocate(tmat(max_stored_rows,max_stored_rows), stat=istat, errmsg=errorMessage)
    call check_alloc("trans_ev_&
    &real&
    &", "tmat", istat, errorMessage)

    allocate(tmp1(max_local_cols*max_stored_rows), stat=istat, errmsg=errorMessage)
    call check_alloc("trans_ev_&
    &real&
    &", "tmp1", istat, errorMessage)

    allocate(tmp2(max_local_cols*max_stored_rows), stat=istat, errmsg=errorMessage)
    call check_alloc("trans_ev_&
    &real&
    &", "tmp2", istat, errorMessage)
  endif

  allocate(h1(max_stored_rows*max_stored_rows), stat=istat, errmsg=errorMessage)
  call check_alloc("trans_ev_&
  &real&
  &", "h1", istat, errorMessage)

  allocate(h2(max_stored_rows*max_stored_rows), stat=istat, errmsg=errorMessage)
  call check_alloc("trans_ev_&
  &real&
  &", "h2", istat, errorMessage)

  allocate(hvb(max_local_rows*nblk), stat=istat, errmsg=errorMessage)
  call check_alloc("trans_ev_&
  &real&
  &", "hvn", istat, errorMessage)

  allocate(hvm(max_local_rows,max_stored_rows), stat=istat, errmsg=errorMessage)
  call check_alloc("trans_ev_&
  &real&
  &", "hvm", istat, errorMessage)

  hvm = 0   ! Must be set to 0 !!!
  hvb = 0   ! Safety only
  blockStep = nblk

  l_cols = local_index(nqc, my_pcol, np_cols, nblk, -1) ! Local columns of q_mat

  nstor = 0
  if (useGPU) then
    hvn_ubnd = 0
  endif

 
  if (useGPU) then
    ! todo: this is used only for copying hmv to device.. it should be possible to go without it
    !allocate(hvm1(max_local_rows*max_stored_rows), stat=istat, errmsg=errorMessage)
    !call check_alloc("trans_ev_&
    !&real&
    !&", "hvm1", istat, errorMessage)
    successGPU = gpu_malloc(tmat_dev, max_stored_rows * max_stored_rows * size_of_datatype)
    call check_alloc_GPU_f("trans_ev", 362,  successGPU)

    successGPU = gpu_malloc(hvm_dev, max_local_rows * max_stored_rows * size_of_datatype)
    call check_alloc_GPU_f("trans_ev", 365,  successGPU)

    successGPU = gpu_malloc(tmp_dev, max_local_cols * max_stored_rows * size_of_datatype)
    call check_alloc_GPU_f("trans_ev", 368,  successGPU)

    num = ldq * matrixCols * size_of_datatype
    successGPU = gpu_malloc(q_dev, num)
    call check_alloc_GPU_f("trans_ev", 372,  successGPU)
  

    successGPU = gpu_memcpy(q_dev, int(loc(q_mat(1,1)),kind=c_intptr_t), &
                  num, gpuMemcpyHostToDevice)
    call check_memcpy_GPU_f("trans_ev", 394,  successGPU)
  endif  ! useGPU

  do istep = 1, na, blockStep


    ics = MAX(istep,3)
    ice = MIN(istep+nblk-1,na)
    if (ice<ics) cycle

    cur_pcol = pcol(istep, nblk, np_cols)

    nb = 0
    do ic = ics, ice

      l_colh = local_index(ic  , my_pcol, np_cols, nblk, -1) ! Column of Householder Vector
      l_rows = local_index(ic-1, my_prow, np_rows, nblk, -1) ! # rows of Householder Vector


      if (my_pcol == cur_pcol) then
        hvb(nb+1:nb+l_rows) = a_mat(1:l_rows,l_colh)
        if (my_prow == prow(ic-1, nblk, np_rows)) then
          hvb(nb+l_rows) = 1.
        endif
      endif

      nb = nb+l_rows
    enddo


    nb = 0
    do ic = ics, ice
      l_rows = local_index(ic-1, my_prow, np_rows, nblk, -1) ! # rows of Householder Vector
      hvm(1:l_rows,nstor+1) = hvb(nb+1:nb+l_rows)
      if (useGPU) then
        hvm_ubnd = l_rows
      endif
      nstor = nstor+1
      nb = nb+l_rows
    enddo

    ! Please note: for smaller matix sizes (na/np_rows<=256), a value of 32 for nstor is enough!
    if (nstor+nblk > max_stored_rows .or. istep+nblk > na .or. (na/np_rows <= 256 .and. nstor >= 32)) then

      ! Calculate scalar products of stored vectors.
      ! This can be done in different ways, we use dsyrk or zherk

      tmat = 0
      call obj%timer%start("blas")
      if (l_rows>0) then
        call SSYRK('U', 'T',   &
                         int(nstor,kind=BLAS_KIND), int(l_rows,kind=BLAS_KIND), ONE, &
                         hvm, int(ubound(hvm,dim=1),kind=BLAS_KIND), ZERO, tmat, int(max_stored_rows,kind=BLAS_KIND))
      endif
      call obj%timer%stop("blas")
      nc = 0
      do n = 1, nstor-1
        h1(nc+1:nc+n) = tmat(1:n,n+1)
        nc = nc+n
      enddo

      if (nc > 0) h2 = h1

      ! Calculate triangular matrix T

      nc = 0
      tmat(1,1) = tau(ice-nstor+1)
      do n = 1, nstor-1
        call obj%timer%start("blas")
        call STRMV('L', 'T' , 'N', int(n,kind=BLAS_KIND), tmat, &
                            int(max_stored_rows,kind=BLAS_KIND), h2(nc+1), 1_BLAS_KIND)
        call obj%timer%stop("blas")

        tmat(n+1,1:n) = &
        -h2(nc+1:nc+n)  &
        *tau(ice-nstor+n+1)

        tmat(n+1,n+1) = tau(ice-nstor+n+1)
        nc = nc+n
      enddo

      if (useGPU) then
        ! todo: is this reshape really neccessary?
        hvm1(1:hvm_ubnd*nstor) = reshape(hvm(1:hvm_ubnd,1:nstor), (/ hvm_ubnd*nstor /))

        !hvm_dev(1:hvm_ubnd*nstor) = hvm1(1:hvm_ubnd*nstor)
        successGPU = gpu_memcpy(hvm_dev, int(loc(hvm1(1)),kind=c_intptr_t),   &
                      hvm_ubnd * nstor * size_of_datatype, gpuMemcpyHostToDevice)

        call check_memcpy_GPU_f("trans_ev", 546,  successGPU)

        !tmat_dev = tmat
        successGPU = gpu_memcpy(tmat_dev, int(loc(tmat(1,1)),kind=c_intptr_t),   &
                      max_stored_rows * max_stored_rows * size_of_datatype, gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("trans_ev", 551,  successGPU)
      endif


      ! Q = Q - V * T * V**T * Q

      if (l_rows > 0) then
        if (useGPU) then
          call obj%timer%start("gpublas")
          gpuHandle = obj%gpu_setup%gpublasHandleArray(0)
          call gpublas_SGEMM('T', 'N',   &
                                   nstor, l_cols, l_rows, ONE, hvm_dev, hvm_ubnd,  &
                                   q_dev, ldq, ZERO, tmp_dev, nstor, gpuHandle)
          call obj%timer%stop("gpublas")
        else ! useGPU

          call obj%timer%start("blas")
          call SGEMM('T', 'N',  &
                              int(nstor,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), &
                              int(l_rows,kind=BLAS_KIND), ONE, hvm, int(ubound(hvm,dim=1),kind=BLAS_KIND), &
                              q_mat, int(ldq,kind=BLAS_KIND), ZERO, tmp1, int(nstor,kind=BLAS_KIND))
          call obj%timer%stop("blas")
        endif ! useGPU

      else !l_rows>0

        if (useGPU) then
          if (gpu_vendor() /= OPENMP_OFFLOAD_GPU) then
            successGPU = gpu_memset(tmp_dev, 0, l_cols * nstor * size_of_datatype)
            call check_memcpy_GPU_f("trans_ev", 590,  successGPU)
          else
            allocate(tmp_debug(l_cols * nstor))
            tmp_debug(:) = 0.
            successGPU = gpu_memcpy(tmp_dev, int(loc(tmp_debug),kind=c_intptr_t), &
                                    l_cols*nstor*size_of_datatype, gpuMemcpyHostToDevice)
            call check_memcpy_GPU_f("trans_ev", 596,  successGPU)
            deallocate(tmp_debug)
          endif
        else
          tmp1(1:l_cols*nstor) = 0
        endif
      endif  !l_rows>0

!     tmp2 = tmp1

      if (l_rows > 0) then
        if (useGPU) then
          call obj%timer%start("gpublas")
          gpuHandle = obj%gpu_setup%gpublasHandleArray(0)
          call gpublas_STRMM('L', 'L', 'N', 'N',     &
                                   nstor, l_cols, ONE, tmat_dev, max_stored_rows,  &
                                   tmp_dev, nstor, gpuHandle)

          call gpublas_SGEMM('N', 'N' ,l_rows ,l_cols ,nstor,  &
                                   -ONE, hvm_dev, hvm_ubnd, tmp_dev, nstor,   &
                                   ONE, q_dev, ldq, gpuHandle)
          call obj%timer%stop("gpublas")
        else !useGPU
          call obj%timer%start("blas")

          call STRMM('L', 'L', 'N', 'N', int(nstor,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND),   &
                              ONE, tmat, int(max_stored_rows,kind=BLAS_KIND), tmp1, int(nstor,kind=BLAS_KIND))
          call SGEMM('N', 'N', int(l_rows,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), &
                              int(nstor,kind=BLAS_KIND), -ONE, hvm, int(ubound(hvm,dim=1),kind=BLAS_KIND), &
                              tmp1, int(nstor,kind=BLAS_KIND), ONE, q_mat, int(ldq,kind=BLAS_KIND))
          call obj%timer%stop("blas")
        endif ! useGPU
      endif  ! l_rows>0
      nstor = 0
    endif  ! (nstor+nblk>max_stored_rows .or. istep+nblk>na .or. (na/np_rows<=256 .and. nstor>=32))

    
  enddo ! istep = 1, na, blockStep

  deallocate(h1, h2, hvb, hvm, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_&     &MATH_DATATYPE&     &: h1, h2, hvb, hvm", 848,  istat,  errorMessage)

  if (useGPU) then

    !q_mat = q_dev
    successGPU = gpu_memcpy(int(loc(q_mat(1,1)),kind=c_intptr_t), &
                  q_dev, ldq * matrixCols * size_of_datatype, gpuMemcpyDeviceToHost)
    call check_memcpy_GPU_f("trans_ev", 864,  successGPU)

    !deallocate(hvm1, stat=istat, errmsg=errorMessage)
    !if (istat .ne. 0) then
    !  print *,"trans_ev_&
    !  &real&
    !  &: error when deallocating hvm1 "//errorMessage
    !  stop 1
    !endif

    !deallocate(q_dev, tmp_dev, hvm_dev, tmat_dev)
    successGPU = gpu_free(q_dev)
    call check_dealloc_GPU_f("trans_ev", 906,  successGPU)

    successGPU = gpu_free(tmp_dev)
    call check_dealloc_GPU_f("trans_ev", 909,  successGPU)

    successGPU = gpu_free(hvm_dev)
    call check_dealloc_GPU_f("trans_ev", 912,  successGPU)

    successGPU = gpu_free(tmat_dev)
    call check_dealloc_GPU_f("trans_ev", 915,  successGPU)
  else ! useGPU
    deallocate(tmat, tmp1, tmp2, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("trans_ev_&     &MATH_DATATYPE&     &: tmat, tmp1, tmp2", 920,  istat,  errorMessage)
  endif ! useGPU

!#if defined (USE_CCL_TRANS_EV) && !defined(WITH_GPU_STREAMS)
!  successGPU = cuda_stream_destroy(obj%gpu_setup%my_stream)
!  if (.not.(successGPU)) then
!    print *,"Cannot destroy gpu stream handle"
!  endif
!#endif

  call obj%timer%stop("trans_ev_&
  &real&
  &" // &
  &"_single" // &
  gpuString )

end subroutine trans_ev_&
&real&
&_&
&single


! now comes a dirty hack:
! the file elpa1_solve_tridi_real_template.F90 must be included twice
! for the legacy and for the new API. In the new API, however, some routines
! must be named "..._impl"

!#include "elpa1_solve_tridi_real_template.F90"




!cannot use "../src/elpa1/../solve_tridi/../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)


subroutine solve_tridi_&
&single_impl &
    ( obj, na, nev, d, e, q, ldq, nblk, matrixCols, mpi_comm_all, mpi_comm_rows, &
                                           mpi_comm_cols, useGPU, wantDebug, success, max_threads )

      use precision
      use elpa_abstract_impl
      use merge_recursive
      use merge_systems
      use elpa_mpi
      use ELPA_utilities
      use distribute_global_column
      use elpa_mpi
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
      integer(kind=ik), intent(in)               :: na, nev, ldq, nblk, matrixCols, &
                                                    mpi_comm_all, mpi_comm_rows, mpi_comm_cols
      real(kind=rk4), intent(inout)    :: d(na), e(na)
      real(kind=rk4), intent(inout)    :: q(ldq,*)
      logical, intent(in)                        :: useGPU, wantDebug
      logical, intent(out)                       :: success

      integer(kind=ik)                           :: i, j, n, np, nc, nev1, l_cols, l_rows
      integer(kind=ik)                           :: my_prow, my_pcol, np_rows, np_cols
      integer(kind=MPI_KIND)                     :: mpierr, my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI
      integer(kind=ik), allocatable              :: limits(:), l_col(:), p_col(:), l_col_bc(:), p_col_bc(:)

      integer(kind=ik)                           :: istat
      character(200)                             :: errorMessage
      character(20)                              :: gpuString
      integer(kind=ik), intent(in)               :: max_threads

      if(useGPU) then
        gpuString = "_gpu"
      else
        gpuString = ""
      endif

      call obj%timer%start("solve_tridi" // "_single" // gpuString)

      call obj%timer%start("mpi_communication")
      call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND) ,my_prowMPI, mpierr)
      call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND) ,np_rowsMPI, mpierr)
      call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND) ,my_pcolMPI, mpierr)
      call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND) ,np_colsMPI, mpierr)

      my_prow = int(my_prowMPI,kind=c_int)
      np_rows = int(np_rowsMPI,kind=c_int)
      my_pcol = int(my_pcolMPI,kind=c_int)
      np_cols = int(np_colsMPI,kind=c_int)

      call obj%timer%stop("mpi_communication")

      success = .true.

      l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a and q
      l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local columns of q

      ! Set Q to 0
      q(1:l_rows, 1:l_cols) = 0.0_rk

      ! Get the limits of the subdivisons, each subdivison has as many cols
      ! as fit on the respective processor column

      allocate(limits(0:np_cols), stat=istat, errmsg=errorMessage)
      call check_allocate_f("solve_tridi: limits", 128,  istat,  errorMessage)

      limits(0) = 0
      do np=0,np_cols-1
        nc = local_index(na, np, np_cols, nblk, -1) ! number of columns on proc column np

        ! Check for the case that a column has have zero width.
        ! This is not supported!
        ! Scalapack supports it but delivers no results for these columns,
        ! which is rather annoying
        if (nc==0) then
          call obj%timer%stop("solve_tridi" // "_single")
          if (wantDebug) write(error_unit,*) 'ELPA1_solve_tridi: ERROR: Problem contains processor column with zero width'
          success = .false.
          return
        endif
        limits(np+1) = limits(np) + nc
      enddo

      ! Subdivide matrix by subtracting rank 1 modifications

      do i=1,np_cols-1
        n = limits(i)
        d(n) = d(n)-abs(e(n))
        d(n+1) = d(n+1)-abs(e(n))
      enddo

      ! Solve sub problems on processsor columns

      nc = limits(my_pcol) ! column after which my problem starts

      if (np_cols>1) then
        nev1 = l_cols ! all eigenvectors are needed
      else
        nev1 = MIN(nev,l_cols)
      endif
      call solve_tridi_col_&
           &single_impl &
             (obj, l_cols, nev1, nc, d(nc+1), e(nc+1), q, ldq, nblk,  &
                        matrixCols, mpi_comm_rows, useGPU, wantDebug, success, max_threads)
      if (.not.(success)) then
        call obj%timer%stop("solve_tridi" // "_single" // gpuString)
        return
      endif
      ! If there is only 1 processor column, we are done

      if (np_cols==1) then
        deallocate(limits, stat=istat, errmsg=errorMessage)
        call check_deallocate_f("solve_tridi: limits", 176,  istat,  errorMessage)

        call obj%timer%stop("solve_tridi" // "_single" // gpuString)
        return
      endif

      ! Set index arrays for Q columns

      ! Dense distribution scheme:

      allocate(l_col(na), stat=istat, errmsg=errorMessage)
      call check_allocate_f("solve_tridi: l_col", 187,  istat,  errorMessage)

      allocate(p_col(na), stat=istat, errmsg=errorMessage)
      call check_allocate_f("solve_tridi: p_col", 190,  istat,  errorMessage)

      n = 0
      do np=0,np_cols-1
        nc = local_index(na, np, np_cols, nblk, -1)
        do i=1,nc
          n = n+1
          l_col(n) = i
          p_col(n) = np
        enddo
      enddo

      ! Block cyclic distribution scheme, only nev columns are set:

      allocate(l_col_bc(na), stat=istat, errmsg=errorMessage)
      call check_allocate_f("solve_tridi: l_col_bc", 205,  istat,  errorMessage)

      allocate(p_col_bc(na), stat=istat, errmsg=errorMessage)
      call check_allocate_f("solve_tridi: p_col_bc", 208,  istat,  errorMessage)

      p_col_bc(:) = -1
      l_col_bc(:) = -1

      do i = 0, na-1, nblk*np_cols
        do j = 0, np_cols-1
          do n = 1, nblk
            if (i+j*nblk+n <= MIN(nev,na)) then
              p_col_bc(i+j*nblk+n) = j
              l_col_bc(i+j*nblk+n) = i/np_cols + n
             endif
           enddo
         enddo
      enddo

      ! Recursively merge sub problems
      call merge_recursive_&
           &single &
           (obj, 0, np_cols, ldq, matrixCols, nblk, &
           l_col, p_col, l_col_bc, p_col_bc, limits, &
           np_cols, na, q, d, e, &
           mpi_comm_all, mpi_comm_rows, mpi_comm_cols,&
           useGPU, wantDebug, success, max_threads)

      if (.not.(success)) then
        call obj%timer%stop("solve_tridi" // "_single" // gpuString)
        return
      endif

      deallocate(limits,l_col,p_col,l_col_bc,p_col_bc, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("solve_tridi: limits, l_col, p_col, l_col_bc, p_col_bc", 239,  istat,  errorMessage)

      call obj%timer%stop("solve_tridi" // "_single" // gpuString)
      return

    end subroutine solve_tridi_&
        &single_impl

    subroutine solve_tridi_col_&
    &single_impl &
      ( obj, na, nev, nqoff, d, e, q, ldq, nblk, matrixCols, mpi_comm_rows, useGPU, wantDebug, success, max_threads )

   ! Solves the symmetric, tridiagonal eigenvalue problem on one processor column
   ! with the divide and conquer method.
   ! Works best if the number of processor rows is a power of 2!
      use precision
      use elpa_abstract_impl
      use elpa_mpi
      use merge_systems
      use ELPA_utilities
      use distribute_global_column
      implicit none
      class(elpa_abstract_impl_t), intent(inout) :: obj

      integer(kind=ik)              :: na, nev, nqoff, ldq, nblk, matrixCols, mpi_comm_rows
      real(kind=rk4)      :: d(na), e(na)
      real(kind=rk4)      :: q(ldq,*)

      integer(kind=ik), parameter   :: min_submatrix_size = 16 ! Minimum size of the submatrices to be used

      real(kind=rk4), allocatable    :: qmat1(:,:), qmat2(:,:)
      integer(kind=ik)              :: i, n, np
      integer(kind=ik)              :: ndiv, noff, nmid, nlen, max_size
      integer(kind=ik)              :: my_prow, np_rows
      integer(kind=MPI_KIND)        :: mpierr, my_prowMPI, np_rowsMPI

      integer(kind=ik), allocatable :: limits(:), l_col(:), p_col_i(:), p_col_o(:)
      logical, intent(in)           :: useGPU, wantDebug
      logical, intent(out)          :: success
      integer(kind=ik)              :: istat
      character(200)                :: errorMessage

      integer(kind=ik), intent(in)  :: max_threads

      integer(kind=MPI_KIND)        :: bcast_request1, bcast_request2
      logical                       :: useNonBlockingCollectivesRows
      integer(kind=c_int)           :: non_blocking_collectives, error

      success = .true.

      call obj%timer%start("solve_tridi_col" // "_single")

      call obj%get("nbc_row_solve_tridi", non_blocking_collectives, error)
      if (error .ne. ELPA_OK) then
        write(error_unit,*) "Problem setting option for non blocking collectives for rows in solve_tridi. Aborting..."
        success = .false.
        return
      endif

      if (non_blocking_collectives .eq. 1) then
        useNonBlockingCollectivesRows = .true.
      else
        useNonBlockingCollectivesRows = .false.
      endif

      call obj%timer%start("mpi_communication")
      call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND), my_prowMPI, mpierr)
      call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND), np_rowsMPI, mpierr)

      my_prow = int(my_prowMPI,kind=c_int)
      np_rows = int(np_rowsMPI,kind=c_int)
      call obj%timer%stop("mpi_communication")
      success = .true.
      ! Calculate the number of subdivisions needed.

      n = na
      ndiv = 1
      do while(2*ndiv<=np_rows .and. n>2*min_submatrix_size)
        n = ((n+3)/4)*2 ! the bigger one of the two halves, we want EVEN boundaries
        ndiv = ndiv*2
      enddo

      ! If there is only 1 processor row and not all eigenvectors are needed
      ! and the matrix size is big enough, then use 2 subdivisions
      ! so that merge_systems is called once and only the needed
      ! eigenvectors are calculated for the final problem.

      if (np_rows==1 .and. nev<na .and. na>2*min_submatrix_size) ndiv = 2

      allocate(limits(0:ndiv), stat=istat, errmsg=errorMessage)
      call check_deallocate_f("solve_tridi_col: limits", 445,  istat,  errorMessage)

      limits(0) = 0
      limits(ndiv) = na

      n = ndiv
      do while(n>1)
        n = n/2 ! n is always a power of 2
        do i=0,ndiv-1,2*n
          ! We want to have even boundaries (for cache line alignments)
          limits(i+n) = limits(i) + ((limits(i+2*n)-limits(i)+3)/4)*2
        enddo
      enddo

      ! Calculate the maximum size of a subproblem

      max_size = 0
      do i=1,ndiv
        max_size = MAX(max_size,limits(i)-limits(i-1))
      enddo

      ! Subdivide matrix by subtracting rank 1 modifications

      do i=1,ndiv-1
        n = limits(i)
        d(n) = d(n)-abs(e(n))
        d(n+1) = d(n+1)-abs(e(n))
      enddo

      if (np_rows==1)    then

        ! For 1 processor row there may be 1 or 2 subdivisions
        do n=0,ndiv-1
          noff = limits(n)        ! Start of subproblem
          nlen = limits(n+1)-noff ! Size of subproblem

          call solve_tridi_single_problem_&
          &single_impl &
                                  (obj, nlen,d(noff+1),e(noff+1), &
                                    q(nqoff+noff+1,noff+1),ubound(q,dim=1), wantDebug, success)

          if (.not.(success)) return
        enddo

      else

        ! Solve sub problems in parallel with solve_tridi_single
        ! There is at maximum 1 subproblem per processor

        allocate(qmat1(max_size,max_size), stat=istat, errmsg=errorMessage)
        call check_deallocate_f("solve_tridi_col: qmat1", 495,  istat,  errorMessage)

        allocate(qmat2(max_size,max_size), stat=istat, errmsg=errorMessage)
        call check_deallocate_f("solve_tridi_col: qmat2", 498,  istat,  errorMessage)

        qmat1 = 0 ! Make sure that all elements are defined

        if (my_prow < ndiv) then

          noff = limits(my_prow)        ! Start of subproblem
          nlen = limits(my_prow+1)-noff ! Size of subproblem
          call solve_tridi_single_problem_&
          &single_impl &
                                    (obj, nlen,d(noff+1),e(noff+1),qmat1, &
                                    ubound(qmat1,dim=1), wantDebug, success)

          if (.not.(success)) return
        endif

        ! Fill eigenvectors in qmat1 into global matrix q

        do np = 0, ndiv-1

          noff = limits(np)
          nlen = limits(np+1)-noff
!          qmat2 = qmat1 ! is this correct
          do i=1,nlen

            call distribute_global_column_&
            &single &
                     (obj, qmat1(1,i), q(1,noff+i), nqoff+noff, nlen, my_prow, np_rows, nblk)
          enddo

        enddo

        deallocate(qmat1, qmat2, stat=istat, errmsg=errorMessage)
        call check_deallocate_f("solve_tridi_col: qmat1, qmat2", 561,  istat,  errorMessage)

      endif

      ! Allocate and set index arrays l_col and p_col

      allocate(l_col(na), p_col_i(na),  p_col_o(na), stat=istat, errmsg=errorMessage)
      call check_deallocate_f("solve_tridi_col: l_col, p_col_i, p_col_o", 568,  istat,  errorMessage)

      do i=1,na
        l_col(i) = i
        p_col_i(i) = 0
        p_col_o(i) = 0
      enddo

      ! Merge subproblems

      n = 1
      do while(n<ndiv) ! if ndiv==1, the problem was solved by single call to solve_tridi_single

        do i=0,ndiv-1,2*n

          noff = limits(i)
          nmid = limits(i+n) - noff
          nlen = limits(i+2*n) - noff

          if (nlen == na) then
            ! Last merge, set p_col_o=-1 for unneeded (output) eigenvectors
            p_col_o(nev+1:na) = -1
          endif
          call merge_systems_&
          &single &
                              (obj, nlen, nmid, d(noff+1), e(noff+nmid), q, ldq, nqoff+noff, nblk, &
                               matrixCols, int(mpi_comm_rows,kind=ik), int(mpi_comm_self,kind=ik), &
                               l_col(noff+1), p_col_i(noff+1), &
                               l_col(noff+1), p_col_o(noff+1), 0, 1, useGPU, wantDebug, success, max_threads)
          if (.not.(success)) return

        enddo

        n = 2*n

      enddo

      deallocate(limits, l_col, p_col_i, p_col_o, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("solve_tridi_col: limits, l_col, p_col_i, p_col_o", 606,  istat,  errorMessage)

      call obj%timer%stop("solve_tridi_col" // "_single")

    end subroutine solve_tridi_col_&
    &single_impl

    subroutine solve_tridi_single_problem_&
    &single_impl &
    (obj, nlen, d, e, q, ldq, wantDebug, success)

   ! Solves the symmetric, tridiagonal eigenvalue problem on a single processor.
   ! Takes precautions if DSTEDC fails or if the eigenvalues are not ordered correctly.
     use precision
     use elpa_abstract_impl
     use elpa_blas_interfaces
     use ELPA_utilities
     implicit none
     class(elpa_abstract_impl_t), intent(inout) :: obj
     integer(kind=ik)                         :: nlen, ldq
     real(kind=rk4)                 :: d(nlen), e(nlen), q(ldq,nlen)

     real(kind=rk4), allocatable    :: work(:), qtmp(:), ds(:), es(:)
     real(kind=rk4)                 :: dtmp

     integer(kind=ik)              :: i, j, lwork, liwork, info
     integer(kind=BLAS_KIND)       :: infoBLAS
     integer(kind=ik), allocatable :: iwork(:)

     logical, intent(in)           :: wantDebug
     logical, intent(out)          :: success
      integer(kind=ik)             :: istat
      character(200)               :: errorMessage

     call obj%timer%start("solve_tridi_single" // "_single")

     success = .true.
     allocate(ds(nlen), es(nlen), stat=istat, errmsg=errorMessage)
     call check_allocate_f("solve_tridi_single: ds, es", 644,  istat,  errorMessage)

     ! Save d and e for the case that dstedc fails

     ds(:) = d(:)
     es(:) = e(:)

     ! First try dstedc, this is normally faster but it may fail sometimes (why???)

     lwork = 1 + 4*nlen + nlen**2
     liwork =  3 + 5*nlen
     allocate(work(lwork), iwork(liwork), stat=istat, errmsg=errorMessage)
     call check_allocate_f("solve_tridi_single: work, iwork", 656,  istat,  errorMessage)
     call obj%timer%start("lapack")
     call SSTEDC('I', int(nlen,kind=BLAS_KIND), d, e, q, int(ldq,kind=BLAS_KIND),    &
                          work, int(lwork,kind=BLAS_KIND), int(iwork,kind=BLAS_KIND), int(liwork,kind=BLAS_KIND), &
                          infoBLAS)
     info = int(infoBLAS,kind=ik)
     call obj%timer%stop("lapack")

     if (info /= 0) then

       ! DSTEDC failed, try DSTEQR. The workspace is enough for DSTEQR.

       write(error_unit,'(a,i8,a)') 'Warning: Lapack routine DSTEDC failed, info= ',info,', Trying DSTEQR!'

       d(:) = ds(:)
       e(:) = es(:)
       call obj%timer%start("lapack")
       call SSTEQR('I', int(nlen,kind=BLAS_KIND), d, e, q, int(ldq,kind=BLAS_KIND), work, infoBLAS )
       info = int(infoBLAS,kind=ik)
       call obj%timer%stop("lapack")

       ! If DSTEQR fails also, we don't know what to do further ...

       if (info /= 0) then
         if (wantDebug) then
           write(error_unit,'(a,i8,a)') 'ELPA1_solve_tridi_single: ERROR: Lapack routine DSTEQR failed, info= ',info,', Aborting!'
         endif
         success = .false.
         return
       endif
     end if

       deallocate(work,iwork,ds,es, stat=istat, errmsg=errorMessage)
       call check_deallocate_f("solve_tridi_single: work, iwork, ds, es", 689,  istat,  errorMessage)

      ! Check if eigenvalues are monotonically increasing
      ! This seems to be not always the case  (in the IBM implementation of dstedc ???)

      do i=1,nlen-1
        if (d(i+1)<d(i)) then
          if (abs(d(i+1) - d(i)) / abs(d(i+1) + d(i)) > 1e-14_rk4) then
            write(error_unit,'(a,i8,2g25.16)') '***WARNING: Monotony error dste**:',i+1,d(i),d(i+1)
          else
            write(error_unit,'(a,i8,2g25.16)') 'Info: Monotony error dste{dc,qr}:',i+1,d(i),d(i+1)
            write(error_unit,'(a)') 'The eigenvalues from a lapack call are not sorted to machine precision.'
            write(error_unit,'(a)') 'In this extent, this is completely harmless.'
            write(error_unit,'(a)') 'Still, we keep this info message just in case.'
          end if
          allocate(qtmp(nlen), stat=istat, errmsg=errorMessage)
          call check_allocate_f("solve_tridi_single: qtmp", 709,  istat,  errorMessage)

          dtmp = d(i+1)
          qtmp(1:nlen) = q(1:nlen,i+1)
          do j=i,1,-1
            if (dtmp<d(j)) then
              d(j+1)        = d(j)
              q(1:nlen,j+1) = q(1:nlen,j)
            else
              exit ! Loop
            endif
          enddo
          d(j+1)        = dtmp
          q(1:nlen,j+1) = qtmp(1:nlen)
          deallocate(qtmp, stat=istat, errmsg=errorMessage)
          call check_deallocate_f("solve_tridi_single: qtmp", 724,  istat,  errorMessage)

       endif
     enddo
     call obj%timer%stop("solve_tridi_single" // "_single")

    end subroutine solve_tridi_single_problem_&
    &single_impl

!#include "elpa1_merge_systems_real_template.F90"






subroutine hh_transform_real_&
&single &
(obj, alpha, xnorm_sq, xf, tau, wantDebug)
  ! Similar to LAPACK routine DLARFP, but uses ||x||**2 instead of x(:)
  ! and returns the factor xf by which x has to be scaled.
  ! It also hasn't the special handling for numbers < 1.d-300 or > 1.d150
  ! since this would be expensive for the parallel implementation.
  use precision
  use elpa_abstract_impl
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

  class(elpa_abstract_impl_t), intent(inout)    :: obj
  logical, intent(in)                           :: wantDebug
  real(kind=rk), intent(inout)       :: alpha
  real(kind=rk), intent(in)          :: xnorm_sq
  real(kind=rk), intent(out)         :: xf, tau

  real(kind=rk)                      :: BETA

  if (wantDebug) call obj%timer%start("hh_transform_&
                   &real&
     	      &" // &
                   &"_single" )


  if ( XNORM_SQ==0.0_rk ) then

    if ( ALPHA>=0.0_rk ) then
      TAU = 0.0_rk
    else
      TAU = 2.0_rk
      ALPHA = -ALPHA
    endif
    XF = 0.0_rk

  else

    BETA = SIGN( SQRT( ALPHA**2 + XNORM_SQ ), ALPHA )
    ALPHA = ALPHA + BETA
    IF ( BETA<0 ) THEN
      BETA = -BETA
      TAU  = -ALPHA / BETA
    ELSE
      ALPHA = XNORM_SQ / ALPHA

      TAU = ALPHA / BETA
      ALPHA = -ALPHA
    END IF
    XF = 1.0_rk/ALPHA
    ALPHA = BETA
  endif

  if (wantDebug) call obj%timer%stop("hh_transform_&
  &real&
  &" // &
  &"_single" )

end subroutine hh_transform_real_&
    &single



! complex double precision

































!> \brief Reduces a distributed symmetric matrix to tridiagonal form (like Scalapack Routine PDSYTRD)
!>
!  Parameters
!
!> \param obj	      object of elpa_type
!> \param na          Order of matrix
!>
!> \param a_mat(matrixRows,matrixCols)    Distributed matrix which should be reduced.
!>              Distribution is like in Scalapack.
!>              Opposed to PDSYTRD, a(:,:) must be set completely (upper and lower half)
!>              a(:,:) is overwritten on exit with the Householder vectors
!>
!> \param matrixRows         Leading dimension of a
!>
!> \param nblk        blocksize of cyclic distribution, must be the same in both directions!
!>
!> \param matrixCols  local columns of matrix
!>
!> \param mpi_comm_rows        MPI-Communicator for rows
!> \param mpi_comm_cols        MPI-Communicator for columns
!>
!> \param d_vec(na)       Diagonal elements (returned), identical on all processors
!>
!> \param e_vec(na)       Off-Diagonal elements (returned), identical on all processors
!>
!> \param tau(na)     Factors for the Householder vectors (returned), needed for back transformation
!>
!> \param useGPU      If true,  GPU version of the subroutine will be used
!> \param wantDebug   if true more debug information
!>
subroutine tridiag_&
  &complex&
  &_&
  &double &
  (obj, na, a_mat, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols, d_vec, e_vec, tau, useGPU, wantDebug, &
   max_threads_in, isSkewsymmetric, success)
  use, intrinsic :: iso_c_binding
  use precision
  use elpa_abstract_impl
  use matrix_plot
  use elpa_omp
  use elpa_blas_interfaces
  use elpa_gpu
  use elpa_gpu_util
  use tridiag_gpu

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
  class(elpa_abstract_impl_t), intent(inout)    :: obj
  integer(kind=ik), intent(in)                  :: na, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols
  logical, intent(in)                           :: useGPU, wantDebug
  logical, intent(in)                           :: isSkewsymmetric

  logical                                       :: useCCL=.false.

  complex(kind=rck), intent(out)          :: tau(na)
  complex(kind=rck), intent(inout)        :: a_mat(matrixRows,*)
  real(kind=rk), intent(out)                    :: d_vec(na)
  real(kind=rk), intent(out)                    :: e_vec(na)
  integer(kind=ik)                              :: max_stored_uv = 32 ! TODO_23_11 - make it tunable instead of hard-coded
  logical,          parameter                   :: mat_vec_as_one_block = .true.

  ! id in processor row and column and total numbers of processor rows and columns
  integer(kind=ik)                              :: my_prow, my_pcol, np_rows, np_cols
  integer(kind=MPI_KIND)                        :: my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI
  integer(kind=MPI_KIND)                        :: mpierr
  integer(kind=ik)                              :: totalblocks, max_loc_block_rows, max_loc_block_cols, max_local_rows, &
                                                   max_local_cols
  ! updated after each istep (in the main cycle) to contain number of
  ! local columns and rows of the remaining part of the matrix
  integer(kind=ik)                              :: l_cols, l_rows
  integer(kind=ik)                              :: n_stored_vecs
  integer(kind=ik)                              :: isOurProcessRowInt ! TODO_23_11 - get rid of it
  logical                                       :: isOurProcessRow, isOurProcessCol, isOurProcessCol_prev


  integer(kind=C_intptr_T)                      :: a_dev, v_row_dev, v_col_dev, u_row_dev, u_col_dev, vu_stored_rows_dev, &
                                                   uv_stored_cols_dev, d_vec_dev, e_vec_dev, tau_dev
  logical                                       :: successGPU

  integer(kind=ik)                              :: istep, i, j, l_col_beg, l_col_end, l_row_beg, l_row_end
  integer(kind=ik)                              :: tile_size, l_rows_per_tile, l_cols_per_tile
  integer(kind=c_intptr_t)                      :: offset_dev

  integer(kind=ik), intent(in)                  :: max_threads_in
  integer(kind=ik)                              :: max_threads

  real(kind=rk)                                 :: vnorm2
  complex(kind=rck)                       :: vav, x, aux1(2), vrl, xf, conjg_tau, dot_prod
  complex(kind=rck), allocatable          :: aux(:)

  integer(kind=c_intptr_t)                      :: aux_dev, aux1_dev, aux_complex_dev, vav_dev, vav_host_or_dev, dot_prod_dev, & 
                                                   xf_dev, a_updated_element_dev, xf_host_or_dev, tau_istep_host_or_dev
  integer(kind=c_intptr_t)                      :: vnorm2_dev, vrl_dev ! alias pointers for aux1_dev(1), aux1_dev(2)

  complex(kind=rck)                             :: aux3(1)

  integer(kind=c_intptr_t)                      :: num
  complex(kind=rck), pointer              :: v_row(:), & ! used to store calculated Householder Vector
                                                   v_col(:)   ! the same Vector, but transposed 
  complex(kind=rck), pointer              :: u_row_debug(:), & ! used to store calculated Householder Vector
                                                   u_col_debug(:)   ! the same Vector, but transposed 
  complex(kind=rck), pointer              :: u_col(:), u_row(:)

  ! the following two matrices store pairs of vectors v and u calculated in each step
  ! at most max_stored_uv Vector pairs are stored, than the matrix A_i is explicitli updated
  ! u and v are stored both in row and Vector forms
  ! pattern: v1,u1,v2,u2,v3,u3,....
  ! todo: It is little bit confusing, I think, that variables _row actually store columns and vice versa
  !complex(kind=rck), pointer             :: vu_stored_rows(:,:)
  complex(kind=rck), allocatable         :: vu_stored_rows(:,:)
  ! pattern: u1,v1,u2,v2,u3,v3,....
  complex(kind=rck), allocatable         :: uv_stored_cols(:,:)


  type(c_ptr)                                   :: v_row_host, v_col_host
  type(c_ptr)                                   :: u_row_host, u_col_host
  !type(c_ptr)                                   :: vu_stored_rows_host, uv_stored_cols_host
  real(kind=rk), allocatable                    :: tmp_real(:)
  integer(kind=ik)                              :: min_tile_size, error
  integer(kind=ik)                              :: istat
  character(200)                                :: errorMessage
  character(20)                                 :: gpuString
  integer(kind=ik)                              :: nblockEnd
  integer(kind=c_intptr_t), parameter           :: size_of_datatype = size_of_&
                                                                      &double&
                                                                      &_&
                                                                      &complex
  integer(kind=c_intptr_t), parameter           :: size_of_datatype_real = size_of_&
                                                                      &double&
                                                                      &_real
  integer(kind=MPI_KIND)                        :: bcast_request1, bcast_request2, bcast_request3
  integer(kind=MPI_KIND)                        :: allreduce_request1, allreduce_request2, allreduce_request3
  integer(kind=MPI_KIND)                        :: allreduce_request4, allreduce_request5, allreduce_request6, &
                                                   allreduce_request7
  logical                                       :: useNonBlockingCollectivesCols
  logical                                       :: useNonBlockingCollectivesRows
  integer(kind=c_int)                           :: non_blocking_collectives_rows, non_blocking_collectives_cols
  logical                                       :: success

  integer(kind=c_intptr_t)                      :: gpuHandle, my_stream
  integer(kind=c_int) :: pointerMode

  integer(kind=ik)                              :: string_length, sm_count


  allocate(aux(2*max_stored_uv), stat=istat, errmsg=errorMessage)

  success = .true.

  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  if (useGPU) then
    max_threads=1
  else
    max_threads=max_threads_in
  endif

  call obj%timer%start("tridiag_&
  &complex&
  &" // &
  "_double" // &
  gpuString )

  call obj%get("nbc_row_elpa1_full_to_tridi", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    write(error_unit,*) "Problem setting option for non blocking collectives for rows in elpa1_tridiag. Aborting..."
    success = .false.
    call obj%timer%stop("tridiag_&
    &complex&
    &" // &
    "_double" // &
    gpuString )
    return
  endif

  call obj%get("nbc_col_elpa1_full_to_tridi", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    write(error_unit,*) "Problem setting option for non blocking collectives for cols in elpa1_tridiag. Aborting..."
    success = .false.
    call obj%timer%stop("tridiag_&
    &complex&
    &" // &
    "_double" // &
    gpuString )
    return
  endif

  if (non_blocking_collectives_rows .eq. 1) then
    useNonBlockingCollectivesRows = .true.
  else
    useNonBlockingCollectivesRows = .false.
  endif

  if (non_blocking_collectives_cols .eq. 1) then
    useNonBlockingCollectivesCols = .true.
  else
    useNonBlockingCollectivesCols = .false.
  endif


  !mpi_comm_all    = obj%mpi_setup%mpi_comm_parent
  !mpi_comm_cols   = obj%mpi_setup%mpi_comm_cols
  !mpi_comm_rows   = obj%mpi_setup%mpi_comm_rows

  my_prow = obj%mpi_setup%myRank_comm_rows
  my_pcol = obj%mpi_setup%myRank_comm_cols

  np_rows = obj%mpi_setup%nRanks_comm_rows
  np_cols = obj%mpi_setup%nRanks_comm_cols


  !if (wantDebug) call obj%timer%start("mpi_communication")
  !call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND), my_prowMPI, mpierr)
  !call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND), np_rowsMPI, mpierr)
  !call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND), my_pcolMPI, mpierr)
  !call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND), np_colsMPI, mpierr)


  !my_prow = int(my_prowMPI, kind=c_int)
  !np_rows = int(np_rowsMPI, kind=c_int)
  !my_pcol = int(my_pcolMPI, kind=c_int)
  !np_cols = int(np_colsMPI, kind=c_int)
  !if (wantDebug) call obj%timer%stop("mpi_communication")

  ! Matrix is split into tiles; work is done only for tiles on the diagonal or above
  ! seems that tile is a square submatrix, consisting by several blocks
  ! it is a smallest possible square submatrix, where blocks being distributed among
  ! processors are "aligned" in both rows and columns
  !  -----------------
  ! | 1 4 | 1 4 | 1 4 | ...
  ! | 2 5 | 2 5 | 2 5 | ...
  ! | 3 6 | 3 6 | 3 6 | ...
  !  ----------------- ...
  ! | 1 4 | 1 4 | 1 4 | ...
  ! | 2 5 | 2 5 | 2 5 | ...
  ! | 3 6 | 3 6 | 3 6 | ...
  !  ----------------- .
  !   : :   : :   : :    .
  !   : :   : :   : :      .
  !
  ! this is a tile, where each number represents block, assigned to a processor with the shown number
  ! size of this small block is nblk
  ! Image is for situation with 6 processors, 3 processor rows and 2 columns
  ! tile_size is thus nblk * 6
  !
  tile_size = nblk*least_common_multiple(np_rows,np_cols) ! minimum global tile size

  ! make tile_size a smallest possible multiple of previously defined tile size, such that it is
  ! larger or equal to min_tile_size
  ! min_tile_size has been originally hardcoded as 128 * max(np_rows, np_cols), so it is now the implicit value
  ! it can, however, be set by the user
  call obj%get("min_tile_size", min_tile_size ,error)
  if (error .ne. ELPA_OK) then
    write(error_unit,*) "Problem setting option for min_tile_size. Aborting..."
    success = .false.
    call obj%timer%stop("tridiag_&
    &complex&
    &" // &
    "_double" // &
    gpuString )
    return
  endif
  if(min_tile_size == 0) then
    ! not set by the user, use the default value
    min_tile_size = 128*max(np_rows, np_cols)
  endif
  tile_size = ((min_tile_size-1)/tile_size+1)*tile_size

  nblockEnd = 3

  l_rows_per_tile = tile_size/np_rows ! local rows of a tile
  l_cols_per_tile = tile_size/np_cols ! local cols of a tile

  totalblocks = (na-1)/nblk + 1
  max_loc_block_rows = (totalblocks-1)/np_rows + 1
  max_loc_block_cols = (totalblocks-1)/np_cols + 1

  ! localy owned submatrix has size at most max_local_rows x max_local_cols at each processor
  max_local_rows = max_loc_block_rows*nblk
  max_local_cols = max_loc_block_cols*nblk

  ! allocate memmory for vectors
  ! todo: It is little bit confusing, I think, that variables _row actually store columns and vice versa
  ! todo: if something has length max_local_rows, it is actually a column, no?
  ! todo: probably one should read it as v_row = Vector v distributed among rows

  allocate(uv_stored_cols(max_local_cols,2*max_stored_uv), stat=istat, errmsg=errorMessage)
  call check_alloc("tridiag_&
       &complex ", "uv_stored_cols", istat, errorMessage)

  allocate(vu_stored_rows(max_local_rows,2*max_stored_uv), stat=istat, errmsg=errorMessage)
  call check_alloc("tridiag_&
       &complex ", "vu_stored_rows", istat, errorMessage)

  if (useGPU) then

    ! allocate v_row 1 element longer to allow store and broadcast tau together with it
    if (gpu_vendor() /= OPENMP_OFFLOAD_GPU .and. gpu_vendor() /= SYCL_GPU) then
      num = (max_local_rows+1) * size_of_datatype
      successGPU = gpu_malloc_host(v_row_host, num)
      call check_host_alloc_GPU_f("tridiag: v_row_host", 444,  successGPU)
      call c_f_pointer(v_row_host,v_row,(/(max_local_rows+1)/))
    else
      allocate(v_row(max_local_rows+1))
    endif

    if (gpu_vendor() /= OPENMP_OFFLOAD_GPU .and. gpu_vendor() /= SYCL_GPU) then
      num = (max_local_cols) * size_of_datatype
      successGPU = gpu_malloc_host(v_col_host,num)
      call check_host_alloc_GPU_f("tridiag: v_col_host", 453,  successGPU)
      call c_f_pointer(v_col_host,v_col,(/(max_local_cols)/))
    else
      allocate(v_col(max_local_cols))
    endif

    if (gpu_vendor() /= OPENMP_OFFLOAD_GPU .and. gpu_vendor() /= SYCL_GPU) then
      num = (max_local_cols) * size_of_datatype
      successGPU = gpu_malloc_host(u_col_host,num)
      call check_host_alloc_GPU_f("tridiag: u_col_host", 462,  successGPU)
      call c_f_pointer(u_col_host,u_col,(/(max_local_cols)/))
    else
      allocate(u_col(max_local_cols))
    endif

    if (gpu_vendor() /= OPENMP_OFFLOAD_GPU .and. gpu_vendor() /= SYCL_GPU) then
      num = (max_local_rows) * size_of_datatype
      successGPU = gpu_malloc_host(u_row_host,num)
      call check_host_alloc_GPU_f("tridiag: u_row_host", 471,  successGPU)
      call c_f_pointer(u_row_host,u_row,(/(max_local_rows)/))
    else
      allocate(u_row(max_local_rows))
    endif

    
    if (gpu_vendor() /= OPENMP_OFFLOAD_GPU .and. gpu_vendor() /= SYCL_GPU) then
      num = (max_local_rows * 2*max_stored_uv) * size_of_datatype
      successGPU = gpu_host_register(int(loc(vu_stored_rows),kind=c_intptr_t), num, gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: vu_stored_rows", 481,  successGPU)

      num = (max_local_cols * 2*max_stored_uv) * size_of_datatype
      successGPU = gpu_host_register(int(loc(uv_stored_cols),kind=c_intptr_t), num, gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: uv_stored_cols", 485,  successGPU)

      num = (1 * 2*max_stored_uv) * size_of_datatype
      successGPU = gpu_host_register(int(loc(aux),kind=c_intptr_t), num, gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: aux", 489,  successGPU)

      num = na * size_of_datatype_real
      successGPU = gpu_host_register(int(loc(d_vec),kind=c_intptr_t), num, gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: d_vec", 493,  successGPU)

      num = na * size_of_datatype_real
      successGPU = gpu_host_register(int(loc(e_vec),kind=c_intptr_t), num, gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: e_vec", 497,  successGPU)

      num = na * size_of_datatype
      successGPU = gpu_host_register(int(loc(tau),kind=c_intptr_t), num, gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: tau", 501,  successGPU)

      num = 2 * size_of_datatype
      successGPU = gpu_host_register(int(loc(aux1),kind=c_intptr_t), num, gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: aux1", 505,  successGPU)

      num = 1 * size_of_datatype
      successGPU = gpu_host_register(int(loc(vav),kind=c_intptr_t), num, gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: vav", 509,  successGPU)

      num = 1 * size_of_datatype
      successGPU = gpu_host_register(int(loc(dot_prod),kind=c_intptr_t), num, gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: dot_prod", 513,  successGPU)

      num = 1 * size_of_datatype
      successGPU = gpu_host_register(int(loc(xf),kind=c_intptr_t), num, gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: xf", 517,  successGPU)
    endif ! gpu_vendor() /= OPENMP_OFFLOAD_GPU .and. gpu_vendor() /= SYCL_GPU
  else ! useGPU

    allocate(v_row(max_local_rows+1), stat=istat, errmsg=errorMessage)
    call check_alloc("tridiag_&
    &complex ", "v_row", istat, errorMessage)

    allocate(v_col(max_local_cols), stat=istat, errmsg=errorMessage)
    call check_alloc("tridiag_&
     &complex ", "v_col", istat, errorMessage)

    allocate(u_col(max_local_cols), stat=istat, errmsg=errorMessage)
    call check_alloc("tridiag_&
    &complex ", "u_col", istat, errorMessage)

    allocate(u_row(max_local_rows), stat=istat, errmsg=errorMessage)
    call check_alloc("tridiag_&
    &complex ", "u_row", istat, errorMessage)
      
  endif ! useGPU


  v_row = 0
  u_row = 0
  v_col = 0
  u_col = 0

  if (useGPU) then
    successGPU = gpu_malloc(v_row_dev, (max_local_rows+1) * size_of_datatype)
    call check_alloc_GPU_f("tridiag: v_row_dev", 556,  successGPU)

    successGPU = gpu_malloc(u_row_dev, max_local_rows * size_of_datatype)

    call check_alloc_GPU_f("tridiag: u_row_dev", 560,  successGPU)

    successGPU = gpu_malloc(v_col_dev, max_local_cols * size_of_datatype)
    call check_alloc_GPU_f("tridiag: v_col_dev", 563,  successGPU)

    successGPU = gpu_malloc(u_col_dev, max_local_cols * size_of_datatype)
    call check_alloc_GPU_f("tridiag: u_col_dev", 566,  successGPU)

    successGPU = gpu_malloc(vu_stored_rows_dev, max_local_rows * 2 * max_stored_uv * size_of_datatype)
    call check_alloc_GPU_f("tridiag: vu_stored_rows_dev", 569,  successGPU)

    successGPU = gpu_malloc(uv_stored_cols_dev, max_local_cols * 2 * max_stored_uv * size_of_datatype)
    call check_alloc_GPU_f("tridiag: uv_stored_cols_dev", 572,  successGPU)

    successGPU = gpu_malloc(d_vec_dev, na * size_of_datatype_real)
    call check_alloc_GPU_f("tridiag: d_vec_dev", 575,  successGPU)

    successGPU = gpu_malloc(e_vec_dev, na * size_of_datatype_real)
    call check_alloc_GPU_f("tridiag: e_vec_dev", 578,  successGPU)

    successGPU = gpu_malloc(tau_dev, na * size_of_datatype)
    call check_alloc_GPU_f("tridiag: tau_dev", 581,  successGPU)

    successGPU = gpu_malloc(aux_dev, 2*max_stored_uv * size_of_datatype)
    call check_alloc_GPU_f("tridiag: aux_dev", 584,  successGPU)

    successGPU = gpu_malloc(aux1_dev, 2 * size_of_datatype)
    call check_alloc_GPU_f("tridiag: aux1_dev", 587,  successGPU)

    successGPU = gpu_malloc(aux_complex_dev, 2 *max_stored_uv* size_of_datatype)
    call check_alloc_GPU_f("tridiag: aux_complex_dev", 590,  successGPU)

    successGPU = gpu_malloc(vav_dev, 1 * size_of_datatype)
    call check_alloc_GPU_f("tridiag: vav_dev", 593,  successGPU)

    successGPU = gpu_malloc(dot_prod_dev, 1 * size_of_datatype)
    call check_alloc_GPU_f("tridiag: dot_prod_dev", 596,  successGPU)

    successGPU = gpu_malloc(xf_dev, 1 * size_of_datatype)
    call check_alloc_GPU_f("tridiag: xf_dev", 599,  successGPU)

  endif !useGPU

  d_vec(:) = 0
  e_vec(:) = 0
  tau(:) = 0

  if (useGPU) then
    successGPU = gpu_memset(d_vec_dev, 0, na * size_of_datatype_real)
    call check_memcpy_GPU_f("tridiag: d_vec_dev", 630,  successGPU)

    successGPU = gpu_memset(e_vec_dev, 0, na * size_of_datatype_real)
    call check_memcpy_GPU_f("tridiag: e_vec_dev", 633,  successGPU)

    successGPU = gpu_memset(tau_dev, 0, na * size_of_datatype)
    call check_memcpy_GPU_f("tridiag: tau_dev", 636,  successGPU)
  endif

  n_stored_vecs = 0

  l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a_mat
  l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a_mat

  if (my_prow == prow(na, nblk, np_rows) .and. my_pcol == pcol(na, nblk, np_cols)) then
      d_vec(na) = real(a_mat(l_rows,l_cols), kind=rk)
  endif

  if (useGPU) then
    ! allocate memory for matrix A on the device and than copy the matrix

    num = matrixRows * matrixCols * size_of_datatype

    successGPU = gpu_malloc(a_dev, num)
    call check_alloc_GPU_f("tridiag: a_dev", 659,  successGPU)


    successGPU = gpu_memcpy(a_dev, int(loc(a_mat(1,1)),kind=c_intptr_t), &
                              num, gpuMemcpyHostToDevice)
    call check_memcpy_GPU_f("tridiag: a_dev", 678,  successGPU)

  endif ! useGPU

  ! main cycle of tridiagonalization
  ! in each step, 1 Householder Vector is calculated
  do istep = na, nblockEnd ,-1


    ! Calculate number of local rows and columns of the still remaining matrix
    ! on the local processor
    l_rows = local_index(istep-1, my_prow, np_rows, nblk, -1)
    l_cols = local_index(istep-1, my_pcol, np_cols, nblk, -1)

    ! Calculate Vector for Householder transformation on all procs
    ! owning column istep

    if (useGPU) then 
      ! copy l_cols + 1 column of a_dev to v_row_dev
      isOurProcessRow      = (my_prow == prow(istep-1, nblk, np_rows))
      isOurProcessCol      = (my_pcol == pcol(istep-1, nblk, np_cols))
      isOurProcessCol_prev = (my_pcol == pcol(istep  , nblk, np_cols)) ! isOurProcessCol from the previous step
      call gpu_copy_and_set_zeros_double_complex(v_row_dev, a_dev, l_rows, l_cols, matrixRows, istep, &
                                            aux1_dev, vav_dev, d_vec_dev, &
                                            isOurProcessRow, isOurProcessCol, isOurProcessCol_prev, &
                                            isSkewsymmetric, useCCL, wantDebug, my_stream)
      if (gpu_vendor() /= SYCL_GPU) successGPU = gpu_DeviceSynchronize()
    endif ! useGPU

    if (my_pcol == pcol(istep, nblk, np_cols)) then

      ! Get Vector to be transformed; distribute last element and norm of
      ! remaining elements to all procs in current column

!             ! copy l_cols + 1 column of A to v_row
      ! if (useGPU) then

! #ifdef WITH_NVTX
!         call nvtxRangePush("memcpy new D-D a_dev(:,l_cols+1)->v_row_dev")
! #endif
!         ! TODO_23_11:  create a dev-dev copy kernel or merge it to another kernel
!         offset_dev = l_cols * matrixRows * size_of_datatype
!         successGPU = gpu_memcpy(v_row_dev, a_dev + offset_dev, (l_rows)* size_of_datatype, gpuMemcpyDeviceToDevice)
!         call check_memcpy_GPU_f("tridiag a_dev 1", 731,  successGPU)

! #ifdef WITH_NVTX
!         call nvtxRangePop()
! #endif

!       else ! useGPU
!         v_row(1:l_rows) = a_mat(1:l_rows,l_cols+1)
!       endif ! useGPU

      ! copy l_cols + 1 column of A to v_row
      if (.not. useGPU) then
        v_row(1:l_rows) = a_mat(1:l_rows,l_cols+1)
      endif ! useGPU

      if (n_stored_vecs > 0 .and. l_rows > 0) then
        if (useGPU) then
          if (wantDebug) call obj%timer%start("gpublas gemv skinny with copying")
          ! v_row_dev = vu_stored_rows_dev * uv_stored_cols_dev(l_cols+1,1:2*n_stored_vecs) + v_row_dev
          gpuHandle = obj%gpu_setup%gpublasHandleArray(0)
          call gpublas_ZGEMV('N', l_rows, 2*n_stored_vecs,  &
                                    ONE, vu_stored_rows_dev, max_local_rows,  &
                                    aux_complex_dev, 1,   &
                                    ONE, v_row_dev, 1, gpuHandle)

          if (wantDebug .and. gpu_vendor() /= SYCL_GPU) successGPU = gpu_DeviceSynchronize()
          if (wantDebug) call obj%timer%stop("gpublas gemv skinny with copying")

        else ! useGPU
          aux(1:2*n_stored_vecs) = conjg(uv_stored_cols(l_cols+1,1:2*n_stored_vecs))

          if (wantDebug) call obj%timer%start("blas")
          ! v_row = vu_stored_rows * uv_stored_cols(l_cols+1,1:2*n_stored_vecs) + v_row
          call ZGEMV('N',   &
                            int(l_rows,kind=BLAS_KIND), int(2*n_stored_vecs,kind=BLAS_KIND), &
                            ONE, vu_stored_rows, int(ubound(vu_stored_rows,dim=1),kind=BLAS_KIND), &
                            aux, 1_BLAS_KIND,  &
                            ONE, v_row, 1_BLAS_KIND)
          if (wantDebug) call obj%timer%stop("blas")
        endif ! useGPU
      endif ! (n_stored_vecs > 0 .and. l_rows > 0)

      if (useGPU) then
        if (my_prow == prow(istep-1, nblk, np_rows)) then
          isOurProcessRowInt = 1
        else
          isOurProcessRowInt = 0
        end if

        my_stream = obj%gpu_setup%my_stream
        call gpu_dot_product_and_assign_double_complex(v_row_dev, l_rows, isOurProcessRowInt, aux1_dev, wantDebug, my_stream)
        if (.not. useCCL) then


          successGPU = gpu_memcpy(int(loc(aux1),kind=c_intptr_t), aux1_dev, 2*size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("tridiag: aux1_dev -> aux1", 819,  successGPU)


        endif ! .not. useCCL 
      endif ! useGPU  


      if (.not. useGPU) then
        if (my_prow == prow(istep-1, nblk, np_rows)) then
          aux1(1) = dot_product(v_row(1:l_rows-1),v_row(1:l_rows-1)) ! = "q"
          aux1(2) = v_row(l_rows) ! = "a_11" (or rather a_nn)
        else
          aux1(1) = dot_product(v_row(1:l_rows),v_row(1:l_rows))
          aux1(2) = 0.
        endif
      endif ! .not. useGPU 


      if (useCCL) then
      else ! useCCL
        vnorm2 = real(aux1(1),kind=rk)
        vrl    = aux1(2)
      endif ! useCCL


      ! Householder transformation
      if (useCCL) then
      else ! useCCL
        call hh_transform_complex_&
                &double &
                (obj, vrl, vnorm2, xf, tau(istep), wantDebug)
      endif ! useCCL
      
      ! vrl is newly computed off-diagonal element of the final tridiagonal matrix
      if (my_prow == prow(istep-1, nblk, np_rows)) then
        if (.not. useCCL) then
          e_vec(istep-1) = real(vrl,kind=rk)
        endif ! .not. useCCL
      endif

      if (.not. useGPU) then
        ! Scale v_row and store Householder Vector for back transformation
        v_row(1:l_rows) = v_row(1:l_rows) * xf

        if (my_prow == prow(istep-1, nblk, np_rows)) then
          v_row(l_rows) = 1.
        endif

        ! store Householder Vector for back transformation
        ! update a_mat
        a_mat(1:l_rows,l_cols+1) = v_row(1:l_rows)
      endif ! .not. useGPU 

      if (.not. useCCL) then
        ! add tau after the end of actuall v_row, to be broadcasted with it
        v_row(l_rows+1) = tau(istep)
      endif

      
      if (useGPU) then
        if (useCCL) then
          xf_host_or_dev = xf_dev
        else 
          xf_host_or_dev = int(loc(xf),kind=c_intptr_t)
        endif       
        
        isOurProcessRow = (my_prow == prow(istep-1, nblk, np_rows))
        call gpu_set_e_vec_scale_set_one_store_v_row_double_complex(e_vec_dev, vrl_dev, a_dev, v_row_dev, tau_dev, xf_host_or_dev, & 
                                                  l_rows, l_cols, matrixRows, istep, isOurProcessRow, useCCL, wantDebug, my_stream)
      endif ! useGPU  


      if (useGPU .and. .not. useCCL) then      
        !v_row_dev -> v_row

        successGPU = gpu_memcpy(int(loc(v_row),kind=c_intptr_t), v_row_dev, l_rows*size_of_datatype, gpuMemcpyDeviceToHost)
        call check_memcpy_GPU_f("tridiag: v_row_dev -> v_row", 995,  successGPU)

      endif ! useGPU .and. .not. useCCL

    endif !(my_pcol == pcol(istep, nblk, np_cols))



    !recover tau, which has been broadcasted together with v_row
    if (.not. useCCL) then
      tau(istep) =  v_row(l_rows+1)
    endif

    ! Transpose Householder Vector v_row -> v_col
    if (useCCL) then
    else ! useCCL
      call elpa_transpose_vectors_&
          &complex&
          &_&
          &double &
                (obj, v_row, ubound(v_row,dim=1), mpi_comm_rows, v_col, ubound(v_col,dim=1), mpi_comm_cols, &
                1, istep-1, 1, nblk, max_threads, .true., success)
      if (wantDebug .and. gpu_vendor() /= SYCL_GPU) successGPU = gpu_DeviceSynchronize()
      if (.not.(success)) then
        write(error_unit,*) "Error in elpa_transpose_vectors. Aborting!"
        return
      endif
    endif ! useCCL

    ! Calculate u = (A + VU**T + UV**T)*v // Dongarra 1987: "y = (A - UV**T - VU**T)*u"

    ! For cache efficiency, we use only the upper half of the matrix tiles for this,
    ! thus the result is partly in u_col(:) and partly in u_row(:)

    if (.not. useGPU) then
      u_col(1:l_cols) = 0
      u_row(1:l_rows) = 0
    endif

    if (useGPU .and. useCCL) then
      successGPU = gpu_memset(u_col_dev, 0, l_cols * size_of_datatype) ! TODO_23_11: omit this, but change gpublas_gemm to u_col_dev=a_dev^T*v_row_dev+0*u_col_dev?
      call check_memcpy_GPU_f("tridiag: u_col_dev", 1101,  successGPU)
    endif

    if (l_rows>0 .and. l_cols>0) then
      if (useGPU) then

        if (.not. mat_vec_as_one_block) then

          successGPU = gpu_memcpy(v_col_dev, int(loc(v_col(1)),kind=c_intptr_t), &
                        l_cols * size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("tridiag: v_col_dev", 1146,  successGPU)


        endif ! .not. mat_vec_as_one_block

        if (.not. useCCL) then

          successGPU = gpu_memcpy(v_row_dev, int(loc(v_row(1)),kind=c_intptr_t), &
                                    l_rows * size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("tridiag: v_row_dev", 1170,  successGPU)


        endif ! .not. useCCL
      endif ! useGPU


      if (.not. useGPU) then
        do i=0, (istep-2)/tile_size ! iteration over tiles
          l_col_beg = i*l_cols_per_tile+1
          l_col_end = min(l_cols,(i+1)*l_cols_per_tile)
          if (l_col_end < l_col_beg) cycle
          do j = 0, i
            l_row_beg = j*l_rows_per_tile+1
            l_row_end = min(l_rows,(j+1)*l_rows_per_tile)
            if (l_row_end < l_row_beg) cycle

            ! multiplication by blocks is efficient only for CPU
            ! for GPU we introduced 2 other ways, either by stripes (more simmilar to the original
            ! CPU implementation) or by one large matrix Vector multiply
            if (wantDebug) call obj%timer%start("blas")
            ! u_col = a_mat*v_row + u_col(=0)
            call ZGEMV('C',  &
                        int(l_row_end-l_row_beg+1,kind=BLAS_KIND), int(l_col_end-l_col_beg+1,kind=BLAS_KIND), &
                        ONE, a_mat(l_row_beg, l_col_beg), int(matrixRows,kind=BLAS_KIND),         &
                        v_row(l_row_beg:max_local_rows+1), 1_BLAS_KIND,                           &
                        ONE, u_col(l_col_beg:max_local_cols), 1_BLAS_KIND)

            if (i/=j) then
              if (isSkewsymmetric) then
                call ZGEMV('N',int(l_row_end-l_row_beg+1,kind=BLAS_KIND), int(l_col_end-l_col_beg+1,kind=BLAS_KIND), &
                                    -ONE, a_mat(l_row_beg,l_col_beg), int(matrixRows,kind=BLAS_KIND),               &
                                    v_col(l_col_beg:max_local_cols), 1_BLAS_KIND, ONE, u_row(l_row_beg:max_local_rows), &
                                    1_BLAS_KIND)

              else
                ! u_row = a_mat*v_col + u_row(=0)
                call ZGEMV('N',int(l_row_end-l_row_beg+1,kind=BLAS_KIND), int(l_col_end-l_col_beg+1,kind=BLAS_KIND),  &
                                    ONE, a_mat(l_row_beg,l_col_beg), int(matrixRows,kind=BLAS_KIND),               &
                                    v_col(l_col_beg:max_local_cols), 1_BLAS_KIND, ONE, u_row(l_row_beg:max_local_rows), &
                                    1_BLAS_KIND)
              endif
            endif
            if (wantDebug) call obj%timer%stop("blas")

          enddo  ! j=0,i
        enddo  ! i=0,(istep-2)/tile_size
      endif ! .not. useGPU

      if (useGPU) then
        if (mat_vec_as_one_block) then
          ! Unlike for CPU, we (for each MPI thread) do just one large mat-vec multiplication
          ! this requires altering of the algorithm when later explicitly updating the matrix
          ! after max_stored_uv is reached : we need to update all tiles, not only those above diagonal
          if (wantDebug) call obj%timer%start("gpublas")
          
          gpuHandle = obj%gpu_setup%gpublasHandleArray(0) ! 0.3-0.5 us

              
          ! u_col_dev = a_dev^T*v_row_dev
          call gpublas_ZGEMV('C', l_rows,l_cols,  &
                                    ONE, a_dev, matrixRows,                   &
                                    v_row_dev , 1,                          &
                                    ZERO, u_col_dev, 1, gpuHandle)
              
          if (wantDebug .and. gpu_vendor() /= SYCL_GPU) successGPU = gpu_DeviceSynchronize()
              
       ! todo: try with non transposed!!!
!                 if(i/=j) then
!                   call gpublas_ZGEMV('N', l_row_end-l_row_beg+1,l_col_end-l_col_beg+1,  &
!                                             ONE, a_dev + offset_dev, matrixRows,                        &
!                                             v_col_dev + (l_col_beg - 1) *                      &
!                                             size_of_datatype, 1,                          &
!                                             ONE, u_row_dev + (l_row_beg - 1) *                 &
!                                             size_of_datatype, 1)
!                 endif
          if (wantDebug) call obj%timer%stop("gpublas")
        else  ! mat_vec_as_one_block
          !perform multiplication by stripes - it is faster than by blocks, since we call cublas with
          !larger matrices. In general, however, this algorithm is very simmilar to the one with CPU
          do i=0,(istep-2)/tile_size
            l_col_beg = i*l_cols_per_tile+1
            l_col_end = min(l_cols,(i+1)*l_cols_per_tile)
            if (l_col_end<l_col_beg) cycle

            l_row_beg = 1
            l_row_end = min(l_rows,(i+1)*l_rows_per_tile)
                  
            offset_dev = ((l_row_beg-1) + (l_col_beg - 1) * matrixRows) * size_of_datatype
  
            gpuHandle = obj%gpu_setup%gpublasHandleArray(0)
            call gpublas_ZGEMV('C', &
                          l_row_end-l_row_beg+1, l_col_end-l_col_beg+1, &
                          ONE, a_dev + offset_dev, matrixRows,  &
                          v_row_dev + (l_row_beg - 1) * size_of_datatype, 1,  &
                          ONE, u_col_dev + (l_col_beg - 1) * size_of_datatype, 1, gpuHandle)
          enddo !i=0,(istep-2)/tile_size

          do i=0,(istep-2)/tile_size
            l_col_beg = i*l_cols_per_tile+1
            l_col_end = min(l_cols,(i+1)*l_cols_per_tile)
            if (l_col_end<l_col_beg) cycle

            l_row_beg = 1
            l_row_end = min(l_rows,i*l_rows_per_tile)
              
            offset_dev = ((l_row_beg-1) + (l_col_beg - 1) * matrixRows) * size_of_datatype

            if (isSkewsymmetric) then
                gpuHandle = obj%gpu_setup%gpublasHandleArray(0)
                call gpublas_ZGEMV('N', l_row_end-l_row_beg+1, l_col_end-l_col_beg+1, &
                            -ONE, a_dev + offset_dev, matrixRows, &
                            v_col_dev + (l_col_beg - 1) * size_of_datatype,1, &
                            ONE, u_row_dev + (l_row_beg - 1) * size_of_datatype, 1, gpuHandle)
            else
                gpuHandle = obj%gpu_setup%gpublasHandleArray(0)
                call gpublas_ZGEMV('N', l_row_end-l_row_beg+1, l_col_end-l_col_beg+1, &
                            ONE, a_dev + offset_dev, matrixRows, &
                            v_col_dev + (l_col_beg - 1) * size_of_datatype,1, &
                            ONE, u_row_dev + (l_row_beg - 1) * size_of_datatype, 1, gpuHandle)
            endif
          enddo ! i=0,(istep-2)/tile_size
        end if ! mat_vec_as_one_block / per stripes


        if (.not. mat_vec_as_one_block) then

          successGPU = gpu_memcpy(int(loc(u_row(1)),kind=c_intptr_t), u_row_dev, l_rows*size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("tridiag: u_row_dev 1", 1378,  successGPU)


        endif ! .not. mat_vec_as_one_block
      endif ! useGPU


      ! second calculate (VU**T + UV**T)*v part of (A + VU**T + UV**T)*v
      if (n_stored_vecs > 0) then
        if (.not. useGPU) then
          if (wantDebug) call obj%timer%start("blas")

      
          ! aux = vu_stored_rows^T*v_row
          call ZGEMV('C',     &
                              int(l_rows,kind=BLAS_KIND), int(2*n_stored_vecs,kind=BLAS_KIND),   &
                              ONE, vu_stored_rows, int(ubound(vu_stored_rows,dim=1),kind=BLAS_KIND),   &
                              v_row,  1_BLAS_KIND, ZERO, aux, 1_BLAS_KIND)
              
          ! u_col = uv_stored_cols*aux + u_col
          call ZGEMV('N', int(l_cols,kind=BLAS_KIND), int(2*n_stored_vecs,kind=BLAS_KIND),   &
                              ONE, uv_stored_cols, int(ubound(uv_stored_cols,dim=1),kind=BLAS_KIND),   &
                              aux, 1_BLAS_KIND, ONE, u_col,  1_BLAS_KIND)
          if (wantDebug) call obj%timer%stop("blas")
        else ! .not. useGPU
          if (wantDebug) call obj%timer%start("gpublas gemv x2 skinny")


          ! aux_dev = vu_stored_rows_dev^T*v_row_dev
          gpuHandle = obj%gpu_setup%gpublasHandleArray(0)
          call gpublas_ZGEMV('C', l_rows, 2*n_stored_vecs, &
                                      ONE, vu_stored_rows_dev, max_local_rows,   &
                                      v_row_dev,  1, ZERO, aux_dev, 1, gpuHandle)
              
          ! u_col_dev = uv_stored_cols_dev*aux_dev + u_col_dev
          call gpublas_ZGEMV('N', l_cols, 2*n_stored_vecs, &
                                      ONE, uv_stored_cols_dev, max_local_cols,   &
                                      aux_dev, 1, ONE, u_col_dev, 1, gpuHandle)

          if (wantDebug .and. gpu_vendor() /= SYCL_GPU) successGPU = gpu_DeviceSynchronize()

          if (wantDebug) call obj%timer%stop("gpublas gemv x2 skinny")  
        endif ! .not. useGPU
      endif ! n_stored_vecs > 0

    endif  ! (l_rows>0 .and. l_cols>0)

    if (useGPU .and. l_cols>0 .and. (.not. useCCL)) then

      successGPU = gpu_memcpy(int(loc(u_col(1)),kind=c_intptr_t), u_col_dev, l_cols*size_of_datatype, gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("tridiag: u_col_dev 1", 1462,  successGPU)

    endif ! useGPU


    if (useGPU .and. .not. useCCL) then 
        if (l_rows==0) then
        u_col(1:l_cols) = 0
      endif
    endif ! useGPU

    ! Sum up all u_row(:) parts along rows and add them to the u_col(:) parts
    ! on the processors containing the diagonal
    ! This is only necessary if u_row has been calculated, i.e. if the
    ! global tile size is smaller than the global remaining matrix
    
    ! in GPU case, u_row_dev=0 up to now, so elpa_reduce_add_vectors is skipped

    if (tile_size < istep-1 .and. .not. useGPU) then
      call elpa_reduce_add_vectors_&
            &complex&
            &_&
            &double &
            (obj, u_row, ubound(u_row,dim=1), mpi_comm_rows, u_col, ubound(u_col,dim=1), &
            mpi_comm_cols, istep-1, 1, nblk, max_threads)
    endif ! (tile_size < istep-1 .and. .not. useGPU)

    ! Sum up all the u_col(:) parts, transpose u_col -> u_row

    if (l_cols > 0) then
    endif ! (l_cols > 0)
    
    ! Transpose Householder Vector u_col -> u_row
    if (useCCL) then
    else ! useCCL
      if (isSkewsymmetric) then
        call elpa_transpose_vectors_ss_&
        &complex&
        &_&
        &double &
        (obj, u_col, ubound(u_col,dim=1), mpi_comm_cols, u_row, ubound(u_row,dim=1), &
          mpi_comm_rows, 1, istep-1, 1, nblk, max_threads, .false., success)
        if (.not.(success)) then
          write(error_unit,*) "Error in elpa_transpose_vectors_ss. Aborting!"
          return
        endif
      else ! isSkewsymmetric
        call elpa_transpose_vectors_&
        &complex&
        &_&
        &double &
        (obj, u_col, ubound(u_col,dim=1), mpi_comm_cols, u_row, ubound(u_row,dim=1), &
        mpi_comm_rows, 1, istep-1, 1, nblk, max_threads, .false., success)
        if (.not.(success)) then
          write(error_unit,*) "Error in elpa_transpose_vectors. Aborting!"
          return
        endif
      endif ! isSkewsymmetric
    endif ! useCCL


    if (useGPU .and. .not. useCCL) then

      successGPU = gpu_memcpy(v_col_dev, int(loc(v_col(1)),kind=c_intptr_t), l_cols*size_of_datatype, gpuMemcpyHostToDevice)
      call check_memcpy_GPU_f("tridiag: v_col_dev", 1607,  successGPU)


      successGPU = gpu_memcpy(u_col_dev, int(loc(u_col(1)),kind=c_intptr_t), l_cols*size_of_datatype, gpuMemcpyHostToDevice)
      call check_memcpy_GPU_f("tridiag: u_col_dev", 1622,  successGPU)

    endif ! (useGPU .and. .not. useCCL)



    if (useGPU .and. .not. useCCL) then

      successGPU = gpu_memcpy(u_row_dev, int(loc(u_row(1)),kind=c_intptr_t), l_rows*size_of_datatype, gpuMemcpyHostToDevice)
      call check_memcpy_GPU_f("tridiag: u_row_dev", 1657,  successGPU)


      successGPU = gpu_memcpy(v_row_dev, int(loc(v_row(1)),kind=c_intptr_t), l_rows*size_of_datatype, gpuMemcpyHostToDevice)
      call check_memcpy_GPU_f("tridiag: v_row_dev", 1671,  successGPU)

    endif ! (useGPU .and. .not. useCCL)


    ! calculate u**T * v (same as v**T * (A + VU**T + UV**T) * v )
    if (.not. useGPU .or. (useGPU .and. .not. useCCL)) then
      vav = 0 ! x=0
      if (l_cols>0) vav = dot_product(v_col(1:l_cols), u_col(1:l_cols))
    endif



    ! store u and v in the matrices U and V
    ! these matrices are stored combined in one here

   if (.not. useCCL) then
      conjg_tau = conjg(tau(istep))
    endif ! .not. useCCL
    
    if (.not. useGPU) then
      if (l_rows > 0) then
        ! update vu_stored_rows
        vu_stored_rows(1:l_rows,2*n_stored_vecs+1) = conjg_tau*v_row(1:l_rows)
        vu_stored_rows(1:l_rows,2*n_stored_vecs+2) = 0.5*conjg_tau*vav*v_row(1:l_rows) - u_row(1:l_rows)
      endif
      if (l_cols > 0) then
        ! update uv_stored_cols
        uv_stored_cols(1:l_cols,2*n_stored_vecs+1) = 0.5*conjg_tau*vav*v_col(1:l_cols) - u_col(1:l_cols)
        uv_stored_cols(1:l_cols,2*n_stored_vecs+2) = conjg_tau*v_col(1:l_cols)
      endif
    endif ! .not. useGPU


    if (useGPU) then

      ! kernel: update vu_stored_rows_dev, uv_stored_cols
      ! then cpu's "store u,v in U,V" can be deleted. But we should take care of dot_prod below, where vu_stored_rows and uv_stored_cols are used
      if (useCCL) then
        vav_host_or_dev = vav_dev
        !tau_istep_host_or_dev = tau_dev + (istep-1)*size_of_datatype
        tau_istep_host_or_dev = v_row_dev + (l_rows+1-1)*size_of_datatype
      else ! useCCL
        vav_host_or_dev = int(loc(vav),kind=c_intptr_t)
        tau_istep_host_or_dev = int(loc(tau(istep)), kind=c_intptr_t)
      endif ! useCCL
      call gpu_store_u_v_in_uv_vu_double_complex(vu_stored_rows_dev, uv_stored_cols_dev, v_row_dev, u_row_dev, &
                                          v_col_dev, u_col_dev, tau_dev, aux_complex_dev, &
                                          vav_host_or_dev, tau_istep_host_or_dev, &
                                          l_rows, l_cols, n_stored_vecs,  max_local_rows, max_local_cols, istep, &
                                          useCCL, wantDebug, my_stream)
    endif ! useGPU


    ! We have calculated another Householder Vector, number of implicitly stored increased
    n_stored_vecs = n_stored_vecs+1

    ! If the limit of max_stored_uv is reached, calculate A + VU**T + UV**T
    if (n_stored_vecs == max_stored_uv .or. istep == 3) then        
      
      if (.not. useGPU .OR. .not. mat_vec_as_one_block) then
        do i = 0, (istep-2)/tile_size
          ! go over tiles above (or on) the diagonal
          l_col_beg = i*l_cols_per_tile+1
          l_col_end = min(l_cols,(i+1)*l_cols_per_tile)
          l_row_beg = 1
          l_row_end = min(l_rows,(i+1)*l_rows_per_tile)
          if (l_col_end<l_col_beg .or. l_row_end<l_row_beg) then
            cycle
          endif

          if (useGPU) then
            if (.not. mat_vec_as_one_block) then
              ! if using mat-vec multiply by stripes, it is enough to update tiles above (or on) the diagonal only
              ! we than use the same calls as for CPU version
              if (wantDebug) call obj%timer%start("gpublas_gemm")
              ! update a_dev
              ! a_dev = vu_stored_rows_dev*uv_stored_cols_dev + a_dev
              gpuHandle = obj%gpu_setup%gpublasHandleArray(0)
              call gpublas_ZGEMM('N', 'C',     &
                                      l_row_end-l_row_beg+1, l_col_end-l_col_beg+1, 2*n_stored_vecs,                      &
                                      ONE, vu_stored_rows_dev + (l_row_beg - 1) *                                         &
                                      size_of_datatype,  &
                                      max_local_rows, uv_stored_cols_dev + (l_col_beg - 1) *                              &
                                      size_of_datatype,  &
                                      max_local_cols, ONE, a_dev + ((l_row_beg - 1) + (l_col_beg - 1) * matrixRows) *     &
                                      size_of_datatype , matrixRows, gpuHandle)
              if (wantDebug .and. gpu_vendor() /= SYCL_GPU) successGPU = gpu_DeviceSynchronize()
              if (wantDebug) call obj%timer%stop("gpublas_gemm")
            endif ! .not. mat_vec_as_one_block
          else ! useGPU
            if (wantDebug) call obj%timer%start("blas_gemm")
            ! update a_mat
            ! a_mat = vu_stored_rows*uv_stored_cols + a_mat
            call ZGEMM('N', 'C',                &
                                int(l_row_end-l_row_beg+1,kind=BLAS_KIND), int(l_col_end-l_col_beg+1,kind=BLAS_KIND), &
                                int(2*n_stored_vecs,kind=BLAS_KIND),    &
                                ONE, vu_stored_rows(l_row_beg:max_local_rows,1:2*max_stored_uv), &
                                int(ubound(vu_stored_rows,dim=1),kind=BLAS_KIND),   &
                                uv_stored_cols(l_col_beg,1), &
                                int(ubound(uv_stored_cols,dim=1),kind=BLAS_KIND),        &
                                ONE, a_mat(l_row_beg,l_col_beg), int(matrixRows,kind=BLAS_KIND))
            if (wantDebug) call obj%timer%stop("blas_gemm")
          endif !useGPU
        enddo ! i = 0, (istep-2)/tile_size

      else !.not. useGPU or .not. mat_vec_as_one_block (i.e. useGPU and mat_vec_as_one_block)

        !update whole (remaining) part of matrix, including tiles below diagonal
        !we can do that in one large cublas call
        if (wantDebug) call obj%timer%start("gpublas_gemm")


        gpuHandle = obj%gpu_setup%gpublasHandleArray(0)
        ! update a_dev
        ! a_dev = vu_stored_rows_dev*uv_stored_cols_dev + a_dev
        call gpublas_ZGEMM('N', 'C', l_rows, l_cols, 2*n_stored_vecs,   &
                                  ONE, vu_stored_rows_dev, max_local_rows, &
                                  uv_stored_cols_dev, max_local_cols,  &
                                  ONE, a_dev, matrixRows, gpuHandle)

        if (wantDebug .and. gpu_vendor() /= SYCL_GPU) successGPU = gpu_DeviceSynchronize()
        if (wantDebug) call obj%timer%stop("gpublas_gemm")
        
        ! copy real(a_dev(l_rows,l_cols)) -> d_vec_dev(istep-1) for correct initial value of d_vec_dev before atomicAdd
        if ((.not. isSkewsymmetric) .and. &
            (my_prow == prow(istep-1, nblk, np_rows) .and. my_pcol == pcol(istep-1, nblk, np_cols))) then
          offset_dev = ((l_rows-1) + (l_cols-1)*matrixRows) * size_of_datatype
          successGPU = gpu_memcpy(d_vec_dev + (istep-2)*size_of_datatype_real, &
                                  a_dev + offset_dev, 1*size_of_datatype_real, gpuMemcpyDeviceToDevice)
          call check_memcpy_GPU_f("tridiag a_dev->d_vec_dev", 1863,  successGPU)
        endif

      endif !.not. useGPU or .not. mat_vec_as_one_block

      n_stored_vecs = 0
    endif ! (n_stored_vecs == max_stored_uv .or. istep == 3)

    if (my_prow == prow(istep-1, nblk, np_rows) .and. my_pcol == pcol(istep-1, nblk, np_cols)) then

      if (useGPU) then
        my_stream = obj%gpu_setup%my_stream
        call gpu_update_matrix_element_add_double_complex(vu_stored_rows_dev, uv_stored_cols_dev, a_dev, d_vec_dev, &
                                            l_rows, l_cols, matrixRows, max_local_rows, max_local_cols, istep, n_stored_vecs, &
                                            isSkewsymmetric, wantDebug, my_stream)

      else ! useGPU
        if (n_stored_vecs > 0) then
          ! update a_mat (only one elememt!)
          dot_prod = dot_product(vu_stored_rows(l_rows,1:2*n_stored_vecs), uv_stored_cols(l_cols,1:2*n_stored_vecs))
          a_mat(l_rows,l_cols) = a_mat(l_rows,l_cols) + dot_prod
        endif
        d_vec(istep-1) = real(a_mat(l_rows,l_cols),kind=rk)
      endif ! useGPU

    endif ! (my_prow == prow(istep-1, nblk, np_rows) .and. my_pcol == pcol(istep-1, nblk, np_cols))


  enddo ! main cycle over istep=na,3,-1


  if (useGPU) then
    ! copy a_dev -> a_mat for backtransformation
    num = matrixRows * matrixCols * size_of_datatype
    successGPU = gpu_memcpy(int(loc(a_mat(1,1)),kind=c_intptr_t), a_dev, num, gpuMemcpyDeviceToHost)
    call check_memcpy_GPU_f("tridiag: a_dev", 1921,  successGPU)

  endif ! useGPU


  ! Store e_vec(1) and d_vec(1)

  if (my_pcol==pcol(2, nblk, np_cols)) then
    if (my_prow==prow(1, nblk, np_rows)) then
      ! We use last l_cols value of loop above
      if (useGPU) then
        successGPU = gpu_memcpy(int(loc(aux3(1)),kind=c_intptr_t), a_dev + (matrixRows * (l_cols - 1)) * size_of_datatype, &
                                1 * size_of_datatype, gpuMemcpyDeviceToHost)
        call check_memcpy_GPU_f("tridiag: a_dev 5", 1944,  successGPU)
        vrl = aux3(1)
      else !useGPU
        vrl = a_mat(1,l_cols)
      endif !useGPU

      call hh_transform_complex_&
            &double &
            (obj, vrl, 0.0_rk, xf, tau(2), wantDebug)
      e_vec(1) = real(vrl,kind=rk)
      a_mat(1,l_cols) = 1. ! for consistency only
    endif ! (my_prow==prow(1, nblk, np_rows))


  endif ! (my_pcol==pcol(2, nblk, np_cols))


  if (my_prow == prow(1, nblk, np_rows) .and. my_pcol == pcol(1, nblk, np_cols))  then
    if (useGPU) then
      successGPU = gpu_memcpy(int(loc(aux3(1)),kind=c_intptr_t), a_dev, &
                             1 * size_of_datatype, gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("tridiag: a_dev 6", 2002,  successGPU)
      d_vec(1) = DREAL(aux3(1))
    else !useGPU
      d_vec(1) = DREAL(a_mat(1,1))
    endif !useGPU
  endif ! (my_prow == prow(1, nblk, np_rows) .and. my_pcol == pcol(1, nblk, np_cols))


  if (useGPU) then
    offset_dev = 1 * size_of_datatype_real
    ! first and last elements of d_vec are treated separately
    successGPU = gpu_memcpy(int(loc(d_vec(2)),kind=c_intptr_t), &
                            d_vec_dev + offset_dev, (na-2) * size_of_datatype_real, gpuMemcpyDeviceToHost)
    call check_memcpy_GPU_f("tridiag: d_vec", 2062,  successGPU)

    if (useCCL) then
      ! e_vec(1) is treated separately
      offset_dev = 1 * size_of_datatype_real
      successGPU = gpu_memcpy(int(loc(e_vec(2)),kind=c_intptr_t), &
                              e_vec_dev + offset_dev, (na-1) * size_of_datatype_real, gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("tridiag: e_vec", 2069,  successGPU)

      ! tau(2) is treated separately, tau(1) is not used
      offset_dev = 2 * size_of_datatype
      successGPU = gpu_memcpy(int(loc(tau(3)),kind=c_intptr_t), &
                              tau_dev + offset_dev, (na-2) * size_of_datatype, gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("tridiag: tau", 2075,  successGPU)
    endif

    ! todo: should we leave a_mat on the device for further use?
    successGPU = gpu_free(a_dev)
    call check_dealloc_GPU_f("tridiag: a_dev 9", 2080,  successGPU)

    successGPU = gpu_free(v_row_dev)
    call check_dealloc_GPU_f("tridiag: v_row_dev", 2083,  successGPU)

    successGPU = gpu_free(u_row_dev)
    call check_dealloc_GPU_f("tridiag: (u_row_dev", 2086,  successGPU)

    successGPU = gpu_free(v_col_dev)
    call check_dealloc_GPU_f("tridiag: v_col_dev", 2089,  successGPU)

    successGPU = gpu_free(u_col_dev)
    call check_dealloc_GPU_f("tridiag: u_col_dev ", 2092,  successGPU)

    successGPU = gpu_free(vu_stored_rows_dev)
    call check_dealloc_GPU_f("tridiag: vu_stored_rows_dev ", 2095,  successGPU)

    successGPU = gpu_free(uv_stored_cols_dev)
    call check_dealloc_GPU_f("tridiag:uv_stored_cols_dev ", 2098,  successGPU)

    successGPU = gpu_free(d_vec_dev)
    call check_dealloc_GPU_f("tridiag: d_vec_dev", 2101,  successGPU)

    successGPU = gpu_free(e_vec_dev)
    call check_dealloc_GPU_f("tridiag: e_vec_dev", 2104,  successGPU)

    successGPU = gpu_free(tau_dev)
    call check_dealloc_GPU_f("tridiag: tau_dev", 2107,  successGPU)

    successGPU = gpu_free(aux_dev)
    call check_dealloc_GPU_f("tridiag: aux_dev", 2110,  successGPU)

    successGPU = gpu_free(aux1_dev)
    call check_dealloc_GPU_f("tridiag: aux1_dev", 2113,  successGPU)

    successGPU = gpu_free(aux_complex_dev)
    call check_dealloc_GPU_f("tridiag: aux_complex_dev", 2116,  successGPU)

    successGPU = gpu_free(vav_dev)
    call check_dealloc_GPU_f("tridiag: vav_dev", 2119,  successGPU)

    successGPU = gpu_free(dot_prod_dev)
    call check_dealloc_GPU_f("tridiag: dot_prod_dev", 2122,  successGPU)

    successGPU = gpu_free(xf_dev)
    call check_dealloc_GPU_f("tridiag: xf_dev", 2125,  successGPU)

  endif ! useGPU

  ! distribute the arrays d_vec and e_vec to all processors

  allocate(tmp_real(na), stat=istat, errmsg=errorMessage)
  call check_allocate_f("tridiag: tmp_real", 2138,  istat,  errorMessage)


  deallocate(tmp_real, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag: tmp_real", 2187,  istat,  errorMessage)

  if (useGPU) then

  else ! useGPU
    deallocate(v_row, v_col, u_row, u_col, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("tridiag: v_row, v_col, u_row, u_col", 2260,  istat,  errorMessage)
  endif ! useGPU

  deallocate(vu_stored_rows, uv_stored_cols, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag: vu_stored_rows, uv_stored_cols", 2264,  istat,  errorMessage)

  deallocate(aux, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag: aux", 2267,  istat,  errorMessage)



  call obj%timer%stop("tridiag_&
  &complex&
  &" // &
  "_double" // &
  gpuString )

end subroutine tridiag_&
&complex&
&_&
&double





!> \brief Transforms the eigenvectors of a tridiagonal matrix back
!>                     to the eigenvectors of the original matrix
!>                     (like Scalapack Routine PDORMTR)
!>
!  Parameters
!
!> \param na          Order of matrix a_mat, number of rows of matrix q_mat
!>
!> \param nqc         Number of columns of matrix q_mat
!>
!> \param a_mat(lda,matrixCols)  Matrix containing the Householder vectors (i.e. matrix a after tridiag_real)
!>                           Distribution is like in Scalapack.
!>
!> \param lda         Leading dimension of a_mat
!>
!> \param tau(na)     Factors of the Householder vectors
!>
!> \param q_mat           On input: Eigenvectors of tridiagonal matrix
!>                    On output: Transformed eigenvectors
!>                    Distribution is like in Scalapack.
!>
!> \param ldq         Leading dimension of q_mat
!>
!> \param nblk        blocksize of cyclic distribution, must be the same in both directions!
!>
!> \param matrixCols  local columns of matrix a_mat and q_mat
!>
!> \param mpi_comm_rows        MPI-Communicator for rows
!>
!> \param mpi_comm_cols        MPI-Communicator for columns
!>
!> \param useGPU      If true,  GPU version of the subroutine will be used
!>


subroutine trans_ev_&
&complex&
&_&
&double &
(obj, na, nqc, a_mat, lda, tau, q_mat, ldq, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols, useGPU, success)
  use, intrinsic :: iso_c_binding
  use precision
  use elpa_abstract_impl
  use elpa_blas_interfaces
  use elpa_gpu
  use elpa_gpu_util

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
  class(elpa_abstract_impl_t), intent(inout)    :: obj
  integer(kind=ik), intent(in)                  :: na, nqc, lda, ldq, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols
  complex(kind=rck), intent(in)           :: tau(na)

  complex(kind=rck), intent(inout)        :: a_mat(lda,*)
  complex(kind=rck), intent(inout)        :: q_mat(ldq,*)
  logical, intent(in)                           :: useGPU
  integer(kind=ik)                              :: max_stored_rows, max_stored_rows_fac

  integer(kind=ik)                              :: my_prow, my_pcol, np_rows, np_cols
  integer(kind=MPI_KIND)                        :: mpierr, my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI
  integer(kind=ik)                              :: totalblocks, max_blocks_row, max_blocks_col, max_local_rows, max_local_cols
  integer(kind=ik)                              :: l_cols, l_rows, l_colh, nstor
  integer(kind=ik)                              :: istep, n, nc, ic, ics, ice, nb, cur_pcol
  integer(kind=ik)                              :: hvn_ubnd, hvm_ubnd

  complex(kind=rck), allocatable          :: hvb(:), hvm(:,:)
  complex(kind=rck), pointer              :: tmp1(:), tmp2(:)
  complex(kind=rck), allocatable          :: h1(:), h2(:), tmp_debug(:)
  complex(kind=rck), pointer              :: tmat(:,:)
  complex(kind=rck), pointer              :: hvm1(:)
  type(c_ptr)                                   :: tmp1_host, tmp2_host
  type(c_ptr)                                   :: hvm1_host, tmat_host

  integer(kind=ik)                              :: istat
  character(200)                                :: errorMessage
  character(20)                                 :: gpuString

  integer(kind=c_intptr_t)                      :: num
  integer(kind=C_intptr_T)                      :: q_dev, tmp_dev, hvm_dev, tmat_dev

  integer(kind=ik)                              :: blockStep
  logical                                       :: successGPU
  integer(kind=c_intptr_t), parameter           :: size_of_datatype = size_of_&
                                                                      &double&
                                                                      &_&
                                                                      &complex
  integer(kind=ik)                              :: error
  integer(kind=MPI_KIND)                        :: bcast_request1, allreduce_request1, allreduce_request2
  logical                                       :: useNonBlockingCollectivesCols
  logical                                       :: useNonBlockingCollectivesRows
  integer(kind=c_int)                           :: non_blocking_collectives_rows, non_blocking_collectives_cols
  logical                                       :: success
  integer(kind=c_intptr_t)                      :: gpuHandle, my_stream


! NCCL requires streams
!#if defined (USE_CCL_TRANS_EV) && !defined(WITH_GPU_STREAMS)
!  successGPU = cuda_stream_create(obj%gpu_setup%my_stream)
!  if (.not.(successGPU)) then
!    print *,"Cannot create gpu stream handle"
!  endif
!  my_stream = obj%gpu_setup%my_stream
!#endif

  success = .true.

  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  call obj%timer%start("trans_ev_&
  &complex&
  &" // &
  &"_double" //&
  gpuString)

  call obj%get("nbc_row_elpa1_tridi_to_full", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    write(error_unit,*) "Problem setting option for non blocking collectives for rows in elpa1_tridi_to_full. Aborting..."
    call obj%timer%stop("trans_ev_&
    &complex&
    &" // &
    &"_double" //&
    gpuString)
    success = .false.
    return
  endif

  call obj%get("nbc_col_elpa1_tridi_to_full", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    write(error_unit,*) "Problem setting option for non blocking collectives for cols in elpa1_tridi_to_full. Aborting..."
    call obj%timer%stop("trans_ev_&
    &complex&
    &" // &
    &"_double" //&
    gpuString)
    success = .false.
    return
  endif

  if (non_blocking_collectives_rows .eq. 1) then
    useNonBlockingCollectivesRows = .true.
  else
    useNonBlockingCollectivesRows = .false.
  endif

  if (non_blocking_collectives_cols .eq. 1) then
    useNonBlockingCollectivesCols = .true.
  else
    useNonBlockingCollectivesCols = .false.
  endif

  my_prow = obj%mpi_setup%myRank_comm_rows
  my_pcol = obj%mpi_setup%myRank_comm_cols

  np_rows = obj%mpi_setup%nRanks_comm_rows
  np_cols = obj%mpi_setup%nRanks_comm_cols


  !call obj%timer%start("mpi_communication")
  !call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND) ,my_prowMPI, mpierr)
  !call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND) ,np_rowsMPI, mpierr)
  !call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND) ,my_pcolMPI, mpierr)
  !call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND) ,np_colsMPI, mpierr)

  !my_prow = int(my_prowMPI, kind=c_int)
  !np_rows = int(np_rowsMPI, kind=c_int)
  !my_pcol = int(my_pcolMPI, kind=c_int)
  !np_cols = int(np_colsMPI, kind=c_int)
  !call obj%timer%stop("mpi_communication")

  call obj%get("max_stored_rows",max_stored_rows_fac, error)

  totalblocks = (na-1)/nblk + 1
  max_blocks_row = (totalblocks-1)/np_rows + 1
  max_blocks_col = ((nqc-1)/nblk)/np_cols + 1  ! Columns of q_mat!

  max_local_rows = max_blocks_row*nblk
  max_local_cols = max_blocks_col*nblk

  max_stored_rows = (max_stored_rows_fac/nblk+1)*nblk
 
  if (.not.(useGPU)) then
    allocate(tmat(max_stored_rows,max_stored_rows), stat=istat, errmsg=errorMessage)
    call check_alloc("trans_ev_&
    &complex&
    &", "tmat", istat, errorMessage)

    allocate(tmp1(max_local_cols*max_stored_rows), stat=istat, errmsg=errorMessage)
    call check_alloc("trans_ev_&
    &complex&
    &", "tmp1", istat, errorMessage)

    allocate(tmp2(max_local_cols*max_stored_rows), stat=istat, errmsg=errorMessage)
    call check_alloc("trans_ev_&
    &complex&
    &", "tmp2", istat, errorMessage)
  endif

  allocate(h1(max_stored_rows*max_stored_rows), stat=istat, errmsg=errorMessage)
  call check_alloc("trans_ev_&
  &complex&
  &", "h1", istat, errorMessage)

  allocate(h2(max_stored_rows*max_stored_rows), stat=istat, errmsg=errorMessage)
  call check_alloc("trans_ev_&
  &complex&
  &", "h2", istat, errorMessage)

  allocate(hvb(max_local_rows*nblk), stat=istat, errmsg=errorMessage)
  call check_alloc("trans_ev_&
  &complex&
  &", "hvn", istat, errorMessage)

  allocate(hvm(max_local_rows,max_stored_rows), stat=istat, errmsg=errorMessage)
  call check_alloc("trans_ev_&
  &complex&
  &", "hvm", istat, errorMessage)

  hvm = 0   ! Must be set to 0 !!!
  hvb = 0   ! Safety only
  blockStep = nblk

  l_cols = local_index(nqc, my_pcol, np_cols, nblk, -1) ! Local columns of q_mat

  nstor = 0
  if (useGPU) then
    hvn_ubnd = 0
  endif

  ! In the complex case tau(2) /= 0
  if (my_prow == prow(1, nblk, np_rows)) then
    q_mat(1,1:l_cols) = q_mat(1,1:l_cols)*(ONE-tau(2))
  endif
 
  if (useGPU) then
    ! todo: this is used only for copying hmv to device.. it should be possible to go without it
    !allocate(hvm1(max_local_rows*max_stored_rows), stat=istat, errmsg=errorMessage)
    !call check_alloc("trans_ev_&
    !&complex&
    !&", "hvm1", istat, errorMessage)
    successGPU = gpu_malloc(tmat_dev, max_stored_rows * max_stored_rows * size_of_datatype)
    call check_alloc_GPU_f("trans_ev", 362,  successGPU)

    successGPU = gpu_malloc(hvm_dev, max_local_rows * max_stored_rows * size_of_datatype)
    call check_alloc_GPU_f("trans_ev", 365,  successGPU)

    successGPU = gpu_malloc(tmp_dev, max_local_cols * max_stored_rows * size_of_datatype)
    call check_alloc_GPU_f("trans_ev", 368,  successGPU)

    num = ldq * matrixCols * size_of_datatype
    successGPU = gpu_malloc(q_dev, num)
    call check_alloc_GPU_f("trans_ev", 372,  successGPU)
  

    successGPU = gpu_memcpy(q_dev, int(loc(q_mat(1,1)),kind=c_intptr_t), &
                  num, gpuMemcpyHostToDevice)
    call check_memcpy_GPU_f("trans_ev", 394,  successGPU)
  endif  ! useGPU

  do istep = 1, na, blockStep


    ics = MAX(istep,3)
    ice = MIN(istep+nblk-1,na)
    if (ice<ics) cycle

    cur_pcol = pcol(istep, nblk, np_cols)

    nb = 0
    do ic = ics, ice

      l_colh = local_index(ic  , my_pcol, np_cols, nblk, -1) ! Column of Householder Vector
      l_rows = local_index(ic-1, my_prow, np_rows, nblk, -1) ! # rows of Householder Vector


      if (my_pcol == cur_pcol) then
        hvb(nb+1:nb+l_rows) = a_mat(1:l_rows,l_colh)
        if (my_prow == prow(ic-1, nblk, np_rows)) then
          hvb(nb+l_rows) = 1.
        endif
      endif

      nb = nb+l_rows
    enddo


    nb = 0
    do ic = ics, ice
      l_rows = local_index(ic-1, my_prow, np_rows, nblk, -1) ! # rows of Householder Vector
      hvm(1:l_rows,nstor+1) = hvb(nb+1:nb+l_rows)
      if (useGPU) then
        hvm_ubnd = l_rows
      endif
      nstor = nstor+1
      nb = nb+l_rows
    enddo

    ! Please note: for smaller matix sizes (na/np_rows<=256), a value of 32 for nstor is enough!
    if (nstor+nblk > max_stored_rows .or. istep+nblk > na .or. (na/np_rows <= 256 .and. nstor >= 32)) then

      ! Calculate scalar products of stored vectors.
      ! This can be done in different ways, we use dsyrk or zherk

      tmat = 0
      call obj%timer%start("blas")
      if (l_rows>0) then
        call ZHERK('U', 'C',   &
                         int(nstor,kind=BLAS_KIND), int(l_rows,kind=BLAS_KIND), ONE, &
                         hvm, int(ubound(hvm,dim=1),kind=BLAS_KIND), ZERO, tmat, int(max_stored_rows,kind=BLAS_KIND))
      endif
      call obj%timer%stop("blas")
      nc = 0
      do n = 1, nstor-1
        h1(nc+1:nc+n) = tmat(1:n,n+1)
        nc = nc+n
      enddo

      if (nc > 0) h2 = h1

      ! Calculate triangular matrix T

      nc = 0
      tmat(1,1) = tau(ice-nstor+1)
      do n = 1, nstor-1
        call obj%timer%start("blas")
        call ZTRMV('L', 'C' , 'N', int(n,kind=BLAS_KIND), tmat, &
                            int(max_stored_rows,kind=BLAS_KIND), h2(nc+1), 1_BLAS_KIND)
        call obj%timer%stop("blas")

        tmat(n+1,1:n) = &
        -conjg(h2(nc+1:nc+n)) &
        *tau(ice-nstor+n+1)

        tmat(n+1,n+1) = tau(ice-nstor+n+1)
        nc = nc+n
      enddo

      if (useGPU) then
        ! todo: is this reshape really neccessary?
        hvm1(1:hvm_ubnd*nstor) = reshape(hvm(1:hvm_ubnd,1:nstor), (/ hvm_ubnd*nstor /))

        !hvm_dev(1:hvm_ubnd*nstor) = hvm1(1:hvm_ubnd*nstor)
        successGPU = gpu_memcpy(hvm_dev, int(loc(hvm1(1)),kind=c_intptr_t),   &
                      hvm_ubnd * nstor * size_of_datatype, gpuMemcpyHostToDevice)

        call check_memcpy_GPU_f("trans_ev", 546,  successGPU)

        !tmat_dev = tmat
        successGPU = gpu_memcpy(tmat_dev, int(loc(tmat(1,1)),kind=c_intptr_t),   &
                      max_stored_rows * max_stored_rows * size_of_datatype, gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("trans_ev", 551,  successGPU)
      endif


      ! Q = Q - V * T * V**T * Q

      if (l_rows > 0) then
        if (useGPU) then
          call obj%timer%start("gpublas")
          gpuHandle = obj%gpu_setup%gpublasHandleArray(0)
          call gpublas_ZGEMM('C', 'N',   &
                                   nstor, l_cols, l_rows, ONE, hvm_dev, hvm_ubnd,  &
                                   q_dev, ldq, ZERO, tmp_dev, nstor, gpuHandle)
          call obj%timer%stop("gpublas")
        else ! useGPU

          call obj%timer%start("blas")
          call ZGEMM('C', 'N',  &
                              int(nstor,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), &
                              int(l_rows,kind=BLAS_KIND), ONE, hvm, int(ubound(hvm,dim=1),kind=BLAS_KIND), &
                              q_mat, int(ldq,kind=BLAS_KIND), ZERO, tmp1, int(nstor,kind=BLAS_KIND))
          call obj%timer%stop("blas")
        endif ! useGPU

      else !l_rows>0

        if (useGPU) then
          if (gpu_vendor() /= OPENMP_OFFLOAD_GPU) then
            successGPU = gpu_memset(tmp_dev, 0, l_cols * nstor * size_of_datatype)
            call check_memcpy_GPU_f("trans_ev", 590,  successGPU)
          else
            allocate(tmp_debug(l_cols * nstor))
            tmp_debug(:) = 0.
            successGPU = gpu_memcpy(tmp_dev, int(loc(tmp_debug),kind=c_intptr_t), &
                                    l_cols*nstor*size_of_datatype, gpuMemcpyHostToDevice)
            call check_memcpy_GPU_f("trans_ev", 596,  successGPU)
            deallocate(tmp_debug)
          endif
        else
          tmp1(1:l_cols*nstor) = 0
        endif
      endif  !l_rows>0

!     tmp2 = tmp1

      if (l_rows > 0) then
        if (useGPU) then
          call obj%timer%start("gpublas")
          gpuHandle = obj%gpu_setup%gpublasHandleArray(0)
          call gpublas_ZTRMM('L', 'L', 'N', 'N',     &
                                   nstor, l_cols, ONE, tmat_dev, max_stored_rows,  &
                                   tmp_dev, nstor, gpuHandle)

          call gpublas_ZGEMM('N', 'N' ,l_rows ,l_cols ,nstor,  &
                                   -ONE, hvm_dev, hvm_ubnd, tmp_dev, nstor,   &
                                   ONE, q_dev, ldq, gpuHandle)
          call obj%timer%stop("gpublas")
        else !useGPU
          call obj%timer%start("blas")

          call ZTRMM('L', 'L', 'N', 'N', int(nstor,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND),   &
                              ONE, tmat, int(max_stored_rows,kind=BLAS_KIND), tmp1, int(nstor,kind=BLAS_KIND))
          call ZGEMM('N', 'N', int(l_rows,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), &
                              int(nstor,kind=BLAS_KIND), -ONE, hvm, int(ubound(hvm,dim=1),kind=BLAS_KIND), &
                              tmp1, int(nstor,kind=BLAS_KIND), ONE, q_mat, int(ldq,kind=BLAS_KIND))
          call obj%timer%stop("blas")
        endif ! useGPU
      endif  ! l_rows>0
      nstor = 0
    endif  ! (nstor+nblk>max_stored_rows .or. istep+nblk>na .or. (na/np_rows<=256 .and. nstor>=32))

    
  enddo ! istep = 1, na, blockStep

  deallocate(h1, h2, hvb, hvm, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_&     &MATH_DATATYPE&     &: h1, h2, hvb, hvm", 848,  istat,  errorMessage)

  if (useGPU) then

    !q_mat = q_dev
    successGPU = gpu_memcpy(int(loc(q_mat(1,1)),kind=c_intptr_t), &
                  q_dev, ldq * matrixCols * size_of_datatype, gpuMemcpyDeviceToHost)
    call check_memcpy_GPU_f("trans_ev", 864,  successGPU)

    !deallocate(hvm1, stat=istat, errmsg=errorMessage)
    !if (istat .ne. 0) then
    !  print *,"trans_ev_&
    !  &complex&
    !  &: error when deallocating hvm1 "//errorMessage
    !  stop 1
    !endif

    !deallocate(q_dev, tmp_dev, hvm_dev, tmat_dev)
    successGPU = gpu_free(q_dev)
    call check_dealloc_GPU_f("trans_ev", 906,  successGPU)

    successGPU = gpu_free(tmp_dev)
    call check_dealloc_GPU_f("trans_ev", 909,  successGPU)

    successGPU = gpu_free(hvm_dev)
    call check_dealloc_GPU_f("trans_ev", 912,  successGPU)

    successGPU = gpu_free(tmat_dev)
    call check_dealloc_GPU_f("trans_ev", 915,  successGPU)
  else ! useGPU
    deallocate(tmat, tmp1, tmp2, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("trans_ev_&     &MATH_DATATYPE&     &: tmat, tmp1, tmp2", 920,  istat,  errorMessage)
  endif ! useGPU

!#if defined (USE_CCL_TRANS_EV) && !defined(WITH_GPU_STREAMS)
!  successGPU = cuda_stream_destroy(obj%gpu_setup%my_stream)
!  if (.not.(successGPU)) then
!    print *,"Cannot destroy gpu stream handle"
!  endif
!#endif

  call obj%timer%stop("trans_ev_&
  &complex&
  &" // &
  &"_double" // &
  gpuString )

end subroutine trans_ev_&
&complex&
&_&
&double







subroutine hh_transform_complex_&
&double &
(obj, alpha, xnorm_sq, xf, tau, wantDebug)
  ! Similar to LAPACK routine ZLARFP, but uses ||x||**2 instead of x(:)
  ! and returns the factor xf by which x has to be scaled.
  ! It also hasn't the special handling for numbers < 1.d-300 or > 1.d150
  ! since this would be expensive for the parallel implementation.
  use precision
  use elpa_abstract_impl
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
  class(elpa_abstract_impl_t), intent(inout)    :: obj
  logical, intent(in)                           :: wantDebug
  complex(kind=ck), intent(inout) :: alpha
  real(kind=rk), intent(in)          :: xnorm_sq
  complex(kind=ck), intent(out)   :: xf, tau
  real(kind=rk)                      :: ALPHR, ALPHI

  real(kind=rk)                      :: BETA

  if (wantDebug) call obj%timer%start("hh_transform_&
                   &complex&
     	      &" // &
                   &"_double" )

  ALPHR = real( ALPHA, kind=rk )
  ALPHI = DIMAG( ALPHA )

  if ( XNORM_SQ==0.0_rk .AND. ALPHI==0.0_rk ) then

    if ( ALPHR>=0.0_rk ) then
      TAU = 0.0_rk
    else
      TAU = 2.0_rk
      ALPHA = -ALPHA
    endif
    XF = 0.0_rk

  else

    BETA = SIGN( SQRT( ALPHR**2 + ALPHI**2 + XNORM_SQ ), ALPHR )
    ALPHA = ALPHA + BETA
    IF ( BETA<0 ) THEN
      BETA = -BETA
      TAU  = -ALPHA / BETA
    ELSE
      ALPHR = ALPHI * (ALPHI/real( ALPHA , kind=rk))
      ALPHR = ALPHR + XNORM_SQ/real( ALPHA, kind=rk )

      TAU = DCMPLX( ALPHR/BETA, -ALPHI/BETA )
      ALPHA = DCMPLX( -ALPHR, ALPHI )
    END IF
    XF = 1.0_rk/ALPHA
    ALPHA = BETA
  endif

  if (wantDebug) call obj%timer%stop("hh_transform_&
  &complex&
  &" // &
  &"_double" )

    end subroutine hh_transform_complex_&
    &double



! complex single precision
































!> \brief Reduces a distributed symmetric matrix to tridiagonal form (like Scalapack Routine PDSYTRD)
!>
!  Parameters
!
!> \param obj	      object of elpa_type
!> \param na          Order of matrix
!>
!> \param a_mat(matrixRows,matrixCols)    Distributed matrix which should be reduced.
!>              Distribution is like in Scalapack.
!>              Opposed to PDSYTRD, a(:,:) must be set completely (upper and lower half)
!>              a(:,:) is overwritten on exit with the Householder vectors
!>
!> \param matrixRows         Leading dimension of a
!>
!> \param nblk        blocksize of cyclic distribution, must be the same in both directions!
!>
!> \param matrixCols  local columns of matrix
!>
!> \param mpi_comm_rows        MPI-Communicator for rows
!> \param mpi_comm_cols        MPI-Communicator for columns
!>
!> \param d_vec(na)       Diagonal elements (returned), identical on all processors
!>
!> \param e_vec(na)       Off-Diagonal elements (returned), identical on all processors
!>
!> \param tau(na)     Factors for the Householder vectors (returned), needed for back transformation
!>
!> \param useGPU      If true,  GPU version of the subroutine will be used
!> \param wantDebug   if true more debug information
!>
subroutine tridiag_&
  &complex&
  &_&
  &single &
  (obj, na, a_mat, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols, d_vec, e_vec, tau, useGPU, wantDebug, &
   max_threads_in, isSkewsymmetric, success)
  use, intrinsic :: iso_c_binding
  use precision
  use elpa_abstract_impl
  use matrix_plot
  use elpa_omp
  use elpa_blas_interfaces
  use elpa_gpu
  use elpa_gpu_util
  use tridiag_gpu

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
  class(elpa_abstract_impl_t), intent(inout)    :: obj
  integer(kind=ik), intent(in)                  :: na, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols
  logical, intent(in)                           :: useGPU, wantDebug
  logical, intent(in)                           :: isSkewsymmetric

  logical                                       :: useCCL=.false.

  complex(kind=rck), intent(out)          :: tau(na)
  complex(kind=rck), intent(inout)        :: a_mat(matrixRows,*)
  real(kind=rk), intent(out)                    :: d_vec(na)
  real(kind=rk), intent(out)                    :: e_vec(na)
  integer(kind=ik)                              :: max_stored_uv = 32 ! TODO_23_11 - make it tunable instead of hard-coded
  logical,          parameter                   :: mat_vec_as_one_block = .true.

  ! id in processor row and column and total numbers of processor rows and columns
  integer(kind=ik)                              :: my_prow, my_pcol, np_rows, np_cols
  integer(kind=MPI_KIND)                        :: my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI
  integer(kind=MPI_KIND)                        :: mpierr
  integer(kind=ik)                              :: totalblocks, max_loc_block_rows, max_loc_block_cols, max_local_rows, &
                                                   max_local_cols
  ! updated after each istep (in the main cycle) to contain number of
  ! local columns and rows of the remaining part of the matrix
  integer(kind=ik)                              :: l_cols, l_rows
  integer(kind=ik)                              :: n_stored_vecs
  integer(kind=ik)                              :: isOurProcessRowInt ! TODO_23_11 - get rid of it
  logical                                       :: isOurProcessRow, isOurProcessCol, isOurProcessCol_prev


  integer(kind=C_intptr_T)                      :: a_dev, v_row_dev, v_col_dev, u_row_dev, u_col_dev, vu_stored_rows_dev, &
                                                   uv_stored_cols_dev, d_vec_dev, e_vec_dev, tau_dev
  logical                                       :: successGPU

  integer(kind=ik)                              :: istep, i, j, l_col_beg, l_col_end, l_row_beg, l_row_end
  integer(kind=ik)                              :: tile_size, l_rows_per_tile, l_cols_per_tile
  integer(kind=c_intptr_t)                      :: offset_dev

  integer(kind=ik), intent(in)                  :: max_threads_in
  integer(kind=ik)                              :: max_threads

  real(kind=rk)                                 :: vnorm2
  complex(kind=rck)                       :: vav, x, aux1(2), vrl, xf, conjg_tau, dot_prod
  complex(kind=rck), allocatable          :: aux(:)

  integer(kind=c_intptr_t)                      :: aux_dev, aux1_dev, aux_complex_dev, vav_dev, vav_host_or_dev, dot_prod_dev, & 
                                                   xf_dev, a_updated_element_dev, xf_host_or_dev, tau_istep_host_or_dev
  integer(kind=c_intptr_t)                      :: vnorm2_dev, vrl_dev ! alias pointers for aux1_dev(1), aux1_dev(2)

  complex(kind=rck)                             :: aux3(1)

  integer(kind=c_intptr_t)                      :: num
  complex(kind=rck), pointer              :: v_row(:), & ! used to store calculated Householder Vector
                                                   v_col(:)   ! the same Vector, but transposed 
  complex(kind=rck), pointer              :: u_row_debug(:), & ! used to store calculated Householder Vector
                                                   u_col_debug(:)   ! the same Vector, but transposed 
  complex(kind=rck), pointer              :: u_col(:), u_row(:)

  ! the following two matrices store pairs of vectors v and u calculated in each step
  ! at most max_stored_uv Vector pairs are stored, than the matrix A_i is explicitli updated
  ! u and v are stored both in row and Vector forms
  ! pattern: v1,u1,v2,u2,v3,u3,....
  ! todo: It is little bit confusing, I think, that variables _row actually store columns and vice versa
  !complex(kind=rck), pointer             :: vu_stored_rows(:,:)
  complex(kind=rck), allocatable         :: vu_stored_rows(:,:)
  ! pattern: u1,v1,u2,v2,u3,v3,....
  complex(kind=rck), allocatable         :: uv_stored_cols(:,:)


  type(c_ptr)                                   :: v_row_host, v_col_host
  type(c_ptr)                                   :: u_row_host, u_col_host
  !type(c_ptr)                                   :: vu_stored_rows_host, uv_stored_cols_host
  real(kind=rk), allocatable                    :: tmp_real(:)
  integer(kind=ik)                              :: min_tile_size, error
  integer(kind=ik)                              :: istat
  character(200)                                :: errorMessage
  character(20)                                 :: gpuString
  integer(kind=ik)                              :: nblockEnd
  integer(kind=c_intptr_t), parameter           :: size_of_datatype = size_of_&
                                                                      &single&
                                                                      &_&
                                                                      &complex
  integer(kind=c_intptr_t), parameter           :: size_of_datatype_real = size_of_&
                                                                      &single&
                                                                      &_real
  integer(kind=MPI_KIND)                        :: bcast_request1, bcast_request2, bcast_request3
  integer(kind=MPI_KIND)                        :: allreduce_request1, allreduce_request2, allreduce_request3
  integer(kind=MPI_KIND)                        :: allreduce_request4, allreduce_request5, allreduce_request6, &
                                                   allreduce_request7
  logical                                       :: useNonBlockingCollectivesCols
  logical                                       :: useNonBlockingCollectivesRows
  integer(kind=c_int)                           :: non_blocking_collectives_rows, non_blocking_collectives_cols
  logical                                       :: success

  integer(kind=c_intptr_t)                      :: gpuHandle, my_stream
  integer(kind=c_int) :: pointerMode

  integer(kind=ik)                              :: string_length, sm_count


  allocate(aux(2*max_stored_uv), stat=istat, errmsg=errorMessage)

  success = .true.

  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  if (useGPU) then
    max_threads=1
  else
    max_threads=max_threads_in
  endif

  call obj%timer%start("tridiag_&
  &complex&
  &" // &
  "_single" // &
  gpuString )

  call obj%get("nbc_row_elpa1_full_to_tridi", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    write(error_unit,*) "Problem setting option for non blocking collectives for rows in elpa1_tridiag. Aborting..."
    success = .false.
    call obj%timer%stop("tridiag_&
    &complex&
    &" // &
    "_single" // &
    gpuString )
    return
  endif

  call obj%get("nbc_col_elpa1_full_to_tridi", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    write(error_unit,*) "Problem setting option for non blocking collectives for cols in elpa1_tridiag. Aborting..."
    success = .false.
    call obj%timer%stop("tridiag_&
    &complex&
    &" // &
    "_single" // &
    gpuString )
    return
  endif

  if (non_blocking_collectives_rows .eq. 1) then
    useNonBlockingCollectivesRows = .true.
  else
    useNonBlockingCollectivesRows = .false.
  endif

  if (non_blocking_collectives_cols .eq. 1) then
    useNonBlockingCollectivesCols = .true.
  else
    useNonBlockingCollectivesCols = .false.
  endif


  !mpi_comm_all    = obj%mpi_setup%mpi_comm_parent
  !mpi_comm_cols   = obj%mpi_setup%mpi_comm_cols
  !mpi_comm_rows   = obj%mpi_setup%mpi_comm_rows

  my_prow = obj%mpi_setup%myRank_comm_rows
  my_pcol = obj%mpi_setup%myRank_comm_cols

  np_rows = obj%mpi_setup%nRanks_comm_rows
  np_cols = obj%mpi_setup%nRanks_comm_cols


  !if (wantDebug) call obj%timer%start("mpi_communication")
  !call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND), my_prowMPI, mpierr)
  !call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND), np_rowsMPI, mpierr)
  !call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND), my_pcolMPI, mpierr)
  !call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND), np_colsMPI, mpierr)


  !my_prow = int(my_prowMPI, kind=c_int)
  !np_rows = int(np_rowsMPI, kind=c_int)
  !my_pcol = int(my_pcolMPI, kind=c_int)
  !np_cols = int(np_colsMPI, kind=c_int)
  !if (wantDebug) call obj%timer%stop("mpi_communication")

  ! Matrix is split into tiles; work is done only for tiles on the diagonal or above
  ! seems that tile is a square submatrix, consisting by several blocks
  ! it is a smallest possible square submatrix, where blocks being distributed among
  ! processors are "aligned" in both rows and columns
  !  -----------------
  ! | 1 4 | 1 4 | 1 4 | ...
  ! | 2 5 | 2 5 | 2 5 | ...
  ! | 3 6 | 3 6 | 3 6 | ...
  !  ----------------- ...
  ! | 1 4 | 1 4 | 1 4 | ...
  ! | 2 5 | 2 5 | 2 5 | ...
  ! | 3 6 | 3 6 | 3 6 | ...
  !  ----------------- .
  !   : :   : :   : :    .
  !   : :   : :   : :      .
  !
  ! this is a tile, where each number represents block, assigned to a processor with the shown number
  ! size of this small block is nblk
  ! Image is for situation with 6 processors, 3 processor rows and 2 columns
  ! tile_size is thus nblk * 6
  !
  tile_size = nblk*least_common_multiple(np_rows,np_cols) ! minimum global tile size

  ! make tile_size a smallest possible multiple of previously defined tile size, such that it is
  ! larger or equal to min_tile_size
  ! min_tile_size has been originally hardcoded as 128 * max(np_rows, np_cols), so it is now the implicit value
  ! it can, however, be set by the user
  call obj%get("min_tile_size", min_tile_size ,error)
  if (error .ne. ELPA_OK) then
    write(error_unit,*) "Problem setting option for min_tile_size. Aborting..."
    success = .false.
    call obj%timer%stop("tridiag_&
    &complex&
    &" // &
    "_single" // &
    gpuString )
    return
  endif
  if(min_tile_size == 0) then
    ! not set by the user, use the default value
    min_tile_size = 128*max(np_rows, np_cols)
  endif
  tile_size = ((min_tile_size-1)/tile_size+1)*tile_size

  nblockEnd = 3

  l_rows_per_tile = tile_size/np_rows ! local rows of a tile
  l_cols_per_tile = tile_size/np_cols ! local cols of a tile

  totalblocks = (na-1)/nblk + 1
  max_loc_block_rows = (totalblocks-1)/np_rows + 1
  max_loc_block_cols = (totalblocks-1)/np_cols + 1

  ! localy owned submatrix has size at most max_local_rows x max_local_cols at each processor
  max_local_rows = max_loc_block_rows*nblk
  max_local_cols = max_loc_block_cols*nblk

  ! allocate memmory for vectors
  ! todo: It is little bit confusing, I think, that variables _row actually store columns and vice versa
  ! todo: if something has length max_local_rows, it is actually a column, no?
  ! todo: probably one should read it as v_row = Vector v distributed among rows

  allocate(uv_stored_cols(max_local_cols,2*max_stored_uv), stat=istat, errmsg=errorMessage)
  call check_alloc("tridiag_&
       &complex ", "uv_stored_cols", istat, errorMessage)

  allocate(vu_stored_rows(max_local_rows,2*max_stored_uv), stat=istat, errmsg=errorMessage)
  call check_alloc("tridiag_&
       &complex ", "vu_stored_rows", istat, errorMessage)

  if (useGPU) then

    ! allocate v_row 1 element longer to allow store and broadcast tau together with it
    if (gpu_vendor() /= OPENMP_OFFLOAD_GPU .and. gpu_vendor() /= SYCL_GPU) then
      num = (max_local_rows+1) * size_of_datatype
      successGPU = gpu_malloc_host(v_row_host, num)
      call check_host_alloc_GPU_f("tridiag: v_row_host", 444,  successGPU)
      call c_f_pointer(v_row_host,v_row,(/(max_local_rows+1)/))
    else
      allocate(v_row(max_local_rows+1))
    endif

    if (gpu_vendor() /= OPENMP_OFFLOAD_GPU .and. gpu_vendor() /= SYCL_GPU) then
      num = (max_local_cols) * size_of_datatype
      successGPU = gpu_malloc_host(v_col_host,num)
      call check_host_alloc_GPU_f("tridiag: v_col_host", 453,  successGPU)
      call c_f_pointer(v_col_host,v_col,(/(max_local_cols)/))
    else
      allocate(v_col(max_local_cols))
    endif

    if (gpu_vendor() /= OPENMP_OFFLOAD_GPU .and. gpu_vendor() /= SYCL_GPU) then
      num = (max_local_cols) * size_of_datatype
      successGPU = gpu_malloc_host(u_col_host,num)
      call check_host_alloc_GPU_f("tridiag: u_col_host", 462,  successGPU)
      call c_f_pointer(u_col_host,u_col,(/(max_local_cols)/))
    else
      allocate(u_col(max_local_cols))
    endif

    if (gpu_vendor() /= OPENMP_OFFLOAD_GPU .and. gpu_vendor() /= SYCL_GPU) then
      num = (max_local_rows) * size_of_datatype
      successGPU = gpu_malloc_host(u_row_host,num)
      call check_host_alloc_GPU_f("tridiag: u_row_host", 471,  successGPU)
      call c_f_pointer(u_row_host,u_row,(/(max_local_rows)/))
    else
      allocate(u_row(max_local_rows))
    endif

    
    if (gpu_vendor() /= OPENMP_OFFLOAD_GPU .and. gpu_vendor() /= SYCL_GPU) then
      num = (max_local_rows * 2*max_stored_uv) * size_of_datatype
      successGPU = gpu_host_register(int(loc(vu_stored_rows),kind=c_intptr_t), num, gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: vu_stored_rows", 481,  successGPU)

      num = (max_local_cols * 2*max_stored_uv) * size_of_datatype
      successGPU = gpu_host_register(int(loc(uv_stored_cols),kind=c_intptr_t), num, gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: uv_stored_cols", 485,  successGPU)

      num = (1 * 2*max_stored_uv) * size_of_datatype
      successGPU = gpu_host_register(int(loc(aux),kind=c_intptr_t), num, gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: aux", 489,  successGPU)

      num = na * size_of_datatype_real
      successGPU = gpu_host_register(int(loc(d_vec),kind=c_intptr_t), num, gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: d_vec", 493,  successGPU)

      num = na * size_of_datatype_real
      successGPU = gpu_host_register(int(loc(e_vec),kind=c_intptr_t), num, gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: e_vec", 497,  successGPU)

      num = na * size_of_datatype
      successGPU = gpu_host_register(int(loc(tau),kind=c_intptr_t), num, gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: tau", 501,  successGPU)

      num = 2 * size_of_datatype
      successGPU = gpu_host_register(int(loc(aux1),kind=c_intptr_t), num, gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: aux1", 505,  successGPU)

      num = 1 * size_of_datatype
      successGPU = gpu_host_register(int(loc(vav),kind=c_intptr_t), num, gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: vav", 509,  successGPU)

      num = 1 * size_of_datatype
      successGPU = gpu_host_register(int(loc(dot_prod),kind=c_intptr_t), num, gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: dot_prod", 513,  successGPU)

      num = 1 * size_of_datatype
      successGPU = gpu_host_register(int(loc(xf),kind=c_intptr_t), num, gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: xf", 517,  successGPU)
    endif ! gpu_vendor() /= OPENMP_OFFLOAD_GPU .and. gpu_vendor() /= SYCL_GPU
  else ! useGPU

    allocate(v_row(max_local_rows+1), stat=istat, errmsg=errorMessage)
    call check_alloc("tridiag_&
    &complex ", "v_row", istat, errorMessage)

    allocate(v_col(max_local_cols), stat=istat, errmsg=errorMessage)
    call check_alloc("tridiag_&
     &complex ", "v_col", istat, errorMessage)

    allocate(u_col(max_local_cols), stat=istat, errmsg=errorMessage)
    call check_alloc("tridiag_&
    &complex ", "u_col", istat, errorMessage)

    allocate(u_row(max_local_rows), stat=istat, errmsg=errorMessage)
    call check_alloc("tridiag_&
    &complex ", "u_row", istat, errorMessage)
      
  endif ! useGPU


  v_row = 0
  u_row = 0
  v_col = 0
  u_col = 0

  if (useGPU) then
    successGPU = gpu_malloc(v_row_dev, (max_local_rows+1) * size_of_datatype)
    call check_alloc_GPU_f("tridiag: v_row_dev", 556,  successGPU)

    successGPU = gpu_malloc(u_row_dev, max_local_rows * size_of_datatype)

    call check_alloc_GPU_f("tridiag: u_row_dev", 560,  successGPU)

    successGPU = gpu_malloc(v_col_dev, max_local_cols * size_of_datatype)
    call check_alloc_GPU_f("tridiag: v_col_dev", 563,  successGPU)

    successGPU = gpu_malloc(u_col_dev, max_local_cols * size_of_datatype)
    call check_alloc_GPU_f("tridiag: u_col_dev", 566,  successGPU)

    successGPU = gpu_malloc(vu_stored_rows_dev, max_local_rows * 2 * max_stored_uv * size_of_datatype)
    call check_alloc_GPU_f("tridiag: vu_stored_rows_dev", 569,  successGPU)

    successGPU = gpu_malloc(uv_stored_cols_dev, max_local_cols * 2 * max_stored_uv * size_of_datatype)
    call check_alloc_GPU_f("tridiag: uv_stored_cols_dev", 572,  successGPU)

    successGPU = gpu_malloc(d_vec_dev, na * size_of_datatype_real)
    call check_alloc_GPU_f("tridiag: d_vec_dev", 575,  successGPU)

    successGPU = gpu_malloc(e_vec_dev, na * size_of_datatype_real)
    call check_alloc_GPU_f("tridiag: e_vec_dev", 578,  successGPU)

    successGPU = gpu_malloc(tau_dev, na * size_of_datatype)
    call check_alloc_GPU_f("tridiag: tau_dev", 581,  successGPU)

    successGPU = gpu_malloc(aux_dev, 2*max_stored_uv * size_of_datatype)
    call check_alloc_GPU_f("tridiag: aux_dev", 584,  successGPU)

    successGPU = gpu_malloc(aux1_dev, 2 * size_of_datatype)
    call check_alloc_GPU_f("tridiag: aux1_dev", 587,  successGPU)

    successGPU = gpu_malloc(aux_complex_dev, 2 *max_stored_uv* size_of_datatype)
    call check_alloc_GPU_f("tridiag: aux_complex_dev", 590,  successGPU)

    successGPU = gpu_malloc(vav_dev, 1 * size_of_datatype)
    call check_alloc_GPU_f("tridiag: vav_dev", 593,  successGPU)

    successGPU = gpu_malloc(dot_prod_dev, 1 * size_of_datatype)
    call check_alloc_GPU_f("tridiag: dot_prod_dev", 596,  successGPU)

    successGPU = gpu_malloc(xf_dev, 1 * size_of_datatype)
    call check_alloc_GPU_f("tridiag: xf_dev", 599,  successGPU)

  endif !useGPU

  d_vec(:) = 0
  e_vec(:) = 0
  tau(:) = 0

  if (useGPU) then
    successGPU = gpu_memset(d_vec_dev, 0, na * size_of_datatype_real)
    call check_memcpy_GPU_f("tridiag: d_vec_dev", 630,  successGPU)

    successGPU = gpu_memset(e_vec_dev, 0, na * size_of_datatype_real)
    call check_memcpy_GPU_f("tridiag: e_vec_dev", 633,  successGPU)

    successGPU = gpu_memset(tau_dev, 0, na * size_of_datatype)
    call check_memcpy_GPU_f("tridiag: tau_dev", 636,  successGPU)
  endif

  n_stored_vecs = 0

  l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a_mat
  l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a_mat

  if (my_prow == prow(na, nblk, np_rows) .and. my_pcol == pcol(na, nblk, np_cols)) then
      d_vec(na) = real(a_mat(l_rows,l_cols), kind=rk)
  endif

  if (useGPU) then
    ! allocate memory for matrix A on the device and than copy the matrix

    num = matrixRows * matrixCols * size_of_datatype

    successGPU = gpu_malloc(a_dev, num)
    call check_alloc_GPU_f("tridiag: a_dev", 659,  successGPU)


    successGPU = gpu_memcpy(a_dev, int(loc(a_mat(1,1)),kind=c_intptr_t), &
                              num, gpuMemcpyHostToDevice)
    call check_memcpy_GPU_f("tridiag: a_dev", 678,  successGPU)

  endif ! useGPU

  ! main cycle of tridiagonalization
  ! in each step, 1 Householder Vector is calculated
  do istep = na, nblockEnd ,-1


    ! Calculate number of local rows and columns of the still remaining matrix
    ! on the local processor
    l_rows = local_index(istep-1, my_prow, np_rows, nblk, -1)
    l_cols = local_index(istep-1, my_pcol, np_cols, nblk, -1)

    ! Calculate Vector for Householder transformation on all procs
    ! owning column istep

    if (useGPU) then 
      ! copy l_cols + 1 column of a_dev to v_row_dev
      isOurProcessRow      = (my_prow == prow(istep-1, nblk, np_rows))
      isOurProcessCol      = (my_pcol == pcol(istep-1, nblk, np_cols))
      isOurProcessCol_prev = (my_pcol == pcol(istep  , nblk, np_cols)) ! isOurProcessCol from the previous step
      call gpu_copy_and_set_zeros_float_complex(v_row_dev, a_dev, l_rows, l_cols, matrixRows, istep, &
                                            aux1_dev, vav_dev, d_vec_dev, &
                                            isOurProcessRow, isOurProcessCol, isOurProcessCol_prev, &
                                            isSkewsymmetric, useCCL, wantDebug, my_stream)
      if (gpu_vendor() /= SYCL_GPU) successGPU = gpu_DeviceSynchronize()
    endif ! useGPU

    if (my_pcol == pcol(istep, nblk, np_cols)) then

      ! Get Vector to be transformed; distribute last element and norm of
      ! remaining elements to all procs in current column

!             ! copy l_cols + 1 column of A to v_row
      ! if (useGPU) then

! #ifdef WITH_NVTX
!         call nvtxRangePush("memcpy new D-D a_dev(:,l_cols+1)->v_row_dev")
! #endif
!         ! TODO_23_11:  create a dev-dev copy kernel or merge it to another kernel
!         offset_dev = l_cols * matrixRows * size_of_datatype
!         successGPU = gpu_memcpy(v_row_dev, a_dev + offset_dev, (l_rows)* size_of_datatype, gpuMemcpyDeviceToDevice)
!         call check_memcpy_GPU_f("tridiag a_dev 1", 731,  successGPU)

! #ifdef WITH_NVTX
!         call nvtxRangePop()
! #endif

!       else ! useGPU
!         v_row(1:l_rows) = a_mat(1:l_rows,l_cols+1)
!       endif ! useGPU

      ! copy l_cols + 1 column of A to v_row
      if (.not. useGPU) then
        v_row(1:l_rows) = a_mat(1:l_rows,l_cols+1)
      endif ! useGPU

      if (n_stored_vecs > 0 .and. l_rows > 0) then
        if (useGPU) then
          if (wantDebug) call obj%timer%start("gpublas gemv skinny with copying")
          ! v_row_dev = vu_stored_rows_dev * uv_stored_cols_dev(l_cols+1,1:2*n_stored_vecs) + v_row_dev
          gpuHandle = obj%gpu_setup%gpublasHandleArray(0)
          call gpublas_CGEMV('N', l_rows, 2*n_stored_vecs,  &
                                    ONE, vu_stored_rows_dev, max_local_rows,  &
                                    aux_complex_dev, 1,   &
                                    ONE, v_row_dev, 1, gpuHandle)

          if (wantDebug .and. gpu_vendor() /= SYCL_GPU) successGPU = gpu_DeviceSynchronize()
          if (wantDebug) call obj%timer%stop("gpublas gemv skinny with copying")

        else ! useGPU
          aux(1:2*n_stored_vecs) = conjg(uv_stored_cols(l_cols+1,1:2*n_stored_vecs))

          if (wantDebug) call obj%timer%start("blas")
          ! v_row = vu_stored_rows * uv_stored_cols(l_cols+1,1:2*n_stored_vecs) + v_row
          call CGEMV('N',   &
                            int(l_rows,kind=BLAS_KIND), int(2*n_stored_vecs,kind=BLAS_KIND), &
                            ONE, vu_stored_rows, int(ubound(vu_stored_rows,dim=1),kind=BLAS_KIND), &
                            aux, 1_BLAS_KIND,  &
                            ONE, v_row, 1_BLAS_KIND)
          if (wantDebug) call obj%timer%stop("blas")
        endif ! useGPU
      endif ! (n_stored_vecs > 0 .and. l_rows > 0)

      if (useGPU) then
        if (my_prow == prow(istep-1, nblk, np_rows)) then
          isOurProcessRowInt = 1
        else
          isOurProcessRowInt = 0
        end if

        my_stream = obj%gpu_setup%my_stream
        call gpu_dot_product_and_assign_float_complex(v_row_dev, l_rows, isOurProcessRowInt, aux1_dev, wantDebug, my_stream)
        if (.not. useCCL) then


          successGPU = gpu_memcpy(int(loc(aux1),kind=c_intptr_t), aux1_dev, 2*size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("tridiag: aux1_dev -> aux1", 819,  successGPU)


        endif ! .not. useCCL 
      endif ! useGPU  


      if (.not. useGPU) then
        if (my_prow == prow(istep-1, nblk, np_rows)) then
          aux1(1) = dot_product(v_row(1:l_rows-1),v_row(1:l_rows-1)) ! = "q"
          aux1(2) = v_row(l_rows) ! = "a_11" (or rather a_nn)
        else
          aux1(1) = dot_product(v_row(1:l_rows),v_row(1:l_rows))
          aux1(2) = 0.
        endif
      endif ! .not. useGPU 


      if (useCCL) then
      else ! useCCL
        vnorm2 = real(aux1(1),kind=rk)
        vrl    = aux1(2)
      endif ! useCCL


      ! Householder transformation
      if (useCCL) then
      else ! useCCL
        call hh_transform_complex_&
                &single &
                (obj, vrl, vnorm2, xf, tau(istep), wantDebug)
      endif ! useCCL
      
      ! vrl is newly computed off-diagonal element of the final tridiagonal matrix
      if (my_prow == prow(istep-1, nblk, np_rows)) then
        if (.not. useCCL) then
          e_vec(istep-1) = real(vrl,kind=rk)
        endif ! .not. useCCL
      endif

      if (.not. useGPU) then
        ! Scale v_row and store Householder Vector for back transformation
        v_row(1:l_rows) = v_row(1:l_rows) * xf

        if (my_prow == prow(istep-1, nblk, np_rows)) then
          v_row(l_rows) = 1.
        endif

        ! store Householder Vector for back transformation
        ! update a_mat
        a_mat(1:l_rows,l_cols+1) = v_row(1:l_rows)
      endif ! .not. useGPU 

      if (.not. useCCL) then
        ! add tau after the end of actuall v_row, to be broadcasted with it
        v_row(l_rows+1) = tau(istep)
      endif

      
      if (useGPU) then
        if (useCCL) then
          xf_host_or_dev = xf_dev
        else 
          xf_host_or_dev = int(loc(xf),kind=c_intptr_t)
        endif       
        
        isOurProcessRow = (my_prow == prow(istep-1, nblk, np_rows))
        call gpu_set_e_vec_scale_set_one_store_v_row_float_complex(e_vec_dev, vrl_dev, a_dev, v_row_dev, tau_dev, xf_host_or_dev, & 
                                                  l_rows, l_cols, matrixRows, istep, isOurProcessRow, useCCL, wantDebug, my_stream)
      endif ! useGPU  


      if (useGPU .and. .not. useCCL) then      
        !v_row_dev -> v_row

        successGPU = gpu_memcpy(int(loc(v_row),kind=c_intptr_t), v_row_dev, l_rows*size_of_datatype, gpuMemcpyDeviceToHost)
        call check_memcpy_GPU_f("tridiag: v_row_dev -> v_row", 995,  successGPU)

      endif ! useGPU .and. .not. useCCL

    endif !(my_pcol == pcol(istep, nblk, np_cols))



    !recover tau, which has been broadcasted together with v_row
    if (.not. useCCL) then
      tau(istep) =  v_row(l_rows+1)
    endif

    ! Transpose Householder Vector v_row -> v_col
    if (useCCL) then
    else ! useCCL
      call elpa_transpose_vectors_&
          &complex&
          &_&
          &single &
                (obj, v_row, ubound(v_row,dim=1), mpi_comm_rows, v_col, ubound(v_col,dim=1), mpi_comm_cols, &
                1, istep-1, 1, nblk, max_threads, .true., success)
      if (wantDebug .and. gpu_vendor() /= SYCL_GPU) successGPU = gpu_DeviceSynchronize()
      if (.not.(success)) then
        write(error_unit,*) "Error in elpa_transpose_vectors. Aborting!"
        return
      endif
    endif ! useCCL

    ! Calculate u = (A + VU**T + UV**T)*v // Dongarra 1987: "y = (A - UV**T - VU**T)*u"

    ! For cache efficiency, we use only the upper half of the matrix tiles for this,
    ! thus the result is partly in u_col(:) and partly in u_row(:)

    if (.not. useGPU) then
      u_col(1:l_cols) = 0
      u_row(1:l_rows) = 0
    endif

    if (useGPU .and. useCCL) then
      successGPU = gpu_memset(u_col_dev, 0, l_cols * size_of_datatype) ! TODO_23_11: omit this, but change gpublas_gemm to u_col_dev=a_dev^T*v_row_dev+0*u_col_dev?
      call check_memcpy_GPU_f("tridiag: u_col_dev", 1101,  successGPU)
    endif

    if (l_rows>0 .and. l_cols>0) then
      if (useGPU) then

        if (.not. mat_vec_as_one_block) then

          successGPU = gpu_memcpy(v_col_dev, int(loc(v_col(1)),kind=c_intptr_t), &
                        l_cols * size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("tridiag: v_col_dev", 1146,  successGPU)


        endif ! .not. mat_vec_as_one_block

        if (.not. useCCL) then

          successGPU = gpu_memcpy(v_row_dev, int(loc(v_row(1)),kind=c_intptr_t), &
                                    l_rows * size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("tridiag: v_row_dev", 1170,  successGPU)


        endif ! .not. useCCL
      endif ! useGPU


      if (.not. useGPU) then
        do i=0, (istep-2)/tile_size ! iteration over tiles
          l_col_beg = i*l_cols_per_tile+1
          l_col_end = min(l_cols,(i+1)*l_cols_per_tile)
          if (l_col_end < l_col_beg) cycle
          do j = 0, i
            l_row_beg = j*l_rows_per_tile+1
            l_row_end = min(l_rows,(j+1)*l_rows_per_tile)
            if (l_row_end < l_row_beg) cycle

            ! multiplication by blocks is efficient only for CPU
            ! for GPU we introduced 2 other ways, either by stripes (more simmilar to the original
            ! CPU implementation) or by one large matrix Vector multiply
            if (wantDebug) call obj%timer%start("blas")
            ! u_col = a_mat*v_row + u_col(=0)
            call CGEMV('C',  &
                        int(l_row_end-l_row_beg+1,kind=BLAS_KIND), int(l_col_end-l_col_beg+1,kind=BLAS_KIND), &
                        ONE, a_mat(l_row_beg, l_col_beg), int(matrixRows,kind=BLAS_KIND),         &
                        v_row(l_row_beg:max_local_rows+1), 1_BLAS_KIND,                           &
                        ONE, u_col(l_col_beg:max_local_cols), 1_BLAS_KIND)

            if (i/=j) then
              if (isSkewsymmetric) then
                call CGEMV('N',int(l_row_end-l_row_beg+1,kind=BLAS_KIND), int(l_col_end-l_col_beg+1,kind=BLAS_KIND), &
                                    -ONE, a_mat(l_row_beg,l_col_beg), int(matrixRows,kind=BLAS_KIND),               &
                                    v_col(l_col_beg:max_local_cols), 1_BLAS_KIND, ONE, u_row(l_row_beg:max_local_rows), &
                                    1_BLAS_KIND)

              else
                ! u_row = a_mat*v_col + u_row(=0)
                call CGEMV('N',int(l_row_end-l_row_beg+1,kind=BLAS_KIND), int(l_col_end-l_col_beg+1,kind=BLAS_KIND),  &
                                    ONE, a_mat(l_row_beg,l_col_beg), int(matrixRows,kind=BLAS_KIND),               &
                                    v_col(l_col_beg:max_local_cols), 1_BLAS_KIND, ONE, u_row(l_row_beg:max_local_rows), &
                                    1_BLAS_KIND)
              endif
            endif
            if (wantDebug) call obj%timer%stop("blas")

          enddo  ! j=0,i
        enddo  ! i=0,(istep-2)/tile_size
      endif ! .not. useGPU

      if (useGPU) then
        if (mat_vec_as_one_block) then
          ! Unlike for CPU, we (for each MPI thread) do just one large mat-vec multiplication
          ! this requires altering of the algorithm when later explicitly updating the matrix
          ! after max_stored_uv is reached : we need to update all tiles, not only those above diagonal
          if (wantDebug) call obj%timer%start("gpublas")
          
          gpuHandle = obj%gpu_setup%gpublasHandleArray(0) ! 0.3-0.5 us

              
          ! u_col_dev = a_dev^T*v_row_dev
          call gpublas_CGEMV('C', l_rows,l_cols,  &
                                    ONE, a_dev, matrixRows,                   &
                                    v_row_dev , 1,                          &
                                    ZERO, u_col_dev, 1, gpuHandle)
              
          if (wantDebug .and. gpu_vendor() /= SYCL_GPU) successGPU = gpu_DeviceSynchronize()
              
       ! todo: try with non transposed!!!
!                 if(i/=j) then
!                   call gpublas_CGEMV('N', l_row_end-l_row_beg+1,l_col_end-l_col_beg+1,  &
!                                             ONE, a_dev + offset_dev, matrixRows,                        &
!                                             v_col_dev + (l_col_beg - 1) *                      &
!                                             size_of_datatype, 1,                          &
!                                             ONE, u_row_dev + (l_row_beg - 1) *                 &
!                                             size_of_datatype, 1)
!                 endif
          if (wantDebug) call obj%timer%stop("gpublas")
        else  ! mat_vec_as_one_block
          !perform multiplication by stripes - it is faster than by blocks, since we call cublas with
          !larger matrices. In general, however, this algorithm is very simmilar to the one with CPU
          do i=0,(istep-2)/tile_size
            l_col_beg = i*l_cols_per_tile+1
            l_col_end = min(l_cols,(i+1)*l_cols_per_tile)
            if (l_col_end<l_col_beg) cycle

            l_row_beg = 1
            l_row_end = min(l_rows,(i+1)*l_rows_per_tile)
                  
            offset_dev = ((l_row_beg-1) + (l_col_beg - 1) * matrixRows) * size_of_datatype
  
            gpuHandle = obj%gpu_setup%gpublasHandleArray(0)
            call gpublas_CGEMV('C', &
                          l_row_end-l_row_beg+1, l_col_end-l_col_beg+1, &
                          ONE, a_dev + offset_dev, matrixRows,  &
                          v_row_dev + (l_row_beg - 1) * size_of_datatype, 1,  &
                          ONE, u_col_dev + (l_col_beg - 1) * size_of_datatype, 1, gpuHandle)
          enddo !i=0,(istep-2)/tile_size

          do i=0,(istep-2)/tile_size
            l_col_beg = i*l_cols_per_tile+1
            l_col_end = min(l_cols,(i+1)*l_cols_per_tile)
            if (l_col_end<l_col_beg) cycle

            l_row_beg = 1
            l_row_end = min(l_rows,i*l_rows_per_tile)
              
            offset_dev = ((l_row_beg-1) + (l_col_beg - 1) * matrixRows) * size_of_datatype

            if (isSkewsymmetric) then
                gpuHandle = obj%gpu_setup%gpublasHandleArray(0)
                call gpublas_CGEMV('N', l_row_end-l_row_beg+1, l_col_end-l_col_beg+1, &
                            -ONE, a_dev + offset_dev, matrixRows, &
                            v_col_dev + (l_col_beg - 1) * size_of_datatype,1, &
                            ONE, u_row_dev + (l_row_beg - 1) * size_of_datatype, 1, gpuHandle)
            else
                gpuHandle = obj%gpu_setup%gpublasHandleArray(0)
                call gpublas_CGEMV('N', l_row_end-l_row_beg+1, l_col_end-l_col_beg+1, &
                            ONE, a_dev + offset_dev, matrixRows, &
                            v_col_dev + (l_col_beg - 1) * size_of_datatype,1, &
                            ONE, u_row_dev + (l_row_beg - 1) * size_of_datatype, 1, gpuHandle)
            endif
          enddo ! i=0,(istep-2)/tile_size
        end if ! mat_vec_as_one_block / per stripes


        if (.not. mat_vec_as_one_block) then

          successGPU = gpu_memcpy(int(loc(u_row(1)),kind=c_intptr_t), u_row_dev, l_rows*size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("tridiag: u_row_dev 1", 1378,  successGPU)


        endif ! .not. mat_vec_as_one_block
      endif ! useGPU


      ! second calculate (VU**T + UV**T)*v part of (A + VU**T + UV**T)*v
      if (n_stored_vecs > 0) then
        if (.not. useGPU) then
          if (wantDebug) call obj%timer%start("blas")

      
          ! aux = vu_stored_rows^T*v_row
          call CGEMV('C',     &
                              int(l_rows,kind=BLAS_KIND), int(2*n_stored_vecs,kind=BLAS_KIND),   &
                              ONE, vu_stored_rows, int(ubound(vu_stored_rows,dim=1),kind=BLAS_KIND),   &
                              v_row,  1_BLAS_KIND, ZERO, aux, 1_BLAS_KIND)
              
          ! u_col = uv_stored_cols*aux + u_col
          call CGEMV('N', int(l_cols,kind=BLAS_KIND), int(2*n_stored_vecs,kind=BLAS_KIND),   &
                              ONE, uv_stored_cols, int(ubound(uv_stored_cols,dim=1),kind=BLAS_KIND),   &
                              aux, 1_BLAS_KIND, ONE, u_col,  1_BLAS_KIND)
          if (wantDebug) call obj%timer%stop("blas")
        else ! .not. useGPU
          if (wantDebug) call obj%timer%start("gpublas gemv x2 skinny")


          ! aux_dev = vu_stored_rows_dev^T*v_row_dev
          gpuHandle = obj%gpu_setup%gpublasHandleArray(0)
          call gpublas_CGEMV('C', l_rows, 2*n_stored_vecs, &
                                      ONE, vu_stored_rows_dev, max_local_rows,   &
                                      v_row_dev,  1, ZERO, aux_dev, 1, gpuHandle)
              
          ! u_col_dev = uv_stored_cols_dev*aux_dev + u_col_dev
          call gpublas_CGEMV('N', l_cols, 2*n_stored_vecs, &
                                      ONE, uv_stored_cols_dev, max_local_cols,   &
                                      aux_dev, 1, ONE, u_col_dev, 1, gpuHandle)

          if (wantDebug .and. gpu_vendor() /= SYCL_GPU) successGPU = gpu_DeviceSynchronize()

          if (wantDebug) call obj%timer%stop("gpublas gemv x2 skinny")  
        endif ! .not. useGPU
      endif ! n_stored_vecs > 0

    endif  ! (l_rows>0 .and. l_cols>0)

    if (useGPU .and. l_cols>0 .and. (.not. useCCL)) then

      successGPU = gpu_memcpy(int(loc(u_col(1)),kind=c_intptr_t), u_col_dev, l_cols*size_of_datatype, gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("tridiag: u_col_dev 1", 1462,  successGPU)

    endif ! useGPU


    if (useGPU .and. .not. useCCL) then 
        if (l_rows==0) then
        u_col(1:l_cols) = 0
      endif
    endif ! useGPU

    ! Sum up all u_row(:) parts along rows and add them to the u_col(:) parts
    ! on the processors containing the diagonal
    ! This is only necessary if u_row has been calculated, i.e. if the
    ! global tile size is smaller than the global remaining matrix
    
    ! in GPU case, u_row_dev=0 up to now, so elpa_reduce_add_vectors is skipped

    if (tile_size < istep-1 .and. .not. useGPU) then
      call elpa_reduce_add_vectors_&
            &complex&
            &_&
            &single &
            (obj, u_row, ubound(u_row,dim=1), mpi_comm_rows, u_col, ubound(u_col,dim=1), &
            mpi_comm_cols, istep-1, 1, nblk, max_threads)
    endif ! (tile_size < istep-1 .and. .not. useGPU)

    ! Sum up all the u_col(:) parts, transpose u_col -> u_row

    if (l_cols > 0) then
    endif ! (l_cols > 0)
    
    ! Transpose Householder Vector u_col -> u_row
    if (useCCL) then
    else ! useCCL
      if (isSkewsymmetric) then
        call elpa_transpose_vectors_ss_&
        &complex&
        &_&
        &single &
        (obj, u_col, ubound(u_col,dim=1), mpi_comm_cols, u_row, ubound(u_row,dim=1), &
          mpi_comm_rows, 1, istep-1, 1, nblk, max_threads, .false., success)
        if (.not.(success)) then
          write(error_unit,*) "Error in elpa_transpose_vectors_ss. Aborting!"
          return
        endif
      else ! isSkewsymmetric
        call elpa_transpose_vectors_&
        &complex&
        &_&
        &single &
        (obj, u_col, ubound(u_col,dim=1), mpi_comm_cols, u_row, ubound(u_row,dim=1), &
        mpi_comm_rows, 1, istep-1, 1, nblk, max_threads, .false., success)
        if (.not.(success)) then
          write(error_unit,*) "Error in elpa_transpose_vectors. Aborting!"
          return
        endif
      endif ! isSkewsymmetric
    endif ! useCCL


    if (useGPU .and. .not. useCCL) then

      successGPU = gpu_memcpy(v_col_dev, int(loc(v_col(1)),kind=c_intptr_t), l_cols*size_of_datatype, gpuMemcpyHostToDevice)
      call check_memcpy_GPU_f("tridiag: v_col_dev", 1607,  successGPU)


      successGPU = gpu_memcpy(u_col_dev, int(loc(u_col(1)),kind=c_intptr_t), l_cols*size_of_datatype, gpuMemcpyHostToDevice)
      call check_memcpy_GPU_f("tridiag: u_col_dev", 1622,  successGPU)

    endif ! (useGPU .and. .not. useCCL)



    if (useGPU .and. .not. useCCL) then

      successGPU = gpu_memcpy(u_row_dev, int(loc(u_row(1)),kind=c_intptr_t), l_rows*size_of_datatype, gpuMemcpyHostToDevice)
      call check_memcpy_GPU_f("tridiag: u_row_dev", 1657,  successGPU)


      successGPU = gpu_memcpy(v_row_dev, int(loc(v_row(1)),kind=c_intptr_t), l_rows*size_of_datatype, gpuMemcpyHostToDevice)
      call check_memcpy_GPU_f("tridiag: v_row_dev", 1671,  successGPU)

    endif ! (useGPU .and. .not. useCCL)


    ! calculate u**T * v (same as v**T * (A + VU**T + UV**T) * v )
    if (.not. useGPU .or. (useGPU .and. .not. useCCL)) then
      vav = 0 ! x=0
      if (l_cols>0) vav = dot_product(v_col(1:l_cols), u_col(1:l_cols))
    endif



    ! store u and v in the matrices U and V
    ! these matrices are stored combined in one here

   if (.not. useCCL) then
      conjg_tau = conjg(tau(istep))
    endif ! .not. useCCL
    
    if (.not. useGPU) then
      if (l_rows > 0) then
        ! update vu_stored_rows
        vu_stored_rows(1:l_rows,2*n_stored_vecs+1) = conjg_tau*v_row(1:l_rows)
        vu_stored_rows(1:l_rows,2*n_stored_vecs+2) = 0.5*conjg_tau*vav*v_row(1:l_rows) - u_row(1:l_rows)
      endif
      if (l_cols > 0) then
        ! update uv_stored_cols
        uv_stored_cols(1:l_cols,2*n_stored_vecs+1) = 0.5*conjg_tau*vav*v_col(1:l_cols) - u_col(1:l_cols)
        uv_stored_cols(1:l_cols,2*n_stored_vecs+2) = conjg_tau*v_col(1:l_cols)
      endif
    endif ! .not. useGPU


    if (useGPU) then

      ! kernel: update vu_stored_rows_dev, uv_stored_cols
      ! then cpu's "store u,v in U,V" can be deleted. But we should take care of dot_prod below, where vu_stored_rows and uv_stored_cols are used
      if (useCCL) then
        vav_host_or_dev = vav_dev
        !tau_istep_host_or_dev = tau_dev + (istep-1)*size_of_datatype
        tau_istep_host_or_dev = v_row_dev + (l_rows+1-1)*size_of_datatype
      else ! useCCL
        vav_host_or_dev = int(loc(vav),kind=c_intptr_t)
        tau_istep_host_or_dev = int(loc(tau(istep)), kind=c_intptr_t)
      endif ! useCCL
      call gpu_store_u_v_in_uv_vu_float_complex(vu_stored_rows_dev, uv_stored_cols_dev, v_row_dev, u_row_dev, &
                                          v_col_dev, u_col_dev, tau_dev, aux_complex_dev, &
                                          vav_host_or_dev, tau_istep_host_or_dev, &
                                          l_rows, l_cols, n_stored_vecs,  max_local_rows, max_local_cols, istep, &
                                          useCCL, wantDebug, my_stream)
    endif ! useGPU


    ! We have calculated another Householder Vector, number of implicitly stored increased
    n_stored_vecs = n_stored_vecs+1

    ! If the limit of max_stored_uv is reached, calculate A + VU**T + UV**T
    if (n_stored_vecs == max_stored_uv .or. istep == 3) then        
      
      if (.not. useGPU .OR. .not. mat_vec_as_one_block) then
        do i = 0, (istep-2)/tile_size
          ! go over tiles above (or on) the diagonal
          l_col_beg = i*l_cols_per_tile+1
          l_col_end = min(l_cols,(i+1)*l_cols_per_tile)
          l_row_beg = 1
          l_row_end = min(l_rows,(i+1)*l_rows_per_tile)
          if (l_col_end<l_col_beg .or. l_row_end<l_row_beg) then
            cycle
          endif

          if (useGPU) then
            if (.not. mat_vec_as_one_block) then
              ! if using mat-vec multiply by stripes, it is enough to update tiles above (or on) the diagonal only
              ! we than use the same calls as for CPU version
              if (wantDebug) call obj%timer%start("gpublas_gemm")
              ! update a_dev
              ! a_dev = vu_stored_rows_dev*uv_stored_cols_dev + a_dev
              gpuHandle = obj%gpu_setup%gpublasHandleArray(0)
              call gpublas_CGEMM('N', 'C',     &
                                      l_row_end-l_row_beg+1, l_col_end-l_col_beg+1, 2*n_stored_vecs,                      &
                                      ONE, vu_stored_rows_dev + (l_row_beg - 1) *                                         &
                                      size_of_datatype,  &
                                      max_local_rows, uv_stored_cols_dev + (l_col_beg - 1) *                              &
                                      size_of_datatype,  &
                                      max_local_cols, ONE, a_dev + ((l_row_beg - 1) + (l_col_beg - 1) * matrixRows) *     &
                                      size_of_datatype , matrixRows, gpuHandle)
              if (wantDebug .and. gpu_vendor() /= SYCL_GPU) successGPU = gpu_DeviceSynchronize()
              if (wantDebug) call obj%timer%stop("gpublas_gemm")
            endif ! .not. mat_vec_as_one_block
          else ! useGPU
            if (wantDebug) call obj%timer%start("blas_gemm")
            ! update a_mat
            ! a_mat = vu_stored_rows*uv_stored_cols + a_mat
            call CGEMM('N', 'C',                &
                                int(l_row_end-l_row_beg+1,kind=BLAS_KIND), int(l_col_end-l_col_beg+1,kind=BLAS_KIND), &
                                int(2*n_stored_vecs,kind=BLAS_KIND),    &
                                ONE, vu_stored_rows(l_row_beg:max_local_rows,1:2*max_stored_uv), &
                                int(ubound(vu_stored_rows,dim=1),kind=BLAS_KIND),   &
                                uv_stored_cols(l_col_beg,1), &
                                int(ubound(uv_stored_cols,dim=1),kind=BLAS_KIND),        &
                                ONE, a_mat(l_row_beg,l_col_beg), int(matrixRows,kind=BLAS_KIND))
            if (wantDebug) call obj%timer%stop("blas_gemm")
          endif !useGPU
        enddo ! i = 0, (istep-2)/tile_size

      else !.not. useGPU or .not. mat_vec_as_one_block (i.e. useGPU and mat_vec_as_one_block)

        !update whole (remaining) part of matrix, including tiles below diagonal
        !we can do that in one large cublas call
        if (wantDebug) call obj%timer%start("gpublas_gemm")


        gpuHandle = obj%gpu_setup%gpublasHandleArray(0)
        ! update a_dev
        ! a_dev = vu_stored_rows_dev*uv_stored_cols_dev + a_dev
        call gpublas_CGEMM('N', 'C', l_rows, l_cols, 2*n_stored_vecs,   &
                                  ONE, vu_stored_rows_dev, max_local_rows, &
                                  uv_stored_cols_dev, max_local_cols,  &
                                  ONE, a_dev, matrixRows, gpuHandle)

        if (wantDebug .and. gpu_vendor() /= SYCL_GPU) successGPU = gpu_DeviceSynchronize()
        if (wantDebug) call obj%timer%stop("gpublas_gemm")
        
        ! copy real(a_dev(l_rows,l_cols)) -> d_vec_dev(istep-1) for correct initial value of d_vec_dev before atomicAdd
        if ((.not. isSkewsymmetric) .and. &
            (my_prow == prow(istep-1, nblk, np_rows) .and. my_pcol == pcol(istep-1, nblk, np_cols))) then
          offset_dev = ((l_rows-1) + (l_cols-1)*matrixRows) * size_of_datatype
          successGPU = gpu_memcpy(d_vec_dev + (istep-2)*size_of_datatype_real, &
                                  a_dev + offset_dev, 1*size_of_datatype_real, gpuMemcpyDeviceToDevice)
          call check_memcpy_GPU_f("tridiag a_dev->d_vec_dev", 1863,  successGPU)
        endif

      endif !.not. useGPU or .not. mat_vec_as_one_block

      n_stored_vecs = 0
    endif ! (n_stored_vecs == max_stored_uv .or. istep == 3)

    if (my_prow == prow(istep-1, nblk, np_rows) .and. my_pcol == pcol(istep-1, nblk, np_cols)) then

      if (useGPU) then
        my_stream = obj%gpu_setup%my_stream
        call gpu_update_matrix_element_add_float_complex(vu_stored_rows_dev, uv_stored_cols_dev, a_dev, d_vec_dev, &
                                            l_rows, l_cols, matrixRows, max_local_rows, max_local_cols, istep, n_stored_vecs, &
                                            isSkewsymmetric, wantDebug, my_stream)

      else ! useGPU
        if (n_stored_vecs > 0) then
          ! update a_mat (only one elememt!)
          dot_prod = dot_product(vu_stored_rows(l_rows,1:2*n_stored_vecs), uv_stored_cols(l_cols,1:2*n_stored_vecs))
          a_mat(l_rows,l_cols) = a_mat(l_rows,l_cols) + dot_prod
        endif
        d_vec(istep-1) = real(a_mat(l_rows,l_cols),kind=rk)
      endif ! useGPU

    endif ! (my_prow == prow(istep-1, nblk, np_rows) .and. my_pcol == pcol(istep-1, nblk, np_cols))


  enddo ! main cycle over istep=na,3,-1


  if (useGPU) then
    ! copy a_dev -> a_mat for backtransformation
    num = matrixRows * matrixCols * size_of_datatype
    successGPU = gpu_memcpy(int(loc(a_mat(1,1)),kind=c_intptr_t), a_dev, num, gpuMemcpyDeviceToHost)
    call check_memcpy_GPU_f("tridiag: a_dev", 1921,  successGPU)

  endif ! useGPU


  ! Store e_vec(1) and d_vec(1)

  if (my_pcol==pcol(2, nblk, np_cols)) then
    if (my_prow==prow(1, nblk, np_rows)) then
      ! We use last l_cols value of loop above
      if (useGPU) then
        successGPU = gpu_memcpy(int(loc(aux3(1)),kind=c_intptr_t), a_dev + (matrixRows * (l_cols - 1)) * size_of_datatype, &
                                1 * size_of_datatype, gpuMemcpyDeviceToHost)
        call check_memcpy_GPU_f("tridiag: a_dev 5", 1944,  successGPU)
        vrl = aux3(1)
      else !useGPU
        vrl = a_mat(1,l_cols)
      endif !useGPU

      call hh_transform_complex_&
            &single &
            (obj, vrl, 0.0_rk, xf, tau(2), wantDebug)
      e_vec(1) = real(vrl,kind=rk)
      a_mat(1,l_cols) = 1. ! for consistency only
    endif ! (my_prow==prow(1, nblk, np_rows))


  endif ! (my_pcol==pcol(2, nblk, np_cols))


  if (my_prow == prow(1, nblk, np_rows) .and. my_pcol == pcol(1, nblk, np_cols))  then
    if (useGPU) then
      successGPU = gpu_memcpy(int(loc(aux3(1)),kind=c_intptr_t), a_dev, &
                             1 * size_of_datatype, gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("tridiag: a_dev 6", 2002,  successGPU)
      d_vec(1) = REAL(aux3(1))
    else !useGPU
      d_vec(1) = REAL(a_mat(1,1))
    endif !useGPU
  endif ! (my_prow == prow(1, nblk, np_rows) .and. my_pcol == pcol(1, nblk, np_cols))


  if (useGPU) then
    offset_dev = 1 * size_of_datatype_real
    ! first and last elements of d_vec are treated separately
    successGPU = gpu_memcpy(int(loc(d_vec(2)),kind=c_intptr_t), &
                            d_vec_dev + offset_dev, (na-2) * size_of_datatype_real, gpuMemcpyDeviceToHost)
    call check_memcpy_GPU_f("tridiag: d_vec", 2062,  successGPU)

    if (useCCL) then
      ! e_vec(1) is treated separately
      offset_dev = 1 * size_of_datatype_real
      successGPU = gpu_memcpy(int(loc(e_vec(2)),kind=c_intptr_t), &
                              e_vec_dev + offset_dev, (na-1) * size_of_datatype_real, gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("tridiag: e_vec", 2069,  successGPU)

      ! tau(2) is treated separately, tau(1) is not used
      offset_dev = 2 * size_of_datatype
      successGPU = gpu_memcpy(int(loc(tau(3)),kind=c_intptr_t), &
                              tau_dev + offset_dev, (na-2) * size_of_datatype, gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("tridiag: tau", 2075,  successGPU)
    endif

    ! todo: should we leave a_mat on the device for further use?
    successGPU = gpu_free(a_dev)
    call check_dealloc_GPU_f("tridiag: a_dev 9", 2080,  successGPU)

    successGPU = gpu_free(v_row_dev)
    call check_dealloc_GPU_f("tridiag: v_row_dev", 2083,  successGPU)

    successGPU = gpu_free(u_row_dev)
    call check_dealloc_GPU_f("tridiag: (u_row_dev", 2086,  successGPU)

    successGPU = gpu_free(v_col_dev)
    call check_dealloc_GPU_f("tridiag: v_col_dev", 2089,  successGPU)

    successGPU = gpu_free(u_col_dev)
    call check_dealloc_GPU_f("tridiag: u_col_dev ", 2092,  successGPU)

    successGPU = gpu_free(vu_stored_rows_dev)
    call check_dealloc_GPU_f("tridiag: vu_stored_rows_dev ", 2095,  successGPU)

    successGPU = gpu_free(uv_stored_cols_dev)
    call check_dealloc_GPU_f("tridiag:uv_stored_cols_dev ", 2098,  successGPU)

    successGPU = gpu_free(d_vec_dev)
    call check_dealloc_GPU_f("tridiag: d_vec_dev", 2101,  successGPU)

    successGPU = gpu_free(e_vec_dev)
    call check_dealloc_GPU_f("tridiag: e_vec_dev", 2104,  successGPU)

    successGPU = gpu_free(tau_dev)
    call check_dealloc_GPU_f("tridiag: tau_dev", 2107,  successGPU)

    successGPU = gpu_free(aux_dev)
    call check_dealloc_GPU_f("tridiag: aux_dev", 2110,  successGPU)

    successGPU = gpu_free(aux1_dev)
    call check_dealloc_GPU_f("tridiag: aux1_dev", 2113,  successGPU)

    successGPU = gpu_free(aux_complex_dev)
    call check_dealloc_GPU_f("tridiag: aux_complex_dev", 2116,  successGPU)

    successGPU = gpu_free(vav_dev)
    call check_dealloc_GPU_f("tridiag: vav_dev", 2119,  successGPU)

    successGPU = gpu_free(dot_prod_dev)
    call check_dealloc_GPU_f("tridiag: dot_prod_dev", 2122,  successGPU)

    successGPU = gpu_free(xf_dev)
    call check_dealloc_GPU_f("tridiag: xf_dev", 2125,  successGPU)

  endif ! useGPU

  ! distribute the arrays d_vec and e_vec to all processors

  allocate(tmp_real(na), stat=istat, errmsg=errorMessage)
  call check_allocate_f("tridiag: tmp_real", 2138,  istat,  errorMessage)


  deallocate(tmp_real, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag: tmp_real", 2187,  istat,  errorMessage)

  if (useGPU) then

  else ! useGPU
    deallocate(v_row, v_col, u_row, u_col, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("tridiag: v_row, v_col, u_row, u_col", 2260,  istat,  errorMessage)
  endif ! useGPU

  deallocate(vu_stored_rows, uv_stored_cols, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag: vu_stored_rows, uv_stored_cols", 2264,  istat,  errorMessage)

  deallocate(aux, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag: aux", 2267,  istat,  errorMessage)



  call obj%timer%stop("tridiag_&
  &complex&
  &" // &
  "_single" // &
  gpuString )

end subroutine tridiag_&
&complex&
&_&
&single





!> \brief Transforms the eigenvectors of a tridiagonal matrix back
!>                     to the eigenvectors of the original matrix
!>                     (like Scalapack Routine PDORMTR)
!>
!  Parameters
!
!> \param na          Order of matrix a_mat, number of rows of matrix q_mat
!>
!> \param nqc         Number of columns of matrix q_mat
!>
!> \param a_mat(lda,matrixCols)  Matrix containing the Householder vectors (i.e. matrix a after tridiag_real)
!>                           Distribution is like in Scalapack.
!>
!> \param lda         Leading dimension of a_mat
!>
!> \param tau(na)     Factors of the Householder vectors
!>
!> \param q_mat           On input: Eigenvectors of tridiagonal matrix
!>                    On output: Transformed eigenvectors
!>                    Distribution is like in Scalapack.
!>
!> \param ldq         Leading dimension of q_mat
!>
!> \param nblk        blocksize of cyclic distribution, must be the same in both directions!
!>
!> \param matrixCols  local columns of matrix a_mat and q_mat
!>
!> \param mpi_comm_rows        MPI-Communicator for rows
!>
!> \param mpi_comm_cols        MPI-Communicator for columns
!>
!> \param useGPU      If true,  GPU version of the subroutine will be used
!>


subroutine trans_ev_&
&complex&
&_&
&single &
(obj, na, nqc, a_mat, lda, tau, q_mat, ldq, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols, useGPU, success)
  use, intrinsic :: iso_c_binding
  use precision
  use elpa_abstract_impl
  use elpa_blas_interfaces
  use elpa_gpu
  use elpa_gpu_util

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
  class(elpa_abstract_impl_t), intent(inout)    :: obj
  integer(kind=ik), intent(in)                  :: na, nqc, lda, ldq, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols
  complex(kind=rck), intent(in)           :: tau(na)

  complex(kind=rck), intent(inout)        :: a_mat(lda,*)
  complex(kind=rck), intent(inout)        :: q_mat(ldq,*)
  logical, intent(in)                           :: useGPU
  integer(kind=ik)                              :: max_stored_rows, max_stored_rows_fac

  integer(kind=ik)                              :: my_prow, my_pcol, np_rows, np_cols
  integer(kind=MPI_KIND)                        :: mpierr, my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI
  integer(kind=ik)                              :: totalblocks, max_blocks_row, max_blocks_col, max_local_rows, max_local_cols
  integer(kind=ik)                              :: l_cols, l_rows, l_colh, nstor
  integer(kind=ik)                              :: istep, n, nc, ic, ics, ice, nb, cur_pcol
  integer(kind=ik)                              :: hvn_ubnd, hvm_ubnd

  complex(kind=rck), allocatable          :: hvb(:), hvm(:,:)
  complex(kind=rck), pointer              :: tmp1(:), tmp2(:)
  complex(kind=rck), allocatable          :: h1(:), h2(:), tmp_debug(:)
  complex(kind=rck), pointer              :: tmat(:,:)
  complex(kind=rck), pointer              :: hvm1(:)
  type(c_ptr)                                   :: tmp1_host, tmp2_host
  type(c_ptr)                                   :: hvm1_host, tmat_host

  integer(kind=ik)                              :: istat
  character(200)                                :: errorMessage
  character(20)                                 :: gpuString

  integer(kind=c_intptr_t)                      :: num
  integer(kind=C_intptr_T)                      :: q_dev, tmp_dev, hvm_dev, tmat_dev

  integer(kind=ik)                              :: blockStep
  logical                                       :: successGPU
  integer(kind=c_intptr_t), parameter           :: size_of_datatype = size_of_&
                                                                      &single&
                                                                      &_&
                                                                      &complex
  integer(kind=ik)                              :: error
  integer(kind=MPI_KIND)                        :: bcast_request1, allreduce_request1, allreduce_request2
  logical                                       :: useNonBlockingCollectivesCols
  logical                                       :: useNonBlockingCollectivesRows
  integer(kind=c_int)                           :: non_blocking_collectives_rows, non_blocking_collectives_cols
  logical                                       :: success
  integer(kind=c_intptr_t)                      :: gpuHandle, my_stream


! NCCL requires streams
!#if defined (USE_CCL_TRANS_EV) && !defined(WITH_GPU_STREAMS)
!  successGPU = cuda_stream_create(obj%gpu_setup%my_stream)
!  if (.not.(successGPU)) then
!    print *,"Cannot create gpu stream handle"
!  endif
!  my_stream = obj%gpu_setup%my_stream
!#endif

  success = .true.

  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  call obj%timer%start("trans_ev_&
  &complex&
  &" // &
  &"_single" //&
  gpuString)

  call obj%get("nbc_row_elpa1_tridi_to_full", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    write(error_unit,*) "Problem setting option for non blocking collectives for rows in elpa1_tridi_to_full. Aborting..."
    call obj%timer%stop("trans_ev_&
    &complex&
    &" // &
    &"_single" //&
    gpuString)
    success = .false.
    return
  endif

  call obj%get("nbc_col_elpa1_tridi_to_full", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    write(error_unit,*) "Problem setting option for non blocking collectives for cols in elpa1_tridi_to_full. Aborting..."
    call obj%timer%stop("trans_ev_&
    &complex&
    &" // &
    &"_single" //&
    gpuString)
    success = .false.
    return
  endif

  if (non_blocking_collectives_rows .eq. 1) then
    useNonBlockingCollectivesRows = .true.
  else
    useNonBlockingCollectivesRows = .false.
  endif

  if (non_blocking_collectives_cols .eq. 1) then
    useNonBlockingCollectivesCols = .true.
  else
    useNonBlockingCollectivesCols = .false.
  endif

  my_prow = obj%mpi_setup%myRank_comm_rows
  my_pcol = obj%mpi_setup%myRank_comm_cols

  np_rows = obj%mpi_setup%nRanks_comm_rows
  np_cols = obj%mpi_setup%nRanks_comm_cols


  !call obj%timer%start("mpi_communication")
  !call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND) ,my_prowMPI, mpierr)
  !call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND) ,np_rowsMPI, mpierr)
  !call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND) ,my_pcolMPI, mpierr)
  !call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND) ,np_colsMPI, mpierr)

  !my_prow = int(my_prowMPI, kind=c_int)
  !np_rows = int(np_rowsMPI, kind=c_int)
  !my_pcol = int(my_pcolMPI, kind=c_int)
  !np_cols = int(np_colsMPI, kind=c_int)
  !call obj%timer%stop("mpi_communication")

  call obj%get("max_stored_rows",max_stored_rows_fac, error)

  totalblocks = (na-1)/nblk + 1
  max_blocks_row = (totalblocks-1)/np_rows + 1
  max_blocks_col = ((nqc-1)/nblk)/np_cols + 1  ! Columns of q_mat!

  max_local_rows = max_blocks_row*nblk
  max_local_cols = max_blocks_col*nblk

  max_stored_rows = (max_stored_rows_fac/nblk+1)*nblk
 
  if (.not.(useGPU)) then
    allocate(tmat(max_stored_rows,max_stored_rows), stat=istat, errmsg=errorMessage)
    call check_alloc("trans_ev_&
    &complex&
    &", "tmat", istat, errorMessage)

    allocate(tmp1(max_local_cols*max_stored_rows), stat=istat, errmsg=errorMessage)
    call check_alloc("trans_ev_&
    &complex&
    &", "tmp1", istat, errorMessage)

    allocate(tmp2(max_local_cols*max_stored_rows), stat=istat, errmsg=errorMessage)
    call check_alloc("trans_ev_&
    &complex&
    &", "tmp2", istat, errorMessage)
  endif

  allocate(h1(max_stored_rows*max_stored_rows), stat=istat, errmsg=errorMessage)
  call check_alloc("trans_ev_&
  &complex&
  &", "h1", istat, errorMessage)

  allocate(h2(max_stored_rows*max_stored_rows), stat=istat, errmsg=errorMessage)
  call check_alloc("trans_ev_&
  &complex&
  &", "h2", istat, errorMessage)

  allocate(hvb(max_local_rows*nblk), stat=istat, errmsg=errorMessage)
  call check_alloc("trans_ev_&
  &complex&
  &", "hvn", istat, errorMessage)

  allocate(hvm(max_local_rows,max_stored_rows), stat=istat, errmsg=errorMessage)
  call check_alloc("trans_ev_&
  &complex&
  &", "hvm", istat, errorMessage)

  hvm = 0   ! Must be set to 0 !!!
  hvb = 0   ! Safety only
  blockStep = nblk

  l_cols = local_index(nqc, my_pcol, np_cols, nblk, -1) ! Local columns of q_mat

  nstor = 0
  if (useGPU) then
    hvn_ubnd = 0
  endif

  ! In the complex case tau(2) /= 0
  if (my_prow == prow(1, nblk, np_rows)) then
    q_mat(1,1:l_cols) = q_mat(1,1:l_cols)*(ONE-tau(2))
  endif
 
  if (useGPU) then
    ! todo: this is used only for copying hmv to device.. it should be possible to go without it
    !allocate(hvm1(max_local_rows*max_stored_rows), stat=istat, errmsg=errorMessage)
    !call check_alloc("trans_ev_&
    !&complex&
    !&", "hvm1", istat, errorMessage)
    successGPU = gpu_malloc(tmat_dev, max_stored_rows * max_stored_rows * size_of_datatype)
    call check_alloc_GPU_f("trans_ev", 362,  successGPU)

    successGPU = gpu_malloc(hvm_dev, max_local_rows * max_stored_rows * size_of_datatype)
    call check_alloc_GPU_f("trans_ev", 365,  successGPU)

    successGPU = gpu_malloc(tmp_dev, max_local_cols * max_stored_rows * size_of_datatype)
    call check_alloc_GPU_f("trans_ev", 368,  successGPU)

    num = ldq * matrixCols * size_of_datatype
    successGPU = gpu_malloc(q_dev, num)
    call check_alloc_GPU_f("trans_ev", 372,  successGPU)
  

    successGPU = gpu_memcpy(q_dev, int(loc(q_mat(1,1)),kind=c_intptr_t), &
                  num, gpuMemcpyHostToDevice)
    call check_memcpy_GPU_f("trans_ev", 394,  successGPU)
  endif  ! useGPU

  do istep = 1, na, blockStep


    ics = MAX(istep,3)
    ice = MIN(istep+nblk-1,na)
    if (ice<ics) cycle

    cur_pcol = pcol(istep, nblk, np_cols)

    nb = 0
    do ic = ics, ice

      l_colh = local_index(ic  , my_pcol, np_cols, nblk, -1) ! Column of Householder Vector
      l_rows = local_index(ic-1, my_prow, np_rows, nblk, -1) ! # rows of Householder Vector


      if (my_pcol == cur_pcol) then
        hvb(nb+1:nb+l_rows) = a_mat(1:l_rows,l_colh)
        if (my_prow == prow(ic-1, nblk, np_rows)) then
          hvb(nb+l_rows) = 1.
        endif
      endif

      nb = nb+l_rows
    enddo


    nb = 0
    do ic = ics, ice
      l_rows = local_index(ic-1, my_prow, np_rows, nblk, -1) ! # rows of Householder Vector
      hvm(1:l_rows,nstor+1) = hvb(nb+1:nb+l_rows)
      if (useGPU) then
        hvm_ubnd = l_rows
      endif
      nstor = nstor+1
      nb = nb+l_rows
    enddo

    ! Please note: for smaller matix sizes (na/np_rows<=256), a value of 32 for nstor is enough!
    if (nstor+nblk > max_stored_rows .or. istep+nblk > na .or. (na/np_rows <= 256 .and. nstor >= 32)) then

      ! Calculate scalar products of stored vectors.
      ! This can be done in different ways, we use dsyrk or zherk

      tmat = 0
      call obj%timer%start("blas")
      if (l_rows>0) then
        call CHERK('U', 'C',   &
                         int(nstor,kind=BLAS_KIND), int(l_rows,kind=BLAS_KIND), ONE, &
                         hvm, int(ubound(hvm,dim=1),kind=BLAS_KIND), ZERO, tmat, int(max_stored_rows,kind=BLAS_KIND))
      endif
      call obj%timer%stop("blas")
      nc = 0
      do n = 1, nstor-1
        h1(nc+1:nc+n) = tmat(1:n,n+1)
        nc = nc+n
      enddo

      if (nc > 0) h2 = h1

      ! Calculate triangular matrix T

      nc = 0
      tmat(1,1) = tau(ice-nstor+1)
      do n = 1, nstor-1
        call obj%timer%start("blas")
        call CTRMV('L', 'C' , 'N', int(n,kind=BLAS_KIND), tmat, &
                            int(max_stored_rows,kind=BLAS_KIND), h2(nc+1), 1_BLAS_KIND)
        call obj%timer%stop("blas")

        tmat(n+1,1:n) = &
        -conjg(h2(nc+1:nc+n)) &
        *tau(ice-nstor+n+1)

        tmat(n+1,n+1) = tau(ice-nstor+n+1)
        nc = nc+n
      enddo

      if (useGPU) then
        ! todo: is this reshape really neccessary?
        hvm1(1:hvm_ubnd*nstor) = reshape(hvm(1:hvm_ubnd,1:nstor), (/ hvm_ubnd*nstor /))

        !hvm_dev(1:hvm_ubnd*nstor) = hvm1(1:hvm_ubnd*nstor)
        successGPU = gpu_memcpy(hvm_dev, int(loc(hvm1(1)),kind=c_intptr_t),   &
                      hvm_ubnd * nstor * size_of_datatype, gpuMemcpyHostToDevice)

        call check_memcpy_GPU_f("trans_ev", 546,  successGPU)

        !tmat_dev = tmat
        successGPU = gpu_memcpy(tmat_dev, int(loc(tmat(1,1)),kind=c_intptr_t),   &
                      max_stored_rows * max_stored_rows * size_of_datatype, gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("trans_ev", 551,  successGPU)
      endif


      ! Q = Q - V * T * V**T * Q

      if (l_rows > 0) then
        if (useGPU) then
          call obj%timer%start("gpublas")
          gpuHandle = obj%gpu_setup%gpublasHandleArray(0)
          call gpublas_CGEMM('C', 'N',   &
                                   nstor, l_cols, l_rows, ONE, hvm_dev, hvm_ubnd,  &
                                   q_dev, ldq, ZERO, tmp_dev, nstor, gpuHandle)
          call obj%timer%stop("gpublas")
        else ! useGPU

          call obj%timer%start("blas")
          call CGEMM('C', 'N',  &
                              int(nstor,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), &
                              int(l_rows,kind=BLAS_KIND), ONE, hvm, int(ubound(hvm,dim=1),kind=BLAS_KIND), &
                              q_mat, int(ldq,kind=BLAS_KIND), ZERO, tmp1, int(nstor,kind=BLAS_KIND))
          call obj%timer%stop("blas")
        endif ! useGPU

      else !l_rows>0

        if (useGPU) then
          if (gpu_vendor() /= OPENMP_OFFLOAD_GPU) then
            successGPU = gpu_memset(tmp_dev, 0, l_cols * nstor * size_of_datatype)
            call check_memcpy_GPU_f("trans_ev", 590,  successGPU)
          else
            allocate(tmp_debug(l_cols * nstor))
            tmp_debug(:) = 0.
            successGPU = gpu_memcpy(tmp_dev, int(loc(tmp_debug),kind=c_intptr_t), &
                                    l_cols*nstor*size_of_datatype, gpuMemcpyHostToDevice)
            call check_memcpy_GPU_f("trans_ev", 596,  successGPU)
            deallocate(tmp_debug)
          endif
        else
          tmp1(1:l_cols*nstor) = 0
        endif
      endif  !l_rows>0

!     tmp2 = tmp1

      if (l_rows > 0) then
        if (useGPU) then
          call obj%timer%start("gpublas")
          gpuHandle = obj%gpu_setup%gpublasHandleArray(0)
          call gpublas_CTRMM('L', 'L', 'N', 'N',     &
                                   nstor, l_cols, ONE, tmat_dev, max_stored_rows,  &
                                   tmp_dev, nstor, gpuHandle)

          call gpublas_CGEMM('N', 'N' ,l_rows ,l_cols ,nstor,  &
                                   -ONE, hvm_dev, hvm_ubnd, tmp_dev, nstor,   &
                                   ONE, q_dev, ldq, gpuHandle)
          call obj%timer%stop("gpublas")
        else !useGPU
          call obj%timer%start("blas")

          call CTRMM('L', 'L', 'N', 'N', int(nstor,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND),   &
                              ONE, tmat, int(max_stored_rows,kind=BLAS_KIND), tmp1, int(nstor,kind=BLAS_KIND))
          call CGEMM('N', 'N', int(l_rows,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), &
                              int(nstor,kind=BLAS_KIND), -ONE, hvm, int(ubound(hvm,dim=1),kind=BLAS_KIND), &
                              tmp1, int(nstor,kind=BLAS_KIND), ONE, q_mat, int(ldq,kind=BLAS_KIND))
          call obj%timer%stop("blas")
        endif ! useGPU
      endif  ! l_rows>0
      nstor = 0
    endif  ! (nstor+nblk>max_stored_rows .or. istep+nblk>na .or. (na/np_rows<=256 .and. nstor>=32))

    
  enddo ! istep = 1, na, blockStep

  deallocate(h1, h2, hvb, hvm, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_&     &MATH_DATATYPE&     &: h1, h2, hvb, hvm", 848,  istat,  errorMessage)

  if (useGPU) then

    !q_mat = q_dev
    successGPU = gpu_memcpy(int(loc(q_mat(1,1)),kind=c_intptr_t), &
                  q_dev, ldq * matrixCols * size_of_datatype, gpuMemcpyDeviceToHost)
    call check_memcpy_GPU_f("trans_ev", 864,  successGPU)

    !deallocate(hvm1, stat=istat, errmsg=errorMessage)
    !if (istat .ne. 0) then
    !  print *,"trans_ev_&
    !  &complex&
    !  &: error when deallocating hvm1 "//errorMessage
    !  stop 1
    !endif

    !deallocate(q_dev, tmp_dev, hvm_dev, tmat_dev)
    successGPU = gpu_free(q_dev)
    call check_dealloc_GPU_f("trans_ev", 906,  successGPU)

    successGPU = gpu_free(tmp_dev)
    call check_dealloc_GPU_f("trans_ev", 909,  successGPU)

    successGPU = gpu_free(hvm_dev)
    call check_dealloc_GPU_f("trans_ev", 912,  successGPU)

    successGPU = gpu_free(tmat_dev)
    call check_dealloc_GPU_f("trans_ev", 915,  successGPU)
  else ! useGPU
    deallocate(tmat, tmp1, tmp2, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("trans_ev_&     &MATH_DATATYPE&     &: tmat, tmp1, tmp2", 920,  istat,  errorMessage)
  endif ! useGPU

!#if defined (USE_CCL_TRANS_EV) && !defined(WITH_GPU_STREAMS)
!  successGPU = cuda_stream_destroy(obj%gpu_setup%my_stream)
!  if (.not.(successGPU)) then
!    print *,"Cannot destroy gpu stream handle"
!  endif
!#endif

  call obj%timer%stop("trans_ev_&
  &complex&
  &" // &
  &"_single" // &
  gpuString )

end subroutine trans_ev_&
&complex&
&_&
&single







subroutine hh_transform_complex_&
&single &
(obj, alpha, xnorm_sq, xf, tau, wantDebug)
  ! Similar to LAPACK routine ZLARFP, but uses ||x||**2 instead of x(:)
  ! and returns the factor xf by which x has to be scaled.
  ! It also hasn't the special handling for numbers < 1.d-300 or > 1.d150
  ! since this would be expensive for the parallel implementation.
  use precision
  use elpa_abstract_impl
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
  class(elpa_abstract_impl_t), intent(inout)    :: obj
  logical, intent(in)                           :: wantDebug
  complex(kind=ck), intent(inout) :: alpha
  real(kind=rk), intent(in)          :: xnorm_sq
  complex(kind=ck), intent(out)   :: xf, tau
  real(kind=rk)                      :: ALPHR, ALPHI

  real(kind=rk)                      :: BETA

  if (wantDebug) call obj%timer%start("hh_transform_&
                   &complex&
     	      &" // &
                   &"_single" )

  ALPHR = real( ALPHA, kind=rk )
  ALPHI = AIMAG( ALPHA )

  if ( XNORM_SQ==0.0_rk .AND. ALPHI==0.0_rk ) then

    if ( ALPHR>=0.0_rk ) then
      TAU = 0.0_rk
    else
      TAU = 2.0_rk
      ALPHA = -ALPHA
    endif
    XF = 0.0_rk

  else

    BETA = SIGN( SQRT( ALPHR**2 + ALPHI**2 + XNORM_SQ ), ALPHR )
    ALPHA = ALPHA + BETA
    IF ( BETA<0 ) THEN
      BETA = -BETA
      TAU  = -ALPHA / BETA
    ELSE
      ALPHR = ALPHI * (ALPHI/real( ALPHA , kind=rk))
      ALPHR = ALPHR + XNORM_SQ/real( ALPHA, kind=rk )

      TAU = CMPLX( ALPHR/BETA, -ALPHI/BETA )
      ALPHA = CMPLX( -ALPHR, ALPHI )
    END IF
    XF = 1.0_rk/ALPHA
    ALPHA = BETA
  endif

  if (wantDebug) call obj%timer%stop("hh_transform_&
  &complex&
  &" // &
  &"_single" )

    end subroutine hh_transform_complex_&
    &single



end module elpa1_compute
