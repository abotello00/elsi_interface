! Copyright (c) 2015-2024, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! Interface to DLA-Future.
!!
module ELSI_DLAF

   use dlaf_fortran, only: dlaf_initialize, dlaf_finalize, &
                           dlaf_create_grid_from_blacs, dlaf_free_grid, &
                           dlaf_pdsyevd, dlaf_pdsygvd, dlaf_pdsygvd_factorized, &
                           dlaf_pzheevd, dlaf_pzhegvd, dlaf_pzhegvd_factorized

   use ELSI_DATATYPE, only: elsi_param_t, elsi_basic_t
   use ELSI_MPI, only: elsi_stop
   use ELSI_OUTPUT, only: elsi_say, elsi_get_time
   use ELSI_PRECISION, only: r8, i4

   implicit none

   private

   public :: elsi_init_dlaf, elsi_cleanup_dlaf
   public :: elsi_solve_dlaf

   interface elsi_solve_dlaf
      module procedure elsi_solve_dlaf_real
      module procedure elsi_solve_dlaf_cmplx
   end interface

contains

!>
!! Initialize DLA-Future.
!!
   subroutine elsi_init_dlaf(ph, bh)

      implicit none

      type(elsi_param_t), intent(inout) :: ph
      type(elsi_basic_t), intent(in) :: bh

      character(len=*), parameter :: caller = "elsi_init_dlaf"

      if (.not. ph%dlaf_started) then
         call dlaf_initialize()
         ph%dlaf_started = .true.
      end if

      if (.not. ph%dlaf_grid_created) then
         call dlaf_create_grid_from_blacs(bh%blacs_ctxt)
         ph%dlaf_grid_created = .true.
      end if

   end subroutine

!>
!! Interface to DLA-Future.
!!
   subroutine elsi_solve_dlaf_real(ph, bh, ham, ovlp, eval, evec)

      implicit none

      type(elsi_param_t), intent(inout) :: ph
      type(elsi_basic_t), intent(in) :: bh
      real(kind=r8), intent(inout) :: ham(bh%n_lrow, bh%n_lcol) ! A
      real(kind=r8), intent(inout) :: ovlp(bh%n_lrow, bh%n_lcol) ! B
      real(kind=r8), intent(out) :: eval(ph%n_basis)
      real(kind=r8), intent(out) :: evec(bh%n_lrow, bh%n_lcol)

      real(kind=r8) :: t0
      real(kind=r8) :: t1
      integer(kind=i4) :: ierr
      character(len=200) :: msg

      character(len=*), parameter :: caller = "elsi_solve_dlaf_real"

      write (msg, "(A)") "Starting DLA-Future eigensolver"
      call elsi_say(bh, msg)

      call elsi_get_time(t0)

      ! Solve
      if (ph%unit_ovlp) then
         call dlaf_pdsyevd("L", ph%n_basis, ham, 1, 1, bh%desc, eval, evec, 1, 1, bh%desc, ierr)
      else
         write (msg, *) "DEBUG", bh%desc, ph%n_basis
         call elsi_say(bh, msg)
         if (ph%dlaf_first) then
            call dlaf_pdsygvd("L", ph%n_basis, ham, 1, 1, bh%desc, ovlp, 1, 1, bh%desc, eval, evec, 1, 1, bh%desc, ierr)
         else
            ! Overlap matrix contains Cholesky factorization
            call dlaf_pdsygvd_factorized("L", ph%n_basis, ham, 1, 1, bh%desc, ovlp, 1, 1, bh%desc, eval, evec, 1, 1, bh%desc, ierr)
         end if
      end if

      if (ierr /= 0) then
         write (msg, "(A)") "DLA-Future eigensolver failed"
         call elsi_stop(bh, msg, caller)
      end if

      evec(:, :) = ham

      call elsi_get_time(t1)

      write (msg, "(A)") "Finished solving eigenproblem"
      call elsi_say(bh, msg)
      write (msg, "(A,F10.3,A)") "| Time :", t1 - t0, " s"
      call elsi_say(bh, msg)

      ph%dlaf_first = .false.

   end subroutine

!>
!! Interface to DLA-Future.
!!
   subroutine elsi_solve_dlaf_cmplx(ph, bh, ham, ovlp, eval, evec)

      implicit none

      type(elsi_param_t), intent(inout) :: ph
      type(elsi_basic_t), intent(in) :: bh
      complex(kind=r8), intent(inout) :: ham(bh%n_lrow, bh%n_lcol) ! A
      complex(kind=r8), intent(inout) :: ovlp(bh%n_lrow, bh%n_lcol) ! B
      real(kind=r8), intent(out) :: eval(ph%n_basis)
      complex(kind=r8), intent(out) :: evec(bh%n_lrow, bh%n_lcol)

      real(kind=r8) :: t0
      real(kind=r8) :: t1
      integer(kind=i4) :: ierr
      character(len=200) :: msg

      character(len=*), parameter :: caller = "elsi_solve_dlaf_complex"

      write (msg, "(A)") "Starting DLA-Future eigensolver"
      call elsi_say(bh, msg)

      call elsi_get_time(t0)

      ! Solve
      if (ph%unit_ovlp) then
         call dlaf_pzheevd("L", ph%n_basis, ham, 1, 1, bh%desc, eval, evec, 1, 1, bh%desc, ierr)
      else
         write (msg, *) "DEBUG", bh%desc, ph%n_basis
         call elsi_say(bh, msg)
         if (ph%dlaf_first) then
            call dlaf_pzhegvd("L", ph%n_basis, ham, 1, 1, bh%desc, ovlp, 1, 1, bh%desc, eval, evec, 1, 1, bh%desc, ierr)
         else
            ! Overlap matrix contains Cholesky factorization
            call dlaf_pzhegvd_factorized("L", ph%n_basis, ham, 1, 1, bh%desc, ovlp, 1, 1, bh%desc, eval, evec, 1, 1, bh%desc, ierr)
         end if
      end if

      if (ierr /= 0) then
         write (msg, "(A)") "DLA-Future eigensolver failed"
         call elsi_stop(bh, msg, caller)
      end if

      evec(:, :) = ham

      call elsi_get_time(t1)

      write (msg, "(A)") "Finished solving eigenproblem"
      call elsi_say(bh, msg)
      write (msg, "(A,F10.3,A)") "| Time :", t1 - t0, " s"
      call elsi_say(bh, msg)

      ph%dlaf_first = .false.

   end subroutine

!>
!! Clean up DLA-Future.
!!
   subroutine elsi_cleanup_dlaf(ph, bh)

      implicit none

      type(elsi_param_t), intent(inout) :: ph
      type(elsi_basic_t), intent(in) :: bh

      character(len=*), parameter :: caller = "elsi_cleanup_dlaf"

      if (ph%dlaf_started) then
         call dlaf_finalize()
         ph%dlaf_started = .false.
      end if

      if (ph%dlaf_grid_created) then
         call dlaf_free_grid(bh%blacs_ctxt)
         ph%dlaf_grid_created = .false.
      end if

   end subroutine

end module ELSI_DLAF
