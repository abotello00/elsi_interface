! Copyright (c) 2015-2024, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! Provide stub routines which are only compiled when the actual DLA-Future is not
!! available.
!!
module ELSI_DLAF

   use ELSI_DATATYPE, only: elsi_param_t, elsi_basic_t
   use ELSI_PRECISION, only: r8

   implicit none

   private

   public :: elsi_init_dlaf
   public :: elsi_cleanup_dlaf
   public :: elsi_solve_dlaf

   interface elsi_solve_dlaf
      module procedure elsi_solve_dlaf_real
      module procedure elsi_solve_dlaf_cmplx
   end interface

contains

   subroutine elsi_init_dlaf(ph, bh)

      implicit none

      type(elsi_param_t) :: ph
      type(elsi_basic_t) :: bh

   end subroutine

   subroutine elsi_solve_dlaf_real(ph, bh, ham, ovlp, eval, evec)

      implicit none

      type(elsi_param_t), intent(inout) :: ph
      type(elsi_basic_t), intent(in) :: bh
      real(kind=r8), intent(inout) :: ham(bh%n_lrow, bh%n_lcol) ! A
      real(kind=r8), intent(inout) :: ovlp(bh%n_lrow, bh%n_lcol) ! B
      real(kind=r8), intent(out) :: eval(ph%n_basis)
      real(kind=r8), intent(out) :: evec(bh%n_lrow, bh%n_lcol)

      write (*, "(A)") "**Error! A DLA-Future stub routine was called"
      stop

   end subroutine

   subroutine elsi_solve_dlaf_cmplx(ph, bh, ham, ovlp, eval, evec)

      implicit none

      type(elsi_param_t), intent(inout) :: ph
      type(elsi_basic_t), intent(in) :: bh
      complex(kind=r8), intent(inout) :: ham(bh%n_lrow, bh%n_lcol) ! A
      complex(kind=r8), intent(inout) :: ovlp(bh%n_lrow, bh%n_lcol) ! B
      real(kind=r8), intent(out) :: eval(ph%n_basis)
      complex(kind=r8), intent(out) :: evec(bh%n_lrow, bh%n_lcol)

      write (*, "(A)") "**Error! A DLA-Future stub routine was called"
      stop

   end subroutine

   subroutine elsi_cleanup_dlaf(ph, bh)

      implicit none

      type(elsi_basic_t) :: bh
      type(elsi_param_t) :: ph

   end subroutine

end module ELSI_DLAF
