! Copyright (c) 2015-2021, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! Perform wavefunction and density matrix extrapolation.
!!
module ELSI_GEO

   use ELSI_CONSTANT, only: PEXSI_CSC,SIESTA_CSC,GENERIC_COO,LT_MAT,UT_MAT
   use ELSI_DATATYPE, only: elsi_handle
   use ELSI_ELPA, only: elsi_update_dm_elpa, elsi_factor_ovlp_elpa
   use ELSI_MALLOC, only: elsi_allocate,elsi_deallocate
   use ELSI_MPI
   use ELSI_NTPOLY, only: elsi_update_dm_ntpoly
   use ELSI_PRECISION, only: r8
   use ELSI_REDIST, only: elsi_generic_to_blacs_hs,elsi_generic_to_ntpoly_hs,&
       elsi_ntpoly_to_generic_dm,elsi_ntpoly_to_siesta_dm,&
       elsi_ntpoly_to_sips_dm,elsi_siesta_to_blacs_hs,elsi_siesta_to_ntpoly_hs,&
       elsi_sips_to_blacs_hs,elsi_sips_to_ntpoly_hs
   use ELSI_UTIL, only: elsi_check_init,elsi_gram_schmidt,elsi_set_full_mat,elsi_check_err
   use ELSI_OUTPUT, only: elsi_say,elsi_get_time

   implicit none

   private

   public :: elsi_orthonormalize_ev_real
   public :: elsi_orthonormalize_ev_complex
   public :: elsi_orthonormalize_ev_real_sparse
   public :: elsi_orthonormalize_ev_complex_sparse
   public :: elsi_extrapolate_dm_real
   public :: elsi_extrapolate_dm_complex
   public :: elsi_extrapolate_dm_real_sparse
   public :: elsi_extrapolate_dm_complex_sparse
   !KL: new routines for extended langrangian MD approach
   public :: elsi_get_ccdm_real
   public :: elsi_get_ccdm_complex
   public :: elsi_get_dm_from_ccdm_real
   public :: elsi_get_dm_from_ccdm_complex
   public :: elsi_residual_ccdm_real
   public :: elsi_residual_ccdm_complex
   public :: elsi_ccdm_scaled_delta_approx_real
   public :: elsi_ccdm_scaled_delta_approx_complex
   public :: elsi_ccdm_dissipation_real
   public :: elsi_ccdm_dissipation_complex
   public :: elsi_init_ccdm_dissipation
   public :: elsi_ccdm_integration_real
   public :: elsi_ccdm_integration_complex
   public :: elsi_trace_glob_dm_real
   public :: elsi_trace_glob_dm_complex
   public :: elsi_norm_residual_real
   public :: elsi_norm_residual_complex
   public :: elsi_edm_from_dm_real
   public :: elsi_edm_from_dm_complex
   public :: elsi_get_inverse_ovlp_real_2
   public :: elsi_get_inverse_ovlp_complex_2
   public :: elsi_inv_ovlp_response_real 
   public :: elsi_inv_ovlp_response_complex


contains

!>
!! Orthonormalize eigenvectors with respect to an overlap matrix.
!!
subroutine elsi_orthonormalize_ev_real(eh,ovlp,evec)

   implicit none

   type(elsi_handle), intent(in) :: eh !< Handle
   real(kind=r8), intent(in) :: ovlp(eh%bh%n_lrow,eh%bh%n_lcol) !< Overlap
   real(kind=r8), intent(inout) :: evec(eh%bh%n_lrow,eh%bh%n_lcol) !< Eigenvectors

   character(len=*), parameter :: caller = "elsi_orthonormalize_ev_real"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   call elsi_gram_schmidt(eh%ph,eh%bh,ovlp,evec)

end subroutine

!>
!! Orthonormalize eigenvectors with respect to an overlap matrix.
!!
subroutine elsi_orthonormalize_ev_complex(eh,ovlp,evec)

   implicit none

   type(elsi_handle), intent(in) :: eh !< Handle
   complex(kind=r8), intent(in) :: ovlp(eh%bh%n_lrow,eh%bh%n_lcol) !< Overlap
   complex(kind=r8), intent(inout) :: evec(eh%bh%n_lrow,eh%bh%n_lcol) !< Eigenvectors

   character(len=*), parameter :: caller = "elsi_orthonormalize_ev_complex"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   call elsi_gram_schmidt(eh%ph,eh%bh,ovlp,evec)

end subroutine

!>
!! Orthonormalize eigenvectors with respect to an overlap matrix.
!!
subroutine elsi_orthonormalize_ev_real_sparse(eh,ovlp,evec)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: ovlp(eh%bh%nnz_l_sp) !< New overlap
   real(kind=r8), intent(inout) :: evec(eh%bh%n_lrow,eh%bh%n_lcol) !< Eigenvectors

   real(kind=r8) :: dummy1(1)
   real(kind=r8) :: dummy2(1,1)
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_orthonormalize_ev_real_sparse"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(.not. allocated(eh%ovlp_real_den)) then
      call elsi_allocate(eh%bh,eh%ovlp_real_den,eh%bh%n_lrow,eh%bh%n_lcol,&
           "ovlp_real_den",caller)
   end if

   eh%ph%unit_ovlp = .true.

   select case(eh%ph%matrix_format)
   case(PEXSI_CSC)
      call elsi_sips_to_blacs_hs(eh%ph,eh%bh,ovlp,dummy1,eh%row_ind_sp1,&
           eh%col_ptr_sp1,eh%ovlp_real_den,dummy2)
   case(SIESTA_CSC)
      call elsi_siesta_to_blacs_hs(eh%ph,eh%bh,ovlp,dummy1,eh%row_ind_sp2,&
           eh%col_ptr_sp2,eh%ovlp_real_den,dummy2)
   case(GENERIC_COO)
      if(.not. allocated(eh%map_den)) then
         call elsi_allocate(eh%bh,eh%map_den,eh%bh%n_lrow,eh%bh%n_lcol,&
              "map_den",caller)
      end if

      call elsi_generic_to_blacs_hs(eh%ph,eh%bh,ovlp,dummy1,eh%row_ind_sp3,&
           eh%col_ind_sp3,eh%ovlp_real_den,dummy2,eh%map_den)
   case default
      write(msg,"(A)") "Unsupported matrix format"
      call elsi_stop(eh%bh,msg,caller)
   end select

   eh%ph%unit_ovlp = .false.

   call elsi_gram_schmidt(eh%ph,eh%bh,eh%ovlp_real_den,evec)

end subroutine

!>
!! Orthonormalize eigenvectors with respect to an overlap matrix.
!!
subroutine elsi_orthonormalize_ev_complex_sparse(eh,ovlp,evec)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   complex(kind=r8), intent(in) :: ovlp(eh%bh%nnz_l_sp) !< New overlap
   complex(kind=r8), intent(inout) :: evec(eh%bh%n_lrow,eh%bh%n_lcol) !< Eigenvectors

   complex(kind=r8) :: dummy1(1)
   complex(kind=r8) :: dummy2(1,1)
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_orthonormalize_ev_complex_sparse"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(.not. allocated(eh%ovlp_cmplx_den)) then
      call elsi_allocate(eh%bh,eh%ovlp_cmplx_den,eh%bh%n_lrow,eh%bh%n_lcol,&
           "ovlp_cmplx_den",caller)
   end if

   eh%ph%unit_ovlp = .true.

   select case(eh%ph%matrix_format)
   case(PEXSI_CSC)
      call elsi_sips_to_blacs_hs(eh%ph,eh%bh,ovlp,dummy1,eh%row_ind_sp1,&
           eh%col_ptr_sp1,eh%ovlp_cmplx_den,dummy2)
   case(SIESTA_CSC)
      call elsi_siesta_to_blacs_hs(eh%ph,eh%bh,ovlp,dummy1,eh%row_ind_sp2,&
           eh%col_ptr_sp2,eh%ovlp_cmplx_den,dummy2)
   case(GENERIC_COO)
      if(.not. allocated(eh%map_den)) then
         call elsi_allocate(eh%bh,eh%map_den,eh%bh%n_lrow,eh%bh%n_lcol,&
              "map_den",caller)
      end if

      call elsi_generic_to_blacs_hs(eh%ph,eh%bh,ovlp,dummy1,eh%row_ind_sp3,&
           eh%col_ind_sp3,eh%ovlp_cmplx_den,dummy2,eh%map_den)
   case default
      write(msg,"(A)") "Unsupported matrix format"
      call elsi_stop(eh%bh,msg,caller)
   end select

   eh%ph%unit_ovlp = .false.

   call elsi_gram_schmidt(eh%ph,eh%bh,eh%ovlp_cmplx_den,evec)

end subroutine

!>
!! Extrapolate density matrix for a new overlap.
!!
subroutine elsi_extrapolate_dm_real(eh,ovlp,dm)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: ovlp(eh%bh%n_lrow,eh%bh%n_lcol) !< New overlap
   real(kind=r8), intent(out) :: dm(eh%bh%n_lrow,eh%bh%n_lcol) !< Density matrix

   real(kind=r8), allocatable :: tmp(:,:)

   character(len=*), parameter :: caller = "elsi_extrapolate_dm_real"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   call elsi_allocate(eh%bh,tmp,eh%bh%n_lrow,eh%bh%n_lcol,"tmp",caller)

   tmp(:,:) = ovlp

   call elsi_update_dm_elpa(eh%ph,eh%bh,eh%ovlp_real_copy,tmp,eh%dm_real_copy,&
        dm)

   call elsi_deallocate(eh%bh,tmp,"tmp")

end subroutine

!>
!! Extrapolate density matrix for a new overlap.
!!
subroutine elsi_extrapolate_dm_complex(eh,ovlp,dm)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   complex(kind=r8), intent(in) :: ovlp(eh%bh%n_lrow,eh%bh%n_lcol) !< New overlap
   complex(kind=r8), intent(out) :: dm(eh%bh%n_lrow,eh%bh%n_lcol) !< Density matrix

   complex(kind=r8), allocatable :: tmp(:,:)

   character(len=*), parameter :: caller = "elsi_extrapolate_dm_complex"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   call elsi_allocate(eh%bh,tmp,eh%bh%n_lrow,eh%bh%n_lcol,"tmp",caller)

   tmp(:,:) = ovlp

   call elsi_update_dm_elpa(eh%ph,eh%bh,eh%ovlp_cmplx_copy,tmp,&
        eh%dm_cmplx_copy,dm)

   call elsi_deallocate(eh%bh,tmp,"tmp")

end subroutine

!>
!! Extrapolate density matrix for a new overlap.
!!
subroutine elsi_extrapolate_dm_real_sparse(eh,ovlp,dm)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: ovlp(eh%bh%nnz_l_sp) !< New overlap
   real(kind=r8), intent(out) :: dm(eh%bh%nnz_l_sp) !< New density matrix

   real(kind=r8) :: dummy1(1)

   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_extrapolate_dm_real_sparse"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   eh%ph%unit_ovlp = .true.

   select case(eh%ph%matrix_format)
   case(PEXSI_CSC)
      call elsi_sips_to_ntpoly_hs(eh%ph,eh%bh,ovlp,dummy1,eh%row_ind_sp1,&
           eh%col_ptr_sp1,eh%nt_ovlp,eh%nt_dm)
   case(SIESTA_CSC)
      call elsi_siesta_to_ntpoly_hs(eh%ph,eh%bh,ovlp,dummy1,eh%row_ind_sp2,&
           eh%col_ptr_sp2,eh%nt_ovlp,eh%nt_dm)
   case(GENERIC_COO)
      call elsi_generic_to_ntpoly_hs(eh%ph,eh%bh,ovlp,dummy1,eh%row_ind_sp3,&
           eh%col_ind_sp3,eh%nt_ovlp,eh%nt_dm,eh%nt_map)
   case default
      write(msg,"(A)") "Unsupported matrix format"
      call elsi_stop(eh%bh,msg,caller)
   end select

   eh%ph%unit_ovlp = .false.

   call elsi_update_dm_ntpoly(eh%ph,eh%bh,eh%nt_ovlp_copy,eh%nt_ovlp,&
        eh%nt_dm_copy,eh%nt_dm)

   select case(eh%ph%matrix_format)
   case(PEXSI_CSC)
      call elsi_ntpoly_to_sips_dm(eh%ph,eh%bh,eh%nt_dm,dm,eh%row_ind_sp1,&
           eh%col_ptr_sp1)
   case(SIESTA_CSC)
      call elsi_ntpoly_to_siesta_dm(eh%bh,eh%nt_dm,dm,eh%row_ind_sp2,&
           eh%col_ptr_sp2)
   case(GENERIC_COO)
      call elsi_ntpoly_to_generic_dm(eh%ph,eh%bh,eh%nt_dm,eh%nt_map,dm,&
           eh%perm_sp3)
   case default
      write(msg,"(A)") "Unsupported matrix format"
      call elsi_stop(eh%bh,msg,caller)
   end select

end subroutine

!>
!! Extrapolate density matrix for a new overlap.
!!
subroutine elsi_extrapolate_dm_complex_sparse(eh,ovlp,dm)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   complex(kind=r8), intent(in) :: ovlp(eh%bh%nnz_l_sp) !< New overlap
   complex(kind=r8), intent(out) :: dm(eh%bh%nnz_l_sp) !< New density matrix

   complex(kind=r8) :: dummy1(1)

   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_extrapolate_dm_complex_sparse"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   eh%ph%unit_ovlp = .true.

   select case(eh%ph%matrix_format)
   case(PEXSI_CSC)
      call elsi_sips_to_ntpoly_hs(eh%ph,eh%bh,ovlp,dummy1,eh%row_ind_sp1,&
           eh%col_ptr_sp1,eh%nt_ovlp,eh%nt_dm)
   case(SIESTA_CSC)
      call elsi_siesta_to_ntpoly_hs(eh%ph,eh%bh,ovlp,dummy1,eh%row_ind_sp2,&
           eh%col_ptr_sp2,eh%nt_ovlp,eh%nt_dm)
   case(GENERIC_COO)
      call elsi_generic_to_ntpoly_hs(eh%ph,eh%bh,ovlp,dummy1,eh%row_ind_sp3,&
           eh%col_ind_sp3,eh%nt_ovlp,eh%nt_dm,eh%nt_map)
   case default
      write(msg,"(A)") "Unsupported matrix format"
      call elsi_stop(eh%bh,msg,caller)
   end select

   eh%ph%unit_ovlp = .false.

   call elsi_update_dm_ntpoly(eh%ph,eh%bh,eh%nt_ovlp_copy,eh%nt_ovlp,&
        eh%nt_dm_copy,eh%nt_dm)

   select case(eh%ph%matrix_format)
   case(PEXSI_CSC)
      call elsi_ntpoly_to_sips_dm(eh%ph,eh%bh,eh%nt_dm,dm,eh%row_ind_sp1,&
           eh%col_ptr_sp1)
   case(SIESTA_CSC)
      call elsi_ntpoly_to_siesta_dm(eh%bh,eh%nt_dm,dm,eh%row_ind_sp2,&
           eh%col_ptr_sp2)
   case(GENERIC_COO)
      call elsi_ntpoly_to_generic_dm(eh%ph,eh%bh,eh%nt_dm,eh%nt_map,dm,&
           eh%perm_sp3)
   case default
      write(msg,"(A)") "Unsupported matrix format"
      call elsi_stop(eh%bh,msg,caller)
   end select

end subroutine


!>
!! Here are some auxiliary routines which are used here for now but can be moved somewhere
!! else for general use.
!!

!>
!! Compute the energy-weighted DM from an initial DM. Uses the relationship
!! EDM = DM * H * DM 
!! that is also used in the normal EDM calculation with the NT-Poly library. (Sparse DM)
!! We need this since we are using an integrated DM as a new starting point in each MD step.
!! The EDM is only used in the calculation of Pulay forces.
subroutine elsi_edm_from_dm_real(eh,dm,ham,edm)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: dm(eh%bh%n_lrow,eh%bh%n_lcol) !< Density matrix
   real(kind=r8), intent(in) :: ham(eh%bh%n_lrow,eh%bh%n_lcol) !< Hamiltonian
   real(kind=r8), intent(out) :: edm(eh%bh%n_lrow,eh%bh%n_lcol) !< energy-weighted Density matrix

   character(len=*), parameter :: caller = "elsi_edm_from_dm_real"

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   real(kind=r8) :: factor

   real(kind=r8), allocatable :: tmp(:,:)

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   call elsi_allocate(eh%bh,tmp,eh%bh%n_lrow,eh%bh%n_lcol,"tmp",caller)

   factor = 1.0_r8/eh%ph%spin_degen ! Needed to scale the edm at the end. Formula for EDM assumes occupation 0-1

   call elsi_get_time(t0) ! Timestamp

   ! Do EDM = DM * H * DM
   call pdgemm("N","N",eh%ph%n_basis,eh%ph%n_basis,eh%ph%n_basis,1.0_r8,dm,1,&
        1,eh%bh%desc,ham,1,1,eh%bh%desc,0.0_r8,tmp,1,1,&
        eh%bh%desc)

   call pdgemm("N","N",eh%ph%n_basis,eh%ph%n_basis,eh%ph%n_basis,1.0_r8,tmp,1,&
        1,eh%bh%desc,dm,1,1,eh%bh%desc,0.0_r8,edm,1,1,&
        eh%bh%desc)

   ! Scale the edm at the end
   edm = factor * edm

   call elsi_get_time(t1) ! Timestamp

   write(msg,"(A)") "Finished EDM from DM calculation."

   call elsi_say(eh%bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(eh%bh,msg)

   call elsi_deallocate(eh%bh,tmp,"tmp")

end subroutine

subroutine elsi_edm_from_dm_complex(eh,dm,ham,edm)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   complex(kind=r8), intent(in) :: dm(eh%bh%n_lrow,eh%bh%n_lcol) !< Density matrix
   complex(kind=r8), intent(in) :: ham(eh%bh%n_lrow,eh%bh%n_lcol) !< Hamiltonian
   complex(kind=r8), intent(out) :: edm(eh%bh%n_lrow,eh%bh%n_lcol) !< energy-weighted Density matrix

   character(len=*), parameter :: caller = "elsi_get_ccdm_real"

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   real(kind=r8) :: factor

   complex(kind=r8), allocatable :: tmp(:,:)

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   call elsi_allocate(eh%bh,tmp,eh%bh%n_lrow,eh%bh%n_lcol,"tmp",caller)

   factor = 1.0_r8/eh%ph%spin_degen ! Needed to scale the edm at the end. Formula for EDM assumes occupation 0-1

   call elsi_get_time(t0) ! Timestamp

   ! Do EDM = DM * H * DM
   call pzgemm("N","N",eh%ph%n_basis,eh%ph%n_basis,eh%ph%n_basis,(1.0_r8,0.0_r8),dm,1,&
        1,eh%bh%desc,ham,1,1,eh%bh%desc,(0.0_r8,0.0_r8),tmp,1,1,&
        eh%bh%desc)

   call pzgemm("N","N",eh%ph%n_basis,eh%ph%n_basis,eh%ph%n_basis,(1.0_r8,0.0_r8),tmp,1,&
        1,eh%bh%desc,dm,1,1,eh%bh%desc,(0.0_r8,0.0_r8),edm,1,1,&
        eh%bh%desc)

   ! Scale the edm at the end
   edm = factor * edm

   call elsi_get_time(t1) ! Timestamp

   write(msg,"(A)") "Finished EDM from DM calculation."

   call elsi_say(eh%bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(eh%bh,msg)

   call elsi_deallocate(eh%bh,tmp,"tmp")

end subroutine

!>
!! Compute the trace of a global matrix with the dimensions of the DM. This can be
!! DM or CCDM etc. Outputs the trace upon running to every process.
!!
subroutine elsi_trace_glob_dm_real(eh,dm,trace)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: dm(eh%bh%n_lrow,eh%bh%n_lcol) !< Density matrix
   real(kind=r8), intent(out) :: trace !< Density matrix

   character(len=*), parameter :: caller = "elsi_trace_glob_dm_real"

   character(len=200) :: msg

   ! Local
   real(kind=r8), external :: pdlatra

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   trace = pdlatra(eh%ph%n_basis,dm,1,1,eh%bh%desc)

end subroutine

subroutine elsi_trace_glob_dm_complex(eh,dm,trace)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   complex(kind=r8), intent(in) :: dm(eh%bh%n_lrow,eh%bh%n_lcol) !< Density matrix
   complex(kind=r8), intent(out) :: trace !< Density matrix

   character(len=*), parameter :: caller = "elsi_trace_glob_dm_complex"

   character(len=200) :: msg

   ! Local
   complex(kind=r8), external :: pzlatra

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   trace = pzlatra(eh%ph%n_basis,dm,1,1,eh%bh%desc)

end subroutine

!>
!! Computes the Frobenius norm of a matrix (sqrt of squared sum of elements) 
!! in the form of the DM. 
!!

subroutine elsi_norm_residual_real(eh,residual,norm)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: residual(eh%bh%n_lrow,eh%bh%n_lcol) !< Residual
   real(kind=r8), intent(out) :: norm !< RMS error of residual

   character(len=*), parameter :: caller = "elsi_norm_residual_real"

   character(len=200) :: msg

   ! Local
   real(kind=r8), external :: pdlange
   real(kind=r8) :: factor

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   norm = pdlange("F",eh%ph%n_basis,eh%ph%n_basis,residual,1,1,eh%bh%desc,0)

   ! Also calculate the RMS error by scaling Froebnius norm , 1/sqrt(n_basis**2) = 1 / n_basis
   factor = eh%ph%n_basis
   norm = norm / factor

end subroutine

subroutine elsi_norm_residual_complex(eh,residual,norm)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   complex(kind=r8), intent(in) :: residual(eh%bh%n_lrow,eh%bh%n_lcol) !< Residual
   real(kind=r8), intent(out) :: norm !< RMS error of residual

   character(len=*), parameter :: caller = "elsi_norm_residual_complex"

   character(len=200) :: msg
   real(kind=r8) :: factor

   ! Local
   real(kind=r8), external :: pzlange

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   norm = pzlange("F",eh%ph%n_basis,eh%ph%n_basis,residual,1,1,eh%bh%desc,0)

   ! Also calculate the RMS error by scaling Froebnius norm , 1/sqrt(n_basis**2) = 1 / n_basis
   factor = eh%ph%n_basis
   norm = norm / factor

end subroutine


!>
!! KL: From this point on the routines are used for the extended langrangian MD
!! approach by Niklasson. For now this is stored here but could be moved somewhere
!! else once the src is split into categories.



!>
!! Get the matrix product of density matrix and overlap matrix. The resulting matrix
!! is sometimes called contra-covariant density matrix (ccdm).
!! ccdm = dm * ovlp
!!
subroutine elsi_get_ccdm_real(eh,ovlp,dm,ccdm)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: ovlp(eh%bh%n_lrow,eh%bh%n_lcol) !< Overlap matrix
   real(kind=r8), intent(in) :: dm(eh%bh%n_lrow,eh%bh%n_lcol) !< Density matrix
   real(kind=r8), intent(out) :: ccdm(eh%bh%n_lrow,eh%bh%n_lcol) !< contra-covariant density matrix

   character(len=*), parameter :: caller = "elsi_get_ccdm_real"

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   call elsi_get_time(t0) ! Timestamp

   !ccdm = matmul(dm,ovlp) ! For now simple fortran matrix-matrix multiplication
   call pdgemm("N","N",eh%ph%n_basis,eh%ph%n_basis,eh%ph%n_basis,1.0_r8,dm,1,&
        1,eh%bh%desc,ovlp,1,1,eh%bh%desc,0.0_r8,ccdm,1,1,&
        eh%bh%desc)

   call elsi_get_time(t1) ! Timestamp

   write(msg,"(A)") "Finished contra-covariant density matrix calculation."

   call elsi_say(eh%bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(eh%bh,msg)

end subroutine

!>
!! Get the matrix product of density matrix and overlap matrix. The resulting matrix
!! is sometimes called contra-covariant density matrix (ccdm).
!! ccdm = dm * ovlp
!!
subroutine elsi_get_ccdm_complex(eh,ovlp,dm,ccdm)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   complex(kind=r8), intent(in) :: ovlp(eh%bh%n_lrow,eh%bh%n_lcol) !< Overlap matrix
   complex(kind=r8), intent(in) :: dm(eh%bh%n_lrow,eh%bh%n_lcol) !< Density matrix
   complex(kind=r8), intent(out) :: ccdm(eh%bh%n_lrow,eh%bh%n_lcol) !< contra-covariant density matrix

   character(len=*), parameter :: caller = "elsi_get_ccdm_complex"

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   call elsi_get_time(t0) ! Timestamp

   !ccdm = matmul(dm,ovlp) ! For now simple fortran matrix-matrix multiplication
   call pzgemm("N","N",eh%ph%n_basis,eh%ph%n_basis,eh%ph%n_basis,(1.0_r8,0.0_r8),dm,1,&
        1,eh%bh%desc,ovlp,1,1,eh%bh%desc,(0.0_r8,0.0_r8),ccdm,1,1,&
        eh%bh%desc)

   call elsi_get_time(t1) ! Timestamp

   write(msg,"(A)") "Finished contra-covariant density matrix calculation."

   call elsi_say(eh%bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(eh%bh,msg)

end subroutine

!>
!! Get the inverse of the overlap matrix 
!! DO I NEED AN ILL_CONDITIONED CHECK HERE?
!!

subroutine elsi_get_inverse_ovlp_real(eh,ovlp,inv_ovlp)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: ovlp(eh%bh%n_lrow,eh%bh%n_lcol) !< Overlap matrix
   real(kind=r8), intent(out) :: inv_ovlp(eh%bh%n_lrow,eh%bh%n_lcol) !< contra-covariant Density matrix

   character(len=*), parameter :: caller = "elsi_get_inverse_ovlp_real"

   real(kind=r8), allocatable :: tmp(:,:), tmp1(:,:)

   integer :: i_k, j
   integer(kind=i4) :: ierr
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   real(kind=r8) :: trace, norm
   character(len=200) :: msg
   logical :: use_scalapack_intern

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   call elsi_get_time(t0) ! Timestamp

   ! Allocate tmp arrays
   call elsi_allocate(eh%bh,tmp,eh%bh%n_lrow,eh%bh%n_lcol,"tmp",caller)
   call elsi_allocate(eh%bh,tmp1,eh%bh%n_lrow,eh%bh%n_lcol,"tmp1",caller)

   tmp(:,:) = ovlp ! Temporary overlap

   ! write(msg,"(A)") "Writing S here as result."
   ! call elsi_say(eh%bh,msg)
   ! do i_k = 1,eh%bh%n_lrow
   !     write(*,*) "ELSI", ( tmp(i_k,j), j=1,eh%bh%n_lcol ), eh%bh%myid
   ! end do

   ! Maybe ill-conditioned check here with ?
   !call elsi_check_ovlp_elpa(eh%ph,eh%bh,tmp,eval_tmp,evec_tmp)

   ! ovlp = U^(-1)
   call elsi_factor_ovlp_elpa(eh%ph,eh%bh,tmp) ! ELPA routine to get inverse
   tmp1(:,:) = tmp

   ! Get ovlp^(-1) = U^(-1)*(U^(-1))^T
   ! Here we are multiplying an upper- times a lower-trianglular matrix. Instead of pdgemm it might be good to use
   ! specific routines for upper, lower matrices instead
   !tmp = matmul(tmp,transpose(tmp)) ! Simple matrix-matrix multiplication for now
   ! call pdgemm("N","T",eh%ph%n_basis,eh%ph%n_basis,eh%ph%n_basis,1.0_r8,tmp,1,&
   !      1,eh%bh%desc,tmp1,1,1,eh%bh%desc,0.0_r8,tmp2,1,1,&
   !      eh%bh%desc)
   ! Try with specific routine, this should be checked in more detail
   ! tmp1 = U^(-T)
   call pdtran(eh%ph%n_basis,eh%ph%n_basis,1.0_r8,tmp,1,1,eh%bh%desc,0.0_r8,tmp1,1,1,&
           eh%bh%desc)
   tmp(:,:) = tmp1 
   ! tmp2 = U^(-1)*(U^(-1))^T
   call eh%ph%elpa_aux%hermitian_multiply("L","L",eh%ph%n_basis,tmp,tmp1,&
           eh%bh%n_lrow,eh%bh%n_lcol,inv_ovlp,eh%bh%n_lrow,eh%bh%n_lcol,ierr)
   call elsi_check_err(eh%bh,"ELPA matrix multiplication",ierr,caller)

   ! ovlp^(-1) should be symmetric in the end again, Several options to ensure that
   ! Option 1 is just setting one triangular half of the matrix
   call elsi_set_full_mat(eh%ph,eh%bh,LT_MAT,inv_ovlp) ! Full matrix
   ! Option 2 is symmetrizing by 0.5 * (ovlp + ovlp^T)
   ! tmp1(:,:) = tmp2
   ! call pdgeadd("T",eh%ph%n_basis,eh%ph%n_basis,1.0_r8,tmp1,1,1,eh%bh%desc,1.0_r8,tmp2,1,1,eh%bh%desc)
   ! tmp2(:,:) = 0.5_r8 * tmp2

   ! Check if overlap times inverse is identity via the trace, Should implement check also for off-diagonal elements
   ! tmp = S * S^{-1}
   tmp = 0.0_r8
   call pdgemm("N","N",eh%ph%n_basis,eh%ph%n_basis,eh%ph%n_basis,1.0_r8,ovlp,1,&
        1,eh%bh%desc,inv_ovlp,1,1,eh%bh%desc,0.0_r8,tmp,1,1,&
        eh%bh%desc)
   ! Get Trace
   call elsi_trace_glob_dm_real(eh,tmp,trace)
   !Check if trace is equal to n_basis as check
   write(msg,"(A,E10.3)") "***Overlap identity via trace:", (abs(real(trace))-eh%ph%n_basis)
   call elsi_say(eh%bh,msg)
   if ((abs(trace)-eh%ph%n_basis)>1e-10_r8) then
        write(msg,"(A,E10.3)") "***Overlap is not correct. Difference to identity via trace:", (abs(trace)-eh%ph%n_basis)
        call elsi_say(eh%bh,msg)
   end if
   ! Check RMS of norm to be 1
   call elsi_norm_residual_real(eh,tmp,norm)
   ! This should be 1
   norm = (norm ** 2) * eh%ph%n_basis
   !Check if trace is equal to n_basis as check
   write(msg,"(A,E10.3)") "***Overlap identity via norm:", (norm-1.0_r8)
   call elsi_say(eh%bh,msg)
   if ((norm-1.0_r8)>1e-10_r8) then
        write(msg,"(A,E10.3)") "***Overlap is not correct. Difference to identity via norm:", (norm-1.0_r8)
        call elsi_say(eh%bh,msg)
   end if

   ! Write the inverse overlap here as test 
   ! write(msg,"(A)") "Writing S^-1 here as result."
   ! call elsi_say(eh%bh,msg)
   ! do i_k = 1,eh%bh%n_lrow
   !     write(*,*) "ELSI", ( tmp2(i_k,j), j=1,eh%bh%n_lcol ), eh%bh%myid
   ! end do

   ! use_scalapack_intern = .true.
   ! ! Try internal scalapack routine for inversion
   ! if (use_scalapack_intern) then
   !    tmp(:,:) = ovlp ! Temporary overlap
   !    call pdgetri(eh%ph%n_basis , tmp , 1 ,  1 ,eh%bh%desc ,IPVT, WORK ,  0  , IWORK , 0  , ierr )
   ! end if 

   ! Deallocate tmp arrays
   call elsi_deallocate(eh%bh,tmp,"tmp")
   call elsi_deallocate(eh%bh,tmp1,"tmp1")

   call elsi_get_time(t1) ! Timestamp

   write(msg,"(A)") "Finished inversion of overlap matrix."

   call elsi_say(eh%bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(eh%bh,msg)


end subroutine

subroutine elsi_get_inverse_ovlp_complex(eh,ovlp,inv_ovlp)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   complex(kind=r8), intent(in) :: ovlp(eh%bh%n_lrow,eh%bh%n_lcol) !< Overlap matrix
   complex(kind=r8), intent(out) :: inv_ovlp(eh%bh%n_lrow,eh%bh%n_lcol) !< Inverse Overlap

   character(len=*), parameter :: caller = "elsi_get_inverse_ovlp_complex"

   complex(kind=r8), allocatable :: tmp(:,:),tmp1(:,:)

   real(kind=r8) :: t0
   integer(kind=i4) :: ierr
   real(kind=r8) :: t1, norm
   complex(kind=r8) :: trace
   character(len=200) :: msg

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   call elsi_get_time(t0) ! Timestamp

   ! Allocate tmp arrays
   call elsi_allocate(eh%bh,tmp,eh%bh%n_lrow,eh%bh%n_lcol,"tmp",caller)
   call elsi_allocate(eh%bh,tmp1,eh%bh%n_lrow,eh%bh%n_lcol,"tmp1",caller)

   tmp(:,:) = ovlp ! Temporary overlap

   ! Maybe ill-conditioned check here with ?
   !call elsi_check_ovlp_elpa(eh%ph,eh%bh,tmp,eval_tmp,evec_tmp)

   ! Alternative to calculate the inverse of S is using scalapack routine psgetri

   ! ovlp = U^(-1)
   call elsi_factor_ovlp_elpa(eh%ph,eh%bh,tmp) ! ELPA routine to get inverse upper triangular
   tmp1(:,:) = tmp

   ! Get ovlp^(-1) = U^(-1)*(U^(-1))^T
   !tmp = matmul(tmp,conjg(transpose(tmp))), serial
   call pzgemm("N","C",eh%ph%n_basis,eh%ph%n_basis,eh%ph%n_basis,(1.0_r8,0.0_r8),tmp,1,&
        1,eh%bh%desc,tmp1,1,1,eh%bh%desc,(0.0_r8,0.0_r8),inv_ovlp,1,1,&
        eh%bh%desc)
   ! Try with specific routine, this should be checked in more detail
   ! tmp1 = U^(-T)
   call pztranc(eh%ph%n_basis,eh%ph%n_basis,(1.0_r8,0.0_r8),tmp,1,1,eh%bh%desc,(0.0_r8,0.0_r8),tmp1,1,1,&
           eh%bh%desc)
   tmp(:,:) = tmp1 
   ! tmp2 = U^(-1)*(U^(-1))^T
   call eh%ph%elpa_aux%hermitian_multiply("L","L",eh%ph%n_basis,tmp,tmp1,&
           eh%bh%n_lrow,eh%bh%n_lcol,inv_ovlp,eh%bh%n_lrow,eh%bh%n_lcol,ierr)
   call elsi_check_err(eh%bh,"ELPA matrix multiplication",ierr,caller)

   ! ovlp^(-1) should be symmetric in the end again, Several options to ensure that
   ! Option 1 is just setting one triangular half of the matrix
   call elsi_set_full_mat(eh%ph,eh%bh,LT_MAT,inv_ovlp) ! Full matrix
   ! Option 2 is symmetrizing by 0.5 * (ovlp + ovlp^T)
   ! tmp1(:,:) = tmp2
   ! call pzgeadd("C",eh%ph%n_basis,eh%ph%n_basis,(1.0_r8,0.0_r8),tmp1,1,1,eh%bh%desc,(1.0_r8,0.0_r8),tmp2,1,1,eh%bh%desc)
   ! tmp2(:,:) = (0.5_r8,0.0_r8) * tmp2

   ! Check if overlap times inverse is identity via the trace, Should implement check also for off-diagonal elements
   ! tmp = S * S^{-1}
   tmp = 0.0_r8
   call pzgemm("N","N",eh%ph%n_basis,eh%ph%n_basis,eh%ph%n_basis,(1.0_r8,0.0_r8),ovlp,1,&
        1,eh%bh%desc,inv_ovlp,1,1,eh%bh%desc,(0.0_r8,0.0_r8),tmp,1,1,&
        eh%bh%desc)
   ! Get Trace
   call elsi_trace_glob_dm_complex(eh,tmp,trace)
   ! Check if trace is equal to n_basis as check
   !write(msg,"(A,E10.3)") "***Overlap identity via trace:", (abs(real(trace))-eh%ph%n_basis)
   !call elsi_say(eh%bh,msg)
   if ((abs(real(trace))-eh%ph%n_basis)>1e-10_r8) then
        write(msg,"(A,E10.3)") "***Overlap is not correct. Difference to identity:", (abs(real(trace))-eh%ph%n_basis)
        call elsi_say(eh%bh,msg)
   end if
   ! Check RMS of norm to be 1
   call elsi_norm_residual_complex(eh,tmp,norm)
   ! This should be 1
   norm = (norm ** 2) * eh%ph%n_basis
   ! Check if trace is equal to n_basis as check
   !write(msg,"(A,E10.3)") "***Overlap identity via norm:", (norm-1.0_r8)
   !call elsi_say(eh%bh,msg)
   if ((norm-1.0_r8)>1e-10_r8) then
        write(msg,"(A,E10.3)") "***Overlap is not correct. Difference to identity via norm:", (norm-1.0_r8)
        call elsi_say(eh%bh,msg)
   end if

   ! Deallocate tmp arrays
   call elsi_deallocate(eh%bh,tmp,"tmp")
   call elsi_deallocate(eh%bh,tmp1,"tmp1")

   call elsi_get_time(t1) ! Timestamp

   write(msg,"(A)") "Finished inversion of overlap matrix."

   call elsi_say(eh%bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(eh%bh,msg)

end subroutine 

!>
!! Version two of the matrix inversion of S
!!
!!
subroutine elsi_get_inverse_ovlp_real_2(eh,ovlp,inv_ovlp)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: ovlp(eh%bh%n_lrow,eh%bh%n_lcol) !< Overlap matrix
   real(kind=r8), intent(out) :: inv_ovlp(eh%bh%n_lrow,eh%bh%n_lcol) !< contra-covariant Density matrix

   character(len=*), parameter :: caller = "elsi_get_inverse_ovlp_real_2"

   real(kind=r8), allocatable :: tmp(:,:), tmp1(:,:)

   integer :: i_k, j
   integer(kind=i4) :: ierr
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   real(kind=r8) :: trace, norm
   character(len=200) :: msg
   logical :: use_scalapack_intern

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   ! Allocate tmp arrays
   call elsi_allocate(eh%bh,tmp,eh%bh%n_lrow,eh%bh%n_lcol,"tmp",caller)
   call elsi_allocate(eh%bh,tmp1,eh%bh%n_lrow,eh%bh%n_lcol,"tmp1",caller)

   call elsi_get_time(t0) ! Timestamp

   ! Tmp overlap
   tmp(:,:) = ovlp
   tmp1(:,:) = 0.0_r8

   ! tmp = U^(-1) , still U matrix
   call elsi_factor_ovlp_elpa(eh%ph,eh%bh,tmp)

   ! tmp^(-1) = U^(-T) , this now a L matrix
   call pdtran(eh%ph%n_basis,eh%ph%n_basis,1.0_r8,tmp,1,1,eh%bh%desc,0.0_r8,tmp1,1,1,&
        eh%bh%desc)

   ! tmp = tmp^(-1) = U^(-T)
   tmp(:,:) = tmp1

   ! ovlp^(-1) = U^(-1) U^(-T)
   call eh%ph%elpa_aux%hermitian_multiply("L","full",eh%ph%n_basis,tmp,tmp1,eh%bh%n_lrow,&
        eh%bh%n_lcol,inv_ovlp,eh%bh%n_lrow,eh%bh%n_lcol,ierr)

   call elsi_check_err(eh%bh,"ELPA matrix multiplication",ierr,caller)

   ! Check if overlap times inverse is identity via the trace, Should implement check also for off-diagonal elements
   ! tmp = S * S^{-1}
   tmp = 0.0_r8
   call pdgemm("N","N",eh%ph%n_basis,eh%ph%n_basis,eh%ph%n_basis,1.0_r8,ovlp,1,&
        1,eh%bh%desc,inv_ovlp,1,1,eh%bh%desc,0.0_r8,tmp,1,1,&
        eh%bh%desc)
   ! Get Trace
   call elsi_trace_glob_dm_real(eh,tmp,trace)
   !Check if trace is equal to n_basis as check
   write(msg,"(A,E10.3)") "***Overlap identity via trace:", (abs(real(trace))-eh%ph%n_basis)
   call elsi_say(eh%bh,msg)
   if ((abs(trace)-eh%ph%n_basis)>1e-10_r8) then
        write(msg,"(A,E10.3)") "***Overlap is not correct. Difference to identity via trace:", (abs(trace)-eh%ph%n_basis)
        call elsi_say(eh%bh,msg)
   end if
   ! Check RMS of norm to be 1
   call elsi_norm_residual_real(eh,tmp,norm)
   ! This should be 1
   norm = (norm ** 2) * eh%ph%n_basis
   !Check if trace is equal to n_basis as check
   write(msg,"(A,E10.3)") "***Overlap identity via norm:", (norm-1.0_r8)
   call elsi_say(eh%bh,msg)
   if ((norm-1.0_r8)>1e-10_r8) then
        write(msg,"(A,E10.3)") "***Overlap is not correct. Difference to identity via norm:", (norm-1.0_r8)
        call elsi_say(eh%bh,msg)
   end if

   ! Deallocate tmp arrays
   call elsi_deallocate(eh%bh,tmp,"tmp")
   call elsi_deallocate(eh%bh,tmp1,"tmp1")

   call elsi_get_time(t1) ! Timestamp

   write(msg,"(A)") "Finished inversion of overlap matrix."

   call elsi_say(eh%bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(eh%bh,msg)

end subroutine

subroutine elsi_get_inverse_ovlp_complex_2(eh,ovlp,inv_ovlp)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   complex(kind=r8), intent(in) :: ovlp(eh%bh%n_lrow,eh%bh%n_lcol) !< Overlap matrix
   complex(kind=r8), intent(out) :: inv_ovlp(eh%bh%n_lrow,eh%bh%n_lcol) !< contra-covariant Density matrix

   character(len=*), parameter :: caller = "elsi_get_inverse_ovlp_complex_2"

   complex(kind=r8), allocatable :: tmp(:,:), tmp1(:,:)

   integer :: i_k, j
   integer(kind=i4) :: ierr
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   real(kind=r8) :: norm
   complex(kind=r8) :: trace
   character(len=200) :: msg
   logical :: use_scalapack_intern

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   ! Allocate tmp arrays
   call elsi_allocate(eh%bh,tmp,eh%bh%n_lrow,eh%bh%n_lcol,"tmp",caller)
   call elsi_allocate(eh%bh,tmp1,eh%bh%n_lrow,eh%bh%n_lcol,"tmp1",caller)

   call elsi_get_time(t0) ! Timestamp

   ! Tmp overlap
   tmp(:,:) = ovlp
   tmp1(:,:) = (0.0_r8,0.0_r8)

   ! tmp = U^(-1) , still U matrix
   call elsi_factor_ovlp_elpa(eh%ph,eh%bh,tmp)

   ! tmp^(-1) = U^(-T) , this now a L matrix
   call pztranc(eh%ph%n_basis,eh%ph%n_basis,(1.0_r8,0.0_r8),tmp,1,1,eh%bh%desc,(0.0_r8,0.0_r8),tmp1,1,1,&
        eh%bh%desc)

   ! tmp = tmp^(-1) = U^(-T)
   tmp(:,:) = tmp1

   ! ovlp^(-1) = U^(-1) U^(-T)
   call eh%ph%elpa_aux%hermitian_multiply("L","full",eh%ph%n_basis,tmp,tmp1,eh%bh%n_lrow,&
        eh%bh%n_lcol,inv_ovlp,eh%bh%n_lrow,eh%bh%n_lcol,ierr)

   call elsi_check_err(eh%bh,"ELPA matrix multiplication",ierr,caller)

   ! Check if overlap times inverse is identity via the trace, Should implement check also for off-diagonal elements
   ! tmp = S * S^{-1}
   tmp = (0.0_r8,0.0_r8)
   call pzgemm("N","N",eh%ph%n_basis,eh%ph%n_basis,eh%ph%n_basis,(1.0_r8,0.0_r8),ovlp,1,&
        1,eh%bh%desc,inv_ovlp,1,1,eh%bh%desc,(0.0_r8,0.0_r8),tmp,1,1,&
        eh%bh%desc)
   ! Get Trace
   call elsi_trace_glob_dm_complex(eh,tmp,trace)
   !Check if trace is equal to n_basis as check
   write(msg,"(A,E10.3)") "***Overlap identity via trace:", (abs(real(trace))-eh%ph%n_basis)
   call elsi_say(eh%bh,msg)
   if ((abs(real(trace))-eh%ph%n_basis)>1e-10_r8) then
        write(msg,"(A,E10.3)") "***Overlap is not correct. Difference to identity via trace:", (abs(real(trace))-eh%ph%n_basis)
        call elsi_say(eh%bh,msg)
   end if
   ! Check RMS of norm to be 1
   call elsi_norm_residual_complex(eh,tmp,norm)
   ! This should be 1
   norm = (norm ** 2) * eh%ph%n_basis
   !Check if trace is equal to n_basis as check
   write(msg,"(A,E10.3)") "***Overlap identity via norm:", (norm-1.0_r8)
   call elsi_say(eh%bh,msg)
   if ((norm-1.0_r8)>1e-10_r8) then
        write(msg,"(A,E10.3)") "***Overlap is not correct. Difference to identity via norm:", (norm-1.0_r8)
        call elsi_say(eh%bh,msg)
   end if

   ! Deallocate tmp arrays
   call elsi_deallocate(eh%bh,tmp,"tmp")
   call elsi_deallocate(eh%bh,tmp1,"tmp1")

   call elsi_get_time(t1) ! Timestamp

   write(msg,"(A)") "Finished inversion of overlap matrix."

   call elsi_say(eh%bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(eh%bh,msg)

end subroutine

!>
!! Calculate the response of the inverse overlap dS^(-1)/dR from the already calculated dS/dR
!!
!!
subroutine elsi_inv_ovlp_response_real(eh,ovlp_response,inv_ovlp,inv_ovlp_response)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: ovlp_response(eh%bh%n_lrow,eh%bh%n_lcol) !< ovlp response dS/dR
   real(kind=r8), intent(in) :: inv_ovlp(eh%bh%n_lrow,eh%bh%n_lcol) !< inverse ovlp
   real(kind=r8), intent(out) :: inv_ovlp_response(eh%bh%n_lrow,eh%bh%n_lcol) !< inv ovlp response dS^(-1)/dR

   character(len=*), parameter :: caller = "elsi_inv_ovlp_response_real"

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   real(kind=r8) :: factor

   real(kind=r8), allocatable :: tmp(:,:)

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   call elsi_allocate(eh%bh,tmp,eh%bh%n_lrow,eh%bh%n_lcol,"tmp",caller)

   factor = -1.0_r8 ! Needed to scale the inv ovlp at the end.

   call elsi_get_time(t0) ! Timestamp

   ! Do dS^(-1)/dR = - inv_ovlp * dS/dR * inv_ovlp
   call pdgemm("N","N",eh%ph%n_basis,eh%ph%n_basis,eh%ph%n_basis,1.0_r8,inv_ovlp,1,&
        1,eh%bh%desc,ovlp_response,1,1,eh%bh%desc,0.0_r8,tmp,1,1,&
        eh%bh%desc)

   call pdgemm("N","N",eh%ph%n_basis,eh%ph%n_basis,eh%ph%n_basis,1.0_r8,tmp,1,&
        1,eh%bh%desc,inv_ovlp,1,1,eh%bh%desc,0.0_r8,inv_ovlp_response,1,1,&
        eh%bh%desc)

   ! Scale the inv_ovlp_response at the end with factor -1
   inv_ovlp_response = factor * inv_ovlp_response

   call elsi_get_time(t1) ! Timestamp

   write(msg,"(A)") "Finished dS^(-1)/dR from dS/dR calculation."

   call elsi_say(eh%bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(eh%bh,msg)

   call elsi_deallocate(eh%bh,tmp,"tmp")

end subroutine

subroutine elsi_inv_ovlp_response_complex(eh,ovlp_response,inv_ovlp,inv_ovlp_response)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   complex(kind=r8), intent(in) :: ovlp_response(eh%bh%n_lrow,eh%bh%n_lcol) !< ovlp response dS/dR
   complex(kind=r8), intent(in) :: inv_ovlp(eh%bh%n_lrow,eh%bh%n_lcol) !< inverse ovlp
   complex(kind=r8), intent(out) :: inv_ovlp_response(eh%bh%n_lrow,eh%bh%n_lcol) !< inv ovlp response dS^(-1)/dR

   character(len=*), parameter :: caller = "elsi_inv_ovlp_response_real"

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   complex(kind=r8) :: factor

   complex(kind=r8), allocatable :: tmp(:,:)

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   call elsi_allocate(eh%bh,tmp,eh%bh%n_lrow,eh%bh%n_lcol,"tmp",caller)

   factor = -(1.0_r8,0.0_r8) ! Needed to scale the inv ovlp at the end.

   call elsi_get_time(t0) ! Timestamp

   ! Do dS^(-1)/dR = - inv_ovlp * dS/dR * inv_ovlp
   call pzgemm("N","N",eh%ph%n_basis,eh%ph%n_basis,eh%ph%n_basis,(1.0_r8,0.0_r8),inv_ovlp,1,&
        1,eh%bh%desc,ovlp_response,1,1,eh%bh%desc,(0.0_r8,0.0_r8),tmp,1,1,&
        eh%bh%desc)

   call pzgemm("N","N",eh%ph%n_basis,eh%ph%n_basis,eh%ph%n_basis,(1.0_r8,0.0_r8),tmp,1,&
        1,eh%bh%desc,inv_ovlp,1,1,eh%bh%desc,(0.0_r8,0.0_r8),inv_ovlp_response,1,1,&
        eh%bh%desc)

   ! Scale the inv_ovlp_response at the end with factor -1
   inv_ovlp_response = factor * inv_ovlp_response

   call elsi_get_time(t1) ! Timestamp

   write(msg,"(A)") "Finished dS^(-1)/dR from dS/dR calculation."

   call elsi_say(eh%bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(eh%bh,msg)

   call elsi_deallocate(eh%bh,tmp,"tmp")

end subroutine

!>
!! Get the density matrix from contra-covariant density matrix. 
!! dm = ccdm * ovlp^(-1)
!! DO I NEED AN ILL_CONDITIONED CHECK HERE?
!!
subroutine elsi_get_dm_from_ccdm_real(eh,ovlp,ccdm,dm)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: ovlp(eh%bh%n_lrow,eh%bh%n_lcol) !< Overlap matrix
   real(kind=r8), intent(in) :: ccdm(eh%bh%n_lrow,eh%bh%n_lcol) !< contra-covariant Density matrix
   real(kind=r8), intent(out) :: dm(eh%bh%n_lrow,eh%bh%n_lcol) !<  New Density matrix

   character(len=*), parameter :: caller = "elsi_get_dm_from_ccdm_real"

   real(kind=r8), allocatable :: tmp(:,:)

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   call elsi_get_time(t0) ! Timestamp

   ! Allocate tmp arrays
   call elsi_allocate(eh%bh,tmp,eh%bh%n_lrow,eh%bh%n_lcol,"tmp",caller)

   ! Get inverse overlap and store in tmp
   call elsi_get_inverse_ovlp_real_2(eh,ovlp,tmp)

   ! dm = ccdm * ovlp
   !dm = matmul(ccdm,tmp) ! Simple matrix-matrix multiplication for now
   call pdgemm("N","N",eh%ph%n_basis,eh%ph%n_basis,eh%ph%n_basis,1.0_r8,ccdm,1,&
        1,eh%bh%desc,tmp,1,1,eh%bh%desc,0.0_r8,dm,1,1,&
        eh%bh%desc)

   ! Deallocate tmp arrays
   call elsi_deallocate(eh%bh,tmp,"tmp")

   call elsi_get_time(t1) ! Timestamp

   write(msg,"(A)") "Finished density matrix calculation from ccdm."

   call elsi_say(eh%bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(eh%bh,msg)

end subroutine

! subroutine elsi_get_dm_from_ccdm_real(eh,ovlp,ccdm,dm)

!    implicit none

!    type(elsi_handle), intent(inout) :: eh !< Handle
!    real(kind=r8), intent(in) :: ovlp(eh%bh%n_lrow,eh%bh%n_lcol) !< Overlap matrix
!    real(kind=r8), intent(in) :: ccdm(eh%bh%n_lrow,eh%bh%n_lcol) !< contra-covariant Density matrix
!    real(kind=r8), intent(out) :: dm(eh%bh%n_lrow,eh%bh%n_lcol) !<  New Density matrix

!    character(len=*), parameter :: caller = "elsi_get_dm_from_ccdm_real"

!    real(kind=r8), allocatable :: tmp(:,:), tmp1(:,:), tmp2(:,:)

!    integer :: i_k, j
!    integer(kind=i4) :: ierr
!    real(kind=r8) :: t0
!    real(kind=r8) :: t1
!    real(kind=r8) :: trace, norm
!    character(len=200) :: msg

!    call elsi_check_init(eh%bh,eh%handle_init,caller)

!    call elsi_get_time(t0) ! Timestamp

!    ! Allocate tmp arrays
!    call elsi_allocate(eh%bh,tmp,eh%bh%n_lrow,eh%bh%n_lcol,"tmp",caller)
!    call elsi_allocate(eh%bh,tmp1,eh%bh%n_lrow,eh%bh%n_lcol,"tmp1",caller)
!    call elsi_allocate(eh%bh,tmp2,eh%bh%n_lrow,eh%bh%n_lcol,"tmp2",caller)

!    tmp(:,:) = ovlp ! Temporary overlap

!    ! write(msg,"(A)") "Writing S here as result."
!    ! call elsi_say(eh%bh,msg)
!    ! do i_k = 1,eh%bh%n_lrow
!    !     write(*,*) "ELSI", ( tmp(i_k,j), j=1,eh%bh%n_lcol ), eh%bh%myid
!    ! end do

!    ! Maybe ill-conditioned check here with ?
!    !call elsi_check_ovlp_elpa(eh%ph,eh%bh,tmp,eval_tmp,evec_tmp)

!    ! ovlp = U^(-1)
!    call elsi_factor_ovlp_elpa(eh%ph,eh%bh,tmp) ! ELPA routine to get inverse
!    tmp1(:,:) = tmp

!    ! Get ovlp^(-1) = U^(-1)*(U^(-1))^T
!    ! Here we are multiplying an upper- times a lower-trianglular matrix. Instead of pdgemm it might be good to use
!    ! specific routines for upper, lower matrices instead
!    !tmp = matmul(tmp,transpose(tmp)) ! Simple matrix-matrix multiplication for now
!    ! call pdgemm("N","T",eh%ph%n_basis,eh%ph%n_basis,eh%ph%n_basis,1.0_r8,tmp,1,&
!    !      1,eh%bh%desc,tmp1,1,1,eh%bh%desc,0.0_r8,tmp2,1,1,&
!    !      eh%bh%desc)
!    ! Try with specific routine, this should be checked in more detail
!    ! tmp1 = U^(-T)
!    call pdtran(eh%ph%n_basis,eh%ph%n_basis,1.0_r8,tmp,1,1,eh%bh%desc,0.0_r8,tmp1,1,1,&
!            eh%bh%desc)
!    tmp(:,:) = tmp1 
!    ! tmp2 = U^(-1)*(U^(-1))^T
!    call eh%ph%elpa_aux%hermitian_multiply("L","L",eh%ph%n_basis,tmp,tmp1,&
!            eh%bh%n_lrow,eh%bh%n_lcol,tmp2,eh%bh%n_lrow,eh%bh%n_lcol,ierr)
!    call elsi_check_err(eh%bh,"ELPA matrix multiplication",ierr,caller)

!    ! ovlp^(-1) should be symmetric in the end again, Several options to ensure that
!    ! Option 1 is just setting one triangular half of the matrix
!    call elsi_set_full_mat(eh%ph,eh%bh,LT_MAT,tmp2) ! Full matrix
!    ! Option 2 is symmetrizing by 0.5 * (ovlp + ovlp^T)
!    ! tmp1(:,:) = tmp2
!    ! call pdgeadd("T",eh%ph%n_basis,eh%ph%n_basis,1.0_r8,tmp1,1,1,eh%bh%desc,1.0_r8,tmp2,1,1,eh%bh%desc)
!    ! tmp2(:,:) = 0.5_r8 * tmp2

!    ! Check if overlap times inverse is identity via the trace, Should implement check also for off-diagonal elements
!    ! tmp = S * S^{-1}
!    tmp = 0.0_r8
!    call pdgemm("N","N",eh%ph%n_basis,eh%ph%n_basis,eh%ph%n_basis,1.0_r8,ovlp,1,&
!         1,eh%bh%desc,tmp2,1,1,eh%bh%desc,0.0_r8,tmp,1,1,&
!         eh%bh%desc)
!    ! Get Trace
!    call elsi_trace_glob_dm_real(eh,tmp,trace)
!    !Check if trace is equal to n_basis as check
!    write(msg,"(A,E10.3)") "***Overlap identity via trace:", (abs(real(trace))-eh%ph%n_basis)
!    call elsi_say(eh%bh,msg)
!    if ((abs(trace)-eh%ph%n_basis)>1e-10_r8) then
!         write(msg,"(A,E10.3)") "***Overlap is not correct. Difference to identity via trace:", (abs(trace)-eh%ph%n_basis)
!         call elsi_say(eh%bh,msg)
!    end if
!    ! Check RMS of norm to be 1
!    call elsi_norm_residual_real(eh,tmp,norm)
!    ! This should be 1
!    norm = (norm ** 2) * eh%ph%n_basis
!    !Check if trace is equal to n_basis as check
!    write(msg,"(A,E10.3)") "***Overlap identity via norm:", (norm-1.0_r8)
!    call elsi_say(eh%bh,msg)
!    if ((norm-1.0_r8)>1e-10_r8) then
!         write(msg,"(A,E10.3)") "***Overlap is not correct. Difference to identity via norm:", (norm-1.0_r8)
!         call elsi_say(eh%bh,msg)
!    end if

!    ! Write the inverse overlap here as test 
!    ! write(msg,"(A)") "Writing S^-1 here as result."
!    ! call elsi_say(eh%bh,msg)
!    ! do i_k = 1,eh%bh%n_lrow
!    !     write(*,*) "ELSI", ( tmp2(i_k,j), j=1,eh%bh%n_lcol ), eh%bh%myid
!    ! end do

!    ! dm = ccdm * ovlp
!    !dm = matmul(ccdm,tmp) ! Simple matrix-matrix multiplication for now
!    call pdgemm("N","N",eh%ph%n_basis,eh%ph%n_basis,eh%ph%n_basis,1.0_r8,ccdm,1,&
!         1,eh%bh%desc,tmp2,1,1,eh%bh%desc,0.0_r8,dm,1,1,&
!         eh%bh%desc)

!    ! Deallocate tmp arrays
!    call elsi_deallocate(eh%bh,tmp,"tmp")
!    call elsi_deallocate(eh%bh,tmp1,"tmp1")
!    call elsi_deallocate(eh%bh,tmp2,"tmp2")

!    call elsi_get_time(t1) ! Timestamp

!    write(msg,"(A)") "Finished density matrix calculation from ccdm."

!    call elsi_say(eh%bh,msg)
!    write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
!    call elsi_say(eh%bh,msg)

! end subroutine


!>
!! Get the density matrix from contra-covariant density matrix. 
!! dm = ccdm * ovlp^(-1)
!! DO I NEED AN ILL_CONDITIONED CHECK HERE?
!!
subroutine elsi_get_dm_from_ccdm_complex(eh,ovlp,ccdm,dm)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   complex(kind=r8), intent(in) :: ovlp(eh%bh%n_lrow,eh%bh%n_lcol) !< Overlap matrix
   complex(kind=r8), intent(in) :: ccdm(eh%bh%n_lrow,eh%bh%n_lcol) !< contra-covariant Density matrix
   complex(kind=r8), intent(out) :: dm(eh%bh%n_lrow,eh%bh%n_lcol) !<  New Density matrix

   character(len=*), parameter :: caller = "elsi_get_dm_from_ccdm_complex"

   complex(kind=r8), allocatable :: tmp(:,:)

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   call elsi_get_time(t0) ! Timestamp

   ! Allocate tmp arrays
   call elsi_allocate(eh%bh,tmp,eh%bh%n_lrow,eh%bh%n_lcol,"tmp",caller)

   ! Get inverse overlap and store in tmp
   call elsi_get_inverse_ovlp_complex_2(eh,ovlp,tmp)

   ! dm = ccdm * ovlp^-1
   !dm = matmul(ccdm,tmp) ! Simple matrix-matrix multiplication for now, serial
   call pzgemm("N","N",eh%ph%n_basis,eh%ph%n_basis,eh%ph%n_basis,(1.0_r8,0.0_r8),ccdm,1,&
        1,eh%bh%desc,tmp,1,1,eh%bh%desc,(0.0_r8,0.0_r8),dm,1,1,&
        eh%bh%desc)

   ! Deallocate tmp arrays
   call elsi_deallocate(eh%bh,tmp,"tmp")

   call elsi_get_time(t1) ! Timestamp

   write(msg,"(A)") "Finished density matrix calculation from ccdm."

   call elsi_say(eh%bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(eh%bh,msg)

end subroutine

! subroutine elsi_get_dm_from_ccdm_complex(eh,ovlp,ccdm,dm)

!    implicit none

!    type(elsi_handle), intent(inout) :: eh !< Handle
!    complex(kind=r8), intent(in) :: ovlp(eh%bh%n_lrow,eh%bh%n_lcol) !< Overlap matrix
!    complex(kind=r8), intent(in) :: ccdm(eh%bh%n_lrow,eh%bh%n_lcol) !< contra-covariant Density matrix
!    complex(kind=r8), intent(out) :: dm(eh%bh%n_lrow,eh%bh%n_lcol) !<  New Density matrix

!    character(len=*), parameter :: caller = "elsi_get_dm_from_ccdm_complex"

!    complex(kind=r8), allocatable :: tmp(:,:),tmp1(:,:), tmp2(:,:)

!    real(kind=r8) :: t0
!    integer(kind=i4) :: ierr
!    real(kind=r8) :: t1, norm
!    complex(kind=r8) :: trace
!    character(len=200) :: msg

!    call elsi_check_init(eh%bh,eh%handle_init,caller)

!    call elsi_get_time(t0) ! Timestamp

!    ! Allocate tmp arrays
!    call elsi_allocate(eh%bh,tmp,eh%bh%n_lrow,eh%bh%n_lcol,"tmp",caller)
!    call elsi_allocate(eh%bh,tmp1,eh%bh%n_lrow,eh%bh%n_lcol,"tmp1",caller)
!    call elsi_allocate(eh%bh,tmp2,eh%bh%n_lrow,eh%bh%n_lcol,"tmp2",caller)

!    tmp(:,:) = ovlp ! Temporary overlap

!    ! Maybe ill-conditioned check here with ?
!    !call elsi_check_ovlp_elpa(eh%ph,eh%bh,tmp,eval_tmp,evec_tmp)

!    ! Alternative to calculate the inverse of S is using scalapack routine psgetri

!    ! ovlp = U^(-1)
!    call elsi_factor_ovlp_elpa(eh%ph,eh%bh,tmp) ! ELPA routine to get inverse
!    tmp1(:,:) = tmp

!    ! Get ovlp^(-1) = U^(-1)*(U^(-1))^T
!    !tmp = matmul(tmp,conjg(transpose(tmp))), serial
!    call pzgemm("N","C",eh%ph%n_basis,eh%ph%n_basis,eh%ph%n_basis,(1.0_r8,0.0_r8),tmp,1,&
!         1,eh%bh%desc,tmp1,1,1,eh%bh%desc,(0.0_r8,0.0_r8),tmp2,1,1,&
!         eh%bh%desc)
!    ! Try with specific routine, this should be checked in more detail
!    ! tmp1 = U^(-T)
!    call pztranc(eh%ph%n_basis,eh%ph%n_basis,(1.0_r8,0.0_r8),tmp,1,1,eh%bh%desc,(0.0_r8,0.0_r8),tmp1,1,1,&
!            eh%bh%desc)
!    tmp(:,:) = tmp1 
!    ! tmp2 = U^(-1)*(U^(-1))^T
!    call eh%ph%elpa_aux%hermitian_multiply("L","L",eh%ph%n_basis,tmp,tmp1,&
!            eh%bh%n_lrow,eh%bh%n_lcol,tmp2,eh%bh%n_lrow,eh%bh%n_lcol,ierr)
!    call elsi_check_err(eh%bh,"ELPA matrix multiplication",ierr,caller)

!    ! ovlp^(-1) should be symmetric in the end again, Several options to ensure that
!    ! Option 1 is just setting one triangular half of the matrix
!    call elsi_set_full_mat(eh%ph,eh%bh,LT_MAT,tmp2) ! Full matrix
!    ! Option 2 is symmetrizing by 0.5 * (ovlp + ovlp^T)
!    ! tmp1(:,:) = tmp2
!    ! call pzgeadd("C",eh%ph%n_basis,eh%ph%n_basis,(1.0_r8,0.0_r8),tmp1,1,1,eh%bh%desc,(1.0_r8,0.0_r8),tmp2,1,1,eh%bh%desc)
!    ! tmp2(:,:) = (0.5_r8,0.0_r8) * tmp2

!    ! Check if overlap times inverse is identity via the trace, Should implement check also for off-diagonal elements
!    ! tmp = S * S^{-1}
!    tmp = 0.0_r8
!    call pzgemm("N","N",eh%ph%n_basis,eh%ph%n_basis,eh%ph%n_basis,(1.0_r8,0.0_r8),ovlp,1,&
!         1,eh%bh%desc,tmp2,1,1,eh%bh%desc,(0.0_r8,0.0_r8),tmp,1,1,&
!         eh%bh%desc)
!    ! Get Trace
!    call elsi_trace_glob_dm_complex(eh,tmp,trace)
!    ! Check if trace is equal to n_basis as check
!    !write(msg,"(A,E10.3)") "***Overlap identity via trace:", (abs(real(trace))-eh%ph%n_basis)
!    !call elsi_say(eh%bh,msg)
!    if ((abs(real(trace))-eh%ph%n_basis)>1e-13_r8) then
!         write(msg,"(A,E10.3)") "***Overlap is not correct. Difference to identity:", (abs(real(trace))-eh%ph%n_basis)
!         call elsi_say(eh%bh,msg)
!    end if
!    ! Check RMS of norm to be 1
!    call elsi_norm_residual_complex(eh,tmp,norm)
!    ! This should be 1
!    norm = (norm ** 2) * eh%ph%n_basis
!    ! Check if trace is equal to n_basis as check
!    !write(msg,"(A,E10.3)") "***Overlap identity via norm:", (norm-1.0_r8)
!    !call elsi_say(eh%bh,msg)
!    if ((norm-1.0_r8)>1e-13_r8) then
!         write(msg,"(A,E10.3)") "***Overlap is not correct. Difference to identity via norm:", (norm-1.0_r8)
!         call elsi_say(eh%bh,msg)
!    end if

!    ! dm = ccdm * ovlp
!    !dm = matmul(ccdm,tmp) ! Simple matrix-matrix multiplication for now, serial
!    call pzgemm("N","N",eh%ph%n_basis,eh%ph%n_basis,eh%ph%n_basis,(1.0_r8,0.0_r8),ccdm,1,&
!         1,eh%bh%desc,tmp2,1,1,eh%bh%desc,(0.0_r8,0.0_r8),dm,1,1,&
!         eh%bh%desc)

!    ! Deallocate tmp arrays
!    call elsi_deallocate(eh%bh,tmp,"tmp")
!    call elsi_deallocate(eh%bh,tmp1,"tmp1")
!    call elsi_deallocate(eh%bh,tmp2,"tmp2")

!    call elsi_get_time(t1) ! Timestamp

!    write(msg,"(A)") "Finished density matrix calculation from ccdm."

!    call elsi_say(eh%bh,msg)
!    write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
!    call elsi_say(eh%bh,msg)

! end subroutine

!>
!! Get the residual between ccdm (D*S) and X. The resulting matrix
!! is called the residual.
!! res = D*S - X
!!
subroutine elsi_residual_ccdm_real(eh,ccdm_current,X_current,res)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: ccdm_current(eh%bh%n_lrow,eh%bh%n_lcol) !< contra-covariant density matrix in current step
   real(kind=r8), intent(in) :: X_current(eh%bh%n_lrow,eh%bh%n_lcol) !< contra-covariant density matrix in previous step
   real(kind=r8), intent(out) :: res(eh%bh%n_lrow,eh%bh%n_lcol) !< residual

   character(len=*), parameter :: caller = "elsi_residual_ccdm_real"

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   integer :: i_k, j

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   call elsi_get_time(t0) ! Timestamp

   res = ccdm_current - X_current ! For now simple fortran matrix-matrix subtraction

   call elsi_get_time(t1) ! Timestamp

   write(msg,"(A)") "Finished residual calculation."

   call elsi_say(eh%bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(eh%bh,msg)

end subroutine

!>
!! Get the residual between ccdm of two steps. The resulting matrix
!! is called the residual.
!! res = D*S - X
!!
subroutine elsi_residual_ccdm_complex(eh,ccdm_current,X_current,res)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   complex(kind=r8), intent(in) :: ccdm_current(eh%bh%n_lrow,eh%bh%n_lcol) !< contra-covariant density matrix in current step
   complex(kind=r8), intent(in) :: X_current(eh%bh%n_lrow,eh%bh%n_lcol) !< contra-covariant density matrix in previous step
   complex(kind=r8), intent(out) :: res(eh%bh%n_lrow,eh%bh%n_lcol) !< residual

   character(len=*), parameter :: caller = "elsi_residual_ccdm_complex"

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   call elsi_get_time(t0) ! Timestamp

   res = ccdm_current - X_current ! For now simple fortran matrix-matrix subtraction

   call elsi_get_time(t1) ! Timestamp

   write(msg,"(A)") "Finished residual calculation."

   call elsi_say(eh%bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(eh%bh,msg)

end subroutine

!>
!! Here come the specific Kernel-based routines that evaluate the second time-derivative
!! ddot(X).
!!
!!

!>
!! Get the second time derivative, accelaration (ccdm_td) of the X in the scaled delta approximation.
!! 
!! dX/dt^2 = c * (DS-X) , c is real number in intervall [0,1]
!!
!!
subroutine elsi_ccdm_scaled_delta_approx_real(eh,res,c,X_td)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: res(eh%bh%n_lrow,eh%bh%n_lcol) !< residual in current step
   real(kind=r8), intent(in) :: c !< real number in [0,1]
   real(kind=r8), intent(out) :: X_td(eh%bh%n_lrow,eh%bh%n_lcol) !< second time deriv. of ccdm

   character(len=*), parameter :: caller = "elsi_ccdm_scaled_delta_approx_real"

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   call elsi_get_time(t0) ! Timestamp

   X_td = c * res! For now simple fortran matrix-scalar multiplication

   call elsi_get_time(t1) ! Timestamp

   write(msg,"(A)") "Finished second time derivative of X calculation in the scaled delta approximation."

   call elsi_say(eh%bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(eh%bh,msg)

end subroutine

!>
!! Get the second time derivative (ccdm_td), accelaration of the X in the scaled delta approximation.
!! 
!! dX/dt^2 = -c * (DS-X) , c is real number in intervall [0,1]
!!
!!
subroutine elsi_ccdm_scaled_delta_approx_complex(eh,res,c,X_td)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   complex(kind=r8), intent(in) :: res(eh%bh%n_lrow,eh%bh%n_lcol) !< residual in current step
   real(kind=r8), intent(in) :: c !< real number in [0,1]
   complex(kind=r8), intent(out) :: X_td(eh%bh%n_lrow,eh%bh%n_lcol) !< second time deriv. of ccdm

   character(len=*), parameter :: caller = "elsi_ccdm_scaled_delta_approx_complex"

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   call elsi_get_time(t0) ! Timestamp

   X_td = c * res! For now simple fortran matrix-scalar multiplication

   call elsi_get_time(t1) ! Timestamp

   write(msg,"(A)") "Finished second time derivative of X calculation in the scaled delta approximation."

   call elsi_say(eh%bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(eh%bh,msg)

end subroutine

!>
!! Routines related to the dissipation of the ccdm.
!!
!!

!>
!! This routine initializes all the relevant variables for the dissipation in the XL-BOMD run.
!! Only needs the number of dissipative terms as input.
!!
!!
subroutine elsi_init_ccdm_dissipation (eh, MD_tstep, k_max, c_k, kappa_XL, alpha_XL, omega_XL)

        implicit none

        type(elsi_handle), intent(inout) :: eh !< Handle
        real(kind=r8), intent(in) :: MD_tstep !< MD timestep 
        integer(kind=i4), intent(in) :: k_max !< Maximum # of past ccdm used in dissipation
        real(kind=r8), intent(out) :: c_k(k_max) !< Constants c_k in the dissipation
        real(kind=r8), intent(out) :: kappa_XL !< Constant kappa
        real(kind=r8), intent(out) :: alpha_XL !< Constant alpha
        real(kind=r8), intent(out) :: omega_XL !< Constant omega = kappa / MD_ts**2

        character(len=*), parameter :: caller = "elsi_init_ccdm_dissipation"

        character(len=200) :: msg

        ! Local variables
        integer :: K_iter

        ! Set the maximum number of dissipative terms
        K_iter = k_max - 1

        ! Define the constants c_K and alpha based on K_iter, MD_tstep is timestep in MD
        ! From Niklasson et al., J. Chem. Phys. 130, 214109 (2009)
        select case(K_iter)
        case (3)
            kappa_XL = 1.69_r8
            alpha_XL = 0.150_r8
            c_K(1) = -2.0_r8
            c_K(2) = 3.0_r8
            c_K(3) = 0.0_r8
            c_K(4) = -1.0_r8
        case (4)
            kappa_XL = 1.75_r8
            alpha_XL = 0.057_r8
            c_K(1) = -3.0_r8
            c_K(2) = 6.0_r8
            c_K(3) = -2.0_r8
            c_K(4) = -2.0_r8
            c_K(5) = 1.0_r8
        case (5)
            kappa_XL = 1.82_r8
            alpha_XL = 0.018_r8
            c_K(1) = -6.0_r8
            c_K(2) = 14.0_r8
            c_K(3) = -8.0_r8
            c_K(4) = -3.0_r8
            c_K(5) = 4.0_r8
            c_K(6) = -1.0_r8
        case (6)
            kappa_XL = 1.84_r8
            alpha_XL = 0.0055_r8
            c_K(1) = -14.0_r8
            c_K(2) = 36.0_r8
            c_K(3) = -27.0_r8
            c_K(4) = -2.0_r8
            c_K(5) = 12.0_r8
            c_K(6) = -6.0_r8
            c_K(7) = 1.0_r8
        case (7)
            kappa_XL = 1.86_r8
            alpha_XL = 0.0016_r8
            c_K(1) = -36.0_r8
            c_K(2) = 99.0_r8
            c_K(3) = -88.0_r8
            c_K(4) = 11.0_r8
            c_K(5) = 32.0_r8
            c_K(6) = -25.0_r8
            c_K(7) = 8.0_r8
            c_K(8) = -1.0_r8
        case default
            write(msg,"(A)") "The specified number of dissipation terms are not available. Select a number of 3 to 7."
            call elsi_say(eh%bh,msg)
        end select

        ! Get the other constant omega
        if (MD_tstep==0.0_r8) then
            omega_XL = 1e200_r8
        else
            omega_XL = kappa_XL / MD_tstep**2
        end if

        write(msg,"(A)") "The relevant constants for XL-BOMD have been set."

        call elsi_say(eh%bh,msg)
        write(msg,"(A,I5)") "| Number of dissipative iterations # :",K_iter
        call elsi_say(eh%bh,msg)


    end subroutine

!>
!! Get the second time derivative (ccdm_td) of the ccdm in the scaled delta approximation.
!! 
!! ccdm_dissipation = alpha*sum_k (c_k * X(t-k*delta t)) , alpha, c_k are constants
!!
!!
subroutine elsi_ccdm_dissipation_real(eh,k_max,alpha,c_k,X_instances,X_dissipation)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: k_max !< Maximum # of past X used in dissipation
   real(kind=r8), intent(in) :: alpha !< Constant
   real(kind=r8), intent(in) :: c_k(k_max) !< Constant
   real(kind=r8), intent(in) :: X_instances(k_max,eh%bh%n_lrow,eh%bh%n_lcol) !< multiple instances of X
   real(kind=r8), intent(out) :: X_dissipation(eh%bh%n_lrow,eh%bh%n_lcol) !< dissipation of c_k

   character(len=*), parameter :: caller = "elsi_ccdm_dissipation_real"

   real(kind=r8) :: t0, t1
   character(len=200) :: msg
   integer(kind=i4) :: j

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   call elsi_get_time(t0) ! Timestamp

   X_dissipation = 0.0_r8

   ! Sum over the ccdm entries
   do j=1,k_max
      X_dissipation(:,:) = X_dissipation(:,:) + c_k(j) * X_instances(j,:,:)
   end do

   ! Scaling by alpha
   X_dissipation = X_dissipation * alpha

   call elsi_get_time(t1) ! Timestamp

   write(msg,"(A)") "Finished calculating X dissipation term."

   call elsi_say(eh%bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(eh%bh,msg)

end subroutine

!>
!! Get the second time derivative (ccdm_td) of the ccdm in the scaled delta approximation.
!! 
!! ccdm_dissipation = alpha*sum_k (c_k * X(t-k*delta t)) , alpha, c_k are constants
!!
!!
subroutine elsi_ccdm_dissipation_complex(eh,k_max,alpha,c_k,X_instances,X_dissipation)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: k_max !< Maximum # of past X used in dissipation
   real(kind=r8), intent(in) :: alpha !< Constant
   real(kind=r8), intent(in) :: c_k(k_max) !< Constant
   complex(kind=r8), intent(in) :: X_instances(k_max,eh%bh%n_lrow,eh%bh%n_lcol) !< multiple instances of X
   complex(kind=r8), intent(out) :: X_dissipation(eh%bh%n_lrow,eh%bh%n_lcol) !< dissipation of c_k

   character(len=*), parameter :: caller = "elsi_ccdm_dissipation_complex"

   real(kind=r8) :: t0, t1
   character(len=200) :: msg
   integer(kind=i4) :: j

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   call elsi_get_time(t0) ! Timestamp

   X_dissipation = (0.0_r8,0.0_r8)

   ! Sum over the ccdm entries
   do j=1,k_max
      X_dissipation(:,:) = X_dissipation(:,:) + c_k(j) * X_instances(j,:,:)
   end do
   ! Scaling by alpha
   X_dissipation = X_dissipation * alpha

   call elsi_get_time(t1) ! Timestamp

   write(msg,"(A)") "Finished calculating X dissipation term."

   call elsi_say(eh%bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(eh%bh,msg)

end subroutine

!>
!! Routine to integrate the Contra-covariant DM X as specified in the 
!! XL-BOMD approach.
!!
!! X(t+delta t) = 
!!   2*X(t) - X(t-delta t) + kappa * dX(t)/dt^2 + 
!!   alpha*sum_k (c_k * X(t-k*delta t))
!!   = 2*X(t) - X(t-delta t) + kappa * integrated_X + dissipative_X
!!
!!
subroutine elsi_ccdm_integration_real(eh,current_index,k_max,md_timestep,ccdm_current,X_array,X_integrated,residual)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: current_index !< Current index / MD step, needed to define the current and previous ccdm
   integer(kind=i4), intent(in) :: k_max !<  # of saved ccdm
   real(kind=r8), intent(in) :: md_timestep !< MD timestep
   real(kind=r8), intent(in) :: ccdm_current(eh%bh%n_lrow,eh%bh%n_lcol) !< multiple instances of ccdm
   real(kind=r8), intent(in) :: X_array(k_max,eh%bh%n_lrow,eh%bh%n_lcol) !< multiple instances of ccdm
   real(kind=r8), intent(out) :: X_integrated(eh%bh%n_lrow,eh%bh%n_lcol) !< integrated ccdm
   real(kind=r8), intent(out) :: residual(eh%bh%n_lrow,eh%bh%n_lcol) !< residual of current and previous step

   character(len=*), parameter :: caller = "elsi_ccdm_integration_real"

   real(kind=r8) :: t0, t1
   real(kind=r8) :: alpha, omega, kappa
   character(len=200) :: msg
   real(kind=r8) :: c_k(k_max)
   integer(kind=i4) :: previous_index

   real(kind=r8), allocatable :: diss_tmp(:,:)
   real(kind=r8), allocatable :: td_tmp(:,:)

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   call elsi_allocate(eh%bh,diss_tmp,eh%bh%n_lrow,eh%bh%n_lcol,"diss_tmp",caller)
   call elsi_allocate(eh%bh,td_tmp,eh%bh%n_lrow,eh%bh%n_lcol,"td_tmp",caller)

   call elsi_get_time(t0) ! Timestamp

   X_integrated = 0.0_r8

   ! Init the relevant constants.
   call elsi_init_ccdm_dissipation(eh,md_timestep,k_max,c_k,kappa,alpha,omega)

   ! Get previous array index based on current one
    previous_index = current_index + 1

   write(msg,"(A,I2,A)") "Current index :",current_index
   call elsi_say(eh%bh,msg)
   write(msg,"(A,I2,A)") "Previous index :",previous_index
   call elsi_say(eh%bh,msg)


   ! Get residual
   call elsi_residual_ccdm_real(eh,ccdm_current,X_array(current_index,:,:),residual)

   ! Calculate the dissipation
   !write(msg,"(A)") "Dissipation turned of for testing."
   !call elsi_say(eh%bh,msg)

   !alpha = 0.0_r8 ! turn of dissipation for testing
   call elsi_ccdm_dissipation_real(eh,k_max,alpha,c_k,X_array,diss_tmp)

   ! Calculate the second time derivative of X, The constant c is here set to 0.6 as a random stable guess, accelaration term
   call elsi_ccdm_scaled_delta_approx_real(eh,residual,0.6_r8,td_tmp)

   ! Sum all contributions to the integrated ccdm
   X_integrated = 2.0_r8 * X_array(current_index,:,:) - X_array(previous_index,:,:) + kappa * td_tmp + diss_tmp

   call elsi_get_time(t1) ! Timestamp

   write(msg,"(A)") "Finished integration of X."

   call elsi_say(eh%bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(eh%bh,msg)

   call elsi_deallocate(eh%bh,diss_tmp,"diss_tmp")
   call elsi_deallocate(eh%bh,td_tmp,"td_tmp")

end subroutine

!>
!! Routine to integrate the Contra-covariant DM X as specified in the 
!! XL-BOMD approach.
!!
!! X(t+delta t) = 
!!   2*X(t) - X(t-delta t) +  kappa * dX(t)/dt^2 +
!!   alpha*sum_k (c_k * X(t-k*delta t))
!!   = 2*X(t) - X(t-delta t) + kappa * integrated_X + dissipative_X
!!
!!
subroutine elsi_ccdm_integration_complex(eh,current_index,k_max,md_timestep,ccdm_current,X_array,X_integrated,residual)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: current_index !< Current index / MD step, needed to define the current and previous ccdm
   integer(kind=i4), intent(in) :: k_max !<  # of saved ccdm
   real(kind=r8), intent(in) :: md_timestep !< MD timestep
   complex(kind=r8), intent(in) :: ccdm_current(eh%bh%n_lrow,eh%bh%n_lcol) !< multiple instances of ccdm
   complex(kind=r8), intent(in) :: X_array(k_max,eh%bh%n_lrow,eh%bh%n_lcol) !< multiple instances of ccdm
   complex(kind=r8), intent(out) :: X_integrated(eh%bh%n_lrow,eh%bh%n_lcol) !< integrated ccdm
   complex(kind=r8), intent(out) :: residual(eh%bh%n_lrow,eh%bh%n_lcol) !< residual of current and previous step

   character(len=*), parameter :: caller = "elsi_ccdm_integration_real"

   real(kind=r8) :: t0, t1
   real(kind=r8) :: alpha, omega, kappa
   character(len=200) :: msg
   real(kind=r8) :: c_k(k_max)
   integer(kind=i4) :: previous_index

   complex(kind=r8), allocatable :: diss_tmp(:,:)
   complex(kind=r8), allocatable :: td_tmp(:,:)

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   call elsi_allocate(eh%bh,diss_tmp,eh%bh%n_lrow,eh%bh%n_lcol,"diss_tmp",caller)
   call elsi_allocate(eh%bh,td_tmp,eh%bh%n_lrow,eh%bh%n_lcol,"td_tmp",caller)

   call elsi_get_time(t0) ! Timestamp

   X_integrated = (0.0_r8,0.0_r8)

   ! Init the relevant constants.
   call elsi_init_ccdm_dissipation(eh,md_timestep,k_max,c_k,kappa,alpha,omega)

   ! Get previous array index based on current one
   if (current_index==1) then
       previous_index = k_max
   else
       previous_index = current_index - 1
   end if

   write(msg,"(A,I2,A)") "Current index :",current_index
   call elsi_say(eh%bh,msg)
   write(msg,"(A,I2,A)") "Previous index :",previous_index
   call elsi_say(eh%bh,msg)


   ! Get residual
   call elsi_residual_ccdm_complex(eh,ccdm_current,X_array(current_index,:,:),residual)

   ! Calculate the dissipation
   !write(msg,"(A)") "Dissipation turned of for testing."
   !call elsi_say(eh%bh,msg)

   !alpha = 0 ! turn of dissipation for testing
   call elsi_ccdm_dissipation_complex(eh,k_max,alpha,c_k,X_array,diss_tmp)

   ! Calculate the second time derivative of X, The constant c is here set to 0.6 as a random stable guess, accelaration term
   call elsi_ccdm_scaled_delta_approx_complex(eh,residual,0.6_r8,td_tmp)

   ! Sum all contributions to the integrated ccdm
   X_integrated = (2.0_r8,0.0_r8) * X_array(current_index,:,:) - X_array(previous_index,:,:) + kappa * td_tmp + diss_tmp

   call elsi_get_time(t1) ! Timestamp

   write(msg,"(A)") "Finished integration of X."

   call elsi_say(eh%bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(eh%bh,msg)

   call elsi_deallocate(eh%bh,diss_tmp,"diss_tmp")
   call elsi_deallocate(eh%bh,td_tmp,"td_tmp")

end subroutine


end module ELSI_GEO
