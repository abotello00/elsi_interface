module ELSI_CHASE
   use chase_diag, only: dchase, zchase
   use ELSI_CONSTANT, only: FC_BASIC,FC_PLUS_V
   use ELSI_DATATYPE, only: elsi_param_t,elsi_basic_t
   use ELSI_ELPA, only: elsi_elpa_tridiag
   use ELSI_MALLOC, only: elsi_allocate,elsi_deallocate
   use ELSI_OUTPUT, only: elsi_say,elsi_get_time
   use ELSI_PRECISION, only: r8,i4,i8
   use ELSI_SORT, only: elsi_heapsort,elsi_permute
  
   use ELSI_LAPACK, only: elsi_factor_ovlp_sp,elsi_reduce_evp_sp,elsi_back_ev_sp,elsi_check_ovlp_sp

   implicit none
   private

   public :: elsi_solve_chase_sp

   interface elsi_solve_chase_sp
       module procedure elsi_solve_chase_real_sp
       module procedure elsi_solve_chase_cmplx_sp
   end interface

contains

subroutine elsi_solve_chase_real_sp(ph,bh,ham,ovlp,eval,evec)
   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   real(kind=r8), intent(inout) :: ham(ph%n_basis,ph%n_basis)
   real(kind=r8), intent(inout) :: ovlp(ph%n_basis,ph%n_basis)
   real(kind=r8), intent(out) :: eval(ph%n_states)
   real(kind=r8), intent(out) :: evec(ph%n_basis,ph%n_states)

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   integer(kind=i4) :: ierr
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_solve_chase_real_sp"
   
   !chase
   integer(kind=i4) :: nev, nex
   integer(kind=i4) :: i, j

   if(.not. ph%unit_ovlp .and. ph%ill_check) then
      call elsi_check_ovlp_sp(ph,bh,ovlp,eval,evec)
   end if
   
   nev = min(ph%n_states, ph%n_good)
   nex = min(max(int(nev * ph%chase_extra_space), ph%chase_min_extra_space), ph%n_good - nev)

   if(.not. ph%chase_started) then
      call elsi_allocate(bh, ph%pre_evec_real, ph%n_good, nev+nex,&
             "pre_evec_real",caller)
   end if
   
   if(.not. ph%chase_started) then
      call elsi_allocate(bh,ph%pre_eval, nev+nex,"pre_eval",caller)
   end if

   ! Transform to standard form
   if(.not. ph%unit_ovlp) then
      if(ph%n_good == ph%n_basis) then
         ! Do Cholesky if not singular
         call elsi_factor_ovlp_sp(ph,bh,ovlp)
      end if

      call elsi_reduce_evp_sp(ph,bh,ham,ovlp,evec)
   end if
   call elsi_get_time(t0)
   ! Explicitly ensure the symmetricity of ham
   ! Required by ChASE
   do j = 1, ph%n_basis
      do i = 1,j
           ham(j, i) = ham(i, j)
      end do
   end do
   
   ! solve
   if(ph%chase_started) then
      call dchase(ham, ph%n_good,ph%pre_evec_real, ph%pre_eval,nev, nex, ph%chase_filter_deg, ph%chase_tol, 'A', 'S')
   else
      call dchase(ham, ph%n_good,ph%pre_evec_real, ph%pre_eval,nev, nex, ph%chase_filter_deg, ph%chase_tol, 'R', 'S')   
   end if

   evec(1:ph%n_good,1:nev+nex) = ph%pre_evec_real(1:ph%n_good,1:nev+nex)
   eval(1:nev+nex) = ph%pre_eval(1:nev+nex)
   
   call elsi_get_time(t1)
   write(msg,"(A)") "Finished solving standard eigenproblem"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

   ! Back-transform eigenvectors
   if(.not. ph%unit_ovlp) then
      call elsi_back_ev_sp(ph,bh,ovlp,evec)
   end if

   ph%chase_started = .true.
   
   ! Switch back to full dimension in case of frozen core
   if(ph%n_basis_c > 0) then
      ph%n_basis = ph%n_basis_v+ph%n_basis_c
      ph%n_states = ph%n_states+ph%n_basis_c
      ph%n_good = ph%n_good+ph%n_basis_c
      ph%n_states_solve = ph%n_states_solve+ph%n_basis_c
   end if
   
end subroutine         

subroutine elsi_solve_chase_cmplx_sp(ph,bh,ham,ovlp,eval,evec)
   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   complex(kind=r8), intent(inout) :: ham(ph%n_basis,ph%n_basis)
   complex(kind=r8), intent(inout) :: ovlp(ph%n_basis,ph%n_basis)
   real(kind=r8), intent(out) :: eval(ph%n_basis)
   complex(kind=r8), intent(out) :: evec(ph%n_basis,ph%n_basis)

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   integer(kind=i4) :: ierr
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_solve_chase_cmplx_sp"

   !chase
   integer(kind=i4) :: nev, nex
   integer(kind=i4) :: i, j

   ! Ill-conditioning check
   if(.not. ph%unit_ovlp .and. ph%ill_check) then
      call elsi_check_ovlp_sp(ph,bh,ovlp,eval,evec)
   end if

   nev = min(ph%n_states, ph%n_good)
   nex = min(max(int(nev * ph%chase_extra_space), ph%chase_min_extra_space), ph%n_good - nev)
  
   if(.not. ph%chase_started) then
      call elsi_allocate(bh, ph%pre_evec_cmplx, ph%n_good, nev+nex,&
             "pre_evec_cmplx",caller)
   end if

   if(.not. ph%chase_started) then
      call elsi_allocate(bh,ph%pre_eval, nev+nex,"pre_eval",caller)
   end if

   ! Transform to standard form
   if(.not. ph%unit_ovlp) then
      if(ph%n_good == ph%n_basis) then
         ! Do Cholesky if not singular
         call elsi_factor_ovlp_sp(ph,bh,ovlp)
      end if

      call elsi_reduce_evp_sp(ph,bh,ham,ovlp,evec)
   end if

   call elsi_get_time(t0)
   ! Explicitly ensure the hermeticity of ham
   ! Required by ChASE
   DO j = 1, ph%n_basis
       DO i = 1,j
           ham(j, i) = ham(i, j)
       END DO
   END DO

   ! solve
   if(ph%chase_started) then
      call zchase(ham, ph%n_good,ph%pre_evec_cmplx, ph%pre_eval,nev, nex, ph%chase_filter_deg, ph%chase_tol, 'A', 'S')
   else
      call zchase(ham, ph%n_good,ph%pre_evec_cmplx, ph%pre_eval,nev, nex, ph%chase_filter_deg, ph%chase_tol, 'R', 'S')
   end if

   evec(1:ph%n_good,1:nev+nex) = ph%pre_evec_cmplx(1:ph%n_good,1:nev+nex)
   eval(1:nev+nex) = ph%pre_eval(1:nev+nex)

   call elsi_get_time(t1)
   write(msg,"(A)") "Finished solving standard eigenproblem"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)
   
   ! Back-transform eigenvectors
   if(.not. ph%unit_ovlp) then
      call elsi_back_ev_sp(ph,bh,ovlp,evec)
   end if

   ph%chase_started = .true.

   ! Switch back to full dimension in case of frozen core
   if(ph%n_basis_c > 0) then
      ph%n_basis = ph%n_basis_v+ph%n_basis_c
      ph%n_states = ph%n_states+ph%n_basis_c
      ph%n_good = ph%n_good+ph%n_basis_c
      ph%n_states_solve = ph%n_states_solve+ph%n_basis_c
   end if
end subroutine

end module ELSI_CHASE
