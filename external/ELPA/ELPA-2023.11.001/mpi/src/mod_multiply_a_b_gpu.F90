












module multiply_a_b_gpu
  use, intrinsic :: iso_c_binding
  use precision


  implicit none

  public
  !This file is auto-generated do NOT edit!

  interface gpu_copy_double_a_aux_bc
    module procedure gpu_copy_double_a_aux_bc_intptr
    module procedure gpu_copy_double_a_aux_bc_cptr
  end interface

  interface gpu_copy_double_tmp2_c
    module procedure gpu_copy_double_tmp2_c_intptr
    module procedure gpu_copy_double_tmp2_c_cptr
  end interface
  
  interface gpu_copy_float_a_aux_bc
    module procedure gpu_copy_float_a_aux_bc_intptr
    module procedure gpu_copy_float_a_aux_bc_cptr
  end interface

  interface gpu_copy_float_tmp2_c
    module procedure gpu_copy_float_tmp2_c_intptr
    module procedure gpu_copy_float_tmp2_c_cptr
  end interface
  
  interface gpu_copy_double_complex_a_aux_bc
    module procedure gpu_copy_double_complex_a_aux_bc_intptr
    module procedure gpu_copy_double_complex_a_aux_bc_cptr
  end interface

  interface gpu_copy_double_complex_tmp2_c
    module procedure gpu_copy_double_complex_tmp2_c_intptr
    module procedure gpu_copy_double_complex_tmp2_c_cptr
  end interface
  
  interface gpu_copy_float_complex_a_aux_bc
    module procedure gpu_copy_float_complex_a_aux_bc_intptr
    module procedure gpu_copy_float_complex_a_aux_bc_cptr
  end interface

  interface gpu_copy_float_complex_tmp2_c
    module procedure gpu_copy_float_complex_tmp2_c_intptr
    module procedure gpu_copy_float_complex_tmp2_c_cptr
  end interface
  

  contains

    subroutine gpu_copy_double_tmp2_c_intptr(tmp2_dev, c_dev, nr_done, nstor, lcs, lce, ldc, ldcCols, my_stream)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: nr_done, nstor, lcs, lce, ldc, ldcCols
      integer(kind=C_intptr_T)        :: tmp2_dev, c_dev
      integer(kind=C_intptr_T)        :: my_stream





    end subroutine

    subroutine gpu_copy_double_tmp2_c_cptr(tmp2_dev, c_dev, nr_done, nstor, lcs, lce, ldc, ldcCols, my_stream)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: nr_done, nstor, lcs, lce, ldc, ldcCols
      integer(kind=C_intptr_T)        :: tmp2_dev
      type(c_ptr)                     :: c_dev
      integer(kind=C_intptr_T)        :: my_stream





    end subroutine

    subroutine gpu_copy_double_a_aux_bc_intptr(a_dev, aux_bc_dev, n_aux_bc, nvals, lrs, lre, noff, nblk, n, l_rows, &
                                        lda, ldaCols, my_stream)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: n_aux_bc, nvals, lrs, lre, noff, nblk, n, l_rows, lda, ldaCols
      integer(kind=C_intptr_T)        :: a_dev, aux_bc_dev
      integer(kind=C_intptr_T)        :: my_stream





    end subroutine

    subroutine gpu_copy_double_a_aux_bc_cptr(a_dev, aux_bc_dev, n_aux_bc, nvals, lrs, lre, noff, nblk, n, l_rows, &
                                        lda, ldaCols, my_stream)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: n_aux_bc, nvals, lrs, lre, noff, nblk, n, l_rows, lda, ldaCols
      type(c_ptr)                     :: a_dev
      integer(kind=C_intptr_T)        :: aux_bc_dev
      integer(kind=C_intptr_T)        :: my_stream





    end subroutine

    subroutine gpu_copy_double_aux_bc_aux_mat(aux_bc_dev, aux_mat_dev, lrs, lre, nstor, n_aux_bc, nvals, &
                                        l_rows, nblk, nblk_mult, my_stream)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: lrs, lre, nstor, n_aux_bc, nvals, l_rows, nblk, nblk_mult
      integer(kind=C_intptr_T)        :: aux_mat_dev, aux_bc_dev
      integer(kind=C_intptr_T)        :: my_stream





    end subroutine

    subroutine gpu_copy_float_tmp2_c_intptr(tmp2_dev, c_dev, nr_done, nstor, lcs, lce, ldc, ldcCols, my_stream)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: nr_done, nstor, lcs, lce, ldc, ldcCols
      integer(kind=C_intptr_T)        :: tmp2_dev, c_dev
      integer(kind=C_intptr_T)        :: my_stream





    end subroutine

    subroutine gpu_copy_float_tmp2_c_cptr(tmp2_dev, c_dev, nr_done, nstor, lcs, lce, ldc, ldcCols, my_stream)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: nr_done, nstor, lcs, lce, ldc, ldcCols
      integer(kind=C_intptr_T)        :: tmp2_dev
      type(c_ptr)                     :: c_dev
      integer(kind=C_intptr_T)        :: my_stream





    end subroutine

    subroutine gpu_copy_float_a_aux_bc_intptr(a_dev, aux_bc_dev, n_aux_bc, nvals, lrs, lre, noff, nblk, n, l_rows, &
                                        lda, ldaCols, my_stream)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: n_aux_bc, nvals, lrs, lre, noff, nblk, n, l_rows, lda, ldaCols
      integer(kind=C_intptr_T)        :: a_dev, aux_bc_dev
      integer(kind=C_intptr_T)        :: my_stream





    end subroutine

    subroutine gpu_copy_float_a_aux_bc_cptr(a_dev, aux_bc_dev, n_aux_bc, nvals, lrs, lre, noff, nblk, n, l_rows, &
                                        lda, ldaCols, my_stream)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: n_aux_bc, nvals, lrs, lre, noff, nblk, n, l_rows, lda, ldaCols
      type(c_ptr)                     :: a_dev
      integer(kind=C_intptr_T)        :: aux_bc_dev
      integer(kind=C_intptr_T)        :: my_stream





    end subroutine

    subroutine gpu_copy_float_aux_bc_aux_mat(aux_bc_dev, aux_mat_dev, lrs, lre, nstor, n_aux_bc, nvals, &
                                        l_rows, nblk, nblk_mult, my_stream)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: lrs, lre, nstor, n_aux_bc, nvals, l_rows, nblk, nblk_mult
      integer(kind=C_intptr_T)        :: aux_mat_dev, aux_bc_dev
      integer(kind=C_intptr_T)        :: my_stream





    end subroutine

    subroutine gpu_copy_double_complex_tmp2_c_intptr(tmp2_dev, c_dev, nr_done, nstor, lcs, lce, ldc, ldcCols, my_stream)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: nr_done, nstor, lcs, lce, ldc, ldcCols
      integer(kind=C_intptr_T)        :: tmp2_dev, c_dev
      integer(kind=C_intptr_T)        :: my_stream





    end subroutine

    subroutine gpu_copy_double_complex_tmp2_c_cptr(tmp2_dev, c_dev, nr_done, nstor, lcs, lce, ldc, ldcCols, my_stream)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: nr_done, nstor, lcs, lce, ldc, ldcCols
      integer(kind=C_intptr_T)        :: tmp2_dev
      type(c_ptr)                     :: c_dev
      integer(kind=C_intptr_T)        :: my_stream





    end subroutine

    subroutine gpu_copy_double_complex_a_aux_bc_intptr(a_dev, aux_bc_dev, n_aux_bc, nvals, lrs, lre, noff, nblk, n, l_rows, &
                                        lda, ldaCols, my_stream)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: n_aux_bc, nvals, lrs, lre, noff, nblk, n, l_rows, lda, ldaCols
      integer(kind=C_intptr_T)        :: a_dev, aux_bc_dev
      integer(kind=C_intptr_T)        :: my_stream





    end subroutine

    subroutine gpu_copy_double_complex_a_aux_bc_cptr(a_dev, aux_bc_dev, n_aux_bc, nvals, lrs, lre, noff, nblk, n, l_rows, &
                                        lda, ldaCols, my_stream)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: n_aux_bc, nvals, lrs, lre, noff, nblk, n, l_rows, lda, ldaCols
      type(c_ptr)                     :: a_dev
      integer(kind=C_intptr_T)        :: aux_bc_dev
      integer(kind=C_intptr_T)        :: my_stream





    end subroutine

    subroutine gpu_copy_double_complex_aux_bc_aux_mat(aux_bc_dev, aux_mat_dev, lrs, lre, nstor, n_aux_bc, nvals, &
                                        l_rows, nblk, nblk_mult, my_stream)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: lrs, lre, nstor, n_aux_bc, nvals, l_rows, nblk, nblk_mult
      integer(kind=C_intptr_T)        :: aux_mat_dev, aux_bc_dev
      integer(kind=C_intptr_T)        :: my_stream





    end subroutine

    subroutine gpu_copy_float_complex_tmp2_c_intptr(tmp2_dev, c_dev, nr_done, nstor, lcs, lce, ldc, ldcCols, my_stream)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: nr_done, nstor, lcs, lce, ldc, ldcCols
      integer(kind=C_intptr_T)        :: tmp2_dev, c_dev
      integer(kind=C_intptr_T)        :: my_stream





    end subroutine

    subroutine gpu_copy_float_complex_tmp2_c_cptr(tmp2_dev, c_dev, nr_done, nstor, lcs, lce, ldc, ldcCols, my_stream)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: nr_done, nstor, lcs, lce, ldc, ldcCols
      integer(kind=C_intptr_T)        :: tmp2_dev
      type(c_ptr)                     :: c_dev
      integer(kind=C_intptr_T)        :: my_stream





    end subroutine

    subroutine gpu_copy_float_complex_a_aux_bc_intptr(a_dev, aux_bc_dev, n_aux_bc, nvals, lrs, lre, noff, nblk, n, l_rows, &
                                        lda, ldaCols, my_stream)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: n_aux_bc, nvals, lrs, lre, noff, nblk, n, l_rows, lda, ldaCols
      integer(kind=C_intptr_T)        :: a_dev, aux_bc_dev
      integer(kind=C_intptr_T)        :: my_stream





    end subroutine

    subroutine gpu_copy_float_complex_a_aux_bc_cptr(a_dev, aux_bc_dev, n_aux_bc, nvals, lrs, lre, noff, nblk, n, l_rows, &
                                        lda, ldaCols, my_stream)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: n_aux_bc, nvals, lrs, lre, noff, nblk, n, l_rows, lda, ldaCols
      type(c_ptr)                     :: a_dev
      integer(kind=C_intptr_T)        :: aux_bc_dev
      integer(kind=C_intptr_T)        :: my_stream





    end subroutine

    subroutine gpu_copy_float_complex_aux_bc_aux_mat(aux_bc_dev, aux_mat_dev, lrs, lre, nstor, n_aux_bc, nvals, &
                                        l_rows, nblk, nblk_mult, my_stream)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: lrs, lre, nstor, n_aux_bc, nvals, l_rows, nblk, nblk_mult
      integer(kind=C_intptr_T)        :: aux_mat_dev, aux_bc_dev
      integer(kind=C_intptr_T)        :: my_stream





    end subroutine

end module
