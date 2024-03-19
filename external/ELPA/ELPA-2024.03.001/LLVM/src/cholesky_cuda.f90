













! This file is auto-generated. Do NOT edit

module cholesky_cuda
  use, intrinsic :: iso_c_binding
  use precision

  implicit none

  public

  interface
    subroutine cuda_check_device_info_c(info_dev,  my_stream)&
                 bind(C, name="cuda_check_device_info_FromC")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_intptr_t), value  :: info_dev
      integer(kind=c_intptr_t), value  :: my_stream
    end subroutine
  end interface

  interface
    subroutine cuda_accumulate_device_info_c(info_abs_dev, info_new_dev, my_stream)&
                 bind(C, name="cuda_accumulate_device_info_FromC")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_intptr_t), value  :: info_abs_dev, info_new_dev
      integer(kind=c_intptr_t), value  :: my_stream
    end subroutine
  end interface

  interface
    subroutine cuda_copy_double_a_tmatc_c(a_dev, tmatc_dev, nblk, matrixRows, l_cols, &
                                                             l_colx, l_row1, my_stream)&
                                                     bind(C, name="cuda_copy_double_a_tmatc_FromC")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=C_intptr_T), value  :: a_dev, tmatc_dev
      integer(kind=c_int), intent(in)  :: nblk, matrixRows, l_cols, l_colx, l_row1
      integer(kind=c_intptr_t), value  :: my_stream
    end subroutine
  end interface

  interface
    subroutine cuda_copy_float_a_tmatc_c(a_dev, tmatc_dev, nblk, matrixRows, l_cols, &
                                                             l_colx, l_row1, my_stream)&
                                                     bind(C, name="cuda_copy_float_a_tmatc_FromC")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=C_intptr_T), value  :: a_dev, tmatc_dev
      integer(kind=c_int), intent(in)  :: nblk, matrixRows, l_cols, l_colx, l_row1
      integer(kind=c_intptr_t), value  :: my_stream
    end subroutine
  end interface

  interface
    subroutine cuda_copy_double_complex_a_tmatc_c(a_dev, tmatc_dev, nblk, matrixRows, l_cols, &
                                                             l_colx, l_row1, my_stream)&
                                                     bind(C, name="cuda_copy_double_complex_a_tmatc_FromC")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=C_intptr_T), value  :: a_dev, tmatc_dev
      integer(kind=c_int), intent(in)  :: nblk, matrixRows, l_cols, l_colx, l_row1
      integer(kind=c_intptr_t), value  :: my_stream
    end subroutine
  end interface

  interface
    subroutine cuda_copy_float_complex_a_tmatc_c(a_dev, tmatc_dev, nblk, matrixRows, l_cols, &
                                                             l_colx, l_row1, my_stream)&
                                                     bind(C, name="cuda_copy_float_complex_a_tmatc_FromC")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=C_intptr_T), value  :: a_dev, tmatc_dev
      integer(kind=c_int), intent(in)  :: nblk, matrixRows, l_cols, l_colx, l_row1
      integer(kind=c_intptr_t), value  :: my_stream
    end subroutine
  end interface


  contains

    subroutine cuda_check_device_info(info_dev, my_stream)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=c_intptr_t)        :: info_dev
      integer(kind=c_intptr_t)        :: my_stream

      call cuda_check_device_info_c(info_dev, my_stream)
    end subroutine

    subroutine cuda_accumulate_device_info(info_abs_dev, info_new_dev, my_stream)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=c_intptr_t)        :: info_abs_dev, info_new_dev
      integer(kind=c_intptr_t)        :: my_stream

      call cuda_accumulate_device_info_c(info_abs_dev, info_new_dev, my_stream)
    end subroutine

    subroutine cuda_copy_double_a_tmatc(a_dev, tmatc_dev, nblk, matrixRows, l_cols, l_colx, &
                                                            l_row1, my_stream)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: nblk, matrixRows, l_cols, l_colx, l_row1
      integer(kind=C_intptr_T)        :: a_dev, tmatc_dev
      integer(kind=c_intptr_t)        :: my_stream

      call cuda_copy_double_a_tmatc_c(a_dev, tmatc_dev, nblk, matrixRows, l_cols, l_colx, &
                                                          l_row1, my_stream)

    end subroutine

    subroutine cuda_copy_float_a_tmatc(a_dev, tmatc_dev, nblk, matrixRows, l_cols, l_colx, &
                                                            l_row1, my_stream)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: nblk, matrixRows, l_cols, l_colx, l_row1
      integer(kind=C_intptr_T)        :: a_dev, tmatc_dev
      integer(kind=c_intptr_t)        :: my_stream

      call cuda_copy_float_a_tmatc_c(a_dev, tmatc_dev, nblk, matrixRows, l_cols, l_colx, &
                                                          l_row1, my_stream)

    end subroutine

    subroutine cuda_copy_double_complex_a_tmatc(a_dev, tmatc_dev, nblk, matrixRows, l_cols, l_colx, &
                                                            l_row1, my_stream)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: nblk, matrixRows, l_cols, l_colx, l_row1
      integer(kind=C_intptr_T)        :: a_dev, tmatc_dev
      integer(kind=c_intptr_t)        :: my_stream

      call cuda_copy_double_complex_a_tmatc_c(a_dev, tmatc_dev, nblk, matrixRows, l_cols, l_colx, &
                                                          l_row1, my_stream)

    end subroutine

    subroutine cuda_copy_float_complex_a_tmatc(a_dev, tmatc_dev, nblk, matrixRows, l_cols, l_colx, &
                                                            l_row1, my_stream)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: nblk, matrixRows, l_cols, l_colx, l_row1
      integer(kind=C_intptr_T)        :: a_dev, tmatc_dev
      integer(kind=c_intptr_t)        :: my_stream

      call cuda_copy_float_complex_a_tmatc_c(a_dev, tmatc_dev, nblk, matrixRows, l_cols, l_colx, &
                                                          l_row1, my_stream)

    end subroutine

end module
