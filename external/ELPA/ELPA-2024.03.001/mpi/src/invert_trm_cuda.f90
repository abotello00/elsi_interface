













! This file is auto-generated. Do NOT edit

module invert_trm_cuda
  use, intrinsic :: iso_c_binding
  use precision

  implicit none

  public

  interface
    subroutine cuda_copy_double_a_tmat2_c(a_dev, tmat2_dev, nblk, matrixRows, l_cols, &
                                                             l_colx, l_row1, nb, my_stream)&
                                                     bind(C, name="cuda_copy_double_a_tmat2_FromC")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=C_intptr_T), value  :: a_dev, tmat2_dev
      integer(kind=c_int), intent(in)  :: nblk, matrixRows, l_cols, l_colx, l_row1, nb
      integer(kind=c_intptr_t), value  :: my_stream
    end subroutine
  end interface

  interface
    subroutine cuda_copy_double_tmp2_tmat2_c(tmp2_dev, tmat2_dev, nblk, l_col1, &
                                                                 nb, my_stream)&
                                                     bind(C, name="cuda_copy_double_tmp2_tmat2_FromC")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=C_intptr_T), value  :: tmp2_dev, tmat2_dev
      integer(kind=c_int), intent(in)  :: nblk, l_col1, nb
      integer(kind=c_intptr_t), value  :: my_stream
    end subroutine
  end interface

  interface
    subroutine cuda_copy_double_a_tmat1_c(a_dev, tmat1_dev, l_rows, matrixRows, nb, l_row1, &
                                                       l_col1, my_stream)&
                                                     bind(C, name="cuda_copy_double_a_tmat1_FromC")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=C_intptr_T), value  :: a_dev, tmat1_dev
      integer(kind=c_int), intent(in)  :: l_rows, matrixRows, nb, l_row1, l_col1
      integer(kind=c_intptr_t), value  :: my_stream
    end subroutine
  end interface

  interface
    subroutine cuda_copy_double_tmp1_tmp2_c(tmp1_dev, tmp2_dev, nblk, nb, &
                                                       my_stream)&
                                                     bind(C, name="cuda_copy_double_tmp1_tmp2_FromC")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=C_intptr_T), value  :: tmp1_dev, tmp2_dev
      integer(kind=c_int), intent(in)  :: nblk, nb
      integer(kind=c_intptr_t), value  :: my_stream
    end subroutine
  end interface

  interface
    subroutine cuda_copy_double_a_tmp1_c(a_dev, tmp1_dev, l_row1, l_col1, matrixRows, &
                                                       nb, my_stream)&
                                                     bind(C, name="cuda_copy_double_a_tmp1_FromC")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=C_intptr_T), value  :: a_dev, tmp1_dev
      integer(kind=c_int), intent(in)  :: l_row1, l_col1, matrixRows, nb
      integer(kind=c_intptr_t), value  :: my_stream
    end subroutine
  end interface

  interface
    subroutine cuda_copy_float_a_tmat2_c(a_dev, tmat2_dev, nblk, matrixRows, l_cols, &
                                                             l_colx, l_row1, nb, my_stream)&
                                                     bind(C, name="cuda_copy_float_a_tmat2_FromC")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=C_intptr_T), value  :: a_dev, tmat2_dev
      integer(kind=c_int), intent(in)  :: nblk, matrixRows, l_cols, l_colx, l_row1, nb
      integer(kind=c_intptr_t), value  :: my_stream
    end subroutine
  end interface

  interface
    subroutine cuda_copy_float_tmp2_tmat2_c(tmp2_dev, tmat2_dev, nblk, l_col1, &
                                                                 nb, my_stream)&
                                                     bind(C, name="cuda_copy_float_tmp2_tmat2_FromC")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=C_intptr_T), value  :: tmp2_dev, tmat2_dev
      integer(kind=c_int), intent(in)  :: nblk, l_col1, nb
      integer(kind=c_intptr_t), value  :: my_stream
    end subroutine
  end interface

  interface
    subroutine cuda_copy_float_a_tmat1_c(a_dev, tmat1_dev, l_rows, matrixRows, nb, l_row1, &
                                                       l_col1, my_stream)&
                                                     bind(C, name="cuda_copy_float_a_tmat1_FromC")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=C_intptr_T), value  :: a_dev, tmat1_dev
      integer(kind=c_int), intent(in)  :: l_rows, matrixRows, nb, l_row1, l_col1
      integer(kind=c_intptr_t), value  :: my_stream
    end subroutine
  end interface

  interface
    subroutine cuda_copy_float_tmp1_tmp2_c(tmp1_dev, tmp2_dev, nblk, nb, &
                                                       my_stream)&
                                                     bind(C, name="cuda_copy_float_tmp1_tmp2_FromC")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=C_intptr_T), value  :: tmp1_dev, tmp2_dev
      integer(kind=c_int), intent(in)  :: nblk, nb
      integer(kind=c_intptr_t), value  :: my_stream
    end subroutine
  end interface

  interface
    subroutine cuda_copy_float_a_tmp1_c(a_dev, tmp1_dev, l_row1, l_col1, matrixRows, &
                                                       nb, my_stream)&
                                                     bind(C, name="cuda_copy_float_a_tmp1_FromC")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=C_intptr_T), value  :: a_dev, tmp1_dev
      integer(kind=c_int), intent(in)  :: l_row1, l_col1, matrixRows, nb
      integer(kind=c_intptr_t), value  :: my_stream
    end subroutine
  end interface

  interface
    subroutine cuda_copy_double_complex_a_tmat2_c(a_dev, tmat2_dev, nblk, matrixRows, l_cols, &
                                                             l_colx, l_row1, nb, my_stream)&
                                                     bind(C, name="cuda_copy_double_complex_a_tmat2_FromC")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=C_intptr_T), value  :: a_dev, tmat2_dev
      integer(kind=c_int), intent(in)  :: nblk, matrixRows, l_cols, l_colx, l_row1, nb
      integer(kind=c_intptr_t), value  :: my_stream
    end subroutine
  end interface

  interface
    subroutine cuda_copy_double_complex_tmp2_tmat2_c(tmp2_dev, tmat2_dev, nblk, l_col1, &
                                                                 nb, my_stream)&
                                                     bind(C, name="cuda_copy_double_complex_tmp2_tmat2_FromC")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=C_intptr_T), value  :: tmp2_dev, tmat2_dev
      integer(kind=c_int), intent(in)  :: nblk, l_col1, nb
      integer(kind=c_intptr_t), value  :: my_stream
    end subroutine
  end interface

  interface
    subroutine cuda_copy_double_complex_a_tmat1_c(a_dev, tmat1_dev, l_rows, matrixRows, nb, l_row1, &
                                                       l_col1, my_stream)&
                                                     bind(C, name="cuda_copy_double_complex_a_tmat1_FromC")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=C_intptr_T), value  :: a_dev, tmat1_dev
      integer(kind=c_int), intent(in)  :: l_rows, matrixRows, nb, l_row1, l_col1
      integer(kind=c_intptr_t), value  :: my_stream
    end subroutine
  end interface

  interface
    subroutine cuda_copy_double_complex_tmp1_tmp2_c(tmp1_dev, tmp2_dev, nblk, nb, &
                                                       my_stream)&
                                                     bind(C, name="cuda_copy_double_complex_tmp1_tmp2_FromC")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=C_intptr_T), value  :: tmp1_dev, tmp2_dev
      integer(kind=c_int), intent(in)  :: nblk, nb
      integer(kind=c_intptr_t), value  :: my_stream
    end subroutine
  end interface

  interface
    subroutine cuda_copy_double_complex_a_tmp1_c(a_dev, tmp1_dev, l_row1, l_col1, matrixRows, &
                                                       nb, my_stream)&
                                                     bind(C, name="cuda_copy_double_complex_a_tmp1_FromC")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=C_intptr_T), value  :: a_dev, tmp1_dev
      integer(kind=c_int), intent(in)  :: l_row1, l_col1, matrixRows, nb
      integer(kind=c_intptr_t), value  :: my_stream
    end subroutine
  end interface

  interface
    subroutine cuda_copy_float_complex_a_tmat2_c(a_dev, tmat2_dev, nblk, matrixRows, l_cols, &
                                                             l_colx, l_row1, nb, my_stream)&
                                                     bind(C, name="cuda_copy_float_complex_a_tmat2_FromC")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=C_intptr_T), value  :: a_dev, tmat2_dev
      integer(kind=c_int), intent(in)  :: nblk, matrixRows, l_cols, l_colx, l_row1, nb
      integer(kind=c_intptr_t), value  :: my_stream
    end subroutine
  end interface

  interface
    subroutine cuda_copy_float_complex_tmp2_tmat2_c(tmp2_dev, tmat2_dev, nblk, l_col1, &
                                                                 nb, my_stream)&
                                                     bind(C, name="cuda_copy_float_complex_tmp2_tmat2_FromC")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=C_intptr_T), value  :: tmp2_dev, tmat2_dev
      integer(kind=c_int), intent(in)  :: nblk, l_col1, nb
      integer(kind=c_intptr_t), value  :: my_stream
    end subroutine
  end interface

  interface
    subroutine cuda_copy_float_complex_a_tmat1_c(a_dev, tmat1_dev, l_rows, matrixRows, nb, l_row1, &
                                                       l_col1, my_stream)&
                                                     bind(C, name="cuda_copy_float_complex_a_tmat1_FromC")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=C_intptr_T), value  :: a_dev, tmat1_dev
      integer(kind=c_int), intent(in)  :: l_rows, matrixRows, nb, l_row1, l_col1
      integer(kind=c_intptr_t), value  :: my_stream
    end subroutine
  end interface

  interface
    subroutine cuda_copy_float_complex_tmp1_tmp2_c(tmp1_dev, tmp2_dev, nblk, nb, &
                                                       my_stream)&
                                                     bind(C, name="cuda_copy_float_complex_tmp1_tmp2_FromC")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=C_intptr_T), value  :: tmp1_dev, tmp2_dev
      integer(kind=c_int), intent(in)  :: nblk, nb
      integer(kind=c_intptr_t), value  :: my_stream
    end subroutine
  end interface

  interface
    subroutine cuda_copy_float_complex_a_tmp1_c(a_dev, tmp1_dev, l_row1, l_col1, matrixRows, &
                                                       nb, my_stream)&
                                                     bind(C, name="cuda_copy_float_complex_a_tmp1_FromC")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=C_intptr_T), value  :: a_dev, tmp1_dev
      integer(kind=c_int), intent(in)  :: l_row1, l_col1, matrixRows, nb
      integer(kind=c_intptr_t), value  :: my_stream
    end subroutine
  end interface

  contains
    subroutine cuda_copy_double_a_tmat2(a_dev, tmat2_dev, nblk, matrixRows, l_cols, l_colx, &
                                                            l_row1, nb, my_stream)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: nblk, matrixRows, l_cols, l_colx, l_row1, nb
      integer(kind=C_intptr_T)        :: a_dev, tmat2_dev
      integer(kind=c_intptr_t)        :: my_stream

      call cuda_copy_double_a_tmat2_c(a_dev, tmat2_dev, nblk, matrixRows, l_cols, l_colx, &
                                                          l_row1, nb, my_stream)

    end subroutine

    subroutine cuda_copy_double_tmp2_tmat2(tmp2_dev, tmat2_dev, nblk, l_col1, &
                                                             nb, my_stream)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: nblk, l_col1, nb
      integer(kind=C_intptr_T)        :: tmp2_dev, tmat2_dev
      integer(kind=c_intptr_t)        :: my_stream

      call cuda_copy_double_tmp2_tmat2_c(tmp2_dev, tmat2_dev, nblk, l_col1, nb, &
                                                          my_stream)

    end subroutine

    subroutine cuda_copy_double_a_tmat1(a_dev, tmat1_dev, l_rows, matrixRows, nb, l_row1, l_col1,&
                                                            my_stream)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: l_rows, matrixRows, nb, l_row1, l_col1
      integer(kind=C_intptr_T)        :: a_dev, tmat1_dev
      integer(kind=c_intptr_t)        :: my_stream

      call cuda_copy_double_a_tmat1_c(a_dev, tmat1_dev, l_rows, matrixRows, nb, l_row1, l_col1, &
                                                          my_stream)

    end subroutine

    subroutine cuda_copy_double_tmp1_tmp2(tmp1_dev, tmp2_dev, nblk, nb,&
                                                            my_stream)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: nblk, nb
      integer(kind=C_intptr_T)        :: tmp1_dev, tmp2_dev
      integer(kind=c_intptr_t)        :: my_stream

      call cuda_copy_double_tmp1_tmp2_c(tmp1_dev, tmp2_dev, nblk, nb, &
                                                            my_stream)

    end subroutine

    subroutine cuda_copy_double_a_tmp1(a_dev, tmp1_dev, l_row1, l_col1, matrixRows, nb,&
                                                            my_stream)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: l_row1, l_col1, matrixRows, nb
      integer(kind=C_intptr_T)        :: a_dev, tmp1_dev
      integer(kind=c_intptr_t)        :: my_stream

      call cuda_copy_double_a_tmp1_c(a_dev, tmp1_dev, l_row1, l_col1, matrixRows, nb, &
                                                            my_stream)

    end subroutine

    subroutine cuda_copy_float_a_tmat2(a_dev, tmat2_dev, nblk, matrixRows, l_cols, l_colx, &
                                                            l_row1, nb, my_stream)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: nblk, matrixRows, l_cols, l_colx, l_row1, nb
      integer(kind=C_intptr_T)        :: a_dev, tmat2_dev
      integer(kind=c_intptr_t)        :: my_stream

      call cuda_copy_float_a_tmat2_c(a_dev, tmat2_dev, nblk, matrixRows, l_cols, l_colx, &
                                                          l_row1, nb, my_stream)

    end subroutine

    subroutine cuda_copy_float_tmp2_tmat2(tmp2_dev, tmat2_dev, nblk, l_col1, &
                                                             nb, my_stream)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: nblk, l_col1, nb
      integer(kind=C_intptr_T)        :: tmp2_dev, tmat2_dev
      integer(kind=c_intptr_t)        :: my_stream

      call cuda_copy_float_tmp2_tmat2_c(tmp2_dev, tmat2_dev, nblk, l_col1, nb, &
                                                          my_stream)

    end subroutine

    subroutine cuda_copy_float_a_tmat1(a_dev, tmat1_dev, l_rows, matrixRows, nb, l_row1, l_col1,&
                                                            my_stream)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: l_rows, matrixRows, nb, l_row1, l_col1
      integer(kind=C_intptr_T)        :: a_dev, tmat1_dev
      integer(kind=c_intptr_t)        :: my_stream

      call cuda_copy_float_a_tmat1_c(a_dev, tmat1_dev, l_rows, matrixRows, nb, l_row1, l_col1, &
                                                          my_stream)

    end subroutine

    subroutine cuda_copy_float_tmp1_tmp2(tmp1_dev, tmp2_dev, nblk, nb,&
                                                            my_stream)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: nblk, nb
      integer(kind=C_intptr_T)        :: tmp1_dev, tmp2_dev
      integer(kind=c_intptr_t)        :: my_stream

      call cuda_copy_float_tmp1_tmp2_c(tmp1_dev, tmp2_dev, nblk, nb, &
                                                            my_stream)

    end subroutine

    subroutine cuda_copy_float_a_tmp1(a_dev, tmp1_dev, l_row1, l_col1, matrixRows, nb,&
                                                            my_stream)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: l_row1, l_col1, matrixRows, nb
      integer(kind=C_intptr_T)        :: a_dev, tmp1_dev
      integer(kind=c_intptr_t)        :: my_stream

      call cuda_copy_float_a_tmp1_c(a_dev, tmp1_dev, l_row1, l_col1, matrixRows, nb, &
                                                            my_stream)

    end subroutine

    subroutine cuda_copy_double_complex_a_tmat2(a_dev, tmat2_dev, nblk, matrixRows, l_cols, l_colx, &
                                                            l_row1, nb, my_stream)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: nblk, matrixRows, l_cols, l_colx, l_row1, nb
      integer(kind=C_intptr_T)        :: a_dev, tmat2_dev
      integer(kind=c_intptr_t)        :: my_stream

      call cuda_copy_double_complex_a_tmat2_c(a_dev, tmat2_dev, nblk, matrixRows, l_cols, l_colx, &
                                                          l_row1, nb, my_stream)

    end subroutine

    subroutine cuda_copy_double_complex_tmp2_tmat2(tmp2_dev, tmat2_dev, nblk, l_col1, &
                                                             nb, my_stream)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: nblk, l_col1, nb
      integer(kind=C_intptr_T)        :: tmp2_dev, tmat2_dev
      integer(kind=c_intptr_t)        :: my_stream

      call cuda_copy_double_complex_tmp2_tmat2_c(tmp2_dev, tmat2_dev, nblk, l_col1, nb, &
                                                          my_stream)

    end subroutine

    subroutine cuda_copy_double_complex_a_tmat1(a_dev, tmat1_dev, l_rows, matrixRows, nb, l_row1, l_col1,&
                                                            my_stream)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: l_rows, matrixRows, nb, l_row1, l_col1
      integer(kind=C_intptr_T)        :: a_dev, tmat1_dev
      integer(kind=c_intptr_t)        :: my_stream

      call cuda_copy_double_complex_a_tmat1_c(a_dev, tmat1_dev, l_rows, matrixRows, nb, l_row1, l_col1, &
                                                          my_stream)

    end subroutine

    subroutine cuda_copy_double_complex_tmp1_tmp2(tmp1_dev, tmp2_dev, nblk, nb,&
                                                            my_stream)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: nblk, nb
      integer(kind=C_intptr_T)        :: tmp1_dev, tmp2_dev
      integer(kind=c_intptr_t)        :: my_stream

      call cuda_copy_double_complex_tmp1_tmp2_c(tmp1_dev, tmp2_dev, nblk, nb, &
                                                            my_stream)

    end subroutine

    subroutine cuda_copy_double_complex_a_tmp1(a_dev, tmp1_dev, l_row1, l_col1, matrixRows, nb,&
                                                            my_stream)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: l_row1, l_col1, matrixRows, nb
      integer(kind=C_intptr_T)        :: a_dev, tmp1_dev
      integer(kind=c_intptr_t)        :: my_stream

      call cuda_copy_double_complex_a_tmp1_c(a_dev, tmp1_dev, l_row1, l_col1, matrixRows, nb, &
                                                            my_stream)

    end subroutine

    subroutine cuda_copy_float_complex_a_tmat2(a_dev, tmat2_dev, nblk, matrixRows, l_cols, l_colx, &
                                                            l_row1, nb, my_stream)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: nblk, matrixRows, l_cols, l_colx, l_row1, nb
      integer(kind=C_intptr_T)        :: a_dev, tmat2_dev
      integer(kind=c_intptr_t)        :: my_stream

      call cuda_copy_float_complex_a_tmat2_c(a_dev, tmat2_dev, nblk, matrixRows, l_cols, l_colx, &
                                                          l_row1, nb, my_stream)

    end subroutine

    subroutine cuda_copy_float_complex_tmp2_tmat2(tmp2_dev, tmat2_dev, nblk, l_col1, &
                                                             nb, my_stream)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: nblk, l_col1, nb
      integer(kind=C_intptr_T)        :: tmp2_dev, tmat2_dev
      integer(kind=c_intptr_t)        :: my_stream

      call cuda_copy_float_complex_tmp2_tmat2_c(tmp2_dev, tmat2_dev, nblk, l_col1, nb, &
                                                          my_stream)

    end subroutine

    subroutine cuda_copy_float_complex_a_tmat1(a_dev, tmat1_dev, l_rows, matrixRows, nb, l_row1, l_col1,&
                                                            my_stream)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: l_rows, matrixRows, nb, l_row1, l_col1
      integer(kind=C_intptr_T)        :: a_dev, tmat1_dev
      integer(kind=c_intptr_t)        :: my_stream

      call cuda_copy_float_complex_a_tmat1_c(a_dev, tmat1_dev, l_rows, matrixRows, nb, l_row1, l_col1, &
                                                          my_stream)

    end subroutine

    subroutine cuda_copy_float_complex_tmp1_tmp2(tmp1_dev, tmp2_dev, nblk, nb,&
                                                            my_stream)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: nblk, nb
      integer(kind=C_intptr_T)        :: tmp1_dev, tmp2_dev
      integer(kind=c_intptr_t)        :: my_stream

      call cuda_copy_float_complex_tmp1_tmp2_c(tmp1_dev, tmp2_dev, nblk, nb, &
                                                            my_stream)

    end subroutine

    subroutine cuda_copy_float_complex_a_tmp1(a_dev, tmp1_dev, l_row1, l_col1, matrixRows, nb,&
                                                            my_stream)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: l_row1, l_col1, matrixRows, nb
      integer(kind=C_intptr_T)        :: a_dev, tmp1_dev
      integer(kind=c_intptr_t)        :: my_stream

      call cuda_copy_float_complex_a_tmp1_c(a_dev, tmp1_dev, l_row1, l_col1, matrixRows, nb, &
                                                            my_stream)

    end subroutine

end module
