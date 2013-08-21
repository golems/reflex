!! -*- mode: F90; -*-
!!
!! Copyright (c) 2013, Georgia Tech Research Corporation
!! All rights reserved.
!!
!! Author(s): Neil T. Dantam <ntd@gatech.edu>
!! Georgia Tech Humanoid Robotics Lab
!! Under Direction of Prof. Mike Stilman <mstilman@cc.gatech.edu>
!!
!!
!! This file is provided under the following "BSD-style" License:
!!
!!
!!   Redistribution and use in source and binary forms, with or
!!   without modification, are permitted provided that the following
!!   conditions are met:
!!
!!   * Redistributions of source code must retain the above copyright
!!     notice, this list of conditions and the following disclaimer.
!!
!!   * Redistributions in binary form must reproduce the above
!!     copyright notice, this list of conditions and the following
!!     disclaimer in the documentation and/or other materials provided
!!     with the distribution.
!!
!!   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
!!   CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
!!   INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
!!   MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!!   DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
!!   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
!!   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
!!   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
!!   USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
!!   AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
!!   LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
!!   ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!!   POSSIBILITY OF SUCH DAMAGE.


module reflex_kin

  use ISO_C_BINDING
  use amino_la
  use amino_tf
  implicit none

contains

  !!! Matrices !!!

  pure subroutine rfx_kin_jac_col_rev( T_abs, axis_rel, pe_abs, J ) &
       bind( C, name="rfx_kin_jac_col_rev" )
    real(C_DOUBLE), intent(in) :: T_abs(3,4), axis_rel(3), pe_abs(3)
    real(C_DOUBLE), intent(out) :: J(6)
    real(C_DOUBLE) :: tmp(3)
    ! rotation
    call aa_tf_9(T_abs(:,1:3), axis_rel, J(4:6))
    ! translation
    tmp = pe_abs - T_abs(:,4)
    call aa_tf_cross( J(3:6), tmp, J(1:3) )
  end subroutine rfx_kin_jac_col_rev

  pure subroutine rfx_kin_tf_jac_rev( n, TT_abs, axis, pe, J, ldJ ) &
       bind( C, name="rfx_kin_tf_jac_rev" )
    integer(C_SIZE_T), intent(in), value :: ldJ, n
    real(C_DOUBLE), intent(in) :: TT_abs(3,4,n), axis(3,n), pe(3)
    real(C_DOUBLE), intent(inout) :: J(ldJ,n)
    integer(C_SIZE_T) :: i
    do i = 1,n
       call rfx_kin_jac_col_rev( TT_abs(:,:,i), axis(:,i), pe, J(1:6,i) )
    end do
  end subroutine rfx_kin_tf_jac_rev

  pure subroutine rfx_kin_tf_chain( n, T0, TT_rel, TT_abs ) &
       bind( C, name="rfx_kin_tf_chain" )
    integer(C_SIZE_T), intent(in), value ::  n
    real(C_DOUBLE), intent(in) :: T0(3,4), TT_rel(3,4,n)
    real(C_DOUBLE), intent(out) :: TT_abs(3,4,n)
    integer(C_SIZE_T) :: i

    call aa_tf_12chain( T0, TT_rel(:,:,1), TT_abs(:,:,1) )

    do i=2,n
       call aa_tf_12chain( TT_abs(:,:,i-1), TT_rel(:,:,i), TT_abs(:,:,i) )
    end do

  end subroutine rfx_kin_tf_chain


  pure subroutine rfx_kin_revchain( n, T0, TT_rel, Te_rel, axis, TT_abs, J, ldJ ) &
       bind( C, name="rfx_kin_revchain" )
    integer(C_SIZE_T), intent(in), value ::  n, ldJ
    real(C_DOUBLE), intent(in) :: T0(3,4), TT_rel(3,4,n), Te_rel(3,4), axis(3,n)
    real(C_DOUBLE), intent(out) :: TT_abs(3,4,n)
    real(C_DOUBLE), intent(inout) :: J(ldJ,n)
    real(C_DOUBLE) :: Te_abs(3,4)
    call rfx_kin_tf_chain( n, T0, TT_rel, TT_abs )
    call aa_tf_12chain( TT_abs(:,:,n), Te_rel, Te_abs )
    call rfx_kin_tf_jac_rev( n, TT_abs, axis, Te_abs(:,4), J, ldJ )
  end subroutine rfx_kin_revchain

  !!! Dual Quaternions !!!

  subroutine rfx_kin_duqu_jac_col_rev( T_abs, axis_rel, pe_abs, J ) &
       bind( C, name="rfx_kin_duqu_jac_col_rev" )
    real(C_DOUBLE), intent(in) :: T_abs(8), axis_rel(3), pe_abs(3)
    real(C_DOUBLE), intent(out) :: J(6)
    real(C_DOUBLE) :: tmp(3)
    ! rotation
    call aa_tf_qrot(T_abs(1:4), axis_rel, J(4:6))
    ! translation
    call aa_tf_duqu_trans( T_abs, tmp )
    tmp = pe_abs - tmp
    call aa_tf_cross( J(3:6), tmp, J(1:3) )
  end subroutine rfx_kin_duqu_jac_col_rev

  subroutine rfx_kin_duqu_jac_rev( n, TT_abs, axis, pe, J, ldJ ) &
       bind( C, name="rfx_kin_duqu_jac_rev" )
    integer(C_SIZE_T), intent(in), value :: ldJ, n
    real(C_DOUBLE), intent(in) :: TT_abs(8,n), axis(3,n), pe(3)
    real(C_DOUBLE), intent(inout) :: J(ldJ,n)
    integer(C_SIZE_T) :: i
    do i = 1,n
       call rfx_kin_duqu_jac_col_rev( TT_abs(:,i), axis(:,i), pe, J(1:6,i) )
    end do
  end subroutine rfx_kin_duqu_jac_rev

  subroutine rfx_kin_duqu_chain( n, T0, TT_rel, TT_abs ) &
       bind( C, name="rfx_kin_duqu_chain" )
    integer(C_SIZE_T), intent(in), value ::  n
    real(C_DOUBLE), intent(in) :: T0(8), TT_rel(8,n)
    real(C_DOUBLE), intent(out) :: TT_abs(8,n)
    integer(C_SIZE_T) :: i

    call aa_tf_duqu_mul( T0, TT_rel(:,1), TT_abs(:,1) )

    do i=2,n
       call aa_tf_duqu_mul( TT_abs(:,i-1), TT_rel(:,i), TT_abs(:,i) )
    end do

  end subroutine rfx_kin_duqu_chain

  subroutine rfx_kin_duqu_revchain( n, T0, TT_rel, Te_rel, axis, TT_abs, J, ldJ ) &
       bind( C, name="rfx_kin_duqu_revchain" )
    integer(C_SIZE_T), intent(in), value ::  n, ldJ
    real(C_DOUBLE), intent(in) :: T0(8), TT_rel(8,n), Te_rel(8), axis(3,n)
    real(C_DOUBLE), intent(out) :: TT_abs(8,n)
    real(C_DOUBLE), intent(inout) :: J(ldJ,n)
    real(C_DOUBLE) :: Te_abs(8), pe_abs(8)
    call rfx_kin_duqu_chain( n, T0, TT_rel, TT_abs )
    call aa_tf_duqu_mul( TT_abs(:,n), Te_rel, Te_abs )
    call aa_tf_duqu_trans( Te_abs, pe_abs )
    call rfx_kin_duqu_jac_rev( n, TT_abs, axis, pe_abs, J, ldJ )
  end subroutine rfx_kin_duqu_revchain


end module reflex_kin
