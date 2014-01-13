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


module reflex
  use ISO_C_BINDING
  use amino_la
  use amino_tf
  implicit none


  Interface
     subroutine rfx_lqg_kbf_gain_work( nx, nz, A, C, V, W, P, K ) &
          bind( C, name="rfx_lqg_kbf_gain_work")
       use ISO_C_BINDING
       integer(C_SIZE_T), intent(in), value :: nx, nz
       real(C_DOUBLE), intent(in) :: A(nx,nx), C(nz,nx), V(nx,nx), W(nz,nz)
       real(C_DOUBLE), intent(out) ::  P(nx,nx), K(nx,nz)
     end subroutine rfx_lqg_kbf_gain_work
  End Interface


  Interface
     subroutine rfx_lqg_kf_predict_cov(n, A, V, P) &
          bind(C,name="rfx_lqg_kf_predict_cov")
       use ISO_C_BINDING
       integer(C_SIZE_T), intent(in), value :: n
       real(C_DOUBLE), intent(in),    dimension(n,n) :: A, V
       real(C_DOUBLE), intent(inout), dimension(n,n) :: P
     end subroutine rfx_lqg_kf_predict_cov
  End Interface

  Interface
     function rfx_lqg_kf_correct_gain(nx, nz, C, P, W, K) result(i) &
          bind(C,name="rfx_lqg_kf_correct_gain")
       use ISO_C_BINDING
       integer(C_SIZE_T), intent(in), value :: nx, nz
       real(C_DOUBLE), intent(in)  :: C(nz, nx), P(nx,nx), W(nz,nz)
       real(C_DOUBLE), intent(out)  :: K(nx,nz)
       integer(C_INT) :: i
     end function rfx_lqg_kf_correct_gain
  End Interface

  Interface
     subroutine rfx_lqg_kf_correct_cov(nx, nz, C, P, K) &
          bind(C,name="rfx_lqg_kf_correct_cov")
       use ISO_C_BINDING
       integer(C_SIZE_T), intent(in), value :: nx, nz
       real(C_DOUBLE), intent(in)  :: C(nz, nx), K(nx,nz)
       real(C_DOUBLE), intent(inout)  :: P(nx,nx)
       integer(C_INT) :: i
     end subroutine rfx_lqg_kf_correct_cov
  End Interface
 contains

#include "kinematics.f90"
#include "lqg_tf_f.f90"

end module reflex
