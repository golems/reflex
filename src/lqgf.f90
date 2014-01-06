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


!!! Kalman Filters !!!

subroutine rfx_lqg_kf_predict_work( n_x, n_u, A, x, B, u, P, V, work ) &
     bind( C, name="rfx_lqg_kf_predict_work" )
  integer(C_SIZE_T), intent(in), value :: n_x, n_u
  real(C_DOUBLE), intent(in) :: A(n_x,n_x), B(n_x, n_u), u(n_u), V(n_x,n_x)
  real(C_DOUBLE), intent(inout) :: x(n_x), P(n_x,n_x)
  real(C_DOUBLE), intent(out) :: work(n_x,n_x)

  call state(work(:,1))
  call covariance(work)

contains

  subroutine state( tmpx )
    real(C_DOUBLE), intent(out) :: tmpx(:)
    ! x = A*x + B*u
    tmpx = matmul(A,x)
    x = tmpx
    tmpx = matmul(B,u)
    x = x + tmpx
  end subroutine state

  subroutine covariance( tmpp )
    real(C_DOUBLE), intent(out) :: tmpp(:,:)
    ! P = A*P*A^T + V
    tmpp = matmul(A,P)
    P = matmul(tmpp,transpose(A))
    P = P + V
  end subroutine covariance

end subroutine rfx_lqg_kf_predict_work

subroutine rfx_lqg_kf_gain_work( n_x, n_z, C, P, W, K, work_xz, work_zz ) &
     bind( C, name="rfx_lqg_kf_gain_work" )
  integer(C_SIZE_T), intent(in), value :: n_x, n_z
  real(C_DOUBLE), intent(in) :: C(n_z,n_x), P(n_x, n_x), W(n_z,n_z)
  real(C_DOUBLE), intent(out) :: work_xz(n_x,n_z), work_zz(n_z,n_z), K(n_x, n_z)
  integer :: i
  !! K = P * C**T * (C * P * C**T + W)**-1
  work_xz = matmul(P,transpose(C))
  work_zz = matmul(C,work_xz)
  work_zz = work_zz + W
  i = aa_la_inv( n_z, work_zz )
  K = matmul(work_xz, work_zz)
end subroutine rfx_lqg_kf_gain_work

subroutine rfx_lqg_kf_innovate_work( n_x, n_z, x, K, z, C, workx, workz ) &
     bind( C, name="rfx_lqg_kf_innovate_work" )
  integer(C_SIZE_T) :: n_x, n_z
  real(C_DOUBLE), intent(inout) :: x(n_x)
  real(C_DOUBLE), intent(in) :: K(n_x, n_z), z(n_z), C(n_x,n_z)
  real(C_DOUBLE), intent(out) :: workz(n_z), workx(n_x)
  !! xh = x + K * (z - C*x)
  workz =  matmul(C,x)
  workz = z - workz
  workx = matmul(K,workz)
  x = x + workx
end subroutine rfx_lqg_kf_innovate_work

subroutine rfx_lqg_kf_covariance_work( n_x, n_z, C, P, K, work_xx2 ) &
     bind( C, name="rfx_lqg_kf_covariance_work" )
  integer(C_SIZE_T), intent(in), value :: n_x, n_z
  real(C_DOUBLE), intent(inout) :: P(n_x,n_x)
  real(C_DOUBLE), intent(in) :: C(n_z, n_x), K(n_x, n_z)
  real(C_DOUBLE), intent(out) :: work_xx2(n_x,n_x,2)
  integer(C_SIZE_T) :: i
  !! P = (I - K*C) * P
  work_xx2(:,:,1) = matmul(K, C)
  work_xx2(:,:,1) = -work_xx2(:,:,1)
  forall (i = 1:n_x)
     work_xx2(i,i,1) = 1 + work_xx2(i,i,1)
  end forall
  work_xx2(:,:,2) = P
  P = matmul( work_xx2(:,:,1), work_xx2(:,:,1) )
end subroutine rfx_lqg_kf_covariance_work


elemental subroutine rfx_lqg_kf_simple(  x, z, p, v, w, k )
  real(C_DOUBLE), intent(in) :: z,v,w
  real(C_DOUBLE), intent(inout) :: x,p
  real(C_DOUBLE), intent(out) :: k

  ! A = C = I

  ! predict
  p = p + v

  ! correct
  k = p / (p+w)
  x = x + k * (z - x)
  p = (1-k) * p
end subroutine rfx_lqg_kf_simple

subroutine rfx_lqg_kf_simple_work( n_x, x, z, P, V, W, K ) &
     bind( C, name="rfx_lqg_kf_simple_work" )
  integer(C_SIZE_T), intent(in), value :: n_x
  real(C_DOUBLE), intent(in) :: z(n_x), W(n_x), V(n_x)
  real(C_DOUBLE), intent(inout) :: x(n_x), P(n_x)
  real(C_DOUBLE), intent(out) :: K(n_x)
  call rfx_lqg_kf_simple( x, z, P, V, W, K )
end subroutine rfx_lqg_kf_simple_work

subroutine rfx_lqg_kf_predict_diag_work( n_x, A, x, P, V ) &
     bind( C, name="rfx_lqg_kf_predict_diag_work" )
  integer(C_SIZE_T), intent(in), value :: n_x
  real(C_DOUBLE), intent(in) :: A(n_x), V(n_x)
  real(C_DOUBLE), intent(inout) :: x(n_x), P(n_x)
  ! state
  x = A*x
  ! covariance
  P = A*P*A + V
end subroutine rfx_lqg_kf_predict_diag_work

subroutine rfx_lqg_kf_correct_diag_work( n_x, z, C, x, P, W, K ) &
     bind( C, name="rfx_lqg_kf_correct_diag_work" )
  integer(C_SIZE_T), intent(in), value :: n_x
  real(C_DOUBLE), intent(in) :: z(n_x), C(n_x), W(n_x)
  real(C_DOUBLE), intent(inout) :: x(n_x), P(n_x)
  real(C_DOUBLE), intent(out) :: K(n_x)

  ! kalman gain
  K = (P * C ) / (C*P*C + W)

  ! state
  x = x + K * (z - C*x)

  ! covariance
  P = (1-K*C) * P
end subroutine rfx_lqg_kf_correct_diag_work
