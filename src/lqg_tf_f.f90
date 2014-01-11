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


#define TF_R 1:4
#define TF_V 5:7
#define TF_DV 8:10
#define TF_W 11:13

#define TF_RV 1:7
#define TF_DX 8:13

function rfx_tf_filter_update_work( dt, XX, UU, ZZ, P, V, W ) result(info) &
     bind( C, name="rfx_tf_filter_update_work" )
  real(C_DOUBLE), intent(in), value :: dt
  real(C_DOUBLE), intent(inout), dimension(13) :: XX
  real(C_DOUBLE), intent(in), dimension(13) :: ZZ, UU
  real(C_DOUBLE), intent(inout), dimension(13,13) :: P
  real(C_DOUBLE), intent(in), dimension(13,13) :: V, W
  integer(C_INT) :: info

  real(C_DOUBLE), dimension(13) :: dx, zh, dxh, XX2
  real(C_DOUBLE), dimension(13,13) ::  A, C, K
  real(C_DOUBLE) :: r(4), R_r(4,4)
  integer :: i


  !! Compute Kalman Gain
  ! A = [0   0   0   0.5R_r ]
  !     [0   0   1   0 ]
  !     [0   0   0   0 ]
  !     [0   0   0   0 ]
  A = real(0,C_DOUBLE)

  r = 0.5 * XX(TF_R)
  call aa_tf_qmatrix_r(r, R_r)
  A(TF_R,TF_W) = R_r(1:4,1:3)

  forall (i=1:3)
     A(4+i,7+i) = real(1,C_DOUBLE)
  end forall

  ! C = 1
  C=real(0,C_DOUBLE)
  forall (i=1:13)
     C(i,i) = real(1,C_DOUBLE)
  end forall


  call rfx_lqg_kbf_gain_work( int(13,C_SIZE_T), int(13,C_SIZE_T), A, C, V, W, P, K )

  !! Update

  !dx = Ax+Bu
  call aa_tf_qvel2diff( XX(TF_R), XX(TF_W), dx(TF_R) )
  dx(TF_V) = XX(TF_DV)
  dx(TF_DX) = real(0,C_DOUBLE)
  dx = dx+UU

  ! innovate: e = z - Cx
  call innovate_r( XX(TF_R), ZZ(TF_R), zh(TF_R) )
  zh(TF_V) = ZZ(TF_V) - XX(TF_V)
  zh(TF_DX) = ZZ(TF_DX) - XX(TF_DX)
  dxh = matmul(K, zh)

  dx = dx + dxh

  !! Integrate
  call aa_tf_qutr_sdiff( XX(TF_RV), dx(TF_RV), dt, XX2(TF_RV) )
  XX2(TF_DX) = XX(TF_DX) + dt*dx(TF_DX)

  XX = XX2
  call aa_tf_qnormalize( XX(TF_R) )

  !! TODO: actually check that things worked out
  info = 0
contains

  subroutine innovate_r( r_est, r_obs, rh )
    real(C_DOUBLE), intent(in), dimension(4) :: r_est, r_obs
    real(C_DOUBLE), intent(out), dimension(4) :: rh
    real(C_DOUBLE), dimension(4) :: r_rel, w_2
    ! rh = log( r_obs * conj(r_est) ) * r_est
    call aa_tf_qmulc( r_obs, r_est, r_rel )
    call aa_tf_qln( r_rel, w_2)
    call aa_tf_qmul( w_2, r_est, rh )
  end subroutine innovate_r

end function rfx_tf_filter_update_work


function rfx_lqg_duqu_predict( dt, S, dx, P, V ) result(info) &
     bind( C, name="rfx_lqg_duqu_predict" )
  real(C_DOUBLE), intent(in), value :: dt
  real(C_DOUBLE), intent(inout) :: S(8), dx(6), P(14,14)
  real(C_DOUBLE), intent(in) :: V(14,14)
  integer(C_INT) :: info

  real(C_DOUBLE) :: A(14,14), AP(14,14)
  real(C_DOUBLE) :: omega(8), omega_exp(8), S_1(8)
  integer(C_INT) :: i

  ! X = [S dx]^T

  ! update state: x_1 = f(x_0)
  ! S_1 = exp(omega*dt/2) * S_0
  ! dx_1 = dx_0
  call aa_tf_duqu_vel2twist(S, dx, omega)
  omega = 0.5*dt*omega
  call aa_tf_duqu_exp(omega, omega_exp)
  call aa_tf_duqu_mul(omega_exp, S, S_1)

  ! Linearize
  ! X_1 = [ [exp(omega*dt/2)] 0 ] (S_0)
  !       [ 0                 1 ] (dx)
  A = real(0,C_DOUBLE)
  call aa_tf_duqu_matrix_l( omega_exp, A(1:8,1:8) )
  forall (i=1:8)
     A(8+i,8+i) = real(1,C_DOUBLE)
  end forall

  ! update covariance: P = F_0 P_0 F_0^T + V_0
  AP = matmul(A,P)
  P = matmul(AP,transpose(A)) + V

  ! store result
  info = 0
  S = S_1
end function rfx_lqg_duqu_predict

function rfx_lqg_duqu_correct( dt, S_est, dx_est, S_obs, dx_obs, P, W ) result(info) &
     bind( C, name="rfx_lqg_duqu_correct" )
  real(C_DOUBLE), intent(in), value :: dt
  real(C_DOUBLE), intent(inout) :: S_est(8), dx_est(6), P(14,14)
  real(C_DOUBLE), intent(inout) :: S_obs(8), dx_obs(6)
  real(C_DOUBLE), intent(in) :: W(14,14)
  integer(C_INT) :: info

  integer(C_INT) :: i

  real (C_DOUBLE) :: y(14), S_rel(8), y_twist(8)
  real(C_DOUBLE) :: H(14,14), K(14,14), Ky(14)
  real(C_DOUBLE) :: PHT(14,14), E(14,14)
  real(C_DOUBLE) :: S_1(8)
  real(C_DOUBLE) :: KH(14,14), P_1(14,14)


  ! y = z-h(x)
  !
  ! Rather than computing y as a difference, get a dual quaternion
  ! derivative as the log of the relative dual quaternion

  call aa_tf_duqu_mulc( S_obs, S_est, S_rel )
  call aa_tf_duqu_minimize( S_rel )
  call aa_tf_duqu_ln( S_rel, y_twist)
  call aa_tf_duqu_mul( y_twist, S_est, y(1:8) ) ! y is a duqu derivative
  y(9:14) = dx_obs - dx_est

  ! H = [ [ln(S_obs S_est^*)]   0 ]
  !     [            0          1 ]
  H = real(0,C_DOUBLE)
  !call aa_tf_duqu_matrix_l( y_twist, H(:,1:8), int(14, C_SIZE_T) )
  ! forall (i = 1:6)
  !    H(8+i,8+i) = real(1,C_DOUBLE)
  ! end forall
  forall (i = 1:14)
     H(i,i) = real(1,C_DOUBLE)
  end forall

  ! E = HPH^T + W
  PHT = matmul(P,transpose(H))
  E = matmul(H,PHT) + W

  ! K = PH^TE^{-1}
  info = aa_la_inv(int(14,C_SIZE_T), E)
  K = matmul(PHT, E)

  ! x = x + Ky
  Ky = matmul(K,y)
  call aa_tf_duqu_sdiff(S_est, Ky(1:8), dt, S_1 )
  S_est = S_1
  call aa_tf_duqu_normalize(S_est)
  dx_est  = dx_est + y(9:14)

  ! P = (I - KH) P
  KH = matmul(K,H)
  KH = -KH
  forall(i=1:14)
     KH(i,i) = KH(i,i)+1
  end forall
  P_1 = matmul(KH,P)
  P = P_1
end function rfx_lqg_duqu_correct

