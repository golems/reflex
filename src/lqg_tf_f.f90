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


function rfx_lqg_duqu_predict( dt, S, dS, P, V ) result(info) &
     bind( C, name="rfx_lqg_duqu_predict" )
  real(C_DOUBLE), intent(in), value :: dt
  real(C_DOUBLE), intent(inout) :: S(8), dS(8), P(16,16)
  real(C_DOUBLE), intent(in) :: V(16,16)
  integer(C_INT) :: info

  real(C_DOUBLE) :: A(16,16)
  real(C_DOUBLE) :: S_1(8)
  !real(C_DOUBLE) :: omega(8), omega_exp(8)
  integer(C_INT) :: i

  call aa_tf_duqu_sdiff( S, dS, dt, S_1 )

  ! Linearize
  ! X_1 = [ 1 dt*dS ] (S_0)
  !       [ 0     1 ] (dS)
  A = real(0,C_DOUBLE)
  forall (i=1:8)
     A(i,i) = real(1,C_DOUBLE)
     A(i,8+i) = dt
  end forall

  ! call aa_tf_duqu_diff2twist(S, dS, omega)
  ! omega = 0.5*dt*omega
  ! call aa_tf_duqu_exp(omega, omega_exp)

  ! call aa_tf_duqu_matrix_l( omega_exp, A(1:8,1:8) )

  forall (i=9:16)
     A(i,i) = real(1,C_DOUBLE)
  end forall

  call rfx_lqg_kf_predict_cov(int(16,C_SIZE_T),  A, V, P )

  ! store result
  info = 0
  S = S_1
end function rfx_lqg_duqu_predict

! function rfx_lqg_duqu_correct( dt, S_est, dx_est, S_obs, dx_obs, P, W ) result(info) &
!      bind( C, name="rfx_lqg_duqu_correct" )
!   real(C_DOUBLE), intent(in), value :: dt
!   real(C_DOUBLE), intent(inout) :: S_est(8), dx_est(6), P(14,14)
!   real(C_DOUBLE), intent(inout) :: S_obs(8), dx_obs(6)
!   real(C_DOUBLE), intent(in) :: W(14,14)
!   integer(C_INT) :: info

!   integer(C_INT) :: i

!   real (C_DOUBLE) :: y(14), S_rel(8), y_twist(8)
!   real(C_DOUBLE) :: H(14,14), K(14,14), Ky(14)
!   real(C_DOUBLE) :: S_1(8)


!   ! y = z-h(x)
!   !
!   ! Rather than computing y as a difference, get a dual quaternion
!   ! derivative as the log of the relative dual quaternion

!   call aa_tf_duqu_mulc( S_obs, S_est, S_rel )
!   call aa_tf_duqu_minimize( S_rel )
!   call aa_tf_duqu_ln( S_rel, y_twist)
!   call aa_tf_duqu_mul( y_twist, S_est, y(1:8) ) ! y is a duqu derivative
!   y(9:14) = dx_obs - dx_est

!   H = real(0,C_DOUBLE)
!   forall (i = 1:14)
!      H(i,i) = real(1,C_DOUBLE)
!   end forall

!   info = rfx_lqg_kf_correct_gain( int(14, C_SIZE_T), int(14, C_SIZE_T), H, P, W, K )

!   ! x = x + Ky
!   Ky = matmul(K,y)
!   print *, "Ky1", Ky(1:8)
!   print *, "Ky2", Ky(9:14)
!   print *, "K"
!   do i=1,14
!      print *,K(i,:)
!   end do
!   call aa_tf_duqu_sdiff(S_est, Ky(1:8), dt, S_1 )
!   S_est = S_1
!   call aa_tf_duqu_normalize(S_est)
!   dx_est  = dx_est + Ky(9:14)

!   ! P = (I - KH) P
!   call rfx_lqg_kf_correct_cov( int(14, C_SIZE_T), int(14,C_SIZE_T), H, P, K )
! end function rfx_lqg_duqu_correct



function rfx_lqg_duqu_correct( dt, S_est, dS_est, S_obs, P, W ) result(info) &
     bind( C, name="rfx_lqg_duqu_correct" )
  real(C_DOUBLE), intent(in), value :: dt
  real(C_DOUBLE), intent(inout) :: S_est(8), dS_est(8), P(16,16)
  real(C_DOUBLE), intent(inout) :: S_obs(8)
  real(C_DOUBLE), intent(in) :: W(8,8)
  integer(C_INT) :: info

  integer(C_INT) :: i

  real (C_DOUBLE) :: y(8), S_rel(8), y_twist(8)
  real(C_DOUBLE) :: H(8,16), K(16,8), Ky(16)
  real(C_DOUBLE) :: S_1(8)

  ! y = z-h(x)
  !
  ! Rather than computing y as a difference, get a dual quaternion
  ! derivative as the log of the relative dual quaternion

  call aa_tf_duqu_mulc( S_obs, S_est, S_rel )
  call aa_tf_duqu_minimize( S_rel )
  call aa_tf_duqu_ln( S_rel, y_twist)
  call aa_tf_duqu_mul( y_twist, S_est, y(1:8) ) ! y is a duqu derivative

  H = real(0,C_DOUBLE)
  forall (i = 1:8)
     H(i,i) = real(1,C_DOUBLE)
  end forall

  info = rfx_lqg_kf_correct_gain( int(16, C_SIZE_T), int(8, C_SIZE_T), H, P, W, K )

  ! x = x + Ky
  Ky = matmul(K,y)
  call aa_tf_duqu_sdiff(S_est, Ky(1:8), dt, S_1 )
  S_est = S_1
  !S_est = S_est + Ky(1:8)
  call aa_tf_duqu_normalize( S_est )

  dS_est = dS_est + Ky(9:16)

  ! P = (I - KH) P
  call rfx_lqg_kf_correct_cov( int(16, C_SIZE_T), int(8,C_SIZE_T), H, P, K )
end function rfx_lqg_duqu_correct


subroutine rfx_lqg_qutr_process_noise( dt, dtheta, dx, E, V ) &
     bind( C, name="rfx_lqg_qutr_process_noise" )
  real(C_DOUBLE), intent(in), value :: dt, dtheta, dx
  real(C_DOUBLE), intent(in) :: E(7)
  real(C_DOUBLE), intent(out) :: V(13,13)
  integer :: i
  real(C_DOUBLE) :: qs(4), Q_R(4,4)
  real(C_DOUBLE) :: qQ(4), vv(3), dvdv(3), dwdw(3)
  real(C_DOUBLE) :: vdv(3), Qdw(4,3)
  V = real(0,C_DOUBLE)

  QQ = dtheta*dt**4 / 4
  forall (i=1:4) V(i,i) = QQ(i)

  vv = dx*dt**4 / 4
  forall (i=1:3) V(i+4,i+4) = vv(i)

  dvdv = dx*dt**2
  forall (i=1:3) V(i+4+3,i+4+3) = dvdv(i)

  dwdw = dtheta*dt**2
  forall (i=1:3) V(i+4+3+3,i+4+3+3) = dwdw(i)

  vdv = dx*dt**2 / 3
  forall (i=1:3)
     V(i+4,i+7) = vdv(i)
     V(i+7,i+4) = vdv(i)
  end forall

  qs = (dt**2)/3*dtheta * E(1:4)
  call aa_tf_qmatrix_r(qs, Q_r)
  Qdw = Q_r(:,1:3)

  V(1:4,11:13) = Qdw
  V(11:13,1:4) = transpose(Qdw)

end subroutine rfx_lqg_qutr_process_noise

function rfx_lqg_qutr_predict( dt, E, dx, P, V ) result(info) &
     bind( C, name="rfx_lqg_qutr_predict" )
  real(C_DOUBLE), intent(in), value :: dt
  real(C_DOUBLE), intent(inout) :: E(7), dx(6), P(13,13)
  real(C_DOUBLE), intent(in) :: V(size(P,1),size(P,1))
  integer(C_INT) :: info

  real(C_DOUBLE) :: A(size(P,1),size(P,1))
  real(C_DOUBLE) :: E_1(7)
  real(C_DOUBLE) :: Q_R(4,4)
  integer(C_SIZE_T) :: i

  integer(C_SIZE_T) :: nx

  nx = size(P,1,C_SIZE_T)

  call aa_tf_qutr_svel( E, dx, dt, E_1 )
  call aa_tf_qnormalize( E(1:4) )
  !call aa_tf_qminimize( E(1:4) )


  ! Linearize
  ! X_1 = [ 1 dt*dS ] (S_0)
  !       [ 0     1 ] (dS)
  A = real(0,C_DOUBLE)

  ! call aa_tf_qdiff2vel(E(1:4), dE(1:4), omega)
  ! omega = 0.5*dt*omega
  ! call aa_tf_qexp(omega, omega_exp)
  ! call aa_tf_qmatrix_l( omega_exp, A(1:4,1:4) )


  forall (i=1:4)
     A(i,i) = real(1,C_DOUBLE)
  end forall
  call aa_tf_qmatrix_r( E(1:4), Q_R)
  A(1:4,11:13) = 0.5*dt*Q_R(:,1:3)

  do i=5,7
     A(i,i) = real(1,C_DOUBLE)
     A(i,7+i-3-1) = dt
  end do

  forall (i=7:nx)
     A(i,i) = real(1,C_DOUBLE)
  end forall



  !call rfx_lqg_kf_predict_cov(int(14,C_SIZE_T),  A, V, P )
  call rfx_lqg_kf_predict_cov(size(A,1,C_SIZE_T),  A, V, P )

  ! print *,"V"
  ! do i=1,nx
  !    print *,real(V(i,:))
  ! end do

  ! print *,"A"
  ! do i=1,nx
  !    print *,real(A(i,:))
  ! end do

  ! print *,"P_p"
  ! do i=1,nx
  !    print *,real(P(i,:))
  ! end do

  ! store result
  info = 0
  E = E_1
end function rfx_lqg_qutr_predict

function rfx_lqg_qutr_correct( dt, E_est, dx_est, E_obs, P, W ) result(info) &
     bind( C, name="rfx_lqg_qutr_correct" )
  real(C_DOUBLE), intent(in), value :: dt
  real(C_DOUBLE), intent(inout) :: E_est(7), dx_est(6), P(13,13)
  real(C_DOUBLE), intent(inout) :: E_obs(7)
  real(C_DOUBLE), intent(in) :: W(7,7)
  integer(C_INT) :: info

  integer(C_SIZE_T) :: i

  real (C_DOUBLE) :: y(7), E_rel(7), omega(4)
  real(C_DOUBLE) :: H(7,13), K(13,7), Ky(13)
  real(C_DOUBLE) :: E_1(7)
  integer(C_SIZE_T) :: nx, nz

  nx = size(P,1,C_SIZE_T)
  nz = size(W,1,C_SIZE_T)

  ! y = z-h(x)
  !
  ! Rather than computing y as a difference, get a dual quaternion
  ! derivative as the log of the relative dual quaternion

  call aa_tf_qmulc( E_obs(1:4), E_est(1:4), E_rel(1:4) )
  call aa_tf_qminimize( E_rel(1:4) )
  call aa_tf_qln( E_rel(1:4), omega)
  call aa_tf_qmul( omega, E_est(1:4), y(1:4) ) ! y is a quat derivative
  y(5:7) = E_obs(5:7) - E_est(5:7)

  H = real(0,C_DOUBLE)
  forall (i = 1:7)
     H(i,i) = real(1,C_DOUBLE)
  end forall

  info = rfx_lqg_kf_correct_gain( nx, nz, H, P, W, K )

  ! x = x + Ky
  Ky = matmul(K,y)

  call aa_tf_qutr_sdiff(E_est, Ky(1:7), dt, E_1 )
  E_est = E_1
  !E_est = E_est + Ky(1:7)
  !call aa_tf_qnormalize( E_est(1:4) )
  !call aa_tf_qminimize( E_est(1:4) )

  dx_est = dx_est + Ky(8:13)

  ! P = (I - KH) P
  call rfx_lqg_kf_correct_cov( nx, nz, H, P, K )

  ! print *,"info", info
  ! print *,"K"
  ! do i=1,nx
  !    print *,real(K(i,:))
  ! end do
  ! print *,"E_est", real(E_est)
  ! print *,"dx_est", real(dx_est)
  ! print *,"E_obs", real(E_obs)
  ! print *,"y", real(y)
  ! print *,"KyE", real(Ky(1:7))
  ! print *,"KydE", real(Ky(8:))

end function rfx_lqg_qutr_correct



! function rfx_lqg_kf_correct_gain(nx, nz, C, P, W, K) result(i) &
!      bind(C,name="rfx_lqg_kf_correct_gain")
!   use ISO_C_BINDING
!   integer(C_SIZE_T), intent(in), value :: nx, nz
!   real(C_DOUBLE), intent(in)  :: C(nz, nx), P(nx,nx), W(nz,nz)
!   real(C_DOUBLE), intent(out)  :: K(nx,nz)
!   integer(C_INT) :: i
!   real(C_DOUBLE) :: PC_T(nx,nz), Kp(nz,nz)
!   PC_T = matmul(P,transpose(C))
!   Kp = matmul(C, PC_T) + W
!   i = aa_la_inv(nz,Kp)
!   !K = real(0,C_DOUBLE)
!   K = matmul(PC_T,Kp)
!   print *,"info", i
! end function rfx_lqg_kf_correct_gain

! subroutine rfx_lqg_kf_correct_cov(nx, nz, C, P, K) &
!      bind(C,name="rfx_lqg_kf_correct_cov")
!   use ISO_C_BINDING
!   integer(C_SIZE_T), intent(in), value :: nx, nz
!   real(C_DOUBLE), intent(in)  :: C(nz, nx), K(nx,nz)
!   real(C_DOUBLE), intent(inout)  :: P(nx,nx)
!   real(C_DOUBLE) :: KC(nx,nx), P0(nx,nx)
!   integer(C_SIZE_T) :: i
!   KC = matmul(K,C)
!   KC = -KC
!   forall (i = 1:nx)
!      KC(i,i) = KC(i,i) + real(1,C_DOUBLE)
!   end forall
!   P0 = P
!   P = matmul(KC,P0)
! end subroutine rfx_lqg_kf_correct_cov





















! function rfx_lqg_duqu_process( dt, x, u, F ) result(info)
!   real(C_DOUBLE), intent(in) :: dt, u(14)
!   real(C_DOUBLE), intent(inout) :: x(14)
!   real(C_DOUBLE), intent(out) :: F(14,14)
!   real(C_DOUBLE) :: S(8)
!   real(C_DOUBLE) :: omega(8), omega_exp(8)
!   integer(C_INT) :: info
!   integer :: i

!   S = x(1:8)

!   call aa_tf_duqu_vel2twist(S, x(9:14), omega)
!   omega = 0.5*dt*omega
!   call aa_tf_duqu_exp(omega, omega_exp)
!   call aa_tf_duqu_mul(omega_exp, S, x(1:8))

!   ! Linearize
!   ! X_1 = [ [exp(omega*dt/2)] 0 ] (S_0)
!   !       [ 0                 1 ] (dx)
!   F = real(0,C_DOUBLE)
!   call aa_tf_duqu_matrix_l( omega_exp, F(1:8,1:8) )
!   forall (i=1:8)
!      F(8+i,8+i) = real(1,C_DOUBLE)
!   end forall

!   info = 0
! end function rfx_lqg_duqu_process

! function rfx_lqg_duqu_measure( dt, x, y, H ) result(info)
!   real(C_DOUBLE), intent(in) :: dt, x(14)
!   real(C_DOUBLE), intent(in) :: y(14), H(14)
!   integer i;

!   H = real(0,C_DOUBLE)
!   forall (i = 1:14)
!      H(i,i) = real(1,C_DOUBLE)
!   end forall

!   y = x
! end function rfx_lqg_duqu_measure
