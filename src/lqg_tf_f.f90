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


 !! interface to inversion function (via lapack) in amino
  ! Interface
  !    Function aa_la_inv( n, A ) result(info)
  !      use ISO_C_BINDING
  !      integer(C_SIZE_T), intent(in), value :: n
  !      real(C_DOUBLE), intent(inout) :: A(n,n)
  !      integer :: info
  !    End Function aa_la_inv
  ! End Interface

#define TF_R 1:4
#define TF_V 5:7
#define TF_DV 8:10
#define TF_W 11:13

subroutine rfx_tf_filter_update_work( dt, XX, UU, ZZ, P, V, W ) &
     bind( C, name="rfx_tf_filter_update_work" )
  ! type :: rfx_tf_state
  !    real(C_DOUBLE) :: r(4)
  !    real(C_DOUBLE) :: v(3)
  !    real(C_DOUBLE) :: dv(3)
  !    real(C_DOUBLE) :: w(3)
  ! end type rfx_tf_state
  real(C_DOUBLE), intent(in), value :: dt
  real(C_DOUBLE), intent(inout), dimension(13) :: XX
  real(C_DOUBLE), intent(in), dimension(13) :: ZZ, UU
  real(C_DOUBLE), intent(inout), dimension(13,13) :: P
  real(C_DOUBLE), intent(in), dimension(13,13) :: V, W

  real(C_DOUBLE), dimension(13) :: dx, zh, dxh, XX2
  real(C_DOUBLE), dimension(13,13) ::  A, C, K
  real(C_DOUBLE) :: r(4), R_r(4,4)
  ! type(rfx_tf_state) :: Xs, Zs, Hs
  integer :: i

  ! call extract(XX, Xs)
  ! call extract(ZZ, Zs)

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


  !print *, "A"
  !print *, A
  !print *, "C"
  !print *, C

  !print *, "call kbf gain"
  call rfx_lqg_kbf_gain_work( int(13,C_SIZE_T), int(13,C_SIZE_T), A, C, V, W, P, K )
  !print *, "got kbf gain", K

  !! Update

  !dx = Ax+Bu
  call aa_tf_qvel2diff( XX(TF_R), XX(TF_W), dx(TF_R) )
  dx(TF_V) = XX(TF_DV)
  dx(TF_DV) = real(0,C_DOUBLE)
  dx(TF_W) = real(0,C_DOUBLE)
  dx = dx+UU

  ! innovate
  call innovate_r( XX(TF_R), ZZ(TF_R), zh(TF_R) )
  zh(TF_V) = ZZ(TF_V) - XX(TF_V)
  zh(TF_DV) = ZZ(TF_DV) - XX(TF_DV)
  zh(TF_W) = ZZ(TF_W) - XX(TF_W)
  dxh = matmul(K, zh)

  dx = dx + dxh

  !! Integrate
  !! TODO: better integration (consider the twist)
  call aa_tf_qsdiff( XX(TF_R), dx(TF_R), dt, XX2(TF_R) )
  XX2(TF_V) = XX(TF_V) + dt*dx(TF_V)
  XX2(TF_DV) = XX(TF_DV) + dt*dx(TF_DV)
  XX2(TF_W) = XX(TF_W) + dt*dx(TF_W)

  XX = XX2
contains

  subroutine innovate_r( r_est, r_obs, rh )
    real(C_DOUBLE), intent(in), dimension(4) :: r_est, r_obs
    real(C_DOUBLE), intent(out), dimension(4) :: rh
    real(C_DOUBLE), dimension(4) :: a, b, c
    call aa_tf_qmulc( r_obs, r_est, a ) ! a = r_est * r_obs
    call aa_tf_qln( a, b)              ! b = log( r_est * r_obs )
    c = b / 2                          ! c = 1/2 * log( r_est * r_obs )
    call aa_tf_qmul(c, r_est, rh )     ! c = 1/2 * log( r_est * r_obs ) * r_est
  end subroutine innovate_r

  ! subroutine extract(x,s)
  !   real(C_DOUBLE) :: x(13)
  !   type(rfx_tf_state) :: s
  !   s%r = x(1:4)
  !   s%v = x(5:7)
  !   s%dv = x(8:10)
  !   s%w = x(11:13)
  ! end subroutine extract

  ! subroutine insert(s,x)
  !   real(C_DOUBLE) :: x(13)
  !   type(rfx_tf_state) :: s
  !   x(1:4) =  s%r
  !   x(5:7) =  s%v
  !   x(8:10) = s%dv
  !   x(11:13) = s%w
  ! end subroutine extract

end subroutine rfx_tf_filter_update_work



