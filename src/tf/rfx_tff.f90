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
!!!!!!!!!!!!!!!
!!! UMEYAMA !!!
!!!!!!!!!!!!!!!

function rfx_tf_numeyama( n, X, ldx, Y, ldy, tf ) result(info) &
     bind(C,name="rfx_tf_numeyama")
  integer(C_SIZE_T), intent(in), value :: n, ldx, ldy
  real(C_DOUBLE), intent(inout) :: X(ldx,n), Y(ldy,n)
  real(C_DOUBLE), intent(out) :: tf(3,4)
  integer(C_INT) :: info
  real(C_DOUBLE) :: U(3,3), Sd(3,3), S(3), Vt(3,3), ux(3), uy(3), sigma(3,3), tmp(3,3), tmpv(3)
  integer(C_SIZE_T) :: i, m
  !! de-mean
  call aa_la_colmean( X(1:3,:), ux )
  call aa_la_colmean( Y(1:3,:), uy )
  forall (i=1:n)
     X(1:3,i) = X(1:3,i) - ux
     Y(1:3,i) = Y(1:3,i) - uy
  end forall

  sigma = 1.0/real(n,C_DOUBLE) * matmul(Y(1:3,:), transpose(X(1:3,:)))

  m = int(3,C_SIZE_T)
  call aa_la_svd(m,m, sigma,m, U,m, S, Vt,m )

  Sd = real(0.0,C_DOUBLE)
  forall (i = 1:3) Sd(i,i) = real(1,C_DOUBLE)
  if( aa_la_det(sigma) < real(0,C_DOUBLE) ) Sd(3,3) = real(-1.0,C_DOUBLE)
  !! TODO: rank defficient


  tmp = matmul(U,Sd)
  tf(:,1:3) = matmul(U, Vt)
  !tf(:,1:3) = tf(:,1:3)

  tmpv = matmul(tf(:,1:3),ux)
  tf(:,4) = uy - tmpv

  info = 0
end function rfx_tf_numeyama



