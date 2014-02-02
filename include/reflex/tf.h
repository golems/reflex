/* -*- mode: C; c-basic-offset: 4 -*- */
/* ex: set shiftwidth=4 tabstop=4 expandtab: */
/*
 * Copyright (c) 2013, Georgia Tech Research Corporation
 * All rights reserved.
 *
 * Author(s): Neil T. Dantam <ntd@gatech.edu>
 * Georgia Tech Humanoid Robotics Lab
 * Under Direction of Prof. Mike Stilman <mstilman@cc.gatech.edu>
 *
 *
 * This file is provided under the following "BSD-style" License:
 *
 *
 *   Redistribution and use in source and binary forms, with or
 *   without modification, are permitted provided that the following
 *   conditions are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
 *   CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
 *   INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *   MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *   DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
 *   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 *   USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *   AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *   LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *   ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *   POSSIBILITY OF SUCH DAMAGE.
 *
 */

#ifndef REFLEX_TF_H
#define REFLEX_TF_H

#ifdef __cplusplus
extern "C" {
#endif //__cplusplus


/** Compute absolute tfs from relative tfs.
 *
 * Each frame must have a greater index than it's parent.
 * Parent value <0 indicates relative frame is in the global frame.
 */
int rfx_tf_abs( size_t n,
                const rfx_tf *tf_rel,
                ssize_t *parents,
                rfx_tf *tf_abs );

/** Compute jacobian for revolute links */
void rfx_tf_rev_jacobian( const double *AA_RESTRICT tf_abs, const double *axes,
                          size_t n, const size_t *AA_RESTRICT indices, const double *AA_RESTRICT pe,
                          double *AA_RESTRICT J, size_t ldJ );


/** Randomly corrupt a transform */
void rfx_tf_corrupt
( double theta_max, double x_max, const double e0[7], double e1[7] );

struct rfx_tf_filter {
    rfx_tf_dx X;  ///< state
    rfx_tf_dx Z;  ///< measurement
    rfx_tf_dx U;  ///< input

    double P[13*13];   ///< covariance
    double V[13*13];   ///< process noise
    double W[13*13];   ///< measurement noise
};


/** Update transform estimate.
 *
 * @pre: current measurement (Z) and input (U) written to F
 * @post: current state estimate (X) written to F
 */
int rfx_tf_filter_update( double dt, struct rfx_tf_filter *F);

int rfx_tf_filter_update_work
( double dt, double *XX, const double *UU, const double *ZZ, double *P, const double *V, const double *W );


struct rfx_lqg_duqu {
    double S[8];       ///< state
    double dx[6];      ///< state
    double P[14*14];   ///< covariance
    double V[14*14];   ///< process noise
    double W[14*14];   ///< measurement noise
};

int rfx_lqg_duqu_predict
( double dt, double *S, double *dS, double *P, const double *V );

int rfx_lqg_duqu_correct
( double dt, double *S_est, double *dS_est,
  const double *S_obs,
  double *P, const double *W );


void rfx_lqg_qutr_process_noise( double dt, double dx, double dtheta,
                                 double *E, double *V );


int rfx_lqg_qutr_process( void *cx, double *x, const double *u, double *F );
int rfx_lqg_qutr_measure( void *cx, const double *x, double *y, double *H );
int rfx_lqg_qutr_innovate( void *cx, const double *x, const double *z, double *y );
int rfx_lqg_qutr_update( void *cx, double *x, const double *Ky );

int rfx_lqg_qutr_predict
( double dt, double *E, double *dE, double *P, const double *V );

int rfx_lqg_qutr_correct
( double dt, double *E_est, double *dE_est,
  const double *E_obs,
  double *P, const double *W );


int rfx_tf_numeyama
( size_t n, double *_X, size_t ldx, double *_Y, size_t ldy, double tf[12] );
int rfx_tf_umeyama
( size_t n, const double *_X, size_t ldx, const double *_Y, size_t ldy, double tf[12] );

void rfx_tf_qlnmedian
( size_t n, const double *u, const double *Q, size_t ldq, double p[4] );

int rfx_tf_dud_median
( size_t n, const double *Ex, size_t ldx, const double *Ey, size_t ldy, double z[7] );

int rfx_tf_dud_rejected_mean
( size_t n, double z_theta, double z_x,
  const double *Ex, size_t ldx, const double *Ey, size_t ldy, double e[7] );

#ifdef __cplusplus
}
#endif //__cplusplus

#endif //REFLEX_TF_H
