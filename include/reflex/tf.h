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

/**
 * Descriptor for transform operations.
 *
 * These descriptors are produced by the frame code generator.
 */
struct rfx_tf_desc {

    /**
     * Number of configurations
     */
    size_t n_config;

    /**
     * Number of frames
     */
    size_t n_frame;

    /**
     * Compute absolute transforms.
     *
     * @param[in]  q       the configuration variables
     * @param[in]  incQ    increment of q array
     * @param[out] E_abs   the absolute frames in quaternion-translation form
     * @param[in]  ldAbs   leading dimension of E_abs
     */
    void (*frames)( const double *AA_RESTRICT q, size_t incQ,
                    double * AA_RESTRICT E_abs, size_t ldAbs );

    /**
     * Array of axes
     */
    const double *config_axes;

    /**
     * Array of frame names
     */
    const char **frame_name;

    /**
     * Array of frame names
     */
    const char **config_name;

    /**
     * Array of frame parent indices
     */
    const ssize_t *frame_parent;
};

/** Compute absolute tfs from relative tfs.
 *
 * Each frame must have a greater index than it's parent.
 * Parent value <0 indicates relative frame is in the global frame.
 */
int rfx_tf_abs( size_t n,
                const rfx_tf *tf_rel,
                ssize_t *parents,
                rfx_tf *tf_abs );

/** Compute Jacobian for kinematic chain with revolute joints
 *
 * \param[in] tf_abs Array of absolute frames in quaternion-translation form
 * \param[in] ldT Leading dimension of T (typically 7)
 * \param[in] axes Array of joint axes
 * \param[in] indices_tf Array of indices for frames
 * \param[in] indices_axis Array of indices for joint axes
 * \param[in] pe translational position of the end-effector
 * \param[out] J Jacobian matrix
 * \param[in] ldJ Leading dimension of Jacobian matrix (typically 6)
 */
void rfx_tf_rev_jacobian( const double *AA_RESTRICT tf_abs, size_t ldT,
                          const double *axes,
                          size_t n, const size_t *AA_RESTRICT indices_tf,  const size_t *AA_RESTRICT indices_axis,
                          const double *AA_RESTRICT pe,
                          double *AA_RESTRICT J, size_t ldJ );

/* Generate random tf */
void rfx_tf_rand( double theta_max, double x_max, double e[7] );

/** Randomly corrupt a transform */
void rfx_tf_corrupt
( double theta_max, double x_max, const double e0[7], double e1[7] );


enum rfx_tf_cor_opts {
    RFX_TF_COR_O_TRANS_MEDIAN = 0x1,
    RFX_TF_COR_O_TRANS_MEAN = 0x2,
    RFX_TF_COR_O_ROT_UMEYAMA = 0x4,
    RFX_TF_COR_O_ROT_DAVENPORT = 0x8,
    RFX_TF_COR_O_ROT_MEDIAN = 0x10
};

/**
 * Compute fit between corresponding transforms
 *
 * @param[in] opts type of fit to perform (bitmask of enum rfx_tf_cor_opts)
 * @param[in] qx quaternion array 0
 * @param[in] ldqx leading dimensions of qx
 * @param[in] qy quaternion array 1
 * @param[in] ldqy leading dimensions of qy
 * @param[in] vx translation array 0
 * @param[in] ldvx leading dimensions of vx
 * @param[in] vy translation array 1
 * @param[in] ldvy leading dimensions of vy
 * @param[out] Z output transform (quaternion, translation)
 */

void rfx_tf_cor( int opts, size_t n,
                 const double *qx, size_t ldqx,
                 const double *vx, size_t ldvx,
                 const double *qy, size_t ldqy,
                 const double *vy, size_t ldvy,
                 double *Z );


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
int rfx_tf_filter_update( double dt, struct rfx_tf_filter *F) AA_DEPRECATED;

int rfx_tf_filter_update_work
( double dt, double *XX, const double *UU, const double *ZZ, double *P, const double *V, const double *W ) AA_DEPRECATED;


struct rfx_lqg_duqu {
    double S[8];       ///< state
    double dx[6];      ///< state
    double P[14*14];   ///< covariance
    double V[14*14];   ///< process noise
    double W[14*14];   ///< measurement noise
};

int rfx_lqg_duqu_predict
( double dt, double *S, double *dS, double *P, const double *V ) AA_DEPRECATED;

int rfx_lqg_duqu_correct
( double dt, double *S_est, double *dS_est,
  const double *S_obs,
  double *P, const double *W ) AA_DEPRECATED;


void rfx_lqg_qutr_process_noise( double dt, double dx, double dtheta,
                                 double *E, double *V );


int rfx_lqg_qutr_process( void *cx, double *x, const double *u, double *F );
int rfx_lqg_qutr_measure( void *cx, const double *x, double *y, double *H );
int rfx_lqg_qutr_innovate( void *cx, const double *x, const double *z, double *y );
int rfx_lqg_qutr_update( void *cx, double *x, const double *Ky );

/**
 * Kalman filter prediction step for quaternion-translation pose
 *
 * @param dt time step
 * @param E pose (quaternion, translation)
 * @param dE pose derivatice (quaternion-derivative, translation-derivative)
 * @param P convarance (13x13)
 * @param V process noise (13x13)
 */
int rfx_lqg_qutr_predict
( double dt, double *E, double *dE, double *P, const double *V );

/**
 * Kalman filter correction step for quaternion-translation pose
 *
 * @param dt time step
 * @param E_est pose estimate
 * @param dE_est pose derivative estimate
 * @param dE_obs pose measurement
 * @param P covariance
 * @param W measurement noise (7x7)
 *
 */
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

void rfx_tf_qangmedian
( size_t n, const double *Q, size_t ldq, double p[4] );

int rfx_tf_dud_median
( size_t n, const double *Ex, size_t ldx, const double *Ey, size_t ldy, double z[7] ) AA_DEPRECATED;

void rfx_tf_median
( size_t n, const double *E, size_t lde, double *z );

int rfx_tf_dud_rejected_mean
( size_t n, double z_theta, double z_x,
  const double *Ex, size_t ldx, const double *Ey, size_t ldy, double e[7] ) AA_DEPRECATED;

/* MAD Filtering */
/* Use Median Absolute Deviation to reject outliers */

int rfx_tf_madqg_predict
( double dt, double *E, double *dE, double *P, const double *V );

int rfx_tf_madqg_correct
( double dt,
  size_t max_delta, double *delta_theta, double *delta_x, size_t *n_delta, size_t *i_delta,
  double *E_est, double *dE_est,
  size_t n_obs, const double *E_obs,
  double *P, const double *W );


int rfx_tf_madqg_correct2
( double dt,
  size_t max_delta, double *delta_theta, double *delta_x, size_t *n_delta, size_t *i_delta,
  double *E_est, double *dE_est,
  size_t n_obs, const double *bEo, const double *cEo,
  double *P, const double *W );

int rfx_tf_madqg_correct_median_window
( double dt,
  size_t max_hist, double *E_obs_hist, size_t *n_hist, size_t *i_hist,
  double *E_est, double *dx_est,
  size_t n_obs, const double *E_obs,
  double *P, const double *W );

#ifdef __cplusplus
}
#endif //__cplusplus

#endif //REFLEX_TF_H
