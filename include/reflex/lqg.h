/* -*- mode: C; c-basic-offset: 4 -*- */
/* ex: set shiftwidth=4 tabstop=4 expandtab: */
/*
 * Copyright (c) 2010-2011, Georgia Tech Research Corporation
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

#ifndef REFLEX_LQG_H
#define REFLEX_LQG_H

#ifdef __cplusplus
extern "C" {
#endif

/** @file lqg.h
 *  @author Neil T. Dantam
 */

/*********************************/
/* Basic Kalman Filter Functions */
/*********************************/

/** Update Kalman Filter covariance.
 *
 * \f[ P = A*P*AR^ + V \f]
 */
void rfx_lqg_kf_predict_cov
( size_t n, const double *A, const double *V, double *P );

/** Compute kalman gain matrix.
 *
 * \f[ K = P * C^T * (C * P * C^T + W)^{-1} \f]
 *
 * @param n_x size of state space
 * @param n_z size of measurement space
 * @param C Measurement model, n_z*n_x
 * @param P Coviariance, n_x*n_x
 * @param W Measurement noise model, n_z*n_z
 * @param K Kalman Gain, n_x*n_z
 */
int rfx_lqg_kf_correct_gain
( size_t n_x, size_t n_z, const double *C, const double *P, const double *W, double *K );

/** Update kalman filter covariance.
 *
 * \f[ P = (I - K*C) * P \f]
 */
void rfx_lqg_kf_correct_cov
( size_t n_x, size_t n_z, const double *C, double *P, double *K );


/**************************/
/* Extended Kalman Filter */
/**************************/

/* The next few typedefs define functions for computing state updates,
 * expected measurements, and linearizing the system */

/** Function to compute process update and linearization.
 *
 * \f[ x \approx F x \f]
 *
 * @param cx context struct
 * @param x state vector
 * @param u input vector
 * @param F linearized process model
 *
 * @pre
 *  - x contains the prior state estimate
 * @post
 *  - x contains the predicted state estimate
 *  - F contains the linearized system model
 */
typedef int (*rfx_lqg_ekf_process_fun)( void *cx, double *x, const double *u, double *F );

/** Function to compute innovation and linearization.
 *
 * \f[ y = h(x) \f]
 *
 * \f[ y \approx Hx \f]
 *
 * @param cx context struct
 * @param x state vector
 * @param y predicted measurement
 * @param F linearized measurement model
 *
 * @post
 *  - y contains the predicted measurement for state x
 *  - H contains the linearized measurement model at x
 *
 */
typedef int (*rfx_lqg_ekf_measure_fun)( void *cx, const double *x, double *y, double *H );

/** Function to compute innovation.
 *
 * \f[ y \approx z - y \f]
 *
 * @param cx context struct
 * @param x state vector
 * @param z actual measurement
 * @param y predicted measurement, innovation
 *
 * @post
 *  - y contains the predicted measurement for state x
 *  - y contains the innovation
 */
typedef int (*rfx_lqg_ekf_innovate_fun)( void *cx, const double *x, const double *z, double *y );

/** Function to compute state update.
 *
 * @param cx context struct
 * @param x state vector
 * @param Ky the product of the Kalman gain and measurement innovation
 *
 * @pre
 *  - x contains the estimated state
 *
 * @post
 *  - x contains the corrected state estimate
 *
 */
typedef int (*rfx_lqg_ekf_update_fun)( void *cx, double *x, const double *Ky );

/** Driver for Extended Kalman Filter prediction step.
 *
 * Compute the updated state estimate x and covariance P
 *
 * \f[ x = f(x,u) \f]
 * \f[ x \approx Fx \f]
 *
 * \f[ P = FPF^T + V \f]
 *
 * @param cx context struct
 * @param n_x size of state vector
 * @param x state vector
 * @param u input vector (size implicitly known by process())
 * @param P covariance
 * @param V process noise
 * @param process The process model function
 *
 * @pre
 *  - x contains the estimated state
 *  - P contains the process covariance
 *
 * @post
 *  - x contains the predicted state estimate
 *  - P contains the updated process covariance
 *
 */
int rfx_lqg_ekf_predict
( void *cx, size_t n_x, double *x, const double *u, double *P, const double *V,
  rfx_lqg_ekf_process_fun process );

/** Driver for Extended Kalman Filter correction step.
 *
 * Compute the updated state estimate x and covariance P
 *
 * \f[ y = z - h(x) \f]
 * \f[ y \approx z - Hx \f]
 *
 * \f[ S = HPH^T + R \f]
 * \f[ K = PH^TS^{-1} \f]
 *
 * \f[ x := x + Ky \f]
 * \f[ P := (I - KH)P \f]
 *
 * @param cx context struct
 * @param n_x size of state vector
 * @param n_z size of measurement vector
 * @param x state vector
 * @param z measurement vector
 * @param P covariance
 * @param W measurement noise
 * @param measure The measurement model function
 * @param innovate Function to compute the innovation
 * @param update Function to update the state estimate
 *
 * @pre
 *  - x contains the estimated state
 *  - P contains the process covariance
 *
 * @post
 *  - x contains the corrected state estimate
 *  - P contains the updated process covariance
 */
int rfx_lqg_ekf_correct
( void *cx, size_t n_x, size_t n_z, double *x, const double *z, double *P, const double *W,
  rfx_lqg_ekf_measure_fun measure, rfx_lqg_ekf_innovate_fun innovate, rfx_lqg_ekf_update_fun update );




/****************************************/
/* Linear Quadratic Gaussian Controller */
/****************************************/


typedef struct {
    size_t n_x;  ///< state size
    size_t n_u;  ///< input size
    size_t n_z;  ///< output size

    double *x;   ///< state estimate,         vector size n_x
    double *u;   ///< computed input,         vector size n_u
    double *z;   ///< measurement,            vector size n_z

    double *A;   ///< process model,          matrix size n_x * n_x, column-major
    double *B;   ///< input model,            matrix size n_x * n_u, column-major
    double *C;   ///< measurement model,      matrix size n_z * n_x, column-major

    double *P;   ///< covariance,             matrix size n_x * n_x, column-major
    double *V;   ///< process noise,          matrix size n_x * n_x, column-major
    double *W;   ///< measurement noise,      matrix size n_z * n_z, column-major

    double *Q;   ///< state cost,             matrix size n_x * n_x, column-major
    double *R;   ///< input cost,             matrix size n_u * n_u, column-major

    double *K;   ///< optimal feedback gain,  matrix size n_x * n_z, column-major
    double *L;   ///< optimal control gain,   matrix size n_u * n_x, column-major

    aa_mem_region_t *reg;
} rfx_lqg_t;


/** Initialize the LQG struct. */
AA_API void rfx_lqg_init( rfx_lqg_t *lqg, size_t n_x, size_t n_u, size_t n_z );
/** Free members of the LQG struct. */
AA_API void rfx_lqg_destroy( rfx_lqg_t *lqg );

/** Kalman-Bucy optimal gain.
 *
 * \f[ \dot{P} = AP + PA^T - PC^TW^{-1}CP + V = 0 \f]
 * \f[ K \leftarrow PC^TW^{-1} \f]
 *
 * Computes optimal K by solving the algebraic Riccati equation.
 *
 * \pre
 *   - lqg.n_x contains the state-space size
 *   - lqg.n_z contains the state-space size
 *   - lqg.A contains the process model matrix, size n_x*n_x
 *   - lqg.C contains the measurement model matrix, size n_z*n_x
 *   - lqg.V contains the process noise model matrix, size n_x*n_x
 *   - lqg.W contains the measurement noise model matrix, size n_z*n_z
 *   - lqg.P contains space for the covariance matrix, size n_x*n_x
 *   - lqg.K contains space for the kalman gain matrix, size n_x*n_z
 * \post
 *   - lqg.K overwritten with the kalman gain
 *   - lqg.P overwritten with the covariance
 */
AA_API void rfx_lqg_kbf_gain( rfx_lqg_t *lqg );

/** Kalman-Bucy optimal gain.
 *
 * \f[ \dot{P} = AP + PA^T - PC^TW^{-1}CP + V = 0 \f]
 * \f[ K \leftarrow PC^TW^{-1} \f]
 *
 * Computes optimal K by solving the algebraic Riccati equation.
 *
 * \pre
 *   - n_x contains the state-space size
 *   - n_z contains the state-space size
 *   - A contains the process model matrix, size n_x*n_x
 *   - C contains the measurement model matrix, size n_z*n_x
 *   - V contains the process noise model matrix, size n_x*n_x
 *   - W contains the measurement noise model matrix, size n_z*n_z
 *   - P contains space for the covariance matrix, size n_x*n_x
 *   - K contains space for the kalman gain matrix, size n_x*n_z
 * \post
 *   - K overwritten with the kalman gain
 *   - P overwritten with the covariance
 */
void rfx_lqg_kbf_gain_work
( size_t n_x, size_t n_z,
  const double *A, const double *C, const double *V, const double *W, double *P, double *K );

/** Kalman-Bucy filter, Euler-integrate
 *
 * \f[ \dot{x} = Ax + Bu + K(z-Cx) \f]
 * \f[ x \leftarrow \Delta t \dot{x} + x \f]
 *
 * \pre
 *   - lqg.n_x contains the state-space size
 *   - lqg.n_u contains the input-space size
 *   - lqg.n_z contains the state-space size
 *   - lqg.A contains the process model matrix, size n_x*n_x
 *   - lqg.B contains the input model matrix, size n_x*n_u
 *   - lqg.C contains the measurement model matrix, size n_z*n_x
 *   - lqg.K contains the kalman gain matrix, size n_x*n_z
 *   - lqg.x contains the prior state vector, size n_x
 *   - lqg.u contains the current input vector, size n_u
 *   - lqg.z contains the current measurement vector, size n_z
 * \post
 *   - lqg.x overwritten with the next state vector
 */
AA_API void rfx_lqg_kbf_step1( rfx_lqg_t *lqg, double dt );

/** Kalman-Bucy filter, Runge-Kutta-4  integrate
 *
 * \f[ \dot{x} = Ax + Bu + K(z-Cx) \f]
 *
 * \pre
 *   - lqg.n_x contains the state-space size
 *   - lqg.n_u contains the input-space size
 *   - lqg.n_z contains the state-space size
 *   - lqg.A contains the process model matrix, size n_x*n_x
 *   - lqg.B contains the input model matrix, size n_x*n_u
 *   - lqg.C contains the measurement model matrix, size n_z*n_x
 *   - lqg.K contains the kalman gain matrix, size n_x*n_z
 *   - lqg.x contains the prior state vector, size n_x
 *   - lqg.u contains the current input vector, size n_u
 *   - lqg.z contains the current measurement vector, size n_z
 * \post
 *   - lqg.x overwritten with the next state vector
 */
AA_API void rfx_lqg_kbf_step4( rfx_lqg_t *lqg, double dt );

/** Discrete time Kalman filter predict step.
 *
 * \f[ x \leftarrow Ax + Bu \f]
 * \f[ P \leftarrow A  P  A^T + V \f]
 *
 * \pre
 *   - lqg.n_x contains the state-space size
 *   - lqg.n_u contains the input-space size
 *   - lqg.A contains the process model matrix, size n_x*n_x
 *   - lqg.B contains the input model matrix, size n_x*n_u
 *   - lqg.V contains the process noise model matrix, size n_x*n_x
 *   - lqg.P contains the covariance matrix, size n_x*n_x
 *   - lqg.x contains prior state vector, size n_x
 *   - lqg.u contains input vector, size n_x
 *   - The stack can grow by at least an additional sizeof(double)*n_x*n_x bytes.
 *
 * \post
 *   - lqg.x overwritten with the predicted state vector
 *   - lqg.P overwritten with the updated covariance matrix
 */
AA_API void rfx_lqg_kf_predict( rfx_lqg_t *lqg );

/** Discrete time Kalman filter correct step.
 *
 * \f[ K \leftarrow P  C^T  (C  P  C^T + W)^{-1} \f]
 * \f[ x \leftarrow x + K  (z - Cx) \f]
 * \f[ P \leftarrow (I - KC)  P \f]
 *
 * \pre
 *   - lqg.n_x contains the state-space size
 *   - lqg.n_u contains the input-space size
 *   - lqg.C contains the measurement model
 *   - lqg.P contains the covariance matrix, size n_x*n_x
 *   - lqg.W contains the measurement noise model, size n_z*n_z
 *   - lqg.x contains the predicted state vector, size n_x
 *   - lqg.z contains the measurement vector, size n_z
 *   - lqg.K is an n_x*n_z matrix
 *   - The stack can grow by at least an additional sizeof(double) * MAX( 2*n_x*n_x, n_z*n_z + n_z*n_x ) bytes.
 *
 * \post
 *   - lqg.x overwritten with the corrected state vector
 *   - lqg.K overwritten with the Kalman gain matrix
 *   - lqg.P overwritten with the updated covariance matrix
 */
AA_API void rfx_lqg_kf_correct( rfx_lqg_t *lqg );


AA_API void rfx_lqg_lqr_gain( rfx_lqg_t *lqg );

AA_API void rfx_lqg_lqr_ctrl( rfx_lqg_t *lqg );

AA_API void rfx_lqg_sys( const void *lqg,
                         double t, const double *AA_RESTRICT x,
                         double *AA_RESTRICT dx );



#ifdef __cplusplus
}
#endif

#endif // REFLEX_LQG_H
