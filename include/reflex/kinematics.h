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

#ifndef REFLEX_KINEMATICS_H
#define REFLEX_KINEMATICS_H

#ifdef __cplusplus
extern "C" {
#endif //__cplusplus

/* ---- MATRICES ---- */

/** Combine array of relative transforms into array of absolute transforms */
void rfx_kin_tf_chain( size_t n, const double T0[12], const double *TT_rel, double * TT_abs );

/** Compute columns of jacobian matrix for revolute links */
void rfx_kin_tf_jac_rev( size_t n, const double *TT_abs, const double *axis, const double pe[3],
                         double *J, size_t ldJ );


/** Compute one column of the Jacobian matrix for a revolute joint */
void rfx_kin_jac_col_rev( const double T_abs[12], const double axis_rel[3], const double pe_abs[3],
                          double *J );

/** Compute FK and Jacobian for a chain of links with revolute joints */
void rfx_kin_revchain( size_t n, const double T0[12], const double *TT_rel, const double Te_rel[12],
                       const double *axis,
                       double * TT_abs, double *J, size_t ldJ );


/* ---- Dual Quaternions ---- */

/** Combine array of relative transforms into array of absolute transforms */
void rfx_kin_duqu_chain( size_t n, const double T0[8], const double *TT_rel, double * TT_abs );

/** Compute columns of jacobian matrix for revolute links */
void rfx_kin_duqu_jac_rev( size_t n, const double *TT_abs, const double *axis, const double pe[3],
                           double *J, size_t ldJ );


/** Compute one column of the Jacobian matrix for a revolute joint */
void rfx_kin_duqu_jac_col_rev( const double T_abs[8], const double axis_rel[3], const double pe_abs[3],
                          double *J );

/** Compute FK and Jacobian for a chain of links with revolute joints */
void rfx_kin_duqu_revchain( size_t n, const double T0[8], const double *TT_rel, const double Te_rel[8],
                            const double *axis,
                            double * TT_abs, double *J, size_t ldJ );


/* ---- Kinematic Solvers ---- */

struct rfx_kin_solve_opts {
    double dt0;          ///< initial timestep

    double theta_tol;    ///< angle error tolerate
    double x_tol;        ///< translation error tolerance
    double dq_tol;       ///< translation error tolerance
    double s2min_dls;    ///< minimum square singular value for damped least squares

    double dx_dt;        ///< scaling for cartesian error
    double dq_dt;        ///< scaling for configuration

    double *q_ref;
};

/** Gradient descent kinematic solver
 */
int rfx_kin_solve( size_t n, const double *q0, const double S1[8],
                   rfx_kin_duqu_fun kin_fun,
                   double *q1,
                   struct rfx_kin_solve_opts *opts );

/** Position error to velocity.
 *
  *   Returns cartesian velocity to correct position error in unit
  *   time
 */
void rfx_kin_duqu_werr( const double S[8], const double S_ref[8],
                        double werr[6] );


/** Position error to scalars
 */
void rfx_kin_duqu_serr( const double S[8], const double S_ref[8],
                       double *theta, double *x );

#ifdef __cplusplus
}
#endif //__cplusplus

#endif //REFLEX_KINEMATICS_H
