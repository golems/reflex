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

#ifndef REFLEX_H
#define REFLEX_H
/** \file reflex.h */


/** \mainpage Reflex: Real-Time Control
 *
 * \section feature Features
 *
 * - Robot Arms
 *   - Workspace (Cartesian) Control
 *     - Damped-Least squares Jacobian Inverse with Nullspace Projection
 *   - Trapezoidal velocity profiles in workspace
 * - Linear Quadratic Gaussian \ref rfx_lqg_t
 *   - Solutions to Algebraic Riccati Equation
 *   - Discrete Time Kalman Filter,
 *     - \ref rfx_lqg_kf_predict
 *     - \ref rfx_lqg_kf_correct
 *   - Kalman-Bucy (continuous time) Filter
 *     - \ref rfx_lqg_kbf_gain
 *     - \ref rfx_lqg_kbf_step1, \ref rfx_lqg_kbf_step4
 *
 * \section wsctrl Workspace Control HOWTO
 *
 * -# Initialization
 *    -# get a \ref rfx_ctrl_t struct
 *    -# call \ref rfx_ctrl_ws_init
 *    -# get a \ref rfx_ctrl_ws_lin_k_t struct
 *    -# call \ref rfx_ctrl_ws_lin_k_init
 *    -# Set gains the the rfx_ctrl_ws_lin_k_t
 *    -# Set limits the the rfx_ctrl_t
 * -# At each time step (probably do this at 1kHz)
 *    -# Update the position q, velocity dq, and Jacobian J in the rfx_ctrl_t
 *    -# Update the reference q_r, x_r, etc in the rfx_ctrl_t
 *    -# Compute desired velocities with \ref rfx_ctrl_ws_lin_vfwd
 *    -# Send the velocities to your arm
 * -# Finalization
 *    -# Call \ref rfx_ctrl_ws_destroy
 *    -# Call \ref rfx_ctrl_ws_lin_k_destroy
 *
 * \author Neil T. Dantam
 * \author Developed at the Georgia Tech Humanoid Robotics Lab
 * \author Under Direction of Professor Mike Stilman
 *
 * \section License
 *
 * Copyright (c) 2010-2011, Georgia Tech Research Corporation.
 * All rights reserved.
 *
 *   Redistribution and use in source and binary forms, with or
 *   without modification, are permitted provided that the following
 *   conditions are met:
 *
 *   - Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *
 *   - Redistributions in binary form must reproduce the above
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


// Symbol:
/*
 * q: configuration
 * x: state/workspace position
 * u: input/command
 * z: sensor
 * dq: joint velocity
 * dq: workspace velocity
 * F: generalized force
 * R: rotation matrix
 * p: quaternion
 * k: gain
*/

#ifdef __cplusplus
extern "C" {
#endif //__cplusplus



typedef struct aa_tf_qv rfx_tf;
typedef struct aa_tf_qv_dx rfx_tf_dx;


// FIXME: controller should return these codes
typedef enum {
    RFX_OK = 0,
    RFX_INVAL,
    RFX_LIMIT_POSITION,
    RFX_LIMIT_POSITION_ERROR,
    RFX_LIMIT_FORCE,
    RFX_LIMIT_MOMENT,
    RFX_LIMIT_FORCE_ERROR,
    RFX_LIMIT_MOMENT_ERROR,
    RFX_LIMIT_CONFIGURATION,
    RFX_LIMIT_CONFIGURATION_ERROR
} rfx_status_t;

AA_API const char* rfx_status_string(rfx_status_t i);

#define RFX_PERROR(s, r) {                                              \
        rfx_status_t _rfx_$_perror_r = r;                               \
        const char *_rfx_$_perror_s = s;                                \
        fprintf(stderr, "%s: %s (%d)\n",                                \
                _rfx_$_perror_s ? _rfx_$_perror_s : "",                 \
                rfx_status_string(_rfx_$_perror_r), _rfx_$_perror_r);   \
    }

/************************/
/* Workspace Controller */
/************************/

typedef int (*rfx_kin_fun) ( const void *cx, const double *q, double *x, double *r, double *J);

typedef int (*rfx_kin_duqu_fun) ( const void *cx, const double *q, double S[8],  double *J);

typedef int (*rfx_ctrlx_fun) ( const void *cx, const double *q, const double *dq,
                               double *x, double *r, double *dx );

struct rfx_ctrlx_state {
    double *q;  ///< joint configuration
    double *dq; ///< joint velocity
    double *S;  ///< pose dual quaternion
    double *dx; ///< workspace velocity
    double *F;  ///< workspace forces
};

/** Workspace control state and reference values.
 *
 */
typedef struct {
    size_t n_q;      ///< size of config space
    struct rfx_ctrlx_state act; ///< actual state
    struct rfx_ctrlx_state ref; ///< reference state
    double *J;       ///< jacobian
    // limits
    double F_max;    ///< maximum linear force magnitude (<=0 to ignore)
    double M_max;    ///< maximum moment magnitude (<=0 to ignore)
    double e_q_max;  ///< maximum joint error (<=0 to ignore)
    double e_x_max;  ///< maximum workspace error (<=0 to ignore)
    double e_F_max;  ///< maximum linear force error (<=0 to ignore)
    double e_M_max;  ///< maximum moment error (<=0 to ignore)
    double *q_min;   ///< minimum joint values  (always checked)
    double *q_max;   ///< maximum joint values (always checked)
    double x_min[3]; ///< minimum workspace position (always checked)
    double x_max[3]; ///< maximum workspace position (always checked)
} rfx_ctrl_t;

typedef rfx_ctrl_t rfx_ctrl_ws_t;

/** Initialize workspace control state structure.
 *  Malloc's arrays for each field
 */
AA_API void rfx_ctrlx_state_init( struct rfx_ctrlx_state *x, size_t n );

/// initialize workspace controller
AA_API void rfx_ctrl_ws_init( rfx_ctrl_ws_t *g, size_t n );
/// destroy workspace controller
AA_API void rfx_ctrl_ws_destroy( rfx_ctrl_ws_t *g );

/** Gains for linear workspace control.
 */
typedef struct {
    size_t n_q;
    double p[6];   ///< position error gains
    double *q;     ///< configuration error gains
    //double v[6]; ///< velocity error gains
    double f[6];   ///< force error gains
    double dls;    ///< damped least squares k
    double s2min;  ///< deadzone damped least squares minimum square singular value
} rfx_ctrl_ws_lin_k_t;


/** Gains for linear jointspace control.
 */
typedef struct {
    size_t n_q;
    double *p;      ///< position error gains
} rfx_ctrlq_lin_k_t;

/// initialize
AA_API void rfx_ctrl_ws_lin_k_init( rfx_ctrl_ws_lin_k_t *k, size_t n_q );
/// destroy
AA_API void rfx_ctrl_ws_lin_k_destroy( rfx_ctrl_ws_lin_k_t *k );

/** Linear Workspace Control.
 *
 * \f[ u = J^*  (  \dot{x}_r - k_p(x - x_r) -  k_f(F - F_r) ) \f]
 *
 * Uses the singularity robust damped-least squares Jacobian inverse to
 * control an arm in workspace.
 *
 * The position error for orientation, \f$ \omega \in \Re^3 \f$, is calculated
 * from the axis-angle form of the relative orientation:
 *
 * \f[ \omega = \vec{a}\theta \f]
 *
 * \param ws The state and reference values
 * \param k The gains
 * \param u The configuration velocity to command, \f$ u \in \Re^{n_q} \f$
 */
AA_API rfx_status_t rfx_ctrl_ws_lin_vfwd( const rfx_ctrl_ws_t *ws, const rfx_ctrl_ws_lin_k_t *k,  double *u );


/** Linear jointspace control
 */
AA_API rfx_status_t rfx_ctrlq_lin_vfwd( const rfx_ctrl_t *g, const rfx_ctrlq_lin_k_t *k,  double *u );


/** Integrate the reference velocity in ws to produce an updated
 * reference position.
 *
 * Uses RK1/Euler integration.
 */
rfx_status_t rfx_ctrl_ws_sdx( rfx_ctrl_ws_t *ws, double dt );



typedef struct rfx_ctrlx_lin {
    rfx_ctrl_t *ctrl;
    rfx_ctrl_ws_lin_k_t *k;
    rfx_kin_fun kin_fun;
    void  *kin_fun_cx;
} rfx_ctrlx_lin_t;

rfx_ctrlx_lin_t *rfx_ctrlx_alloc( aa_mem_region_t *reg, size_t n_q, rfx_kin_fun kin_fun, void *kin_fun_cx );

AA_API rfx_status_t rfx_ctrlx_lin_vfwd( const rfx_ctrlx_lin_t *ctrl, const double *q,
                                        double *u );


/**************************/
/* A Simple PD Controller */
/* /\**************************\/ */

/* typedef struct { */
/*         size_t n_x;   ///< dimensionality */
/*         double *x;    ///< position */
/*         double *dx;   ///< velocity */
/*         double *x_r;  ///< reference position */
/*         double *dx_r; ///< reference velocity */
/*         double *k_p;  ///< position gain */
/*         double *k_d;  ///< velocity gain */
/* } ctrl_pd_t; */

/* /\* */
/*  * u = k * (x-x_r) + k*(dx-dx_r) */
/*  *\/ */
/* void ctrl_pd( const ctrl_pd_t *A, size_t n_u, double *u ); */

#include "reflex/lqg.h"
#include "reflex/kinematics.h"
#include "reflex/trajq.h"
#include "reflex/trajx.h"
#include "reflex/body.h"
#include "reflex/tf.h"

#ifdef __cplusplus
}
#endif //__cplusplus

#endif //REFLEX_H
