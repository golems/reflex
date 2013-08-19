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
 *   - lqg.K contains space for the kalman gain matrix, size n_x*n_z
 * \post
 *   - lqg.K overwritten with the kalman gain
 */
AA_API void rfx_lqg_kbf_gain( rfx_lqg_t *lqg );

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


/****************/
/* Trajectories */
/****************/

/*--- Configuration Space Trajectories ---*/

struct rfx_trajq_vtab {
    int (*generate)(void *cx);
    int (*get_q)(void *cx, double t, double *q);
    int (*get_dq)(void *cx, double t, double *dq);
    int (*get_ddq)(void *cx, double t, double *ddq);
};

struct rfx_trajq_data {
    double t;
    double q[1];
};

typedef struct rfx_trajq {
    struct rfx_trajq_vtab *vtab;
    aa_mem_region_t *reg;
    size_t n_q;
    size_t n_t;
    double t_i, t_f;
    double *q_i, *q_f;
    aa_mem_rlist_t *point;
} rfx_trajq_t;

/* Initialize struct, performing future alloctions out of reg */
struct rfx_trajq *rfx_trajq_alloc( aa_mem_region_t *reg, size_t n_q );

void rfx_trajq_init( struct rfx_trajq *traj, aa_mem_region_t *reg, size_t n_q );

/* Add point q to trajectory.
 * q is copied to cx's internally managed memory */
void rfx_trajq_add( struct rfx_trajq *cx, double t, double *q );

void rfx_trajq_plot( struct rfx_trajq *cx, double dt );

static inline int rfx_trajq_get_q( struct rfx_trajq *cx, double t, double *q) {
    return cx->vtab->get_q( cx, t, q );
}

static inline int rfx_trajq_get_dq( struct rfx_trajq *cx, double t, double *dq) {
    return cx->vtab->get_dq( cx, t, dq );
}

static inline int rfx_trajq_get_ddq( struct rfx_trajq *cx, double t, double *ddq) {
    return cx->vtab->get_ddq( cx, t, ddq );
}

static inline int rfx_trajq_generate( struct rfx_trajq *cx ) {
    return cx->vtab->generate( cx );
}

typedef struct rfx_trajq_trapvel {
    struct rfx_trajq traj;
    double *dq_max;   ///< max velocity (specified)
    double *ddq_max;  ///< max accelerations (specified)
    double t_b;       ///< blend time
    double *dq_r;     ///< reference velocity (generated)
    double *ddq_r;    ///< reference acceleration (generated)
} rfx_trajq_trapvel_t;

void rfx_trajq_trapvel_init( struct rfx_trajq_trapvel *cx, aa_mem_region_t *reg, size_t n_q );

/*--- Cartesian Space Trajectories ---*/

/* TODO: Separate structure types for trajectory via point container
 * and segment list.  Then, can use a generic segment list struct for
 * any trajectory type.  Also gives some possibility to parallelize
 * generation and tracking, or even pass generated trajectories via
 * IPC (would need to fixup vtable links for that somehow).
 */

struct rfx_trajx;

struct rfx_trajx_vtab {
    int (*generate)(struct rfx_trajx *cx);
    void (*add)(struct rfx_trajx *cx, double t, double x[3], double r[4]);
    int (*get_x)(struct rfx_trajx *cx, double t, double x[3], double r[4]);
    int (*get_dx)(struct rfx_trajx *cx, double t, double dx[6]);
    int (*get_ddx)(struct rfx_trajx *cx, double t, double ddx[6]);
};

struct rfx_trajx_point {
    double t;
    double x[3];
    double r[4];
};

typedef struct rfx_trajx {
    struct rfx_trajx_vtab *vtab;
    struct aa_mem_region *reg;
    struct rfx_trajq *trajq;

    struct rfx_trajx_point *pt_i;
    struct rfx_trajx_point *pt_f;

    size_t n_p;

    /* Poses as vector and quaternion */
    aa_mem_rlist_t *point;
} rfx_trajx_t;

static inline void rfx_trajx_add( struct rfx_trajx *cx, double t, double x[3], double r[4]) {
    cx->vtab->add( cx, t, x, r );
}


static inline void rfx_trajx_add_duqu( struct rfx_trajx *cx, double t, const double S[8] ) {
    double r[4], x[3];
    aa_tf_duqu2qv(S, r, x );
    cx->vtab->add( cx, t, x, r );
}

static inline int rfx_trajx_get_x( struct rfx_trajx *cx, double t, double x[3], double r[4]) {
    return cx->vtab->get_x( cx, t, x, r );
}

static inline int rfx_trajx_get_x_duqu( struct rfx_trajx *cx, double t, double S[8] ) {
    double x[3], r[4];
    int i = cx->vtab->get_x( cx, t, x, r );
    aa_tf_qv2duqu( r, x, S );
    return i;
}

static inline int rfx_trajx_get_dx( struct rfx_trajx *cx, double t, double dx[6]) {
    return cx->vtab->get_dx( cx, t, dx );
}

static inline int rfx_trajx_get_ddx( struct rfx_trajx *cx, double t, double ddx[6]) {
    return cx->vtab->get_ddx( cx, t, ddx );
}


int rfx_trajx_set_ctrl( struct rfx_trajx *cx, double t, rfx_ctrlx_lin_t *ctrlx );

static inline int rfx_trajx_generate( struct rfx_trajx *cx ) {
    return cx->vtab->generate( cx );
}


void rfx_trajx_destroy( struct rfx_trajx *cx );

/*-- Rotation Vector orientations --*/
/* Initialize cartesian trajectory generator
 *
 * @pre: trajq has been initialized with length of 6
 */
void rfx_trajx_rv_init( struct rfx_trajx *cx, aa_mem_region_t *reg );



void rfx_trajx_slerp_init( struct rfx_trajx *cx, aa_mem_region_t *reg );


/*--- Cartesian Segments ---*/

struct rfx_trajx_seg {
    struct rfx_trajx_vtab *vtab;
    double t_i, t_f;
};

typedef struct rfx_trajx_seg_lerp {
    struct rfx_trajx_seg seg;
    double x_i[3], x_f[3];
    double r_i[4], r_f[4];
    double tau_i, tau_f;
    double dt;
} rfx_trajx_seg_lerp_t;

/* t_i and t_f are the times when we actually start (used to pick the segment to execute)
 * tau_f and tau_f are the times used to compute the interpolation
 */

struct rfx_trajx_seg *
rfx_trajx_seg_lerp_slerp_alloc( aa_mem_region_t *reg, double t_i, double t_f,
                                double tau_i, double x_i[3], double r_i[4],
                                double tau_f, double x_f[3], double r_f[4] );

typedef struct rfx_trajx_via {
    struct rfx_trajx trajx;
    struct rfx_trajx_seg **seg;
    size_t n_seg;
} rfx_trajx_via_t;


void rfx_trajx_via_init( struct rfx_trajx_via *cx, aa_mem_region_t *reg );

typedef struct rfx_trajx_plot_opts {
    int to_file;                ///< write plot to file
    rfx_ctrlx_lin_t * ctrlx;    ///< plot joint values with this controller
    double *q_0;                ///< initial joint position
} rfx_trajx_plot_opts_t;

void rfx_trajx_plot( struct rfx_trajx *cx, double dt, const struct rfx_trajx_plot_opts *opts );


/*--- Cartesian Parabolic Blends ---*/

typedef struct rfx_trajx_seg_lerp_rv {
    struct rfx_trajx_seg seg;
    double x_i[6], x_f[6];
    double tau_i, tau_f;
    double dt;
} rfx_trajx_seg_lerp_rv_t;

struct rfx_trajx_seg *
rfx_trajx_seg_lerp_rv_alloc( aa_mem_region_t *reg, double t_i, double t_f,
                             double tau_i, double x_i[6],
                             double tau_f, double x_f[6] );

typedef struct rfx_trajx_seg_blend_rv {
    struct rfx_trajx_seg seg;
    double x_i[6], dx_i[6], ddx[6];
    double tau_i;
} rfx_trajx_seg_blend_rv_t;
struct rfx_trajx_seg *
rfx_trajx_seg_blend_rv_alloc( aa_mem_region_t *reg, double t_i, double t_f,
                              double tau_i, double x_i[6], double dx_i[6], double ddx[6] );


typedef struct rfx_trajx_parablend {
    struct rfx_trajx_via via;
    double t_b;
} rfx_trajx_parablend_t;


void rfx_trajx_parablend_init( struct rfx_trajx_parablend *cx, aa_mem_region_t *reg, double t_b );

/*--- Cartesian Sphereical Parabolic Blends ---*/

typedef struct rfx_trajx_seg_blend_q {
    struct rfx_trajx_seg seg;
    double x0[3];
    double dx0[3];
    double ddx[3];

    double r_i[4];
    double r_j[4];
    double r_k[4];

    double t_b;
    double t_i;
    double t_j;
    double t_k;

    double t_ij;

    double tau_i;

    double ddu_ij;
    double ddu_jk;
    double ddu_j;
} rfx_trajx_seg_blend_q_t;

struct rfx_trajx_seg *
rfx_trajx_seg_blend_q_alloc( aa_mem_region_t *reg, double t_0, double t_1,
                             double t_b,
                             double t_i, double x_i[3], double r_i[4],
                             double t_j, double x_j[3], double r_j[4],
                             double t_k, double x_k[3], double r_k[4] );

typedef struct rfx_trajx_parablend rfx_trajx_splend_t;

void rfx_trajx_splend_init( rfx_trajx_splend_t *cx, aa_mem_region_t *reg, double t_b );

#include "reflex/kinematics.h"

#ifdef __cplusplus
}
#endif //__cplusplus

#endif //REFLEX_H
