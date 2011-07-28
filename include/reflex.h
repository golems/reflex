/* -*- mode: C; c-basic-offset: 4  -*- */
/* ex: set shiftwidth=4 expandtab: */
/*
 * Copyright (c) 2010, Georgia Tech Research Corporation
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *     * Redistributions of source code must retain the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer.
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials
 *       provided with the distribution.
 *     * Neither the name of the Georgia Tech Research Corporation nor
 *       the names of its contributors may be used to endorse or
 *       promote products derived from this software without specific
 *       prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY GEORGIA TECH RESEARCH CORPORATION ''AS
 * IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL GEORGIA
 * TECH RESEARCH CORPORATION BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

#ifndef REFLEX_H
#define REFLEX_H
/** \file reflex.h */

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

/** Workspace control state and reference values.
 *
 */
typedef struct {
    size_t n_q;      ///< size of config space
    // actual
    double *q;       ///< actual configuration
    double *dq;      ///< actual config velocity
    double *J;       ///< jacobian
    double x[3];     ///< actual workspace position
    double r[4];     ///< actual workspace orientation quaternion
    double F[6];     ///< actual workspace forces
    // reference
    double *q_r;     ///< reference confguration position
    double *dq_r;    ///< reference confguration velocity
    double x_r[3];   ///< reference workspace position
    double r_r[4];   ///< reference workspace orientation quaternion
    double dx_r[6];  ///< reference workspace velocity
    double F_r[6];   ///< reference workspace forces
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

AA_API void rfx_ctrl_ws_init( rfx_ctrl_ws_t *g, size_t n );
AA_API void rfx_ctrl_ws_destroy( rfx_ctrl_ws_t *g );

/** Gains for linear workspace control.
 */
typedef struct {
    size_t n_q;
    double p[6]; ///< position error gains
    double *q;   ///< configuration error gains
    //double v[6]; ///< velocity error gains
    double f[6]; ///< force error gains
    double dls;  ///< damped least squares k
} rfx_ctrl_ws_lin_k_t;

AA_API void rfx_ctrl_ws_lin_k_init( rfx_ctrl_ws_lin_k_t *k, size_t n_q );
AA_API void rfx_ctrl_ws_lin_k_destroy( rfx_ctrl_ws_lin_k_t *k );

/** Linear Workspace Control.
 *
 * \f[ u = J^*  (  \dot{x}_r - k_p(x - x_r) -  k_f(F - F_r) ) \f]
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

/****************************************/
/* Linear Quadratic Gaussian Controller */
/****************************************/

typedef struct {
    size_t n_x;
    size_t n_u;
    size_t n_z;
    double *A;   ///< process model
    double *B;   ///< input model
    double *C;   ///< measurement model
    double *V;   ///< process noise
    double *W;   ///< measurement noise
    double *Q;   ///< state cost
    double *R;   ///< input cost

    double *x;   ///< state estimate
    double *z;   ///< measurement
    double *u;   ///< computed input
} rfx_clqg_t;

AA_API void rfx_clqg_init( rfx_clqg_t *lqg, size_t n_x, size_t n_u, size_t n_z );
AA_API void rfx_clqg_destroy( rfx_clqg_t *lqg, size_t n_x, size_t n_u, size_t n_z );
AA_API void rfc_clqg_apply( rfx_clqg_t *lqg );


AA_API void rfx_lqg_observe_(
    const int *n_x, const int *n_u, const int *n_z,
    double *A, double *B, double *C,
    double *x, double *u, double *z,
    double *K, double *dx, double *zwork );

/** dx = Ax + Bu + K(z-Cx) */
AA_API void rfx_lqg_observe(
    size_t n_x, size_t n_u, size_t n_z,
    const double *A, const double *B, const double *C,
    const double *x, const double *u, const double *z,
    const double *K,
    double *dx, double *zwork );

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


/*************************/
/* LQG Filter/Controller */
/*************************/

typedef struct {
    size_t n_x;
    size_t n_z;
    size_t n_u;
} rfx_lqg_t;

#endif //REFLEX_H
