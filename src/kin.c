/* -*- mode: C; c-basic-offset: 4 -*- */
/* ex: set shiftwidth=4 tabstop=4 expandtab: */
/*
 * Copyright (c) 2010, Georgia Tech Research Corporation
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

#include <amino.h>
#include "reflex.h"

struct kin_solve_cx {
    size_t n;
    struct rfx_kin_solve_opts *opts;
    rfx_kin_duqu_fun kin_fun;
    const double *S1;
};


void rfx_kin_duqu_werr( const double S[8], const double S_ref[8], double werr[6] ) {
    double twist[8], de[8];
    aa_tf_duqu_mulc( S, S_ref, de );  // de = d*conj(d_r)
    aa_tf_duqu_minimize(de);
    aa_tf_duqu_ln( de, twist );     // twist = log( de )
    aa_tf_duqu_twist2vel( S, twist, werr );
}


void rfx_kin_duqu_serr( const double S[8], const double S_ref[8],
                       double *theta, double *x ) {
    double S_rel[8];
    aa_tf_duqu_cmul( S, S_ref, S_rel ); // relative dual quaternion
    aa_tf_duqu_minimize(S_rel);
    //printf("srel: "); aa_dump_vec(stdout, S_rel, 8 );

    // quaternion angle
    *theta = atan2( sqrt(S_rel[0]*S_rel[0] + S_rel[1]*S_rel[1] + S_rel[2]*S_rel[2]),
                    S_rel[3] );

    // translation
    double xe[3];
    aa_tf_duqu_trans( S_rel, xe );
    *x = sqrt( xe[0]*xe[0] + xe[1]*xe[1] + xe[2]*xe[2] );
}


static void kin_solve_sys( const void *vcx,
                           double t, const double *AA_RESTRICT q,
                           double *AA_RESTRICT dq ) {
    (void) t; // time invariant
    const struct kin_solve_cx *cx = (const struct kin_solve_cx*)vcx;

    // compute kinematics
    double S[8];
    double J[6*cx->n];
    double J_star[6*cx->n];
    cx->kin_fun( NULL, q, S, J );

    // position error
    double w_e[6];
    rfx_kin_duqu_werr( S, cx->S1, w_e );
    for( size_t i = 0; i < 6; i ++ ) w_e[i] *= -cx->opts->dx_dt;


    // damped least squares
    aa_la_dzdpinv( 6, cx->n, cx->opts->s2min_dls, J, J_star );
    if( cx->opts->q_ref ) {
        // nullspace projection
        double dqnull[cx->n];
        for( size_t i = 0; i < cx->n; i ++ )  {
            dqnull[i] =  -.003*( q[i] - cx->opts->q_ref[i] );
        }
        aa_la_xlsnp( 6, cx->n, J, J_star, w_e, dqnull, dq );
    } else {
        aa_la_mvmul(7,6,J_star,w_e,dq);
    }
}


/* Levenberg Marquedt */
int rfx_kin_solve( size_t n, const double *q0, const double S1[8],
                   rfx_kin_duqu_fun kin_fun,
                   double *q1,
                   struct rfx_kin_solve_opts *opts ) {

    struct kin_solve_cx cx;
    cx.n = n;
    cx.opts = opts;
    cx.S1 = S1;
    cx.kin_fun = kin_fun;

    AA_MEM_CPY(q1, q0, n);
    double dq_norm = 0;

    int iters = 0;

    do {
        iters++;
        // check error
        double S[8], theta_err, x_err;
        kin_fun( NULL, q1, S, NULL );
        rfx_kin_duqu_serr( S, S1, &theta_err, &x_err );

        //printf("err: theta: %f, x: %f, dqn: %f\n", theta_err, x_err, dq_norm);

        if( fabs(theta_err) < opts->theta_tol &&
            fabs(x_err) < opts->x_tol &&
            dq_norm < opts->dq_tol )
        {
            break;
        }

        dq_norm = 0;

        // integrate
        double dq[n];
        kin_solve_sys( &cx, 0, q1, dq );
        for( size_t i = 0; i < n; i ++ ) {
            dq_norm += dq[i] * dq[i];
            q1[i] += dq[i];
        }

        /* double q[7]; */
        /* aa_odestep_rk4( 7, kin_solve_sys, &cx, */
        /*                 0, opts->dt0, */
        /*                 q1, q ); */
        /* for( size_t i = 0; i < 7; i ++ ) { */
        /*     dq_norm += (q1[i]-q[i]) * (q1[i]-q[i]); */
        /*     q1[i] = q[i]; */
        /* } */

        //printf("q: "); aa_dump_vec(stdout, q1, 7 );
        //AA_MEM_CPY(q, q1, n);
    } while(1);

    //printf(" iter: %d\n", iters );
    //printf(" norm: %f\n", aa_la_dot( 7, q1, q1 ) );
    return 0;
}
