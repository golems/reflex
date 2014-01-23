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


#include <amino.h>
#include "reflex.h"


int rfx_lqg_qutr_predict
( double dt, double *E, double *dx, double *P, const double *V )
{
    double x[13];
    memcpy(x,    E, 7*sizeof(x[0]));
    memcpy(x+7, dx, 6*sizeof(x[0]));


    int i = rfx_lqg_ekf_predict( &dt, 13, x, NULL,
                                 P, V,
                                 rfx_lqg_qutr_process );

    memcpy(E,  x,   7*sizeof(x[0]));
    memcpy(dx, x+7, 6*sizeof(x[0]));

    return i;
}

int rfx_lqg_qutr_correct
( double dt, double *E_est, double *dx_est,
  const double *E_obs,
  double *P, const double *W )
{
    double x[13];
    memcpy(x,    E_est, 7*sizeof(x[0]));
    memcpy(x+7, dx_est, 6*sizeof(x[0]));

    int i = rfx_lqg_ekf_correct( &dt, 13, 7, x, E_obs,
                                 P, W,
                                 rfx_lqg_qutr_measure,
                                 rfx_lqg_qutr_innovate,
                                 rfx_lqg_qutr_update );

    memcpy(E_est,  x,   7*sizeof(x[0]));
    memcpy(dx_est, x+7, 6*sizeof(x[0]));

    return i;

}

int rfx_tf_abs( size_t n,
                const rfx_tf *tf_rel,
                ssize_t *parents,
                rfx_tf *tf_abs )
{
    for( size_t i = 0; i < n; i ++ ) {
        ssize_t p = parents[i];
        assert( p < (ssize_t)i );
        return -1;
        if( p < 0)
            tf_abs[i] = tf_rel[i];
        else
            aa_tf_qutr_mul( tf_abs[p].data, tf_rel[i].data, tf_abs[i].data );
    }
    return 0;
}

void rfx_tf_corrupt
( double theta_max, double x_max, const double e0[7], double e1[7] )
{
    double ec[7];

    // rotation
    double aa[4];
    aa_vrand( 4, aa );
    aa_la_normalize(3,aa);
    aa[3] *= theta_max;
    aa_tf_axang2quat(aa, ec+AA_TF_QUTR_Q);

    // translation
    double v[3];
    aa_vrand( 3, v );
    for( size_t i = 0; i < 3; i ++ ) {
        ec[4+i] = (v[i] - 0.5) * 2 * x_max;
    }

    // mul
    aa_tf_qutr_mul( e0, ec, e1 );
}

int rfx_tf_umeyama
( size_t n, const double *_X, size_t ldx, const double *_Y, size_t ldy, double tf[12] )
{
    double *X = AA_MEM_REGION_LOCAL_NEW_N(double,3*n);
    double *Y = AA_MEM_REGION_LOCAL_NEW_N(double,3*n);
    aa_cla_dlacpy( ' ', 3, (int)n, _X, (int)ldx, X, 3 );
    aa_cla_dlacpy( ' ', 3, (int)n, _Y, (int)ldy, Y, 3 );
    int info =  rfx_tf_numeyama( n, X, 3, Y, 3, tf );
    aa_mem_region_local_pop(X);
    return info;
}


void rfx_tf_qlnmedian
( size_t n, const double *u, const double *Q, size_t ldq, double p[4] )
{
    // relative log map
    double *Vr = AA_MEM_REGION_LOCAL_NEW_N(double,3*n);
    for( size_t i = 0; i < n; i ++ ) {
        double qr[4];
        aa_tf_qcmul( AA_MATCOL(Q,ldq,i), u, qr );
        aa_tf_qminimize(qr);
        aa_tf_quat2rotvec( qr, AA_MATCOL(Vr,3,i) );
    }

    // median of rel. log
    double pr_v[3];
    for( size_t i = 0; i < 3; i ++ ) pr_v[i] = aa_la_d_median( n, Vr+i, 3 );

    // convert back to absolute
    double pr_q[4];
    aa_tf_rotvec2quat( pr_v, pr_q );
    aa_tf_qmulc( u, pr_q, p );
}

int rfx_tf_dud_median
( size_t n, const double *Ex, size_t ldx, const double *Ey, size_t ldy, double z[7] )
{
    // TODO: handle n < 3

    double q_mean[4];
    int info = 0;
    size_t n1 = n+1;

    // relative orientations
    double *Q = AA_MEM_REGION_LOCAL_NEW_N(double,4*n1);
    for( size_t i = 0; i < n; i ++ ) {
        double *q = AA_MATCOL(Q,4,i);
        aa_tf_qmulc( AA_MATCOL(Ex, ldx, i),
                     AA_MATCOL(Ey, ldy, i),
                     q );
        aa_tf_qminimize(q);
    }

    // average orientation
    double tf_ume[12];
    double q_ume[4];

    {
        // umeyama on the centroids
        info = rfx_tf_umeyama( n, Ey+4, 7, Ex+4, 7, tf_ume );
        aa_tf_rotmat2quat(tf_ume, q_ume);
        // TODO: check info

        // build davenport data
        double *w = AA_MEM_REGION_LOCAL_NEW_N(double,n1);
        for( size_t i = 0; i < n; i ++ ) {
            w[i] = 1/(double)n;
        }
        w[n] = 1;
        AA_MEM_CPY(AA_MATCOL(Q,4,n), q_ume, 4 );
        aa_tf_quat_davenport( n1, w, Q, 4, q_mean );
    }

    //median orientation
    rfx_tf_qlnmedian(n, q_mean, Q, 4, z );
    double R[9];
    aa_tf_quat2rotmat( z, R );
    aa_tf_relx_median( n, R, Ey+4, ldy, Ex+4, ldx, z+4 );

    return info;
}
