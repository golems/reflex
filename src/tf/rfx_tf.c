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

int rfx_tf_madqg_predict
( double dt, double *E, double *dE, double *P, const double *V )
{
    return rfx_lqg_qutr_predict( dt, E, dE, P, V );
}

int rfx_tf_madqg_correct2
( double dt,
  size_t max_delta, double *delta_theta, double *delta_x, size_t *n_delta, size_t *i_delta,
  double *E_est, double *dx_est,
  size_t n_obs, const double *bEo, const double *cEo,
  double *P, const double *W )
{
    // check params
    if( *n_delta > max_delta || *i_delta > max_delta || *i_delta > *n_delta ) return -1;

    // TODO: Use the centrally-located orientation to transform each
    //       measurement as in rfx_tf_cor
    double bEc[7*n_obs];
    for( size_t i = 0; i < n_obs; i ++ ) {
        size_t j = 7*i;
        aa_tf_qutr_mulc( bEo+j, cEo+j, bEc+j );
    }

    return rfx_tf_madqg_correct( dt,
                                 max_delta, delta_theta, delta_x, n_delta, i_delta,
                                 E_est, dx_est,
                                 n_obs, bEc,
                                 P, W );
}

int rfx_tf_madqg_correct
( double dt,
  size_t max_delta, double *delta_theta, double *delta_x, size_t *n_delta, size_t *i_delta,
  double *E_est, double *dx_est,
  size_t n_obs, const double *E_obs,
  double *P, const double *W )
{
    //aa_dump_mat( stdout, W, 7, 7 );
    /* Find rels */
    const double *q_est = E_est;
    const double *x_est = q_est + 4;

    double q_use[4*n_obs];
    size_t n_q_use = 0;

    double x_use[3*n_obs];
    size_t n_x_use = 0;


    double dq_median;
    double dx_median;
    /* Maybe reject */
    if( *n_delta ) {
        dq_median = aa_la_d_median( *n_delta, delta_theta, 1 );
        dx_median = aa_la_d_median( *n_delta, delta_x, 1 );
    } else {
        dq_median = 1e9;
        dx_median = 1e9;
    }

    {
        // get medians
        // Check observatiosn
        for( size_t i = 0; i < n_obs; i ++ ) {
            const double *q_obs = E_obs + 7*i;
            const double *x_obs = q_obs + 4;
            // check angle
            double dq = aa_tf_qangle_rel( q_est, q_obs );
            if ( dq <= dq_median ) {
                AA_MEM_CPY( q_use + 4*n_q_use, q_obs, 4 );
                n_q_use++;
            }
            // check translation
            double dx = sqrt( aa_la_ssd(3, x_est, x_obs) );
            if( dx <= dx_median ) {
                AA_MEM_CPY( x_use + 3*n_x_use, x_obs, 3 );
                n_x_use++;
            }
            /* Insert into window */
            delta_theta[*i_delta] = dq;
            delta_x[*i_delta] = dx;
            *i_delta = (*i_delta + 1) % max_delta;
            assert( *i_delta < max_delta );
        }
    }

    *n_delta += n_obs;
    if( *n_delta > max_delta ) *n_delta = max_delta;

    // Average non-rejected observations
    /* TODO: be smarter about partial updates */
    if( n_q_use || n_x_use ) {
        double Z[7];

        double *q_obs = Z;
        if( n_q_use ) {
            aa_tf_quat_davenport( n_q_use, NULL, q_use, 4, q_obs );
            aa_tf_qminimize(q_obs);
        } else {
            AA_MEM_ZERO( dx_est+3, 3 );
            AA_MEM_CPY( q_obs, q_est, 4 );
        }

        double *x_obs = Z+4;
        if( n_x_use ) {
            aa_la_d_colmean( 3, n_x_use, x_use, 3, x_obs );
        } else {
            AA_MEM_CPY( x_obs, x_est, 3 );
            AA_MEM_ZERO( dx_est, 3 );
        }


        int r = rfx_lqg_qutr_correct( dt, E_est, dx_est, Z, P, W );
    } else {
        AA_MEM_ZERO( dx_est, 6 );
    }

    return 0;
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

void rfx_tf_rand( double theta_max, double x_max, double e[7] )
{
    // rotation
    double aa[4];
    aa_vrand( 4, aa );
    aa_la_normalize(3,aa);
    aa[3] *= theta_max;
    aa_tf_axang2quat(aa, e+AA_TF_QUTR_Q);

    // translation
    double v[3];
    aa_vrand( 3, v );
    for( size_t i = 0; i < 3; i ++ ) {
        e[4+i] = (v[i] - 0.5) * 2 * x_max;
    }
}

void rfx_tf_corrupt
( double theta_max, double x_max, const double e0[7], double e1[7] )
{
    double ec[7];
    rfx_tf_rand(theta_max,x_max, ec);

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

// Find TF from correspondences
void rfx_tf_cor( int opts, size_t n,
                 const double *qx, size_t ldqx,
                 const double *vx, size_t ldvx,
                 const double *qy, size_t ldqy,
                 const double *vy, size_t ldvy,
                 double *Z )
{

    double *top = AA_MEM_REGION_LOCAL_NEW_N( double, 1 );

    /*-- Orientation --*/
    double q_fit[3][4];
    size_t n_fit = 0;

    if( opts & RFX_TF_COR_O_ROT_UMEYAMA ) {
        double tf[12];
        rfx_tf_umeyama( n, vx, ldvx, vy, ldvy, tf );
        assert(aa_tf_isrotmat(tf));
        aa_la_transpose(3,tf);
        aa_tf_rotmat2quat( tf, q_fit[n_fit++] );
    }

    double *Qrel = AA_MEM_REGION_LOCAL_NEW_N( double, 4*n );
    for( size_t i = 0; i < n; i ++ ) {
        double *q = AA_MATCOL(Qrel,4,i);
        aa_tf_qmulc( AA_MATCOL(qx,ldqx,i),
                         AA_MATCOL(qy,ldqy,i),
                     q );
        //aa_tf_qminimize(q);
    }

    if( opts & RFX_TF_COR_O_ROT_DAVENPORT )
        aa_tf_quat_davenport( n, NULL, Qrel, 4, q_fit[n_fit++] );

    if( opts & RFX_TF_COR_O_ROT_MEDIAN )
        rfx_tf_qangmedian ( n, Qrel, 4, q_fit[n_fit++] );

    aa_tf_quat_davenport(n_fit, NULL, q_fit[0], 4, Z);

    /*-- deviation -- */
    /* double *angles = AA_MEM_REGION_LOCAL_NEW_N( double, n ); */
    /* for( size_t i = 0; i < n; i ++ ) { */
    /*     double *q = AA_MATCOL(Qrel,4,i); */
    /*     angles[i] = aa_tf_qangle_rel(Z, q); */
    /* } */
    /* double astd = aa_la_d_vecstd( n, angles, 1, 0 ); */
    /* double amad = aa_la_d_median( n, angles, 1 ); */
    /* double amax = aa_la_max( n, angles ); */

    /* fprintf(global_output, */
    /*         //"# %s\n" */
    /*         "# angle std: %f\n" */
    /*         "# angle mad: %f\n" */
    /*         "# angle max: %f\n", */
    /*         //comment, */
    /*         astd, amad, amax */
    /*     ); */

    /*-- Translation --*/
    double R[9];
    aa_tf_quat2rotmat(Z,R);
    assert( aa_tf_isrotmat(R) );

    // rels
    double *Zt = AA_MEM_REGION_LOCAL_NEW_N(double, 3*n);
    aa_cla_dlacpy( ' ', 3, (int)n,
                   vx, (int)ldvx,
                   Zt, 3 );

    for( size_t j = 0; j < n; j ++ ) {
        double yp[3] = {0};
        aa_tf_9rot( R, AA_MATCOL(vy,ldvy,j), yp );
        for( size_t i = 0; i < 3; i++ )
            AA_MATREF(Zt, 3, i, j) -=  yp[i];
        assert( aa_tf_isrotmat(R) );
    }

    assert( aa_tf_isrotmat(R) );
    if( opts & RFX_TF_COR_O_TRANS_MEAN ) {
        aa_la_d_colmean( 3, n, Zt, 3, Z+4 );
        //aa_tf_relx_mean( n, R, vy,ldvy, vx,ldvx, Z+4 );
    } else if (opts & RFX_TF_COR_O_TRANS_MEDIAN ) {
        //aa_tf_relx_median( n, R, vy,ldvy, vx,ldvx, Z+4 );
        for( size_t i = 0; i < 3; i ++ )
            (Z+4)[i] = aa_la_d_median( n, Zt+i, 3 );
    }


    aa_mem_region_local_pop(top);
}

int rfx_tf_filter_update( double dt, struct rfx_tf_filter *F)
{
    return rfx_tf_filter_update_work( dt,
                                      F->X.data, F->U.data, F->Z.data,
                                      F->P, F->V, F->W );
}
