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

// huber
// cost =  norm( log( z^* x ) )


void rfx_tf_qangmedian
( size_t n, const double *Q, size_t ldq, double p[4] )
{
    if( 0 == n ) {
        AA_MEM_CPY(p, aa_tf_quat_ident, 4 );
        return;
    } else if (1 == n) {
        aa_tf_qslerp( .5, AA_MATCOL(Q,ldq,0), AA_MATCOL(Q,ldq,1), p );
        return;
    } // else do stuff


    // TODO: parallelize, because this is slow

    double *sum_dist = AA_MEM_REGION_LOCAL_NEW_N(double,n);
    AA_MEM_ZERO(sum_dist, n);

    for( size_t i = 0; i < n; i ++ ) {
        const double *qi = AA_MATCOL(Q,ldq,i);
        for( size_t j = 0; j < i; j ++ ) {
            const double *qj = AA_MATCOL(Q,ldq,j);
            double dist = aa_tf_qangle_rel(qi, qj);
            sum_dist[i] += dist;
            sum_dist[j] += dist;
        }
    }

    size_t i_min = aa_fminloc( n, sum_dist );
    AA_MEM_CPY( p, AA_MATCOL(Q, ldq, i_min), 4 );
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
    aa_tf_qminimize(p);
}

void rfx_tf_dud_qrel
( size_t n, const double *Qx, size_t ldx, const double *Qy, size_t ldy, double *Q, size_t ldq )
{
    // relative orientations
    for( size_t i = 0; i < n; i ++ ) {
        double *q = AA_MATCOL(Q,ldq,i);
        aa_tf_qmulc( AA_MATCOL(Qx, ldx, i),
                     AA_MATCOL(Qy, ldy, i),
                     q );
        aa_tf_qminimize(q);
    }
}

int rfx_tf_dud_qmean
( size_t n, const double *Ex, size_t ldx, const double *Ey, size_t ldy, double *Q, size_t ldq, double q_mean[4] )
{
    int info = 0;
    size_t n1 = n+1;

    // relatitivies
    rfx_tf_dud_qrel(n, Ex, ldx, Ey, ldy, Q, ldq );

    // size checks
    if( 1 == n ) {
        AA_MEM_CPY(q_mean, Q, 4);
        return 0;
    } else if( 2 == n ) {
        aa_tf_qslerp( .5, AA_MATCOL(Q,ldq,0), AA_MATCOL(Q,ldq,1), q_mean );
        return 0;
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

    aa_tf_qminimize(q_mean);
    return info;
}

int rfx_tf_dud_median
( size_t n, const double *Ex, size_t ldx, const double *Ey, size_t ldy, double z[7] )
{
    // TODO: handle n < 3
    int info = 0;
    size_t n1 = n+1;

    double *Q = AA_MEM_REGION_LOCAL_NEW_N(double,4*n1);


    //double q_mean[4];
    //rfx_tf_dud_qmean( n, Ex, ldx, Ey, ldy, Q, 4, q_mean );
    //rfx_tf_qlnmedian(n, q_mean, Q, 4, z );

    rfx_tf_dud_qrel(n, Ex, ldx, Ey, ldy, Q, 4 );
    aa_tick("med: ");
    rfx_tf_qangmedian( n, Q, 4, z );
    aa_tock();

    //median translation
    double R[9];
    aa_tf_quat2rotmat( z, R );
    aa_tf_relx_median( n, R, Ey+4, ldy, Ex+4, ldx, z+4 );

    aa_mem_region_local_pop(Q);
    return info;
}

size_t rfx_tf_dud_reject
( size_t n, double z_theta, double z_x, const double *e_mu,
  double *Ex, size_t ldx, double *Ey, size_t ldy )
{
    // compute distances
    double *d_ang = AA_MEM_REGION_LOCAL_NEW_N( double, n );
    double *d_x = AA_MEM_REGION_LOCAL_NEW_N( double, n );
    double R[9];
    aa_tf_quat2rotmat( e_mu, R );
    for( size_t i = 0; i < n; i ++ ) {
        // angle
        double *ex = AA_MATCOL(Ex,ldx,i);
        double *ey = AA_MATCOL(Ey,ldy,i);
        double q[4];
        aa_tf_qmulc( ex, ey,  q );
        aa_tf_qminimize(q);
        d_ang[i] = fabs(aa_tf_qangle_rel(e_mu, q));
        // trans
        // dx = || x - Ry ||
        double yp[3];
        aa_tf_9rot( R, ey+4, yp );
        for( size_t k = 0; k < 3; k++ ) yp[k] = ex[4+k] - yp[k];
        d_x[i] = sqrt( aa_la_ssd(3, e_mu+4, yp) );
    }
    double mad_ang = aa_la_d_median( n, d_ang, 1 );
    double mad_x = aa_la_d_median( n, d_x, 1 );

    size_t j = 0;
    for( size_t i = 0; i < n; i ++ ) {
        // TODO: div by zero
        if( d_ang[i] / mad_ang < z_theta &&
            d_x[i] / mad_x < z_x )
        {
            if( i != j ) {
                AA_MEM_CPY( AA_MATCOL(Ex,ldx,j), AA_MATCOL(Ex,ldx,i), 7 );
                AA_MEM_CPY( AA_MATCOL(Ey,ldx,j), AA_MATCOL(Ey,ldy,i), 7 );
            }
            j++;
        }
    }
    printf("rejected: %lu\n", n - j );
    return j;
}

int rfx_tf_dud_rejected_mean
( size_t n, double z_theta, double z_x, const double *Ex, size_t ldx, const double *Ey, size_t ldy, double e[7] )
{
    //TODO: handle n < 3
    double *top_ptr = AA_MEM_REGION_LOCAL_NEW_N( double, 1 );

    double e_median[7];
    rfx_tf_dud_median( n, Ex, ldx, Ey, ldy, e_median );

    // reject outliers
    double *Exp = AA_MEM_REGION_LOCAL_NEW_N( double, n*7 );
    double *Eyp = AA_MEM_REGION_LOCAL_NEW_N( double, n*7 );
    aa_cla_dlacpy( 0, 7, (int)n, Ex, (int)ldx, Exp, 7 );
    aa_cla_dlacpy( 0, 7, (int)n, Ey, (int)ldy, Eyp, 7 );
    size_t j = rfx_tf_dud_reject( n, z_theta, z_x, e_median,
                                  Exp, 7, Eyp, 7 );

    // TODO: check size
    // re-average
    double *Q = AA_MEM_REGION_LOCAL_NEW_N( double, j*4 );
    rfx_tf_dud_qmean( j, Exp, 7, Eyp, 7, Q, 4, e );
    double R[9];
    aa_tf_quat2rotmat( e, R );
    aa_tf_relx_mean( j, R, Eyp+4, 7, Exp+4, 7, e+4 );

    //printf("reject: %lu\n", n-j);

    aa_mem_region_local_pop(top_ptr);
    return 0;
 }
