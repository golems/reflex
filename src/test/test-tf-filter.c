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
#include <cblas.h>
#include "reflex.h"


void x_fopen( const char *prefix, const char *suffix, size_t n, FILE *f[] ) {
    for( size_t i = 0; i < n; i ++ ) {
        char buf[512];
        sprintf(buf,"%s%lu%s", prefix, i, suffix);
        f[i] = fopen(buf, "w");
    }
}

void x_write( double t, size_t n, double e[], FILE *f[7] ) {
    for( size_t i = 0; i < n; i ++ ) {
        fprintf(f[i], "%f %f\n", t, e[i] );
        //fprintf(stdout, "%f %f\n", t, e[i] );
    }
}

void e_corrupt( double theta_max, double x_max, const double e0[7], double e1[7] )
{
    // rotation
    double aa[4];
    aa_vrand( 4, aa );
    aa_la_normalize(3,aa);
    aa[3] *= theta_max;

    // translation
    double v[3];
    aa_vrand( 3, v );
    for( size_t i = 0; i < 3; i ++ ) {
        v[i] *= x_max;
    }

    // convert
    double ec[7];
    aa_tf_axang2quat(aa, ec+AA_TF_QUTR_Q);
    AA_MEM_CPY(ec+AA_TF_QUTR_V, v, 3 );

    // mul
    aa_tf_qutr_mul( e0, ec, e1 );


}

rfx_tf_dx XX_true;
rfx_tf_dx XX_est;
rfx_tf_dx ZZ;
rfx_tf_dx UU;

double P[14*14] = {0};
double V[14*14] = {0};
double W[8*8] = {0};

double Pa[13*13] = {0};
double Wa[13*13] = {0};
double Va[13*13] = {0};


int main(void)
{

    FILE *fout_z[7];
    FILE *fout_x_true[7];
    FILE *fout_x_est[7];

    FILE *fout_dx_true[6];
    FILE *fout_dx_est[6];

    // state
    memset(&XX_true,0,sizeof(XX_true));
    memset(&XX_est,0,sizeof(XX_est));
    memset(&ZZ,0,sizeof(ZZ));
    memset(&UU,0,sizeof(UU));

    aa_la_diag( 14, P, 1.0 );
    aa_la_diag( 8, W, 1.0 );
    aa_la_diag( 14, V, 1.0e-2 );


    aa_la_diag( 13, Pa, 1.0e-2 );
    aa_la_diag( 13, Wa, 1.0 );
    aa_la_diag( 13, Va, 1.0e2 );

    // files

    x_fopen("xtrue", ".dat", 7,fout_x_true);
    x_fopen("z", ".dat", 7,fout_z);
    x_fopen("xest", ".dat", 7,fout_x_est);
    x_fopen("dx_true", ".dat", 6,fout_dx_true);
    x_fopen("dx_est", ".dat", 6,fout_dx_est);

    memcpy(XX_true.tf.r.data, aa_tf_quat_ident, 4*sizeof(XX_true.tf.r.data[0]));
    memcpy(XX_est.tf.r.data, aa_tf_quat_ident, 4*sizeof(XX_true.tf.r.data[0]));

    double dt = .01;
    for( double t = 0.0; t < 2.0; t += dt ) {
        fprintf(stderr, "t: %f, ", t );
        // pristine velocity
        double dx[6] = {cos(t*M_PI), 0, 0,
                        0, 0, sin(t*M_PI)};
        AA_MEM_CPY(XX_true.dx.data, dx, 6);

        // integrate true
        {
            double e_next[7];
            aa_tf_qutr_svel(XX_true.tf.data, XX_true.dx.data, dt, e_next );
            memcpy( XX_true.tf.data, e_next, sizeof(e_next) );
        }


        // corrupt measurement
        e_corrupt( 10*M_PI/180, 7e-2, XX_true.tf.data, ZZ.tf.data );

        // filter
        aa_tick("filter time: ");
        double S[8];
        double Sz[8];
        aa_tf_qutr2duqu( XX_est.tf.data, S);
        aa_tf_qutr2duqu( ZZ.tf.data, Sz);
        rfx_lqg_duqu_predict( dt, S, XX_est.dx.data, P, V );
        rfx_lqg_duqu_correct( 1,
                              S, XX_est.dx.data,
                              Sz, /*ZZ.dx.data, */
                              P, W );

        printf("foo\n");
        aa_tf_duqu2qutr( S, XX_est.tf.data);
        //rfx_tf_filter_update_work( dt, XX_est.tf.data, UU.tf.data, ZZ.tf.data, Pa, Va, Wa );
        aa_tock();
        printf("--\n");
        aa_dump_mat( stdout, P, 14, 14 );

        // print
        x_write( t, 7, XX_true.tf.data, fout_x_true );
        x_write( t, 6, XX_true.dx.data, fout_dx_true );
        x_write( t, 7, XX_est.tf.data, fout_x_est );
        x_write( t, 6, XX_est.dx.data, fout_dx_est );
        x_write( t, 7, ZZ.tf.data, fout_z );
    }


}
