/* -*- mode: C; c-basic-offset: 4 -*- */
/* ex: set shiftwidth=4 tabstop=4 expandtab: */
/*
 * Copyright (c) 2011, Georgia Tech Research Corporation
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


void spring()  {
    FILE *f_x0 = fopen("x0.dat","w");
    FILE *f_x1 = fopen("x1.dat","w");
    //FILE *f_nx0 = fopen("x0.dat","w");
    //FILE *f_nx1 = fopen("x1.dat","w");

    double A[4] = {0, -1, 1, 0};
    double z2[2] = {0,0};
    aa_sys_affine_t sys = {.n = 2, .A = A, .D = z2};
    double x1[2] = {0,1}, x0[2];
    double dt=.1;
    for( double t = 0; t < 10; t+=dt ) {
        memcpy(x0, x1, sizeof(x1));
        aa_odestep_rk4( 2, (aa_sys_fun*)aa_sys_affine, &sys,
                        t, dt,
                        x0, x1 );
        fprintf(f_x0, "%f %f\n", t, x1[0]);
        fprintf(f_x1, "%f %f\n", t, x1[1]);
    }
    fclose(f_x0);
    fclose(f_x1);
}


void init_lqg(rfx_lqg_t *lqg) {
    rfx_lqg_init(lqg, 2, 1, 1 );
    double A[4] = {0,-1, 1,-.5};
    AA_MEM_CPY( lqg->A, A, 4 );
    double B[2] = {0,1};
    AA_MEM_CPY( lqg->B, B, 2 );
    double C[2] = {1,0};
    AA_MEM_CPY( lqg->C, C, 2 );

    double P[4] = {1,0, 0,1};
    AA_MEM_CPY( lqg->P, P, 4 );
    double V[4] = {1,0, 0,1};
    AA_MEM_CPY( lqg->V, V, 4 );
    double W[1] = {1};
    AA_MEM_CPY( lqg->W, W, 1 );
}

void kf() {
    rfx_lqg_t kf_lqg;
    rfx_lqg_t kbf_lqg;
    rfx_lqg_t sys_lqg;
    init_lqg( &kf_lqg );
    init_lqg( &kbf_lqg );
    init_lqg( &sys_lqg );

    // spring sequence test
    AA_MEM_SET(kf_lqg.x, 0, kf_lqg.n_x);
    AA_MEM_SET(kf_lqg.u, 0, kf_lqg.n_u);
    AA_MEM_SET(kf_lqg.z, 0, kf_lqg.n_z);
    double P[4] = {1e0,0,0,1e0};
    AA_MEM_CPY(kf_lqg.P, P, kf_lqg.n_x*kf_lqg.n_x);
    FILE *f_tx0 = fopen("x0.dat","w");
    FILE *f_tx1 = fopen("x1.dat","w");
    FILE *f_z0 = fopen("z0.dat","w");
    FILE *f_kfx0 = fopen("kf_x0.dat","w");
    FILE *f_kfx1 = fopen("kf_x1.dat","w");
    FILE *f_kbfx0 = fopen("kbf_x0.dat","w");
    FILE *f_kbfx1 = fopen("kbf_x1.dat","w");


    sys_lqg.x[0]=.87;
    sys_lqg.x[1]=-0.125;

    double dt=.01;

    // euler approximation for the discrete kalman filter
    aa_la_scal(4, dt, kf_lqg.A);
    kf_lqg.A[0] += 1;
    kf_lqg.A[3] += 1;
    kf_lqg.W[0] = 100;

    // optimal kalman-bucy gain
    rfx_lqg_kbf_gain(&kbf_lqg);

    for( double t = 0; t < 20; t+=dt ) {
        // get some noise
        double zg[4];
        aa_stat_box_muller( aa_frand(), aa_frand(), &zg[0], &zg[1] );
        aa_stat_box_muller( aa_frand(), aa_frand(), &zg[2], &zg[3] );

        // integrate simulation
        double x1[sys_lqg.n_x];
        aa_odestep_rk4( sys_lqg.n_x, rfx_lqg_sys, &sys_lqg,
                        t, dt,
                        sys_lqg.x, x1 );

        // process noise
        sys_lqg.x[0] = x1[0] + aa_stat_z2x( zg[0], 0, .002 );
        sys_lqg.x[1] = x1[1] + aa_stat_z2x( zg[1], 0, .008 );

        // measurement noise
        kf_lqg.z[0] = sys_lqg.x[0] + aa_stat_z2x( zg[2], 0, .02 );
        kbf_lqg.z[0] = sys_lqg.x[0] + aa_stat_z2x( zg[2], 0, .02 );

        // discrete kf
        rfx_lqg_kf_predict(&kf_lqg);
        rfx_lqg_kf_correct(&kf_lqg);

        // kalman-bucy
        rfx_lqg_kbf_step4(&kbf_lqg, dt);

        // output
        fprintf(f_kfx0, "%f %f\n", t, kf_lqg.x[0]);
        fprintf(f_kfx1, "%f %f\n", t, kf_lqg.x[1]);
        fprintf(f_kbfx0, "%f %f\n", t, kbf_lqg.x[0]);
        fprintf(f_kbfx1, "%f %f\n", t, kbf_lqg.x[1]);
        fprintf(f_z0, "%f %f\n", t, kf_lqg.z[0]);
        fprintf(f_tx0, "%f %f\n", t, sys_lqg.x[0]);
        fprintf(f_tx1, "%f %f\n", t, sys_lqg.x[1]);
    }
    fclose(f_tx0);
    fclose(f_tx1);
    fclose(f_z0);
    fclose(f_kbfx0);
    fclose(f_kbfx1);
}


void kbf() {
    rfx_lqg_t lqg;
    init_lqg( &lqg );

    // test gain
    rfx_lqg_kbf_gain(&lqg);
    double K[2] = {1,0};
    assert( aa_veq( 2, lqg.K, K, .001 ) );

    /* memcpy(lqg.P, (double[]){1,0,0,1}, sizeof(double)*lqg.n_x*lqg.n_x); */
    /* rfx_lqg_kf_predict(&lqg); */
    /* assert( aa_veq( 4, lqg.P, (double[]){2,0,0,2}, .001 ) ); */

    /* // test correct */
    /* rfx_lqg_kf_correct(&lqg); */
    /* assert( aa_veq( 2, lqg.K, (double[]){.666667,0}, .001 ) ); */
    /* assert( aa_veq( 2, lqg.x, (double[]){.733333,1}, .001 ) ); */
    /* assert( aa_veq( 4, lqg.P, (double[]){.666667,0,0,2}, .001 ) ); */


    /* // spring sequence test */
    /* memcpy(lqg.x, (double[]){0,0}, sizeof(double)*lqg.n_x); */
    /* memcpy(lqg.u, (double[]){0}, sizeof(double)*lqg.n_u); */
    /* memcpy(lqg.z, (double[]){0}, sizeof(double)*lqg.n_z); */
    /* memcpy(lqg.P, (double[]){1e0,0,0,1e0}, sizeof(double)*lqg.n_x*lqg.n_x); */
    /* lqg.A = (double[]){0,-1,1,-.5}; */

    /* lqg.W = (double[]){10}; */
    /* lqg.V = (double[]){1,0,0,1}; */

    /* FILE *f_tx0 = fopen("kf_true_x0.dat","w"); */
    /* FILE *f_tx1 = fopen("kf_true_x1.dat","w"); */
    /* FILE *f_zx0 = fopen("kf_sense_x0.dat","w"); */
    /* FILE *f_fx0 = fopen("kf_filt_x0.dat","w"); */
    /* FILE *f_fx1 = fopen("kf_filt_x1.dat","w"); */
    /* //FILE *f_nx0 = fopen("x0.dat","w"); */
    /* //FILE *f_nx1 = fopen("x1.dat","w"); */



    /* aa_sys_affine_t sys = {.n = 2, .A = lqg.A, .D = (double[]){0,0}}; */
    /* double x1[2] = {.87,-.125}, x0[2]; */
    /* double dt=.025; */
    /* lqg.A = (double*)alloca(sizeof(double)*4); */
    /* memcpy(lqg.A, sys.A, 4*sizeof(double)); */
    /* aa_la_scal(4, dt, lqg.A); */
    /* lqg.A[0] += 1; */
    /* lqg.A[3] += 1; */
    /* for( double t = 0; t < 20; t+=dt ) { */
    /*     memcpy(x0, x1, sizeof(x1)); */
    /*     double zg[2]; */
    /*     aa_box_muller( aa_frand(), aa_frand(), &zg[0], &zg[1] ); */
    /*     aa_rk4_step( 2, (aa_sys_fun*)aa_sys_affine, &sys, */
    /*                  t, dt, */
    /*                  x0, x1 ); */
    /*     x1[0] += aa_z2x( zg[0], 0, .01 ); */
    /*     rfx_lqg_kf_predict(&lqg); */
    /*     lqg.z[0] = x1[0] + aa_z2x( zg[1], 0, .02 ); */
    /*     rfx_lqg_kf_correct(&lqg); */
    /*     fprintf(f_fx0, "%f %f\n", t, lqg.x[0]); */
    /*     fprintf(f_fx1, "%f %f\n", t, lqg.x[1]); */
    /*     fprintf(f_zx0, "%f %f\n", t, lqg.z[0]); */
    /*     fprintf(f_tx0, "%f %f\n", t, x1[0]); */
    /*     fprintf(f_tx1, "%f %f\n", t, x1[1]); */
    /* } */
    /* fclose(f_tx0); */
    /* fclose(f_tx1); */
    /* fclose(f_fx0); */
    /* fclose(f_zx0); */
    /* fclose(f_fx1); */

}

int main( int argc, char **argv ) {
    (void)argc;
    (void)argv;

    srand((unsigned int)time(NULL)); // might break in 2038

    //spring();
    kf();
    kbf();
}
