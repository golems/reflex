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
    aa_sys_affine_t sys = {.n = 2, .A = A, .D = (double[]){0,0}};
    double x1[2] = {0,1}, x0[2];
    double dt=.1;
    for( double t = 0; t < 10; t+=dt ) {
        memcpy(x0, x1, sizeof(x1));
        aa_rk4_step( 2, (aa_sys_fun*)aa_sys_affine, &sys,
                     t, dt,
                     x0, x1 );
        fprintf(f_x0, "%f %f\n", t, x1[0]);
        fprintf(f_x1, "%f %f\n", t, x1[1]);
    }
    fclose(f_x0);
    fclose(f_x1);
}


void kf() {
    rfx_lqg_t lqg;
    memset(&lqg,0,sizeof(lqg));
    lqg.n_x = 2;
    lqg.n_u = 1;
    lqg.n_z = 1;

    lqg.x = (double*)alloca(sizeof(double)*2);
    lqg.u = (double*)alloca(sizeof(double)*1);
    lqg.z = (double*)alloca(sizeof(double)*1);

    lqg.A = (double[]){0,-1,1,0};
    lqg.B = (double[]){0,1};
    lqg.C = (double[]){1,0};

    lqg.P = (double[]){1,0,0,1};
    lqg.V = (double[]){1,0,0,1};
    lqg.W = (double[]){1};

    lqg.Q = NULL;
    lqg.R = NULL;

    lqg.K = (double*)alloca(sizeof(double)*2);
    lqg.L = NULL;

    // test predict
    memcpy(lqg.x, (double[]){1,0}, sizeof(double)*lqg.n_x);
    memcpy(lqg.u, (double[]){2}, sizeof(double)*lqg.n_u);
    memcpy(lqg.z, (double[]){1.1}, sizeof(double)*lqg.n_z);
    rfx_lqg_kf_predict(&lqg);
    assert( aa_veq( 2, lqg.x, (double[]){0,1}, .001 ) );
    assert( aa_veq( 4, lqg.P, (double[]){2,0,0,2}, .001 ) );

    // test correct
    rfx_lqg_kf_correct(&lqg);
    assert( aa_veq( 2, lqg.K, (double[]){.666667,0}, .001 ) );
    assert( aa_veq( 2, lqg.x, (double[]){.733333,1}, .001 ) );
    assert( aa_veq( 4, lqg.P, (double[]){.666667,0,0,2}, .001 ) );
}

int main( int argc, char **argv ) {
    (void)argc;
    (void)argv;


    {
    double A[] = {1,2,3,4};
    double B[] = {1,2};
    double C[] = {2,1,2,1};
    double u[] = {5};
    double x[] = {4,2};
    double z[] = {8,2};
    double K[] = {8,2,1,8};
    double dx[2], zwork[2];

    rfx_lqg_observe(2,1,2, A,B,C, x,u,z, K, dx, zwork);

    assert( aa_veq( 2, dx, AA_FAR(-21,-14), .00001 ) );
    }

    spring();
    kf();
}





