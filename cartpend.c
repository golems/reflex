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



/*         .     /
 *         .    /
 *         . th/
 *         .  /
 *         . /
 *      ___./__
 *     |       |
 *     |   M   |    -> x
 *     |_______|
 *       o   o
 */

// X: {x,dx,phi,dphi}





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


void init_lqg(rfx_lqg_t *lqg) {
    rfx_lqg_init(lqg, 4, 1, 2 );
    const double M = .5;       // car mass
    const double m = 0.2;      // pendulum mass
    const double b = 0.1;      // cart friction
    const double g = 9.8;      // gravity
    const double l = 0.3;      // length to pendulum CoM
    const double i = 0.006;    // moment of interia of pendulum

    double p = i*(M+m)+M*m*l*l;
    aa_la_transpose2( 4, 4,
                      AA_FAR( 0, 1,               0,             0,
                              0, -(i+m*l*l)*b/p, (m*m*g*l*l)/p,  0,
                              0, 0,               0,             1,
                              0, -(m*l*b)/p,      m*g*l*(M+m)/p, 0 ),
                      lqg->A );

    memcpy( lqg->B,
            AA_FAR( 0, (i+m*l*l)/p, 0, m*l/p ),
            sizeof(double)*4 );

    aa_la_transpose2( 4, 2,
                      AA_FAR( 1, 0, 0, 0,
                              0, 0, 1, 0),
                      lqg->C );

    aa_la_ident( lqg->n_x, lqg->P );
    aa_la_ident( lqg->n_x, lqg->V );
    aa_la_ident( lqg->n_x, lqg->Q );
    aa_la_ident( lqg->n_z, lqg->W );
    aa_la_ident( lqg->n_u, lqg->R );
}

int main( int argc, char **argv ) {
    (void)argc;
    (void)argv;

    srand((unsigned int)time(NULL)); // might break in 2038

    rfx_lqg_t sys; // simulated system
    rfx_lqg_t lqg; // controller and observer

    init_lqg(&lqg);
    init_lqg(&sys);

    printf("A\n");
    printf("----\n");
    aa_dump_mat(stdout, lqg.A, lqg.n_x,lqg.n_x);
    printf("B\n");
    printf("----\n");
    aa_dump_mat(stdout, lqg.B, lqg.n_x,lqg.n_u);
    printf("C\n");
    printf("----\n");
    aa_dump_mat(stdout, lqg.C, lqg.n_z,lqg.n_x);
    printf("P\n");
    printf("----\n");
    aa_dump_mat(stdout, lqg.P, lqg.n_x,lqg.n_x);
    printf("V\n");
    printf("----\n");
    aa_dump_mat(stdout, lqg.V, lqg.n_x,lqg.n_x);
    printf("Q\n");
    printf("----\n");
    aa_dump_mat(stdout, lqg.Q, lqg.n_x,lqg.n_x);
    printf("W\n");
    printf("----\n");
    aa_dump_mat(stdout, lqg.W, lqg.n_z,lqg.n_z);
    printf("R\n");
    printf("----\n");
    aa_dump_mat(stdout, lqg.R, lqg.n_u,lqg.n_u);


    //aa_tick("lqr: ",NULL);
    rfx_lqg_lqr_gain( &lqg );
    //aa_tock();


    //aa_tick("kbf: ",NULL);
    rfx_lqg_kbf_gain( &lqg );
    //aa_tock();

    printf("L\n");
    printf("----\n");
    aa_dump_mat(stdout, lqg.L, lqg.n_u,lqg.n_x);


    printf("K\n");
    printf("----\n");
    aa_dump_mat(stdout, lqg.K, lqg.n_x,lqg.n_z);

    assert( aa_veq(4, lqg.L,
                   AA_FAR( -1.0000,  -2.0408,  20.3672,   3.9302),
                   .001) );


    FILE *f_x[4];
    for( size_t i = 0; i < 4; i ++ ) {
        char buf[64];
        sprintf(buf, "cp_x%d.dat", i);
        f_x[i] = fopen(buf, "w");
    }
    FILE *f_u = fopen("cp_u.dat", "w");

    double dt=.01;
    sys.x[2] = .1;


    for( double t = 0; t < 5; t+=dt ) {
        // simulate
        double x1[sys.n_x];
        aa_rk4_step( sys.n_x, rfx_lqg_sys, &sys,
                     t, dt,
                     sys.x, x1 );
        memcpy(sys.x, x1, sizeof(x1));

        // measure
        cblas_dgemv( CblasColMajor, CblasNoTrans,
                     (int)lqg.n_z, (int)lqg.n_x,
                     1.0, lqg.C, (int)lqg.n_z,
                     sys.x, 1,
                     0.0, lqg.z, 1 );
        // observe
        rfx_lqg_kbf_step4(&lqg, dt);

        // control
        rfx_lqg_lqr_ctrl( &lqg );
        memcpy( sys.u, lqg.u, sizeof(double)*lqg.n_u );

        // output
        for( size_t i = 0; i < 4; i ++ ) {
            fprintf(f_x[i], "%f %f\n", t, sys.x[i]);
        }
        fprintf(f_u, "%f %f\n", t, sys.u[0]);
    }

    for( size_t i = 0; i < 4; i ++ ) {
        fclose(f_x[i]);
    }
    fclose(f_u);
}





