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
#include <inttypes.h>
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




void init_lqg(rfx_lqg_t *lqg) {
    rfx_lqg_init(lqg, 4, 1, 2 );
    // define some parameters
    const double M = .5;       // car mass
    const double m = 0.2;      // pendulum mass
    const double b = 0.1;      // cart friction
    const double g = 9.8;      // gravity
    const double l = 0.3;      // length to pendulum CoM
    const double i = 0.006;    // moment of interia of pendulum

    // denominator
    double p = i*(M+m)+M*m*l*l;

    // A: Process Model
    aa_la_transpose2( 4, 4,
                      (double[]){ 0, 1,               0,             0,
                                  0, -(i+m*l*l)*b/p, (m*m*g*l*l)/p,  0,
                                  0, 0,               0,             1,
                                  0, -(m*l*b)/p,      m*g*l*(M+m)/p, 0 },
                      lqg->A );

    // B: Input Model
    memcpy( lqg->B,
            (double[]){ 0, (i+m*l*l)/p, 0, m*l/p },
            sizeof(double)*4 );

    // C: Measurement Model
    aa_la_transpose2( 4, 2,
                      (double[]){ 1, 0, 0, 0,
                                  0, 0, 1, 0},
                      lqg->C );

    // Covariance, Noise, Cost as identity matrices
    aa_la_ident( lqg->n_x, lqg->P );
    aa_la_ident( lqg->n_x, lqg->V );
    aa_la_ident( lqg->n_x, lqg->Q );
    aa_la_ident( lqg->n_z, lqg->W );
    aa_la_ident( lqg->n_u, lqg->R );
}


void print_matrices( rfx_lqg_t *lqg ) {
    printf("A\n");
    printf("----\n");
    aa_dump_mat(stdout, lqg->A, lqg->n_x,lqg->n_x);
    printf("B\n");
    printf("----\n");
    aa_dump_mat(stdout, lqg->B, lqg->n_x,lqg->n_u);
    printf("C\n");
    printf("----\n");

    aa_dump_mat(stdout, lqg->C, lqg->n_z,lqg->n_x);
    printf("P\n");
    printf("----\n");
    aa_dump_mat(stdout, lqg->P, lqg->n_x,lqg->n_x);
    printf("V\n");
    printf("----\n");
    aa_dump_mat(stdout, lqg->V, lqg->n_x,lqg->n_x);
    printf("Q\n");
    printf("----\n");
    aa_dump_mat(stdout, lqg->Q, lqg->n_x,lqg->n_x);
    printf("W\n");
    printf("----\n");
    aa_dump_mat(stdout, lqg->W, lqg->n_z,lqg->n_z);
    printf("R\n");
    printf("----\n");
    aa_dump_mat(stdout, lqg->R, lqg->n_u,lqg->n_u);

    printf("L\n");
    printf("----\n");
    aa_dump_mat(stdout, lqg->L, lqg->n_u,lqg->n_x);


    printf("K\n");
    printf("----\n");
    aa_dump_mat(stdout, lqg->K, lqg->n_x,lqg->n_z);
}

int main( int argc, char **argv ) {
    (void)argc;
    (void)argv;

    srand((unsigned int)time(NULL)); // might break in 2038

    /* --- controller and observer --- */
    rfx_lqg_t lqg;
    init_lqg(&lqg);

    /* --- compute optimal gains --- */
    // Control LQR gain
    aa_tick("lqr: ",NULL);
    rfx_lqg_lqr_gain( &lqg );
    aa_tock();

    // Estimation Kalman-Bucy gain
    aa_tick("kbf: ",NULL);
    rfx_lqg_kbf_gain( &lqg );
    aa_tock();

    // check the the LQR gain is what we expect
    assert( aa_veq(4, lqg.L,
                   AA_FAR( -1.0000,  -2.0408,  20.3672,   3.9302),
                   .001) );


    /* open output files */
    FILE *f_x[4];
    FILE *f_xh[4];
    for( size_t i = 0; i < 4; i ++ ) {
        char buf[64];
        sprintf(buf, "cp_x%"PRIuPTR".dat", i);
        f_x[i] = fopen(buf, "w");
        sprintf(buf, "cp_xh%"PRIuPTR".dat", i);
        f_xh[i] = fopen(buf, "w");
    }
    FILE *f_u = fopen("cp_u.dat", "w");

    /* initial conditions */
    double x[4] = {0.1,0,0.0,0};

    // Step Size
    double dt=.010;

    /* --- Simulation and Control Loop --- */
    for( double t = 0; t < 10; t+=dt ) {
        // simulate
        {
            double x1[lqg.n_x];
            aa_lsim_rk4step( lqg.n_x, lqg.n_u,
                             dt, lqg.A, lqg.B,
                             x, lqg.u, x1 );
            memcpy(x, x1, sizeof(x1));
        }

        // measure
        cblas_dgemv( CblasColMajor, CblasNoTrans,
                     (int)lqg.n_z, (int)lqg.n_x,
                     1.0, lqg.C, (int)lqg.n_z,
                     x, 1,
                     0.0, lqg.z, 1 );

        // estimate
        rfx_lqg_kbf_step4(&lqg, dt);

        // control
        rfx_lqg_lqr_ctrl( &lqg );

        // output
        for( size_t i = 0; i < 4; i ++ ) {
            fprintf(f_x[i], "%f %f\n", t, x[i]);
            fprintf(f_xh[i], "%f %f\n", t, lqg.x[i]);
        }
        fprintf(f_u, "%f %f\n", t, lqg.u[0]);
    }

    /* close output files */
    for( size_t i = 0; i < 4; i ++ ) {
        fclose(f_x[i]);
        fclose(f_xh[i]);
    }
    fclose(f_u);
}





