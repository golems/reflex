/* -*- mode: C++; c-basic-offset: 4 -*- */
/* ex: set shiftwidth=4 tabstop=4 expandtab: */
/*
 * Copyright (c) -2013, Georgia Tech Research Corporation
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


void plotq () {
    aa_mem_region_t reg;
    aa_mem_region_init( &reg, 1024*64 );


    double q0[6] =  {0.296500, 0.000000, 0.178000, 0.000000, 3.141593, 0.000000};
    double q1[6] =  {0.380000, 0.003357, 0.089718, 0.482354, 3.045461, 0.989530};

    size_t n_q = sizeof(q0)/sizeof(q0[0]);
    rfx_trajq_trapvel_t traj;
    rfx_trajq_trapvel_init( &traj, &reg, n_q );

    for( size_t i = 0; i < n_q; i ++ ) {
        traj.dq_max[i] = 10.0;
        traj.ddq_max[i] = 10.0;
    }

    rfx_trajq_add( &traj.traj, 0, q0 );
    rfx_trajq_add( &traj.traj, 10, q1 );

    printf("q0: ");
    aa_dump_vec( stdout, traj.traj.q_i, traj.traj.n_q );
    printf("q1: ");
    aa_dump_vec( stdout, traj.traj.q_f, traj.traj.n_q );

    rfx_trajq_generate( &traj.traj );

    printf("times: %f, %f\n", traj.traj.t_i, traj.traj.t_f );
    printf("blend: %f\n", traj.t_b );
    printf("dq_r: " );
    aa_dump_vec( stdout, traj.dq_r, traj.traj.n_q );
    printf("ddq_r: " );
    aa_dump_vec( stdout, traj.ddq_r, traj.traj.n_q );

    fprintf( stderr, "%lu\n",
             traj.traj.n_q );

    printf("q0: ");
    aa_dump_vec( stdout, traj.traj.q_i, traj.traj.n_q );
    printf("q1: ");
    aa_dump_vec( stdout, traj.traj.q_f, traj.traj.n_q );

    rfx_trajq_plot( &traj.traj, .001 );

    aa_mem_region_destroy( &reg );
}



void plot_qseg () {
    aa_mem_region_t reg;
    aa_mem_region_init( &reg, 1024*64 );


    double q0[6] =  {0.296500, 0.000000, 0.178000, 0.000000, 3.141593, 0.000000};
    double q1[6] =  {0.380000, 0.003357, 0.089718, 0.482354, 3.045461, 0.989530};
    double q2[6] =  {1, 2, 3, 4, 5, 6};

    size_t n_q = sizeof(q0)/sizeof(q0[0]);

    rfx_trajq_points_t *points = rfx_trajq_points_alloc( &reg, n_q );
    rfx_trajq_points_add( points, 0, q0 );
    rfx_trajq_points_add( points, 10, q1 );
    rfx_trajq_points_add( points, 20, q2 );

    printf("n points: %lu\n", points->n_t );
    printf("points 0: %lx\n", points->point->head );
    printf("points 1: %lx\n", points->point->head->next );

    rfx_trajq_seg_list_t *seglist = rfx_trajq_gen_pblend_tm1( &reg, points, 1 );

    rfx_trajq_seg_plot( seglist, .001 );

    aa_mem_region_destroy( &reg );
}

void plotx2() {
    aa_mem_region_t reg;
    aa_mem_region_init( &reg, 1024*32 );

    aa_mem_region_release( &reg );

    /* rfx_trajq_trapvel_t trajq; */
    /* rfx_trajq_trapvel_init( &trajq, &reg, 6 ); */
    /* for( size_t i = 0; i < 6; i ++ ) { */
    /*     trajq.dq_max[i] = 10.0; */
    /*     trajq.ddq_max[i] = 10.0; */
    /* } */


    double x0[3] = {0.384311,0.004544,0.089198};
    double r0[4] = {0.148775, 0.939373, 0.305129, -0.048381};
    double x1[3] = {0.380000, 0.053318, 0.020952};
    double r1[4] = {0.148778, 0.939347, 0.305212, -0.048341};

    rfx_trajx_t trajx;
    rfx_trajx_rv_init( &trajx, &reg );

    rfx_trajx_add( &trajx, 0, x0, r0 );
    rfx_trajx_add( &trajx, 10, x1, r1 );

    rfx_trajx_generate( &trajx );

    rfx_trajq_plot( trajx.trajq, .001 );

    aa_mem_region_destroy( &reg );
}

void plot_viax() {
    aa_mem_region_t reg;
    aa_mem_region_init( &reg, 1024*32 );

    /* rfx_trajx_via_t T; */
    /* rfx_trajx_via_init( &T, &reg ); */
    rfx_trajx_parablend_t T;
    rfx_trajx_splend_init( &T, &reg, 1 );
    //rfx_trajx_parablend_init( &T, &reg, 1 );

    rfx_trajx_t *pT = (rfx_trajx_t*)&T;

    double X[4][3] = { {0,0,0}, {1,0,0}, {1,1,0}, {1,1,1} };
    double E[4][3] = { {0,0,0}, {M_PI_2,0,0}, {M_PI_2,M_PI_2,0}, {M_PI_2,M_PI_2,M_PI_2} };
    double R[4][4];
    for( size_t i = 0; i < sizeof(E)/sizeof(E[0]); i ++ ) {
        aa_tf_eulerzyx2quat( E[i][0], E[i][1], E[i][2], R[i] );
        rfx_trajx_add( pT, 5.0*(double)i, X[i], R[i] );
    }

    rfx_trajx_generate( pT );
    rfx_trajx_plot( pT, .001, NULL );

    aa_mem_region_destroy( &reg );
}

int main( void ) {
    /* /rfx_trajq_plot( &traj.traj, 0.01 ); */
    //plotx2();
    //plot_viax();
    plot_qseg();
    //plotq();


}
