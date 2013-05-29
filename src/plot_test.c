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

    rfx_trajq_add( &traj.traj, 0, 0, q0 );
    rfx_trajq_add( &traj.traj, 1, 10, q1 );

    rfx_trajq_generate( &traj.traj );
    printf("times: %f, %f\n", traj.traj.T[0], traj.traj.T[1] );
    printf("blend: %f\n", traj.t_b );
    printf("dq_r: " );
    aa_dump_vec( stdout, traj.dq_r, traj.traj.n_q );
    printf("ddq_r: " );
    aa_dump_vec( stdout, traj.ddq_r, traj.traj.n_q );

    fprintf( stderr, "%lu\n",
             traj.traj.n_q );

    printf("q0: ");
    aa_dump_vec( stdout, traj.traj.Q, traj.traj.n_q );
    printf("q1: ");
    aa_dump_vec( stdout, traj.traj.Q+traj.traj.n_q, traj.traj.n_q );

}


void plotx() {
    aa_mem_region_t reg;
    aa_mem_region_init( &reg, 1024*32 );

    rfx_trajq_trapvel_t trajq;
    rfx_trajq_trapvel_init( &trajq, &reg, 6 );
    for( size_t i = 0; i < 6; i ++ ) {
        trajq.dq_max[i] = 10.0;
        trajq.ddq_max[i] = 10.0;
    }


    double x0[3] = {0.384311,0.004544,0.089198};
    double r0[4] = {0.148775, 0.939373, 0.305129, -0.048381};
    double x1[3] = {0.380000, 0.053318, 0.020952};
    double r1[4] = {0.148778, 0.939347, 0.305212, -0.048341};

    rfx_trajx_t trajx;
    rfx_trajx_rv_init( &trajx, &trajq.traj );

    rfx_trajx_add( &trajx, 0, 0, x0, r0 );
    rfx_trajx_add( &trajx, 1, 10, x1, r1 );

    rfx_trajx_generate( &trajx );

    rfx_trajq_plot( trajx.trajq, .001 );
}

int main( void ) {
    //rfx_trajq_plot( &traj.traj, 0.01 );
    plotx();


}
