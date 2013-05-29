/* -*- mode: C++; c-basic-offset: 4 -*- */
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


void rfx_trajq_init( struct rfx_trajq *cx, aa_mem_region_t *reg, size_t n_q, size_t n_t ) {
    cx->n_q = n_q;
    cx->n_t = n_t;
    cx->reg = reg;

    cx->T = (double*)aa_mem_region_alloc( cx->reg, cx->n_t * sizeof(cx->T[0]) );
    cx->Q = (double*)aa_mem_region_alloc( cx->reg, cx->n_t * cx->n_q * sizeof(cx->Q[0]) );
}

void rfx_trajq_destroy( struct rfx_trajq *cx ) {
    aa_mem_region_pop(cx->reg, cx->T);
}

void rfx_trajq_add( struct rfx_trajq *cx, size_t i, double t, double *q ) {
    cx->T[i] = t;
    memcpy( &cx->Q[i*cx->n_q], q, cx->n_q*sizeof(q[0]) );
}

int rfx_trapvel_generate( size_t n, double t_f,
                          const double *x_i, const double *x_f,
                          const double *dx_max, const double *ddx_max,
                          double *ptb, double *dx_r, double *ddx_r ) {
    double t3=t_f;
    double tb = t3/2;
    //fprintf(stderr, "t3: %f\n", t3 );

    double x[n];
    aa_la_vsub(n, x_f, x_i, x);
    //aa_dump_vec(stderr, x, 6 );

    // try triangular profile
    int is_tri = 1;
    for( size_t i = 0; i < n; i++ ) {
        dx_r[i] = x[i] / tb;
        // check for velocity limit
        if( fabs(dx_r[i]) > fabs(dx_max[i]) ) {
            is_tri = 0;
            break;
        }
        // check for acceleration limit
        ddx_r[i] = dx_r[i] / tb;
        if( fabs(ddx_r[i]) > fabs(ddx_max[i]) ) {
            fprintf(stderr, "accel limit\n");
            return -1;
        }
    }
    if( is_tri ) {
        //fprintf(stderr, "tri\n");
        *ptb = tb;
        return 0;
    }
    //fprintf(stderr, "trap\n");

    // needs to be trapezoid
    // find longest acceptable blend time
    for( size_t i = 0; i < n; i++ ) {
        double t = t3 - x[i]/dx_max[i];
        tb = AA_MIN(tb,t);
    }
    // calc dx, ddx
    double t2 = t3 - tb;
    for( size_t i = 0; i < n; i++ ) {
        dx_r[i] = x[i]/t2;
        ddx_r[i] = dx_r[i]/tb;
        // check a
        if( fabs(ddx_r[i]) > fabs(ddx_max[i]) ||
            fabs(dx_r[i])  > fabs(dx_max[i]) ) {
            return -1;
        }
    }
    *ptb = tb;
    return 0;
}



int q_trapvel_generate( void *vcx ) {
    struct rfx_trajq_trapvel *cx = (struct rfx_trajq_trapvel*)vcx;
    double delta_t = cx->traj.T[1] - cx->traj.T[0];

    rfx_trapvel_generate( cx->traj.n_q, delta_t,
                          cx->traj.Q, cx->traj.Q + cx->traj.n_q,
                          cx->dq_max, cx->ddq_max,
                          &cx->t_b, cx->dq_r, cx->ddq_r );
    return 0;
}

int q_trapvel_get_q( void *vcx, double t, double *q ) {
    struct rfx_trajq_trapvel *cx = (struct rfx_trajq_trapvel*)vcx;
    const double *q0 = cx->traj.Q;
    const double *q1 = cx->traj.Q+cx->traj.n_q;
    if( t < cx->traj.T[0] ) {
        // before t0
        memcpy( q, q0, cx->traj.n_q*sizeof(q[0]) ) ;
    }
    else if( t > cx->traj.T[1] ) {
        // after t1
        memcpy( q, q1, cx->traj.n_q*sizeof(q[0]) );
    } else {
        //normal
        double t1 = cx->traj.T[0] + cx->t_b;
        double t2 = cx->traj.T[1] - cx->t_b; // switching timesreturn
        if( t < t1 ) {
            double tt = t - cx->traj.T[0];
            for( size_t i = 0; i < cx->traj.n_q; i ++ )
                q[i] = q0[i] + 0.5 * cx->ddq_r[i] * tt * tt;
        } else if (t < t2 ) {
            double tt = t - t1;
            for( size_t i = 0; i < cx->traj.n_q; i ++ )
                q[i] = q0[0] + 0.5*cx->dq_r[i]*cx->t_b + cx->dq_r[i]*tt;
        } else if (t < cx->traj.T[1] ) {
            double tt = cx->traj.T[1] - t;
            for( size_t i = 0; i < cx->traj.n_q; i ++ )
                q[i] = q1[i] - .5*cx->ddq_r[i]*tt*tt;
        } else {
            assert(0);
        }
    }

    return 0;
}

int q_trapvel_get_dq( void *vcx, double t, double *dq ) {
    struct rfx_trajq_trapvel *cx = (struct rfx_trajq_trapvel*)vcx;
    double t_i = cx->traj.T[0], t_f = cx->traj.T[1];
    if( t <= t_i || t >= t_f ) {
        memset( dq, 0, cx->traj.n_q*sizeof(dq[0]) );
    } else if( t <= t_i + cx->t_b ) {
        // accelerating blend
        for( size_t i = 0; i < cx->traj.n_q; i ++ )
            dq[i] = (t-t_i) * cx->ddq_r[i];
    } else if (t <= t_f - cx->t_b ) {
        // constant velocity
        memcpy( dq, cx->dq_r, cx->traj.n_q*sizeof(dq[0]) );
    } else if (t <= t_f ) {
        // deccelerating blend
        double t2 = t_f - cx->t_b;
        for( size_t i = 0; i < cx->traj.n_q; i ++ )
            dq[i] = cx->dq_r[i] - (t-t2)*cx->ddq_r[i]; // steady state - decelleration
    } else {
        // bogus
        assert(0);
    }
    return 0;
}

int q_trapvel_get_ddq( void *vcx, double t, double *ddq ) {
    struct rfx_trajq_trapvel *cx = (struct rfx_trajq_trapvel*)vcx;
    double t_i = cx->traj.T[0], t_f = cx->traj.T[1];
    if( t <= t_i || t >= t_f ) {
        memset( ddq, 0, cx->traj.n_q * sizeof(ddq[0]) );
    } else if (t <= t_i + cx->t_b ) {
        // accelerating blend
        memcpy( ddq, cx->ddq_r, cx->traj.n_q * sizeof(ddq[0]) );
    } else if (t <= t_f - cx->t_b ) {
        // constant velocity
        memset( ddq, 0, cx->traj.n_q * sizeof(ddq[0]) );
    } else if (t <= t_f ) {
        // deccelerating blend
        for( size_t i = 0; i < cx->traj.n_q; i ++ )
            ddq[i] = - cx->ddq_r[i];
    } else {
        assert(0);
    }

    return 0;
}

static struct rfx_trajq_vtab vtab_q_trapvel = {
    q_trapvel_generate,
    q_trapvel_get_q,
    q_trapvel_get_dq,
    q_trapvel_get_ddq
};

void rfx_trajq_trapvel_init( struct rfx_trajq_trapvel *cx, aa_mem_region_t *reg, size_t n_q ) {
    rfx_trajq_init( &cx->traj, reg, n_q, 2 );
    cx->dq_r    = (double*)aa_mem_region_alloc( reg, n_q*sizeof(cx->dq_r[0]) );
    cx->ddq_r   = (double*)aa_mem_region_alloc( reg, n_q*sizeof(cx->ddq_r[0]) );
    cx->dq_max  = (double*)aa_mem_region_alloc( reg, n_q*sizeof(cx->dq_max[0]) );
    cx->ddq_max = (double*)aa_mem_region_alloc( reg, n_q*sizeof(cx->ddq_max[0]) );
    cx->traj.vtab = &vtab_q_trapvel;

}
