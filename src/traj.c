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

void rfx_trajq_init( struct rfx_trajq *traj, aa_mem_region_t *reg, size_t n_q ) {
    traj->n_q = n_q;
    traj->reg = reg;
    traj->n_t = 0;
    traj->point = aa_mem_rlist_alloc( reg );
}

struct rfx_trajq *rfx_trajq_alloc( aa_mem_region_t *reg, size_t n_q ) {
    struct rfx_trajq* traj = AA_MEM_REGION_NEW( reg, struct rfx_trajq );
    rfx_trajq_init( traj, reg, n_q );
    return traj;
}

/* void rfx_trajq_alloc( struct rfx_trajq *cx, aa_mem_region_t *reg, size_t n_q, size_t n_t ) { */
/*     cx->n_q = n_q; */
/*     cx->n_t = n_t; */

/*     cx->T = (double*)aa_mem_region_alloc( cx->reg, cx->n_t * sizeof(cx->T[0]) ); */
/*     cx->Q = (double*)aa_mem_region_alloc( cx->reg, cx->n_t * cx->n_q * sizeof(cx->Q[0]) ); */
/* } */

/* void rfx_trajq_destroy( struct rfx_trajq *cx ) { */
/*     aa_mem_region_pop(cx->reg, cx->T); */
/* } */

void rfx_trajq_add( struct rfx_trajq *cx,  double t, double *q ) {
    double *qq = AA_MEM_REGION_NEW_N( cx->reg, double, cx->n_q+1 );
    qq[0] = t;
    AA_MEM_CPY(qq+1, q, cx->n_q);
    aa_mem_rlist_enqueue_ptr( cx->point, qq );

    if( 0 ==  cx->n_t ) {
        cx->t_i = t;
        cx->q_i = qq+1;
    }
    cx->t_f = t;
    cx->q_f = qq+1;
    cx->n_t++;
    printf("addq: "); aa_dump_vec( stdout, qq, cx->n_q+1 );
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

    if( 2 != cx->traj.n_t ) {
        abort();
        return -1;
    }
    /* double *X_i = (double*)cx->traj.point->head->data; */
    /* double *X_f = (double*)cx->traj.point->head->next->data; */
    /* cx->t_i = X_i[0]; */
    /* cx->t_f = X_f[0]; */
    /* cx->q_i = X_i+1; */
    /* cx->q_f = X_f+1; */

    double delta_t = cx->traj.t_f - cx->traj.t_i;

    rfx_trapvel_generate( cx->traj.n_q, delta_t,
                          cx->traj.q_i, cx->traj.q_f,
                          cx->dq_max, cx->ddq_max,
                          &cx->t_b, cx->dq_r, cx->ddq_r );
    return 0;
}

int q_trapvel_get_q( void *vcx, double t, double *q ) {
    struct rfx_trajq_trapvel *cx = (struct rfx_trajq_trapvel*)vcx;
    const double *q0 = cx->traj.q_i;
    const double *q1 = cx->traj.q_f;
    if( t < cx->traj.t_i ) {
        // before t0
        memcpy( q, q0, cx->traj.n_q*sizeof(q[0]) ) ;
    }
    else if( t > cx->traj.t_f ) {
        // after t1
        memcpy( q, q1, cx->traj.n_q*sizeof(q[0]) );
    } else {
        //normal
        double t1 = cx->traj.t_i + cx->t_b;
        double t2 = cx->traj.t_f - cx->t_b; // switching timesreturn
        if( t < t1 ) {
            double tt = t - cx->traj.t_i;
            for( size_t i = 0; i < cx->traj.n_q; i ++ )
                q[i] = q0[i] + 0.5 * cx->ddq_r[i] * tt * tt;
        } else if (t < t2 ) {
            double tt = t - t1;
            for( size_t i = 0; i < cx->traj.n_q; i ++ )
                q[i] = q0[0] + 0.5*cx->dq_r[i]*cx->t_b + cx->dq_r[i]*tt;
        } else if (t < cx->traj.t_f ) {
            double tt = cx->traj.t_f - t;
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
    double t_i = cx->traj.t_i, t_f = cx->traj.t_f;
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
    double t_i = cx->traj.t_i, t_f = cx->traj.t_f;
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
    .generate = q_trapvel_generate,
    .get_q = q_trapvel_get_q,
    .get_dq = q_trapvel_get_dq,
    .get_ddq = q_trapvel_get_ddq
};

void rfx_trajq_trapvel_init( struct rfx_trajq_trapvel *cx, aa_mem_region_t *reg, size_t n_q ) {
    rfx_trajq_init( &cx->traj, reg, n_q );
    cx->dq_r    = AA_MEM_REGION_NEW_N( reg, double, n_q );
    cx->ddq_r   = AA_MEM_REGION_NEW_N( reg, double, n_q );
    cx->dq_max  = AA_MEM_REGION_NEW_N( reg, double, n_q );
    cx->ddq_max = AA_MEM_REGION_NEW_N( reg, double, n_q );
    cx->traj.vtab = &vtab_q_trapvel;

}



/*************/
/* WORKSPACE */
/*************/

static void x_add( struct rfx_trajx *cx, double t, double x[3], double r[4] ) {
    struct rfx_trajx_point *pt = AA_MEM_REGION_NEW( cx->trajq->reg, struct rfx_trajx_point );
    pt->t = t;
    AA_MEM_CPY( pt->x, x, 3 );
    AA_MEM_CPY( pt->r, r, 4 );
    aa_mem_rlist_enqueue_ptr( cx->point, pt );
}


/*-- Rotation Vector --*/

static int x_rv_generate( struct rfx_trajx *cx ) {

    cx->pt_i = (struct rfx_trajx_point*) cx->point->head->data;
    cx->pt_f = (struct rfx_trajx_point*) cx->point->head->next->data;

    double rp[3];
    for( aa_mem_cons_t *cons = cx->point->head; cons; cons = cons->next ) {
        double xp[6];
        struct rfx_trajx_point *pt = (struct rfx_trajx_point*)cons->data;
        AA_MEM_CPY( xp, pt->x, 3 );
        if( cons == cx->point->head ) {
            aa_tf_quat2rotvec( pt->r, xp+3 );
        } else {
            aa_tf_quat2rotvec_near( pt->r, rp, xp+3 );
        }
        AA_MEM_CPY(rp, xp, 3);
        rfx_trajq_add( cx->trajq, pt->t, xp );
    }

    int i = cx->trajq->vtab->generate(cx->trajq);

    return i;
}


static int x_rv_get_x( struct rfx_trajx *cx, double t, double x[3], double r[4] ) {
    double xp[6];
    int i = rfx_trajq_get_q( cx->trajq, t, xp );
    if( !i ) {
        memcpy( x, xp, 3*sizeof(x[0]) );
        aa_tf_rotvec2quat( xp+3, r );
    }
    return i;
}
static int x_rv_get_dx( struct rfx_trajx *cx, double t, double dx[6] ) {
    return rfx_trajq_get_dq( cx->trajq, t, dx );
}
static int x_rv_get_ddx( struct rfx_trajx *cx, double t, double ddx[6] ) {
    return rfx_trajq_get_ddq( cx->trajq, t, ddx );
}

static struct rfx_trajx_vtab vtab_x_rv = {
    .generate = x_rv_generate,
    .add = x_add,
    .get_x = x_rv_get_x,
    .get_dx = x_rv_get_dx,
    .get_ddx = x_rv_get_ddx,
};

void rfx_trajx_destroy( struct rfx_trajx *cx ) {
    aa_mem_region_pop(cx->trajq->reg, cx->trajq);
}

void rfx_trajx_rv_init( struct rfx_trajx *traj, aa_mem_region_t *reg ) {
    struct rfx_trajq_trapvel *trajq = (struct rfx_trajq_trapvel*) aa_mem_region_alloc( reg, sizeof(*trajq) );
    traj->trajq = &trajq->traj;
    rfx_trajq_trapvel_init( trajq, reg, 6 );
    traj->vtab = &vtab_x_rv;

    traj->point = aa_mem_rlist_alloc( reg );

    for( size_t i = 0; i < 6; i ++ ) {
        trajq->dq_max[i] = 1.0; // TODO: pick better
        trajq->ddq_max[i] = 1.0;
    }

}

/*-- SLERP --*/

static int x_slerp_generate( struct rfx_trajx *cx ) {

    cx->pt_i = (struct rfx_trajx_point*) cx->point->head->data;
    cx->pt_f = (struct rfx_trajx_point*) cx->point->head->next->data;

    int i = 0;
    for( aa_mem_cons_t *cons = cx->point->head; cons; cons = cons->next ) {
        struct rfx_trajx_point *pt = (struct rfx_trajx_point*)cons->data;
        double xp[4];
        AA_MEM_CPY( xp, pt->x, 3 );

        if( 0 == i ) {
            xp[3] = 0;
        } else if( 1 == i ) {
            xp[3] = 1;
        } else {
            abort();
        }
        rfx_trajq_add( cx->trajq, pt->t, xp );
        i++;
    }

    return cx->trajq->vtab->generate(cx->trajq);
}

static int x_slerp_get_x( struct rfx_trajx *cx, double t, double x[3], double r[4] ) {
    double xp[4];
    int i = rfx_trajq_get_q( cx->trajq, t, xp );


    //if( !i ) {
        memcpy( x, xp, 3*sizeof(x[0]) );
        aa_tf_qslerp( xp[3], cx->pt_i->r, cx->pt_f->r, r );
    //}
        //printf("slerp %f: ", xp[3]); aa_dump_vec( stdout, r, 4 );
    return i;
}

static int x_slerp_get_dx( struct rfx_trajx *cx, double t, double dx[6] ) {
    double xp[4];
    int i = rfx_trajq_get_dq( cx->trajq, t, xp );

    memcpy( dx, xp, 3*sizeof(xp[0]) );

    double tau = xp[3];
    double delta_t = cx->pt_f->t - cx->pt_i->t;
    double dtau_dt = 1.0 / delta_t;
    double dr_dtau[4], r[4], dr_dt[4];
    aa_tf_qslerp( xp[3], cx->pt_i->r, cx->pt_f->r, r );
    aa_tf_qslerpdiff( tau, cx->pt_i->r, cx->pt_f->r, dr_dtau );
    for( size_t j = 0; j < 4; j ++ ) dr_dt[j] = dr_dtau[j] * dtau_dt;
    aa_tf_qdiff2vel( r, dr_dt, dx+3 );

    return i;
}
static int x_slerp_get_ddx( struct rfx_trajx *cx, double t, double ddx[6] ) {
    (void)cx; (void)t; (void) ddx;
    /* double xp[4]; */
    /* int i = rfx_trajq_get_ddq( cx->trajq, t, xp ); */
    /* double tau = xp[3]; */
    /* memcpy( ddx, xp, 3*sizeof(xp[0]) ); */
    abort();
    return 0;
}

static struct rfx_trajx_vtab vtab_x_slerp = {
    .generate = x_slerp_generate,
    .add = x_add,
    .get_x = x_slerp_get_x,
    .get_dx = x_slerp_get_dx,
    .get_ddx = x_slerp_get_ddx,
};


void rfx_trajx_slerp_init( struct rfx_trajx *traj, aa_mem_region_t *reg ) {
    struct rfx_trajq_trapvel *trajq = (struct rfx_trajq_trapvel*) aa_mem_region_alloc( reg, sizeof(*trajq) );
    traj->trajq = &trajq->traj;
    rfx_trajq_trapvel_init( trajq, reg, 4 );
    traj->vtab = &vtab_x_slerp;

    traj->point = aa_mem_rlist_alloc( reg );

    for( size_t i = 0; i < 3; i ++ ) {
        trajq->dq_max[i] = 1.0; // TODO: pick better
        trajq->ddq_max[i] = 1.0;
    }
    trajq->dq_max[3] = 50;
    trajq->ddq_max[3] = 50;
}
