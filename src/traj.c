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


int rfx_trajx_set_ctrl( struct rfx_trajx *cx, double t, rfx_ctrlx_lin_t *ctrlx ) {
    int r;

    double x_r[3], r_r[4];
    if( (r = rfx_trajx_get_x( cx, t, x_r,r_r )) )
        return r;
    aa_tf_qv2duqu( r_r, x_r, ctrlx->ctrl->ref.S );

    r = rfx_trajx_get_dx( cx, t, ctrlx->ctrl->ref.dx );

    return r;
}

static void x_add( struct rfx_trajx *cx, double t, double x[3], double r[4] ) {
    struct rfx_trajx_point *pt = AA_MEM_REGION_NEW( cx->reg, struct rfx_trajx_point );
    pt->t = t;
    AA_MEM_CPY( pt->x, x, 3 );
    AA_MEM_CPY( pt->r, r, 4 );
    cx->pt_f = pt;
    if( NULL == cx->point->head ) cx->pt_i = pt;
    aa_mem_rlist_enqueue_ptr( cx->point, pt );
    cx->n_p++;
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
    memset(traj,0,sizeof(*traj));
    struct rfx_trajq_trapvel *trajq = (struct rfx_trajq_trapvel*) aa_mem_region_alloc( reg, sizeof(*trajq) );
    traj->reg = reg;
    traj->trajq = &trajq->traj;
    rfx_trajq_trapvel_init( trajq, reg, 6 );
    traj->vtab = &vtab_x_rv;

    traj->point = aa_mem_rlist_alloc( reg );

    for( size_t i = 0; i < 6; i ++ ) {
        trajq->dq_max[i] = 10.0; // TODO: pick better
        trajq->ddq_max[i] = 10.0;
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
    memset(traj,0,sizeof(*traj));
    struct rfx_trajq_trapvel *trajq = (struct rfx_trajq_trapvel*) aa_mem_region_alloc( reg, sizeof(*trajq) );
    traj->trajq = &trajq->traj;
    traj->reg = reg;

    rfx_trajq_trapvel_init( trajq, reg, 4 );
    traj->vtab = &vtab_x_slerp;

    traj->point = aa_mem_rlist_alloc( reg );

    for( size_t i = 0; i < 3; i ++ ) {
        trajq->dq_max[i] = 10.0; // TODO: pick better
        trajq->ddq_max[i] = 10.0;
    }
    trajq->dq_max[3] = 50;
    trajq->ddq_max[3] = 50;
}

/*-- VIA --*/

static int x_seg_lerp_get_x( struct rfx_trajx *cx, double t, double x[3], double r[4] ) {
    rfx_trajx_seg_lerp_t *S = (rfx_trajx_seg_lerp_t*)cx;
    double u = (t - S->tau_i) / S->dt;
    aa_la_d_lerp( 3, u,
                  S->x_i, 1,
                  S->x_f, 1,
                  x, 1 );
    aa_tf_qslerp( u, S->r_i, S->r_f, r );
    return 0;
}

static int x_seg_lerp_get_dx( struct rfx_trajx *cx, double t, double dx[6] ) {
    rfx_trajx_seg_lerp_t *S = (rfx_trajx_seg_lerp_t*)cx;
    /* translational */
    for( size_t i = 0; i < 3; i ++ ) {
        dx[i] = (S->x_f[i] - S->x_i[i]) / S->dt;
    }

    /* rotational */
    double u = (t - S->tau_i) / S->dt;
    double r[4], dr[4];
    aa_tf_qslerp( u, S->r_i, S->r_f, r );
    aa_tf_qslerpdiff( u, S->r_i, S->r_f, dr );
    // dr/dt = dr/du * du/dt
    for( size_t i = 0; i < 4; i ++ ) dr[i] /= S->dt;
    aa_tf_qdiff2vel( r, dr, dx+3 );
    return 0;
}

static struct rfx_trajx_vtab x_seg_lerp_vtab = {
    .get_x = x_seg_lerp_get_x,
    .get_dx = x_seg_lerp_get_dx
};

struct rfx_trajx_seg *
rfx_trajx_seg_lerp_slerp_alloc( aa_mem_region_t *reg, double t_i, double t_f,
                                double tau_i, double x_i[3], double r_i[4],
                                double tau_f, double x_f[3], double r_f[4] ) {
    rfx_trajx_seg_lerp_t *S = AA_MEM_REGION_NEW( reg, rfx_trajx_seg_lerp_t );
    S->seg.vtab = &x_seg_lerp_vtab;
    S->seg.t_i = t_i;
    S->seg.t_f = t_f;
    S->tau_i = tau_i;
    S->tau_f = tau_f;
    S->dt = tau_f - tau_i;
    AA_MEM_CPY( S->x_i, x_i, 3 );
    AA_MEM_CPY( S->x_f, x_f, 3 );
    AA_MEM_CPY( S->r_i, r_i, 4 );
    AA_MEM_CPY( S->r_f, r_f, 4 );
    return &S->seg;
}



static int x_via_generate( struct rfx_trajx *cx ) {
    cx->pt_i = (struct rfx_trajx_point*) cx->point->head->data;
    struct rfx_trajx_via *traj = (struct rfx_trajx_via*)cx;

    // allocate arrays
    traj->n_seg = cx->n_p - 1;
    traj->seg = AA_MEM_REGION_NEW_N( cx->reg, struct rfx_trajx_seg*, traj->n_seg );

    // fill arrays
    size_t i = 0;
    for( aa_mem_cons_t *pcons = cx->point->head; pcons; pcons = pcons->next )
    {
        struct rfx_trajx_point *pt = (struct rfx_trajx_point*)pcons->data;
        struct rfx_trajx_point *pt_next = NULL;
        if( NULL == pcons->next ) cx->pt_f = pt;
        else {
            pt_next = (struct rfx_trajx_point*)pcons->next->data;
            traj->seg[i] = rfx_trajx_seg_lerp_slerp_alloc( traj->trajx.reg, pt->t, pt_next->t,
                                                           pt->t, pt->x, pt->r,
                                                           pt_next->t, pt_next->x, pt_next->r );
        }
        i++;
    }
    return 0;
}

static int x_via_compar( const void *a, const void *b ) {
    double t = *(double*)a;
    struct rfx_trajx_seg *B = *(struct rfx_trajx_seg**)b;
    if( t < B->t_i ) return -1;
    else if( t > B->t_f ) return 1;
    else {
        return 0;
    }
}

static struct rfx_trajx_seg  *
x_via_search( struct rfx_trajx *cx, double t ) {
    rfx_trajx_via_t *traj = (rfx_trajx_via_t*)cx;
    if( t <= cx->pt_i->t ) return traj->seg[0];
    else if (t >= cx->pt_f->t) return traj->seg[traj->n_seg-1];
    /* else search */
    struct rfx_trajx_seg *B = *(struct rfx_trajx_seg**)bsearch( &t, traj->seg,
                                                                traj->n_seg, sizeof(traj->seg[0]),
                                                                x_via_compar );
    return B;
}

static int x_via_get_x( struct rfx_trajx *cx, double t, double x[3], double r[4] ) {
    struct rfx_trajx_seg *s = x_via_search(cx,t);
    return rfx_trajx_get_x( (rfx_trajx_t*)s, t, x, r );
}

static int x_via_get_dx( struct rfx_trajx *cx, double t, double dx[6] ) {
    struct rfx_trajx_seg *s = x_via_search(cx,t);
    return rfx_trajx_get_dx( (rfx_trajx_t*)s, t, dx );
}

static int x_via_get_ddx( struct rfx_trajx *cx, double t, double ddx[6] ) {
    struct rfx_trajx_seg *s = x_via_search(cx,t);
    return rfx_trajx_get_ddx( (rfx_trajx_t*)s, t, ddx );
}

static void x_via_add( struct rfx_trajx *cx, double t, double x[3], double r[4] ) {
    x_add( cx, t, x, r );
}

static struct rfx_trajx_vtab vtab_x_via = {
    .generate = x_via_generate,
    .add = x_via_add,
    .get_x = x_via_get_x,
    .get_dx = x_via_get_dx,
    .get_ddx = x_via_get_ddx,
};

void rfx_trajx_via_init( struct rfx_trajx_via *traj, aa_mem_region_t *reg ) {
    memset(traj,0,sizeof(*traj));
    traj->trajx.vtab = &vtab_x_via;
    traj->trajx.reg = reg;
    traj->trajx.point = aa_mem_rlist_alloc( reg );

    traj->trajx.n_p = 0;
}

/*-- Parabolic Blends --*/

// Rotation Vector Linear Segment
static int x_seg_lerp_rv_get_x( struct rfx_trajx *cx, double t, double x[3], double r[4] ) {
    rfx_trajx_seg_lerp_rv_t *S = (rfx_trajx_seg_lerp_rv_t*)cx;
    double u = (t - S->tau_i) / S->dt;
    double xp[6];
    aa_la_d_lerp( 6, u,
                  S->x_i, 1,
                  S->x_f, 1,
                  xp, 1 );
    AA_MEM_CPY(x, xp, 3);
    aa_tf_rotvec2quat( xp+3, r );
    return 0;
}

static int x_seg_lerp_rv_get_dx( struct rfx_trajx *cx, double t, double dx[6] ) {
    rfx_trajx_seg_lerp_rv_t *S = (rfx_trajx_seg_lerp_rv_t*)cx;
    double dv[3], v[3];
    for( size_t i = 0; i < 3; i ++ ) {
        dx[i] = (S->x_f[i] - S->x_i[i]) / S->dt;
        dv[i] = (S->x_f[i+3] - S->x_i[i+3]) / S->dt;
        v[i]  = S->x_i[i+3] + dv[i]*(t-S->tau_i);
    }

    aa_tf_rotvec_diff2vel( v, dv, dx+3 );

    return 0;
}

static struct rfx_trajx_vtab x_seg_lerp_rv_vtab = {
    .get_x = x_seg_lerp_rv_get_x,
    .get_dx = x_seg_lerp_rv_get_dx
};


struct rfx_trajx_seg *
rfx_trajx_seg_lerp_rv_alloc( aa_mem_region_t *reg, double t_i, double t_f,
                             double tau_i, double x_i[6],
                             double tau_f, double x_f[6] ) {
    rfx_trajx_seg_lerp_rv_t *S = AA_MEM_REGION_NEW( reg, rfx_trajx_seg_lerp_rv_t );
    S->seg.vtab = &x_seg_lerp_rv_vtab;
    S->seg.t_i = t_i;
    S->seg.t_f = t_f;
    S->tau_i = tau_i;
    S->tau_f = tau_f;
    S->dt = tau_f - tau_i;
    AA_MEM_CPY( S->x_i, x_i, 6 );
    AA_MEM_CPY( S->x_f, x_f, 6 );
    return &S->seg;
}

// Rotation Vector Parabolic Blend Segment
static int x_seg_blend_rv_get_x( struct rfx_trajx *cx, double t, double x[3], double r[4] ) {
    rfx_trajx_seg_blend_rv_t *S = (rfx_trajx_seg_blend_rv_t*)cx;
    double dt = t - S->tau_i;
    double xp[6];
    for( size_t i = 0; i < 6; i ++ ) {
        xp[i] = S->x_i[i] + dt*S->dx_i[i] + 0.5*dt*dt*S->ddx[i];
    }
    AA_MEM_CPY(x, xp, 3);
    aa_tf_rotvec2quat( xp+3, r );
    return 0;
}

static int x_seg_blend_rv_get_dx( struct rfx_trajx *cx, double t, double dx[6] ) {
    rfx_trajx_seg_blend_rv_t *S = (rfx_trajx_seg_blend_rv_t*)cx;
    double dt = t - S->tau_i;
    double dv[3], v[3];
    for( size_t i = 0; i < 3; i ++ ) {
        dv[i] = S->dx_i[i+3] + dt*S->ddx[i+3];
        dx[i] = S->dx_i[i] + dt*S->ddx[i];
        v[i] = S->x_i[i+3] + dt*S->dx_i[i+3] + 0.5*dt*dt*S->ddx[i+3];
    }
    aa_tf_rotvec_diff2vel( v, dv, dx+3 );
    return 0;
}

static struct rfx_trajx_vtab x_seg_blend_rv_vtab = {
    .get_x = x_seg_blend_rv_get_x,
    .get_dx = x_seg_blend_rv_get_dx
};


struct rfx_trajx_seg *
rfx_trajx_seg_blend_rv_alloc( aa_mem_region_t *reg, double t_i, double t_f,
                              double tau_i,
                              double x_i[6], double dx_i[6], double ddx[6] ) {
    rfx_trajx_seg_blend_rv_t *S = AA_MEM_REGION_NEW( reg, rfx_trajx_seg_blend_rv_t );
    S->seg.vtab = &x_seg_blend_rv_vtab;
    S->seg.t_i = t_i;
    S->seg.t_f = t_f;
    S->tau_i = tau_i;
    AA_MEM_CPY( S->x_i, x_i, 6 );
    AA_MEM_CPY( S->dx_i, dx_i, 6 );
    AA_MEM_CPY( S->ddx, ddx, 6 );
    return &S->seg;
}

// static void point2vector( struct rfx_trajx_point *pt, double x_last[6], double x[6] ) {
//     AA_MEM_CPY( x, pt->x, 3);
//     aa_tf_quat2rotvec_near(pt->r, x_last+3, x+3 );
// }

// static int x_parablend_generate( struct rfx_trajx *cx ) {
//     cx->pt_i = (struct rfx_trajx_point*) cx->point->head->data;
//     struct rfx_trajx_parablend *traj = (struct rfx_trajx_parablend*)cx;

//     // allocate arrays
//     traj->via.n_seg = cx->n_p * 2 - 1;
//     traj->via.seg = AA_MEM_REGION_NEW_N( cx->reg, struct rfx_trajx_seg*, traj->via.n_seg );

//     // add virtual via points
//     {
//         struct rfx_trajx_point *pt_i = AA_MEM_REGION_NEW_CPY( cx->reg, cx->pt_i, struct rfx_trajx_point );
//         struct rfx_trajx_point *pt_f = AA_MEM_REGION_NEW_CPY( cx->reg, cx->pt_f, struct rfx_trajx_point );
//         aa_mem_rlist_push_ptr( cx->point, pt_i );
//         aa_mem_rlist_enqueue_ptr( cx->point, pt_f );
//         cx->pt_i->t += traj->t_b/2;
//         cx->pt_f->t -= traj->t_b/2;
//         cx->pt_i = pt_i;
//         cx->pt_f = pt_f;
//     }

//     // fill arrays
//     size_t i = 0;
//     double dx_p[6] = {0}, x_p[6] = {0};
//     {
//         struct rfx_trajx_point *pt0 = (struct rfx_trajx_point*)cx->point->head->data;
//         AA_MEM_CPY( x_p, pt0->x, 3);
//         aa_tf_quat2rotvec(pt0->r, x_p+3 );
//     }
//     for( aa_mem_cons_t *pcons = cx->point->head->next; pcons && pcons->next; pcons = pcons->next )
//     {
//         double t_b = traj->t_b;
//         /* current point */
//         struct rfx_trajx_point *pt = (struct rfx_trajx_point*)pcons->data;
//         double x[6];
//         point2vector( pt, x_p, x );

//         /* next point point */
//         struct rfx_trajx_point *pt_next = (struct rfx_trajx_point*)pcons->next->data;
//         double x_n[6];
//         point2vector( pt_next, x, x_n );


//         /* compute velocity and acceleration */
//         double ddx[6], dx[6] = {0};
//         for( size_t j = 0; j < 6; j ++ ) {
//             dx[j] = (x_n[j] - x[j]) / (pt_next->t - pt->t);
//             ddx[j] = (dx[j] - dx_p[j]) / t_b;
//         }
//         /* add blend about current point */
//         traj->via.seg[i++] = rfx_trajx_seg_blend_rv_alloc( cx->reg, pt->t-t_b/2, pt->t+t_b/2,
//                                                            pt->t-t_b/2, x_p, dx_p, ddx );
//         /* add linear to next point */
//         if( pcons->next->next ) {
//             traj->via.seg[i++] = rfx_trajx_seg_lerp_rv_alloc( cx->reg, pt->t+t_b/2, pt_next->t-t_b/2,
//                                                               pt->t, x,
//                                                               pt_next->t, x_n );
//             /* previous is end of current linear segment */
//             double r[4];
//             rfx_trajx_get_x( (struct rfx_trajx*)traj->via.seg[i-1], pt_next->t-t_b/2, x_p, r );
//             aa_tf_quat2rotvec_near(r, x+3, x_p+3 );
//             /* copy current to prev */
//             AA_MEM_CPY(dx_p, dx, 6);
//         }
//     }
//     return 0;
// }


// static struct rfx_trajx_vtab vtab_x_parablend = {
//     .generate = x_parablend_generate,
//     .add = x_via_add,
//     .get_x = x_via_get_x,
//     .get_dx = x_via_get_dx,
//     .get_ddx = x_via_get_ddx,
// };

// void rfx_trajx_parablend_init( struct rfx_trajx_parablend *cx, aa_mem_region_t *reg, double t_b ) {
//     memset(cx,0,sizeof(*cx));
//     rfx_trajx_via_init(&cx->via, reg);
//     cx->t_b = t_b;
//     cx->via.trajx.vtab = &vtab_x_parablend;
// }


/*--- Cartesian Sphereical Parabolic Blends ---*/

static void x_seg_blend_q_get_u( struct rfx_trajx *cx, double t,
                                 double *pdt,
                                 double *du_ij, double *u_ij,
                                 double *du_jk, double *u_jk,
                                 double *du_j, double *u_j )
{
    rfx_trajx_seg_blend_q_t *S = (rfx_trajx_seg_blend_q_t*)cx;
    double t_ij = S->t_ij;
    double t_b = S->t_b;
    double ddu_ij = S->ddu_ij;
    double ddu_jk = S->ddu_jk;

    double dt = t - (S->t_j - t_b/2);

    *du_ij = (1.0 / t_ij) + (dt * ddu_ij);
    *u_ij = (t_ij - t_b/2) / t_ij  +  dt/t_ij + 0.5*ddu_ij * dt*dt;

    *du_jk = dt * ddu_jk;
    *u_jk = 0.5*ddu_jk * dt*dt;

    if (t < S->t_j) {
        *du_j = S->ddu_j * dt;
        *u_j = 0.5 * S->ddu_j * dt*dt;
    } else {
        double dtt = (S->t_j + t_b/2) - t;
        *du_j = S->ddu_j * dtt;
        *u_j = 1 - 0.5 * S->ddu_j * dtt*dtt;
    }

    *pdt = dt;

}

static int x_seg_blend_q_get_q( struct rfx_trajx *cx, double t, double r[4], double dr[4] ) {
    rfx_trajx_seg_blend_q_t *S = (rfx_trajx_seg_blend_q_t*)cx;
    double dt, du_ij, u_ij, du_jk, u_jk, du_j, u_j;
    x_seg_blend_q_get_u( cx, t,
                         &dt,
                         &du_ij, &u_ij,
                         &du_jk, &u_jk,
                         &du_j, &u_j );

    aa_tf_qslerp3diff( u_ij, du_ij, S->r_i, S->r_j,
                       u_jk, du_jk, S->r_j, S->r_k,
                       u_j, du_j, r, dr );
    return 0;
}


// Rotation Vector Parabolic Blend Segment
static int x_seg_blend_q_get_x( struct rfx_trajx *cx, double t, double x[3], double r[4] ) {
    rfx_trajx_seg_blend_q_t *S = (rfx_trajx_seg_blend_q_t*)cx;
    double dt = t - S->tau_i;

    /* translation */
    for( size_t i = 0; i < 3; i ++ ) {
        x[i] = S->x0[i] + dt*S->dx0[i] + 0.5*dt*dt*S->ddx[i];
    }

    /* rotation */
    double  dq[4];
    x_seg_blend_q_get_q( cx, t, r, dq );
    return 0;
}

static int x_seg_blend_q_get_dx( struct rfx_trajx *cx, double t, double dx[6] ) {
    rfx_trajx_seg_blend_q_t *S = (rfx_trajx_seg_blend_q_t*)cx;
    double dt = t - (S->t_j-S->t_b/2);

    /* translation */
    for( size_t i = 0; i < 3; i ++ ) {
        dx[i] = S->dx0[i] + dt*S->ddx[i];
    }

    /* rotation */
    double  q[4], dq[4];
    x_seg_blend_q_get_q( cx, t, q, dq );
    aa_tf_qdiff2vel( q, dq, dx+3 );

    return 0;
}

static struct rfx_trajx_vtab x_seg_blend_q_vtab = {
    .get_x = x_seg_blend_q_get_x,
    .get_dx = x_seg_blend_q_get_dx
};


struct rfx_trajx_seg *
rfx_trajx_seg_blend_q_alloc( aa_mem_region_t *reg, double t_0, double t_1,
                             double t_b,
                             double t_i, double x_i[3], double r_i[4],
                             double t_j, double x_j[3], double r_j[4],
                             double t_k, double x_k[3], double r_k[4] )
{
    rfx_trajx_seg_blend_q_t *S = AA_MEM_REGION_NEW( reg, rfx_trajx_seg_blend_q_t );
    S->seg.vtab = &x_seg_blend_q_vtab;
    S->seg.t_i = t_0;
    S->seg.t_f = t_1;

    S->t_i = t_i;
    S->t_j = t_j;
    S->t_k = t_k;
    S->t_b = t_b;

    AA_MEM_CPY( S->r_i, r_i, 4 );
    AA_MEM_CPY( S->r_j, r_j, 4 );
    AA_MEM_CPY( S->r_k, r_k, 4 );

    S->t_ij = t_j-t_i;
    S->tau_i = t_j - t_b/2;
    S->ddu_ij =  - 1.0 / (S->t_ij*t_b);
    S->ddu_jk =  1.0 / ( (t_k-t_j)*t_b );
    S->ddu_j = 4 / (t_b*t_b);

    /* translational parameters */
    aa_la_d_lerp( 3, (S->t_ij-t_b/2)/S->t_ij,
                  x_i, 1,
                  x_j, 1,
                  S->x0, 1 );
    double dx1[3];
    for( size_t i = 0; i < 3; i ++ ) {
        S->dx0[i] = (x_j[i] - x_i[i]) / S->t_ij;
        dx1[i] = (x_k[i] - x_j[i]) / (t_k - t_j);
        S->ddx[i] = (dx1[i] - S->dx0[i]) / t_b;
    }
    return &S->seg;
}

// static int x_splend_generate( struct rfx_trajx *cx ) {
//     rfx_trajx_splend_t *traj = (rfx_trajx_splend_t*)cx;

//     // allocate arrays
//     traj->via.n_seg = cx->n_p * 2 - 1;
//     traj->via.seg = AA_MEM_REGION_NEW_N( cx->reg, struct rfx_trajx_seg*, traj->via.n_seg );

//     // add virtual via points
//     {
//         cx->pt_i = (struct rfx_trajx_point*) cx->point->head->data;
//         struct rfx_trajx_point *pt_i = AA_MEM_REGION_NEW_CPY( cx->reg, cx->pt_i, struct rfx_trajx_point );
//         struct rfx_trajx_point *pt_f = AA_MEM_REGION_NEW_CPY( cx->reg, cx->pt_f, struct rfx_trajx_point );
//         aa_mem_rlist_push_ptr( cx->point, pt_i );
//         aa_mem_rlist_enqueue_ptr( cx->point, pt_f );
//         cx->pt_i->t += traj->t_b/2;
//         cx->pt_f->t -= traj->t_b/2;
//         cx->pt_i = pt_i;
//         cx->pt_f = pt_f;
//     }

//     // fill arrays
//     size_t i = 0;
//     for( aa_mem_cons_t *pcons = cx->point->head; pcons->next->next; pcons = pcons->next )
//     {
//         double t_b = traj->t_b;
//         /* current point */
//         struct rfx_trajx_point *pt_i = (struct rfx_trajx_point*)pcons->data;
//         struct rfx_trajx_point *pt_j = (struct rfx_trajx_point*)pcons->next->data;
//         struct rfx_trajx_point *pt_k = (struct rfx_trajx_point*)pcons->next->next->data;
//         /* add blend about current point */
//         traj->via.seg[i++] = rfx_trajx_seg_blend_q_alloc( cx->reg, pt_j->t-t_b/2, pt_j->t+t_b/2,
//                                                           t_b,
//                                                           pt_i->t, pt_i->x, pt_i->r,
//                                                           pt_j->t, pt_j->x, pt_j->r,
//                                                           pt_k->t, pt_k->x, pt_k->r );
//         /* add linear to next point */
//         if( pcons->next->next->next ) {
//             double t0 = pt_j->t + t_b/2;
//             double t1 = pt_k->t - t_b/2;
//             // FIXME: this is wrong for the first one
//             traj->via.seg[i++] = rfx_trajx_seg_lerp_slerp_alloc( cx->reg, t0, t1,
//                                                                  pt_j->t, pt_j->x, pt_j->r,
//                                                                  pt_k->t, pt_k->x, pt_k->r );
//             double x[3], r[4];
//             rfx_trajx_get_x( (rfx_trajx_t*)traj->via.seg[i-1], t0, x, r );
//         }
//     }
//     return 0;
// }


// static struct rfx_trajx_vtab vtab_x_splend = {
//     .generate = x_splend_generate,
//     .add = x_via_add,
//     .get_x = x_via_get_x,
//     .get_dx = x_via_get_dx,
//     .get_ddx = x_via_get_ddx,
// };

// void rfx_trajx_splend_init( struct rfx_trajx_parablend *cx, aa_mem_region_t *reg, double t_b ) {
//     memset(cx,0,sizeof(*cx));
//     rfx_trajx_via_init(&cx->via, reg);
//     cx->t_b = t_b;
//     cx->via.trajx.vtab = &vtab_x_splend;
// }
