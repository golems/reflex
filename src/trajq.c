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

rfx_trajq_points_t *
rfx_trajq_points_alloc( aa_mem_region_t *reg, size_t n_q ) {
    rfx_trajq_points_t *x = AA_MEM_REGION_NEW( reg, rfx_trajq_points_t );
    memset(x,0,sizeof(*x));
    x->reg = reg;
    x->n_q = n_q;
    x->point = aa_mem_rlist_alloc( reg );
    return x;
}

void rfx_trajq_points_add( rfx_trajq_points_t *pts, double t, double *q ) {
    //printf("add point: %f\n", t );
    double *qq = AA_MEM_REGION_NEW_N( pts->reg, double, pts->n_q+1 );
    qq[0] = t;
    AA_MEM_CPY(qq+1, q, pts->n_q);
    aa_mem_rlist_enqueue_ptr( pts->point, qq );

    if( 0 == pts->n_t ) {
        pts->t_i = t;
        pts->q_i = qq+1;
    }

    pts->t_f = t;
    pts->q_f = qq+1;
    pts->n_t++;
}


/*--- Segment List ---*/

static struct rfx_trajq_seg*
seg_list_lookup( struct rfx_trajq_seg_list *list, double t ) {
    //printf("list: 0x%x\n", list );
    //printf("last_seg: 0x%x\n", list->last_seg );
    //printf("list->seg: 0x%x\n", list->seg );
    //printf("list->seg->head: 0x%x\n", list->seg->head );
    // init cached segment
    // TODO: handle empty list?
    if( ! list->last_seg ) list->last_seg = list->seg->head;

    // check cached segment
    rfx_trajq_seg_t *s_l = (struct rfx_trajq_seg*) list->last_seg->data;
    if( t >= s_l->t_i && t <= s_l->t_f )
        return s_l;

    // check next segment
    aa_mem_cons_t *next = list->last_seg->next;
    if( next ) {
        rfx_trajq_seg_t *s_n = (struct rfx_trajq_seg*) next->data;
        if( t >= s_n->t_i && t <= s_n->t_f ) {
            list->last_seg = next;
            return s_n;
        }
    } else if ( t > s_l->t_f ) {
        return s_l;
    }

    // linear search
    aa_mem_cons_t *pcons = list->seg->head;
    while( pcons->next ) {
        rfx_trajq_seg_t *s_n = (struct rfx_trajq_seg*) pcons->data;
        if( t <= s_n->t_f ) break;
        pcons = pcons->next;
    }
    return (struct rfx_trajq_seg*)pcons->data;
}

static int
seg_list_get_q( void *list, double t, double *q ) {
    return rfx_trajq_seg_get_q( seg_list_lookup((struct rfx_trajq_seg_list*)list,t),
                                t, q );
}

static int
seg_list_get_dq( void *list, double t, double *q, double *dq ) {
    return rfx_trajq_seg_get_dq( seg_list_lookup((struct rfx_trajq_seg_list*)list,t),
                                 t, q, dq );
}

static int
seg_list_get_ddq( void *list, double t, double *q, double *dq, double *ddq ) {
    return rfx_trajq_seg_get_ddq( seg_list_lookup((struct rfx_trajq_seg_list*)list,t),
                                  t, q, dq, ddq );
}

static struct rfx_trajq_seg_vtab seglist_vtab = {
    .get_q   = seg_list_get_q,
    .get_dq  = seg_list_get_dq,
    .get_ddq = seg_list_get_ddq
};

struct rfx_trajq_seg_list *
rfx_trajq_seg_list_alloc( aa_mem_region_t *reg ) {
    rfx_trajq_seg_list_t *x = AA_MEM_REGION_NEW( reg, rfx_trajq_seg_list_t );
    memset(x,0,sizeof(*x));
    x->reg = reg;
    x->seg = aa_mem_rlist_alloc( reg );
    x->vtab = &seglist_vtab;
    return x;
}

void rfx_trajq_seg_list_add( rfx_trajq_seg_list_t *seglist, rfx_trajq_seg_t *seg ) {
    aa_mem_rlist_enqueue_ptr( seglist->seg, seg );

    if( 0 == seglist->n_t ) {
        seglist->t_i = seg->t_i;
    }

    seglist->t_f = seg->t_f;
    seglist->n_t++;
}



/*--- Parabolic Blends ---*/
rfx_trajq_seg_list_t *
rfx_trajq_gen_pblend_tm1( aa_mem_region_t *reg, rfx_trajq_points_t *points, double t_blend ) {
    struct rfx_trajq_seg_list *list = rfx_trajq_seg_list_alloc(reg);
    size_t n_q = list->n_q = points->n_q;
    double tb_2 = t_blend/2;
    struct rfx_trajq_point *pt_prev = NULL;

    double dq_prev[n_q];
    memset(dq_prev, 0, sizeof(dq_prev));

    struct rfx_trajq_point *pt_virt0 =  (struct rfx_trajq_point*) alloca(sizeof(double)*(1+n_q));
    struct rfx_trajq_point *pt_virt1 =  (struct rfx_trajq_point*) alloca(sizeof(double)*(1+n_q));

    for( aa_mem_cons_t *pcons = points->point->head; pcons; pcons = pcons->next ) {
        struct rfx_trajq_point *pt;
        struct rfx_trajq_point *pt_next;
        pt = (struct rfx_trajq_point*)pcons->data;

        //printf("checking point: %f\n", pt->t );

        // initial point
        if( pt_prev ) {
            pt = (struct rfx_trajq_point*)pcons->data;
        } else {
            pt_prev = (struct rfx_trajq_point*)pcons->data;
            pt = pt_virt0;
            pt->t = pt_prev->t+tb_2;
            AA_MEM_CPY( pt->q, pt_prev->q, n_q );
        }

        // next point
        if( NULL == pcons->next ) {
            // final point
            pt_next = pt;
            AA_MEM_CPY( pt_virt1->q, pt_next->q, n_q );
            pt_virt1->t = pt_next->t - t_blend;
            pt = pt_virt1;
        } else {
            pt_next = (struct rfx_trajq_point*)pcons->next->data;
            if( NULL == pcons->next->next ) {
                // penultimate point
                AA_MEM_CPY( pt_virt1->q, pt_next->q, n_q );
                pt_virt1->t = pt_next->t - tb_2;
                pt_next = pt_virt1;
            }
        }

        // velocity and acceleration
        double dq_next[n_q];
        double ddq[n_q];
        for( size_t i = 0; i < n_q; i++ ) {
            dq_next[i] = (pt_next->q[i] - pt->q[i]) / (pt_next->t - pt->t);
            ddq[i] = (dq_next[i] - dq_prev[i]) / t_blend;
        }

        //printf("[%f, %f, %f]\n", pt_prev->t, pt->t, pt_next->t);
        //printf("dq_next: "); aa_dump_vec( stdout, dq_next, n_q );
        //printf("ddq: "); aa_dump_vec( stdout, ddq, n_q );

        int r = 0;

        // acceleration segment
        if( 0 == list->n_t ) {
            // initial point
            rfx_trajq_seg_list_add( list,
                                    (rfx_trajq_seg_t *)
                                    rfx_trajq_seg_2dq_alloc( reg, n_q,
                                                             pt_prev->t, pt_prev->q,
                                                             dq_prev, ddq, pt->t+tb_2 ));
        } else if (pcons->next) {
            // interior point
            r = rfx_trajq_seg_2dq_link( list, ddq, pt->t+tb_2);
        } else {
            r = rfx_trajq_seg_2dq_link( list, ddq, pt->t + t_blend);
        }
        if( r ) return NULL;

        // linear segment
        if( pcons->next )
            if( (r = rfx_trajq_seg_dq_link( list, dq_next, pt_next->t-tb_2 )) )
                return NULL;


        pt_prev = pt;
        memcpy( dq_prev, dq_next, sizeof(dq_prev) );
    }

    //printf("segs: %lu\n", list->n_t);
    return list;
}


        /* if( NULL == pt_prev ) { */
        /*     /\* initial segment *\/ */
        /*     for( size_t i = 0; i < n_q; i++ ) { */
        /*         dq_next[i] = (pt_next->q[i] - pt->q[i]) / (pt_next->t - pt->t - 0.5*t_blend); */
        /*         ddq[i] = dq_next[i]  / t_blend; */
        /*     } */


/*--- Segments ---*/
static inline struct rfx_trajq_seg *
rfx_trajq_seg_init( struct rfx_trajq_seg *seg, struct rfx_trajq_seg_vtab *vtab,
                    size_t n_q, double t_i, double t_f ) {
    //printf("init n_q: %lu\n", n_q );
    seg->vtab = vtab;
    seg->n_q = n_q;
    seg->t_i = t_i;
    seg->t_f = t_f;
    return seg;
}

// Allocate trajectory segment of type with n_p additional space
#define RFX_TRAJQ_SEG_ALLOC( reg, type, vtab, n_q, t_i, t_f, n_p )      \
    ( (type*) rfx_trajq_seg_init( (struct rfx_trajq_seg *)              \
                                  aa_mem_region_alloc( (reg), sizeof(type) + (n_p)), \
                                  (vtab), (n_q), (t_i), (t_f) ) )



/*--- Constant Velocity ---*/

static inline void dq_parm( rfx_trajq_seg_dq_t *cx, double t,
                            size_t *n_q, double **p_q_i, double **p_dq, double *dt ) {
    *n_q = cx->parent.n_q;
    *p_q_i = cx->p;
    *p_dq = cx->p + *n_q;
    *dt = t - cx->parent.t_i;
}

static int dq_get_q( void *vcx, double t, double *q ) {
    rfx_trajq_seg_dq_t *cx = (rfx_trajq_seg_dq_t*) vcx;
    size_t n_q ;
    double *p_q_i, *p_dq, dt;
    dq_parm( cx, t, &n_q, &p_q_i, &p_dq, &dt );
    for( size_t i = 0; i < n_q; i++ ) {
        q[i] = p_q_i[i] + dt * p_dq[i];
    }
    return 0;
}

static int dq_get_dq( void *vcx, double t, double *q, double *dq ) {
    rfx_trajq_seg_dq_t *cx = (rfx_trajq_seg_dq_t*) vcx;
    size_t n_q ;
    double *p_q_i, *p_dq, dt;
    dq_parm( cx, t, &n_q, &p_q_i, &p_dq, &dt );
    for( size_t i = 0; i < n_q; i++ ) {
        q[i] = p_q_i[i] + dt * p_dq[i];
        dq[i] = p_dq[i];
    }
    return 0;
}

static int dq_get_ddq( void *vcx, double t, double *q, double *dq, double *ddq ) {
    rfx_trajq_seg_dq_t *cx = (rfx_trajq_seg_dq_t*) vcx;
    size_t n_q ;
    double *p_q_i, *p_dq, dt;
    dq_parm( cx, t, &n_q, &p_q_i, &p_dq, &dt );
    for( size_t i = 0; i < n_q; i++ ) {
        q[i] = p_q_i[i] + dt * p_dq[i];
        dq[i] = p_dq[i];
        ddq[i] = 0;
    }
    return 0;
}

static struct rfx_trajq_seg_vtab seg_dq_vtab = {
    .get_q   = dq_get_q,
    .get_dq  = dq_get_dq,
    .get_ddq = dq_get_ddq
};

rfx_trajq_seg_dq_t *
rfx_trajq_seg_dq_alloc( aa_mem_region_t *reg, size_t n_q,
                        double t_i, double *q_i, double *dq,
                        double t_f ) {
    rfx_trajq_seg_dq_t *x = RFX_TRAJQ_SEG_ALLOC( reg, rfx_trajq_seg_dq_t, &seg_dq_vtab,
                                                 n_q, t_i, t_f, 2*n_q*sizeof(double) );
    double *p_q_i = x->p;
    double *p_dq = x->p + n_q;
    AA_MEM_CPY(p_q_i, q_i, n_q );
    AA_MEM_CPY(p_dq, dq, n_q );
    return x;
}

rfx_trajq_seg_dq_t *
rfx_trajq_seg_dq_alloc2( aa_mem_region_t *reg, size_t n_q,
                         double t_i, double *q_i,
                         double t_f, double *q_f ) {
    // compute velocity
    double dq[n_q];
    for( size_t i = 0; i < n_q; i ++ ) {
        dq[i] = (q_f[i] - q_i[i]) / (t_f - t_i);
    }
    // allocate
    return rfx_trajq_seg_dq_alloc( reg, n_q, t_i, q_i, dq, t_f );
}

int rfx_trajq_seg_dq_link( rfx_trajq_seg_list_t *seglist, double *dq, double t_f ) {
    // initial state
    double t_i = seglist->t_f;
    if( t_f <= t_i ) return -1;
    double q_i[seglist->n_q];
    rfx_trajq_seg_list_get_q( seglist, t_i, q_i );
    // add segment
    rfx_trajq_seg_list_add( seglist,
                            (rfx_trajq_seg_t*)rfx_trajq_seg_dq_alloc( seglist->reg, seglist->n_q,
                                                                      t_i, q_i, dq, t_f ) );
    return 0;
}

/*--- Constant Acceleration ---*/

static inline void x_2dq_parm( rfx_trajq_seg_dq_t *cx, double t,
                               size_t *n_q, double **p_q_i, double **p_dq, double **p_ddq,
                               double *dt, double *dt2_2 ) {
    *n_q = cx->parent.n_q;
    *p_q_i = cx->p;
    *p_dq = cx->p + *n_q;
    *p_ddq = cx->p + *n_q + *n_q;
    *dt = t - cx->parent.t_i;
    *dt2_2 = (*dt)*(*dt) / 2;
}


static int x_2dq_get_q( void *vcx, double t, double *q ) {
    rfx_trajq_seg_2dq_t *cx = (rfx_trajq_seg_dq_t*) vcx;
    size_t n_q ;
    double *p_q_i, *p_dq, *p_ddq, dt, dt2_2;
    x_2dq_parm( cx, t, &n_q, &p_q_i, &p_dq, &p_ddq, &dt, &dt2_2 );

    for( size_t i = 0; i < n_q; i++ ) {
        q[i] = p_q_i[i] + dt * p_dq[i] + dt2_2 * p_ddq[i] ;
    }
    return 0;
}

static int x_2dq_get_dq( void *vcx, double t, double *q, double *dq ) {
    rfx_trajq_seg_2dq_t *cx = (rfx_trajq_seg_dq_t*) vcx;
    size_t n_q ;
    double *p_q_i, *p_dq, *p_ddq, dt, dt2_2;
    x_2dq_parm( cx, t, &n_q, &p_q_i, &p_dq, &p_ddq, &dt, &dt2_2 );

    for( size_t i = 0; i < n_q; i++ ) {
        q[i] = p_q_i[i] + dt * p_dq[i] + dt2_2 * p_ddq[i] ;
        dq[i] = p_dq[i] + dt*p_ddq[i];
    }
    return 0;
}

static int x_2dq_get_ddq( void *vcx, double t, double *q, double *dq, double *ddq ) {
    rfx_trajq_seg_2dq_t *cx = (rfx_trajq_seg_dq_t*) vcx;
    size_t n_q ;
    double *p_q_i, *p_dq, *p_ddq, dt, dt2_2;
    x_2dq_parm( cx, t, &n_q, &p_q_i, &p_dq, &p_ddq, &dt, &dt2_2 );

    //printf("n_q: %lu\n", cx->n_q );
    //aa_dump_vec( stdout, p_q_i, cx->n_q );
    //aa_dump_vec( stdout, p_dq, cx->n_q );
    //aa_dump_vec( stdout, p_ddq, cx->n_q );

    for( size_t i = 0; i < n_q; i++ ) {
        q[i] = p_q_i[i] + dt * p_dq[i] + dt2_2 * p_ddq[i] ;
        dq[i] = p_dq[i] + dt*p_ddq[i];
        ddq[i] = p_ddq[i];
    }
    return 0;
}

static struct rfx_trajq_seg_vtab seg_2dq_vtab = {
    .get_q   = x_2dq_get_q,
    .get_dq  = x_2dq_get_dq,
    .get_ddq = x_2dq_get_ddq
};


rfx_trajq_seg_2dq_t *
rfx_trajq_seg_2dq_alloc( aa_mem_region_t *reg, size_t n_q,
                         double t_i, double *q_i,
                         double *dq, double *ddq,
                         double t_f ) {
    rfx_trajq_seg_2dq_t *x = RFX_TRAJQ_SEG_ALLOC( reg, rfx_trajq_seg_2dq_t, &seg_2dq_vtab,
                                                  n_q, t_i, t_f, 3*n_q*sizeof(double) );
    double *p_q_i = x->p;
    double *p_dq =  x->p + n_q;
    double *p_ddq = x->p + n_q + n_q;
    AA_MEM_CPY(p_q_i, q_i, n_q );
    AA_MEM_CPY(p_dq,  dq,  n_q );
    AA_MEM_CPY(p_ddq, ddq, n_q );
    return x;
}

rfx_trajq_seg_2dq_t *
rfx_trajq_seg_2dq_alloc2( aa_mem_region_t *reg, size_t n_q,
                          double t_i, double *q_i, double *dq_i,
                          double t_f, double *q_f );


int rfx_trajq_seg_2dq_link( rfx_trajq_seg_list_t *seglist, double *ddq, double t_f ) {
    // initial state
    double t_i = seglist->t_f;
    if( t_f <= t_i ) return -1;
    double q_i[seglist->n_q];
    double dq_i[seglist->n_q];
    rfx_trajq_seg_list_get_dq( seglist, t_i, q_i, dq_i );
    // add segment
    rfx_trajq_seg_list_add( seglist,
                            (rfx_trajq_seg_t*)rfx_trajq_seg_2dq_alloc( seglist->reg, seglist->n_q,
                                                                       t_i, q_i, dq_i, ddq, t_f ) );
    return 0;
}

/*--- Constant Jerk ---*/

static struct rfx_trajq_seg_vtab seg_3dq_vtab = {
    .get_q   = NULL,
    .get_dq  = NULL,
    .get_ddq = NULL
};

struct rfx_trajq_seg_3dq *
rfx_trajq_seg_3dq_alloc( aa_mem_region_t *reg, size_t n_q,
                         double t_i, double *q_i,
                         double *dq, double *ddq, double *dddq,
                         double t_f ) {
    struct rfx_trajq_seg_3dq *x = RFX_TRAJQ_SEG_ALLOC( reg, struct rfx_trajq_seg_3dq,
                                                       &seg_3dq_vtab, n_q, t_i, t_f, 3*n_q*sizeof(double) );
    x->q_i = q_i;
    x->dq_i = dq;
    x->ddq_i = ddq;
    x->dddq = dddq;
    // TODO: vtable
    return x;
}


struct rfx_trajq_seg_3dq *
rfx_trajq_seg_3dq_alloc2( aa_mem_region_t *reg, size_t n_q,
                          double t_i, double *q_i,
                          double *dq, double *ddq,
                          double t_f, double *dq_f, double *ddq_f );


/*--- Dense waypoint segments ---*/

typedef struct rfx_trajq_seg_dense {
    rfx_trajq_seg_t parent;
    size_t n_p;
    double dt;
    double *Q;
} rfx_trajq_seg_dense_t;


static int dense_get_ddq( void *vcx, double t, double *q, double *dq, double *ddq ) {
    rfx_trajq_seg_dense_t *cx = (rfx_trajq_seg_dense_t*) vcx;
    size_t n_q = cx->parent.n_q;
    size_t n_p = cx->n_p;
    double dt = cx->dt;
    double *Q = cx->Q;
    double t_i = cx->parent.t_i;

    if( t <= cx->parent.t_i ) {
        AA_MEM_CPY( q, cx->Q, n_q );
        AA_MEM_SET(dq, 0, n_q );
        AA_MEM_SET(ddq, 0, n_q );
    } else  if( t >= cx->parent.t_i ) {
        AA_MEM_CPY( q, cx->Q + n_q*(n_p-1), n_q );
        AA_MEM_SET(dq, 0, n_q );
        AA_MEM_SET(ddq, 0, n_q );
    } else {
        size_t j = (size_t) ( (t - (double)cx->parent.t_i) / dt );
        double dj = (double)j;
        /* get velocities of previous and next points */
        double dx0[n_q], dx1[n_q];
        double *X0 = Q+n_q*j;
        double *X1 = Q+n_q*(j+1);
        if( 0 == j ) {
            memset(dx0, 0, sizeof(dx0));
        } else  {
            aa_la_quadterp_dx( n_q,
                               t_i+dt*(dj-1), Q+n_q*(j-1),
                               t_i+dt*dj,     X0,
                               t_i+dt*(dj+1), X1,
                               t_i+dt*dj,     dx0 );

        }
        if( n_p-2 <= j ) {
            memset(dx1, 0, sizeof(dx0));
        } else {
            aa_la_quadterp_dx( n_q,
                               t_i+dt*dj,     X0,
                               t_i+dt*(dj+1), X1,
                               t_i+dt*(dj+2), Q+n_q*(j+2),
                               t_i+dt*(dj+2), dx0 );
        }
        /* Compute interpolation parameters */
        double a2[n_q], a3[n_q];
        aa_la_d_3spline_param( n_q, dt, /* assumes x0 is at t=0 */
                               Q+n_q*j,     1,  dx0, 1,
                               Q+n_q*(j+1), 1, dx1, 1,
                               a2, a3 );
        /* Interpolate */
        aa_la_d_3spline( n_q, t - t_i + dj*dt,
                         X0, 1, dx0, 1, a2, a3,
                         q, 1, dq, 1, ddq, 1 );

    }


    return 0;
}

static int dense_get_q( void *vcx, double t, double *q ) {
    rfx_trajq_seg_dense_t *cx = (rfx_trajq_seg_dense_t*) vcx;
    size_t n_q = cx->parent.n_q;
    double dq[n_q];
    double ddq[n_q];
    dense_get_ddq( vcx, t, q, dq, ddq );

    return 0;
}

static int dense_get_dq( void *vcx, double t, double *q, double *dq ) {
    rfx_trajq_seg_dense_t *cx = (rfx_trajq_seg_dense_t*) vcx;
    size_t n_q = cx->parent.n_q;
    double ddq[n_q];
    dense_get_ddq( vcx, t, q, dq, ddq );

    return 0;
}


static struct rfx_trajq_seg_vtab seg_dense_vtab = {
    .get_q   = dense_get_q,
    .get_dq  = dense_get_dq,
    .get_ddq = dense_get_ddq
};

struct rfx_trajq_seg*
rfx_trajq_seg_dense_alloc( aa_mem_region_t *reg, size_t n_q, size_t n_p,
                           double t_i, double dt,
                           double *Q )
{

    struct rfx_trajq_seg_dense *x =
        RFX_TRAJQ_SEG_ALLOC( reg, struct rfx_trajq_seg_dense,
                             &seg_dense_vtab, n_q, t_i, t_i+(double)n_p*dt, 2*sizeof(double)+sizeof(size_t) );
    x->n_p = n_p;
    x->dt = dt;
    x->Q = Q;
    return &x->parent;
}
