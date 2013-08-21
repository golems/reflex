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

struct rfx_trajq_seg_list *
rfx_trajq_seg_list_alloc( aa_mem_region_t *reg ) {
    rfx_trajq_seg_list_t *x = AA_MEM_REGION_NEW( reg, rfx_trajq_seg_list_t );
    memset(x,0,sizeof(*x));
    x->reg = reg;
    x->seg = aa_mem_rlist_alloc( reg );
    // TODO: vtable
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

static inline struct rfx_trajq_seg *
rfx_trajq_seg_init( struct rfx_trajq_seg *seg, size_t n_q,
                    double t_i, double t_f ) {
    seg->n_q = n_q;
    seg->t_i = t_i;
    seg->t_f = t_f;
    return seg;
}

#define RFX_TRAJQ_SEG_ALLOC( reg, type, n_q, t_i, t_f )                 \
    ( (type*) rfx_trajq_seg_init( (struct rfx_trajq_seg *) AA_MEM_REGION_NEW( (reg), type), \
                                  (n_q), (t_i), (t_f) ) )


struct rfx_trajq_seg_dq *
rfx_trajq_seg_dq_alloc( aa_mem_region_t *reg, size_t n_q,
                        double t_i, double *q_i, double *dq,
                        double t_f ) {
    struct rfx_trajq_seg_dq *x = RFX_TRAJQ_SEG_ALLOC( reg, struct rfx_trajq_seg_dq, n_q, t_i, t_f );
    x->q_i = q_i;
    x->dq = dq;
    // TODO: vtable
    return x;
}

struct rfx_trajq_seg_dq *
rfx_trajq_seg_dq_alloc2( aa_mem_region_t *reg, size_t n_q,
                         double t_i, double *q_i,
                         double t_f, double *q_f ) {
    struct rfx_trajq_seg_dq *x = RFX_TRAJQ_SEG_ALLOC( reg, struct rfx_trajq_seg_dq, n_q, t_i, t_f );
    x->q_i = q_i;

    /* Allocate velocity buffer after the struct, so popping the
       struct pops this buffer */
    double *dq = AA_MEM_REGION_NEW_N( reg, double, n_q );

    // computer velocity
    for( size_t i = 0; i < n_q; i ++ ) {
        dq[i] = (q_f[i] - q_i[i]) / (t_f - t_i);
    }

    x->dq = dq;
    return x;
}

struct rfx_trajq_seg_2dq *
rfx_trajq_seg_2dq_alloc( aa_mem_region_t *reg, size_t n_q,
                         double t_i, double *q_i,
                         double *dq, double *ddq,
                         double t_f ) {
    struct rfx_trajq_seg_2dq *x = RFX_TRAJQ_SEG_ALLOC( reg, struct rfx_trajq_seg_2dq, n_q, t_i, t_f );
    x->q_i = q_i;
    x->dq_i = dq;
    x->ddq = ddq;
    // TODO: vtable
    return x;
}

struct rfx_trajq_seg_2dq *
rfx_trajq_seg_2dq_alloc2( aa_mem_region_t *reg, size_t n_q,
                          double t_i, double *q_i, double *dq_i,
                          double t_f, double *q_f );

struct rfx_trajq_seg_3dq *
rfx_trajq_seg_3dq_alloc( aa_mem_region_t *reg, size_t n_q,
                         double t_i, double *q_i,
                         double *dq, double *ddq, double *dddq,
                         double t_f ) {
    struct rfx_trajq_seg_3dq *x = RFX_TRAJQ_SEG_ALLOC( reg, struct rfx_trajq_seg_3dq, n_q, t_i, t_f );
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
