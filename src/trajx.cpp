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

#include <amino.hpp>
#include "reflex.h"

using namespace amino;

struct trajx_point {
    double t;
    double tb;
    DualQuat S;
    trajx_point() {}
    trajx_point( double a_t, double a_tb, const double a_S[8] ) :
        t(a_t),
        tb(a_tb),
        S(DualQuat::from_duqu(a_S))
    { }
};

/*--- Point Lists ---*/
struct rfx_trajx_point_list {
    RegionList<trajx_point>::type list;
    rfx_trajx_point_list( aa_mem_region_t *reg ) :
    list(reg)
    {}
};


struct rfx_trajx_point_list *
rfx_trajx_point_list_alloc( aa_mem_region_t *reg )
{
    return new(reg) rfx_trajx_point_list(reg);
}

int
rfx_trajx_point_list_add_qv( struct rfx_trajx_point_list *list, double t, const double r[4], const double v[3] )
{
    rfx_trajx_point_list_addb_qv( list, t, -1, r, v );
}
int
rfx_trajx_point_list_add_duqu( struct rfx_trajx_point_list *list, double t, const double S[8] )
{
    rfx_trajx_point_list_addb_duqu( list, t, -1, S );
}
int
rfx_trajx_point_list_add_tfmat( struct rfx_trajx_point_list *list, double t, const double T[12] )
{
    rfx_trajx_point_list_addb_tfmat( list, t, -1, T );
}


int
rfx_trajx_point_list_addb_qv( struct rfx_trajx_point_list *list,
                              double t, double t_blend, const double r[4], const double v[3] )
{
    DualQuat S( DualQuat::from_qv(r,v) );
    rfx_trajx_point_list_addb_duqu( list, t, t_blend, S.data );
}
int
rfx_trajx_point_list_addb_duqu( struct rfx_trajx_point_list *list,
                                double t, double t_blend, const double S[8] )
{
    list->list.push_back( trajx_point(t, t_blend, S) );
}
int
rfx_trajx_point_list_addb_tfmat( struct rfx_trajx_point_list *list,
                                 double t, double t_blend, const double T[12] )
{
    double S[8];
    aa_tf_tfmat2duqu( T, S );
    rfx_trajx_point_list_addb_duqu( list, t, t_blend, S );
}

/*--- Segment Lists ---*/

struct rfx_trajx_seg_list {
    RegionList<rfx_trajx_seg*>::type list;
    double t_i, t_f;
    RegionList<rfx_trajx_seg*>::iterator itr;
    rfx_trajx_seg_list( aa_mem_region_t *reg ) :
    list(reg)
    {}
};

struct rfx_trajx_seg_list *
rfx_trajx_seg_list_alloc( aa_mem_region_t *reg )
{
    return new(reg) rfx_trajx_seg_list(reg);
}

double
rfx_trajx_seg_list_get_t_i( struct rfx_trajx_seg_list *seg )
{
    return seg->t_i;
}
double
rfx_trajx_seg_list_get_t_f( struct rfx_trajx_seg_list *seg )
{
    return seg->t_f;
}


int
rfx_trajx_seg_list_add( struct rfx_trajx_seg_list *seg_list,
                        struct rfx_trajx_seg *seg )
{
    assert( seg->t_i < seg->t_f );
    if( seg_list->list.empty() ) {
        /* add to empty list */
        seg_list->t_i = seg->t_i;
    } else if (seg_list->t_i > seg->t_f ||
               seg->t_f < seg->t_i ||
               seg->t_f < seg_list->t_f )
    {
        /* Incorrect ordering */
        assert(0);
        return RFX_INVAL;
    }
    seg_list->list.push_back(seg);
    seg_list->itr = seg_list->list.end();
    seg_list->itr--;
    seg_list->t_f = seg->t_f;

    return 0;
}

static rfx_trajx_seg *
seg_search( struct rfx_trajx_seg_list *seglist, double t )
{
    if( t <= seglist->t_i ) return * seglist->list.begin();
    if( t >= seglist->t_f ) return * seglist->list.rbegin();

    /* check current segment */
    if( t >= (*seglist->itr)->t_i &&
        t <= (*seglist->itr)->t_f )
    {
        return *seglist->itr;
    }

    /* check next segment */
    if( seglist->itr != seglist->list.end() ) {
        RegionList<rfx_trajx_seg*>::iterator next = amino::next(seglist->itr);
        if( next != seglist->list.end() &&
            t >= (*next)->t_i &&
            t <= (*next)->t_f )
        {
            seglist->itr = next;
            return *next;
        }
    }

    /* linear search froms start */
    for( RegionList<rfx_trajx_seg*>::iterator itr = seglist->list.begin();
         itr != seglist->list.end();
         itr ++ )
    {
        if( t >= (*itr)->t_i &&
            t <= (*itr)->t_f )
        {
            seglist->itr = itr;
            return *itr;
        }
    }
    /* nothing found */
    return NULL;
}

// TODO: better types for segments

int
rfx_trajx_seg_list_get_x_duqu( struct rfx_trajx_seg_list *seglist,
                               double t, double S[8] )
{
    rfx_trajx_seg *seg = seg_search( seglist, t );
    if( NULL == seg ) return RFX_INVAL;
    double x[3], r[4];
    seg->vtab->get_x( (rfx_trajx_t*)seg, t, x, r );
    aa_tf_qv2duqu( r, x, S );
    return RFX_OK;
}

int
rfx_trajx_seg_list_get_dx_duqu( struct rfx_trajx_seg_list *seglist,
                                double t,  double S[8], double dx[6] )
{
    rfx_trajx_seg *seg = seg_search( seglist, t );
    if( NULL == seg ) return RFX_INVAL;
    double x[3], r[4];
    seg->vtab->get_x( (rfx_trajx_t*)seg, t, x, r );
    seg->vtab->get_dx( (rfx_trajx_t*)seg, t, dx );
    aa_tf_qv2duqu( r, x, S );
    return RFX_OK;
}

int
rfx_trajx_seg_list_get_ddx_duqu( struct rfx_trajx_seg_list *seglist,
                                 double t, double S[8], double dx[6], double ddx[6] )
{
    rfx_trajx_seg *seg = seg_search( seglist, t );
    if( NULL == seg ) return RFX_INVAL;
    double x[3], r[4];
    seg->vtab->get_x( (rfx_trajx_t*)seg, t, x, r );
    seg->vtab->get_dx( (rfx_trajx_t*)seg, t, dx );
    seg->vtab->get_ddx( (rfx_trajx_t*)seg, t, ddx );
    aa_tf_qv2duqu( r, x, S );
    return RFX_OK;
}


int
rfx_trajx_seg_list_get_x_vl( struct rfx_trajx_seg_list *seg,
                             double t, double x[6] )
{
    double S[8];
    int i = rfx_trajx_seg_list_get_x_duqu( seg, t, S );
    if( RFX_OK == i ) {
        double r[4];
        aa_tf_duqu2qv( S, r, x );
        aa_tf_quat2rotvec(r, x+3 );
    }
    return i;
}


int
rfx_trajx_seg_list_get_x_qv( struct rfx_trajx_seg_list *seg,
                             double t, double r[4], double x[3] )
{
    double S[8];
    int i = rfx_trajx_seg_list_get_x_duqu( seg, t, S );
    if( RFX_OK == i ) {
        aa_tf_duqu2qv( S, r, x );
    }
    return i;
}
int
rfx_trajx_seg_list_get_dx_qv( struct rfx_trajx_seg_list *seg,
                              double t, double r[4], double x[3], double dx[6] )
{
    double S[8];
    int i = rfx_trajx_seg_list_get_dx_duqu( seg, t, S, dx );
    if( RFX_OK == i ) {
        aa_tf_duqu2qv( S, r, x );
    }
    return i;
}
int
rfx_trajx_seg_list_get_ddx_qv( struct rfx_trajx_seg_list *seg,
                               double t, double r[4], double x[3], double dx[6], double ddx[6] )
{
    double S[8];
    int i = rfx_trajx_seg_list_get_ddx_duqu( seg, t, S, dx, ddx );
    if( RFX_OK == i ) {
        aa_tf_duqu2qv( S, r, x );
    }
    return i;
}


int
rfx_trajx_seg_list_get_x_tfmat( struct rfx_trajx_seg_list *seg,
                                double t, double T[12] )
{
    double S[8];
    int i = rfx_trajx_seg_list_get_x_duqu( seg, t, S );
    if( RFX_OK == i ) {
        aa_tf_duqu2tfmat( S, T );
    }
    return i;
}
int
rfx_trajx_seg_list_get_dx_tfmat( struct rfx_trajx_seg_list *seg,
                                 double t, double T[12], double dx[6] )
{
    double S[8];
    int i = rfx_trajx_seg_list_get_x_duqu( seg, t, S );
    if( RFX_OK == i ) {
        aa_tf_duqu2tfmat( S, T );
    }
    return i;
}
int
rfx_trajx_seg_list_get_ddx_tfmat( struct rfx_trajx_seg_list *seg,
                                  double t,double T[12], double dx[6], double ddx[6] )
{
    double S[8];
    int i = rfx_trajx_seg_list_get_ddx_duqu( seg, t, S, dx, ddx );
    if( RFX_OK == i ) {
        aa_tf_duqu2tfmat( S, T );
    }
    return i;
}

/*--- Generators ---*/

static void point2vector( const struct trajx_point *pt, double x_last[6], double x[6] )
{
    double r[4];
    aa_tf_duqu2qv( pt->S.data, r, x );
    if( x_last )
        aa_tf_quat2rotvec_near(r, x_last+3, x+3 );
    else
        aa_tf_quat2rotvec(r, x+3 );
}

void virtpoints( RegionList<trajx_point>::type *list,
                 RegionList<trajx_point>::iterator itr_j,
                 struct trajx_point *pt_i, struct trajx_point *pt_j, struct trajx_point *pt_k )

{

    auto itr_0 = list->begin();
    auto itr_1 = amino::next(itr_0);
    auto itr_n = list->end();
    auto itr_k = amino::next(itr_j);
    assert( itr_1 != itr_n ); /* more than one point */
    assert( itr_j != itr_n ); /* not at end */
    /*-- pt_i --*/
    if( itr_j == itr_0 ) {
        *pt_i = *itr_0;
    } else if( itr_j == itr_1 ) {
        *pt_i = *itr_0;
        pt_i->t += pt_i->tb/2;
    } else {
        auto itr_i = amino::prev(itr_j);
        *pt_i = *itr_i;
    }
    /*-- pt_j --*/
    if( itr_j == itr_0 ) {
        *pt_j = *itr_0;
        pt_j->t += pt_j->tb/2;
    } else if( itr_k == itr_n ) {
        *pt_j = *itr_j;
        pt_j->t -= pt_j->tb/2;
    } else {
        *pt_j = *itr_j;
    }
    /*-- pt_k --*/
    if( itr_k == itr_n ) {
        *pt_k = *itr_j;
    } else if( amino::next(itr_k) == itr_n ) {
        *pt_k = *itr_k;
        pt_k->t -= pt_k->tb/2;
    } else {
        *pt_k = *itr_k;
    }
}

struct rfx_trajx_seg_list *
rfx_trajx_parablend_generate( struct rfx_trajx_point_list *points, aa_mem_region_t *reg )
{
    RegionList<trajx_point>::type &plist = points->list;

    if( plist.empty() ) return NULL;
    if( 1 == plist.size() ) return NULL; /* TODO: return constant segment */

    struct rfx_trajx_seg_list * seg_list = rfx_trajx_seg_list_alloc( reg );

    double dx_p[6] = {0};
    auto itr_j = plist.begin();
    seg_list->t_i = itr_j->t;
    seg_list->t_f = itr_j->t;
    while ( plist.end() != itr_j ) {
        struct trajx_point pt_i, pt_j, pt_k;
        virtpoints( &plist, itr_j, &pt_i, &pt_j, &pt_k );
        double x_i[6], x_j[6], x_k[6];
        point2vector( &pt_i, NULL, x_i );
        point2vector( &pt_j, x_i, x_j );
        point2vector( &pt_k, x_j, x_k );

        /* get start position */
        double x_p[6];
        if( itr_j == plist.begin() ) {
            memcpy( x_p, x_i, sizeof(x_p) );
            //memset( dx_p, 0, sizeof(x_p) );
        } else {
            double  r[4];
            //int i = rfx_trajx_seg_list_get_dx_qv( seg_list, seg_list->t_f, r, x_p, dx_p );
            int i = rfx_trajx_seg_list_get_x_qv( seg_list, seg_list->t_f, r, x_p );
            aa_tf_quat2rotvec_near(r, x_j+3, x_p+3 );
        }

        /* compute velocity and acceleration */
        double ddx[6], dx[6];
        for( size_t j = 0; j < 6; j ++ ) {
            dx[j] = (x_k[j] - x_j[j]) / (pt_k.t - pt_j.t);
            ddx[j] = (dx[j] - dx_p[j]) / pt_j.tb;
        }
        /* blend */
        //assert( seg_list->t_f == pt_j.t - pt_j.tb/2 );
        if( rfx_trajx_seg_list_add( seg_list,
                                    rfx_trajx_seg_blend_rv_alloc( reg,
                                                                  seg_list->t_f,
                                                                  seg_list->t_f + pt_j.tb,
                                                                  seg_list->t_f,
                                                                  x_p, dx_p, ddx ) ) )
        {
            return NULL;
        }

        /* linear */
        if( plist.end() != amino::next(itr_j) ) {
            if( rfx_trajx_seg_list_add( seg_list,
                                        rfx_trajx_seg_lerp_rv_alloc( reg,
                                                                     seg_list->t_f,
                                                                     pt_k.t-pt_k.tb/2,
                                                                     pt_j.t, x_j,
                                                                     pt_k.t, x_k ) ) )
            {
                return NULL;
            }
            memcpy(dx_p, dx, sizeof(dx_p));

        }
        itr_j++;
    }

    return seg_list;
}


struct rfx_trajx_seg_list *
rfx_trajx_splend_generate( struct rfx_trajx_point_list *points, aa_mem_region_t *reg )
{
    RegionList<trajx_point>::type &plist = points->list;

    if( plist.empty() ) return NULL;
    if( 1 == plist.size() ) return NULL; /* TODO: return constant segment */

    struct rfx_trajx_seg_list * seg_list = rfx_trajx_seg_list_alloc( reg );

    auto itr_j = plist.begin();
    while ( plist.end() != itr_j ) {
        struct trajx_point pt_i, pt_j, pt_k;
        virtpoints( &plist, itr_j, &pt_i, &pt_j, &pt_k );
        /* Convert Dual Quaternions to quaternion, vector */
        Vec3 v_i(pt_i.S), v_j(pt_j.S), v_k(pt_k.S);
        /* Blend About Current point */
        if( rfx_trajx_seg_list_add(seg_list,
                                   rfx_trajx_seg_blend_q_alloc(reg,
                                                               pt_j.t-pt_j.tb/2,
                                                               pt_j.t+pt_j.tb/2,
                                                               pt_j.tb,
                                                               pt_i.t, v_i.data, pt_i.S.real.data,
                                                               pt_j.t, v_j.data, pt_j.S.real.data,
                                                               pt_k.t, v_k.data, pt_k.S.real.data)) )
        {
            return NULL;
        }

        /* Linear to next point */
        if( plist.end() != amino::next(itr_j) ) {
            double t0 = pt_j.t + pt_j.tb/2;
            double t1 = pt_k.t - pt_k.tb/2;

            if( rfx_trajx_seg_list_add( seg_list,
                                        rfx_trajx_seg_lerp_slerp_alloc( reg, t0, t1,
                                                                        pt_j.t, v_j.data, pt_j.S.real.data,
                                                                        pt_k.t, v_k.data, pt_k.S.real.data)) )
            {
                return NULL;
            }

        }

        itr_j++;
    }

    return seg_list;
}
