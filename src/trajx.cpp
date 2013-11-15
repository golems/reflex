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
    double S[8];
    trajx_point( double a_t, double a_tb, const double a_S[8] ) :
        t(a_t),
        tb(a_tb)
    {
        memcpy(S, a_S, sizeof(S) );
    }
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
rfx_trajx_point_list_add_xr( struct rfx_trajx_point_list *list, double t, const double x[3], const double r[4] )
{
    rfx_trajx_point_list_addb_xr( list, t, -1, x, r );
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
rfx_trajx_point_list_addb_xr( struct rfx_trajx_point_list *list,
                              double t, double t_blend, const double x[3], const double r[4] )
{
    double S[8];
    aa_tf_qv2duqu( r, x, S );
    rfx_trajx_point_list_addb_duqu( list, t, t_blend, S );
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

int
rfx_trajx_point_list_add( struct rfx_trajx_seg_list *seg_list,
                          struct rfx_trajx_seg *seg )
{
    if( seg_list->list.empty() ) {
        /* add to empty list */
        seg_list->t_i = seg->t_i;
    } else if (seg_list->t_i < seg->t_f ||
               seg->t_f < seg->t_i ||
               seg->t_f < seg_list->t_f )
    {
        /* Incorrect ordering */
        return RFX_INVAL;
    }
    seg_list->list.push_back(seg);
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
        RegionList<rfx_trajx_seg*>::iterator next = seglist->itr;
        next++;
        if( t >= (*next)->t_i &&
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
}

/*--- Generators ---*/

// static void point2vector( const struct trajx_point *pt, double x_last[6], double x[6] )
// {
//     double r[4];
//     aa_tf_duqu2qv( pt->S, r, x );
//     aa_tf_quat2rotvec_near(r, x_last+3, x+3 );
// }

// static void
//  RegionList<rfx_trajx_seg*>::type &slist

// void parablend_add_blend( aa_mem_region_t *reg,
//                           RegionList<rfx_trajx_seg*>::type *slist,
//                           const struct trajx_point *pt,
//                           const struct trajx_point *pt_next,
//                           double x[6], double dx[6] )
// {
//     double dx_p[6], x_p[6];
//     memcpy(x_p, x, sizeof(x_p));
//     memcpy(dx_p, dx, sizeof(dx_p));

//     rfx_trajx_get_dx( (struct rfx_trajx*)*slist.rbegin(), itr_next->t-t_b/2, x_p, r );

//     double t_b = pt->tb;

//     /* current point */
//     point2vector( pt, x_p, x );

//     /* next point point */
//     double x_n[6];
//     point2vector( pt, x, x_n );

//     /* compute velocity and acceleration */
//     double ddx[6];
//     for( size_t j = 0; j < 6; j ++ ) {
//         dx[j] = (x_n[j] - x[j]) / (pt_next->t - pt->t);
//         ddx[j] = (dx[j] - dx_p[j]) / t_b;
//     }

//     /* add blend about current point */
//     slist->push_back( rfx_trajx_seg_blend_rv_alloc(reg,
//                                                    pt->t - pt->tb/2,
//                                                    pt->t + pt->tb/2,
//                                                    pt->t - pt->tb/2,
//                                                    x_p, dx_p, ddx) );

// }


// void parablend_add_linear( aa_mem_region_t *reg,
//                            RegionList<rfx_trajx_seg*>::type *slist,
//                            const struct trajx_point *pt,
//                            const struct trajx_point *pt_next,
//                            double x[6], double x_n[6] )
// {
//     slist->push_back( rfx_trajx_seg_lerp_rv_alloc(reg,
//                                                   pt->t + pt->tb/2,
//                                                   pt_next->t - pt_next->tb/2,
//                                                   pt->t, x,
//                                                   pt_next->t, x_n) );


// }


// struct rfx_trajx_seg_list *
// rfx_trajx_parablend_generate( struct rfx_trajx_point_list *points, aa_mem_region_t *reg )
// {
//     RegionList<trajx_point>::type &plist = points->list;

//     if( plist.empty() ) return NULL;
//     if( 1 == plist.size() ) return NULL; /* TODO: return constant segment */

//     struct rfx_trajx_seg_list * seg_list = rfx_trajx_seg_list_alloc( reg );
//     RegionList<rfx_trajx_seg*>::type &slist = seg_list->list;

//     struct trajx_point virt_i( *plist.begin() );
//     struct trajx_point virt_f( *plist.begin() );

//     RegionList<trajx_point>::iterator itr = plist.begin();
//     RegionList<trajx_point>::iterator itr_next(itr); itr_next++;

//     double dx_p[6] = {0}, x_p[6] = {0};
//     {
//         double r[4];
//         aa_tf_duqu2qv( plist.begin()->S, r, x_p );
//         aa_tf_quat2rotvec( r, x_p+3 );
//     }

//     for(; itr_next != plist.end(); itr++, itr_next++ ) {
//         double t_b = itr->tb;
//         /* current point */
//         double x[6];
//         point2vector( &(*itr), x_p, x );

//         /* next point point */
//         double x_n[6];
//         point2vector( &(*itr_next), x, x_n );

//         /* compute velocity and acceleration */
//         double ddx[6], dx[6] = {0};
//         for( size_t j = 0; j < 6; j ++ ) {
//             dx[j] = (x_n[j] - x[j]) / (itr_next->t - itr->t);
//             ddx[j] = (dx[j] - dx_p[j]) / t_b;
//         }

//         /* add blend about current point */
//         slist.push_back( rfx_trajx_seg_blend_rv_alloc(reg, itr->t-t_b/2, itr->t+t_b/2,
//                                                       itr->t-t_b/2, x_p, dx_p, ddx) );

//         /* add linear to next point */
//         slist.push_back( rfx_trajx_seg_lerp_rv_alloc(reg, itr->t+t_b/2, itr_next->t-t_b/2,
//                                                      itr->t, x,
//                                                      itr_next->t, x_n) );
//         /* previous is end of current linear segment */
//         double r[4];
//         rfx_trajx_get_x( (struct rfx_trajx*)*slist.rbegin(), itr_next->t-t_b/2, x_p, r );
//         aa_tf_quat2rotvec_near(r, x+3, x_p+3 );

//         /* copy current to prev */
//         AA_MEM_CPY(dx_p, dx, 6);
//     }
//    /* add terminal blend */

// }

// struct rfx_trajx_seg_list *
// rfx_trajx_splend_generate( struct rfx_trajx_point_list *points, aa_mem_region_t *reg )
// {

// }
