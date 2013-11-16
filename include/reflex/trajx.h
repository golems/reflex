/* -*- mode: C; c-basic-offset: 4 -*- */
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

#ifndef REFLEX_TRAJX_H
#define REFLEX_TRAJX_H

#ifdef __cplusplus
extern "C" {
#endif //__cplusplus


/*--- Base Segment Type ---*/
/* struct rfx_trajx_seg; */
/* int */
/* rfx_trajx_seg_get_x_duqu( struct rfx_trajx_seg *seg, double S[8] ); */
/* int */
/* rfx_trajx_seg_get_dx_duqu( struct rfx_trajx_seg *seg, double S[8], double dx[6] ); */
/* int */
/* rfx_trajx_seg_get_ddx_duqu( struct rfx_trajx_seg *seg, double S[8], double dx[6], double ddx[6] ); */

/*--- Point Lists ---*/

/** Opaque type for lists of points */
struct rfx_trajx_point_list;

/** Allocate a point list from region
 *
 * All points added to list will be allocated out of region as well
 */
struct rfx_trajx_point_list *
rfx_trajx_point_list_alloc( aa_mem_region_t *region );

/** Add a point given as vector and quaternion to the list.
 *
 * Arguments are copied into memory allocated of the the region for
 * list.
 */
int
rfx_trajx_point_list_add_xr( struct rfx_trajx_point_list *list,
                             double t, const double x[3], const double r[4] );
/** Add a point given as dual quaternion to the list.
 *
 * Arguments are copied into memory allocated of the the region for
 * list.
 */
int
rfx_trajx_point_list_add_duqu( struct rfx_trajx_point_list *list,
                               double t, const double S[8] );

/** Add a point given as transformation matrix to the list.
 *
 * Arguments are copied into memory allocated of the the region for
 * list.
 *
 * T is column major, and the fourth column of the matrix is ommitted.
 */
int
rfx_trajx_point_list_add_tfmat( struct rfx_trajx_point_list *list,
                             double t, const double T[12] );


/** Add a blend point given as vector and quaternion to the list.
 *
 * Arguments are copied into memory allocated of the the region for
 * list.
 */
int
rfx_trajx_point_list_addb_xr( struct rfx_trajx_point_list *list,
                              double t, double t_blend, const double x[3], const double r[4] );
/** Add a blend point given as dual quaternion to the list.
 *
 * Arguments are copied into memory allocated of the the region for
 * list.
 */
int
rfx_trajx_point_list_addb_duqu( struct rfx_trajx_point_list *list,
                                double t, double t_blend, const double S[8] );

/** Add a blend point given as transformation matrix to the list.
 *
 * Arguments are copied into memory allocated of the the region for
 * list.
 *
 * T is column major, and the fourth column of the matrix is ommitted.
 */
int
rfx_trajx_point_list_addb_tfmat( struct rfx_trajx_point_list *list,
                                 double t, double t_blend, const double T[12] );


/*--- Segment Lists ---*/

/** Opaque type for list of segments */
struct rfx_trajx_seg_list;

/** Allocate a segment list from region
 *
 * All points added to list will be allocated out of region as well
 */
struct rfx_trajx_seg_list *
rfx_trajx_seg_list_alloc( aa_mem_region_t *region );

/** Add a segment to the list
 *
 * List storage allocated from previously region for the list.
 */
int
rfx_trajx_seg_list_add( struct rfx_trajx_seg_list seg_list,
                        struct rfx_trajx_seg *seg );

/** Get trajectory pose as dual quaternion */
int
rfx_trajx_seg_list_get_x_duqu( struct rfx_trajx_seg_list *seg,
                               double t, double S[8] );
/** Get trajectory pose as dual quaternion and velocity */
int
rfx_trajx_seg_list_get_dx_duqu( struct rfx_trajx_seg_list *seg,
                                double t, double S[8], double dx[6] );

/** Get trajectory pose as dual quaternion and velocity and acceleration */
int
rfx_trajx_seg_list_get_ddx_duqu( struct rfx_trajx_seg_list *seg,
                                 double t, double S[8], double dx[6], double ddx[6] );


/** Get trajectory pose as quaternion and vector */
int
rfx_trajx_seg_list_get_x_qv( struct rfx_trajx_seg_list *seg,
                             double t, double r[4], double x[3] );
/** Get trajectory pose as quaternion and vector, and get velocity */
int
rfx_trajx_seg_list_get_dx_qv( struct rfx_trajx_seg_list *seg,
                              double t, double r[4], double x[3], double dx[6] );
/** Get trajectory pose as quaternion and vector, and get velocity and acceleration */
int
rfx_trajx_seg_list_get_ddx_qv( struct rfx_trajx_seg_list *seg,
                               double t, double r[4], double x[3], double dx[6], double ddx[6] );

/** Get trajectory pose [translation log(quaternion)] */
int
rfx_trajx_seg_list_get_x_vl( struct rfx_trajx_seg_list *seg,
                             double t, double x[6] );

/** Get trajectory initial time */
double
rfx_trajx_seg_list_get_t_i( struct rfx_trajx_seg_list *seg );

/** Get trajectory final time */
double
rfx_trajx_seg_list_get_t_f( struct rfx_trajx_seg_list *seg );


/** Plot trajectory */
void
rfx_trajx_seglist_plot( struct rfx_trajx_seg_list *cx, double dt,
                        const struct rfx_trajx_plot_opts *xopts ) ;

/*--- Segment Types ---*/

/** Generate a parabolic blend trajectory
 *
 * The trajectory segment list and segments will be allocated out of
 * region.  This may be the same or different from the region used for
 * points.
 */
struct rfx_trajx_seg_list *
rfx_trajx_parablend_generate( struct rfx_trajx_point_list *points, aa_mem_region_t *region );

/** Generate a parabolic blend trajectory
 *
 * The trajectory segment list and segments will be allocated out of
 * region.  This may be the same or different from the region used for
 * points.
 */
struct rfx_trajx_seg_list *
rfx_trajx_splend_generate( struct rfx_trajx_point_list *points, aa_mem_region_t *region );


#ifdef __cplusplus
}
#endif //__cplusplus

#endif //REFLEX_TRAJX_H
