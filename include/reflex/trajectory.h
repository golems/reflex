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

#ifndef REFLEX_TRAJECTORY_H
#define REFLEX_TRAJECTORY_H

#ifdef __cplusplus
extern "C" {
#endif //__cplusplus

/*--- Point List ---*/

struct rfx_trajq_point ;
struct rfx_trajq_points;

struct rfx_trajq_points *
rfx_trajq_points_alloc( aa_mem_region_t *reg, size_t n_q );

void rfx_trajq_points_add( struct rfx_trajq_points *traj, double t, double *q );

/*--- Trajectory Segments ---*/
struct rfx_trajq_seg;

int
rfx_trajq_seg_get_q( struct rfx_trajq_seg *seg, double t, double *q );
int
rfx_trajq_seg_get_dq( struct rfx_trajq_seg *seg, double t, double *q, double *dq );
int
rfx_trajq_seg_get_ddq( struct rfx_trajq_seg *seg, double t, double *q, double *dq, double *ddq );


struct rfx_trajq_seg_list;

/* Initialize struct, performing future alloctions out of reg */
struct rfx_trajq_seg_list *
rfx_trajq_seg_list_alloc( aa_mem_region_t *reg );

/* Add a segment */
void rfx_trajq_seg_list_add( struct rfx_trajq_seg_list *seglist, struct rfx_trajq_seg*seg );


int
rfx_trajq_seg_list_get_q( struct rfx_trajq_seg_list *seglist, double t, double *q ) ;
int
rfx_trajq_seg_list_get_dq(  struct rfx_trajq_seg_list *seglist, double t, double *q, double *dq ) ;
int
rfx_trajq_seg_list_get_ddq( struct rfx_trajq_seg_list *seglist, double t, double *q, double *dq, double *ddq ) ;

double
rfx_trajq_seg_list_get_t_i( struct rfx_trajq_seg_list *seglist );
double
rfx_trajq_seg_list_get_t_f( struct rfx_trajq_seg_list *seglist );
size_t
rfx_trajq_seg_list_get_n_q( struct rfx_trajq_seg_list *seglist );

void rfx_trajq_seg_plot( struct rfx_trajq_seg_list *cx, double dt );



/*--- Trajectory Generators ---*/

/** Generate a parabolic blend trajectory through all points.
 *
 * Uses a single, given blend time to compute accelerations
 */
struct rfx_trajq_seg_list *
rfx_trajq_gen_pblend_tm1( aa_mem_region_t *reg, struct rfx_trajq_points *points, double t_blend );


/*--- LERP Segments (Constant Velocity) ---*/

/** Allocate segment to LERP from q_i with given velocity */
struct rfx_trajq_seg *
rfx_trajq_seg_dq_alloc( aa_mem_region_t *reg, size_t n_q,
                        double t_i, double *q_i, double *dq,
                        double t_f );

/** Allocate segment to LERP between q_i and q_f */
struct rfx_trajq_seg *
rfx_trajq_seg_dq_alloc2( aa_mem_region_t *reg, size_t n_q,
                         double t_i, double *q_i,
                         double t_f, double *q_f );


/** Add a constant velocity segment to end of trajectory. */
int rfx_trajq_seg_dq_link( struct rfx_trajq_seg_list *seglist, double *dq, double t_f );

/*--- Constant Acceleration Segments ---*/
typedef struct rfx_trajq_seg_param rfx_trajq_seg_2dq_t;

/** Allocate a constant acceleration segment */
struct rfx_trajq_seg *
rfx_trajq_seg_2dq_alloc( aa_mem_region_t *reg, size_t n_q,
                         double t_i, double *q_i,
                         double *dq, double *ddq,
                         double t_f );

/** Allocate a constant acceleration segment */
struct rfx_trajq_seg *
rfx_trajq_seg_2dq_alloc2( aa_mem_region_t *reg,
                          size_t n_q,
                          double t_i, double *q_i, double *dq_i,
                          double t_f, double *q_f );

/** Add a constant acceleration segment to end of trajectory. */
int rfx_trajq_seg_2dq_link( struct rfx_trajq_seg_list *seglist, double *ddq, double t_f );

/*--- Constant Jerk Segments ---*/

/** Allocate a constant acceleration segment */
struct rfx_trajq_seg *
rfx_trajq_seg_3dq_alloc( aa_mem_region_t *reg, size_t n_q,
                         double t_i, double *q_i,
                         double *dq, double *ddq, double *dddq,
                         double t_f );


/** Allocate a constant acceleration segment */
struct rfx_trajq_seg *
rfx_trajq_seg_3dq_alloc2( aa_mem_region_t *reg, size_t n_q,
                          double t_i, double *q_i,
                          double *dq, double *ddq,
                          double t_f, double *dq_f, double *ddq_f );



/*--- Dense waypoint segments ---*/

/** Allocate segment for dense waypoints.
 *
 * Does NOT create a copy of Q, but rather stores reference to Q.
 */
struct rfx_trajq_seg*
rfx_trajq_seg_dense_alloc( aa_mem_region_t *reg, size_t n_q, size_t n_p,
                           double t_i, double dt,
                           double *Q );

#ifdef __cplusplus
}
#endif //__cplusplus

#endif //REFLEX_TRAJECTORY_H
