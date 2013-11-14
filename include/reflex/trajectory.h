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

struct rfx_trajq_point {
    double t;       ///< time
    double q[1];    ///< points
};

/** List of via points for a trajectory */
typedef struct rfx_trajq_points {
    aa_mem_region_t *reg;   ///< memory region for allocation
    size_t n_q;             ///< number of joints
    size_t n_t;             ///< number of points
    double t_i, t_f;        ///< initial and final time
    double *q_i, *q_f;      ///< inital and final joint configuration
    aa_mem_rlist_t *point;  ///< list of points
} rfx_trajq_points_t;

rfx_trajq_points_t *
rfx_trajq_points_alloc( aa_mem_region_t *reg, size_t n_q );

void rfx_trajq_points_add( rfx_trajq_points_t *traj, double t, double *q );

/*--- Trajectory Segments ---*/

struct rfx_trajq_seg_vtab {
    int (*get_q)(void *cx, double t, double *q);
    int (*get_dq)(void *cx, double t, double *q, double *dq);
    int (*get_ddq)(void *cx, double t, double *q, double *dq, double *ddq);
};


/** Segment of a configuration space trajectory
 */
typedef struct rfx_trajq_seg {
    struct rfx_trajq_seg_vtab *vtab;
    double t_i; ///< initial time for this segment
    double t_f; ///< final time for this segment
    size_t n_q; ///< number of joints
} rfx_trajq_seg_t;



static inline int
rfx_trajq_seg_get_q( rfx_trajq_seg_t *seg, double t, double *q ) {
    return seg->vtab->get_q( seg, t, q );
}
static inline int
rfx_trajq_seg_get_dq(  rfx_trajq_seg_t *seg, double t, double *q, double *dq ) {
    return seg->vtab->get_dq( seg, t, q, dq );
}
static inline int
rfx_trajq_seg_get_ddq( rfx_trajq_seg_t *seg, double t, double *q, double *dq, double *ddq ) {
    return seg->vtab->get_ddq( seg, t, q, dq, ddq );
}


/** List of trajectory segments.
 */
typedef struct rfx_trajq_seg_list {
    struct rfx_trajq_seg_vtab *vtab;
    double t_i;             ///< initial time for this segment
    double t_f;             ///< final time for this segment
    size_t n_q;             ///< number of joints
    size_t n_t;             ///< number of segments
    aa_mem_region_t *reg;   ///< memory region for allocation
    aa_mem_rlist_t *seg;    ///< list of segments
    struct aa_mem_cons *last_seg; ///< pointer to last referenced segment
} rfx_trajq_seg_list_t;

/* Initialize struct, performing future alloctions out of reg */
struct rfx_trajq_seg_list *
rfx_trajq_seg_list_alloc( aa_mem_region_t *reg );

/* Add a segment */
void rfx_trajq_seg_list_add( rfx_trajq_seg_list_t *seglist, rfx_trajq_seg_t *seg );


static inline int
rfx_trajq_seg_list_get_q( rfx_trajq_seg_list_t *seglist, double t, double *q ) {
    return seglist->vtab->get_q( seglist, t, q );
}
static inline int
rfx_trajq_seg_list_get_dq(  rfx_trajq_seg_list_t *seglist, double t, double *q, double *dq ) {
    return seglist->vtab->get_dq( seglist, t, q, dq );
}
static inline int
rfx_trajq_seg_list_get_ddq( rfx_trajq_seg_list_t *seglist, double t, double *q, double *dq, double *ddq ) {
    return seglist->vtab->get_ddq( seglist, t, q, dq, ddq );
}

void rfx_trajq_seg_plot( struct rfx_trajq_seg_list *cx, double dt );

struct rfx_trajq_seg_param {
    rfx_trajq_seg_t parent;
    double p[0];
};


/*--- Trajectory Generators ---*/

/** Generate a parabolic blend trajectory through all points.
 *
 * Uses a single, given blend time to compute accelerations
 */
rfx_trajq_seg_list_t *
rfx_trajq_gen_pblend_tm1( aa_mem_region_t *reg, rfx_trajq_points_t *points, double t_blend );


/*--- LERP Segments (Constant Velocity) ---*/

/** Linearly interpolate between initial and final configuration
 */
typedef struct rfx_trajq_seg_param rfx_trajq_seg_dq_t;

/** Allocate segment to LERP from q_i with given velocity */
rfx_trajq_seg_dq_t *
rfx_trajq_seg_dq_alloc( aa_mem_region_t *reg, size_t n_q,
                        double t_i, double *q_i, double *dq,
                        double t_f );

/** Allocate segment to LERP between q_i and q_f */
rfx_trajq_seg_dq_t *
rfx_trajq_seg_dq_alloc2( aa_mem_region_t *reg, size_t n_q,
                         double t_i, double *q_i,
                         double t_f, double *q_f );


/** Add a constant velocity segment to end of trajectory. */
int rfx_trajq_seg_dq_link( rfx_trajq_seg_list_t *seglist, double *dq, double t_f );

/*--- Constant Acceleration Segments ---*/
typedef struct rfx_trajq_seg_param rfx_trajq_seg_2dq_t;

/** Allocate a constant acceleration segment */
rfx_trajq_seg_2dq_t *
rfx_trajq_seg_2dq_alloc( aa_mem_region_t *reg, size_t n_q,
                         double t_i, double *q_i,
                         double *dq, double *ddq,
                         double t_f );

/** Allocate a constant acceleration segment */
rfx_trajq_seg_2dq_t *
rfx_trajq_seg_2dq_alloc2( aa_mem_region_t *reg,
                          size_t n_q,
                          double t_i, double *q_i, double *dq_i,
                          double t_f, double *q_f );

/** Add a constant acceleration segment to end of trajectory. */
int rfx_trajq_seg_2dq_link( rfx_trajq_seg_list_t *seglist, double *ddq, double t_f );

/*--- Constant Jerk Segments ---*/
struct rfx_trajq_seg_3dq {
    rfx_trajq_seg_t parent;
    size_t n_q;
    double *q_i;
    double *dq_i;
    double *ddq_i;
    double *dddq;
};

/** Allocate a constant acceleration segment */
struct rfx_trajq_seg_3dq *
rfx_trajq_seg_3dq_alloc( aa_mem_region_t *reg, size_t n_q,
                         double t_i, double *q_i,
                         double *dq, double *ddq, double *dddq,
                         double t_f );


/** Allocate a constant acceleration segment */
struct rfx_trajq_seg_3dq *
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
