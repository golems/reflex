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

// TODO: remove deprecated code

/*--- Cartesian Space Trajectories ---*/

/* TODO: Separate structure types for trajectory via point container
 * and segment list.  Then, can use a generic segment list struct for
 * any trajectory type.  Also gives some possibility to parallelize
 * generation and tracking, or even pass generated trajectories via
 * IPC (would need to fixup vtable links for that somehow).
 */

struct rfx_trajx;

struct rfx_trajx_vtab {
    int (*generate)(struct rfx_trajx *cx);
    void (*add)(struct rfx_trajx *cx, double t, double x[3], double r[4]);
    int (*get_x)(struct rfx_trajx *cx, double t, double x[3], double r[4]);
    int (*get_dx)(struct rfx_trajx *cx, double t, double dx[6]);
    int (*get_ddx)(struct rfx_trajx *cx, double t, double ddx[6]);
};

struct rfx_trajx_point {
    double t;
    double x[3];
    double r[4];
};

typedef struct rfx_trajx {
    struct rfx_trajx_vtab *vtab;
    struct aa_mem_region *reg;
    struct rfx_trajq *trajq;

    struct rfx_trajx_point *pt_i;
    struct rfx_trajx_point *pt_f;

    size_t n_p;

    /* Poses as vector and quaternion */
    aa_mem_rlist_t *point;
} rfx_trajx_t;

static inline void rfx_trajx_add( struct rfx_trajx *cx, double t, double x[3], double r[4]) {
    cx->vtab->add( cx, t, x, r );
}


static inline void rfx_trajx_add_duqu( struct rfx_trajx *cx, double t, const double S[8] ) {
    double r[4], x[3];
    aa_tf_duqu2qv(S, r, x );
    cx->vtab->add( cx, t, x, r );
}

static inline int rfx_trajx_get_x( struct rfx_trajx *cx, double t, double x[3], double r[4]) {
    return cx->vtab->get_x( cx, t, x, r );
}

static inline int rfx_trajx_get_x_duqu( struct rfx_trajx *cx, double t, double S[8] ) {
    double x[3], r[4];
    int i = cx->vtab->get_x( cx, t, x, r );
    aa_tf_qv2duqu( r, x, S );
    return i;
}

static inline int rfx_trajx_get_dx( struct rfx_trajx *cx, double t, double dx[6]) {
    return cx->vtab->get_dx( cx, t, dx );
}

static inline int rfx_trajx_get_ddx( struct rfx_trajx *cx, double t, double ddx[6]) {
    return cx->vtab->get_ddx( cx, t, ddx );
}


int rfx_trajx_set_ctrl( struct rfx_trajx *cx, double t, rfx_ctrlx_lin_t *ctrlx );

static inline int rfx_trajx_generate( struct rfx_trajx *cx ) {
    return cx->vtab->generate( cx );
}


void rfx_trajx_destroy( struct rfx_trajx *cx );

/*-- Rotation Vector orientations --*/
/* Initialize cartesian trajectory generator
 *
 * @pre: trajq has been initialized with length of 6
 */
void rfx_trajx_rv_init( struct rfx_trajx *cx, aa_mem_region_t *reg );



void rfx_trajx_slerp_init( struct rfx_trajx *cx, aa_mem_region_t *reg );


/*--- Cartesian Segments ---*/

struct rfx_trajx_seg {
    struct rfx_trajx_vtab *vtab;
    double t_i, t_f;
};

typedef struct rfx_trajx_seg_lerp {
    struct rfx_trajx_seg seg;
    double x_i[3], x_f[3];
    double r_i[4], r_f[4];
    double tau_i, tau_f;
    double dt;
} rfx_trajx_seg_lerp_t;

/* t_i and t_f are the times when we actually start (used to pick the segment to execute)
 * tau_f and tau_f are the times used to compute the interpolation
 */

struct rfx_trajx_seg *
rfx_trajx_seg_lerp_slerp_alloc( aa_mem_region_t *reg, double t_i, double t_f,
                                double tau_i, double x_i[3], double r_i[4],
                                double tau_f, double x_f[3], double r_f[4] );

typedef struct rfx_trajx_via {
    struct rfx_trajx trajx;
    struct rfx_trajx_seg **seg;
    size_t n_seg;
} rfx_trajx_via_t;


void rfx_trajx_via_init( struct rfx_trajx_via *cx, aa_mem_region_t *reg );

typedef struct rfx_trajx_plot_opts {
    int to_file;                ///< write plot to file
    rfx_ctrlx_fun ctrlx_fun;    ///< plot joint values with this controller
    void *ctrlx_cx;
    size_t n_q;
    double *q_0;                ///< initial joint position
} rfx_trajx_plot_opts_t;

void rfx_trajx_plot( struct rfx_trajx *cx, double dt, const struct rfx_trajx_plot_opts *opts );


/*--- Cartesian Parabolic Blends ---*/

typedef struct rfx_trajx_seg_lerp_rv {
    struct rfx_trajx_seg seg;
    double x_i[6], x_f[6];
    double tau_i, tau_f;
    double dt;
} rfx_trajx_seg_lerp_rv_t;

struct rfx_trajx_seg *
rfx_trajx_seg_lerp_rv_alloc( aa_mem_region_t *reg, double t_i, double t_f,
                             double tau_i, double x_i[6],
                             double tau_f, double x_f[6] );

typedef struct rfx_trajx_seg_blend_rv {
    struct rfx_trajx_seg seg;
    double x_i[6], dx_i[6], ddx[6];
    double tau_i;
} rfx_trajx_seg_blend_rv_t;
struct rfx_trajx_seg *
rfx_trajx_seg_blend_rv_alloc( aa_mem_region_t *reg, double t_i, double t_f,
                              double tau_i, double x_i[6], double dx_i[6], double ddx[6] );


typedef struct rfx_trajx_parablend {
    struct rfx_trajx_via via;
    double t_b;
} rfx_trajx_parablend_t;


void rfx_trajx_parablend_init( struct rfx_trajx_parablend *cx, aa_mem_region_t *reg, double t_b );

/*--- Cartesian Sphereical Parabolic Blends ---*/

typedef struct rfx_trajx_seg_blend_q {
    struct rfx_trajx_seg seg;
    double x0[3];
    double dx0[3];
    double ddx[3];

    double r_i[4];
    double r_j[4];
    double r_k[4];

    double t_b;
    double t_i;
    double t_j;
    double t_k;

    double t_ij;

    double tau_i;

    double ddu_ij;
    double ddu_jk;
    double ddu_j;
} rfx_trajx_seg_blend_q_t;

struct rfx_trajx_seg *
rfx_trajx_seg_blend_q_alloc( aa_mem_region_t *reg, double t_0, double t_1,
                             double t_b,
                             double t_i, double x_i[3], double r_i[4],
                             double t_j, double x_j[3], double r_j[4],
                             double t_k, double x_k[3], double r_k[4] );

typedef struct rfx_trajx_parablend rfx_trajx_splend_t;

void rfx_trajx_splend_init( rfx_trajx_splend_t *cx, aa_mem_region_t *reg, double t_b );




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
rfx_trajx_point_list_add_qv( struct rfx_trajx_point_list *list,
                             double t, const double r[4], const double v[3] );
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


/** Add a blend point given as quaternion and vector to the list.
 *
 * Arguments are copied into memory allocated of the the region for
 * list.
 */
int
rfx_trajx_point_list_addb_qv( struct rfx_trajx_point_list *list,
                              double t, double t_blend, const double r[4], const double v[3] );

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
rfx_trajx_seg_list_add( struct rfx_trajx_seg_list *seg_list,
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

/** Get trajectory pose as transformation matrix */
int
rfx_trajx_seg_list_get_x_tfmat( struct rfx_trajx_seg_list *seg,
                                double t, double T[12] );
/** Get trajectory pose as transformation matrix, and get velocity */
int
rfx_trajx_seg_list_get_dx_tfmat( struct rfx_trajx_seg_list *seg,
                                 double t, double T[12], double dx[6] );
/** Get trajectory pose as transformation matrix, and get velocity and acceleration */
int
rfx_trajx_seg_list_get_ddx_tfmat( struct rfx_trajx_seg_list *seg,
                                  double t,double T[12], double dx[6], double ddx[6] );


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
rfx_trajx_seg_list_plot( struct rfx_trajx_seg_list *cx, double dt,
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
