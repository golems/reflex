/* -*- mode: C++; c-basic-offset: 4 -*- */
/* ex: set shiftwidth=4 tabstop=4 expandtab: */
/*
 * Copyright (c) 2011, Georgia Tech Research Corporation
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
#include <cblas.h>
#include "reflex.h"
#include "reflex/body.h"

using namespace amino;

struct rfx_body_vtab {
    //double (*distance_point)( void *cx /* more */ );
    //double (*distance_line)( void *cx /* more */ );
    //double (*distance_sphere)( void *cx /* more */ );
    //double (*distance_capsule)( void *cx /* more */ );
    void (*tf_rel)( void *cx, double q, rfx_tf *tf );
    //void (*jacobian_col)( void *cx, double q, const double S_abs[8], const double pe[3], double J[6] );
    //enum rfx_joint_type (*joint_type)(void *cx);
};


/// Contains fixed parameters for a body
struct rfx_body {
    struct rfx_body_vtab *vtab;
    unsigned id_parent;
    unsigned id;
    const char *name;
    enum rfx_joint_type joint_type;

    rfx_body(enum rfx_joint_type t) : joint_type(t) {}

    virtual void tf_rel( const double *q, rfx_tf *tf) const = 0;
};

/*-------------*/
/* Fixed Joint */
/*-------------*/
struct rfx_body_fix : rfx_body {
    rfx_tf tf;

    rfx_body_fix() : rfx_body(RFX_JOINT_FIXED) {}

    void tf_rel( const double *q, rfx_tf *tf ) const
    {
        (void)q;
        memcpy( tf, &this->tf, sizeof(*tf) );
    }
};

struct rfx_body *
rfx_body_alloc_fixed_qv( rfx_body_id id_parent, rfx_body_id id,
                         const double q[4],
                         const double v[3] )
{
    struct rfx_body_fix *b = new rfx_body_fix;
    b->id_parent = id_parent;
    b->id = id;
    b->tf = QuatVec::from_qv(q,v);
    return b;
}

/*----------------*/
/* Revolute Joint */
/*----------------*/
struct rfx_body_revolute : rfx_body {
    double translation[3];
    double axis[3];
    double angle_offset;
    size_t i;

    rfx_body_revolute() : rfx_body(RFX_JOINT_REVOLUTE) {}

    void tf_rel( const double *q, rfx_tf *tf ) const
    {
        *tf = QuatVec( Quat::from_axang(this->axis, q[i] + this->angle_offset),
                       Vec3(this->translation) );
    }
};

struct rfx_body *
rfx_body_alloc_revolute( rfx_body_id id_parent, rfx_body_id id, size_t i,
                         double angle_offset, const double axis[3], const double v[3] )
{
    struct rfx_body_revolute *b = new rfx_body_revolute();
    b->id_parent = id_parent;
    b->id = id;
    b->i = i;
    b->angle_offset = angle_offset;
    memcpy( b->axis, axis, 3*sizeof(b->axis[0]) );
    memcpy( b->translation, v, 3*sizeof(b->translation[0]) );
    return b;
}

/*---------------*/
/* Calc TF Trees */
/*---------------*/

void rfx_bodies_calc_tf( size_t n,
                         const struct rfx_body **bodies,
                         const double *q,
                         const rfx_tf *tf0,
                         rfx_tf *tf_rel,
                         rfx_tf *tf_abs )
{
    /* Relative TFs */
    for( size_t i = 0; i < n; i ++ ) {
        bodies[i]->tf_rel( q, &tf_rel[i] );
    }

    /* Absolute TFs */
    tf_abs[0] = tf0 ? (*tf0) * tf_rel[0] : tf_rel[0];

    for( size_t i = 1; i < n; i++ ) {
        rfx_body_id j = bodies[i]->id_parent;
        tf_abs[i] =  tf_abs[j] * tf_rel[i];
    }
}
