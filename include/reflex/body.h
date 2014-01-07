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

#ifndef REFLEX_BODY_H
#define REFLEX_BODY_H
/** \file reflex.h */


/** \mainpage Reflex: Real-Time Control
 *
 * \section feature Features
 *
 * - Robot Arms
 *   - Workspace (Cartesian) Control
 *     - Damped-Least squares Jacobian Inverse with Nullspace Projection
 *   - Trapezoidal velocity profiles in workspace
 * - Linear Quadratic Gaussian \ref rfx_lqg_t
 *   - Solutions to Algebraic Riccati Equation
 *   - Discrete Time Kalman Filter,
 *     - \ref rfx_lqg_kf_predict
 *     - \ref rfx_lqg_kf_correct
 *   - Kalman-Bucy (continuous time) Filter
 *     - \ref rfx_lqg_kbf_gain
 *     - \ref rfx_lqg_kbf_step1, \ref rfx_lqg_kbf_step4
 *
 * \section wsctrl Workspace Control HOWTO
 *
 * -# Initialization
 *    -# get a \ref rfx_ctrl_t struct
 *    -# call \ref rfx_ctrl_ws_init
 *    -# get a \ref rfx_ctrl_ws_lin_k_t struct
 *    -# call \ref rfx_ctrl_ws_lin_k_init
 *    -# Set gains the the rfx_ctrl_ws_lin_k_t
 *    -# Set limits the the rfx_ctrl_t
 * -# At each time step (probably do this at 1kHz)
 *    -# Update the position q, velocity dq, and Jacobian J in the rfx_ctrl_t
 *    -# Update the reference q_r, x_r, etc in the rfx_ctrl_t
 *    -# Compute desired velocities with \ref rfx_ctrl_ws_lin_vfwd
 *    -# Send the velocities to your arm
 * -# Finalization
 *    -# Call \ref rfx_ctrl_ws_destroy
 *    -# Call \ref rfx_ctrl_ws_lin_k_destroy
 *
 * \author Neil T. Dantam
 * \author Developed at the Georgia Tech Humanoid Robotics Lab
 * \author Under Direction of Professor Mike Stilman
 *
 * \section License
 *
 * Copyright (c) 2010-2011, Georgia Tech Research Corporation.
 * All rights reserved.
 *
 *   Redistribution and use in source and binary forms, with or
 *   without modification, are permitted provided that the following
 *   conditions are met:
 *
 *   - Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *
 *   - Redistributions in binary form must reproduce the above
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



#ifdef __cplusplus
extern "C" {
#endif //__cplusplus


/** Identifier for a body.
 *
 * Constraints:
 * - Ordering: child frames must have IDs greater than all
 *   parent frames.
 */
typedef unsigned rfx_body_id;

#define RFX_BODY_INDEX_ROOT ((rfx_body_id)-1)

enum rfx_joint_type {
    RFX_JOINT_FIXED,
    RFX_JOINT_REVOLUTE,
    RFX_JOINT_PRISMATIC
};


/** Opaque struct for a body (rigid member) */
struct rfx_body;

/** Allocate a body at a fixed transform from parent */
struct rfx_body *
rfx_body_alloc_fixed_qv( rfx_body_id id_parent, rfx_body_id id,
                         const double q[4],
                         const double v[3] );

/** Allocate a body with rotating joint from parent */
struct rfx_body *
rfx_body_alloc_revolute( rfx_body_id id_parent, rfx_body_id id,
                         size_t i, double angle_offset,
                         const double axis[3], const double v[3] );

/* Lists of bodies.
 * Must be sorted by ID.
 */

/** Calculate all relative and absolute transforms */
int rfx_bodies_calc_tf( size_t n,
                        const struct rfx_body **bodies,
                        const double *q,
                        const rfx_tf *tf0,
                        rfx_tf *tf_rel,
                        rfx_tf *tf_abs );


/* Clone a contiguous block of bodies.
 *
 *  \pre body[old_id0] must be the root of all cloned bodies
 */
int rfx_bodies_clone( size_t n,
                      struct rfx_body **bodies,
                      rfx_body_id old_id0, rfx_body_id old_id1, size_t old_i,
                      rfx_body_id new_parent,
                      rfx_body_id new_id0, rfx_body_id new_id1,  size_t new_i );


/** Compute jacobian for a subset of the bodies */
int rfx_bodies_jacobian( size_t n,
                         const struct rfx_body **bodies,
                         const rfx_tf *tf_rel,
                         const rfx_tf *tf_abs,
                         size_t n_indices, const size_t *indices,
                         double *J, size_t ldJ );

#ifdef __cplusplus
}
#endif //__cplusplus

#endif //REFLEX_H
