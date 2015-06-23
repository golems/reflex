/* -*- mode: C; c-basic-offset: 4 -*- */
/* ex: set shiftwidth=4 tabstop=4 expandtab: */
/*
 * Copyright (c) 2010-2011, Georgia Tech Research Corporation
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

#ifndef REFLEX_H
#define REFLEX_H
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


// Symbol:
/*
 * q: configuration
 * x: state/workspace position
 * u: input/command
 * z: sensor
 * dq: joint velocity
 * dq: workspace velocity
 * F: generalized force
 * R: rotation matrix
 * p: quaternion
 * k: gain
*/

#ifdef __cplusplus
extern "C" {
#endif //__cplusplus

typedef struct aa_tf_qv rfx_tf;
typedef struct aa_tf_qv_dx rfx_tf_dx;


// FIXME: controller should return these codes
typedef enum {
    RFX_OK = 0,
    RFX_INVAL,
    RFX_LIMIT_POSITION,
    RFX_LIMIT_POSITION_ERROR,
    RFX_LIMIT_FORCE,
    RFX_LIMIT_MOMENT,
    RFX_LIMIT_FORCE_ERROR,
    RFX_LIMIT_MOMENT_ERROR,
    RFX_LIMIT_CONFIGURATION,
    RFX_LIMIT_CONFIGURATION_ERROR
} rfx_status_t;

AA_API const char* rfx_status_string(rfx_status_t i);

#define RFX_PERROR(s, r) {                                              \
        rfx_status_t _rfx_$_perror_r = r;                               \
        const char *_rfx_$_perror_s = s;                                \
        fprintf(stderr, "%s: %s (%d)\n",                                \
                _rfx_$_perror_s ? _rfx_$_perror_s : "",                 \
                rfx_status_string(_rfx_$_perror_r), _rfx_$_perror_r);   \
    }

#define RFX_SIGN(x) ((x > 0) ? 1 : ((x < 0) ? -1 : 0))

#include "reflex/control.h"
#include "reflex/lqg.h"
#include "reflex/kinematics.h"
#include "reflex/trajq.h"
#include "reflex/trajx.h"
#include "reflex/body.h"
#include "reflex/tf.h"

#ifdef __cplusplus
}
#endif //__cplusplus

#endif //REFLEX_H
