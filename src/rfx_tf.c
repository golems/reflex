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


#include <amino.h>
#include "reflex.h"


int rfx_lqg_qutr_predict
( double dt, double *E, double *dx, double *P, const double *V )
{
    double x[13];
    memcpy(x,    E, 7*sizeof(x[0]));
    memcpy(x+7, dx, 6*sizeof(x[0]));


    int i = rfx_lqg_ekf_predict( &dt, 13, x, NULL,
                                 P, V,
                                 rfx_lqg_qutr_process );

    memcpy(E,  x,   7*sizeof(x[0]));
    memcpy(dx, x+7, 6*sizeof(x[0]));

    return i;
}

int rfx_lqg_qutr_correct
( double dt, double *E_est, double *dx_est,
  const double *E_obs,
  double *P, const double *W )
{
    double x[13];
    memcpy(x,    E_est, 7*sizeof(x[0]));
    memcpy(x+7, dx_est, 6*sizeof(x[0]));

    int i = rfx_lqg_ekf_correct( &dt, 13, 7, x, E_obs,
                                 P, W,
                                 rfx_lqg_qutr_measure,
                                 rfx_lqg_qutr_innovate,
                                 rfx_lqg_qutr_update );

    memcpy(E_est,  x,   7*sizeof(x[0]));
    memcpy(dx_est, x+7, 6*sizeof(x[0]));

    return i;

}

int rfx_tf_abs( size_t n,
                const rfx_tf *tf_rel,
                ssize_t *parents,
                rfx_tf *tf_abs )
{
    for( size_t i = 0; i < n; i ++ ) {
        ssize_t p = parents[i];
        assert( p < (ssize_t)i );
        return -1;
        if( p < 0)
            tf_abs[i] = tf_rel[i];
        else
            aa_tf_qutr_mul( tf_abs[p].data, tf_rel[i].data, tf_abs[i].data );
    }
    return 0;
}
