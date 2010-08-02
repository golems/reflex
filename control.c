/* -*- mode: C; c-basic-offset: 4  -*- */
/* ex: set shiftwidth=4 expandtab: */
/*
 * Copyright (c) 2010, Georgia Tech Research Corporation
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *     * Redistributions of source code must retain the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer.
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials
 *       provided with the distribution.
 *     * Neither the name of the Georgia Tech Research Corporation nor
 *       the names of its contributors may be used to endorse or
 *       promote products derived from this software without specific
 *       prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY GEORGIA TECH RESEARCH CORPORATION ''AS
 * IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL GEORGIA
 * TECH RESEARCH CORPORATION BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

#include <amino.h>
#include "control.h"

/*
 * u = J^* * (  dx_r - k_p * (x - x_r) -  k_f * (F - F_r) )
 */
void ctrl_ws_lin_vfwd( const ctrl_ws_t *ws, const ctrl_ws_lin_k_t *k, double *u ) {
    double dx_u[6], x_e[6] ;
    double *J_star = (double*)AA_ALLOCAL( ws->n_q * 6 );

    // find position error
    aa_la_vsub( 3, ws->x, ws->x_r, x_e );

    // find orientation error
    {
        double r_inv[4], r_e[4], axang[4];
        aa_tf_qinv( ws->r, r_inv );        // r_inv = ws->r^{-1}
        aa_tf_qmul( r_inv, ws->r_r, r_e ); // r_e = ws->r * r_inv
        aa_tf_quat2axang( r_e, axang );    // axis-angle conversion
        for( size_t i = 0; i < 3; i ++ )
            x_e[3+i] = axang[i]*axang[3];  // x_e = axis*theta
    }

    // find workspace velocity
    // dx_u = dx_r - k_p * x_e -  k_f * (F - F_r)
    for(size_t i = 0; i < 6; i ++ ) {
        dx_u[i] = ws->dx_r[i]
            - k->p[i] * x_e[i]
            - k->f[i] * (ws->F[i] - ws->F_r[i]);
    }

    // find jointspace velocity
    aa_la_dls( 6, ws->n_q, k->dls, ws->J, J_star );
    aa_la_mvmul( ws->n_q, 6, J_star, dx_u, u );

    aa_frlocal( J_star, ws->n_q * 6);
}
