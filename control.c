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
#include "reflex.h"


void rfx_ctrl_ws_init( rfx_ctrl_ws_t *g, size_t n ) {
    g->n_q = n;
    g->q = AA_NEW0_AR( double, n );
    g->dq = AA_NEW0_AR( double, n );
    g->J = AA_NEW0_AR( double, n*6 );
    g->q_r = AA_NEW0_AR( double, n );
    g->dq_r = AA_NEW0_AR( double, n );
    g->F_max = 0;
    g->M_max = 0;
}

void rfx_ctrl_ws_destroy( rfx_ctrl_ws_t *g ) {
    free(g->q);
    free(g->dq);
    free(g->J);
    free(g->q_r);
    free(g->dq_r);
}

AA_API void rfx_ctrl_ws_lin_k_init( rfx_ctrl_ws_lin_k_t *k, size_t n_q  ) {
    k->n_q = n_q;
    k->q = AA_NEW0_AR( double, k->n_q );
}
AA_API void rfx_ctrl_ws_lin_k_destroy( rfx_ctrl_ws_lin_k_t *k ) {
    free(k->q);
}

/*
 * u = J^* * (  dx_r - k_p * (x - x_r) -  k_f * (F - F_r) )
 */
void rfx_ctrl_ws_lin_vfwd( const rfx_ctrl_ws_t *ws, const rfx_ctrl_ws_lin_k_t *k, double *u ) {
    double dx_u[6], x_e[6] ;
    double dq_r[ws->n_q];

    assert( ws->n_q == k->n_q );

    // check force limits
    {
        double FF = aa_la_dot( 3, ws->F, ws->F );
        double MM = aa_la_dot( 3, ws->F+3, ws->F+3 );
        if( ( ws->F_max > 0 && FF > ws->F_max*ws->F_max ) ||
            ( ws->M_max > 0 && MM > ws->M_max*ws->M_max ) ) {
            aa_fzero( u, ws->n_q );
            return;
        }
    }

    // find position error
    aa_la_vsub( 3, ws->x, ws->x_r, x_e );

    // find orientation error
    {
        double r_e[4];
        aa_tf_qrel(ws->r, ws->r_r, r_e);
        aa_tf_quat2rotvec_near( r_e, AA_FAR(0,0,0), x_e+3 );    // axis-angle conversion
    }

    // find workspace velocity
    // dx_u = dx_r - k_p * x_e -  k_f * (F - F_r)
    for(size_t i = 0; i < 6; i ++ ) {
        dx_u[i] = ws->dx_r[i]
            - k->p[i] * x_e[i]
            - k->f[i] * (ws->F[i] - ws->F_r[i]);
    }
    // jointspace reference velocity
    for( size_t i = 0; i < ws->n_q; i ++ ) {
        dq_r[i] = -k->q[i] * (ws->q[i] - ws->q_r[i]);// + ws->dq_r[i];
    }
    // find jointspace velocity
    aa_la_dlsnp( 6, ws->n_q, k->dls, ws->J, dx_u, dq_r, u );

}
