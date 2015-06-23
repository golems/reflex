/* -*- mode: C; c-basic-offset: 4 -*- */
/* ex: set shiftwidth=4 tabstop=4 expandtab: */
/*
 * Copyright (c) 2010, Georgia Tech Research Corporation
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

const char* rfx_status_string(rfx_status_t i) {
    switch(i) {
    case RFX_OK: return "OK";
    case RFX_INVAL: return "INVAL";
    case RFX_LIMIT_POSITION: return "LIMIT_POSITION";
    case RFX_LIMIT_POSITION_ERROR: return "LIMIT_POSITION_ERROR";
    case RFX_LIMIT_FORCE: return "LIMIT_FORCE";
    case RFX_LIMIT_MOMENT: return "LIMIT_MOMENT";
    case RFX_LIMIT_FORCE_ERROR: return "LIMIT_FORCE_ERROR";
    case RFX_LIMIT_MOMENT_ERROR: return "LIMIT_MOMENT_ERROR";
    case RFX_LIMIT_CONFIGURATION: return "LIMIT_CONFIGURATION";
    case RFX_LIMIT_CONFIGURATION_ERROR: return "LIMIT_CONFIGURATION_ERROR";
    }
    return "unknown";
}


void rfx_ctrlx_state_init( struct rfx_ctrlx_state *x, size_t n ) {
    x->q  = AA_NEW0_AR( double, n );
    x->dq = AA_NEW0_AR( double, n );
    x->S  = AA_NEW0_AR( double, 8 );
    x->dx = AA_NEW0_AR( double, 6 );
    x->F  = AA_NEW0_AR( double, 6 );
    AA_MEM_CPY( x->S, aa_tf_duqu_ident, 8 );
}


void rfx_ctrlx_state_init_region( struct rfx_ctrlx_state *x, aa_mem_region_t *reg, size_t n ) {
    x->q  = AA_MEM_REGION_NEW_N( reg, double, n );
    AA_MEM_SET( x->q, 0, n );
    x->dq = AA_MEM_REGION_NEW_N( reg, double, n );
    AA_MEM_SET( x->dq, 0, n );
    x->S  = AA_MEM_REGION_NEW_N( reg, double, 8 );
    AA_MEM_CPY( x->S, aa_tf_duqu_ident, 8 );
    x->dx = AA_MEM_REGION_NEW_N( reg, double, 6 );
    AA_MEM_SET( x->dx, 0, 6 );
    x->F  = AA_MEM_REGION_NEW_N( reg, double, 6 );
    AA_MEM_SET( x->F, 0, 6 );
}

void rfx_ctrl_ws_init( rfx_ctrl_ws_t *g, size_t n ) {
    memset( g, 0, sizeof(*g) );
    g->n_q = n;

    rfx_ctrlx_state_init( &g->ref, n );
    rfx_ctrlx_state_init( &g->act, n );

    g->J = AA_NEW0_AR( double, n*6 );

    g->q_min = AA_NEW0_AR( double, n );
    g->q_max = AA_NEW0_AR( double, n );
}

void rfx_ctrlx_state_destroy( struct rfx_ctrlx_state *x ) {
    free(x->q);
    free(x->dq);
    free(x->S);
    free(x->dx);
    free(x->F);
}

void rfx_ctrl_ws_destroy( rfx_ctrl_ws_t *g ) {
    rfx_ctrlx_state_destroy(&g->ref);
    rfx_ctrlx_state_destroy(&g->act);

    free(g->J);
    free(g->q_min);
    free(g->q_max);
}

AA_API void rfx_ctrl_ws_lin_k_init( rfx_ctrl_ws_lin_k_t *k, size_t n_q  ) {
    k->n_q = n_q;
    k->q = AA_NEW0_AR( double, k->n_q );
}
AA_API void rfx_ctrl_ws_lin_k_destroy( rfx_ctrl_ws_lin_k_t *k ) {
    free(k->q);
}

// FIXME: check directions for all limits
// FIXME: add hard limits
static rfx_status_t check_limit( const rfx_ctrl_t *g, const double dx[3] ) {
    // F_max
    if( (g->F_max > 0) /* has limit */ &&
        (aa_la_dot( 3, g->act.F, g->act.F ) > g->F_max*g->F_max) /* magnitude check */ &&
        0 < aa_la_dot( 3, g->act.F, dx ) /*direction check*/ )
    {
        return RFX_LIMIT_FORCE;
    }
    // M_max
    if( (g->M_max > 0) &&
        (aa_la_dot( 3, g->act.F+3, g->act.F+3 ) > g->M_max*g->M_max) )
        return RFX_LIMIT_MOMENT;
    // q_min, q_max
    for( size_t i = 0; i < g->n_q; i++ ) {
        if( (g->act.q[i] < g->q_min[i]) || (g->act.q[i] > g->q_max[i] ) )
            return RFX_LIMIT_CONFIGURATION;
    }
    // x_min, x_max
    double x[3], x_r[3];
    aa_tf_duqu_trans( g->act.S, x );
    aa_tf_duqu_trans( g->ref.S, x_r );
    for( size_t i = 0; i < 3; i++ ) {
        if( (x[i] < g->x_min[i]) || (x[i] > g->x_max[i] ) )
            return RFX_LIMIT_POSITION;
    }
    // e_q_max
    if( (g->e_q_max > 0) &&
        ( aa_la_ssd(g->n_q, g->act.q, g->ref.q) > g->e_q_max * g->e_q_max ) )
        return RFX_LIMIT_CONFIGURATION_ERROR;
    // e_x_max
    if( (g->e_x_max > 0) &&
        ( aa_la_ssd(3, x, x_r) > g->e_x_max * g->e_x_max ) )
        return RFX_LIMIT_POSITION_ERROR;
    // e_F_max
    if( (g->e_F_max > 0) &&
        ( aa_la_ssd(3, g->act.F, g->ref.F) > g->e_F_max * g->e_F_max ) )
        return RFX_LIMIT_FORCE_ERROR;;
    // e_M_max
    if( (g->e_M_max > 0) &&
        ( aa_la_ssd(3, g->act.F+3, g->ref.F+3) > g->e_M_max * g->e_M_max ) )
        return RFX_LIMIT_MOMENT_ERROR;

    return RFX_OK;
}

/*
 * u = J^* * (  dx_r - k_p * (x - x_r) -  k_f * (F - F_r) )
 */
rfx_status_t rfx_ctrl_ws_lin_vfwd( const rfx_ctrl_ws_t *ws, const rfx_ctrl_ws_lin_k_t *k, double *u ) {
    double dx_u[6], x_e[6];
    double dq_r[ws->n_q];

    assert( ws->n_q == k->n_q );

    AA_MEM_ZERO(u, ws->n_q );

    // find position error
    /* aa_la_vsub( 3, ws->x, ws->x_r, x_e ); */

    /* // find orientation error */
    /* { */
    /*     double r_e[4], rln[4]; */
    /*     aa_tf_qmulc(ws->r, ws->r_r, r_e); */
    /*     // this is really the quaternion logarithm! */
    /*     aa_tf_qminimize(r_e); */
    /*     aa_tf_quat2rotvec( r_e, x_e+3 );    // axis-angle conversion */
    /* } */

    // relative dual quaternion -> twist -> velocity
    {
        double twist[8], de[8];
        aa_tf_duqu_mulc( ws->act.S, ws->ref.S, de );  // de = d*conj(d_r)
        aa_tf_duqu_minimize(de);
        aa_tf_duqu_ln( de, twist );     // twist = log( de )
        aa_tf_duqu_twist2vel( ws->act.S, twist, x_e );
    }

    //printf("xe: "); aa_dump_vec(stdout, x_e, 6);
    //printf("xed: "); aa_dump_vec(stdout, x_e_d, 6);
    //if( aa_la_ssd(6, x_e_s, x_e_d) > 1e-1 ) {
    //    abort();
    //}

    // find workspace velocity
    // dx_u = dx_r - k_p * x_e -  k_f * (F - F_r)
    for(size_t i = 0; i < 6; i ++ ) {
        dx_u[i] = ws->ref.dx[i]
            - k->p[i] * x_e[i]
            - k->f[i] * (ws->act.F[i] - ws->ref.F[i]);
    }


    // check force limits
    {
        rfx_status_t r = check_limit( ws, dx_u );
        if( RFX_OK != r ) {
            AA_MEM_ZERO( u, ws->n_q );
            return r;
        }
    }


    // jointspace reference velocity
    for( size_t i = 0; i < ws->n_q; i ++ ) {
        dq_r[i] = -k->q[i] * (ws->act.q[i] - ws->ref.q[i]);// + ws->dq_r[i];
    }

    // find damped inverse
    double J_star[6*ws->n_q];  // q is probally small, assume this fits on the stack

    // Compute a damped pseudo inverse
    if( k->s2min > 0 ) {
        aa_la_dzdpinv( 6, ws->n_q, k->s2min, ws->J, J_star );
    } else  {
        aa_la_dpinv( 6, ws->n_q, k->dls, ws->J, J_star );
    }

    // damped least squares with null-space projection
    aa_la_xlsnp( 6, ws->n_q, ws->J, J_star, dx_u, dq_r, u );

    /* aa_la_dlsnp( 6, ws->n_q, k->dls, ws->J, dx_u, dq_r, u ); */

    return RFX_OK;
}

void rfx_ctrl_ws_lin_opt_calc( const rfx_ctrl_ws_t *ws,
                               const rfx_ctrl_ws_lin_k_t *k_t, double dt,
                               double *tJ_sdt, double *tddx_r,
                               double *Ndq_rn ) {

    // Calculate error in workspace position
    double x_e[6];
    {
        double twist[8], de[8];
        aa_tf_duqu_mulc( ws->act.S, ws->ref.S, de );   // de = d*conj(d_r)
        aa_tf_duqu_minimize( de );
        aa_tf_duqu_ln( de, twist );                    // twist = log( de )
        aa_tf_duqu_twist2vel( ws->act.S, twist, x_e );
    }

    // Calculate sign matrix and ~ddx_r
    double M[6 * 6];
    AA_MEM_ZERO( M, 6 * 6 );
    {
        size_t idx;
        for( idx = 0; idx < 6; idx++ ) {
            // Find current acceleration
            double dx_r = ws->ref.dx[idx] - k_t->p[idx] * x_e[idx];
            tddx_r[idx] = (dx_r - ws->act.dx[idx]) / dt;

            // Create bookkeeping matrix of sign and make vector positive
            double s = (double) RFX_SIGN( tddx_r[idx] );
            AA_MATREF( M, 6, idx, idx ) = s;
            tddx_r[idx] *= s;
        }
    }

    double J_s[6 * ws->n_q];
    AA_MEM_ZERO( J_s, 6 * ws->n_q );
    AA_MEM_ZERO( tJ_sdt, 6 * ws->n_q );

    // Compute a damped pseudo inverse
    if( k_t->s2min > 0 )
        aa_la_dzdpinv( 6, ws->n_q, k_t->s2min, ws->J, J_s );
    else
        aa_la_dpinv( 6, ws->n_q, k_t->dls, ws->J, J_s );

    // Calculate the sign preserving J_s and multiply in the constant dt.
    cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
                 (int) ws->n_q, 6, 6,
                 dt, J_s, (int) ws->n_q, M, 6,
                 0, tJ_sdt, (int) ws->n_q );

    // Get change in nullspace velocity to zero joints
    double dq_rn[ws->n_q];
    {
        size_t idx;
        for( idx = 0; idx < ws->n_q; idx++ ) {
            double tzero = ( ws->q_max[idx] + ws->q_min[idx] ) / 2;
            dq_rn[idx] = -ws->act.dq[idx] +
                    (tzero - ws->act.q[idx]) / (ws->q_max[idx] - tzero);
        }
    }

    // Find N dq_rn
    double zero[6];
    AA_MEM_ZERO( zero, 6 );

    aa_la_xlsnp( 6, ws->n_q, ws->J, J_s, zero, dq_rn, Ndq_rn );
}

void rfx_ctrl_ws_lin_opt_cons( const rfx_ctrl_ws_t *ws, double dt,
                               double *c_max, double *c_min ) {
    size_t idx;
    for( idx = 0; idx < ws->n_q; idx++ ) {
        // Velocity constraint
        double vel_max = ws->dq_max[idx] - ws->act.dq[idx];
        double vel_min = ws->dq_min[idx] - ws->act.dq[idx];

        // Acceleration constraint
        double acc_max = ws->ddq_max[idx] * dt;
        double acc_min = ws->ddq_min[idx] * dt;

        // Position constraint
        double pos_max =
                (ws->q_max[idx] - ws->act.q[idx] - (ws->ddq_min[idx] * dt * dt) / 2.)
                / dt - ws->act.dq[idx];

        double pos_min =
                (ws->q_min[idx] - ws->act.q[idx] - (ws->ddq_max[idx] * dt * dt) / 2.)
                / dt - ws->act.dq[idx];

        c_min[idx] = fmax(vel_min, fmax(pos_min, acc_min));
        c_max[idx] = fmin(vel_max, fmin(pos_max, acc_max));
    }
}

AA_API void rfx_ctrl_ws_lin_opt_lp_init( const rfx_ctrl_ws_t *ws, lprec **lp ) {
    *lp = make_lp( (int) ws->n_q * 2, 6 + 1 );

    set_maxim(*lp);
    set_verbose(*lp, IMPORTANT);
}

AA_API rfx_status_t rfx_ctrl_ws_lin_opt( const rfx_ctrl_ws_t *ws,
                                         const rfx_ctrl_ws_lin_k_t *k_t,
                                         lprec *lp, double k_max, double C_u,
                                         double dt, double *u ) {
    AA_MEM_ZERO( u, ws->n_q );

    if (dt <= 0)
        return RFX_OK;

    double tJ_sdt[ws->n_q * 6];
    double tddx_r[6];
    double Ndq_rn[ws->n_q];
    rfx_ctrl_ws_lin_opt_calc( ws, k_t, dt, tJ_sdt, tddx_r, Ndq_rn );

    // Find q constraints
    double c_max[ws->n_q], c_min[ws->n_q];
    rfx_ctrl_ws_lin_opt_cons( ws, dt, c_max, c_min );

    // Set up LP problem
    {
        double cons[6 + 1 + 1];

        {
            // Set objective function
            size_t idx;
            for( idx = 1; idx <= 6; idx++ )
                cons[idx] = tddx_r[idx - 1];
            cons[idx] = C_u;

            set_row( lp, 0, cons );
        }

        {
            // Set constraint of tddx_u <= tddx_r
            int idx;
            for( idx = 1; idx <= 6; idx++ )
                set_bounds( lp, idx, 0, tddx_r[idx - 1] );

            // Set constraint of k <= k_max
            set_bounds( lp, 7, 0, k_max );
        }

        // Set all joint constraints
        {
            size_t idx;
            for( idx = 0; idx < ws->n_q; idx++) {
                // Set ~J*dt coeffecient on ~ddx_u
                size_t jdx;
                for( jdx = 1; jdx <= 6; jdx++ )
                    cons[jdx] = AA_MATREF( tJ_sdt, ws->n_q, idx, jdx - 1);

                // Set Ndq_rn coeffecient on k
                cons[jdx] = Ndq_rn[idx];

                int cdx = (int) idx + 1;
                int ddx = cdx + (int) ws->n_q;

                // Set coeffecients
                set_row( lp, cdx, cons );
                set_row( lp, ddx, cons );

                // Set bounding values
                set_rh( lp, cdx, c_max[idx] );
                set_rh( lp, ddx, c_min[idx] );

                // Set constraint type
                set_constr_type( lp, cdx, LE );
                set_constr_type( lp, ddx, GE );
            }
        }
    }

    // Solve problem
    int ret = solve( lp );
    if( ret != 0 ) {
        printf("LP Solver returned %d\n", ret);
        printf("C_max: ");
        aa_dump_vec(stderr, c_max, ws->n_q);
        printf("C_min: ");
        aa_dump_vec(stderr, c_min, ws->n_q);
        return RFX_INVAL;
    }

    // Get return value
    double tddx_uk[6 + 1];
    AA_MEM_ZERO( tddx_uk, 6 + 1 );
    get_variables( lp, tddx_uk );

    {
        double xys = 0;
        double xs = 0;
        double ys = 0;

        {
            size_t idx;
            for( idx = 0; idx < 6; idx++ ) {
                xys += tddx_uk[idx] * tddx_r[idx];

                double xt = fabs(tddx_uk[idx]);
                xs += xt * xt;

                double yt = fabs(tddx_r[idx]);
                ys += yt * yt;
            }
        }

        xys = fabs(xys);
        xys *= xys;

        if( xys > xs * ys ) {
            // Not linearly dependent
        }
    }

    double Ddq_u[ws->n_q];
    AA_MEM_ZERO( Ddq_u, ws->n_q );
    aa_la_mvmul( ws->n_q, 6, tJ_sdt, tddx_uk, Ddq_u );

    {
        size_t idx;
        for( idx = 0; idx < ws->n_q; idx++ )
            u[idx] = ws->act.dq[idx] + Ddq_u[idx] + Ndq_rn[idx] * tddx_uk[6];
    }

    return RFX_OK;
}

rfx_status_t rfx_ctrl_ws_sdx( rfx_ctrl_ws_t *ws, double dt ) {

    // translation
    /* double r1_split[4], v1_split[3]; */
    /* AA_MEM_CPY( v1_split, ws->x_r, 3 ); */
    /* aa_la_axpy( 3, dt, ws->dx_r, v1_split ); */
    /* // rotation */
    /* aa_tf_qsvel( ws->r_r, ws->dx_r+3, dt, r1_split ); */
    /* aa_tf_qnormalize( r1_split ); */

    double S1[8];
    aa_tf_duqu_svel( ws->ref.S, ws->ref.dx, dt, S1 );
    AA_MEM_CPY( ws->ref.S, S1, 8 );

    /* printf("old:  %f\t%f\t%f\t%f\t|\t%f\t%f\t%f\n", */
    /*        ws->r_r[0], ws->r_r[1], ws->r_r[2], ws->r_r[3], */
    /*        ws->x_r[0], ws->x_r[1], ws->x_r[2] ); */
    /* printf("splt: %f\t%f\t%f\t%f\t|\t%f\t%f\t%f\n", */
    /*        r1_split[0], r1_split[1], r1_split[2], r1_split[3], */
    /*        v1_split[0], v1_split[1], v1_split[2] ); */
    /* printf("dual: %f\t%f\t%f\t%f\t|\t%f\t%f\t%f\n", */
    /*        S1[0], S1[1], S1[2], S1[3], */
    /*        v1_duqu[0], v1_duqu[1], v1_duqu[2] ); */

    return RFX_OK;
}


rfx_status_t rfx_ctrlq_lin_vfwd( const rfx_ctrl_t *g, const rfx_ctrlq_lin_k_t *k,  double *u ) {
    for( size_t i = 0; i < g->n_q; i ++ ) {
        u[i] = g->ref.dq[i]
            - k->p[i] * (g->act.q[i] - g->ref.q[i]);
    }
    return RFX_OK;
}


rfx_ctrlx_lin_t *rfx_ctrlx_lin_alloc( aa_mem_region_t *reg, size_t n_q, rfx_kin_fun kin_fun, void *kin_fun_cx ) {
    rfx_ctrlx_lin_t *p = AA_MEM_REGION_NEW( reg, rfx_ctrlx_lin_t );
    p->ctrl = AA_MEM_REGION_NEW( reg, rfx_ctrl_t );
    p->k = AA_MEM_REGION_NEW( reg, rfx_ctrl_ws_lin_k_t  );

    memset( p->ctrl, 0, sizeof(*p->ctrl) );
    memset( p->k, 0, sizeof(*p->k) );

    p->ctrl->n_q = n_q;
    p->kin_fun = kin_fun;
    p->kin_fun_cx = kin_fun_cx;

    rfx_ctrlx_state_init_region( &p->ctrl->act, reg, n_q );
    rfx_ctrlx_state_init_region( &p->ctrl->ref, reg, n_q );

    p->ctrl->q_min = AA_MEM_REGION_NEW_N( reg, double, n_q );
    p->ctrl->q_max = AA_MEM_REGION_NEW_N( reg, double, n_q );

    AA_MEM_SET( p->ctrl->q_min, 0, n_q );
    AA_MEM_SET( p->ctrl->q_max, 0, n_q );

    p->ctrl->J = AA_MEM_REGION_NEW_N( reg, double, n_q*6 );

    p->k->n_q = n_q;
    p->k->q = AA_MEM_REGION_NEW_N( reg, double, n_q );
    AA_MEM_SET( p->k->q, 0, n_q );


    return p;
}

AA_API rfx_status_t rfx_ctrlx_lin_vfwd( const rfx_ctrlx_lin_t *ctrl, const double *q,
                                        double *u ) {
    AA_MEM_CPY( ctrl->ctrl->act.q, q, ctrl->ctrl->n_q );

    double E[7];
    ctrl->kin_fun( ctrl->kin_fun_cx, q, E, ctrl->ctrl->J );
    aa_tf_qutr2duqu( E, ctrl->ctrl->act.S );
    return rfx_ctrl_ws_lin_vfwd( ctrl->ctrl, ctrl->k, u );
}

int rfx_ctrlx_fun_lin_vfwd ( void *cx,
                             const double *q_a, const double *dq_a,
                             const double *E_r, const double *dx_r,
                             double *dq_r )
{
    const rfx_ctrlx_lin_t *ctrlx = (const rfx_ctrlx_lin_t *)cx;
    // actual
    rfx_ctrl_t *ctrl = ctrlx->ctrl;
    size_t nq = ctrl->n_q;
    AA_MEM_CPY(ctrl->act.q, q_a, nq);
    AA_MEM_CPY(ctrl->act.dq, dq_a, nq);
    double E[7];
    ctrlx->kin_fun( ctrlx->kin_fun_cx, q_a, E, ctrl->J );
    aa_tf_qutr2duqu( E, ctrl->act.S );

    // reference
    aa_tf_qutr2duqu( E_r, ctrl->ref.S );
    AA_MEM_CPY(ctrl->ref.dx, dx_r, 6 );

    // result
    return rfx_ctrl_ws_lin_vfwd( ctrlx->ctrl, ctrlx->k, dq_r );
}

AA_API void
rfx_ctrl_update_act( rfx_ctrl_t *G,
                     const size_t *idx_q, const size_t *idx_frame, int idx_frame_ee,
                     const double *axes,
                     const double *q, size_t incQ,
                     const double *dq, size_t incdQ,
                     const double *E_abs, size_t ldE )
{

    /* Copy joint state */
    for( size_t i = 0; i < G->n_q; i++ ) {
        size_t j = idx_q[i] * incQ;
        size_t k = idx_q[i] * incdQ;
        G->act.q[i] = q[j];
        G->act.dq[i] = dq[k];
    }

    /* compute end-effector dual quaternion */
    const double *E_ee = E_abs + 7*idx_frame_ee;
    aa_tf_qutr2duqu( E_ee, G->act.S );

    const double *pe = & E_ee[AA_TF_QUTR_T];

    /* Compute Jacobian */
    rfx_tf_rev_jacobian( E_abs, ldE,
                         axes,
                         G->n_q,
                         idx_frame, idx_q,
                         pe, G->J, 6 );

}
