/* -*- mode: C++; c-basic-offset: 4 -*- */
/* ex: set shiftwidth=4 tabstop=4 expandtab: */
/*
 * Copyright (c) 2010-2013, Georgia Tech Research Corporation
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

void rfx_trajq_seg_plot( struct rfx_trajq_seg_list *cx, double dt ) {
    double t_i = rfx_trajq_seg_list_get_t_i(cx);
    double t_f = rfx_trajq_seg_list_get_t_f(cx);


    size_t n = (size_t) ( (t_f - t_i) / dt );
    size_t n_q = rfx_trajq_seg_list_get_n_q(cx) ;
    double T[n]; // time
    // matrices for points in the trajectory.
    // each column is a point.
    // the rows are the time series for each axis.
    double Q[n*n_q];     //pos
    double dQ[n*n_q];    // vel
    double ddQ[n*n_q];   // acceleration
    double Qi[n*n_q];    // integrated vel (pos)
    /* double ddX[n*6];  // acc */
    /* double sddX[n*6]; // integrated acc (vel) */

    // get actuals
    {
        double t;
        size_t i;
        for( i = 0, t = t_i; t < t_f && i < n; i++, t+=dt ) {
            T[i] = t;
            rfx_trajq_seg_list_get_ddq( cx, t, Q + i*n_q, dQ + i*n_q,  ddQ + i*n_q );
        }
    }
    // integrate
    {
        //aa_fzero( sdX, 6 );
        //aa_fzero( sddX, 6 );
        AA_MEM_CPY( Qi, Q, n_q );

        for( size_t k = 1; k < n; k ++ ) {
            //qi(i) = qi(i-1) + dt*dq(i)
            for( size_t i = 0; i < n_q; i ++ ) {
                Qi[ k*n_q + i ] = Qi[ (k-1)*n_q + i ] + dt*dQ[ k*n_q + i ];
            }
        }
    }

    { // plot position
        aa_plot_opts_t opts = {0};
        opts.title="Position";
        opts.ylabel="position";
        opts.xlabel="time";
        aa_plot_row_series( n_q, n, T, Q,
                            & opts );

    }

    { // plot position
        aa_plot_opts_t opts = {0};
        opts.title="Position (integrated velocity)";
        opts.ylabel="position";
        opts.xlabel="time";
        aa_plot_row_series( n_q, n, T, Qi,
                            & opts );

    }

    { // plot velocity position
        aa_plot_opts_t opts = {0};
        opts.title="Velocity";
        opts.ylabel="velocity";
        opts.xlabel="time";
        aa_plot_row_series( n_q, n, T, dQ,
                            & opts );

    }
}


void rfx_trajq_plot( struct rfx_trajq *cx, double dt ) {
    double t_i = cx->t_i;
    double t_f = cx->t_f;


    size_t n = (size_t) ( (t_f - t_i) / dt );
    size_t n_q = cx->n_q;
    double T[n]; // time
    // matrices for points in the trajectory.
    // each column is a point.
    // the rows are the time series for each axis.
    double Q[n*n_q];     //pos
    double dQ[n*n_q];    // vel
    double ddQ[n*n_q];   // acceleration
    double Qi[n*n_q];    // integrated vel (pos)
    /* double ddX[n*6];  // acc */
    /* double sddX[n*6]; // integrated acc (vel) */

    // get actuals
    {
        double t;
        size_t i;
        for( i = 0, t = t_i; t < t_f && i < n; i++, t+=dt ) {
            T[i] = t;
            rfx_trajq_get_q( cx, t, Q + i*n_q );
            rfx_trajq_get_dq( cx, t, dQ + i*n_q );
            rfx_trajq_get_ddq( cx, t, ddQ + i*n_q );
        }
    }
    // integrate
    {
        //aa_fzero( sdX, 6 );
        //aa_fzero( sddX, 6 );
        AA_MEM_CPY( Qi, cx->q_i, n_q );

        for( size_t k = 1; k < n; k ++ ) {
            //qi(i) = qi(i-1) + dt*dq(i)
            for( size_t i = 0; i < n_q; i ++ ) {
                Qi[ k*n_q + i ] = Qi[ (k-1)*n_q + i ] + dt*dQ[ k*n_q + i ];
            }
        }
    }

    { // plot position
        aa_plot_opts_t opts = {0};
        opts.title="Position";
        opts.ylabel="position";
        opts.xlabel="time";
        aa_plot_row_series( n_q, n, T, Q,
                            & opts );

    }

    { // plot position
        aa_plot_opts_t opts = {0};
        opts.title="Position (integrated velocity)";
        opts.ylabel="position";
        opts.xlabel="time";
        aa_plot_row_series( n_q, n, T, Qi,
                            & opts );

    }

    { // plot velocity position
        aa_plot_opts_t opts = {0};
        opts.title="Velocity";
        opts.ylabel="velocity";
        opts.xlabel="time";
        aa_plot_row_series( n_q, n, T, dQ,
                            & opts );

    }

/*     { // plot workspace velocity */
/*         aa_plot_opts_t opts; */
/*         opts.title="Velocity"; */
/*         opts.ylabel="velocity"; */
/*         opts.xlabel="time"; */
/*         static const char *c[] = {"dx", "dy", "dz", "dr_x", "dr_y", "dr_z"}; */
/*         opts.axis_label = c; */
/*         aa_plot_row_series( 6, n, T, dX, */
/*                             & opts ); */

/*     } */
/*     { // plot integrated workspace velocity */
/*         aa_plot_opts_t opts; */
/*         opts.title="Integrated Velocity (Position)"; */
/*         opts.ylabel="position"; */
/*         opts.xlabel="time"; */
/*         static const char * c[] = {"x", "y", "z", "r_x", "r_y", "r_z"}; */
/*         opts.axis_label = c; */
/*         aa_plot_row_series( 6, n, T, sdX, */
/*                             & opts ); */

/*     } */
/*     { // plot workspace acceleration */
/*         aa_plot_opts_t opts; */
/*         opts.title="Acceleration"; */
/*         opts.ylabel="acceleration"; */
/*         opts.xlabel="time"; */
/*         static const char *c[] = {"ddx", "ddy", "ddz", "ddr_x", "ddr_y", "ddr_z"}; */
/*         opts.axis_label = c; */
/*         aa_plot_row_series( 6, n, T, ddX, */
/*                             & opts ); */

/*     } */
/*     { // plot workspace velocity */
/*         aa_plot_opts_t opts; */
/*         opts.title="Integrated Acceleration (Velocity)"; */
/*         opts.ylabel="velocity"; */
/*         opts.xlabel="time"; */
/*         static const char *c[] = {"dx", "dy", "dz", "dr_x", "dr_y", "dr_z"}; */
/*         opts.axis_label = c; */
/*         aa_plot_row_series( 6, n, T, sddX, */
/*                             & opts ); */

/*     } */
/*     { // 3-d plot workspace position */
/*         int r; */
/*         FILE *g = popen("gnuplot -persist", "w"); */

/*         assert(g); */

/*         fprintf(g, "set xlabel 'X'\n"); */
/*         fprintf(g, "set ylabel 'Y'\n"); */
/*         fprintf(g, "set zlabel 'Z'\n"); */
/*         fprintf(g, "set title 'Workspace Path'\n"); */
/*         fprintf(g, "splot '-' with points title 'Path'"); */
/*         fprintf(g, "\n"); */
/*         for(size_t i = 0; i < n; i++ ) { */
/*             fprintf(g, "%f %f %f\n", */
/*                     AA_MATREF(X,6,0,i), */
/*                     AA_MATREF(X,6,1,i), */
/*                     AA_MATREF(X,6,2,i) ); */
/*         } */
/*         fprintf(g, "e\n"); */
/*         fflush(g); */
/*         r = pclose(g); */
/*         assert( -1 != r ); */
/*     } */
/* } */

}

static struct rfx_trajx_plot_opts default_trajx_plot_opts = {0};


void rfx_trajx_plot( struct rfx_trajx *cx, double dt, const struct rfx_trajx_plot_opts *xopts ) {

    double t_i = cx->pt_i->t;
    double t_f = cx->pt_f->t;

    if( NULL == xopts ) xopts = &default_trajx_plot_opts;

    size_t n = (size_t) ( (t_f - t_i) / dt );
    double T[n]; // time
    // matrices for points in the trajectory.
    // each column is a point.
    // the rows are the time series for each axis.
    double X[n*3];     // translation
    double Q[n*4];     // Orientation
    double dXt[n*3];   // translational vel
    double dXr[n*3];   // rotational vel

    double sdQ[n*4];   // Integrated Quatenion Derivative (orientation)
    double sdX[n*3];   // Integrated translation

    double dQ[n*4];    // Quaternion Derivative

    /* double ddX[n*6];  // acc */
    /* double sddX[n*6]; // integrated acc (vel) */

    // get actuals
    {
        double t;
        size_t i;
        for( i = 0, t = t_i; t < t_f && i < n; i++, t+=dt ) {
            T[i] = t;
            rfx_trajx_get_x( cx, t, X + 3*i, Q+4*i );

            double dx[6];
            rfx_trajx_get_dx( cx, t, dx );
            AA_MEM_CPY( dXt + 3*i, dx, 3 );
            AA_MEM_CPY( dXr + 3*i, dx+3, 3 );
            aa_tf_qvel2diff( Q+4*i, dx+3, dQ+4*i );
        }
    }
    // integrate
    {
        //aa_fzero( sdX, 6 );
        //aa_fzero( sddX, 6 );
        AA_MEM_CPY( sdQ, cx->pt_i->r, 4 );
        AA_MEM_CPY( sdX, cx->pt_i->x, 3 );

        for( size_t k = 1; k < n; k ++ ) {
            // orientation
            for( size_t i = 0; i < 4; i ++ ) {
                sdQ[ k*4 + i ] = sdQ[ (k-1)*4 + i ] + dt*dQ[ k*4 + i ];
            }

            // translation
            for( size_t i = 0; i < 3; i ++ ) {
                sdX[ k*3 + i ] = sdX[ (k-1)*3 + i ] + dt*dXt[ k*3 + i ];
            }
        }
    }


    static const char *XYZW[] = {"x", "y", "z", "w"};
    { // plot position
        aa_plot_opts_t opts = {0};
        opts.title="Translation";
        opts.ylabel="translation (m)";
        opts.xlabel="time (s)";
        opts.axis_label = XYZW;
        if( xopts->to_file ) opts.script_file = "x.gnuplot";
        aa_plot_row_series( 3, n, T, X,
                            & opts );

    }

    { // plot orientation
        aa_plot_opts_t opts = {0};
        opts.title="Orientation";
        opts.ylabel="orientation";
        opts.xlabel="time (s)";
        opts.axis_label = XYZW;
        if( xopts->to_file ) opts.script_file = "q.gnuplot";
        aa_plot_row_series( 4, n, T, Q,
                            & opts );

    }


    { // plot integrated quaternion derivative
        aa_plot_opts_t opts = {0};
        opts.title="Orientation (integrated derivative)";
        opts.ylabel="quaternion";
        opts.xlabel="time (s)";
        opts.axis_label = XYZW;
        if( xopts->to_file ) opts.script_file = "sum_dq.gnuplot";
        aa_plot_row_series( 4, n, T, sdQ,
                            & opts );

    }

    { // plot integrated translational velocity
        aa_plot_opts_t opts = {0};
        opts.title="Position (integrated velocity)";
        opts.ylabel="position (m)";
        opts.xlabel="time (s)";
        opts.axis_label = XYZW;
        if( xopts->to_file ) opts.script_file = "sum_dx.gnuplot";
        aa_plot_row_series( 3, n, T, sdX,
                            & opts );

    }

    { // plot translational velocity
        aa_plot_opts_t opts = {0};
        opts.title="Translational Velocity";
        opts.ylabel="velocity (m/s)";
        opts.xlabel="time (s)";
        opts.axis_label = XYZW;
        if( xopts->to_file ) opts.script_file = "dx.gnuplot";
        aa_plot_row_series( 3, n, T, dXt,
                            & opts );

    }

    { // plot rotational velocity
        aa_plot_opts_t opts = {0};
        opts.title="Rotational Velocity";
        opts.ylabel="velocity (rad/s)";
        opts.xlabel="time (s)";
        opts.axis_label = XYZW;
        if( xopts->to_file ) opts.script_file = "omega.gnuplot";
        aa_plot_row_series( 3, n, T, dXr,
                            & opts );

    }
    { // plot quaternion derivative
        aa_plot_opts_t opts = {0};
        opts.title="Orientation Derivative";
        opts.ylabel="quaternion derivative";
        opts.xlabel="time (s)";
        opts.axis_label = XYZW;
        if( xopts->to_file ) opts.script_file = "dq.gnuplot";
        aa_plot_row_series( 4, n, T, dQ,
                            & opts );

    }

    // joints

    if ( xopts && xopts->ctrlx && xopts->q_0 ) {
        size_t n_q = xopts->ctrlx->ctrl->n_q;
        double phi[ xopts->ctrlx->ctrl->n_q * (n+1) ];
        double dphi[ n_q * n ];
        AA_MEM_CPY( phi, xopts->q_0, n_q );

        for( size_t i = 0; i < n; i ++ ) {
            rfx_trajx_set_ctrl( cx, T[i], xopts->ctrlx );
            rfx_ctrlx_lin_vfwd( xopts->ctrlx, phi + i*n_q, dphi+i*n_q );
            // integrate
            for( size_t j = 0; j < n_q; j ++ ) {
                phi[ (i+1)*n_q + j ] = phi[i*n_q + j] + dt * dphi[i*n_q + j];
            }
        }
        char lbl[n_q][32];
        const char *plbl[n_q];
        for( size_t j = 0; j < n_q; j ++ ) {
            sprintf( lbl[j], "%lu", j );
            plbl[j] = lbl[j];
        }

        { // plot joint position
            aa_plot_opts_t opts = {0};
            opts.title="Joint Position";
            opts.ylabel="position (rad)";
            opts.xlabel="time (s)";
            opts.axis_label = plbl;
            if( xopts->to_file ) opts.script_file = "phi.gnuplot";
            aa_plot_row_series( n_q, n, T, phi,
                                & opts );

        }
        { // plot joint velocity
            aa_plot_opts_t opts = {0};
            opts.title="Joint Velocity";
            opts.ylabel="velocity (rad/s)";
            opts.xlabel="time (s)";
            opts.axis_label = plbl;
            if( xopts->to_file ) opts.script_file = "dphi.gnuplot";
            aa_plot_row_series( n_q, n, T, dphi,
                                & opts );

        }
    }

}




void rfx_trajx_seglist_plot( struct rfx_trajx_seg_list *cx, double dt, const struct rfx_trajx_plot_opts *xopts ) {

    double t_i = rfx_trajx_seg_list_get_t_i(cx);
    double t_f = rfx_trajx_seg_list_get_t_f(cx);

    if( NULL == xopts ) xopts = &default_trajx_plot_opts;

    size_t n = (size_t) ( (t_f - t_i) / dt );

    double T[n]; // time
    // matrices for points in the trajectory.
    // each column is a point.
    // the rows are the time series for each axis.
    double X[n*3];     // translation
    double Q[n*4];     // Orientation
    double dXt[n*3];   // translational vel
    double dXr[n*3];   // rotational vel

    double sdQ[n*4];   // Integrated Quatenion Derivative (orientation)
    double sdX[n*3];   // Integrated translation

    double dQ[n*4];    // Quaternion Derivative

    /* double ddX[n*6];  // acc */
    /* double sddX[n*6]; // integrated acc (vel) */

    // get actuals
    {
        double t;
        size_t i;
        for( i = 0, t = t_i; t < t_f && i < n; i++, t+=dt ) {
            T[i] = t;
            double dx[6], S[8];
            rfx_trajx_seg_list_get_dx_duqu( cx, t, S, dx );
            aa_tf_duqu2qv( S, Q+4*i, X+3*i );

            AA_MEM_CPY( dXt + 3*i, dx, 3 );
            AA_MEM_CPY( dXr + 3*i, dx+3, 3 );
            aa_tf_qvel2diff( Q+4*i, dx+3, dQ+4*i );
        }
    }
    // integrate
    {
        //aa_fzero( sdX, 6 );
        //aa_fzero( sddX, 6 );
        AA_MEM_CPY( sdQ, Q, 4 );
        AA_MEM_CPY( sdX, X, 3 );

        for( size_t k = 1; k < n; k ++ ) {
            // orientation
            for( size_t i = 0; i < 4; i ++ ) {
                sdQ[ k*4 + i ] = sdQ[ (k-1)*4 + i ] + dt*dQ[ k*4 + i ];
            }
            aa_tf_qnormalize( sdQ + k*4 );

            // translation
            for( size_t i = 0; i < 3; i ++ ) {
                sdX[ k*3 + i ] = sdX[ (k-1)*3 + i ] + dt*dXt[ k*3 + i ];
            }
        }
    }


    static const char *XYZW[] = {"x", "y", "z", "w"};
    { // plot position
        aa_plot_opts_t opts = {0};
        opts.title="Translation";
        opts.ylabel="translation (m)";
        opts.xlabel="time (s)";
        opts.axis_label = XYZW;
        if( xopts->to_file ) opts.script_file = "x.gnuplot";
        aa_plot_row_series( 3, n, T, X,
                            & opts );

    }

    { // plot orientation
        aa_plot_opts_t opts = {0};
        opts.title="Orientation";
        opts.ylabel="orientation";
        opts.xlabel="time (s)";
        opts.axis_label = XYZW;
        if( xopts->to_file ) opts.script_file = "q.gnuplot";
        aa_plot_row_series( 4, n, T, Q,
                            & opts );

    }


    { // plot integrated quaternion derivative
        aa_plot_opts_t opts = {0};
        opts.title="Orientation (integrated derivative)";
        opts.ylabel="quaternion";
        opts.xlabel="time (s)";
        opts.axis_label = XYZW;
        if( xopts->to_file ) opts.script_file = "sum_dq.gnuplot";
        aa_plot_row_series( 4, n, T, sdQ,
                            & opts );

    }

    { // plot integrated translational velocity
        aa_plot_opts_t opts = {0};
        opts.title="Position (integrated velocity)";
        opts.ylabel="position (m)";
        opts.xlabel="time (s)";
        opts.axis_label = XYZW;
        if( xopts->to_file ) opts.script_file = "sum_dx.gnuplot";
        aa_plot_row_series( 3, n, T, sdX,
                            & opts );

    }

    { // plot translational velocity
        aa_plot_opts_t opts = {0};
        opts.title="Translational Velocity";
        opts.ylabel="velocity (m/s)";
        opts.xlabel="time (s)";
        opts.axis_label = XYZW;
        if( xopts->to_file ) opts.script_file = "dx.gnuplot";
        aa_plot_row_series( 3, n, T, dXt,
                            & opts );

    }

    { // plot rotational velocity
        aa_plot_opts_t opts = {0};
        opts.title="Rotational Velocity";
        opts.ylabel="velocity (rad/s)";
        opts.xlabel="time (s)";
        opts.axis_label = XYZW;
        if( xopts->to_file ) opts.script_file = "omega.gnuplot";
        aa_plot_row_series( 3, n, T, dXr,
                            & opts );

    }
    { // plot quaternion derivative
        aa_plot_opts_t opts = {0};
        opts.title="Orientation Derivative";
        opts.ylabel="quaternion derivative";
        opts.xlabel="time (s)";
        opts.axis_label = XYZW;
        if( xopts->to_file ) opts.script_file = "dq.gnuplot";
        aa_plot_row_series( 4, n, T, dQ,
                            & opts );

    }

    // // joints

    // if ( xopts && xopts->ctrlx && xopts->q_0 ) {
    //     size_t n_q = xopts->ctrlx->ctrl->n_q;
    //     double phi[ xopts->ctrlx->ctrl->n_q * (n+1) ];
    //     double dphi[ n_q * n ];
    //     AA_MEM_CPY( phi, xopts->q_0, n_q );

    //     for( size_t i = 0; i < n; i ++ ) {
    //         rfx_trajx_set_ctrl( cx, T[i], xopts->ctrlx );
    //         rfx_ctrlx_lin_vfwd( xopts->ctrlx, phi + i*n_q, dphi+i*n_q );
    //         // integrate
    //         for( size_t j = 0; j < n_q; j ++ ) {
    //             phi[ (i+1)*n_q + j ] = phi[i*n_q + j] + dt * dphi[i*n_q + j];
    //         }
    //     }
    //     char lbl[n_q][32];
    //     const char *plbl[n_q];
    //     for( size_t j = 0; j < n_q; j ++ ) {
    //         sprintf( lbl[j], "%lu", j );
    //         plbl[j] = lbl[j];
    //     }

    //     { // plot joint position
    //         aa_plot_opts_t opts = {0};
    //         opts.title="Joint Position";
    //         opts.ylabel="position (rad)";
    //         opts.xlabel="time (s)";
    //         opts.axis_label = plbl;
    //         if( xopts->to_file ) opts.script_file = "phi.gnuplot";
    //         aa_plot_row_series( n_q, n, T, phi,
    //                             & opts );

    //     }
    //     { // plot joint velocity
    //         aa_plot_opts_t opts = {0};
    //         opts.title="Joint Velocity";
    //         opts.ylabel="velocity (rad/s)";
    //         opts.xlabel="time (s)";
    //         opts.axis_label = plbl;
    //         if( xopts->to_file ) opts.script_file = "dphi.gnuplot";
    //         aa_plot_row_series( n_q, n, T, dphi,
    //                             & opts );

    //     }
    // }

}
