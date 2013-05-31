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

void rfx_trajq_plot( struct rfx_trajq *cx, double dt ) {
    (void)dt;
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
