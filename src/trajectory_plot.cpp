/* -*- mode: C++; c-basic-offset: 4 -*- */
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
#include <amino.h>
#include "reflex.hpp"

using namespace reflex;
using namespace std;

void WorkspaceTrajectory::plot( double dt ) {
    (void)dt;

    size_t n = (size_t) ( (this->t_n - this->t_0)/dt );
    double T[n]; // time
    // matrices for points in the trajectory.
    // each column is a point.
    // the rows are the time series for each axis.
    double X[n*6];     // pos
    double dX[n*6];    // vel
    double sdX[n*6];   // integrated vel (pos)
    double ddX[n*6];   // acc
    double sddX[n*6];  // integrated acc (vel)

    // compute positions
    {
        double t;
        size_t i;
        for( i = 0, t = this->t_0; t< this->t_n && i < n; i++, t+=dt ) {
            T[i] = t;
            double *Xp = X+i*6;
            double r[4];
            this->get_x(t, Xp, r );
            aa_tf_quat2rotvec( r, Xp + 3 );
            this->get_dx(t, dX+6*i );
            this->get_ddx(t, ddX+6*i );
        }
    }
    // integrate
    {
        aa_fzero( sdX, 6 );
        aa_fzero( sddX, 6 );

        for( size_t i = 1; i < n; i ++ ) {
            size_t j = 6*i;
            size_t k = 6*(i-1);
            // using basic euler-integration
            // vel
            aa_la_axpy3( 6, dt, dX+k, sdX+k,  sdX+j );
            // acc
            aa_la_axpy3( 6, dt, ddX+k, sdX+k, sddX+j );
        }
    }

    { // plot workspace position
        aa_plot_opts_t opts;
        opts.title="Position";
        opts.ylabel="position";
        opts.xlabel="time";
        static const char *c[] = {"x", "y", "z", "r_x", "r_y", "r_z"};
        opts.axis_label = c;
        aa_plot_row_series( 6, n, T, X,
                            & opts );

    }
    { // plot workspace velocity
        aa_plot_opts_t opts;
        opts.title="Velocity";
        opts.ylabel="velocity";
        opts.xlabel="time";
        static const char *c[] = {"dx", "dy", "dz", "dr_x", "dr_y", "dr_z"};
        opts.axis_label = c;
        aa_plot_row_series( 6, n, T, dX,
                            & opts );

    }
    { // plot integrated workspace velocity
        aa_plot_opts_t opts;
        opts.title="Integrated Velocity (Position)";
        opts.ylabel="position";
        opts.xlabel="time";
        static const char * c[] = {"x", "y", "z", "r_x", "r_y", "r_z"};
        opts.axis_label = c;
        aa_plot_row_series( 6, n, T, sdX,
                            & opts );

    }
    { // plot workspace acceleration
        aa_plot_opts_t opts;
        opts.title="Acceleration";
        opts.ylabel="acceleration";
        opts.xlabel="time";
        static const char *c[] = {"ddx", "ddy", "ddz", "ddr_x", "ddr_y", "ddr_z"};
        opts.axis_label = c;
        aa_plot_row_series( 6, n, T, ddX,
                            & opts );

    }
    { // plot workspace velocity
        aa_plot_opts_t opts;
        opts.title="Integrated Acceleration (Velocity)";
        opts.ylabel="velocity";
        opts.xlabel="time";
        static const char *c[] = {"dx", "dy", "dz", "dr_x", "dr_y", "dr_z"};
        opts.axis_label = c;
        aa_plot_row_series( 6, n, T, sddX,
                            & opts );

    }
    { // 3-d plot workspace position
        int r;
        FILE *g = popen("gnuplot -persist", "w");

        assert(g);

        fprintf(g, "set xlabel 'X'\n");
        fprintf(g, "set ylabel 'Y'\n");
        fprintf(g, "set zlabel 'Z'\n");
        fprintf(g, "set title 'Workspace Path'\n");
        fprintf(g, "splot '-' with points title 'Path'");
        fprintf(g, "\n");
        for(size_t i = 0; i < n; i++ ) {
            fprintf(g, "%f %f %f\n",
                    AA_MATREF(X,6,0,i),
                    AA_MATREF(X,6,1,i),
                    AA_MATREF(X,6,2,i) );
        }
        fprintf(g, "e\n");
        fflush(g);
        r = pclose(g);
        assert( -1 != r );
    }
}

void WorkspaceTrajectory::jacobian_plot( double dt,
                                         size_t m, const double *q0,
                                         void fkfun(const double *q,
                                                    double R[9], double v[3] ),
                                         void jfun(const double *q, double *J),
                                         double k_dls ) {
    size_t n = (size_t) ( (this->t_n - this->t_0)/dt );
    double Q[m*n], dQ[m*n], T[n], X0[6*n], dX0[6*n], X[6*n], dX[6*n];
    AA_ZERO_AR(Q);
    AA_ZERO_AR(dQ);
    AA_ZERO_AR(T);
    AA_ZERO_AR(X0);
    AA_ZERO_AR(dX0);
    AA_ZERO_AR(X);
    AA_ZERO_AR(dX);
    // compute joint path
    {
        double t;
        size_t i;
        aa_fcpy(Q, q0, m);
        for( i = 0, t = t_0; i < n-1; i ++, t+=dt ) {
            double J[6*m];
            size_t j = m*i, k = m*(i+1);
            T[i] = t;
            // get joint velocities
            jfun(Q+j, J);  // compute jacobian
            this->get_dx(t, dX0 + 6*i);
            aa_la_dls( 6, m, k_dls, J, dX0 + 6*i, dQ+j );  // compute dq
            // integrate
            aa_la_axpy3( m, dt, dQ+j, Q+j, Q+k );
            // original position
            double rq[4];
            this->get_x(t,X0+6*i, rq);
            aa_tf_quat2rotvec( rq, X0 + 6*i + 3 );
        }
    }
    // compute end-effector path
    {
        for( size_t i = 0; i < n; i++ ) {
            double J[6*m];
            double R[9], v[3];
            size_t iq = m*i;
            size_t ix = 6*i;
            // fk
            fkfun( Q+iq, R, v );
            aa_fcpy( X+ix, v, 3 );
            aa_tf_rotmat2rotvec( R, X+ix+3 );
            // back out workspace velocity
            jfun(Q+iq, J);  // compute jacobian
            aa_la_mvmul( 6, m, J, dQ+iq, dX+ix );
        }
    }

    // plot workspace velocity
    {
        aa_plot_opts_t opts;
        opts.title="WS Velocity: Trajectory";
        opts.ylabel="velocity";
        opts.xlabel="time";
        static const char *c[] = {"dx", "dy", "dz", "dr_x", "dr_y", "dr_z"};
        opts.axis_label = c;
        aa_plot_row_series( 6, n, T, dX0,
                            & opts );
    }

    // plot joint position
    {
        aa_plot_opts_t opts;
        opts.title="JS Positions: Velocity Integral";
        opts.ylabel="position";
        opts.xlabel="time";
        char *axes[m];
        for( size_t i = 0; i < m; i ++ ) {
            axes[i] = (char*)alloca(16);
            snprintf(axes[i], 16, "q_%i",  i );
        }
        opts.axis_label = (const char**)axes;
        aa_plot_row_series( m, n - 1 , T, Q, //FIXME: last position entry is broken
                            & opts );
    }

    // plot joint velocity
    {
        aa_plot_opts_t opts;
        opts.title="Joint Velocity: J^* dx";
        opts.ylabel="velocity";
        opts.xlabel="time";
        char *axes[m];
        for( size_t i = 0; i < m; i ++ ) {
            axes[i] = (char*)alloca(16);
            snprintf(axes[i], 16, "dq_%i",  i );
        }
        opts.axis_label = (const char**)axes;
        aa_plot_row_series( m, n, T, dQ,
                            & opts );
    }
    { // plot workspace position
        aa_plot_opts_t opts;
        opts.title="WS Position: Trajectory";
        opts.ylabel="position";
        opts.xlabel="time";
        static const char *c[] = {"x", "y", "z", "r_x", "r_y", "r_z"};
        opts.axis_label = c;
        aa_plot_row_series( 6, n - 1, T, X0,
                            & opts );

    }
    { // plot workspace position
        aa_plot_opts_t opts;
        opts.title="WS Position: FK";
        opts.ylabel="position";
        opts.xlabel="time";
        static const char *c[] = {"x", "y", "z", "r_x", "r_y", "r_z"};
        opts.axis_label = c;
        aa_plot_row_series( 6, n - 1, T, X,
                            & opts );

    }
    { // plot workspace velocity
        aa_plot_opts_t opts;
        opts.title="Workspace Velocity: dx = J dq";
        opts.ylabel="velocity";
        opts.xlabel="time";
        static const char *c[] = {"dx", "dy", "dz", "dr_x", "dr_y", "dr_z"};
        opts.axis_label = c;
        aa_plot_row_series( 6, n, T, dX,
                            & opts );

    }
    { // 3-d plot workspace position
        int r;
        FILE *g = popen("gnuplot -persist", "w");

        assert(g);

        fprintf(g, "set xlabel 'X'\n");
        fprintf(g, "set ylabel 'Y'\n");
        fprintf(g, "set zlabel 'Z'\n");
        fprintf(g, "set title 'Workspace Path: FK'\n");
        fprintf(g, "splot '-' with points title 'Path'");
        fprintf(g, "\n");
        for(size_t i = 0; i < n-1; i++ ) {
            fprintf(g, "%f %f %f\n",
                    AA_MATREF(X,6,0,i),
                    AA_MATREF(X,6,1,i),
                    AA_MATREF(X,6,2,i) );
        }
        fprintf(g, "e\n");
        fflush(g);
        r = pclose(g);
        assert( -1 != r );
    }

}

void WorkspaceTrajectory::wsctrl_plot( double dt,
                                       size_t m, const double *q0,
                                       void fkfun(const double *q,
                                                  double R[9], double v[3] ),
                                       void jfun(const double *q, double *J),
                                       double k_dls, double k_p[6] ) {
    // init controller
    rfx_ctrl_ws_t G;
    rfx_ctrl_ws_lin_k_t kappa;
    aa_fcpy( kappa.p, k_p, 6 );
    aa_fzero( kappa.f, 6 );
    kappa.dls = k_dls;
    rfx_ctrl_ws_init(&G, m);

    size_t n = (size_t) ( (this->t_n - this->t_0 + 2)/dt );
    //n+=50;
    double Q[m*n], dQ[m*n], T[n], X0[6*n], dX0[6*n], X[6*n], dX[6*n];
    n--;
    AA_ZERO_AR(Q);
    AA_ZERO_AR(dQ);
    AA_ZERO_AR(T);
    AA_ZERO_AR(X0);
    AA_ZERO_AR(dX0);
    AA_ZERO_AR(X);
    AA_ZERO_AR(dX);
    // compute path
    {
        double t;
        size_t i;
        aa_fcpy(Q, q0, m);
        for( i = 0, t = t_0; i < n; i ++, t+=dt ) {
            size_t j = m*i, k = m*(i+1);
            T[i] = t;
            // fill in controller state
            aa_fcpy(G.q, Q+j, m);  // joint config
            jfun(G.q, G.J);  // jacobian
            double R[9];
            fkfun(G.q, R, G.x);
            aa_tf_rotmat2quat(R, G.r);  // workspace position
            // fill in controller references
            this->get_x( t, G.x_r, G.r_r );
            this->get_dx( t, G.dx_r );
            // run controller
            rfx_ctrl_ws_lin_vfwd( &G, &kappa, dQ+j );
            // integrate
            aa_la_axpy3( m, dt, dQ+j, Q+j, Q+k );
            // fill in stored state
            aa_fcpy(X0+6*i, G.x_r, 3);
            aa_fcpy(X+6*i, G.x, 3);
            aa_fcpy(dX0+6*i, G.dx_r, 6);
            if( i > 0 ) {
                aa_tf_quat2rotvec_near( G.r_r, X0 + 6*(i-1) + 3, X0 + 6*i + 3 );
                aa_tf_quat2rotvec_near( G.r,    X + 6*(i-1) + 3,  X + 6*i + 3 );
            } else {
                aa_tf_quat2rotvec( G.r_r, X0 + 6*i + 3 );
                aa_tf_quat2rotvec( G.r, X + 6*i + 3 );
            }
            aa_la_mvmul( 6, m, G.J, dQ+j, dX+6*i );
        }
    }

    // plot workspace velocity
    {
        aa_plot_opts_t opts;
        opts.title="WS Velocity: Trajectory";
        opts.ylabel="velocity";
        opts.xlabel="time";
        static const char *c[] = {"dx", "dy", "dz", "dr_x", "dr_y", "dr_z"};
        opts.axis_label = c;
        aa_plot_row_series( 6, n, T, dX0,
                            & opts );
    }

    // plot joint position
    {
        aa_plot_opts_t opts;
        opts.title="JS Positions: Velocity Integral";
        opts.ylabel="position";
        opts.xlabel="time";
        char *axes[m];
        for( size_t i = 0; i < m; i ++ ) {
            axes[i] = (char*)alloca(16);
            snprintf(axes[i], 16, "q_%i",  i );
        }
        opts.axis_label = (const char**)axes;
        aa_plot_row_series( m, n  , T, Q, //FIXME: last position entry is broken
                            & opts );
    }

    // plot joint velocity
    {
        aa_plot_opts_t opts;
        opts.title="Joint Velocity: J^* dx";
        opts.ylabel="velocity";
        opts.xlabel="time";
        char *axes[m];
        for( size_t i = 0; i < m; i ++ ) {
            axes[i] = (char*)alloca(16);
            snprintf(axes[i], 16, "dq_%i",  i );
        }
        opts.axis_label = (const char**)axes;
        aa_plot_row_series( m, n-1, T, dQ,
                            & opts );
    }
    { // plot workspace position
        aa_plot_opts_t opts;
        opts.title="WS Position: Trajectory";
        opts.ylabel="position";
        opts.xlabel="time";
        static const char *c[] = {"x", "y", "z", "r_x", "r_y", "r_z"};
        opts.axis_label = c;
        aa_plot_row_series( 6, n , T, X0,
                            & opts );

    }
    { // plot workspace position
        aa_plot_opts_t opts;
        opts.title="WS Position: FK";
        opts.ylabel="position";
        opts.xlabel="time";
        static const char *c[] = {"x", "y", "z", "r_x", "r_y", "r_z"};
        opts.axis_label = c;
        aa_plot_row_series( 6, n , T, X,
                            & opts );

    }
    { // plot workspace velocity
        aa_plot_opts_t opts;
        opts.title="Workspace Velocity: dx = J dq";
        opts.ylabel="velocity";
        opts.xlabel="time";
        static const char *c[] = {"dx", "dy", "dz", "dr_x", "dr_y", "dr_z"};
        opts.axis_label = c;
        aa_plot_row_series( 6, n, T, dX,
                            & opts );

    }
    { // 3-d plot workspace position
        int r;
        FILE *g = popen("gnuplot -persist", "w");

        assert(g);

        fprintf(g, "set xlabel 'X'\n");
        fprintf(g, "set ylabel 'Y'\n");
        fprintf(g, "set zlabel 'Z'\n");
        fprintf(g, "set title 'Workspace Path: FK'\n");
        fprintf(g, "splot '-' with points title 'Path'");
        fprintf(g, "\n");
        for(size_t i = 0; i < n; i++ ) {
            fprintf(g, "%f %f %f\n",
                    AA_MATREF(X,6,0,i),
                    AA_MATREF(X,6,1,i),
                    AA_MATREF(X,6,2,i) );
        }
        fprintf(g, "e\n");
        fflush(g);
        r = pclose(g);
        assert( -1 != r );
    }
    rfx_ctrl_ws_destroy(&G);
}
