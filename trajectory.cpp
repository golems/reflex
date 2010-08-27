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
#include "reflex.hpp"

using namespace reflex;
using namespace std;

Trajectory::Trajectory() {}
Trajectory::~Trajectory() {}

WorkspaceTrajectory::WorkspaceTrajectory() {}
WorkspaceTrajectory::~WorkspaceTrajectory() {}

static int trapvel_generate( size_t n, double t_f,
                             const double *x_i, const double *x_f,
                             const double *dx_max, const double *ddx_max,
                             double *tb, double *dx_r, double *ddx_r );

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
        opts.axis_label = (const char*[]){"x", "y", "z", "r_x", "r_y", "r_z"};
        aa_plot_row_series( 6, n, T, X,
                            & opts );

    }
    { // plot workspace velocity
        aa_plot_opts_t opts;
        opts.title="Velocity";
        opts.ylabel="velocity";
        opts.xlabel="time";
        opts.axis_label = (const char*[]){"dx", "dy", "dz", "dr_x", "dr_y", "dr_z"};
        aa_plot_row_series( 6, n, T, dX,
                            & opts );

    }
    { // plot integrated workspace velocity
        aa_plot_opts_t opts;
        opts.title="Integrated Velocity (Position)";
        opts.ylabel="position";
        opts.xlabel="time";
        opts.axis_label = (const char*[]){"x", "y", "z", "r_x", "r_y", "r_z"};
        aa_plot_row_series( 6, n, T, sdX,
                            & opts );

    }
    { // plot workspace acceleration
        aa_plot_opts_t opts;
        opts.title="Acceleration";
        opts.ylabel="acceleration";
        opts.xlabel="time";
        opts.axis_label = (const char*[]){"ddx", "ddy", "ddz", "ddr_x", "ddr_y", "ddr_z"};
        aa_plot_row_series( 6, n, T, ddX,
                            & opts );

    }
    { // plot workspace velocity
        aa_plot_opts_t opts;
        opts.title="Integrated Acceleration (Velocity)";
        opts.ylabel="velocity";
        opts.xlabel="time";
        opts.axis_label = (const char*[]){"dx", "dy", "dz", "dr_x", "dr_y", "dr_z"};
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
                                         void jfun(const double *q, double *J) ) {
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
            aa_la_dls( 6, m, .001, J, dX0 + 6*i, dQ+j );  // compute dq
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
        opts.axis_label = (const char*[]){"dx", "dy", "dz", "dr_x", "dr_y", "dr_z"};
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
        opts.axis_label = (const char*[]){"x", "y", "z", "r_x", "r_y", "r_z"};
        aa_plot_row_series( 6, n - 1, T, X0,
                            & opts );

    }
    { // plot workspace position
        aa_plot_opts_t opts;
        opts.title="WS Position: FK";
        opts.ylabel="position";
        opts.xlabel="time";
        opts.axis_label = (const char*[]){"x", "y", "z", "r_x", "r_y", "r_z"};
        aa_plot_row_series( 6, n - 1, T, X,
                            & opts );

    }
    { // plot workspace velocity
        aa_plot_opts_t opts;
        opts.title="Workspace Velocity: dx = J dq";
        opts.ylabel="velocity";
        opts.xlabel="time";
        opts.axis_label = (const char*[]){"dx", "dy", "dz", "dr_x", "dr_y", "dr_z"};
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


TrapvelWS::TrapvelWS()
{
    t_0 = 0;
    t_n = 0;
    AA_ZERO_AR(x_0);
    AA_ZERO_AR(x_n);

    AA_SET_AR(dx_max, 1.0);
    AA_SET_AR(ddx_max, 1.0);
}

TrapvelWS::~TrapvelWS()
{}

int TrapvelWS::validate() {
    // sanity checks
    if( ! aa_feq( this->t_0, 0, 0) ) return -1;
    if( this->tb <= 0 ) return -1;
    if( this->t_n < this->tb ) return -1;
    for( size_t i = 0; i < 6; i ++ ) {
        if( ! aa_feq( this->tb * ddx_r[i], dx_r[i], .001 ) ) {
            return -1;
        }
    }
    // ok
    return 0;
}


int TrapvelWS::generate() {
    int i = trapvel_generate( 6, this->t_n,
                              this->x_0, this->x_n,
                              this->dx_max, this->ddx_max,
                              &this->tb, this->dx_r, this->ddx_r );
    //fprintf(stderr, "tb: %f\n", tb );
    return i;
}

/* static int trapvel_generate( size_t n, double t_f,  */
/*                              const double *x_i, const double *x_f,  */
/*                              const double *dx_max, const double *ddx_max, */
/*                              double *t1, double *t2,  */
/*                              double *dx_r, double *ddx_r ) { */
/*     double t3=t_f; */
/*     for( size_t i = 0; i < n; i++ ) { */
/*         double d = x_f[i] - x_i[i]; */
/*         // try triangular profile */
/*         t1[i] = t3/2; */
/*         dx_r[i] = d / t1[i]; */
/*         // check velocity */
/*         if( abs(dx_r[i]) <= abs(dx_max[i]) ) { */
/*             // triangle ok */
/*             t2[i] = t1[i]; */
/*             ddx_r[i] = dx_r[i] / t1[i]; */
/*         } else { */
/*             // trapezoid */
/*             dx_r[i] = dx_max[i]; */
/*             double tc = 2*d/dx_r[i] - t3; */
/*             t1[i] = (t3-tc)/2; */
/*             t2[i] = t3 - t1[i]; */
/*             ddx_r[i] = dx_max[i]/t1[i]; */
/*             // check acceleration */
/*             if( abs(ddx_r[i]) > ddx_max[i] ) return -1; */
/*         } */
/*     } */
/*     return 0; */
/* } */

static int trapvel_generate( size_t n, double t_f,
                             const double *x_i, const double *x_f,
                             const double *dx_max, const double *ddx_max,
                             double *ptb, double *dx_r, double *ddx_r ) {
    double t3=t_f;
    double tb = t3/2;
    //fprintf(stderr, "t3: %f\n", t3 );

    double x[n];
    aa_la_vsub(n, x_f, x_i, x);
    //aa_dump_vec(stderr, x, 6 );

    // try triangular profile
    int is_tri = 1;
    for( size_t i = 0; i < n; i++ ) {
        dx_r[i] = x[i] / tb;
        // check for velocity limit
        if( abs(dx_r[i]) > abs(dx_max[i]) ) {
            is_tri = 0;
            break;
        }
        // check for acceleration limit
        ddx_r[i] = dx_r[i] / tb;
        if( abs(ddx_r[i]) > abs(ddx_max[i]) )
            return -1;
    }
    if( is_tri ) {
        //fprintf(stderr, "tri\n");
        *ptb = tb;
        return 0;
    }
    //fprintf(stderr, "trap\n");

    // needs to be trapezoid
    // find longest acceptable blend time
    for( size_t i = 0; i < n; i++ ) {
        double t = t3 - x[i]/dx_max[i];
        tb = AA_MIN(tb,t);
    }
    //fprintf(stderr, "found tb: %f\n", tb );
    // calc dx, ddx
    double t2 = t3 - tb;
    for( size_t i = 0; i < n; i++ ) {
        dx_r[i] = x[i]/t2;
        ddx_r[i] = dx_r[i]/tb;
        // check a
        if( abs(ddx_r[i]) > abs(ddx_max[i]) ||
            abs(dx_r[i])  > abs(dx_max[i]) ) {
            return -1;
        }
    }
    *ptb = tb;
    return 0;
}

int TrapvelWS::get_x( double t, double x[3], double r[4]) {
    double *xp, xps[6];  // temp storage

    if( t < this->t_0 ) {
        // before t0
        xp = this->x_0;
    }
    else if( t > this->t_n ) {
        // after t1
        xp = this->x_n;
    } else {
        //normal
        xp = xps;
        double t1 = tb, t2 = t_n-tb;
        for( size_t i = 0; i < 6; i ++ ) {
            if( t < t1 ) {
                double tt = t - t_0;
                xp[i] = x_0[i] + 0.5 * ddx_r[i] * tt * tt;
            } else if (t < t2 ) {
                double tt = t - t1;
                xp[i] = x_0[0] + 0.5*dx_r[i]*t1 + dx_r[i]*tt;
            } else if (t < t_n ) {
                double tt = t_n - t;
                xp[i] = x_n[i] - .5*ddx_r[i]*tt*tt;
            } else {
                assert(0);
            }
            assert( aa_isfok(xp[i]) );
        }
    }

    // convert back to vector+quaternion form
    aa_fcpy( x, xp, 3 );
    aa_tf_rotvec2quat(xp+3, r);
    return 0;
}

int TrapvelWS::get_dx( double t, double dx[6]) {
    if( t < this->t_0 || t > this->t_n ) {
        aa_fset( dx, 0, 6 );
    } else if( t <= this->t_0 + this->tb ) {
        // accelerating blend
        aa_la_smul( 6, t - this->t_0, ddx_r, dx );
    } else if (t <= this->t_n - this->tb ) {
        // constant velocity
        aa_fcpy( dx, dx_r, 6 );
    } else if (t <= this->t_n ) {
        // deccelerating blend
        double t2 = this->t_n - this->tb;
        aa_la_smul( 6, -(t-t2), ddx_r, dx ); // amount of decellertion
        aa_la_vinc(6, dx_r, dx );  // steady state velocity
    } else {
        // bogus
        assert(0);
    }
    return 0;
}

int TrapvelWS::get_ddx( double t, double ddx[6]) {
    if( t < this->t_0 || t > this->t_n ) {
        aa_fset( ddx, 0, 6 );
    } else if (t <= this->t_0 + this->tb ) {
        // accelerating blend
        aa_fcpy( ddx, this->ddx_r, 6 );
    } else if (t <= this->t_n - this->tb ) {
        // constant velocity
        aa_fset( ddx, 0, 6 );
    } else if (t <= this->t_n ) {
        // deccelerating blend
        aa_la_smul(6, -1, this->ddx_r, ddx );
    } else {
        // bogus
        assert(0);
    }
    return 0;
}

int TrapvelWS::add( double t, const double x[3], const double r[4]) {
    double *xp;
    if( aa_feq(t,0,0) ) {
        xp = this->x_0;
        this->t_0 = t;
    } else if (t > 0) {
        xp = this->x_n;
        this->t_n = t;
    } else {
        return 1;
    }
    aa_fcpy( xp, x, 3 );
    aa_tf_quat2rotvec( r, xp+3 );

    assert( aa_isfok(0) );
    for( size_t i = 0; i < 6; i++ )
        assert( aa_isfok( xp[i] ) );
    //aa_dump_vec( stderr, xp, 6 );
    return 0;
}
