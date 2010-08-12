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
    double X[n][3];
    double T[n];
    size_t i;

    // compute positions
    {
        double t;
        for( i = 0, t = this->t_0; t< this->t_n && i < n; i++, t+=dt ) {
            T[i] = t;
            double r[4];
            this->get_x(t, X[i], r );
        }
    }

    { // plot workspace position
        int r;
        FILE *g = popen("gnuplot -persist", "w");

        assert(g);

        fprintf(g, "set xlabel 'time (s)'\n");
        fprintf(g, "set ylabel 'pos (m)'\n");
        fprintf(g, "set title 'Workspace Positions'\n");
        fprintf(g, "plot '-' with lines title 'WS X'");
        fprintf(g, ", '-' with lines title 'WS Y'");
        fprintf(g, ", '-' with lines title 'WS Z'");
        fprintf(g, "\n");
        // print X
        for(i = 0; i < n; i++ ) {
            fprintf(g, "%f %f\n", T[i], X[i][0] );
        }
        fprintf(g, "e\n");
        // print Y
        for(i = 0; i < n; i++ ) {
            fprintf(g, "%f %f\n", T[i], X[i][1] );
        }
        fprintf(g, "e\n");
        // print Z
        for(i = 0; i < n; i++ ) {
            fprintf(g, "%f %f\n", T[i], X[i][2] );
        }
        fprintf(g, "e\n");
        fflush(g);
        r = pclose(g);
        assert( -1 != r );
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
        for(i = 0; i < n; i++ ) {
            fprintf(g, "%f %f %f\n",
                    X[i][0], X[i][1], X[i][2] );
        }
        fprintf(g, "e\n");
        fflush(g);
        r = pclose(g);
        assert( -1 != r );
    }
}

TrapvelWSTrajectory::TrapvelWSTrajectory()
{
    t_0 = 0;
    t_n = 0;
    aa_fset(x_0,0,3);
    aa_fset(x_n,0,3);
    aa_fset(r_0,0,4);
    aa_fset(r_n,0,4);

    aa_fset(dx_max,1.0,3);
    aa_fset(ddx_max,1.0,3);
}

TrapvelWSTrajectory::~TrapvelWSTrajectory()
{}

int TrapvelWSTrajectory::validate() {
    return 0;
}


int TrapvelWSTrajectory::generate() {
    int i = trapvel_generate( 3, this->t_n,
                              this->x_0, this->x_n,
                              this->dx_max, this->ddx_max,
                              &this->tb, this->dx_r, this->ddx_r );
    fprintf(stderr, "tb: %f\n", tb );
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

    double x[n];
    aa_la_vsub(n, x_f, x_i, x);

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
        fprintf(stderr, "tri\n");
        *ptb = tb;
        return 0;
    }
    fprintf(stderr, "trap\n");

    // needs to be trapezoid
    // find longest acceptable blend time
    for( size_t i = 0; i < n; i++ ) {
        double t = t3 - x[i]/dx_max[i];
        tb = AA_MIN(tb,t);
    }
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

int TrapvelWSTrajectory::get_x( double t, double x[3], double r[4]) {
    // before t0
    if( t < this->t_0 ) {
        aa_fcpy(x, this->x_0, 3 );
        return 0;
    }
    // after t1
    if( t > this->t_n ) {
        aa_fcpy(x, this->x_n, 3 );
        return 0;
    }
    double t1 = tb, t2 = t_n-tb;
    //normal
    for( size_t i = 0; i < 3; i ++ ) {
        if( t < t1 ) {
            double tt = t - t_0;
            x[i] = x_0[i] + 0.5 * ddx_r[i] * tt * tt;
        } else if (t < t2 ) {
            double tt = t - t1;
            x[i] = x_0[0] + 0.5*dx_r[i]*t1 + dx_r[i]*tt;
        } else if (t < t_n ) {
            double tt = t_n - t;
            x[i] = x_n[i] - .5*ddx_r[i]*tt*tt;
        } else {
            assert(0);
        }
    }
    aa_fset( r, 0, 4 );
    return 0;
}

int TrapvelWSTrajectory::get_dx( double t, double dx[6]) {
    (void)t;
    aa_fset( dx, 0, 6 );
    return 0;
}

int TrapvelWSTrajectory::get_ddx( double t, double ddx[6]) {
    (void)t;
    aa_fset( ddx, 0, 6 );
    return 0;
}

int TrapvelWSTrajectory::add( double t, const double x[3], const double r[4]) {
    if( aa_feq(t,0,0) ) {
        aa_fcpy( this->x_0, x, 3 );
        aa_fcpy( this->r_0, r, 3 );
    } else if (t > 0) {
        this->t_n = t;
        aa_fcpy( this->x_n, x, 3 );
        aa_fcpy( this->r_n, r, 3 );
    } else {
        return 1;
    }
    return 0;
}
