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

Trajectory::Trajectory() {}
Trajectory::~Trajectory() {}

WorkspaceTrajectory::WorkspaceTrajectory() {}
WorkspaceTrajectory::~WorkspaceTrajectory() {}

extern "C" int rfx_trapvel_generate( size_t n, double t_f,
                                     const double *x_i, const double *x_f,
                                     const double *dx_max, const double *ddx_max,
                                     double *tb, double *dx_r, double *ddx_r );



TrapvelWS::TrapvelWS()
{
    // t_0 = 0;
    // t_n = 0;
    // AA_ZERO_AR(x_0);
    // AA_ZERO_AR(x_n);

    aa_mem_region_init(&reg, 1024*32);
    rfx_trajq_trapvel_init( &this->traj, &reg, 6 );
    for( size_t i = 0; i < traj.traj.n_q; i ++ ) {
        traj.dq_max[i] = 10.0;
        traj.ddq_max[i] = 10.0;
    }
}

TrapvelWS::~TrapvelWS()
{
    aa_mem_region_destroy(&reg);
}

int TrapvelWS::validate() {
    // sanity checks
    // if( ! aa_feq( this->t_0, 0, 0) ) return -1;
    // if( this->tb <= 0 ) return -1;
    // if( this->t_n < this->tb ) return -1;
    // for( size_t i = 0; i < 6; i ++ ) {
    //     if( ! aa_feq( this->tb * ddx_r[i], dx_r[i], .001 ) ) {
    //         return -1;
    //     }
    // }
    // ok
    return 0;
}


int TrapvelWS::generate() {
    double rv[3];
    // aa_fcpy(rv, traj.traj.Q + 6 + 3, 3);
    // aa_tf_rotvec_near( rv, traj.traj.Q+3, traj.traj.Q+6+3 );
    // int i = rfx_trapvel_generate( 6, this->t_n,
    //                               this->x_0, this->x_n,
    //                               this->dx_max, this->ddx_max,
    //                               &this->tb, this->dx_r, this->ddx_r );
    int i = rfx_trajq_generate(&traj.traj);
    return i;

    //fprintf(stderr, "tb: %f\n", tb );
    //return i;
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


int TrapvelWS::get_x( double t, double x[3], double r[4]) {
    // double *xp, xps[6];  // temp storage
    double  xp[6];  // temp storage
    rfx_trajq_get_q( &traj.traj, t, xp );

    // if( t < this->t_0 ) {
    //     // before t0
    //     xp = this->x_0;
    // }
    // else if( t > this->t_n ) {
    //     // after t1
    //     xp = this->x_n;
    // } else {
    //     //normal
    //     xp = xps;
    //     double t1 = tb, t2 = t_n-tb;
    //     for( size_t i = 0; i < 6; i ++ ) {
    //         if( t < t1 ) {
    //             double tt = t - t_0;
    //             xp[i] = x_0[i] + 0.5 * ddx_r[i] * tt * tt;
    //         } else if (t < t2 ) {
    //             double tt = t - t1;
    //             xp[i] = x_0[0] + 0.5*dx_r[i]*t1 + dx_r[i]*tt;
    //         } else if (t < t_n ) {
    //             double tt = t_n - t;
    //             xp[i] = x_n[i] - .5*ddx_r[i]*tt*tt;
    //         } else {
    //             assert(0);
    //         }
    //         assert( aa_isfok(xp[i]) );
    //     }
    // }

    // convert back to vector+quaternion form
    aa_fcpy( x, xp, 3 );
    aa_tf_rotvec2quat(xp+3, r);
    return 0;
}

int TrapvelWS::get_dx( double t, double dx[6]) {
    return rfx_trajq_get_dq( &traj.traj, t, dx );
    // if( t < this->t_0 || t > this->t_n ) {
    //     aa_fset( dx, 0, 6 );
    // } else if( t <= this->t_0 + this->tb ) {
    //     // accelerating blend
    //     aa_la_smul( 6, t - this->t_0, ddx_r, dx );
    // } else if (t <= this->t_n - this->tb ) {
    //     // constant velocity
    //     aa_fcpy( dx, dx_r, 6 );
    // } else if (t <= this->t_n ) {
    //     // deccelerating blend
    //     double t2 = this->t_n - this->tb;
    //     aa_la_smul( 6, -(t-t2), ddx_r, dx ); // amount of decellertion
    //     aa_la_vinc(6, dx_r, dx );  // steady state velocity
    // } else {
    //     // bogus
    //     assert(0);
    // }
    // return 0;
}

int TrapvelWS::get_ddx( double t, double ddx[6]) {
    return rfx_trajq_get_ddq( &traj.traj, t, ddx );
    // if( t < this->t_0 || t > this->t_n ) {
    //     aa_fset( ddx, 0, 6 );
    // } else if (t <= this->t_0 + this->tb ) {
    //     // accelerating blend
    //     aa_fcpy( ddx, this->ddx_r, 6 );
    // } else if (t <= this->t_n - this->tb ) {
    //     // constant velocity
    //     aa_fset( ddx, 0, 6 );
    // } else if (t <= this->t_n ) {
    //     // deccelerating blend
    //     aa_la_smul(6, -1, this->ddx_r, ddx );
    // } else {
    //     // bogus
    //     assert(0);
    // }
    // return 0;
}

int TrapvelWS::add( double t, const double x[3], const double r[4]) {
    double xp[6];
    size_t i = 0;
    if( aa_feq(t,0,0) ) {
        //xp = this->x_0;
        //this->t_0 = t;
        i = 0;
    } else if (t > 0) {
        i = 1;
        //xp = this->x_n;
        //this->t_n = t;
    } else {
        return 1;
    }
    aa_fcpy( xp, x, 3 );
    aa_tf_quat2rotvec( r, xp+3 );

    assert( aa_isfok(0) );
    for( size_t i = 0; i < 6; i++ )
        assert( aa_isfok( xp[i] ) );
    //aa_dump_vec( stderr, xp, 6 );
    //return 0;
    // rfx_trajq_add( &traj.traj, i, t, xp );
}


/*
ParaBlendWS::ParaBlendWS()
{
    t_0 = 0;
    t_n = 0;
    aa_fzero( this->ddx, 6 );
}

ParaBlendWS::~ParaBlendWS()
{}


int ParaBlendWS::validate() {
    return 0;
}

int ParaBlendWS::generate() {

    for(  std::map< double, T >::iterator itr = this->points.begin();
          itr != this->points.end();
          itr++ ) {
        fprintf(stderr, "%f: ", itr->first );
        aa_dump_vec(stderr, itr->second.x, 6 );
    }
    return 0;
}

int ParaBlendWS::get_x(double t, double x[3], double r[4]) {
    return 0;
}

int ParaBlendWS::get_dx(double t, double x[6]) {
    return 0;
}

int ParaBlendWS::get_ddx(double t, double x[6]) {
    return 0;
}

int ParaBlendWS::add(double t, const double x[3], const double r[4]) {
    double (&X)[6] = this->points[t].x;
    aa_fcpy(X,x,3);
    aa_tf_quat2rotvec( r, X+3 );
    return 0;
}

void ParaBlendWS::set_ddx(double _ddx[6]) {
    aa_fcpy(this->ddx, _ddx, 6 );
}

*/
