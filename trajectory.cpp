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

TrapvelWSTrajectory::TrapvelWSTrajectory():
    t_0(0), t_n(0)
{
    aa_fset(x_0,0,3);
    aa_fset(x_n,0,3);
    aa_fset(r_0,0,4);
    aa_fset(r_n,0,4);
}

TrapvelWSTrajectory::~TrapvelWSTrajectory()
{}

int TrapvelWSTrajectory::validate() {
    return 0;
}

int TrapvelWSTrajectory::generate() {
}

int TrapvelWSTrajectory::get_x( double x[3], double r[4], double t ) {

}

int TrapvelWSTrajectory::get_dx( double dx[6], double t ) {

}

int TrapvelWSTrajectory::get_ddx( double ddx[6], double t ) {

}

int TrapvelWSTrajectory::add(const double x[3], const double r[4], double t) {
    if( 0 == t ) {
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
