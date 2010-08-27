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

#ifndef REFLEX_HPP
#define REFLEX_HPP

/** \file reflex.hpp */
#include "reflex.h"

namespace reflex {
    class Trajectory {
    public:
        Trajectory();
        virtual ~Trajectory();

        virtual int validate() = 0;
        virtual int generate() = 0;
    protected:
        double t_0, t_n;
    };

    class WorkspaceTrajectory : public Trajectory {
    public:
        WorkspaceTrajectory();
        virtual ~WorkspaceTrajectory();
        virtual int get_x( double t, double x[3], double r[4] ) = 0;
        virtual int get_dx( double t, double dx[6]) = 0;
        virtual int get_ddx( double t, double ddx[6]) = 0;
        virtual int add(double t, const double x[3], const double r[4] ) = 0;
        void plot(double dt);
        void jacobian_plot(double dt,
                           size_t n, const double *q0,
                           void fkfun(const double *q, double R[9], double v[3] ),
                           void jfun(const double *q, double *J),
                           double k_dls
            );
        void wsctrl_plot(double dt,
                         size_t n, const double *q0,
                         void fkfun(const double *q, double R[9], double v[3] ),
                         void jfun(const double *q, double *J),
                         double k_dls, double k_p[6]
            );
    };

    class TrapvelWS : public WorkspaceTrajectory {
    public:
        TrapvelWS();
        virtual ~TrapvelWS();
        virtual int validate();
        virtual int generate();
        virtual int get_x( double t, double x[3], double r[4] );
        virtual int get_dx( double t, double dx[6] );
        virtual int get_ddx( double t, double ddx[6] );
        virtual int add(double t, const double x[3], const double r[4] );
    protected:
        double x_0[6], x_n[6];
        double dx_max[6], ddx_max[6];
        double tb;
        double dx_r[6], ddx_r[6];
    };
}

#endif //REFLEX_H
