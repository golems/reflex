/*
 * Copyright (c) 2013, Georgia Tech Research Corporation
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

package org.golems.reflex;

import org.golems.amino.Lib;

class Test
{
    public static void main( String[] args ) {
        System.out.println("Run Java Test");


        // Get memory region
        long region = org.golems.amino.Lib.mem_region_create(4096);

        // Create points list
        long points = org.golems.reflex.Lib.trajx_point_list_alloc(region);

        // add points
        double q0[] = new double[]{0,0,0,1};
        double x0[] = new double[]{0,0,0};

        double q1[] = new double[]{0,0,1,0};
        double x1[] = new double[]{0,0,1};

        org.golems.reflex.Lib.trajx_point_list_addb_qv( points, 0, 1,
                                                        q0, x0 );
        org.golems.reflex.Lib.trajx_point_list_addb_qv( points, 10, 1,
                                                        q1, x1 );

        // generate trajectory
        long segments = org.golems.reflex.Lib.trajx_splend_generate( points, region );

        // plot trajectory
        org.golems.reflex.Lib.trajx_seg_list_plot( segments, .01 );

        // print trajectory
        double q[] = new double[4];
        double v[] = new double[3];
        double dx[] = new double[6];
        double ddx[] = new double[6];
        for( double t = 0; t < 10; t += .05 ) {

            //org.golems.reflex.Lib.trajx_seg_list_get_x_qv( segments, t, q, v );
            org.golems.reflex.Lib.trajx_seg_list_get_dx_qv( segments, t, q, v, dx );

            System.out.format("%f : [%f, %f, %f, %f] | [%f, %f, %f]\n",
                              t, q[0], q[1], q[2], q[3],
                              v[0], v[1], v[2]);


        }

        // release memory when finished
        org.golems.amino.Lib.mem_region_release(region);

        // Destroy memory region
        org.golems.amino.Lib.mem_region_destroy(region);


    }

}

