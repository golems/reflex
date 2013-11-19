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

/** Wrapper class for libreflex */
public class Lib
{

    static{
        System.loadLibrary("reflexj");
    }

    /* Native Function Wrappers */

    /** Allocate a point list out of mem_region.
     *
     * All points added to list will be allocated out of the provided
     * memory region.
     *
     * @param mem_region Handle to an org.golems.amino memory region
     *
     * @return A handle to the point list
     */
    public static native long
    trajx_point_list_alloc( long mem_region );

    /** Add blend point to the list.
     *
     * @param point_list Handle to point list from trajx_point_list_alloc()
     * @param t Time to reach this point
     * @param tb Blend time for this point
     * @param q Unit quaternion orientation, in xyzw order
     * @param v Translation, xyz
     */
    public static native int
    trajx_point_list_addb_qv( long point_list, double t, double tb,
                              double[] q, double[] v );

    /** Generate spherical blend trajectory.
     *
     * @param point_list List of way points in the trajectory.
     * @param mem_region Handle to an org.golems.amino memory region
     * to allocation the trajectory segments out of.
     *
     * @return Handle to the trajectory segment list.
     */
    public static native long
    trajx_splend_generate( long point_list, long mem_region );


    /** Get trajectory pose at given time.
     *
     * @param seg_list Handle to trajectory segments.
     * @param q Output for orientation is unit quaternion, xyzw order
     * @param v Output for translation as array xyz
     */
    public static native int
    trajx_seg_list_get_x_qv( long seg_list, double t,
                             double[] q, double[] v );

    /** Get trajectory pose and velocity at given time.
     *
     * @param seg_list Handle to trajectory segments.
     * @param q Output for orientation is unit quaternion, xyzw order
     * @param v Output for translation as array xyz
     * @param dx Output for velocity, xyz translation followed by xyz rotation
     */
    public static native int
    trajx_seg_list_get_dx_qv( long seg_list, double t,
                              double[] q, double[] v,
                              double[] dx);

    /** Get trajectory pose, velocity, and acceleration at given time.
     *
     * @param seg_list Handle to trajectory segments.
     * @param q Output for orientation is unit quaternion, xyzw order
     * @param v Output for translation as array xyz
     * @param dx Output for velocity, xyz translation followed by xyz rotation
     * @param ddx Output for orientation, xyz translation followed by xyz rotation
     */
    public static native int
    trajx_seg_list_get_ddx_qv( long seg_list, double t,
                               double[] q, double[] v,
                               double[] dx, double ddx[] );


    /** Plot the trajectory */
    public static native int
    trajx_seg_list_plot( long seg_list, double dt );

}
