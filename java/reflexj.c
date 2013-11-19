/* -*- mode: C; c-basic-offset: 4 -*- */
/* ex: set shiftwidth=4 tabstop=4 expandtab: */
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


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <amino.h>
#include "reflex.h"
#include "org_golems_reflex_Lib.h"

/*
 * Class:     org_golems_reflex_Lib
 * Method:    trajx_point_list_alloc
 * Signature: (J)J
 */
JNIEXPORT jlong JNICALL Java_org_golems_reflex_Lib_trajx_1point_1list_1alloc
(JNIEnv *env, jclass clazz, jlong region)
{
    (void)env; (void)clazz;
    return (intptr_t)rfx_trajx_point_list_alloc( (struct aa_mem_region*)region );
}

/*
 * Class:     org_golems_reflex_Lib
 * Method:    trajx_point_list_addb_qv
 * Signature: (J[D[D)I
 */
JNIEXPORT jint JNICALL Java_org_golems_reflex_Lib_trajx_1point_1list_1addb_1qv
(JNIEnv *env, jclass clazz, jlong point_list, jdouble t, jdouble tb, jdoubleArray jq, jdoubleArray jv)
{
    (void)clazz;
    // TODO: check sizes
    double *q = (*env)->GetDoubleArrayElements(env,jq,0);
    double *v = (*env)->GetDoubleArrayElements(env,jv,0);

    rfx_trajx_point_list_addb_qv( (struct rfx_trajx_point_list*)point_list,
                                  t, tb, q, v );

    (*env)->ReleaseDoubleArrayElements(env,jq,q,0);
    (*env)->ReleaseDoubleArrayElements(env,jv,v,0);

    return RFX_OK;
}

/*
 * Class:     org_golems_reflex_Lib
 * Method:    trajx_splend_generate
 * Signature: (J)J
 */
JNIEXPORT jlong JNICALL Java_org_golems_reflex_Lib_trajx_1splend_1generate
(JNIEnv *env, jclass clazz, jlong points, jlong region)
{
    (void)env; (void)clazz;

    return (intptr_t)rfx_trajx_splend_generate( (struct rfx_trajx_point_list*)points,
                                                (struct aa_mem_region*)region );
}

/*
 * Class:     org_golems_reflex_Lib
 * Method:    trajx_seg_list_get_x_qv
 * Signature: (JD[D[D)I
 */
JNIEXPORT jint JNICALL Java_org_golems_reflex_Lib_trajx_1seg_1list_1get_1x_1qv
(JNIEnv *env, jclass clazz, jlong seg, jdouble t,
 jdoubleArray jq, jdoubleArray jv)
{
    (void)clazz;
    double q[4], v[3];
    int i = rfx_trajx_seg_list_get_x_qv( (struct rfx_trajx_seg_list*)seg, t,
                                         q, v );
    (*env)->SetDoubleArrayRegion(env, jq, 0, 4, q );
    (*env)->SetDoubleArrayRegion(env, jv, 0, 3, v );

    return i;
}

/*
 * Class:     org_golems_reflex_Lib
 * Method:    trajx_seg_list_get_dx_qv
 * Signature: (JD[D[D[D)I
 */
JNIEXPORT jint JNICALL Java_org_golems_reflex_Lib_trajx_1seg_1list_1get_1dx_1qv
(JNIEnv *env, jclass clazz, jlong seg, jdouble t,
 jdoubleArray jq, jdoubleArray jv, jdoubleArray jdx)
{
    (void)clazz;
    double q[4], v[3], dx[6];
    int i = rfx_trajx_seg_list_get_dx_qv( (struct rfx_trajx_seg_list*)seg, t,
                                          q, v, dx );
    (*env)->SetDoubleArrayRegion(env, jq, 0, 4, q );
    (*env)->SetDoubleArrayRegion(env, jv, 0, 3, v );
    (*env)->SetDoubleArrayRegion(env, jdx, 0, 6, dx );

    return i;
}

/*
 * Class:     org_golems_reflex_Lib
 * Method:    trajx_seg_list_get_ddx_qv
 * Signature: (JD[D[D[D[D)I
 */
JNIEXPORT jint JNICALL Java_org_golems_reflex_Lib_trajx_1seg_1list_1get_1ddx_1qv
(JNIEnv *env, jclass clazz, jlong seg, jdouble t,
 jdoubleArray jq, jdoubleArray jv, jdoubleArray jdx, jdoubleArray jddx)
{
    (void)clazz;
    double q[4], v[3], dx[6], ddx[6];
    int i = rfx_trajx_seg_list_get_ddx_qv( (struct rfx_trajx_seg_list*)seg, t,
                                           q, v, dx, ddx );
    (*env)->SetDoubleArrayRegion(env, jq, 0, 4, q );
    (*env)->SetDoubleArrayRegion(env, jv, 0, 3, v );
    (*env)->SetDoubleArrayRegion(env, jdx, 0, 6, dx );
    (*env)->SetDoubleArrayRegion(env, jddx, 0, 6, ddx );

    return i;
}

/*
 * Class:     org_golems_reflex_Lib
 * Method:    trajx_seg_list_plot
 * Signature: (JD)I
 */
JNIEXPORT jint JNICALL Java_org_golems_reflex_Lib_trajx_1seg_1list_1plot
  (JNIEnv *env, jclass clazz, jlong seg, jdouble dt)
{
    (void)env; (void)clazz;
    rfx_trajx_seg_list_plot( (struct rfx_trajx_seg_list*)seg, dt, NULL );
    return RFX_OK;

}
