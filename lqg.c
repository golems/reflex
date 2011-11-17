/* -*- mode: C; c-basic-offset: 4 -*- */
/* ex: set shiftwidth=4 tabstop=4 expandtab: */
/*
 * Copyright (c) 2011, Georgia Tech Research Corporation
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
#include <cblas.h>
#include "reflex.h"



AA_API void rfx_lqg_init( rfx_lqg_t *lqg, size_t n_x, size_t n_u, size_t n_z ) {
    memset( lqg, 0, sizeof(lqg) );
    lqg->n_x = n_x;
    lqg->n_u = n_u;
    lqg->n_z = n_z;

    lqg->x = AA_NEW0_AR( double, lqg->n_x );
    lqg->u = AA_NEW0_AR( double, lqg->n_u );
    lqg->z = AA_NEW0_AR( double, lqg->n_z );

    lqg->A = AA_NEW0_AR( double, lqg->n_x * lqg->n_x );
    lqg->B = AA_NEW0_AR( double, lqg->n_x * lqg->n_u );
    lqg->C = AA_NEW0_AR( double, lqg->n_z * lqg->n_x );

    lqg->P = AA_NEW0_AR( double, lqg->n_x * lqg->n_x );
    lqg->V = AA_NEW0_AR( double, lqg->n_x * lqg->n_x );
    lqg->W = AA_NEW0_AR( double, lqg->n_z * lqg->n_z );

    lqg->Q = AA_NEW0_AR( double, lqg->n_x * lqg->n_x );
    lqg->R = AA_NEW0_AR( double, lqg->n_u * lqg->n_u );

    lqg->K = AA_NEW0_AR( double, lqg->n_u * lqg->n_x );
    lqg->L = AA_NEW0_AR( double, lqg->n_x * lqg->n_z );

}
AA_API void rfx_lqg_destroy( rfx_lqg_t *lqg ) {

    free(lqg->x);
    free(lqg->u);
    free(lqg->z);

    free(lqg->A);
    free(lqg->B);
    free(lqg->C);

    free(lqg->Q);
    free(lqg->R);

    free(lqg->L);
    free(lqg->K);


}



/*
 * X = op(A)*op(B)*op(C) + beta*X
 *
 * m*n = (m*k) * (k*n)
 *
 * X:     m*n
 * op(A): m*k
 * op(B): k*p
 * op(C): p*n
 */
static inline void matmul3( size_t m, size_t n, size_t k, size_t p,
                            enum CBLAS_ORDER transA,
                            enum CBLAS_ORDER transB,
                            enum CBLAS_ORDER transC,
                            const double *A, size_t lda,
                            const double *B, size_t ldb,
                            const double *C, size_t ldc,
                            double beta, double *X, size_t ldx)
{
    int mi = (int)m;
    int ni = (int)n;
    int ki = (int)k;
    int pi = (int)p;
    double T[m*p];
    // T := A*B
    cblas_dgemm( CblasColMajor, transA, transB,
                 mi, pi, ki,
                 1.0, A, (int)lda,
                 B, (int)ldb,
                 0.0, T, mi );

    // X := (A*B) * C + beta*X
    cblas_dgemm( CblasColMajor, CblasNoTrans, transC,
                 mi, ni, pi,
                 1.0, T, mi,
                 C, (int)ldc,
                 beta, X, (int)ldx );
}

AA_API void rfx_lqg_kf_predict( rfx_lqg_t *lqg ) {
    int m = (int)lqg->n_x;
    // x = A*x + B*u
    {
        double x[lqg->n_x];
        aa_lsim_dstep( lqg->n_x, lqg->n_u,
                       lqg->A, lqg->B,
                       lqg->x, lqg->u,
                       x );
        memcpy(lqg->x, x, sizeof(x));
    }
    // P = A * P * A**T + V
    {
        double T[lqg->n_x * lqg->n_x];
        // T := A*P
        cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
                     m, m, m,
                     1.0, lqg->A, m,
                     lqg->P, m,
                     0.0, T, m );
        // P := (A*P) * A**T + V
        memcpy( lqg->P, lqg->V, sizeof(T) );
        cblas_dgemm( CblasColMajor, CblasNoTrans, CblasTrans,
                     m, m, m,
                     1.0, T, m,
                     lqg->A, m,
                     1.0, lqg->P, m );
    }

}
AA_API void rfx_lqg_kf_correct( rfx_lqg_t *lqg ) {

    // K = P * C**T * (C * P * C**T + W)**-1
    {
        double Kp[lqg->n_z*lqg->n_z];
        // Kp := C * P * C**T + W
        memcpy(Kp, lqg->W, sizeof(Kp) );
        matmul3( lqg->n_z, lqg->n_z, lqg->n_x, lqg->n_x,
                 CblasNoTrans, CblasNoTrans, CblasTrans,
                 lqg->C, lqg->n_z,
                 lqg->P, lqg->n_x,
                 lqg->C, lqg->n_x,
                 1.0, Kp, lqg->n_z );
        aa_la_inv(lqg->n_z, Kp);

        // K := P * C**T * Kp
        matmul3( lqg->n_x, lqg->n_z, lqg->n_x, lqg->n_z,
                 CblasNoTrans, CblasTrans, CblasNoTrans,
                 lqg->P, lqg->n_x,
                 lqg->C, lqg->n_x,
                 Kp, lqg->n_z,
                 0.0, lqg->K, lqg->n_x );
    }
    // x = x + K * (z - C*x)
    {
        // T := z - C*x
        double T[lqg->n_z];
        memcpy(T, lqg->z, sizeof(T));
        cblas_dgemv( CblasColMajor, CblasNoTrans,
                     (int)lqg->n_z, (int)lqg->n_x,
                     -1.0, lqg->C, (int)lqg->n_z,
                     lqg->x, 1,
                     1.0, T, 1 );
        // x = K * (z - C*x) + x
        cblas_dgemv( CblasColMajor, CblasNoTrans,
                     (int)lqg->n_x, (int)lqg->n_z,
                     1.0, lqg->K, (int)lqg->n_x,
                     T, 1,
                     1.0, lqg->x, 1 );
    }

    // P = (I - K*C) * P
    {
        double KC[lqg->n_x * lqg->n_x];
        cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
                     (int)lqg->n_x, (int)lqg->n_x, (int)lqg->n_z,
                     -1.0, lqg->K, (int)lqg->n_x,
                     lqg->C, (int)lqg->n_z,
                     0.0, KC, (int)lqg->n_x );
        // KC += I
        for( size_t i = 0; i < lqg->n_x; i ++ )
            AA_MATREF(KC, lqg->n_x, i, i) += 1;

        double PT[lqg->n_x*lqg->n_x];
        memcpy(PT, lqg->P, sizeof(PT));
        cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
                     (int)lqg->n_x, (int)lqg->n_x, (int)lqg->n_x,
                     1.0, KC, (int)lqg->n_x,
                     PT, (int)lqg->n_x,
                     0.0, lqg->P, (int)lqg->n_x );
    }

}

// kalman-bucy gain
void rfx_lqg_kbf_gain( rfx_lqg_t *lqg ) {

    // dP = A*P + P*A' - P*C'*W^{-1}*C*P + V
    // K = P * C' * W^{-1}

    double Winv[lqg->n_z*lqg->n_z];
    memcpy(Winv, lqg->W, sizeof(Winv));
    aa_la_inv( lqg->n_z, Winv);

    double ric_A[lqg->n_x*lqg->n_x];
    aa_la_transpose2(lqg->n_x, lqg->n_x, lqg->A, ric_A );

    double ric_B[lqg->n_x*lqg->n_x];

    // ric_B := C' * W^{-1} * C
    matmul3( lqg->n_x, lqg->n_x, lqg->n_z, lqg->n_z,
             CblasTrans, CblasNoTrans, CblasNoTrans,
             lqg->C, lqg->n_x,
             Winv, lqg->n_z,
             lqg->C, lqg->n_z,
             0.0, ric_B, lqg->n_x );


    // solve ARE with dP = 0, result is P
    aa_la_care_laub( lqg->n_x, lqg->n_x, lqg->n_x,
                     ric_A, ric_B, lqg->V, lqg->P );

    // K = P * C' * W^{-1}
    matmul3( lqg->n_x, lqg->n_z, lqg->n_x, lqg->n_z,
             CblasNoTrans, CblasTrans, CblasNoTrans,
             lqg->P, lqg->n_x,
             lqg->C, lqg->n_x,
             Winv, lqg->n_z,
             0.0, lqg->K, lqg->n_x );

}

void rfx_lqg_kbf_step1( rfx_lqg_t *lqg, double dt ) {
    // dx = A*x + B*u + K(z-Cx)
    double dx[lqg->n_x];

    // dx := A*x + B*u
    aa_lsim_dstep( lqg->n_x, lqg->n_u,
                   lqg->A, lqg->B,
                   lqg->x, lqg->u,
                   dx );
    // dx := K * (z - C*x) + dx
    {
        // T := z - C*x
        double T[lqg->n_z];
        memcpy(T, lqg->z, sizeof(T));
        cblas_dgemv( CblasColMajor, CblasNoTrans,
                     (int)lqg->n_z, (int)lqg->n_x,
                     -1.0, lqg->C, (int)lqg->n_z,
                     lqg->x, 1,
                     1.0, T, 1 );
        cblas_dgemv( CblasColMajor, CblasNoTrans,
                     (int)lqg->n_x, (int)lqg->n_z,
                     1.0, lqg->K, (int)lqg->n_x,
                     T, 1,
                     1.0, lqg->x, 1 );
    }

    // x := dt*dx + x
    cblas_daxpy( (int)lqg->n_x, dt, dx, 1, lqg->x, 1 );
}

//AA_API void rfx_lqg_observe_euler( rfx_lqg_t *lqg, double dt, aa_region_t *reg ) {
    // --- calculate kalman gain ---


    // --- predict from previous input ---
    // dx = A*x + B*u + K * ( z - C*xh )
    //double *dx = (double*)aa_region_alloc(reg, sizeof(double)*lqg->n_x);
    // zz := z
    //double *zz = (double*)aa_region_alloc(reg, sizeof(double)*lqg->n_z);


    // x = x + dx * dt
    //aa_la_axpy( lqg->n_x, dt, dx, lqg->x );
//}

//AA_API void rfx_lqg_ctrl( rfx_lqg_t *lqg, double dt, aa_region_t *reg ) {
    // compute optimal gain
    //  -dS = A'*S + S*A - S*B*R^{-1}*B' + Q
    // solve ARE, result is S
    // L = R^{-1} * B' * S  : (nu*nu) * (nu*nx) * (nx*nx)

    // compute current input
    // u = -Lx
//}


AA_API void rfx_lqg_observe
( size_t n_x, size_t n_u, size_t n_z,
  const double *A, const double *B, const double *C,
  const double *x, const double *u, const double *z,
  const double *K,
  double *dx, double *zwork )
{
    //zwork := z
    aa_fcpy(zwork, z, n_z);
    // zz := -C*x + zz
    cblas_dgemv( CblasColMajor, CblasNoTrans,
                 (int)n_z, (int)n_x,
                 -1.0, C, (int)n_z,
                 x, 1, 1.0, zwork, 1 );
    // dx := K * zwork
    cblas_dgemv( CblasColMajor, CblasNoTrans,
                 (int)n_x, (int)n_z,
                 1.0, K, (int)n_x,
                 zwork, 1, 0.0, dx, 1 );
    // dx := B*u + dx
    cblas_dgemv( CblasColMajor, CblasNoTrans,
                 (int)n_x, (int)n_u,
                 1.0, B, (int)n_x,
                 u, 1, 1.0, dx, 1 );
    // dx := A*x + dx
    cblas_dgemv( CblasColMajor, CblasNoTrans,
                 (int)n_x, (int)n_x,
                 1.0, A, (int)n_x,
                 x, 1, 1.0, dx, 1 );
}
