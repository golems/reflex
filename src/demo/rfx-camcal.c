/* -*- mode: C++; c-basic-offset: 4 -*- */
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

#include <amino.h>
#include <getopt.h>
#include "reflex.h"


const char *opt_file_cam = NULL;
const char *opt_file_fk = NULL;
const char *opt_file_out = NULL;
size_t opt_test = 0;
int opt_verbosity = 0;
double opt_d_theta = 0;
double opt_d_x = 0;

double opt_zmax_theta = 1;
double opt_zmax_x = 1;


static ssize_t
read_tfs( const char *name, double **A );

static void
write_tfs(const char *comment, const char *name, size_t n, double *A );


static void
write_tf(FILE *tf, const char *comment, double *E );


static FILE*
tryfopen(const char *name );

// static size_t
// reject( size_t n_in, double zmax_theta, double zmax_x, double *EE, const double *E_mean );

void gentest( void );

FILE *global_output = NULL;

// static double qangle( const double *q ) {
//     double aa[4], qr[4];
//     AA_MEM_CPY( qr, q, 4 );
//     aa_tf_qminimize(qr);
//     aa_tf_quat2axang( qr, aa );
//     return aa[3];
// }

// static double relangle( const double *q1, const double *q2 ) {
//     double qr[4];
//     aa_tf_qmulc( q1, q2, qr );
//     return qangle( qr );
// }

// static void tf_dist( const double *E1, const double *E2, double *dtheta, double *dx ) {
//     *dtheta = relangle( E1, E2 );
//     *dx = sqrt(aa_la_ssd( 3, E1+4, E2+4));
// }

//static void eavg( size_t n, const double *EE, double *E_avg );

// static void eavg( size_t n, const double *EE, double *E_avg )
// {

//     double *w = AA_MEM_REGION_LOCAL_NEW_N( double, n );
//     for( size_t i = 0; i < n; i++ ) w[i] = 1/(double)n;

//     aa_tf_qutr_wavg( n, w, EE, 7, E_avg );
//     aa_tf_qminimize(E_avg);


//     if( opt_verbosity ) {

//         for( size_t i = 0; i < n; i ++ ) {
//             const double *e = EE + 7*i;
//             double angle, dist;
//             tf_dist( e, E_avg, &angle, &dist );
//             printf("rel. tf %lu (dp=%f,dx=%f): ", i, angle, dist);
//             aa_dump_vec( stdout, e, 7 );
//         }
//     }
//     aa_mem_region_local_pop(w);
// }


static void iterate( size_t n,
                     double *E_fk, double *E_cam, double *E_avg ) {
    double *w = AA_MEM_REGION_LOCAL_NEW_N( double, n );
    double *Q = AA_MEM_REGION_LOCAL_NEW_N( double, 4*n );

    for( size_t i = 0; i < n; i ++ ) {
        size_t j = 7*i;
        aa_tf_qmulc( E_fk+j, E_cam+j, Q+4*i );
        aa_tf_qminimize(Q+4*i);
        w[i] = 1/(double)n;
    }

    aa_tf_quat_davenport( n, w, Q, 4, E_avg );

    double R[9];
    aa_tf_quat2rotmat(E_avg, R );
    aa_tf_relx_mean( n, R, E_cam+4, 7, E_fk+4, 7, E_avg+4 );
    write_tf( global_output, "Davenport, mean translation", E_avg );

    aa_tf_relx_median( n, R, E_cam+4, 7, E_fk+4, 7, E_avg+4 );
    write_tf( global_output, "Davenport, median translation", E_avg );

    aa_mem_region_local_pop(w);

        //}
}

int main( int argc, char **argv )
{
    /* Parse */
    for( int c; -1 != (c = getopt(argc, argv, "c:k:o:t:d:x:v?")); ) {
        switch(c) {
        case 'v':
            opt_verbosity++;
            break;
        case 'c':
            opt_file_cam = optarg;
            break;
        case 'k':
            opt_file_fk = optarg;
            break;
        case 'o':
            opt_file_out = optarg;
            break;
        case 't':
            opt_test = (size_t)atoi(optarg);
            break;
        case 'd':
            opt_d_theta = atof(optarg)*M_PI/180;
            break;
        case 'x':
            opt_d_x = atof(optarg);
            break;
        case '?':   /* help     */
            puts( "Usage: rfx-camcal -k FK_POSE_FILE -c CAM_POSE_FILE \n"
                  "Calibrate a camera from list of kinematics and camera transforms"
                  "\n"
                  "Options:\n"
                  "  -k FILENAME,                         Forward Kinematics Pose File\n"
                  "  -c FILENAME,                         Camera Marker Pose File\n"
                  "  -o FILENAME,                         Output file\n"
                  "  -t POSE_COUNT,                       Generate test data\n"
                  "  -d DEGREES,                          Corrupt test data rotation by max DEGREES\n"
                  "  -x VALUE,                            Corrupt test translation by max VALUE\n"
                  "  -v,                                  Be verbose\n"
                  "  -?,                                  Program help text\n"
                  "\n"
                  "Files:\n"
                  "  transforms.dat                       Input files are lists transforms, one per line,\n"
                  "                                       given as a quaternion followed by a translation.\n"
                  "                                       Values are separated by spaces.\n"
                  "                                       The quaternion is in xyzw order.\n"
                  "\n"
                  "Examples:\n"
                  "  rfx-camcal -k fk-tf.dat -c cam-tf.dat     Compute average relative transform\n"
                  "\n"
                  "  rfx-camcal -k fk.dat -c cam.dat -t 10     Generate random test poses\n"
                  "\n"
                  "Report bugs to <ntd@gatech.edu>"
                );
            exit(EXIT_SUCCESS);
            break;
        default:
            printf("Unknown argument: `%s'\n", optarg);
            exit(EXIT_FAILURE);
        }
    }

    /* Maybe generate test data */
    if( opt_test ) {
        gentest();
        return 0;
    }

    /* Open */
    if( !opt_file_cam ) {
        fprintf(stderr, "Need to specify camera transform file\n");
        exit(EXIT_FAILURE);
    }
    if( !opt_file_fk ) {
        fprintf(stderr, "Need to specify kinematics transform file\n");
        exit(EXIT_FAILURE);
    }

    global_output = tryfopen( opt_file_out );


    /* Read points */
    double *E_cam=NULL, *E_fk=NULL;
    ssize_t lines_cam = read_tfs( opt_file_cam, &E_cam );
    ssize_t lines_fk = read_tfs( opt_file_fk, &E_fk );
    if( lines_cam != lines_fk ) {
        fprintf(stderr, "Differing line count between `%s' and `%s'\n",
                opt_file_cam, opt_file_fk );
        exit(EXIT_FAILURE);
    }


    /* Compute Rels */
    size_t count = (size_t)lines_cam;

    /* Umeyama */
    double tf[12], EU[7];
    rfx_tf_umeyama( count, E_cam+4, 7, E_fk+4, 7, tf );
    assert(aa_tf_isrotmat(tf));
    aa_tf_tfmat2qutr( tf, EU );
    aa_tf_qminimize(EU);
    {
        write_tf( global_output, "Umeyama, mean translation", EU );
        aa_tf_relx_median( count, tf, E_cam+4, 7, E_fk+4, 7, EU+4 );
        write_tf( global_output, "Umeyama, median translation", EU );
    }

    //write_tfs( "Umeyama Registration", opt_file_out, 1, EU );

    /* Quaternion Average */
    double E_avg[7];
    iterate( count, E_fk, E_cam, E_avg );

    //write_tfs( "Registration", opt_file_out, k, E_avg );
}

static ssize_t
read_tfs( const char *name, double **A )
{

    FILE *f = fopen( name, "r" );
    if( NULL == f ) {
        fprintf(stderr, "Could not open `%s'\n", name );
        exit(EXIT_FAILURE);
    }
    size_t elts = 0;
    ssize_t lines = aa_io_fread_matrix_heap( f, 7, A, &elts );
    if( lines < 0 ) {
        fprintf(stderr, "Error in file `%s' on line %ld.\n", name, lines );
        exit(EXIT_FAILURE);
    }

    return lines;
}

void gentest( void )
{
    if( opt_verbosity) {
        fprintf(stderr, "Generating test data.\n");
    }

    srand((unsigned int)time(NULL)); // might break in 2038

    double E_true[7];
    double *E_cam = (double*)malloc( sizeof(double) * 7 * opt_test );
    double *E_fk = (double*)malloc( sizeof(double) * 7 * opt_test );

    aa_tf_qutr_rand( E_true );
    aa_tf_qminimize(E_true);

    for( size_t i = 0; i < opt_test; i ++ ) {
        size_t j = 7*i;
        aa_tf_qutr_rand( E_fk+j );
        aa_tf_qutr_cmul( E_true, E_fk+j, E_cam+j );

        if( 0 < opt_d_theta || 0 < opt_d_x ) {
            double tmp[7];
            rfx_tf_corrupt( opt_d_theta, opt_d_x, E_cam+j, tmp );
            memcpy( E_cam+j, tmp, 7*sizeof(tmp[0]) );
        }
   }

    write_tfs( "Camera Pose Estimates", opt_file_cam, opt_test, E_cam );
    write_tfs( "Forward Kinematics Poses", opt_file_fk, opt_test, E_fk );
    write_tfs( "True Registration", opt_file_out, 1, E_true );
}


static FILE*
tryfopen(const char *name)
{
    FILE *f = NULL;
    if( NULL == name || 0 == strcmp(name,"-") ) {
        f = stdout;
    } else {
        f = fopen(name, "w");
    }

    if( NULL == f ) {
        fprintf(stderr, "Could not open `%s'\n", name);
        exit(EXIT_FAILURE);
    }
    return f;
}

static void
write_tfs(const char *comment, const char *name, size_t n, double *A )
{
    if( opt_verbosity) {
        fprintf(stderr, "Writing output to `%s'.\n",
                name ? name : "STDOUT");
    }

    FILE *f = tryfopen( name );

    fprintf(f,"# %s\n"
            "\n"
            "# quat_x quat_y quat_z quat_w trans_x trans_y trans_z\n",
            comment);

    for( size_t i = 0; i < n; i ++ ) {
        double *E = A+7*i;
        write_tf( f, NULL, E );
    }
}


static void
write_tf(FILE *f, const char *comment, double *E )
{
    if( comment ) fprintf(f, "# %s\n", comment );
    fprintf(f, "%f %f %f %f %f %f %f\n",
            E[0], E[1], E[2], E[3],
            E[4], E[5], E[6] );
}

// static void
// tf_std( size_t n, const double *EE_in, const double *E_mean,
//         double *dtheta, double *dx,
//         double *theta_std, double *x_std )
// {
//     *theta_std = 0;
//     *x_std = 0;
//     for( size_t i = 0; i < n; i++ ) {
//         size_t j = 7*i;
//         // angle
//         dtheta[i] = relangle( &EE_in[j], E_mean );
//         *theta_std += dtheta[i]*dtheta[i];
//         // translation
//         double x2 = aa_la_ssd( 3, &EE_in[j+3], E_mean+3 );
//         *x_std += x2;
//         dx[i] = sqrt(x2);
//     }

//     *theta_std = sqrt( *theta_std / (double)n );
//     *x_std = sqrt( *x_std / (double)n );
// }


// static size_t
// reject( size_t n_in, double zmax_theta, double zmax_x,
//         double *EE, const double *E_mean )
// {
//     // compute variance
//     double *dx = AA_MEM_REGION_LOCAL_NEW_N( double, n_in );
//     double *dtheta = AA_MEM_REGION_LOCAL_NEW_N( double, n_in );

//     double theta_std, x_std, x_mean_dist = 0;
//     tf_std( n_in, EE, E_mean, dtheta, dx, &theta_std, &x_std );
//     printf( "std_theta: %f\n"
//             "std_x:     %f\n",
//             theta_std, x_std );

//     for( size_t i = 0; i < n_in; i ++ ) x_mean_dist += dx[i] / (double)n_in;
//     printf("x mean dist: %f\n", x_mean_dist);


//     // do it
//     size_t i_out = 0;
//     for( size_t i = 0; i < n_in; i ++ ) {
//         if( dtheta[i] / theta_std < zmax_theta ||
//             dx[i] / x_std < zmax_x )
//         {
//             if( i != i_out ) AA_MEM_CPY( &EE[7*i_out], &EE[7*i], 7 );
//             i_out++;
//         }
//     }

//     return i_out;
// }
