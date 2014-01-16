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


static ssize_t
read_tfs( const char *name, double **A );

static void
write_tfs(const char *comment, const char *name, size_t n, double *A );

void gentest( void );

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
    double *E_rel = (double*)malloc( sizeof(double) * 7 * count );
    double *w = (double*)malloc( sizeof(double) * count );
    for( size_t i = 0; i < count; i ++ ) {
        size_t j = 7*i;
        aa_tf_qutr_mulc( E_cam+j, E_fk+j, E_rel+j );
        w[i] = 1.0 / (double)count;
    }

    /* Quaternion Average */
    double E_avg[7];
    aa_tf_qutr_wavg( count, w, E_rel, 7, E_avg );
    aa_tf_qminimize(E_avg);

    /* print it */
    write_tfs("Average Relative Transform", opt_file_out, 1, E_avg);
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
        aa_tf_qutr_mul( E_true, E_fk+j, E_cam+j );

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

static void
write_tfs(const char *comment, const char *name, size_t n, double *A )
{
    if( opt_verbosity) {
        fprintf(stderr, "Writing output to `%s'.\n",
                name ? name : "STDOUT");
    }

    FILE *f;
    if( NULL == name || 0 == strcmp(name,"-") ) {
        f = stdout;
    } else {
        f = fopen(name, "w");
        if( NULL == f ) {
            fprintf(stderr, "Could not open `%s'\n", name);
            exit(EXIT_FAILURE);
        }
    }

    fprintf(f,"# %s\n"
            "\n"
            "# quat_x quat_y quat_z quat_w trans_x trans_y trans_z\n",
            comment);

    for( size_t i = 0; i < n; i ++ ) {
        double *E = A+7*i;
        fprintf(f, "%f %f %f %f %f %f %f\n",
                E[0], E[1], E[2], E[3],
                E[4], E[5], E[6] );
    }
}
