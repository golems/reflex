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

#include <amino.h>
#include <getopt.h>
#include "reflex.h"



const char *opt_file_in = NULL;
const char *opt_file_out = NULL;
const char *opt_type = "splend";
int   opt_vel = 0;
int opt_verbosity = 0;
double opt_dt = 0.01;


static void
read_points( struct rfx_trajx_point_list*, FILE *);

static void
write_traj ( struct rfx_trajx_seg_list*, FILE *);


static void posarg( char *arg, int i )
{
    (void) i;
    if( opt_file_in ) {
        fprintf(stderr, "Can only handle one input file\n");
    }
    opt_file_in = arg;
}

int main( int argc, char **argv )
{
    /* Parse */
    int i = 0;
    for( int c; -1 != (c = getopt(argc, argv, "t:o:dv?")); ) {
        switch(c) {
        case 'v':
            opt_verbosity++;
            break;
        case 'd':
            opt_vel = 1;
            break;
        case 't':
            opt_dt = atof(optarg);
            if( opt_dt <= 0 ) {
                fprintf(stderr, "Invalid time step: `%s'\n", optarg );
                exit(EXIT_FAILURE);
            }
            break;
        case 'o':
            opt_file_out = optarg;
            break;
        case '?':   /* help     */
            puts( "Usage: rfx-trajgen [OPTIONS...] input-file\n"
                  "Convert a sequence of waypoints to a dense trajectory\n"
                  "\n"
                  "Options:\n"
                  "  -t TIMESTEP,                 timestep (seconds)\n"
                  "  -d,                          Include velocities in output\n"
                  "  -o FILE,                     output file\n"
                  "  -v,                          Be verbose\n"
                  "  -?,                          Program help text\n"
                  "\n"
                  "\n"
                  "Report bugs to <ntd@gatech.edu>"
                );
            exit(EXIT_SUCCESS);
            break;
        default:
            posarg( optarg, i++ );
        }
    }
    while( optind < argc ) {
        posarg(argv[optind++], i++);
    }

    /* Open */
    FILE *in = opt_file_in ? fopen( opt_file_in, "r" ) : stdin;
    FILE *out = opt_file_out ? fopen( opt_file_out, "w" ) : stdout;
    struct aa_mem_region reg;
    aa_mem_region_init( &reg, 1024 * 64 );

    /* Read */
    struct rfx_trajx_point_list *points = rfx_trajx_point_list_alloc( &reg );
    read_points(points,in);

    /* Generate */
    struct rfx_trajx_seg_list *segs = rfx_trajx_splend_generate( points, &reg );

    /* Write */
    write_traj(segs, out);

}

static void
read_points( struct rfx_trajx_point_list *points, FILE *in)
{
    double t, tb;
    rfx_tf tf;

    if( opt_verbosity ) fprintf(stderr, "WAYPOINTS\n");

    struct aa_mem_region reg;
    aa_mem_region_init(&reg, 1);

    while( !feof(in) ) {
        aa_mem_region_release(&reg);

        char *line = aa_io_getline(in, &reg);
        if( NULL == line ) continue;

        char *line2 = aa_io_skipblank(line);
        if('\0' == *line2 || AA_IO_ISCOMMENT(*line2)) continue;

        double x[9];
        size_t i = aa_io_parsevector( line2, 9, x, 1, NULL );

        if( 9 != i ) {
            fprintf(stderr, "Invalid line: `%s'\n", line );
            exit(EXIT_FAILURE);
        }

        t      = x[0];
        tb     = x[1];
        tf.r.x = x[2];
        tf.r.y = x[3];
        tf.r.z = x[4];
        tf.r.w = x[5];
        tf.v.x = x[6];
        tf.v.y = x[7];
        tf.v.z = x[8];

        if( opt_verbosity ) {
            fprintf(stderr, "%f / %f: ", t, tb);
            aa_dump_vec( stderr, tf.data, 7 );
        }

        // TODO: check times
        rfx_trajx_point_list_addb_qv( points, t, tb,
                                      tf.r.data, tf.v.data );
    }

    aa_mem_region_destroy(&reg);
}

static void
write_traj ( struct rfx_trajx_seg_list *segs, FILE *out)
{
    double ti = rfx_trajx_seg_list_get_t_i( segs );
    double tf = rfx_trajx_seg_list_get_t_f( segs );

    for( double t = ti; t < tf; t+= opt_dt ) {
        rfx_tf_dx tf_dx;
        rfx_trajx_seg_list_get_dx_qv( segs,
                                      t, tf_dx.tf.r.data, tf_dx.tf.v.data,
                                      tf_dx.dx.data );
        if( opt_vel ) {
            fprintf( out, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf  %lf  %lf  %lf  %lf\n",
                     t,
                     tf_dx.tf.r.x,
                     tf_dx.tf.r.y,
                     tf_dx.tf.r.z,
                     tf_dx.tf.r.w,
                     tf_dx.tf.v.x,
                     tf_dx.tf.v.y,
                     tf_dx.tf.v.z,
                     tf_dx.dx.dv[0],
                     tf_dx.dx.dv[1],
                     tf_dx.dx.dv[2],
                     tf_dx.dx.omega[0],
                     tf_dx.dx.omega[1],
                     tf_dx.dx.omega[2]
                );
        } else {
            /* print positions */
            fprintf( out, "%lf %lf %lf %lf %lf %lf %lf %lf\n",
                     t,
                     tf_dx.tf.r.x,
                     tf_dx.tf.r.y,
                     tf_dx.tf.r.z,
                     tf_dx.tf.r.w,
                     tf_dx.tf.v.x,
                     tf_dx.tf.v.y,
                     tf_dx.tf.v.z );
        }
    }
}
