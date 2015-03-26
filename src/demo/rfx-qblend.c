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
// TODO: Is this right?
const char *opt_type = "splend";
// ?? Include velocity in output flag
int   opt_vel = 0;
int opt_verbosity = 0;
// TODO: We probably don't need this
double opt_dt = 0.01;
size_t n_q = 7; // default number of joints

static void
read_points( struct rfx_trajq_points*, FILE *);


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
    for( int c; -1 != (c = getopt(argc, argv, "q:t:o:dv?")); ) {
        switch(c) {
        case 'q':
            n_q = (size_t)atoi(optarg);
            if( n_q <= 0 ) {
                fprintf(stderr, "Invalid number of joints: `%lu'\n", n_q );
                exit(EXIT_FAILURE);
            }
            break;
        case 'v':
            opt_verbosity++;
            break;
        case 'd':
            opt_vel = 1;
            break;
        // ?? Not quite sure what t is at the moment. Might need to delete later
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
                  "  -q NUMJOINTS                 number of joints (default 7)\n"
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
    struct aa_mem_region reg;
    aa_mem_region_init( &reg, 1024 * 64 ); // ?? Why allocate by this amount?

    /* Read */
    struct rfx_trajq_points *points = rfx_trajq_points_alloc( &reg, n_q );
    read_points(points,in);

    /* Generate */
    double v_max = 1;
    double a_max = 0.010;
    struct rfx_trajq_seg_list *segs = rfx_trajq_gen_pblend_max( &reg, points, v_max, a_max );

    /* Plot */
    rfx_trajq_seg_plot(segs, 0.01);

    aa_mem_region_destroy(&reg);
}

static void
read_points( struct rfx_trajq_points *points, FILE *in)
{
    if( opt_verbosity ) fprintf(stderr, "WAYPOINTS\n");

    struct aa_mem_region reg;
    aa_mem_region_init(&reg, 1);

    while( !feof(in) ) {
        aa_mem_region_release(&reg);

        char *line = aa_io_getline(in, &reg);
        if( NULL == line ) continue;

        char *line2 = aa_io_skipblank(line);
        if('\0' == *line2 || AA_IO_ISCOMMENT(*line2)) continue;

        // TODO: Automatically determine array size
        double x[n_q];
        size_t i = aa_io_parsevector( line2, n_q, x, 1, NULL );

        if( n_q != i ) {
            fprintf(stderr, "Invalid line: `%s'\n", line );
            exit(EXIT_FAILURE);
        }

        if( opt_verbosity ) {
            aa_dump_vec( stderr, x, n_q );
        }

        rfx_trajq_points_add( points, 0.0, x );
    }


    aa_mem_region_destroy(&reg);
}
