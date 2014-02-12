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


const char *opt_file_cam = NULL;
const char *opt_file_fk = NULL;
const char *opt_file_out = NULL;
const char *opt_file_id = NULL;
size_t opt_test = 0;
int opt_verbosity = 0;
double opt_d_theta = 0;
double opt_d_x = 0;

double opt_zmax_theta = 1;
double opt_zmax_x = 1;

static ssize_t
read_mat( const char *name, size_t n, double **A );

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

struct tf_cor {
    double id;
    double X[7];
    double Y[7];
};

static void tf_cor( int opts, size_t n,
                    const double *qx, size_t ldqx,
                    const double *vx, size_t ldvx,
                    const double *qy, size_t ldqy,
                    const double *vy, size_t ldvy,
                    double *Z );

static void tf_cor_em( int opts, size_t k, size_t n,
                       struct tf_cor *cor,
                       double *Z );

#define TF_COR_LD (sizeof(struct tf_cor)/sizeof(double))

static int tf_cor_compar( const void *_a, const void *_b ) {
    struct tf_cor *a = (struct tf_cor*)_a;
    struct tf_cor *b = (struct tf_cor*)_b;
    if( a->id < b->id ) return -1;
    if( a->id > b->id ) return 1;
    return 0;
}

enum tf_cor_opts {
    TF_COR_O_TRANS_MEDIAN = 0x1,
    TF_COR_O_TRANS_MEAN = 0x2,
    TF_COR_O_ROT_UMEYAMA = 0x4,
    TF_COR_O_ROT_DAVENPORT = 0x8,
    TF_COR_O_ROT_MEDIAN = 0x10
};


// Find TF from correspondences
static void tf_cor( int opts, size_t n,
                    const double *qx, size_t ldqx,
                    const double *vx, size_t ldvx,
                    const double *qy, size_t ldqy,
                    const double *vy, size_t ldvy,
                    double *Z )
{

    double *top = AA_MEM_REGION_LOCAL_NEW_N( double, 1 );

    /*-- Orientation --*/
    double q_fit[3][4];
    size_t n_fit = 0;

    if( opts & TF_COR_O_ROT_UMEYAMA ) {
        double tf[12];
        rfx_tf_umeyama( n, vx, ldvx, vy, ldvy, tf );
        assert(aa_tf_isrotmat(tf));
        aa_la_transpose(3,tf);
        aa_tf_rotmat2quat( tf, q_fit[n_fit++] );
    }

    double *Qrel = AA_MEM_REGION_LOCAL_NEW_N( double, 4*n );
    for( size_t i = 0; i < n; i ++ ) {
        double *q = AA_MATCOL(Qrel,4,i);
        aa_tf_qmulc( AA_MATCOL(qx,ldqx,i),
                         AA_MATCOL(qy,ldqy,i),
                     q );
        aa_tf_qminimize(q);
    }

    if( opts & TF_COR_O_ROT_DAVENPORT )
        aa_tf_quat_davenport( n, NULL, Qrel, 4, q_fit[n_fit++] );

    if( opts & TF_COR_O_ROT_MEDIAN )
        rfx_tf_qangmedian ( n, Qrel, 4, q_fit[n_fit++] );

    aa_tf_quat_davenport(n_fit, NULL, q_fit[0], 4, Z);

    /*-- deviation -- */
    /* double *angles = AA_MEM_REGION_LOCAL_NEW_N( double, n ); */
    /* for( size_t i = 0; i < n; i ++ ) { */
    /*     double *q = AA_MATCOL(Qrel,4,i); */
    /*     angles[i] = aa_tf_qangle_rel(Z, q); */
    /* } */
    /* double astd = aa_la_d_vecstd( n, angles, 1, 0 ); */
    /* double amad = aa_la_d_median( n, angles, 1 ); */
    /* double amax = aa_la_max( n, angles ); */

    /* fprintf(global_output, */
    /*         //"# %s\n" */
    /*         "# angle std: %f\n" */
    /*         "# angle mad: %f\n" */
    /*         "# angle max: %f\n", */
    /*         //comment, */
    /*         astd, amad, amax */
    /*     ); */

    /*-- Translation --*/
    double R[9];
    aa_tf_quat2rotmat(Z,R);
    assert( aa_tf_isrotmat(R) );

    // rels
    double *Zt = AA_MEM_REGION_LOCAL_NEW_N(double, 3*n);
    aa_cla_dlacpy( ' ', 3, (int)n,
                   vx, (int)ldvx,
                   Zt, 3 );

    for( size_t j = 0; j < n; j ++ ) {
        double yp[3] = {0};
        aa_tf_9rot( R, AA_MATCOL(vy,ldvy,j), yp );
        for( size_t i = 0; i < 3; i++ )
            AA_MATREF(Zt, 3, i, j) -=  yp[i];
        assert( aa_tf_isrotmat(R) );
    }

    assert( aa_tf_isrotmat(R) );
    if( opts & TF_COR_O_TRANS_MEAN ) {
        aa_la_d_colmean( 3, n, Zt, 3, Z+4 );
        //aa_tf_relx_mean( n, R, vy,ldvy, vx,ldvx, Z+4 );
    } else if (opts & TF_COR_O_TRANS_MEDIAN ) {
        //aa_tf_relx_median( n, R, vy,ldvy, vx,ldvx, Z+4 );
        for( size_t i = 0; i < 3; i ++ )
            (Z+4)[i] = aa_la_d_median( n, Zt+i, 3 );
    }


    aa_mem_region_local_pop(top);
}

static void tf_cor_em( int opts, size_t k, size_t n,
                       struct tf_cor *cor,
                       double *Z )
{
    void *top = aa_mem_region_local_alloc(1);
    // count markers
    size_t m = 1;
    for( size_t i = 1; i < n; i ++ ) {
        if( !aa_feq(cor[i].id, cor[i-1].id, 0) ) {
            m++;
        }
    }

    // find end indices
    size_t i_end[m];
    size_t i_beg[m];
    size_t *i_mark = AA_MEM_REGION_LOCAL_NEW_N( size_t, n );
    i_beg[0] = 0;
    for( size_t i=0, j=0; i < n-1; i ++ ) {
        i_mark[i] = j;
        if( !aa_feq(cor[i].id, cor[i+1].id, 0) ) {
            i_end[j] = i+1;
            j++;
            i_beg[j] = i+1;
        }
    }
    i_end[m-1] = n;
    i_mark[n-1] = m-1;


    // initialize offsets
    rfx_tf *tf_off = AA_MEM_REGION_LOCAL_NEW_N( rfx_tf, m );
    for( size_t i = 0; i < m; i ++ ) {
        AA_MEM_CPY( tf_off[i].r.data, aa_tf_quat_ident, 4 );
        AA_MEM_SET( tf_off[i].v.data, 0, 3 );
    }

    rfx_tf *X_conj = AA_MEM_REGION_LOCAL_NEW_N( rfx_tf, n );
    for( size_t i = 0; i < n; i ++ ) {
        aa_tf_qutr_conj(cor[i].X, X_conj[i].data);
    }

    // iterate
    rfx_tf *X_p = AA_MEM_REGION_LOCAL_NEW_N( rfx_tf, n );
    rfx_tf *Y_p = AA_MEM_REGION_LOCAL_NEW_N( rfx_tf, n );
    for( size_t j = 0; j < k; j ++ ) {
        // apply offsets
        for( size_t i = 0; i < n; i ++ )
            aa_tf_qutr_mul(cor[i].X, tf_off[ i_mark[i] ].data, X_p[i].data );
        // compute registration
        tf_cor( opts, n,
                X_p[0].r.data, 7,
                X_p[0].v.data, 7,
                cor[0].Y, TF_COR_LD,
                cor[0].Y+4, TF_COR_LD,
                Z );

        printf("em %lu:  ", j);
        aa_dump_vec(stdout, Z, 7 );

        // update offsets
        // T_fk*T_off = Z*T_cam  =>  T_off * (Z*T_cam)^* = T_fk^*
        for( size_t i = 0; i < n; i ++ ) {
            double tmp[7];
            aa_tf_qutr_mul(Z, cor[i].Y, tmp);
            aa_tf_qutr_conj(tmp, Y_p[i].data);
        }

        for( size_t i = 0; i < m; i ++ ) {
            size_t beg = i_beg[i];
            tf_cor( opts, i_end[i] - beg,
                    X_conj[beg].r.data, 7,
                    X_conj[beg].v.data, 7,
                    Y_p[beg].r.data, 7,
                    Y_p[beg].v.data, 7,
                    tf_off[i].data );

            printf("  ");
            aa_dump_vec(stdout, tf_off[i].data, 7 );
        }

    }
    aa_mem_region_local_pop(top);
}


int main( int argc, char **argv )
{
    /* Parse */
    for( int c; -1 != (c = getopt(argc, argv, "c:k:i:o:t:d:x:v?")); ) {
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
        case 'i':
            opt_file_id = optarg;
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
                  "  -i ID-FILE,                          Frame id file\n"
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
    double *E_cam=NULL, *E_fk=NULL, *ids=NULL;
    //return 0;
    ssize_t lines_cam = read_mat( opt_file_cam, 7, &E_cam );
    ssize_t lines_fk =  read_mat( opt_file_fk, 7, &E_fk );
    ssize_t lines_id =  read_mat( opt_file_id, 1, &ids );
    if( lines_cam != lines_fk || lines_id != lines_cam ) {
        fprintf(stderr, "Differing line count between `%s', `%s', and `%s'\n",
                opt_file_cam, opt_file_fk, opt_file_id );
        exit(EXIT_FAILURE);
    }

    size_t count = (size_t)lines_cam;

    struct tf_cor *cor = AA_NEW_AR( struct tf_cor, count );
    AA_MEM_ZERO(cor, count);
    aa_cla_dlacpy( ' ', 7, (int)count, E_fk, 7, cor[0].X, TF_COR_LD );
    aa_cla_dlacpy( ' ', 7, (int)count, E_cam, 7, cor[0].Y, TF_COR_LD );
    aa_cla_dlacpy( ' ', 1, (int)count, ids, 1, &cor[0].id, TF_COR_LD );
    aa_aheap_sort( cor, count, sizeof(*cor), &tf_cor_compar );

    // tf correspondences
    double E[7];
    tf_cor( TF_COR_O_ROT_UMEYAMA | TF_COR_O_TRANS_MEAN,
            count,
            cor[0].X, TF_COR_LD,
            cor[0].X+4, TF_COR_LD,
            cor[0].Y, TF_COR_LD,
            cor[0].Y+4, TF_COR_LD,
            E );
    write_tf( global_output, "Umeyama Mean", E );

    tf_cor( TF_COR_O_ROT_DAVENPORT | TF_COR_O_TRANS_MEAN,
            count,
            cor[0].X, TF_COR_LD,
            cor[0].X+4, TF_COR_LD,
            cor[0].Y, TF_COR_LD,
            cor[0].Y+4, TF_COR_LD, E );
    write_tf( global_output, "Davenport Mean", E );

    /* tf_cor( TF_COR_O_ROT_MEDIAN | TF_COR_O_TRANS_MEDIAN, */
    /*         count, */
    /*         cor[0].X, TF_COR_LD, */
    /*         cor[0].X+4, TF_COR_LD, */
    /*         cor[0].Y, TF_COR_LD, */
    /*         cor[0].Y+4, TF_COR_LD, E ); */
    /* write_tf( global_output, "Median", E ); */

    tf_cor( TF_COR_O_ROT_UMEYAMA | TF_COR_O_TRANS_MEDIAN,
            count,
            cor[0].X, TF_COR_LD,
            cor[0].X+4, TF_COR_LD,
            cor[0].Y, TF_COR_LD,
            cor[0].Y+4, TF_COR_LD, E );
    write_tf( global_output, "Umeyama Median", E );

    tf_cor_em( TF_COR_O_ROT_UMEYAMA | TF_COR_O_ROT_DAVENPORT | TF_COR_O_TRANS_MEAN,
               50, count,
               cor, E );
    write_tf( global_output, "EM", E );
}


static ssize_t
read_mat( const char *name, size_t n, double **A )
{

    FILE *f = fopen( name, "r" );
    if( NULL == f ) {
        fprintf(stderr, "Could not open `%s'\n", name );
        exit(EXIT_FAILURE);
    }
    size_t elts = 0;
    ssize_t lines = aa_io_fread_matrix_heap( f, n, A, &elts );
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
    double E_off[7];

    double *E_cam = (double*)malloc( sizeof(double) * 7 * opt_test );
    double *E_fk = (double*)malloc( sizeof(double) * 7 * opt_test );

    aa_tf_qutr_rand( E_true );
    aa_tf_qminimize(E_true);
    rfx_tf_rand( 5.0*M_PI/180, .25e-2, E_off );

    for( size_t i = 0; i < opt_test; i ++ ) {
        size_t j = 7*i;
        {
            // random FK
            aa_tf_qutr_rand( E_fk+j );
            // marker fk = fk * offset
            double mk_fk[7];
            aa_tf_qutr_mul(E_fk+j, E_off, mk_fk );
            aa_tf_qutr_cmul( E_true, mk_fk, E_cam+j );
        }

        if( 0 < opt_d_theta || 0 < opt_d_x ) {
            double tmp[7];
            rfx_tf_corrupt( opt_d_theta, opt_d_x, E_cam+j, tmp );
            memcpy( E_cam+j, tmp, 7*sizeof(tmp[0]) );
        }
   }

    write_tfs( "Camera Pose Estimates", opt_file_cam, opt_test, E_cam );
    write_tfs( "Forward Kinematics Poses", opt_file_fk, opt_test, E_fk );
    double E_out[2][7];
    AA_MEM_CPY(E_out[0], E_off, 7);
    AA_MEM_CPY(E_out[1], E_true, 7);
    write_tfs( "True Registration", opt_file_out, 2, E_out[0] );

    FILE *f_id = fopen( opt_file_id, "w" );
    if( NULL == f_id ) {
        fprintf(stderr, "Could not open `%s'\n", opt_file_id );
        exit(EXIT_FAILURE);
    }
    for( size_t i = 0; i < opt_test; i ++ )
        fprintf(f_id,"0\n");
    fclose(f_id);
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
    aa_tf_qminimize(E);
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
