#include <amino.h>
//#include <schkin.h>
#include "reflex.hpp"
/*
void trapvel() {
    reflex::TrapvelWS T;
    T.add( 0, (double[3]){0,0,0}, (double[4]){0,0,0,1} );
    T.add( 4.0, (double[3]){1,2,3}, (double[4]){.1,.2,.3,1} );
    int i = T.generate();
    printf("status: %d\n", i );
    i = T.validate();
    printf("validity: %d\n", i );
    if( 0 == i )
        T.plot(.01);
}

void jacobian() {
    // initial configuration
    double q0[7];
    AA_SET_AR(q0,M_PI/4);

    // initial ws pos
    double R0[9], v0[3], rq0[4];
    schkin_lwa3_fk(q0, NULL, NULL, NULL, NULL, R0, v0);
    aa_tf_rotmat2quat(R0, rq0);

    // final ws pos
    double v1[3], rq1[4];
    aa_fcpy(rq1, rq0, 4);
    aa_fcpy(v1, v0, 3);
    v1[2] += .3;
    v1[0] -= .1;
    v1[1] -= .2;

    aa_dump_vec(stderr, v0, 3);
    aa_dump_vec(stderr, v1, 3);

    // trajectory
    reflex::TrapvelWS T;
    T.add( 0, v0, rq0 );
    T.add( 10, v1, rq1 );
    int i = T.generate();
    printf("trapvel generation: %d\n", i );
    i = T.validate();
    printf("trapvel validity: %d\n", i );
    //if( 0 == i )
        //T.jacobian_plot(.005, 7, q0, schkin_lwa3_fk, schkin_lwa3_jacobian, .001);

}


void ws() {
    // initial configuration
    double q0[7];
    AA_SET_AR(q0,M_PI/4);

    // initial ws pos
    double R0[9], v0[3], rq0[4];
    schkin_lwa3_fk(q0, NULL, NULL, NULL, NULL, R0, v0);
    aa_tf_rotmat2quat(R0, rq0);

    // final ws pos
    double v1[3], rq1[4] = {0,0,0,1};
    aa_fcpy(rq1, rq0, 4);
    aa_fcpy(v1, v0, 3);
    v1[2] += .3;
    v1[0] -= .1;
    v1[1] -= .2;

    aa_dump_vec(stderr, v0, 3);
    aa_dump_vec(stderr, v1, 3);

    // trajectory
    reflex::TrapvelWS T;
    T.add( 0, v0, rq0 );
    T.add( 10, v1, rq1 );
    int i = T.generate();
    printf("trapvel generation: %d\n", i );
    i = T.validate();
    printf("trapvel validity: %d\n", i );
    double kp[6];
    aa_fset(kp,3,3);
    aa_fset(kp+3,1,1);
    //AA_SET_AR(kp, I);
    //if( 0 == i )
        //T.wsctrl_plot(.005, 7, q0, schkin_lwa3_fk, schkin_lwa3_jacobian, .001,
                      //kp);

}
*/
/*
void parablend() {
    reflex::ParaBlendWS T;
    T.add( 0, AA_FAR(0,0,0), AA_TF_QUAT_IDENT );
    T.add( 1, AA_FAR(.2,.2,.2), AA_TF_QUAT_IDENT );
    T.add( 2, AA_FAR(.4,.4,.4), AA_TF_QUAT_IDENT );
    T.generate();
}
*/

int main( int argc, char **argv ) {
    (void) argc;
    (void) argv;
    //trapvel();
    //jacobian();
    //ws();
    //parablend();
    return 0;

}
