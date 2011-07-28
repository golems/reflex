
#include <amino.h>
#include "reflex.h"



int main( int argc, char **argv ) {
    (void)argc;
    (void)argv;


    {
    double A[] = {1,2,3,4};
    double B[] = {1,2};
    double C[] = {2,1,2,1};
    double u[] = {5};
    double x[] = {4,2};
    double z[] = {8,2};
    double K[] = {8,2,1,8};
    double dx[2], zwork[2];

    rfx_lqg_observe(2,1,2, A,B,C, x,u,z, K, dx, zwork);

    assert( aa_veq( 2, dx, AA_FAR(-21,-14), .00001 ) );
    }
}





