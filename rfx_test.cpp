#include <amino.h>
#include "reflex.hpp"


int main( int argc, char **argv ) {
    reflex::TrapvelWSTrajectory T;
    T.add( 0, (double[3]){0,0,0}, (double[4]){0,0,0,0} );
    T.add( 4.4, (double[3]){1,2,3}, (double[4]){0,0,0,0} );
    int i = T.generate();
    printf("status: %d\n", i );
    T.plot(.01);
    return 0;
}
