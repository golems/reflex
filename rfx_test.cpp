#include <amino.h>
#include "reflex.hpp"


int main( int argc, char **argv ) {
    (void) argc;
    (void) argv;
    reflex::TrapvelWS T;
    T.add( 0, (double[3]){0,0,0}, (double[4]){0,0,0,1} );
    T.add( 8.4, (double[3]){1,2,3}, (double[4]){1,1,1,1} );
    int i = T.generate();
    printf("status: %d\n", i );
    i = T.validate();
    printf("validity: %d\n", i );
    if( 0 == i )
        T.plot(.01);
    return 0;
}
