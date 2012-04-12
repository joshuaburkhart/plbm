/*
distributed under the terms of the GNU General Public License
Copyright 2012 Joshua Burkhart
*/

#include <stdio.h>

int main(int argc,int argv[]) {
    int nan = 0 / argv[1]; // user is expected to supply a 0
    int i;
    for(i = -1; i < 2; i++) {
        if(nan > i) {
            printf("nan is > %i\n",i);
        }
        if(nan < i) {
            printf("nan is < %i\n",i);
        }
        if(nan == i) {
            printf("nan is == %i\n",i);
        }
        if(nan != i) {
            printf("nan is != %i\n",i);
        }
    }
    return 0;
}
