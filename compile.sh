#!/bin/bash

cd pkg/src
rm --verbose *.o *.so
gcc -std=gnu99 -I/usr/share/R/include -DNDEBUG -fpic -fopenmp -O3 -Wall -pipe  -g  -c *.c 
#gcc -std=gnu99 -I/usr/share/R/include -DNDEBUG -fpic -O2 -Wall -pipe  -g  -c *.c 
gcc -std=gnu99 -shared -o gower.so *.o -L/usr/lib/R/lib -lR
cd ../../

