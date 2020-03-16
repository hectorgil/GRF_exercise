#!/bin/bash
gcc -std=c99  main_code_power_spectrum.c -O3 -lm -lfftw3 -o file.out
./file.out
