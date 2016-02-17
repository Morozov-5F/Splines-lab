//
//  main.c
//  Splines
//
//  Created by Evgeny Morozov on 09.11.15.
//  Copyright (c) 2015 Voronezh State University. All rights reserved.
//

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "spline.h"

#define N 100

#define A 0.0
#define B 1.0
#define SPLINE_N 10

double f(double x)
{
    return x * sin(50 * x);
}

double df(double x)
{
    return sin(50 * x) + 50 * cos(50 * x) * x;
}

double d2f(double x)
{
    return 6 * x;
}

int main(int argc, const char * argv[])
{
    spline_t spline;
    spline_init(&spline);
    spline_create_grid(&spline, A, B, SPLINE_N, GRID_UNIFORM);
    spline_init_with_function(&spline, f);
    spline_solve(&spline, df(A), df(B), BORDERS_TYPE_ONE);
    
    unsigned i = 0;
    double step = (B - A) / (N);
    freopen("output.txt", "w", stdout);
    for (i = 0; i <= N; ++i)
        printf("%f %f %f\n", i * step, spline_s(&spline, i * step), f(i * step));
    
#if defined(DEBUG) || defined(_DEBUG)
    dump_spline_struct(&spline, stderr);
#endif
    return 0;
}
