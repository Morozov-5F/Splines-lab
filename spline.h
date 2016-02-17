//
//  spline.h
//  Splines
//
//  Created by Evgeny Morozov on 09.11.15.
//  Copyright (c) 2015 Voronezh State University. All rights reserved.
//

#ifndef __SPLINE_H__
#define __SPLINE_H__

#include <stdlib.h>
#include <stdio.h>

typedef struct splineDef
{
    double * m;
    double * f;
    double * x;
    size_t   n;
} spline_t;

typedef enum
{
    GRID_UNIFORM,
    GRID_LOGARITHMIC,
} grid_params_t;

typedef enum
{
    BORDERS_TYPE_ONE,
    BORDERS_TYPE_TWO,
    BORDERS_CYCLIC,
} border_params_t;

void spline_init(spline_t * spline);
void spline_create_grid(spline_t * spline, double a, double b, size_t num, grid_params_t parameters);
void spline_init_with_function(spline_t * spline, double (*f)(double x));
void spline_solve(spline_t * spline, double bc_l, double bc_r, border_params_t params);
void spline_destroy(spline_t * spline);

double spline_s(const spline_t * spline, double x);

#if defined(DEBUG) || defined (_DEBUG)

void dump_spline_struct(const spline_t * spline, FILE * stream);

#endif

#endif /* defined(__SPLINE_H__) */
