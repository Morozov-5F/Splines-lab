//
//  spline.c
//  Splines
//
//  Created by Evgeny Morozov on 09.11.15.
//  Copyright (c) 2015 Voronezh State University. All rights reserved.
//

#include "spline.h"

#include <stdio.h>
#include <math.h>

#if defined(DEBUG) || defined (_DEBUG)

void dump_spline_struct(const spline_t * spline, FILE * stream)
{
    fprintf(stream, "\nSpline dump: ");
    fprintf(stream, "size = %zu", spline->n);
    size_t i;
    fprintf(stream, "\nGRID: \n");
    for (i = 0; i <= spline->n; ++i)
    {
        fprintf(stream, "%f ", spline->x[i]);
    }
    fprintf(stream, "\nf_i: \n");
    for (i = 0; i <= spline->n; ++i)
    {
        fprintf(stream, "%f ", spline->f[i]);
    }
    fprintf(stream, "\nm_i: \n");
    for (i = 0; i <= spline->n; ++i)
    {
        fprintf(stream, "%f ", spline->m[i]);
    }
    fprintf(stream, "\n");
}

#endif

static inline double phi_1 (double t)
{
    return (1 - t) * (1 - t) * ( 1 + 2 * t);
}

static inline double phi_2 (double t)
{
    return t * t * (3 - 2 * t);
}

static inline double phi_3 (double t)
{
    return t * (1 - t) * (1 - t);
}

static inline double phi_4 (double t)
{
    return - t * t * (1 - t);
}


static void binary_search(double * a, double x, size_t n, size_t * index)
{
    if (a[0] > x)
    {
        *index = 0;
        return;
    }
    else if (a[n] < x)
    {
        *index = n;
        return;
    }
    size_t left = 0, right = n + 1;
    
    while (right - left > 1)
    {
        size_t m = (left + right) / 2;
        if (a[m] > x)
            right = m;
        else
            left = m;
    }
    
    *index = left;
    return;
}

static void solve_tridiagonal(double * a, double * b, double * c, double * d, double * res, size_t n)
{
    double v = -b[0] / a[0], u = d[0] / a[0], buf;
    a[0] = v;
    b[0] = u;
    size_t i;
    /* Straight */
    for (i = 1; i <= n - 1; ++i)
    {
        buf = (c[i] * a[i - 1] + a[i]);
        v = - b[i] / buf;
        u = (d[i] - c[i] * b[i - 1]) / buf;
        
        a[i] = v;
        b[i] = u;
    }
    u = (d[n] - c[n] * b[n - 1]) / (c[n] * a[n - 1] + a[n]);
    res[n] = u;
    fprintf(stderr, "%f %f %f\n", res[n], a[n - 1], b[n - 1]);
    /* Reverse */
    for (i = n; i-- != 0;)
    {
        fprintf(stderr, "%f %f %f\n", res[i + 1], a[i], b[i]);
        res[i] = a[i] * res[i + 1] + b[i];
    }
}

void spline_init(spline_t * spline)
{
    spline->f = NULL;
    spline->m = NULL;
    spline->x = NULL;
    spline->n = 0;
}

void spline_destroy(spline_t * spline)
{
    if (spline->f)
        free(spline->f);
    if (spline->m)
        free(spline->m);
    if (spline->x)
        free(spline->x);
    spline->n = 0;
}

void spline_create_grid(spline_t * spline, double a, double b, size_t num, grid_params_t parameters)
{
    if (spline == NULL)
    {
        return;
    }

    if (spline->x != NULL)
    {
        return;
    }
    
    spline->x = (double *)calloc(num + 1, sizeof(double));
    spline->n = num;
    unsigned i;
    
    double step;
    for (i = 0, step = (b - a) / (num); i <= num; ++i)
    {
        spline->x[i] = a + step * i;
    }
}

void spline_init_with_function(spline_t * spline, double (*f)(double x))
{
    if (spline == NULL)
    {
        return;
    }
    
    if (spline->x == NULL)
    {
        return;
    }

    spline->f = (double *)calloc(spline->n + 1, sizeof(double));
    unsigned i;
    for (i = 0; i <= spline->n; ++i)
    {
        spline->f[i] = f(spline->x[i]);
    }
}


void spline_solve(spline_t * spline, double bc_l, double bc_r, border_params_t params)
{
    if (spline == NULL)
    {
        return;
    }
    
    if (spline->x == NULL || spline->f == NULL)
    {
        return;
    }
    
    size_t n = spline->n;
    spline->m = (double *)calloc(n + 1, sizeof(double));
    double *a, *b, *c, *d;
    a = (double *)calloc(n + 1, sizeof(double));
    b = (double *)calloc(n + 1, sizeof(double));
    c = (double *)calloc(n + 1, sizeof(double));
    d = (double *)calloc(n + 1, sizeof(double));

    switch (params)
    {
        case BORDERS_TYPE_ONE:
            a[0] = 1; d[0] = bc_l;
            a[n] = 1; d[n] = bc_r;
            break;
            
        case BORDERS_TYPE_TWO:
            a[0] = 2; b[0] = 1;
            d[0] = (3 * spline->f[1] - spline->f[0]) / (spline->x[1] - spline->x[0]) - bc_l * (spline->x[1] - spline->x[0]) / 2;
           
            c[n] = 1; a[spline->n] = 2;
            d[n] = (3 * spline->f[n] - spline->f[n - 1]) / (spline->x[n] - spline->x[n - 1]) + bc_r * (spline->x[n] - spline->x[n - 1]) / 2;
            break;
            
        default:
            return;
            break;
    }
    
    unsigned i;
    for (i = 1; i <= n - 1; ++i)
    {
        double h_prev = spline->x[i] - spline->x[i - 1], h_curr = spline->x[i + 1] - spline->x[i];
        double mu_i = h_prev / (h_prev + h_curr);
        
        c[i] = 1 - mu_i;
        a[i] = 2;
        b[i] = mu_i;
        d[i] = 3 * mu_i * (spline->f[i + 1] - spline->f[i]) / h_curr + 3 * (1 - mu_i) * (spline->f[i] - spline->f[i - 1]) / h_prev;
    }
    
    solve_tridiagonal(a, b, c, d, spline->m, n);
    
    free(a);
    free(b);
    free(c);
    free(d);
}

double spline_s(const spline_t * spline, double x)
{
    if (spline == NULL)
        return NAN;
    
    if (spline->f == NULL || spline->x == NULL || spline->m == NULL)
        return NAN;
    
    double res, t;
    size_t i = 0;
    
    binary_search(spline->x, x, spline->n, &i);
    
    double h_i = spline->x[i + 1] - spline->x[i];
    
    t = (x - spline->x[i]) / h_i;

#if defined(DEBUG) || (_DEBUG)
    fprintf(stderr, "x = %f, i = %zu, t = %f\n", x, i, t);
#endif
    
    res = phi_1(t) * spline->f[i] + phi_2(t) * spline->f[i + 1] + phi_3(t) * h_i * spline->m[i] + phi_4(t) * h_i * spline->m[i + 1];
    return res;
}