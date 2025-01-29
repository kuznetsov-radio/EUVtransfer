#pragma once

#ifdef LINUX
#define finite isfinite
#else
#define finite _finite
#endif

inline double sqr(double x)
{
 return x*x;
}

inline double max(double a, double b)
{
 return (a>b) ? a : b;
}

inline double min(double a, double b)
{
 return (a<b) ? a : b;
}

double IntTabulated(double *x, double *y, int N);

class Spline
{
 int N;
 double *x_arr, *y_arr, *y2_arr;
 public:
 Spline(int _N, double *x, double *y);
 ~Spline();
 void Interpolate(double x, double *y, double *y1);
};