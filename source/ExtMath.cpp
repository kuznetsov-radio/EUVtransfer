#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "ExtMath.h"

double IntTabulated(double *x, double *y, int N)
{
 double s=0;
 for (int i=1; i<N; i++) s+=0.5*(y[i-1]+y[i])*(x[i]-x[i-1]);
 return s;
}

void spline_init(double *x, double *y, int n, double yp1, double ypn, double *y2)
{
 int i, k;
 double p, qn, sig, un, *u;
 u=(double*)malloc(sizeof(double)*n);
 if (!finite(yp1)) y2[0]=u[0]=0.0; 
 else 
 { 
  y2[0]=-0.5;
  u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
 }
 for (i=1; i<n-1; i++) 
 { 
  sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
  p=sig*y2[i-1]+2.0;
  y2[i]=(sig-1.0)/p;
  u[i]=(y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])/(x[i]-x[i-1]);
  u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
 }
 if (!finite(ypn)) qn=un=0.0;
 else 
 { 
  qn=0.5;
  un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
 }
 y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
 for (k=n-2; k>=0; k--) y2[k]=y2[k]*y2[k+1]+u[k]; 
 free(u);
}

void spline_short(double *x, double *y, int n, double *y2)
{
 double dxl, dxr;
 double K1[3];

 dxl=x[1]-x[0];
 dxr=x[2]-x[1];
 K1[0]=-(2.0*dxl+dxr)/dxl/(dxl+dxr);
 K1[1]=(dxr+dxl)/(dxr*dxl);
 K1[2]=-dxl/dxr/(dxr+dxl);
 double y1l=K1[0]*y[0]+K1[1]*y[1]+K1[2]*y[2];

 dxl=x[n-2]-x[n-3];
 dxr=x[n-1]-x[n-2];
 K1[0]=dxr/dxl/(dxr+dxl);
 K1[1]=-(dxr+dxl)/(dxr*dxl);
 K1[2]=(2.0*dxr+dxl)/dxr/(dxr+dxl);
 double y1r=K1[0]*y[n-3]+K1[1]*y[n-2]+K1[2]*y[n-1];

 spline_init(x, y, n, y1l, y1r, y2);
}

Spline :: Spline(int _N, double *x, double *y)
{
 N=_N;

 x_arr=(double*)malloc(sizeof(double)*N);
 y_arr=(double*)malloc(sizeof(double)*N);
 y2_arr=(double*)malloc(sizeof(double)*N);

 for (int i=0; i<N; i++)
 {
  x_arr[i]=x[i];
  y_arr[i]=y[i];
 }

 spline_short(x_arr, y_arr, N, y2_arr);
}

Spline :: ~Spline()
{
 free(x_arr);
 free(y_arr);
 free(y2_arr);
}

void spline_interp(double *xa, double *ya, double *y2a, int n, double x, double *y, double *y1)
{
 int klo, khi, k;
 double h, b, a;

 if ((x<xa[0]) || (x>xa[n-1]))
 {
  if (y) *y=0;
  if (y1) *y1=0;
 }
 else
 {
  klo=0; 
  khi=n-1;
  while (khi-klo>1) 
  {
   k=(khi+klo)>>1;
   if (xa[k]>x) khi=k;
   else klo=k;
  } 
 
  h=xa[khi]-xa[klo];
  a=(xa[khi]-x)/h; 
  b=(x-xa[klo])/h; 

  if (y) *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
  if (y1) *y1=(ya[khi]-ya[klo])/h+((1.0-3.0*a*a)*y2a[klo]+(3.0*b*b-1)*y2a[khi])*h/6.0;
 }
}

void Spline :: Interpolate(double x, double *y, double *y1)
{
 spline_interp(x_arr, y_arr, y2_arr, N, x, y, y1);
}

int value_locate(double *a, int n, double x)
{
 int asc=a[n-1]>a[0];

 if (asc ? x<a[0] : x>a[0]) return -1;
 if (asc ? x>=a[n-1] : x<=a[n-1]) return n-1;

 int j, j1, l;
 j=0; 
 j1=n-1; 
 while (j1-j>1) 
 { 
  l=(j+j1)>>1; 
  if (asc ? a[l]>x : a[l]<x) j1=l; 
  else j=l; 
 } 
 return j;
} 