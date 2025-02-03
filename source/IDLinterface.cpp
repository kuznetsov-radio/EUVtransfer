#include "IDLinterface.h"
#include "EUVmain.h"

#ifndef LINUX
#include <ppl.h>
#else
#include <omp.h>
#endif

#ifndef LINUX
extern "C" __declspec(dllexport) int GET_EUV(int argc, void **argv)
#else
extern "C" int GET_EUV(int argc, void **argv)
#endif
{
 int res=0;

 if (argc==9)
 {
  int *Lparms=(int*)argv[0];
  double *Rparms=(double*)argv[1];
  double *Parms=(double*)argv[2];
  double *logTe_rsp=(double*)argv[3];
  double *rsp=(double*)argv[4];
  double *logTe_DEM=(double*)argv[5];
  double *DEM_cor=(double*)argv[6];
  double *DEM_tr=(double*)argv[7];
  double *flux=(double*)argv[8];

  res=EUVtransfer(Lparms, Rparms, Parms, logTe_rsp, rsp, logTe_DEM, DEM_cor, DEM_tr, flux);
 }
 else res=-1;

 return res;
}

#ifndef LINUX
extern "C" __declspec(dllexport) int GET_GX_EUV(int argc, void **argv)
#else
extern "C" int GET_GX_EUV(int argc, void **argv)
#endif
{
 int res=0;

 if (argc==11)
 {
  int *Lparms=(int*)argv[0];
  double *Rparms=(double*)argv[1];
  double *Parms=(double*)argv[2];
  float *logTe_rsp=(float*)argv[3];
  double *rsp=(double*)argv[4];
  float *Qrun=(float*)argv[5];
  float *Lrun=(float*)argv[6];
  float *logTe_DEM=(float*)argv[7];
  float *DEM_cor_run=(float*)argv[8];
  float *DEM_tr_run=(float*)argv[9];
  double *flux=(double*)argv[10];

  res=EUVtransferGX(Lparms, Rparms, Parms, logTe_rsp, rsp, Qrun, Lrun, logTe_DEM, DEM_cor_run, DEM_tr_run, flux);
 }
 else res=-1;

 return res;
}

#ifndef LINUX
extern "C" __declspec(dllexport) int GET_EUV_SLICE(int argc, void** argv)
#else
extern "C" int GET_EUV_SLICE(int argc, void** argv)
#endif
{
 int res=0;

 if (argc==9)
 {
  int *Lparms_M=(int*)argv[0];
  double *Rparms_M=(double*)argv[1];
  double *Parms_M=(double*)argv[2];
  double *logTe_rsp=(double*)argv[3];
  double *rsp=(double*)argv[4];
  double *logTe_DEM=(double*)argv[5];
  double *DEM_cor_M=(double*)argv[6];
  double *DEM_tr_M=(double*)argv[7];
  double *flux_M=(double*)argv[8];

  int Npix=Lparms_M[0];
  int Nz=Lparms_M[1];
  int Nch=Lparms_M[2];
  int NT_DEM=Lparms_M[4];

  #ifndef LINUX
  concurrency::parallel_for(0, Npix, [&](int pix)
  #else
  #pragma omp parallel for
  for(int pix=0; pix<Npix; pix++)
  #endif
  {
   void *ARGV[9];

   ARGV[0]=(void*)(Lparms_M+1);
   ARGV[1]=(void*)(Rparms_M+pix*RpSize);
   ARGV[2]=(void*)(Parms_M+pix*Nz*ParmSize);
   ARGV[3]=(void*)logTe_rsp;
   ARGV[4]=(void*)rsp;
   ARGV[5]=(void*)logTe_DEM;
   ARGV[6]=(void*)(DEM_cor_M+pix*Nz*NT_DEM);
   ARGV[7]=(void*)(DEM_tr_M+pix*NT_DEM);
   ARGV[8]=(void*)(flux_M+pix*Nch*fluxSize);

   GET_EUV(9, ARGV);
  #ifndef LINUX
  });
  #else
  }
  #endif
 }
 else res=-1;

 return res;
}

#ifndef LINUX
extern "C" __declspec(dllexport) int GET_GX_EUV_SLICE(int argc, void** argv)
#else
extern "C" int GET_GX_EUV_SLICE(int argc, void** argv)
#endif
{
 int res=0;

 if (argc==11)
 {
  int *Lparms_M=(int*)argv[0];
  double *Rparms_M=(double*)argv[1];
  double *Parms_M=(double*)argv[2];
  float *logTe_rsp=(float*)argv[3];
  double *rsp=(double*)argv[4];
  float *Qrun=(float*)argv[5];
  float *Lrun=(float*)argv[6];
  float *logTe_DEM=(float*)argv[7];
  float *DEM_cor_run=(float*)argv[8];
  float *DEM_tr_run=(float*)argv[9];
  double *flux_M=(double*)argv[10];

  int Npix=Lparms_M[0];
  int Nz=Lparms_M[1];
  int Nch=Lparms_M[2];

  #ifndef LINUX
  concurrency::parallel_for(0, Npix, [&](int pix)
  #else
  #pragma omp parallel for
  for(int pix=0; pix<Npix; pix++)
  #endif
  {
   void *ARGV[11];

   ARGV[0]=(void*)(Lparms_M+1);
   ARGV[1]=(void*)(Rparms_M+pix*RpSizeGX);
   ARGV[2]=(void*)(Parms_M+pix*Nz*ParmSizeGX);
   ARGV[3]=(void*)logTe_rsp;
   ARGV[4]=(void*)rsp;
   ARGV[5]=(void*)Qrun;
   ARGV[6]=(void*)Lrun;
   ARGV[7]=(void*)logTe_DEM;
   ARGV[8]=(void*)DEM_cor_run;
   ARGV[9]=(void*)DEM_tr_run;
   ARGV[10]=(void*)(flux_M+pix*Nch*fluxSize);

   GET_GX_EUV(11, ARGV);
  #ifndef LINUX
  });
  #else
  }
  #endif
 }
 else res=-1;

 return res;
}