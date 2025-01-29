#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include "ExtMath.h"
#include "IDLinterface.h"
#include "Messages.h"

int EUVtransfer(int *Lparms, double *Rparms, double *Parms, 
	            double *logTe_rsp, double *rsp, 
				double *logTe_DEM, double *DEM_cor, double *DEM_tr, 
	            double *flux)
{
 int Nz=Lparms[0];
 int Nch=Lparms[1];
 int NT_rsp=Lparms[2];
 int NT_DEM=Lparms[3];

 double dS_map=Rparms[0];
 double dS_rsp=Rparms[1];
 double TRfactor=Rparms[2];

 int *pix_on=(int*)malloc(Nz*sizeof(int));

 for (int k=0; k<Nz; k++) 
 {
  int DEM_on=(int)Parms[D2(ParmSize, 3, k)];
  
  if (DEM_on==0)
  {
   double n0=Parms[D2(ParmSize, 2, k)];
   pix_on[k]=(n0>0);
  }
  else
  {
   pix_on[k]=0;

   for (int j=0; j<NT_DEM; j++) if (DEM_cor[D2(NT_DEM, j, k)]>0)
   {
    pix_on[k]=1;
    break;
   }
  }
 }

 double Tmin=max(logTe_rsp[0], logTe_DEM[0]);
 double Tmax=min(logTe_rsp[NT_rsp-1], logTe_DEM[NT_DEM-1]);

 int i1=-1;
 int i2=0;
 for (int i=0; i<NT_DEM; i++)
 {
  if ((logTe_DEM[i]>=Tmin) && (i1<0)) i1=i;
  if (logTe_DEM[i]<=Tmax) i2=i;
 }
 int NT_rsp_local=i2-i1+1;
 int DEM_idx=i1;
 double *logTe_rsp_local=(double*)malloc(NT_rsp_local*sizeof(double));
 double *Te_rsp_local=(double*)malloc(NT_rsp_local*sizeof(double));
 for (int i=0; i<NT_rsp_local; i++)
 {
  logTe_rsp_local[i]=logTe_DEM[DEM_idx+i];
  Te_rsp_local[i]=pow(10.0, logTe_rsp_local[i]);
 }

 Spline **rsp_spl_arr=(Spline**)malloc(Nch*sizeof(Spline*));
 double *rsp_local=(double*)malloc(NT_rsp_local*Nch*sizeof(double));
 for (int j=0; j<Nch; j++)
 {
  rsp_spl_arr[j]=new Spline(NT_rsp, logTe_rsp, rsp+j*NT_rsp);
  for (int i=0; i<NT_rsp_local; i++) 
   rsp_spl_arr[j]->Interpolate(logTe_rsp_local[i], rsp_local+i+j*NT_rsp_local, 0);
 }

 double *EUV_integrand=(double*)malloc(NT_rsp_local*sizeof(double));

 memset(flux, 0, sizeof(double)*Nch*fluxSize);

 for (int k=0; k<Nz; k++) if (pix_on[k])
 {
  double dz=Parms[D2(ParmSize, 0, k)];
  double T0=Parms[D2(ParmSize, 1, k)];
  double n0=Parms[D2(ParmSize, 2, k)];
  int DEM_on=(int)Parms[D2(ParmSize, 3, k)];

  double logTe_iso=DEM_on ? 0 : log10(T0);

  for (int l=0; l<Nch; l++)
  {
   double fl=0;

   if (DEM_on)
   {
	for (int m=0; m<NT_rsp_local; m++) 
	 EUV_integrand[m]=rsp_local[D2(NT_rsp_local, m, l)]*DEM_cor[D2(NT_DEM, DEM_idx+m, k)]*Te_rsp_local[m];
	fl=IntTabulated(logTe_rsp_local, EUV_integrand, NT_rsp_local)*log(10.0)*dz;
   }
   else
   {
	double r=0;
	rsp_spl_arr[l]->Interpolate(logTe_iso, &r, 0);
	fl=sqr(n0)*dz*r;
   }

   flux[D2(fluxSize, 0, l)]+=fl;
   flux[D2(fluxSize, 1, l)]+=fl;
  }
 }

 if ((DEM_tr!=0) && (TRfactor>0)) for (int l=0; l<Nch; l++)
 {
  for (int m=0; m<NT_rsp_local; m++) EUV_integrand[m]=rsp_local[D2(NT_rsp_local, m, l)]*DEM_tr[DEM_idx+m]*Te_rsp_local[m];
  double fl=IntTabulated(logTe_rsp_local, EUV_integrand, NT_rsp_local)*log(10.0)*TRfactor;
  flux[D2(fluxSize, 1, l)]+=fl;
 }

 for (int l=0; l<Nch*fluxSize; l++) flux[l]*=(dS_map/dS_rsp);

 free(EUV_integrand);
 free(rsp_local);
 for (int j=0; j<Nch; j++) delete rsp_spl_arr[j];
 free(rsp_spl_arr);
 free(logTe_rsp_local);
 free(Te_rsp_local);
 free(pix_on);

 return 0;
}