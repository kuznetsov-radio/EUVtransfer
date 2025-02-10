#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <float.h>
#include "ExtMath.h"
#include "IDLinterface.h"

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
   double T0=Parms[D2(ParmSize, 1, k)];
   double n0=Parms[D2(ParmSize, 2, k)];
   pix_on[k]=(T0>0 && n0>0);
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
  }
 }

 if ((DEM_tr!=0) && (TRfactor>0)) for (int l=0; l<Nch; l++)
 {
  for (int m=0; m<NT_rsp_local; m++) EUV_integrand[m]=rsp_local[D2(NT_rsp_local, m, l)]*DEM_tr[DEM_idx+m]*Te_rsp_local[m];
  double fl=IntTabulated(logTe_rsp_local, EUV_integrand, NT_rsp_local)*log(10.0)*TRfactor;
  flux[D2(fluxSize, 1, l)]=fl;
 }

 for (int m=0; m<2; m++) for (int l=0; l<Nch; l++) flux[D2(fluxSize, m, l)]*=(dS_map/dS_rsp);

 free(EUV_integrand);
 free(rsp_local);
 for (int j=0; j<Nch; j++) delete rsp_spl_arr[j];
 free(rsp_spl_arr);
 free(logTe_rsp_local);
 free(Te_rsp_local);
 free(pix_on);

 return 0;
}

int InterpolateEBTEL(int NQ, int NL, int NT, double Q, double L, double *QgridL, double *LgridL, float *DEM_run, double *DEM)
{
 for (int l=0; l<NT; l++) DEM[l]=0;

 int res=0;

 double Qlog=(finite(Q) && Q>0) ? log(Q) : -1000;
 double Llog=(finite(L) && L>0) ? log(L) : -1000;

 int Lind=value_locate(LgridL, NL, Llog);

 if (Lind>=0 && Lind<(NL-1))
 {
  int Qind1=value_locate(QgridL+Lind*NQ, NQ, Qlog);
  int Qind2=value_locate(QgridL+(Lind+1)*NQ, NQ, Qlog);

  if (Qind1>=0 && Qind1<(NQ-1) && Qind2>=0 && Qind2<(NQ-1))
  {
   res=1;

   double dL=(Llog-LgridL[Lind])/(LgridL[Lind+1]-LgridL[Lind]);
   double dQ1=(Qlog-QgridL[D2(NQ, Qind1, Lind)])/(QgridL[D2(NQ, Qind1+1, Lind)]-QgridL[D2(NQ, Qind1, Lind)]);
   double dQ2=(Qlog-QgridL[D2(NQ, Qind2, Lind+1)])/(QgridL[D2(NQ, Qind2+1, Lind+1)]-QgridL[D2(NQ, Qind2, Lind+1)]);

   for (int l=0; l<NT; l++)
    DEM[l]=DEM_run[D3(NT, NQ, l, Qind1, Lind)]*(1.0-dL)*(1.0-dQ1)+
           DEM_run[D3(NT, NQ, l, Qind1+1, Lind)]*(1.0-dL)*dQ1+
           DEM_run[D3(NT, NQ, l, Qind2, Lind+1)]*dL*(1.0-dQ2)+
           DEM_run[D3(NT, NQ, l, Qind2+1, Lind+1)]*dL*dQ2;
  }
 }

 return res;
}

int EUVtransferGX(int *Lparms, double *Rparms, double *Parms, 
	              float *logTe_rsp, double *rsp, 
				  float *Qrun, float *Lrun, float *logTe_DEM, float *DEM_cor_run, float *DEM_tr_run, 
	              double *flux)
{
 int Nz=Lparms[0];
 int Nchannels=Lparms[1];
 int NT_rsp=Lparms[2];
 int NQ=Lparms[3];
 int NL=Lparms[4];
 int NT_DEM=Lparms[5];

 double dS_map=Rparms[0];
 double dS_rsp=Rparms[1];

 double *LgridL=(double*)malloc(sizeof(double)*NL);
 for (int j=0; j<NL; j++) LgridL[j]=log((double)Lrun[D2(NQ, 0, j)]);

 double *QgridL=(double*)malloc(sizeof(double)*NQ*NL);
 for (int i=0; i<NQ*NL; i++) QgridL[i]=log((double)Qrun[i]);

 double *logTe_rsp_short=(double*)malloc(sizeof(double)*NT_rsp);
 for (int i=0; i<NT_rsp; i++) logTe_rsp_short[i]=logTe_rsp[i];

 double *logTe_DEM_short=(double*)malloc(sizeof(double)*NT_DEM);
 for (int i=0; i<NT_DEM; i++) logTe_DEM_short[i]=logTe_DEM[i];

 double *Parms_short=(double*)malloc(sizeof(double)*Nz*ParmSize);
 double *DEM_cor_short=(double*)malloc(sizeof(double)*Nz*NT_DEM);
 double *DEM_tr_short=(double*)malloc(sizeof(double)*NT_DEM);

 double TRfactor=0;
 int k_first=0;
 int k_last=-1;
 int EUVTR_on=0;

 int done=0;
 for (int i=Nz-1; i>=0; i--) if (!done)
 {
  k_first=i;

  Parms_short[D2(ParmSize, 0, i)]=Parms[D2(ParmSizeGX, 0, i)]; //dz
  Parms_short[D2(ParmSize, 1, i)]=Parms[D2(ParmSizeGX, 1, i)]; //T0
  Parms_short[D2(ParmSize, 2, i)]=Parms[D2(ParmSizeGX, 2, i)]; //n0

  int VoxID=(int)Parms[D2(ParmSizeGX, 3, i)];
  if ((VoxID & 4)!=0) //corona
  {
   int DEM_valid=InterpolateEBTEL(NQ, NL, NT_DEM,
	                              Parms[D2(ParmSizeGX, 4, i)], Parms[D2(ParmSizeGX, 5, i)], 
	                              QgridL, LgridL, DEM_cor_run, 
	                              DEM_cor_short+i*NT_DEM);
   Parms_short[D2(ParmSize, 3, i)]=DEM_valid;
  }
  else Parms_short[D2(ParmSize, 3, i)]=0;

  if (k_last<0) if (Parms_short[D2(ParmSize, 3, i)]==0)
  {
   if (Parms_short[D2(ParmSize, 1, i)]>0 && Parms_short[D2(ParmSize, 2, i)]>0) k_last=i;
  }
  else
  {
   for (int l=0; l<NT_DEM; l++) if (DEM_cor_short[D2(NT_DEM, l, i)]>0)
   {
    k_last=i;
    break;
   }
  }

  if ((VoxID & 2)!=0) //TR
  {
   int DEM_valid=InterpolateEBTEL(NQ, NL, NT_DEM,
	                              Parms[D2(ParmSizeGX, 4, i)], Parms[D2(ParmSizeGX, 5, i)], 
	                              QgridL, LgridL, DEM_tr_run, 
	                              DEM_tr_short);
   TRfactor=DEM_valid ? Parms[D2(ParmSizeGX, 6, i)] : 0;
   done=1;

   if ((VoxID & 8)!=0) EUVTR_on=1;
  }
 }

 double Rparms_short[3]={dS_map, dS_rsp, TRfactor};

 if (k_last>=0)
 {
  int Lparms_short[4]={k_last-k_first+1, Nchannels, NT_rsp, NT_DEM}; 

  EUVtransfer(Lparms_short, Rparms_short, Parms_short+k_first*ParmSize, 
	          logTe_rsp_short, rsp, 
	          logTe_DEM_short, DEM_cor_short+k_first*NT_DEM, DEM_tr_short, flux);
 }
 else memset(flux, 0, sizeof(double)*Nchannels*fluxSize);

 if (EUVTR_on) for (int l=0; l<Nchannels; l++) flux[D2(fluxSize, 2, l)]=1;

 free(Parms_short);
 free(DEM_cor_short);
 free(DEM_tr_short);
 free(logTe_DEM_short);
 free(logTe_rsp_short);
 free(QgridL);
 free(LgridL);

 return 0;
}