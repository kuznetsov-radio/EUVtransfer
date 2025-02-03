#pragma once

int EUVtransfer(int *Lparms, double *Rparms, double *Parms, 
	            double *logTe_rsp, double *rsp, 
				double *logTe_DEM, double *DEM_cor, double *DEM_tr, 
	            double *flux);
int EUVtransferGX(int *Lparms, double *Rparms, double *Parms, 
	              float *logTe_rsp, double *rsp, 
				  float *Qrun, float *Lrun, float *logTe_DEM, float *DEM_cor_run, float *DEM_tr_run, 
	              double *flux);