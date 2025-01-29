#include <stdio.h>

const char* arr1[]={
 "      N_vox;       1   ;int;     data;                    Number of voxels",
 " N_channels;       7   ;int;     data;              Number of EUV channels",
 " N_temp_rsp;     101   ;int;     data;  No. of temperatures (resp. matrix)",
 " N_temp_DEM;     451   ;int;     data;     No. of temperatures (DEM array)"
};

#define N1 4

const char* arr1s[]={
 "      N_pix;       1   ;int;     data;                    Number of pixels",
 "      N_vox;       1   ;int;     data;                    Number of voxels",
 " N_channels;       7   ;int;     data;              Number of EUV channels",
 " N_temp_rsp;     101   ;int;     data;  No. of temperatures (resp. matrix)",
 " N_temp_DEM;     451   ;int;     data;     No. of temperatures (DEM array)"
};

#define N1s 5

const char* arr2[]={
 "     dS_map;     1.0   ;area;    data;                 Source/pixel area",
 "     dS_rsp;    0.36   ;area;    data;   Default instrumental pixel area",
 "   TRfactor;     1.0   ;factor;  data; Transition region geometry factor"
};

#define N2 3

const char* arr3[]={
 "         dR;   1E+09   ;cm;      data;                Source/voxel depth",
 "        T_0;   1E+06   ;K;       data;                Plasma temperature",
 "        n_0;   1E+09   ;cm^{-3}; data;                  Electron density",
 "     DEM_on;       0   ;1/0;     data;                        DEM on/off",
 "   reserved;       0   ;int;     data;                          reserved",
 "   reserved;       0   ;int;     data;                          reserved"
};

#define N3 6

void WriteParms(const char **arr, const char *fname, int N, int add)
{
 FILE *F=fopen(fname, add ? "a" : "w");
 if (F)
 {
  for (int i=0; i<N; i++) fprintf(F, "%s\n", arr[i]);
  fclose(F);
 }
}

#ifndef LINUX
extern "C" __declspec(dllexport) float GET_PARMS(int argc, void **argv)
#else
extern "C" float GET_PARMS(int argc, void **argv)
#endif
{
 WriteParms(arr1, "Long_input.txt",  N1, 0);
 WriteParms(arr2, "Real_input.txt",  N2, 0);
 WriteParms(arr3, "Parms_input.txt", N3, 0);
 return 0;
}

#ifndef LINUX
extern "C" __declspec(dllexport) float GET_PARMS_SLICE(int argc, void **argv)
#else
extern "C" float GET_PARMS_SLICE(int argc, void **argv)
#endif
{
 WriteParms(arr1s, "Long_input.txt",  N1s, 0);
 WriteParms(arr2,  "Real_input.txt",  N2, 0);
 WriteParms(arr3,  "Parms_input.txt", N3, 0);
 return 0;
}