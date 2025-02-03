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

//-------------------------------------

const char* arr1gx[]={
 "      N_vox;       1   ;int;     data;                    Number of voxels",
 " N_channels;       7   ;int;     data;              Number of EUV channels",
 " N_temp_rsp;     101   ;int;     data;  No. of temperatures (resp. matrix)",
 "  N_Q_EBTEL;      60   ;int;     data;   No. of Q values in the EBTEL grid",
 "  N_L_EBTEL;      49   ;int;     data;   No. of L values in the EBTEL grid",
 " N_temp_DEM;     451   ;int;     data;     No. of temperatures (DEM array)"
};

#define N1gx 6

const char* arr1sgx[]={
 "      N_pix;       1   ;int;     data;                    Number of pixels",
 "      N_vox;       1   ;int;     data;                    Number of voxels",
 " N_channels;       7   ;int;     data;              Number of EUV channels",
 " N_temp_rsp;     101   ;int;     data;  No. of temperatures (resp. matrix)",
 "  N_Q_EBTEL;      60   ;int;     data;   No. of Q values in the EBTEL grid",
 "  N_L_EBTEL;      49   ;int;     data;   No. of L values in the EBTEL grid",
 " N_temp_DEM;     451   ;int;     data;     No. of temperatures (DEM array)"
};

#define N1sgx 7

const char* arr2gx[]={
 "     dS_map;     1.0   ;area;    data;                 Source/pixel area",
 "     dS_rsp;    0.36   ;area;    data;   Default instrumental pixel area"
};

#define N2gx 2

const char* arr3gx[]={
 "         dR;   1E+09   ;cm;      data;                Source/voxel depth",
 "        T_0;   1E+06   ;K;       data;                Plasma temperature",
 "        n_0;   1E+09   ;cm^{-3}; data;                  Electron density",
 "         ID;       0   ;0/1/2/3; data;            Corona and/or TR flags",
 "          Q;   1E-02   ;relative;data;                    Heating rate Q",
 "          L;   1E+08   ;cm;      data;                     Loop length L",
 "   TRfactor;     1.0   ;factor;  data; Transition region geometry factor",
 "   reserved;       0   ;int;     data;                          reserved",
 "   reserved;       0   ;int;     data;                          reserved",
 "   reserved;       0   ;int;     data;                          reserved"
};

#define N3gx 10

//---------------------------------

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

#ifndef LINUX
extern "C" __declspec(dllexport) float GET_GX_PARMS(int argc, void **argv)
#else
extern "C" float GET_GX_PARMS(int argc, void **argv)
#endif
{
 WriteParms(arr1gx, "Long_input.txt",  N1gx, 0);
 WriteParms(arr2gx, "Real_input.txt",  N2gx, 0);
 WriteParms(arr3gx, "Parms_input.txt", N3gx, 0);
 return 0;
}

#ifndef LINUX
extern "C" __declspec(dllexport) float GET_GX_PARMS_SLICE(int argc, void **argv)
#else
extern "C" float GET_GX_PARMS_SLICE(int argc, void **argv)
#endif
{
 WriteParms(arr1sgx, "Long_input.txt",  N1sgx, 0);
 WriteParms(arr2gx,  "Real_input.txt",  N2gx, 0);
 WriteParms(arr3gx,  "Parms_input.txt", N3gx, 0);
 return 0;
}