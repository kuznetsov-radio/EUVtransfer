pro Example_EUV_GX_MultiThreads
 ;load an EUV response matrix: 
 r=gx_euv_response('2024-05-14 02:00:00', 'aia', /evenorm, /chiantifix)
 Nchannels=n_elements(r.channels)
 NT_rsp=n_elements(r.logte)
 logTe_rsp=float(r.logte)
 response=double(r.all)
 dS_rsp=double(r.pix_arcsec)^2
 
 ;load a DEM:
 restore, 'C:\ssw\packages\gx_simulator\euv\ebtel\ebtel_scale=0.2_alpha=-2.5.sav'
 s=size(DEM_cor_run, /dimensions)
 NT_DEM=s[0]
 NQ=s[1]
 NL=s[2]
 
 ;some values for the test:
 Npix=4L
 Nz=10L
 dS_map=2d0^2
 TRfactor=0.5d0
 T0=1d6
 n0=1d8
 
 Lparms_M=[Npix, Nz, Nchannels, NT_rsp, NQ, NL, NT_DEM]
 
 Rparms=[dS_map, dS_rsp]
 
 Rparms_M=dblarr(2, Npix)
 for i=0, Npix-1 do Rparms_M[*, i]=Rparms
 
 Parms=dblarr(10, Nz)
 Parms[0, *]=1000d5 ;voxel length
 Parms[1, *]=T0 ;plasma temperature
 Parms[2, *]=n0 ;plasma density
 Parms[3, 0 : 4]=4 ;corona
 Parms[3, 0]=4+2   ;corona & TR
 Parms[4, 0 : 4]=0.0035 ;Q (some typical value)
 Parms[5, 0 : 4]=2.5d8  ;L (some typical value)
 Parms[6, *]=TRfactor ;TRfactor is defined everywhere, but is used only in one voxel
 
 Parms_M=dblarr(10, Nz, Npix)
 for i=0, Npix-1 do begin
  Parms_M[*, *, i]=Parms
  Parms_M[4, *, i]=Parms[4, *]*10d0^(0.1*i)
  Parms_M[5, *, i]=Parms[5, *]*10d0^(0.1*i)
 endfor
 
 Parms_M[3, 0, 2 : 3]+=8 ;the EUVTR mask is set in last two LOSs (#2 and #3)
 
 flux_M=dblarr(3, Nchannels, Npix)
 
 res=call_external('EUVtransfer_64.dll', 'GET_GX_EUV_SLICE', $
                   Lparms_M, Rparms_M, Parms_M, logTe_rsp, response, $
                   Qrun, Lrun, logtDEM, DEM_cor_run, DEM_tr_run, flux_M, /unload)
                   
 for i=0, Npix-1 do print, transpose(flux_M[*, *, i])     
end