pro Example_EUV_GX_SingleThread
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
 Nz=10L
 dS_map=2d0^2
 TRfactor=0.5d0
 T0=1d6
 n0=1d8
 
 Lparms=[Nz, Nchannels, NT_rsp, NQ, NL, NT_DEM]
 
 Rparms=[dS_map, dS_rsp]
 
 Parms=dblarr(10, Nz)
 Parms[0, *]=1000d5 ;voxel length
 Parms[1, *]=T0 ;plasma temperature
 Parms[2, *]=n0 ;plasma density
 Parms[3, 5 : 8]=1 ;corona
 Parms[3, 9]=1+2   ;corona & TR
 Parms[4, 5 : 9]=0.0035 ;Q (some typical value)
 Parms[5, 5 : 9]=2.5d8  ;L (some typical value)
 Parms[6, *]=TRfactor ;TRfactor is defined everywhere, but is used only in one voxel
 
 flux=dblarr(2, Nchannels)
 
 res=call_external('EUVtransfer_64.dll', 'GET_GX_EUV', $
                   Lparms, Rparms, Parms, logTe_rsp, response, $
                   Qrun, Lrun, logtDEM, DEM_cor_run, DEM_tr_run, flux, /unload)
                   
 print, transpose(flux)     
end