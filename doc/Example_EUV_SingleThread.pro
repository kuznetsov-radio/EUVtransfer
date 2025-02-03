pro Example_EUV_SingleThread
 ;load an EUV response matrix: 
 r=gx_euv_response('2024-05-14 02:00:00', 'aia', /evenorm, /chiantifix)
 Nchannels=n_elements(r.channels)
 NT_rsp=n_elements(r.logte)
 logTe_rsp=double(r.logte)
 response=double(r.all)
 dS_rsp=double(r.pix_arcsec)^2
 
 ;load a DEM:
 restore, 'C:\ssw\packages\gx_simulator\euv\ebtel\ebtel_scale=0.2_alpha=-2.5.sav'
 NT_DEM=n_elements(logtdem)
 logTe_DEM=double(logtdem)
 DEM_cor=double(reform(dem_cor_run[*, 30, 25])) ;some typical DEM
 DEM_tr=double(reform(dem_tr_run[*, 30, 25])) ;some typical DEM
 
 ;some values for the test:
 Nz=10L
 dS_map=2d0^2
 TRfactor=0.5d0
 T0=1d6
 n0=1d8
 
 Lparms=[Nz, Nchannels, NT_rsp, NT_DEM]
 
 Rparms=[dS_map, dS_rsp, TRfactor]
 
 Parms=dblarr(6, Nz)
 Parms[0, *]=1000d5 ;voxel length
 Parms[1, *]=T0 ;plasma temperature
 Parms[2, *]=n0 ;plasma density
 Parms[3, 0 : 4]=1; DEM_on in first 5 voxels
 
 DEM_cor_arr=dblarr(NT_DEM, Nz)
 for i=0, 4 do DEM_cor_arr[*, i]=DEM_cor ;DEM is defined in first 5 voxels
 
 flux=dblarr(2, Nchannels)
 
 res=call_external('EUVtransfer_64.dll', 'GET_EUV', $
                   Lparms, Rparms, Parms, logTe_rsp, response, logTe_DEM, DEM_cor_arr, DEM_tr, flux)
                   
 print, transpose(flux)     
end