FUNCTION Poly_fit_cal,coe,indep_variable
  order=n_elements(coe)
  Fit_results=indep_variable & Fit_results[*]=0
  For i=0, order-1 DO Begin
     temp=coe[i]*indep_variable[*]^i
     Fit_results=fit_results+temp
  endfor
  Return,fit_results
END

;FFT of a field
FUNCTION FFT_series,field,seconds,Nspec,Point_Num,lapping_index,zero_pad=zp


;zp = 0
;if keyword_set(zp) then begin
;	;Zero-pad these arrays to speed up FFT
;	fac = 1
;	while 2L^fac lt n_elements(field) do fac++
;	addarr = fltarr(2L^fac - n_elements(field))  ;array of zeros
;;stop
;	field = [field,addarr]
;	Point_Num = 2L^fac
;endif


;Number of Overlap spectrums.
  New_Nspec=lapping_index*(Nspec-1)+1
;Point_num=n_elements(field)
  nfft=LONG64(Point_num/Nspec)
  RAW_spec=make_array(nfft,New_Nspec,/dcomplex)
  spec=make_array(nfft/2.0+1,New_Nspec,/dcomplex)
  XXX=make_array(nfft,New_Nspec,/double)
  for ii=0,New_Nspec-1  do Begin
     xxx[*,ii]=Field[(ii*nfft/lapping_index):(ii*nfft/lapping_index)+nfft-1]
  ENDFOR



  for nn=0,New_Nspec-1 do begin

     RAW_spec[*,nn]=fft(hanning(nfft)*(xxx[*,nn]-mean(xxx[*,nn])))
     spec[*,nn]=Raw_spec[0:nfft/2.0,nn]
  EndFOR
  Return,spec
END

;Phase spectrum of two arbitary components
FUNCTION Phase_spectra,FFT_field_1,FFT_field_2,seconds,Nspec,point_num_1,point_num_2,pow_1,pow_2,lapping_index
  New_Nspec=lapping_index*(Nspec-1)+1
  nfft=LONG64(MIN([point_num_1,point_num_2])/float(nspec))
  Phase=make_array(nfft/2.0+1,/double,value=0)
  G_1_2=make_array(nfft/2.0+1,/dcomplex,value=0)
  for nn=0,New_Nspec-1 do begin
     G_1_2=G_1_2+(2.0*seconds/nspec)*(FFT_field_1[*,nn]*Conj(FFT_field_2[*,nn]))/Float(New_Nspec)
  Endfor
  Phase=ATAN(Imaginary((G_1_2)/(pow_1*pow_2)),double((G_1_2)/(pow_1*pow_2)))
  return,Phase
END

;Average Power Spectrum Density,PSD
FUNCTION Pow_Spectra,FFT_field,seconds,Nspec,point_num,lapping_index
  New_Nspec=lapping_index*(Nspec-1)+1
  Nfft=LONG64(Point_num/Nspec)
  Pspecaverage=make_array(nfft/2.0+1,/double,value=0)
  for mm=0,New_Nspec-1 do Pspecaverage=Pspecaverage+ABS(FFT_field[*,mm])^2/float(New_Nspec)
  Pspecaverage=Pspecaverage*2*seconds/float(nspec)
  Return,Pspecaverage
END

;Coherence spectrum of two arbitary components,cross correlation function
FUNCTION Coherence_spectra,FFT_field_1,FFT_field_2,seconds,Nspec,point_num_1,point_num_2,pow_1,pow_2,lapping_index
  New_Nspec=lapping_index*(Nspec-1)+1
  nfft=LONG64(MIN([point_num_1,point_num_2])/float(nspec))
  Correlation=make_array(nfft/2.0+1,/double,value=0)
  G_1_2=make_array(nfft/2.0+1,/dcomplex,value=0)
  for nn=0,New_Nspec-1 do begin
     G_1_2=G_1_2+(2*seconds/nspec)*(FFT_field_1[*,nn]*Conj(FFT_field_2[*,nn]))/Float(New_Nspec)
  Endfor
  Correlation=ABS(G_1_2)/(pow_1*pow_2)^0.5
  return,correlation
END

;***********************Main function ***********************************************************************************************

;Example, dynamic_cross_spec_tplot,'rbspb_efw_e-spinfit-mgse_e12_spinfit_mgse',1,'rbspb_efw_e-spinfit-mgse_e12_spinfit_mgse',2,T1,T2,1000.0,100.,1500.0,new_name='rbspb_Ey_Ez'
;window: time window used for evaluation of fourier spectra. lag: the time lag between slicing window.. coherence time: the interval over which coherence is evaluated. must larger than window.
Pro dynamic_cross_spec_tplot,var1_name,var1_y_dim,var2_name,var2_y_dim,T1,T2,window,lag,coherence_time,new_name=new_name

  get_data, var1_name, data = d
  d_temp=time_trim([[d.x],[d.y[*,var1_y_dim]]],time_double(T1),time_double(T2))
;d_temp = tsample(var1_name,time_double([T1,T2]),times=tms)

;time_array = tms

  time_array = d_temp[*,0]
;values1 = reform(d_temp)
;values1_interp = values1

  values1 = reform(d_temp[*,1:n_elements(var1_y_dim)])
  values1_interp=interpol(values1,time_array,time_array)
  E_num=Long64(n_elements(time_array))

;***Original method. This isn't very robust and can be way wrong.
;  data_rate=1.0/mean(time_array[2:12]-time_array[1:11],/nan)

;***New method...much more robust.
  boogoo = rbsp_sample_rate(time_array,out_med_avg=medavg,average=avg)
  data_rate = medavg[0]



;get_data, var2_name, data = d
;d_temp = tsample(var2_name,time_double([T1,T2]),times=tms)
;time_array = tms
;;time_array = d_temp[*,0]
;values2 = reform(d_temp)
;values2_interp = values1


  get_data, var2_name, data = d2
  d2_temp=time_trim([[d2.x],[d2.y[*,var2_y_dim]]],time_double(T1),time_double(T2))
  values2 = reform(d2_temp[*,1:n_elements(var2_y_dim)])
  values2_interp=interpol(values2,d2_temp[*,0],time_array)


  Nspec=1
  lapping_index=1
;calculate frequency coordinate
  Enfft=window*data_rate


;zero_pad = 0
;if keyword_set(zero_pad) then begin
;	fac = 1
;	while 2L^fac lt Enfft do fac++
;	Enfft = 2L^fac
;endif



  E_FREQUENCE_coordinate=Nspec*findgen(Enfft/2+1)/window



  New_time_array=!values.f_nan
  dynamic_FFT_1=make_array(n_elements(E_FREQUENCE_coordinate),/double) & dynamic_FFT_1[*]=!values.f_nan
  dynamic_FFT_2=make_array(n_elements(E_FREQUENCE_coordinate),/double) & dynamic_FFT_2[*]=!values.f_nan

  dynamic_phase=make_array(n_elements(E_FREQUENCE_coordinate),/double) & dynamic_phase[*]=!values.f_nan

;; ii is the number of sliding window.
  coherence_num=(coherence_time-window)/lag/2.0


  print,'***********'
  print,'Enum = ',E_num

  ii=Long64(1)
  While ((ii-1)*(lag)*data_rate+(window)*data_rate LT E_num) Do Begin
     FFT_E_1=FFT_series(values1_interp[(ii-1)*(lag)*data_rate:(ii-1)*(lag)*data_rate+(window)*data_rate-1],window,Nspec,window*data_rate,lapping_index)
     FFT_E_2=FFT_series(values2_interp[(ii-1)*(lag)*data_rate:(ii-1)*(lag)*data_rate+(window)*data_rate-1],window,Nspec,window*data_rate,lapping_index)
     Phase_FFT_E_12=ATAN(Imaginary(FFT_E_1/FFT_E_2),double(FFT_E_1/FFT_E_2))*180.0/3.1415

     New_time_array=[New_time_array,mean(time_array[(ii-1)*(lag)*data_rate:(ii-1)*(lag)*data_rate+(window)*data_rate-1],/nan)]
     dynamic_FFT_1=[[dynamic_FFT_1],[FFT_E_1]]
     dynamic_FFT_2=[[dynamic_FFT_2],[FFT_E_2]]

;dynamic_FFT_1=[dynamic_FFT_1,FFT_E_1]
;dynamic_FFT_2=[dynamic_FFT_2,FFT_E_2]



     if ii mod 100 eq 0 then print,'...% finished', 100.*n_elements(dynamic_FFT_1)/E_num
;help,dynamic_FFT_1
     dynamic_phase=[[dynamic_phase],[Phase_FFT_E_12]]
     ii=ii+1
  ENDwhile

  store_data,'dynamic_FFT_1',data={x:new_time_array,y:transpose(real_part(dynamic_FFT_1)),V:E_FREQUENCE_coordinate,spec:1},lim=lim_spec
  store_data,'dynamic_FFT_2',data={x:new_time_array,y:transpose(real_part(dynamic_FFT_2)),V:E_FREQUENCE_coordinate,spec:1},lim=lim_spec


  dynamic_coherence=make_array(n_elements(E_FREQUENCE_coordinate),n_elements(New_time_array),/double) & dynamic_coherence[*]=!values.f_nan
  phase_sum=make_array(n_elements(E_FREQUENCE_coordinate),n_elements(New_time_array),/double)
  power_sum=make_array(n_elements(E_FREQUENCE_coordinate),n_elements(New_time_array),/double)

  For jj=coherence_num, n_elements(New_time_array)-1-coherence_num Do Begin
     For kk=jj-coherence_num, jj+coherence_num Do phase_sum[*,jj]=phase_sum[*,jj]+dynamic_FFT_1[*,kk]*Conj(dynamic_FFT_2[*,kk])
     For kk=jj-coherence_num, jj+coherence_num Do power_sum[*,jj]=power_sum[*,jj]+ABS(dynamic_FFT_1[*,kk])*Abs(dynamic_FFT_2[*,kk])
     dynamic_coherence[*,jj]= ABS(phase_sum[*,jj])/power_sum[*,jj]
  ENDFOR

  lim_spec={YSTYLE:1,PANEL_SIZE:1,XMINOR:5,XTICKLEN:0.04,$
            ZSTYLE:1,ZLOG:0}
  name=var1_name
  if not keyword_set(new_name) then new_name = name+'_spectrum'

;Store_data,new_name+'_power',data={x:new_time_array,y:transpose(dynamic_power_1),V:E_FREQUENCE_coordinate,spec:1},lim=lim_spec
  Store_data,new_name+'_phase',data={x:new_time_array,y:transpose(dynamic_phase),V:E_FREQUENCE_coordinate,spec:1},lim=lim_spec
  Store_data,new_name+'_coherence',data={x:new_time_array,y:transpose(dynamic_coherence),V:E_FREQUENCE_coordinate,spec:1},lim=lim_spec
  zlim,new_name+'_coherence',0,1
  zlim,new_name+'_phase',-180,180,0

  loadct,39
  device,decompose=0
  A=STRMID(STRTRIM(window,2),0,4)
  B=STRMID(STRTRIM(lag,2),0,4)
  C=STRMID(STRTRIM(data_rate,2),0,5)
  D=STRMID(STRTRIM(coherence_time,2),0,5)

;options,new_name+'_phase','title','window='+A+',lag='+B+',data_rate='+C+'/s'
;options,new_name+'_coherence','title','window='+A+',lag='+B+',data_rate='+C+'/s'+', coherence '+D+'s'

END
