;+
; NAME: cross_spec_tplot.pro
; SYNTAX: 
; PURPOSE: Returns the cross-phase,coherence and power spectra of two
;          tplot variables over the interval from T1 to T2.
; INPUT: var1_name, var2_name -> names of the two input tplot
;          timeseries variables
;        T1, T2 --> start and end times used for calculation
;        sub_interval --> number of sub intervals
;        overlap_index --> number of the overlap index (3-4 usually
;        good)
;
; OUTPUT: tplot variables of cross-phase, coherence, power spectra
; KEYWORDS: 
; HISTORY: Written by Lei Dai (c. 2012). Modified by Aaron W Breneman
; VERSION: 
;   $LastChangedBy: $
;   $LastChangedDate: $
;   $LastChangedRevision: $
;   $URL: $
;-


;;------------------------------------------------------------------------------
;; FFT of a field
;;------------------------------------------------------------------------------

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


;;Number of Overlap spectrums.
  New_Nspec=lapping_index*(Nspec-1)+1
;;Point_num=n_elements(field)
  nfft=LONG64(Point_num/Nspec)
  RAW_spec=make_array(nfft,New_Nspec,/dcomplex)
  spec=make_array(nfft/2.0+1,New_Nspec,/dcomplex)
  XXX=make_array(nfft,New_Nspec,/double)

  for ii=0,New_Nspec-1  do Begin
     xxx[*,ii]=Field[(ii*nfft/lapping_index):(ii*nfft/lapping_index)+nfft-1]
  endfor

  for nn=0,New_Nspec-1 do begin
     RAW_spec[*,nn]=fft(hanning(nfft)*(xxx[*,nn]-mean(xxx[*,nn])))
     spec[*,nn]=Raw_spec[0:nfft/2.0,nn]
  endfor

  Return,spec
end

;;------------------------------------------------------------------------------
;; Phase spectrum of two arbitary components
;;------------------------------------------------------------------------------

FUNCTION Phase_spectra,FFT_field_1,FFT_field_2,seconds,Nspec,point_num_1,$
                       point_num_2,pow_1,pow_2,lapping_index

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

;;------------------------------------------------------------------------------
;;Average Power Spectrum Density (PSD)
;;------------------------------------------------------------------------------

function Pow_Spectra,FFT_field,seconds,Nspec,point_num,lapping_index
  New_Nspec=lapping_index*(Nspec-1)+1
  Nfft=LONG64(Point_num/Nspec)
  Pspecaverage=make_array(nfft/2.0+1,/double,value=0)
  for mm=0,New_Nspec-1 do $
     Pspecaverage=Pspecaverage+ABS(FFT_field[*,mm])^2/float(New_Nspec)
  Pspecaverage=Pspecaverage*2*seconds/float(nspec)
  Return,Pspecaverage
end

;;------------------------------------------------------------------------------
;;Coherence spectrum of two arbitary components,cross correlation function
;;------------------------------------------------------------------------------

FUNCTION Coherence_spectra,FFT_field_1,FFT_field_2,seconds,Nspec,point_num_1,$
                           point_num_2,pow_1,pow_2,lapping_index

  New_Nspec=lapping_index*(Nspec-1)+1
  nfft=LONG64(MIN([point_num_1,point_num_2])/float(nspec))
  Correlation=make_array(nfft/2.0+1,/double,value=0)
  G_1_2=make_array(nfft/2.0+1,/dcomplex,value=0)

  for nn=0,New_Nspec-1 do begin
     G_1_2=G_1_2+(2*seconds/nspec)*(FFT_field_1[*,nn]*Conj(FFT_field_2[*,nn]))/Float(New_Nspec)
  Endfor

;Alternative defination of coherence
  Correlation=ABS(G_1_2)/(pow_1*pow_2)^0.5 
;Correlation=ABS(G_1_2)^2/(pow_1*pow_2)
  return,correlation
END



;;------------------------------------------------------------------------------
;; Main function
;;------------------------------------------------------------------------------

Function cross_spec_tplot,var1_name,var1_y_dim, var2_name, var2_y_dim, T1,T2,$
                          sub_interval=sub_interval, overlap_index=overlap_index

  if not keyword_set(sub_interval) then sub_interval = 1.0
  if not keyword_set(overlap_index) then overlap_index = 1.0
  Nspec=sub_interval
  lapping_index=overlap_index

  get_data, var1_name, data = d
  d_temp=time_trim([[d.x],[d.y[*,var1_y_dim]]],time_double(T1),time_double(T2))
  times1 = d_temp[*,0]
  values1 = reform(d_temp[*,1:n_elements(var1_y_dim)])

  get_data, var2_name, data = d
  times2 = d.X
  values2 = reform(d.Y[*,var2_y_dim])

  values1_interp=interpol(values1,times1,times1)
  Values2_interp=interpol(values2,times2,times1)
  times=times1
  E_num=Long64(n_elements(times))
  Enfft=E_num/sub_interval
  Seconds=double(max(times,/nan)-Min(times,/nan)) 
  E_FREQUENCE_coordinate=Nspec*findgen(Enfft/2+1)/seconds

  FFT_E_1=FFT_series(values1_interp,seconds,Nspec,E_num,lapping_index)
  Aver_Pow_E_1=Pow_Spectra(FFT_E_1,seconds,Nspec,E_num,lapping_index)
  FFT_E_2=FFT_series(values2_interp,seconds,Nspec,E_num,lapping_index)
  Aver_Pow_E_2=Pow_Spectra(FFT_E_2,seconds,Nspec,E_num,lapping_index)
  COHERENCE_E1_E2=COHERENCE_spectra(FFT_E_1,FFT_E_2,seconds,Nspec,E_num,E_num,Aver_Pow_E_1,Aver_Pow_E_2,lapping_index)
  PHASE_E1_E2=Phase_spectra(FFT_E_1,FFT_E_2,seconds,Nspec,E_num,E_num,Aver_Pow_E_1,Aver_Pow_E_2,lapping_index)
  return,[[E_FREQUENCE_coordinate],[PHASE_E1_E2],[COHERENCE_E1_E2],[Aver_Pow_E_1],[Aver_Pow_E_2]]

end








