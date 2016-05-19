


path = '/Users/aaronbreneman/Desktop/Research/OTHER/Stuff_for_other_people/Thaller_Scott/'


restore,filename=path + 'F10_pt_7_2012_2016.dat'
store_data,'F10.7',data={x:times,y:ftps}

restore,filename=path + 'rbspa_outbnd_pp.dat'
store_data,'rbspa!Coutbound!Cplasmapause',data={x:times,y:obpp}

restore,filename=path + 'rbspa_inbnd_pp.dat'
store_data,'rbspa!Cinbound!Cplasmapause',data={x:times,y:ibpp}

rbsp_detrend,['rbspa!Coutbound!Cplasmapause','rbspa!Cinbound!Cplasmapause'],86400.*3


;works
;ndays = 200.
;window = ndays*86400.
;lag = window/4.
;coherence_time = window*2

;works
ndays = 250.
window = ndays*86400.
lag = window/4.
coherence_time = window*2



     v1 = 'F10.7'
     v2 = 'rbspa!Cinbound!Cplasmapause'

     TT1 = time_double('2012-10-03')
     TT2 = time_double('2015-11-04')
     dynamic_cross_spec_tplot,v1,0,v2,0,TT1,TT2,window,lag,coherence_time,$
                              new_name='chorus_mb'  

     cormin = 0.6
     get_data,'chorus_mb_coherence',data=coh
     ;;get rid of first element which has t=0
     coh = {x:coh.x[1:n_elements(coh.x)-1],y:coh.y[1:n_elements(coh.x)-1,*],v:coh.v,spec:1}
     goo = where(coh.y lt cormin)
     if goo[0] ne -1 then coh.y[goo] = !values.f_nan
     store_data,'chorus_mb_coherence',data={x:coh.x,y:coh.y,v:1/coh.v/86400.}
     options,'chorus_mb_coherence','ytitle','Chorus, MB!Ccoherence!C[period (days)]'
     copy_data,'chorus_mb_coherence','chorus_mb_coherence_log'
     options,'chorus_mb_coherence*','spec',1


	get_data,'dynamic_FFT_1',data=fft1
     fft1 = {x:fft1.x[1:n_elements(fft1.x)-1],y:fft1.y[1:n_elements(fft1.x)-1,*],v:fft1.v,spec:1}
     store_data,'fft1',data={x:fft1.x,y:abs(fft1.y),v:1/fft1.v/86400.}
     options,'fft1','spec',1
	get_data,'dynamic_FFT_2',data=fft1
     fft1 = {x:fft1.x[1:n_elements(fft1.x)-1],y:fft1.y[1:n_elements(fft1.x)-1,*],v:fft1.v,spec:1}
     store_data,'fft2',data={x:fft1.x,y:abs(fft1.y),v:1/fft1.v/86400.}
     options,'fft2','spec',1

	zlim,'fft1',0.1,10,1
	zlim,'fft2',0.001,0.1,1

	linev = replicate(27.,n_elements(coh.x))
	store_data,'line',data={x:coh.x,y:linev}
	store_data,'coh',data=['chorus_mb_coherence','line']
	store_data,'fft11',data=['fft1','line']
	store_data,'fft22',data=['fft2','line']
	ylim,['fft11','fft22','coh'],0,60,0
	zlim,'coh',0.6,0.9,0
	ylim,'rbspa!Coutbound!Cplasmapause_smoothed',3,6
	ylim,'rbspa!Cinbound!Cplasmapause_smoothed',3,6

	tplot,['F10.7','fft11',$
		'rbspa!Coutbound!Cplasmapause_smoothed','fft22',$
			'coh']


	tplot,['F10.7','fft11',$
		'rbspa!Cinbound!Cplasmapause_smoothed','fft22',$
			'coh']






