pro extract_nexrad,dataDir,outputDir,timeWildCard,limits,utmZone,StandardTimeZone,decompress=decompress


;++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;Example of parameter values
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++

;dataDir='~/rawData/nexradHourlyPrecipitation/CBRFC/stageIII/' ;'dataGila/nexrad\ precipitation\ hourly/CBRFC/stageIII/'
;outputDir='/tmp/outputRain/'

;timeWildCard='2020-0[789]'

;limits=[595887.5, 3592657,871943.5,3818251.0]

;utmZone=12
;StandardTimeZone=-7

;decompress=0

;++++++++++++++++++++++++++++++++++++++++++++++++++++++++

close,1
cd,dataDir
dirs_time=file_search(timeWildCard,/test_directory)

;++++++++++++++++++++++++++++
;DECOMPRESS TAR FILES

if not(keyword_set(decompress)) then goto,jumpDecompress

;++++++++++++++++++++++++++++
;DECOMPRESS TAR FILES

dirs_time=file_search(/test_directory)
for j=0,n_elements(dirs_time)-1 do begin
	cd,dirs_time[j]
	files_to_uncomp=file_search('*.tar')
	if files_to_uncomp[0] ne '' then begin
		for i=0,n_elements(files_to_uncomp)-1 do begin
			spawn,'tar xvf '+files_to_uncomp[i]
			spawn,'rm '+files_to_uncomp[i]
		endfor
	endif

	spawn,'gunzip *'
	files_to_uncomp=file_search('*.tar*')
	if files_to_uncomp[0] ne '' then begin
		for i=0,n_elements(files_to_uncomp)-1 do begin
 		       	spawn,'tar xvf '+files_to_uncomp[i]
        		spawn,'rm '+files_to_uncomp[i]
		endfor
	endif
	spawn,'gunzip *'
	cd,'..'
endfor
cd,'../..'

;++++++++++++++++++++++++++++

jumpDecompress:

window,1,xs=600,ys=600,xp=0,yp=70
window,10,xs=900,ys=300

totalVolumes=0.0

ncol=0L
nrow=0.0D
	
xll=0.0D
yll=0.0D

for j=0,n_elements(dirs_time)-1 do begin        
	cd,dirs_time[j]
	files_to_checkout=file_search('*')

	for i=0,n_elements(files_to_checkout)-1 do begin
				
		breakName=strsplit(files_to_checkout[i],'xmrgz_CB',/extract)
		if n_elements(breakName) eq 1 then begin
			fileMM=strmid(breakName[0],0,2)
			fileDD=strmid(breakName[0],2,2)
			if strlen(breakName[0]) eq 8 then begin
				fileYY=fix(strmid(breakName[0],4,2))+1900
				fileHH=fix(strmid(breakName[0],6,2))
			endif else begin
				fileYY=fix(strmid(breakName[0],4,4))
				fileHH=fix(strmid(breakName[0],8,2))
			endelse
		endif else begin
			fileMM=strmid(breakName[0],0,2)
			fileDD=strmid(breakName[0],2,2)
			fileYY=fix(strmid(breakName[0],4,4))
			fileHH=breakName[1]
		endelse

		print,files_to_checkout[i],"	",fileMM,"-",fileDD,"-",fileYY,"-",fileHH

		fileDate=julday(fileMM,fileDD,fileYY,fileHH)
		fileCurrentTimeZoneDate=julday(fileMM,fileDD,fileYY,fileHH+standardTimeZone)
		CALDAT, fileCurrentTimeZoneDate, staMonth, staDay, staYear, staHour, staMinute, staSecond

		openr,1,files_to_checkout[i],/xdr
		if not eof(1) then begin
			infoFile=fstat(1)
			close,1
			openr,1,files_to_checkout[i]
		        tmp1=intarr(infoFile.size/2)
		        readu,1,tmp1
		        close,1
	
			testOrder=tmp1[0] ne 16
	
			if testOrder then begin
				byteorder,tmp1,/SWAP_IF_LITTLE_ENDIAN
			endif
	
			ncol_h=tmp1[6+testOrder]
	       	        nrow_h=double(tmp1[8+testOrder])
	
       		        xll_h=double(tmp1[2+testOrder])
       	        	yll_h=double(tmp1[4+testOrder])
	
			data=(reform(tmp1[n_elements(tmp1)-(ncol_h+4)*nrow_h:*],ncol_h+4,nrow_h))[2:ncol_h+1,*]>0
		endif else begin
			data=intarr(ncol,nrow)
			close,1
		endelse
		data=data/100.0

		if (j eq 0 and i eq 0) or ncol_h ne ncol or nrow_h ne nrow or xll_h ne xll or yll_h ne yll then begin

			ncol=tmp1[6+testOrder]
	       	        nrow=double(tmp1[8+testOrder])
	
       		        xll=double(tmp1[2+testOrder])
       	        	yll=double(tmp1[4+testOrder])
	
			lat0 = 60.0*!dtor
			lon0 = -105.0*!dtor;
			r0=6372000.0D
		
			hrapx=lindgen(1,ncol*nrow) mod ncol+xll
			hrapy=lindgen(1,ncol*nrow)/ncol+yll
		
			stereoX = 4762.5D * (hrapx - 401D);
			stereoY = 4762.5D * (hrapy - 1601D);
		      
			R = sqrt( (stereoX*stereoX) + (stereoY*stereoY) );
			
			ThetaComponent = (CV_COORD( FROM_RECT=[stereoX,stereoY] , /TO_POLAR))[0,*]
			
			lonprime = -90.0 - (-105.0) - ThetaComponent/!dtor;
			latprime = 90.0 - 2*( atan( R/(r0*(1+sin(lat0))) )/!dtor )
	
			utmArray=[-lonprime,latprime]
		
			openw,1,'/tmp/coordsToConvert.txt'
			printf,1,utmArray
			close,1

			spawn,'proj +proj=utm  +zone='+strtrim(utmZone,2)+' +ellps=WGS84 /tmp/coordsToConvert.txt > /tmp/coordsConverted.txt'

			openr,1,'/tmp/coordsConverted.txt'
			readf,1,utmArray
			close,1

			spawn,'rm /tmp/coordsConverted.txt'

			!p.multi=0

			window,10,xs=900,ys=300
			!p.multi=[0,3,1]
			plot,hrapx,hrapy,/ynozero,title='HARP',/isotropic
			plot,-lonprime,latprime,/ynozero,title='Lat-Long',/isotropic
			plot,utmArray[0,*],utmArray[1,*],/ynozero,title='UTM',/isotropic
			!p.multi=0

		endif
		
		;Toufique
		;XregInterest=[327633-10000,377961+10000]
		;YregInterest=[3945451-10000,4000078+10000]

		;Sevilleta
		;XregInterest=[345711-10000,375711+10000]
		;YregInterest=[3795105-10000,3795105+10000]

		;Generic
		XregInterest=[limits[0],limits[2]]
		YregInterest=[limits[1],limits[3]]

		elemInSubSquare=where(utmArray[0,*] ge XregInterest[0] and utmArray[0,*] le XregInterest[1] and utmArray[1,*] ge YregInterest[0] and utmArray[1,*] le YregInterest[1])
		bigRegData=data[elemInSubSquare]

		outputFileName=strmid(strtrim(staMonth/100.0,2),2,2)+strmid(strtrim(staDay/100.0,2),2,2)+strtrim(staYear,2)+strmid(strtrim(staHour/100.,2),2,2)

		if max(bigRegData) gt 0 then begin
			wset,1
			!p.multi=0

			contour,data[elemInSubSquare],utmArray[0,elemInSubSquare],utmArray[1,elemInSubSquare],/fill,$
				/irregular,title=files_to_checkout[i],nlevels=32,$
				xrange=[XregInterest[0],XregInterest[1]],yrange=[YregInterest[0],YregInterest[1]],$
				/isotropic,xstyle=1,ystyle=1
			oplot,utmArray[0,elemInSubSquare],utmArray[1,elemInSubSquare],thick=2,psym=3,color=150
			
			XEndSide=XregInterest[0]+long((XregInterest[1]-XregInterest[0])/4000)*4000
			YEndSide=YregInterest[0]+long((YregInterest[1]-YregInterest[0])/4000)*4000

			out_nCols=long((-XregInterest[0]+XEndSide)/4000)+1
			out_nRows=long((-YregInterest[0]+YEndSide)/4000)+1

			TRIANGULATE, utmArray[0,elemInSubSquare],utmArray[1,elemInSubSquare], tr, b  
			dataOutput=TRIGRID(utmArray[0,elemInSubSquare],utmArray[1,elemInSubSquare],data[elemInSubSquare], tr,[4000,4000],[XregInterest[0],YregInterest[0],XEndSide,YEndSide])

			openw,1,outputDir+'/p'+outputFileName+'.txt'

			printf,1,"ncols		"+strtrim(out_nCols,2)
			printf,1,"nrows		"+strtrim(out_nRows,2)
			printf,1,"xllcorner	"+strtrim(XregInterest[0],2)
			printf,1,"yllcorner	"+strtrim(YregInterest[0],2)
			printf,1,"cellsize	"+strtrim(4000,2)
			printf,1,"NODATA_value	-9999"
			printf,1,reverse(dataOutput,2),format='('+strtrim(out_nCols,2)+'F10.4)'
			
			close,1
			print,files_to_checkout[i],'>>>>',max(dataOutput)

			!p.multi=0

		endif else begin
			wset,1
			!p.multi=0
			plot,utmArray[0,elemInSubSquare],utmArray[1,elemInSubSquare],psym=3,/ynozero,color=150,$
				xrange=[XregInterest[0],XregInterest[1]],yrange=[YregInterest[0],YregInterest[1]],$
				/isotropic,xstyle=1,ystyle=1,title=files_to_checkout[i]
			!p.multi=0

			XEndSide=XregInterest[0]+long((XregInterest[1]-XregInterest[0])/4000)*4000
			YEndSide=YregInterest[0]+long((YregInterest[1]-YregInterest[0])/4000)*4000

			out_nCols=long((-XregInterest[0]+XEndSide)/4000)+1
			out_nRows=long((-YregInterest[0]+YEndSide)/4000)+1
			
			openw,1,outputDir+'/p'+outputFileName+'.txt'

			printf,1,"ncols		"+strtrim(out_nCols,2)
			printf,1,"nrows		"+strtrim(out_nRows,2)
			printf,1,"xllcorner	"+strtrim(XregInterest[0],2)
			printf,1,"yllcorner	"+strtrim(YregInterest[0],2)
			printf,1,"cellsize	"+strtrim(4000,2)
			printf,1,"NODATA_value	-9999"
			printf,1,intarr(out_nCols,out_nRows),format='('+strtrim(out_nCols,2)+'I5)'
			
			close,1

		endelse

	endfor
	;print,totalVolumes
	cd,'..'
endfor

end
