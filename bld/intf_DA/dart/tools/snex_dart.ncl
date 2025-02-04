;Convert SNEX snow depths to DART observations
;Text Format
; Obs_type, Lat, Lon, Hgt, Year, Month, Day, H, M, S, Value, Obs_Error
;======================================================
; Prabhakar Shrestha
; E-mail: ps98@illinois.edu
;=======================================================

load "./functions/read_snex.ncl"
load "./functions/read_xml.ncl"

begin

; User Settings
;-----------------------------------------------------------------------------------------

 ssid = (/ 255640259, 255640260, 255640261, 255640262, 255640267/)
 obs_type = 1
 obs_hgt  = 0.0 ;m
 obs_err  = 0.02 ;m
 ;Filter parameters
 year     = 2016
 hour     = 0.

;;; Crop index of DA
 yyyy = 2016
 mm   = 10
 dd   = 21
 hh   = 0
 ndays = 162
 nfreq = 24

; Read Data
;-----------------------------------------------------------------------------------------
;For water year, extend to next year
 year1 = year + 1

 sid  = ssid(id)
 diri = "/home/ps98/data/NASA_SnowEX_5000004096036"
 filname =  systemfunc("ls " + diri + "/" + sid + "/*.csv")
 xml     = systemfunc("ls " + diri + "/" + sid + "/*.xml")

 lonlat = read_xml(xml)
 if (lonlat(1).lt.0) then
   lonlat(1) = 360. + lonlat(1)
 end if

 col11  = read_snex(filname,20)  ;Snow Depth
 time   = col11@cal

 ;Filter data for choosen time of assimilation
 sind = ind(time(0,:).ge.year .and. time(0,:).le.year1 .and. time(4,:).eq.hour .and. col11.ne.-9999)
 tmp1 = col11(sind)
 tmp2 = time(:,sind)
 delete([/col11,time/])
 col11 = (/tmp1/)
 time  = tmp2

 ndim   = dimsizes(col11)

 col1 = new(ndim,"integer")
 col2 = new(ndim,"float")
 col3 = new(ndim,"float")
 col4 = new(ndim,"float")  
 col5 = new(ndim,"integer")
 col6 = new(ndim,"integer")
 col7 = new(ndim,"integer")
 col8 = new(ndim,"integer")
 col9 = new(ndim,"integer")
 col10 = new(ndim,"integer")
 col12 = new(ndim,"float")
 ;
 col1 = obs_type
 col2 = lonlat(0)
 col3 = lonlat(1)
 col4 = obs_hgt
 col5 = time(0,:)
 col6 = time(1,:)
 col7 = time(2,:)
 col8 = time(3,:)
 col9 = time(4,:)
 col10 = time(5,:)
 col12 = obs_err

; alist   = [/col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12/]

 stind = ind(col5.eq.yyyy .and. col6.eq.mm .and. col7.eq.dd .and. col8.eq.hh)
 etind = stind + ndays*nfreq

 do ic = stind, etind, nfreq
   fname   = "data_snex_" + sprinti("%0.4i",col5(ic)) + sprinti("%0.2i",col6(ic)) + \
                sprinti("%0.2i",col7(ic)) + sprinti("%0.2i",col8(ic)) + \
                sprinti("%0.2i",col9(ic)) + sprinti("%0.2i",col10(ic))
   alist   = [/col1(ic),col2(ic),col3(ic),col4(ic),col5(ic),col6(ic),col7(ic),col8(ic),col9(ic),col10(ic),col11(ic),col12(ic)/]
   write_table(fname, "w", alist, "%i %f %f %f %i %i %i %i %i %i %f %f")
 end do

 
end
