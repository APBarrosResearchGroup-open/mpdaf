;Generated MSHM forcing using HRRR data and NLDAS Albedo
; HRRR and NLDAS data are first interpolated to the same regular
; geographic grid with fixed dx and dy (kms)

function gapfill(var,xo)
begin
 nx = dimsizes(xo)
 xi = (/xo(0)-1,xo(nx-1)+1/)
 ndx = dimsizes(var)
 do ip = 0, ndx(1)-1
   fi = (/var(xi(0),ip),var(xi(1),ip)/)
   var(xo(0):xo(nx-1),ip) = linint1_Wrap (xi, fi, False, xo,0)
 end do
 return(var)
end

function extract_nldas(fil)
begin
 nc_fil = addfile(fil,"r")
 lat = nc_fil->lat
 lon = nc_fil->lon
 ny = dimsizes(lat)
 nx = dimsizes(lon)
 np = ny*nx
 var = new((/np/),"float")
 var = ndtooned(nc_fil->Band1)
 var = var/100. ;% to fraction
 return(var)
end

function extract_hrrr(fil,vname)
begin
  nc_fil = addfile(fil,"r")
  lat = nc_fil->lat
  lon = nc_fil->lon
  ny = dimsizes(lat)
  nx = dimsizes(lon)
  np = ny*nx
  var = new((/np/),"float")
 
  if (vname.eq."ws") then
    uu        = ndtooned(nc_fil->uu(0,:,:))   
    vv        = ndtooned(nc_fil->vv(0,:,:))
    var       = sqrt(uu^2 + vv^2)    
  else
    var = ndtooned(nc_fil->$(vname)$)
  end if
  return(var)
end
;;;;;

begin
 diri=getenv("DATA")
 scratch=getenv("SCRATCH")
 intpdir = scratch + "/HRRR_90m"
 ;;
 wyr = 2020
 mnths = (/10,11,12,1,2,3,4,5,6,7,8,9/)
 band = "X"

 ntime = 8784; 5832
 ctr=-1
 mctr=-1
 do ii = 0, dimsizes(mnths)-1
 mm = mnths(ii)
 if (ii.le.2) then
   yyyy = wyr-1
 else
   yyyy = wyr
 end if
 yyyymm = sprinti("%0.4i",yyyy) + sprinti("%0.2i",mm)
 fdir = intpdir +  "/interp_" + yyyymm
 do dd = 1,days_in_month(yyyy,mm)
 yyyymmdd = yyyymm + sprinti("%0.2i",dd)
 do hh = 0,23
   ctr = ctr + 1
   yyyymmddhh = yyyymmdd + sprinti("%0.2i",hh)
   ;;fname = fdir + "/interp_" + yyyymmddhh + "_crop.nc"
   fname = fdir + "/interp_" + yyyymmddhh + ".nc"
   fname_nldas = scratch + "/NLDAS_90m/"  + "/interp_ALBEDO/" + "interp_NLDAS_MOS0125_H.A" + yyyymmdd + "." + sprinti("%0.2i",hh) + "00.020.nc.SUB.nc4" 
   print(" " + fname)
   print(" " + fname_nldas)
   if (ctr.eq.0) then
     val = extract_hrrr(fname,"t2m")  ;K
     ndim = dimsizes(val)
     abh2 = new((/ntime,ndim/),"float")
     tkh2 = new((/ntime,ndim/),"float")
     sph2 = new((/ntime,ndim/),"float")
     pph2 = new((/ntime,ndim/),"float")
     prh2 = new((/ntime,ndim/),"float")
     wsh2 = new((/ntime,ndim/),"float")
     lwh2 = new((/ntime,ndim/),"float")
     swh2 = new((/ntime,ndim/),"float")
   end if

   if isfilepresent(fname) then  ;Some HRRR files are missing 
     tkh2(ctr,:) = extract_hrrr(fname,"t2m")  ;K
     sph2(ctr,:) = extract_hrrr(fname,"sp")  ;kg/kg
     wsh2(ctr,:) = extract_hrrr(fname,"ws")  ;m/s
     pph2(ctr,:) = extract_hrrr(fname,"pp")  ;Pa
     prh2(ctr,:) = extract_hrrr(fname,"pr")  ;m/s
     lwh2(ctr,:) = extract_hrrr(fname,"dl")  ;W/m2
     swh2(ctr,:) = extract_hrrr(fname,"dx")  ;W/m2
   else

     tkh2(ctr,:) = -999.
     sph2(ctr,:) = -999.
     wsh2(ctr,:) = -999.
     pph2(ctr,:) = -999.
     prh2(ctr,:) = 0.
     lwh2(ctr,:) = -999.
     swh2(ctr,:) = -999.
     mctr = mctr + 1
   end if

   if isfilepresent(fname_nldas) then
     abh2(ctr,:) = extract_nldas(fname_nldas) ;-
   else
     abh2(ctr,:) = -999.
   end if

 end do
 end do
 end do

 print("Total data : " + ctr)
 print("Missing HRRR files : " + mctr)
 print(" " )

 print("Gap filling ...")
 ;Gap filling and interpolate to half-hourly data
 mind = ind(tkh2(:,0).eq.-999.)
 print(mind)

 xo = ispan(1265,1266,1)
 tkh2 = gapfill(tkh2,xo) 
 sph2 = gapfill(sph2,xo)
 pph2 = gapfill(pph2,xo)
 lwh2 = gapfill(lwh2,xo)
 swh2 = gapfill(swh2,xo)
 wsh2 = gapfill(wsh2,xo)
 delete(xo)

 delete(mind)
 mind = ind(abh2(:,0).eq.-999.)
 print(mind)

 abh2 = gapfill(abh2,2422)
 abh2 = gapfill(abh2,2488)
 abh2 = gapfill(abh2,2525)
 abh2 = gapfill(abh2,4555)
 abh2 = gapfill(abh2,4687)
 abh2 = gapfill(abh2,7310)

 print(" ")
 print("Albedo Min : " + min(abh2) + " Max : " + max(abh2))
 print("Tair Min : " + min(tkh2) + " Max : " + max(tkh2))
 print("Qair Min : " + min(sph2) + " Max : " + max(sph2))
 print("P  Min : " + min(pph2) + " Max : " + max(pph2))
 print("Rain Min : " + min(prh2) + " Max : " + max(prh2))
 print("LW Min : " + min(lwh2) + " Max : " + max(lwh2))
 print("SW Min : " + min(swh2) + " Max : " + max(swh2))
 print("UX Min : " + min(wsh2) + " Max : " + max(wsh2))

 totrain = dim_sum_n(prh2,0)*3600.
 print("Total Precipitation " + min(totrain) + " " + max(totrain) + " m")

 print("Half-hourly linear interpolation ...")
 xi = ispan(0,ntime-1,1)
 xo = ispan(0,ntime*2-2,1)*0.5
 abh = new((/ntime*2-1,ndim/),"float")
 tkh = new((/ntime*2-1,ndim/),"float")
 sph = new((/ntime*2-1,ndim/),"float")
 pph = new((/ntime*2-1,ndim/),"float")
 prh = new((/ntime*2-1,ndim/),"float")
 wsh = new((/ntime*2-1,ndim/),"float")
 lwh = new((/ntime*2-1,ndim/),"float")
 swh = new((/ntime*2-1,ndim/),"float")

 do ip = 0, ndim-1
   abh(:,ip) = linint1_Wrap (xi,abh2(:,ip), False, xo, 0)
   tkh(:,ip) = linint1_Wrap (xi,tkh2(:,ip), False, xo, 0)
   sph(:,ip) = linint1_Wrap (xi,sph2(:,ip), False, xo, 0)
   pph(:,ip) = linint1_Wrap (xi,pph2(:,ip), False, xo, 0)
   prh(:,ip) = linint1_Wrap (xi,prh2(:,ip), False, xo, 0)
   wsh(:,ip) = linint1_Wrap (xi,wsh2(:,ip), False, xo, 0)
   lwh(:,ip) = linint1_Wrap (xi,lwh2(:,ip), False, xo, 0)
   swh(:,ip) = linint1_Wrap (xi,swh2(:,ip), False, xo, 0)
 end do

;write_matrix ("A-" + band + ".csv", abh2) ;-
; write_matrix("example1.txt", "w", alist, "%s")

 opt = True
 fmt = ndim + "f15.8"
 opt@fout = "A-" + band + ".csv"
 write_matrix (abh, fmt, opt) ; -
 opt@fout = "P-" + band + ".csv"
 write_matrix (pph, fmt, opt) ; Pa
 opt@fout = "T-" + band + ".csv"
 write_matrix (tkh, fmt, opt) ; T
 opt@fout = "J-" + band + ".csv"
 write_matrix (prh, fmt, opt) ; kg/m2/s
 opt@fout = "F-" + band + ".csv"
 write_matrix (wsh, fmt, opt) ; m/s
 opt@fout = "B-" + band + ".csv"
 write_matrix (sph, fmt, opt) ; kg/kg
 opt@fout = "CC-" + band + ".csv"
 write_matrix (lwh, fmt, opt) ; W/m2
 opt@fout = "DC-" + band + ".csv"
 write_matrix (swh, fmt, opt) ; W/m2

end
