;Process HRRR 3km grib2 data downloaded from AWS
;Interpolates data to specified resolution (dlx,dly) in regular geographic
;co-ordinates from curvilinear grids
;This interpolation will not work for downscaling 
;P. Shrestha

function getdgeo(dsx,dsy)
begin
 llat = get_pi("float")*39./180.
 dlat = dsy*(1./119.574)
 dlon = dsx*(1./(111.320*cos(llat)))
 ;print(dlat + " " + dlon)
 dgeo = (/dlat,dlon/)
 return(dgeo)
end

function extract_interp(grb_fil,lat,lon,varname,latx,lonx,opt)
begin
  ;Read data
  varc  = grb_fil->$(varname)$
  lat2d = grb_fil->gridlat_0 
  lon2d = grb_fil->gridlon_0

  ;print(" " + varc@units)
  ;Interpolate using inverse distance squared weighting method
  varr = rcm2rgrid_Wrap(lat2d,lon2d,varc,lat,lon,0)
  ;
  if (opt) then
    itmp = area_conserve_remap_Wrap (lon,lat, varr , lonx, latx, False) 
    delete(varr)
    varr = itmp
  end if

  return(varr)
end

begin
 diri=getenv("DATA")
 scratch=getenv("SCRATCH")

 yyyy = 2020
 mm   = 10 

 days   = days_in_month(yyyy,mm)
 dly    = 3.0/33. ;;; 3.0  ;km
 dlx    = 3.0/33. ;;; 3.0  ;km

 ;;;;;
 do dd = 1,days

 ddiri=diri + "/aws/hrrr/" + yyyy + "" + sprinti("%0.2i",mm) + "" + sprinti("%0.2i",dd)

 fils = systemfunc("ls " + ddiri + "/*.grib2")
 nfils = dimsizes(fils)

 dgeo3km = getdgeo(3.,3.)
 ;3 km setup
 dgeo = getdgeo(3.,3.)
 dlat = dgeo(0)
 dlon = dgeo(1)
 nlat = tointeger((39.2-38.8)/dlat) + 2
 nlon = tointeger((108.39-107.525)/dlon) + 2 
 lat  = ispan(0,nlat-1,1)*dlat + 38.8
 lon  = ispan(0,nlon-1,1)*dlon -108.39

 edgelat = (/lat(8)-dlat/2., lat(9)+dlat/2./)
 edgelon = (/lon(6)-dlon/2., lon(10)+dlon/2./)
 print(edgelat + " " + edgelon)

 ;Cropping extent from 3 km domain
 ;39.0007 39.0258 251.818 251.957
 ;dlat= 0.0250893 dlon= 0.0346832

 ;90 m setup
 dgeo = getdgeo(dlx,dly)
 dlatx = dgeo(0)
 dlonx = dgeo(1)
 swlat = lat(8) - dlat/2. + dlatx/2.
 swlon = (lon(6) - dlon/2. + dlonx/2.) 
 nlat = 2*33 ;65  (8-9) 
 nlon = 5*33 ;165 (6-10)
 latx  = ispan(0,nlat-1,1)*dlatx + swlat 
 lonx  = ispan(0,nlon-1,1)*dlonx + swlon

 edgelat = (/min(latx)-dlatx/2., max(latx)+dlatx/2./)
 edgelon = (/min(lonx)-dlonx/2., max(lonx)+dlonx/2./)
 print(edgelat + " " + edgelon)

 print(swlat + " " + swlon)
 print(nlat + " " + nlon)
 print(dlatx + " " + dlonx)

 do nf = 0, nfils-1
   print("Processing " + fils(nf))
   tstr = str_split(fils(nf),"/")
   tstr2 = str_split(tstr(6),".")
   tstr3 = stringtochar(tstr2(1))
   tstr4 = chartostring(tstr3(1:2))

   grb_fil = addfile(fils(nf),"r")
   ;varnames = getfilevarnames(grb_fil)
   fopt = True ;downscale
   t2m = extract_interp(grb_fil,lat,lon,"TMP_P0_L103_GLC0",latx,lonx,fopt)    ;K
   sp  = extract_interp(grb_fil,lat,lon,"SPFH_P0_L103_GLC0",latx,lonx,fopt)   ;kg/kg
   uu  = extract_interp(grb_fil,lat,lon,"UGRD_P0_L103_GLC0",latx,lonx,fopt)   ;m/s
   vv  = extract_interp(grb_fil,lat,lon,"VGRD_P0_L103_GLC0",latx,lonx,fopt)   ;m/s
   pp  = extract_interp(grb_fil,lat,lon,"PRES_P0_L1_GLC0",latx,lonx,fopt)     ;Pa
   pr  = extract_interp(grb_fil,lat,lon,"PRATE_P0_L1_GLC0",latx,lonx,fopt)   ;kg/m2/s
   dl  = extract_interp(grb_fil,lat,lon,"DLWRF_P0_L1_GLC0",latx,lonx,fopt)    ;W/m2
   dx  = extract_interp(grb_fil,lat,lon,"DSWRF_P0_L1_GLC0",latx,lonx,fopt)    ;W/m2

   oupfile = scratch + "/interp_" + tstr(5) + "" + tstr4 + ".nc"
   system("rm " + oupfile)
   fout = addfile(oupfile,"c")
   fout->t2m = t2m
   fout->sp  = sp
   fout->uu  = uu
   fout->vv  = vv
   fout->pp  = pp
   fout->pr  = pr
   fout->dl  = dl
   fout->dx  = dx
   if (fopt) then
    fout->lat = latx
    fout->lon = lonx
   else
    fout->lat = lat
    fout->lon = lon
   end if
 end do
 delete(fils)
 end do

 ;Files
  outdir = scratch + "/interp_" + yyyy + "" + sprinti("%0.2i",mm) 
  system("rm -rf " + outdir)
  system("mkdir " + outdir)
  system("mv " + scratch + "/*.nc " + outdir)
end
