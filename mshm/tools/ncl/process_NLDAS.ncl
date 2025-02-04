function getdgeo(dsx,dsy)
begin
 llat = get_pi("float")*39./180.
 dlat = dsy*(1./119.574)
 dlon = dsx*(1./(111.320*cos(llat)))
 ;print(dlat + " " + dlon)
 dgeo = (/dlat,dlon/)
 return(dgeo)
end

begin
 diri=getenv("DATA")
 scratch=getenv("SCRATCH")

 dly    = 3./33. ; 3.0  ;km
 dlx    = 3./33. ; 3.0  ;km

 ;;;;;

 diri=diri + "/ALBEDO/2020/"

 fils = systemfunc("ls " + diri + "/*.nc4")
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

 ;
 fopt = True
 ;Downscale

 do nf = 0 , nfils-1
   fname = fils(nf)
   tmp   = str_split(fname,"/")

   nc_in = addfile(fname,"r")
   alb   = nc_in->Albedo
   if (nf.eq.0) then
     ilat  = nc_in->lat
     ilon  = nc_in->lon
     system("mkdir " + scratch + "/interp_ALBEDO")
   end if
   if (fopt) then
     ialb  = area_conserve_remap_Wrap (ilon,ilat, alb(0,:,:) , lonx, latx, False)
   else
     ialb  = area_conserve_remap_Wrap (ilon,ilat, alb(0,:,:) , lon, lat, False)
   end if

   foutname = scratch + "/interp_ALBEDO/interp_" + tmp(5) 
   fout = addfile(foutname,"c")
   fout->Band1 = ialb
   if (fopt) then
    fout->lat = latx
    fout->lon = lonx
   else
    fout->lat = lat
    fout->lon = lon
   end if
 end do

end
