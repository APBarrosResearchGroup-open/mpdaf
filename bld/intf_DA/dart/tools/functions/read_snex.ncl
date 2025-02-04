function read_snex(csvfil,icol)

begin

  var1d = asciiread(csvfil,-1,"string")
  str2d = str_split_csv(var1d, ",", 0)
  ndim  = dimsizes(str2d)
  var2d = tofloat(str2d(1:ndim(0)-1,1:ndim(1)-1))
  tstr  = str2d(1:ndim(0)-1,0)
  hdr   = str2d(0,1:ndim(1)-1)

  ;2022-08-03 19:00:00
  units  = "hours since 1900-01-01 00:00:00"
  time   = new(ndim(0)-1,"double")
  doy    = new(ndim(0)-1,"double")
  cal    = new((/6,ndim(0)-1/),"integer")

  do it  = 0, ndim(0)-2
    ymdhms = tstr(it)
    tmp = str_split(ymdhms," ")
    ymd = tointeger(str_split(tmp(0),"-"))
    hms = tointeger(str_split(tmp(1),":"))

    time(it)  = cd_inv_calendar(ymd(0),ymd(1),ymd(2),hms(0),hms(1),hms(2),units, 0)
    doy(it)   = day_of_year(ymd(0),ymd(1),ymd(2))*1. + (hms(0) + hms(1)/60. + hms(2)/3600.)/24.
    cal(:,it) = (/ymd(0),ymd(1),ymd(2),hms(0),hms(1),hms(2)/)
  end do

  ;print(" " + hdr)
  print("Returning ..." + hdr(icol))
  var    = (/var2d(:,icol)/)
  var@long_name = hdr(icol)
  var!0  = "time"
  var&time = time
  var@_FillValue = -9999.
  var@doy  = doy
  var@cal  = cal
  return(var)
end
