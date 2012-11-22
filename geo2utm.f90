!!$--------!---------!---------!---------!---------!---------!---------%
!!$
!!$     File Name: GEO2UTM.F90
!!$
!!$     Author:    Rui Faria Pereira, Megajoule Inovação Lda.
!!$
!!$     URL:       www.megajoule.pt
!!$
!!$     Version:   1.1
!!$
!!$     Date:      11.August.2012
!!$
!!$--------!---------!---------!---------!---------!---------!---------%
  program geo2utm

  real :: xlat,xlong,eastings,northings
  integer :: utm_zone
  character (len=1) :: utm_zonestr


  write(*,*) 'Input latitude [decimal degrees, WGS84 Datum]: '
  read(*,*) xlat

  write(*,*) 'Input longitude [decimal degrees, WGS84 Datum]: '
  read(*,*) xlong

  write(*,*) 'Input UTM zone: '
  read(*,*) utm_zone

!!$ Coordinate transformation LAT/LON -> UTM WGS84
  call geo2utm_wgs84(xlat,xlong,eastings,northings,utm_zone_aux,utm_zonestr)
  write(*,*) 'LAT/LON      : ',xlat,xlon
  write(*,*) 'XUTM/YUTM [m]: ',eastings,northings

  end program

!!$###############################
!!$###############################
subroutine geo2utm_wgs84(lat,long,eastings,northings,utm_zone,utm_zonestr) 

  implicit none

  real :: lat,long,eastings,northings,latt,longg
  real :: a,b,k0,e,ee,n,rho,nu,p,centmerid,A0,B0,C0,D0,E0,S
  real :: Ki,Kii,Kiii,Kiv,Kv,pi
  integer :: utm_zone
  character(len=1) :: utm_zonestr

  pi = 4.*atan(1.0)

  longg=long*pi/180
  latt=lat*pi/180

  !! Cálculo da zona UTM
  if (long.lt.0) then
     utm_zone=int(((180+long)/6)+1) ! Arredondar por defeito
  else
     utm_zone=int(long/6)+31        ! Arredondar por defeito
  endif

  !! Cálculo do Meridiano central
  centmerid=6.*utm_zone-183.

  !! Constantes do Datum, no caso WGS84. 
  !! Para outros datums basta mudar estes dois valores
  a = 6378137         !Raio Equatorial em metros
  b = 6356752.3142    !Raio Polar em metros

  k0 = 0.9996         !Factor de escala
  e = (1-b**2./a**2.)**0.5 !Excentricidade da terra numa secção eliptica
  ee = (e*a/b)**2. 
  n = (a-b)/(a+b) 
  rho = a*(1.-e**2.)/(1.-e**2.*sin(latt)**2.)**(3./2.)
  nu = a/((1.-(e*sin(latt))**2.)**(1./2.))
  p=(longg-centmerid*pi/180.)

  !! Cálculo do comprimento do arco do meridiano S
  A0=a*(1.-n+(5.*n*n/4.)*(1.-n)+(81.*n**4./64.)*(1.-n))
  B0=(3.*a*n/2.)*(1.-n-(7.*n*n/8.)*(1.-n)+55.*n**4./64.)
  C0=(15.*a*n*n/16.)*(1.-n+(3.*n*n/4.)*(1.-n))
  D0=(35.*a*n**3./48.)*(1.-n+11.*n*n/16.)
  E0=(315.*a*n**4./51.)*(1.-n)

  S=A0*latt-B0*sin(2.*latt)+C0*sin(4.*latt)-D0*sin(6.*latt)+E0*sin(8.*latt)

  !! Cálculo dos coeficientes UTM
  Ki=S*k0
  Kii=(nu*sin(latt)*cos(latt)*k0)/2.
  Kiii=((nu*sin(latt)*cos(latt)**3.)/24.)*(5.-tan(latt)**2.+9.*ee*cos(latt)**2.+ &
       4.*ee**2.*cos(latt)**4.)*k0
  Kiv=nu*cos(latt)*k0
  Kv=(cos(latt))**3.*(nu/6.)*(1.-tan(latt)**2.+ee*cos(latt)**2.)*k0

  !! Coordenadas em UTM
  if (lat.lt.0) then
     northings=10E6+(Ki+Kii*p*p+Kiii*p**4.)
     utm_zonestr='S'
  else
     northings=(Ki+Kii*p*p+Kiii*p**4.)
     utm_zonestr='N'
  endif
  eastings=500000.+(Kiv*p+Kv*p**3.)

end subroutine geo2utm_wgs84
