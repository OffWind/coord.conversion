!!$--------!---------!---------!---------!---------!---------!---------%
!!$
!!$     File Name: UTM2GEO.F90
!!$
!!$     Author:    Rui Faria Pereira, Megajoule Inovação Lda.
!!$
!!$     URL:       www.megajoule.pt
!!$
!!$     Version:   1.2
!!$
!!$     Date:      22.August.2012
!!$
!!$--------!---------!---------!---------!---------!---------!---------%
program utm2geo

  real :: xlat,xlong,eastings,northings
  integer :: utm_zone

  write(*,*) 'Input eastings, XUTM [m, WGS84 Datum]: '
  read(*,*) eastings

  write(*,*) 'Input northings, XUTM [m, WGS84 Datum]: '
  read(*,*) northings

  write(*,*) 'Input UTM zone: '
  read(*,*) utm_zone

!!$ Coordinate transformation LAT/LON -> UTM WGS84
  call utm2geo_wgs84(eastings,northings,xlon,xlat,utm_zone)
  write(*,*) 'XUTM/YUTM [m]: ',eastings,northings
  write(*,*) 'LAT/LON      : ',xlat,xlon

end program utm2geo

!!$###############################
!!$###############################
      subroutine utm2geo_wgs84(xutm,yutm,lon,lat,zoneutm) 

        implicit none

        real :: a,b,k0,e,ee,arc,mu,e1,J1,J2,J3,J4,fp,C1,T1,N1,R1,D,X
        real :: Q1,Q2,Q3,Q4,Q5,Q6,Q7
        real :: xutm,yutm,lon,lat,xutmm,yutmm,pi
        integer :: zoneutm,centmerid

        pi = 4.0 * atan(1.0)

        if (zoneutm.lt.0) yutm = 10000000.0 - yutm

        xutmm = xutm
        yutmm = yutm

        X=500000.0-xutmm

        centmerid = 6*abs(zoneutm)-183

!!$ Constantes do Datum, no caso WGS84, 
        a=6378137.0            !!$ Raio Equatorial em metros
        b=6356752.3142         !!% Raio Polar em metros
        k0=0.9996              !!% Factor de escala

!!$ Excentricidade da terra numa secÁ„o eliptica
        e=(1.0-b**2.0/a**2.0)**0.5      
        ee=(e*a/b)**2.0

!!$ Calculate the Meridional Arc
        arc=yutmm/k0   

!!$ Calculate Footprint Latitude
        mu=arc/(a*(1.0-e**2.0/4.0-3.0*e**4.0/64.0-5.0*e**6.0/256.0))
        e1=(1.0-(1.0-e**2.0)**(1.0/2.0))/(1.0+(1.0-e**2.0)**(1.0/2.0))
        J1=3.0*e1/2.0-27.0*e1**3.0/32.0
        J2=21.0*e1**2.0/16.0-55*e1**4.0/32.0
        J3=151.0*e1**3.0/96.0
        J4=1097.0*e1**4.0/512.0

        fp=mu+J1*sin(2.0*mu)+J2*sin(4.0*mu)+J3*sin(6.0*mu)+J4*sin(8.0*mu)

!!$ Constants for formulas
        C1=ee*cos(fp)**2.0
        T1=tan(fp)**2.0
        N1=a/(1.0-(e*sin(fp))**2.0)**(1.0/2.0)
        R1=a*(1.0-e**2.0)/(1-(e*sin(fp))**2.0)**(3.0/2.0)
        D=X/(N1*k0)

!!$ Coefficients for Calculating Latitude
        Q1=N1*tan(fp)/R1
        Q2=D**2.0/2.0
        Q3=(5.0+3.0*T1+10.0*C1-4.0*C1**2.0-9.0*ee)*D**4.0/24.0
        Q4=(61.0+90.0*T1+298.0*C1+45.0*T1**2.0-252.0*ee-3.0*C1**2.0)*D**6.0/720.0

!!$latitude final em GEO
        lat=(fp-Q1*(Q2-Q3+Q4))*180.0/pi
        if (zoneutm.lt.0) lat = -lat 

!!$ Coefficients for Calculating Longitude
        Q5=D
        Q6=(1.0+2.0*T1+C1)*D**3.0/6.0
        Q7=(5.0-2.0*C1+28.0*T1-3.0*C1**2.0+8.0*ee+24.0*T1**2.0)*D**5.0/120.0

!!$longitude final em GEO
        lon=centmerid-(Q5-Q6+Q7)/cos(fp)*180.0/pi 

      end subroutine utm2geo_wgs84
