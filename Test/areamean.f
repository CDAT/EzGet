c f77 -e areamean.f -L$PCMDI/lib -lezget -lcdms -ldrs -L/usr/local/netcdf-2.3.2/libsrc -lnetcdf -o areamean
c
c Macintosh OS X
c gfortran -g -save-temps -e areamean.f -arch x86_64 -L$HOME/drs/ezget -L$HOME/cdat/libcdms/lib -L$HOME/drs/libdrs/lib -lezget -lcdms -ldrs -lnetcdf -o areamean
c f77 -e -g areamean.f -L/zooks0/ktaylor/pcmdi/util/ezget -lezgetdebug -L$PCMDI/lib -lcdms -ldrs -L/usr/local/netcdf-2.3.2/libsrc -lnetcdf -o areamean
c  IBM
c xlf areamean.f  -L/u/williams/ezget -lezget -L/u/williams/devel/cdms/lib -lcdms -L/u/williams/drs/lib -ldrs -L/u/williams/netcdf/libsrc -lnetcdf -o areamean
c   HP
c fort77  +U77 areamean.f -L/afs/nersc.gov/u/ahs/u4160/ezget -lezget -L/afs/nersc.gov/software/pub/drs/lib -lcdms -ldrs -lm -L/afs/nersc.gov/projects/graphics/netcdf/2.3.2pl4/lib -lnetcdf -o areamean
c   sun/SOLARIS
c f77 -e areamean.f -L/pcmdi/cirrus/ktaylor/ezget/solaris -lezget -L/pcmdi/cirrus/drs_solaris -lcdms -ldrs -L/pcmdi/cirrus/netcdf-2.3.2/libsrc -lnetcdf -o areamean
c   SGI
c f77 areamean.f -L/pcmdi/ktaylor/pcmdi/util/ezget/sgi -lezget -L/usr/local/cdms/lib -lcdms -L/usr/local/drs/drslib -ldrs -L/usr/local/lib -lnetcdf -o areamean

       program areamean
      
c   This program
c
c    ** retrieves surface temperature data for the region of
c       North America north of 23.5 N latitude.  It then uses the 
c       regridding capability of EzGet to compute the area average of a 
c       single "grid-cell" covering the region of North America north of
c       23.5 N.  Up to 120 months of data can be extracted at once and 
c       the area means are returned to this program as a vector, one
c       element for each month of data extracted.
c   
c    ** creates a vector of "weights" (one for each month)
c       proportional to the area of the region over which the means have
c       been computed.  The elements of this vector should be 
c       identical (except possibly if grid cells are missing data). 

c   This program is in many ways similar to examples 1, 2 and 3, and 
c   further explanation can be found in the comments appearing in 
c   those examples (programs extract, getregn, and regrding).

c  >>> Note that nlon and nlat can be declared as small as 1, because 
c  >>> these dimensions reduce to a single grid cell after regridding.

      parameter (nmon=120, nlon=1, nlat=1, n4=0)

      real amean(nmon, nlon, nlat), wtsmask(nmon, nlon, nlat) 
      integer lons, lats, mons, i4
c     -------------------------------------------
c     Initialize EzGet and define "missing" data value:

      call initget
      call defmisc('input missing value', 'real', 1.0e20) 

c     -------------------------------------------
c     Define variable 1 as observed surface temperature data ('tas')
c     and define its domain:

      call defvar(1, 'tas', '/scratch/staff/lisa/amip_obs/nasa-amip_t') 

      call defdim(1, 1, 'time', 'unit', 'as saved', 0.0, 0.0, 0.0)
      call defdim(1, 2, 'longitude', 'width', 'range', 0.0, 0.0, 360.0)
      call defdim(1, 3, 'latitude', 'cosine', 'range', 0.0, 0.0,   0.0)

c     >>>    Note that in the 1st call to defdim above, the last 4 
c     >>>    arguments are:
c     >>>    'as saved', 0.0, 0.0, 0.0:
c     >>>    'as saved' specifies that all months should be retrieved in
c     >>>    the order that they were originally stored in the file.
c     >>>    (In this case the last 3 arguments [0.0, 0.0, 0.0] are 
c     >>>    ignored.) 

c     -------------------------------------------
c     Define variable 2 as the geography data needed and specify that
c     it be used to select North American data only for variable 1:

      call defvar(2, 'sftbyrgn', '/amipsp/drs/sftbyrgn/sftbyrgn_gla') 
      call defgeog(1, 'in', 2, 'North America')

c     -------------------------------------------
c     Define target grid to which data should be mapped. 

      call defregrd(1, 'uniform', 0, 'area-weighted', 
     &                        1, 56.75, 66.5,  1, -125.0, 150.0)

c     where the target grid has been specified as a single grid cell 
c         centered at 125 W longitude, 56.75 N latitude and has 
c         latitude-longitude dimensions of 66.5 x 150. 

c     -------------------------------------------
c     Extract variable 1, which, because of regridding to a single cell, 
c     contains the area-weighted mean for each month. 
 
      mons = 0
      lons = 1
      lats = 1
      i4 = 0
      call getdata(1,  nmon,nlon,nlat,n4,  mons,lons,lats,i4,
     &             wtsmask, amean)

c     -------------------------------------------
c     Close all files opened by EzGet.

      call closeget

c     -------------------------------------------
c     Write data retrieved:

      write(*,'("North American Mean (North of 23.5 N)" //
     &        " month      area      mean" / )') 
    
      write(*,'(i5, f12.7, f10.3)') 
     &                  (m, wtsmask(m,1,1), amean(m,1,1), m=1,mons)

      end
