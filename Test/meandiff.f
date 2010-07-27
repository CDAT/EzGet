c f77 -e meandiff.f -L$PCMDI/lib -lezget -lcdms -ldrs -L/usr/local/netcdf-2.3.2/libsrc -lnetcdf -o meandiff
c f77 -e -g meandiff.f -L/zooks0/ktaylor/pcmdi/util/ezget -lezgetdebug -L$PCMDI/lib -lcdms -ldrs -L/usr/local/netcdf-2.3.2/libsrc -lnetcdf -o meandiff
c  IBM
c xlf meandiff.f  -L/u/williams/ezget -lezget -L/u/williams/devel/cdms/lib -lcdms -L/u/williams/drs/lib -ldrs -L/u/williams/netcdf/libsrc -lnetcdf -o meandiff
c   HP
c fort77  +U77 meandiff.f -L/afs/nersc.gov/u/ahs/u4160/ezget -lezget -L/afs/nersc.gov/software/pub/drs/lib -lcdms -ldrs -lm -L/afs/nersc.gov/projects/graphics/netcdf/2.3.2pl4/lib -lnetcdf -o meandiff
c   sun/SOLARIS
c f77 -e meandiff.f -L/pcmdi/cirrus/ktaylor/ezget/solaris -lezget -L/pcmdi/cirrus/drs_solaris -lcdms -ldrs -L/pcmdi/cirrus/netcdf-2.3.2/libsrc -lnetcdf -o meandiff
c   SGI
c f77 meandiff.f -L/pcmdi/ktaylor/pcmdi/util/ezget/sgi -lezget -L/usr/local/cdms/lib -lcdms -L/usr/local/drs/drslib -ldrs -L/usr/local/lib -lnetcdf -o meandiff


      program meandiff
      
c   This program
c
c    ** retrieves observed and model-simulated monthly surface air 
c       temperature for the year 1979 (first year of AMIP run).
c    ** determines which grid cells are not missing any data and then
c       for these grid cells
c    ** computes the global average, annual mean difference between the
c       observed and model-simulated temperatures.
c
c       The data contributing to the averages are area-weighted and also
c       weighted by the number of days in each month.

      parameter (nlon=100, nlat=50, nmon=12, n4=0)

      real datamodl(nlon,nlat,nmon), wtsmodl(nlon,nlat,nmon),
     &     dataobs(nlon,nlat,nmon), wtsobs(nlon,nlat,nmon),
     &     wtobs(nlon,nlat)
      double precision asum, wtsum
      integer lons, lats, mons, i4, i, j, m, mm

c     -------------------------------------------
c     Initialize EzGet and define "missing" data value:

      call initget
      call defmisc('input missing value', 'real', 1.0e20) 

c     -------------------------------------------
c     Define variable 1 as the observed surface air temperature and 
c     indicate path/filename where data are stored.

      call defvar(1, 'tas', '/scratch/staff/lisa/amip_obs/nasa-amip_t')

c     -------------------------------------------
c     Define domain for variable 1.

      call defdim(1, 1, 'longitude', 'width', 'range',
     &                                             -180.0, 180.0, 360.0)
      call defdim(1, 2, 'latitude', 'cosine', 'range', -90.0, 90.0, 0.0)
      call defdim(1, 3, 'time',    'month', 'nearest', 109., 120., 0.0)

c     -------------------------------------------
c     Extract variable 1 (observed field). 

      lons = 0
      lats = 0
      mons = 0
      i4 = 0
      call getdata(1,  nlon,nlat,nmon,n4,  lons,lats,mons,i4,
     &             wtsobs, dataobs)

c     -------------------------------------------
c     Define variable 2 as the modeled surface air temperature and 
c     indicate path/filename where data are stored.

      call defvar(2, 'tas', '/amipsp/drs/tas/tas_bmr')

c     -------------------------------------------
c     Define domain for variable 2.

      call defdim(2, 1, 'longitude', 'width', 'range',
     &                                             0.0, 0.0, 360.0)
      call defdim(2, 2, 'latitude', 'gaussian', 'range',
     &                                              0.0, 0.0, 0.0)
      call defdim(2, 3, 'time',    'month', 'nearest', 109., 120., 0.0)

c     -------------------------------------------
c      Regrid model output to observed grid

      call defregrd(2, 'to', 1, 'area-weighted', 0,0.0,0.0, 0,0.0,0.0)

c     -------------------------------------------
c     Extract variable 2 (modeled field). 

      call getdata(2,  nlon,nlat,nmon,n4,  lons,lats,mons,i4,
     &             wtsmodl, datamodl)


c       -----------------------------------------
c       Compute annual mean at each grid cell, but only if all 12 months
c       of data are available for both the model and the observations.

c       If data is available, weight by the area of the grid cell and
c       the number of days in a month. 

        do 300 j=1,lats
          do 200 i=1,lons

            mm = 0
            do 100 m=1,mons
              if (wtsmodl(i,j,m)*wtsobs(i,j,m) .gt. 0.0) mm = mm + 1
  100       continue

            if (mm .eq. 12) then
              wtobs(i,j) = wtsobs(i,j,1)
            else
              wtobs(i,j) = 0.0
            endif

  200     continue
  300   continue

c     -------------------------------------------
c     Compute global average, annual mean difference

      asum = 0.0
      wtsum = 0.0

      do 600 m=1,nmon
        do 500 j=1,nlat
          do 400 i=1,nlon

            asum = asum + wtobs(i,j)*(datamodl(i,j,m)-dataobs(i,j,m))
            wtsum = wtsum + wtobs(i,j)

  400     continue
  500   continue
  600 continue

c     -----------------------------------------
c     Write area-weighted annual mean difference between model and 
c     observed fields.

      if (wtsum .gt. 0.0) then 

        asum = asum/wtsum

        write(*,'("annually-averaged fraction of globe with data: ",
     &            f12.7)') wtsum
        write(*,'("global, annual mean difference: ", f12.3)') asum

      else

        write(*,'( 
     &  " data missing everywhere for at least 1 month of the year")')

      endif

c     -------------------------------------------
c     Close all files opened by EzGet.

      call closeget

      end
