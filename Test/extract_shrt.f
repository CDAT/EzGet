c f77 -e extract_shrt.f -L$PCMDI/lib -lezget -lcdms -ldrs -L/usr/local/netcdf-2.3.2/libsrc -lnetcdf -o extract_shrt
c f77 -e -g extract_shrt.f -L/zooks0/ktaylor/pcmdi/util/ezget -lezgetdebug -L$PCMDI/lib -lcdms -ldrs -L/usr/local/netcdf-2.3.2/libsrc -lnetcdf -o extract_shrt
c  IBM
c xlf extract_shrt.f  -L/u/williams/ezget -lezget -L/u/williams/devel/cdms/lib -lcdms -L/u/williams/drs/lib -ldrs -L/u/williams/netcdf/libsrc -lnetcdf -o extract_shrt
c   HP
c fort77  +U77 extract_shrt.f -L/afs/nersc.gov/u/ahs/u4160/ezget -lezget -L/afs/nersc.gov/software/pub/drs/lib -lcdms -ldrs -lm -L/afs/nersc.gov/projects/graphics/netcdf/2.3.2pl4/lib -lnetcdf -o extract_shrt
c   sun/SOLARIS
c f77 -e extract_shrt.f -L/pcmdi/cirrus/ktaylor/ezget/solaris -lezget -L/pcmdi/cirrus/drs_solaris -lcdms -ldrs -L/pcmdi/cirrus/netcdf-2.3.2/libsrc -lnetcdf -o extract_shrt
c   SGI
c f77 extract_shrt.f -L/pcmdi/ktaylor/pcmdi/util/ezget/sgi -lezget -L/usr/local/cdms/lib -lcdms -L/usr/local/drs/drslib -ldrs -L/usr/local/lib -lnetcdf -o extract_shrt

      program extract

c   This program 
c     ** retrieves the first 12 months of global data for a variable 
c        (named 'tas' in this example) that is stored in a file (named 
c        '/scratch/staff/lisa/amip_obs/nasa-amip_t').
c     ** creates an array of "weights" with elements proportional to
c        grid cell area (except for grid cells with missing data
c        where the "weight" is set to 0.0).
c     ** obtains the longitude and latitude coordinates and the 
c        length of each dimension retrieved.
c     ** prints out a portion of the retrieved data.


      parameter (nlon=100, nlat=50, nmon=12, n4=0)

      real adata(nlon,nlat,nmon), wtsmask(nlon,nlat,nmon),
     &      alon(nlon), alat(nlat)
      integer lons, lats, mons, i4

c     -------------------------------------------
c     Initialize EzGet:

      call initget

c     -------------------------------------------
c     Define "missing" data value. 

      call defmisc('input missing value', 'real', 1.0e20) 

c     -------------------------------------------
c     Define variable 1 as 'tas' and indicate path/filename where 
c     it is stored.  In subsequent calls to EzGet (e.g., defdim
c     getdata, and getcoord), this field will be referred to by
c     the index assigned here (i.e., 1).

      call defvar(1, 'tas', '/scratch/staff/lisa/amip_obs/nasa-amip_t') 

c     -------------------------------------------
c     Define domain for variable 1 and define order that EzGet will 
c     retrieve data.  

      call defdim(1, 1, 'longitude', 'width', 'nearest',
     &                                             -180.0, 180.0, 360.0)
      call defdim(1, 2, 'latitude', 'cosine', 'range', 
     &                                               90.0, -90.0, 0.0)
      call defdim(1, 3, 'time',    'unit', 'nearest', 1.0, 12.0, 0.0)

c     -------------------------------------------
c     Extract variable 1 from file and create missing data mask.
c     Also return actual dimensions of retrieved data.

      lons = 0
      lats = 0
      mons = 12
      i4 = 0
      call getdata(1,  nlon,nlat,nmon,n4,  lons,lats,mons,i4,
     &             wtsmask, adata)

c     -------------------------------------------
c     Retrieve longitude and latitude coordinates of variable 1 (in 
c     the same order as the retrieved data).

      call getcoord(1, 1, alon)
      call getcoord(1, 2, alat)

c     -------------------------------------------
c     Close all files opened by EzGet.

      call closeget

c     -------------------------------------------
c     Write data retrieved (last month only):

      write(*,'("longitudes = " / (10f8.3))') (alon(i), i=1,lons)

      write(*,'(/ "month = ", i3)') mons

      do 100 j=1,lats
        write(*,'(/"latitude = ", f8.3)') alat(j)
        write(*,'(10f8.3)') (adata(i,j,mons), i=1,lons)
  100 continue

      end
