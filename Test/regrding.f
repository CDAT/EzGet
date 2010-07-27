c f77 -e regrding.f -L$PCMDI/lib -lezget -lcdms -ldrs -L/usr/local/netcdf-2.3.2/libsrc -lnetcdf -o regrding
c f77 -e -g regrding.f -L/zooks0/ktaylor/pcmdi/util/ezget -lezgetdebug -L$PCMDI/lib -lcdms -ldrs -L/usr/local/netcdf-2.3.2/libsrc -lnetcdf -o regrding
c  IBM
c xlf regrding.f  -L/u/williams/ezget -lezget -L/u/williams/devel/cdms/lib -lcdms -L/u/williams/drs/lib -ldrs -L/u/williams/netcdf/libsrc -lnetcdf -o regrding
c   HP
c fort77  +U77 regrding.f -L/afs/nersc.gov/u/ahs/u4160/ezget -lezget -L/afs/nersc.gov/software/pub/drs/lib -lcdms -ldrs -lm -L/afs/nersc.gov/projects/graphics/netcdf/2.3.2pl4/lib -lnetcdf -o regrding
c   sun/SOLARIS
c f77 -e regrding.f -L/pcmdi/cirrus/ktaylor/ezget/solaris -lezget -L/pcmdi/cirrus/drs_solaris -lcdms -ldrs -L/pcmdi/cirrus/netcdf-2.3.2/libsrc -lnetcdf -o regrding
c   SGI
c f77 regrding.f -L/pcmdi/ktaylor/pcmdi/util/ezget/sgi -lezget -L/usr/local/cdms/lib -lcdms -L/usr/local/drs/drslib -ldrs -L/usr/local/lib -lnetcdf -o regrding

       program regrding
      
c   This program 
c     ** retrieves global surface temperature data and maps
c        it to a 10 degree by 15 degree latitude-longitude "target"  
c        grid using an area-weighted averaging scheme.  
c     ** creates an array of "weights", with elements set
c        proportional to the sum of the areas of the original (i.e.,
c        source) grid cells that contribute to each target cell. 
c     ** obtains the longitude and latitude coordinates
c        of the target grid and the length of each dimension of the 
c        target grid.

c   This program is in many ways similar to example 1, and 
c   further explanation can be found in the comments appearing in 
c   that example (program extract).

c   Note that in order to apply the area-weighting regridding algorithm,
c   the user must specify what type of grid the original data were
c   stored on (e.g., gaussian, evenly-spaced, etc.)

c   The size of the arrays, adata and wtsmask, only need to be large
c   enough to accommodate the region extracted (i.e. after regridding), 
c   so nlon .ge. 360/15 = 24 and nlat .ge. 180/10 = 18. 

      parameter (nlon=24, nlat=18, nmon=12, n4=0)

      real adata(nlon,nlat,nmon), wtsmask(nlon,nlat,nmon),
     &      alon(nlon), alat(nlat)
      integer lons, lats, mons, i4
c     -------------------------------------------
c     Initialize EzGet and define "missing" data value:

      call initget
      call defmisc('input missing value', 'real', 1.0e20) 

c     -------------------------------------------
c     Define variable 1 as 'tas' and indicate path/filename where 
c     it is stored.

      call defvar(1, 'tas', 
     &             '/afs/nersc.gov/u/ahs/u4160/ezget/nasa-amip_t') 
c     &             '/u/williams/ezget/nasa-amip_t') 
c     &             '/pcmdi/ktaylor/pcmdi/util/ezget/nasa-amip_t') 

c     -------------------------------------------
c     Define domain for variable 1 and define order that EzGet will 
c     retrieve data.   Also indicate what type of grid the source
c     data were stored on. 

      call defdim(1, 1, 'longitude', 'width', 'range', 0.0, 0.0, 360.0)
      call defdim(1, 2, 'latitude', 'cosine', 'range', 0.0, 0.0,   0.0)
      call defdim(1, 3, 'time',    'unit', 'nearest',  1.0, 12.0, 0.0)
c
c      When regridding data, the "cycle" for longitude should always 
c      be specified (set to 360.0 in the above example).
c      The domain limits specified for longitude and latitude will be 
c      overridden by the arguments in the call to subroutine defregrd
c      below.

c     -------------------------------------------
c     Define target grid to which data should be mapped. 

      call defregrd(1, 'uniform', 0, 'area-weighted', 18, 85.0, -10.0,
     &                                               24, -172.5, 15.0)

c     where the arguments in the subroutine call indicate the following:
c
c      1:    The first argument has the value 1 and indicates that the  
c            the regridding will be applied to variable 1. 
c      'uniform': indicates that the target grid will be rectangular 
c            grid of evenly-spaced latitude and longitude cells.
c      0:    This argument is ignored because the 2nd argument was set 
c            to 'uniform'.   
c      'area-weighted': indicates that an area-weighted mapping scheme
c            should be used.
c      18, 85.0, -10.0:
c            These 3 arguments define the latitude grid that will be 
c            created.  In this case a grid with 18 latitude cells 
c            is created with the first grid cell centered at 85.0 
c            degrees north and proceeding southward in increments of
c            10 degrees.
c      24, -172.5, 15.0:
c            These 3 arguments define the longitude grid, indicating 
c            that there will be 24 longitude grid cells, with the first
c            grid cell centered at -172.5 degrees west and proceeding 
c            eastward in increments of 15 degrees.
c
c     -------------------------------------------
c     Extract variable 1 from file and map to target grid defined above.

      lons = 24
      lats = 18
      mons = 12
      i4 = 0
      call getdata(1,  nlon,nlat,nmon,n4,  lons,lats,mons,i4,
     &             wtsmask, adata)

c     Note: we have defined the expected longitude and latitude
c     dimension to be 24 and 18, respectively (as specified by lons and 
c     lats) because we know the data will be regridded to the grid 
c     defined by the call to defregrd, which specifies a global 10 x 15 
c     degree latitude-longitude grid.  It is always good practice to 
c     define the expected dimensions if they are known.

c     -------------------------------------------
c     Retrieve longitude and latitude coordinates of variable 1.

      call getcoord(1, 1, alon)
      call getcoord(1, 2, alat)

c     Note:  According to the parameters defining the target grid,
c          alon should contain -172.5, -157.5, . . . 172.5 and
c          alat should contain  85.0, 75.0. . . -85.0

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

