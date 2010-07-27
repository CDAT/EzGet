c f77 -e getregn.f -L$PCMDI/lib -lezget -lcdms -ldrs -L/usr/local/netcdf-2.3.2/libsrc -lnetcdf -o getregn
c f77 -e -g getregn.f -L/zooks0/ktaylor/pcmdi/util/ezget -lezgetdebug -L$PCMDI/lib -lcdms -ldrs -L/usr/local/netcdf-2.3.2/libsrc -lnetcdf -o getregn
c  IBM
c xlf getregn.f  -L/u/williams/ezget -lezget -L/u/williams/devel/cdms/lib -lcdms -L/u/williams/drs/lib -ldrs -L/u/williams/netcdf/libsrc -lnetcdf -o getregn
c   HP
c fort77  +U77 getregn.f -L/afs/nersc.gov/u/ahs/u4160/ezget -lezget -L/afs/nersc.gov/software/pub/drs/lib -lcdms -ldrs -lm -L/afs/nersc.gov/projects/graphics/netcdf/2.3.2pl4/lib -lnetcdf -o getregn
c   sun/SOLARIS
c f77 -e getregn.f -L/pcmdi/cirrus/ktaylor/ezget/solaris -lezget -L/pcmdi/cirrus/drs_solaris -lcdms -ldrs -L/pcmdi/cirrus/netcdf-2.3.2/libsrc -lnetcdf -o getregn
c   SGI
c f77 getregn.f -L/pcmdi/ktaylor/pcmdi/util/ezget/sgi -lezget -L/usr/local/cdms/lib -lcdms -L/usr/local/drs/drslib -ldrs -L/usr/local/lib -lnetcdf -o getregn

       program getregn
      
c   This program
c    ** retrieves surface temperature data for the region of
c       North America north of 23.5 N latitude.  (Note: data are 
c       extracted only for grid cells whose coordinates are not outside 
c       the range 23.5 N to 90.0 N and 190 W to 40 W.)
c    ** creates an array of "weights", with elements set proportional
c       to the area of the grid cells falling within the selected 
c       geographical region (except for those grid cells that contain 
c       "missing" data or which lie outside North America, in which case
c       the weight is set to 0.0). 
c    ** obtains the longitude and latitude coordinates and
c       the length of each dimension retrieved.

c   This program is in many ways similar to example 1, and 
c   further explanation can be found in the comments appearing in 
c   that example (program extract).

c   Note that the region is selected in two ways.  
c   First, all data for grid cells in the region 23.5 N to 90.0 N and
c   190 W and 40 W are extracted, and then all grid-cells that are
c   outside the North American boundaries are masked out.  The size of 
c   the arrays, adata and wtsmask, only need to be large enough to
c   accommodate the region extracted (before the geography mask is 
c   applied).

      parameter (nlon=35, nlat=25, nmon=15, n4=0)

      real adata(nlon,nlat,nmon), wtsmask(nlon,nlat,nmon),
     &      alon(nlon), alat(nlat)
      integer lons, lats, mons, i4
c     -------------------------------------------
c     Initialize EzGet and define "missing" data value:

      call initget
      call defmisc('input missing value', 'real', 1.0e20) 

c     -------------------------------------------
c     Define variable 1 as 'tas' and define its domain and the order 
c     that EzGet will retrieve data. 

      call defvar(1, 'tas', '/scratch/staff/lisa/amip_obs/nasa-amip_t') 

      call defdim(1, 1, 'longitude', 'width', 'range',
     &                                             -190.0, -40.0, 360.0)

c     >>> Note that the longitude range specified above is large enough 
c     >>> to accommodate all North American grid cells.

      call defdim(1, 2, 'latitude', 'cosine', 'range', 23.5, 90.0, 0.0)
      call defdim(1, 3, 'time',    'unit', 'nearest', 1.0, 12.0, 0.0)

c     -------------------------------------------
c    >>>  Define variable 2 as 'sftbyrgn' and indicate path/filename 
c    >>>  where it is stored.  This variable should contain a "geography
c    >>>  mask" that is compatible with the EzGet convention for 
c    >>>  identifying different geographical regions (see documentation 
c    >>>  of subroutine defgeog), and it should be on the same grid as 
c    >>>  variable 1 defined above.  In subsequent calls to EzGet (e.g., 
c    >>>  defgeog) this field will be referred to by the index assigned
c    >>>  here (i.e., 2).

      call defvar(2, 'sftbyrgn', '/amipsp/drs/sftbyrgn/sftbyrgn_gla') 

c     -------------------------------------------
c    >>>  Control geographical masking of retrieved data.

      call defgeog(1, 'in', 2, 'North America')

c    >>>  where the arguments indicate the following:
c
c    >>>   1:    The first argument has the value 1 and indicates that 
c    >>>     the geography mask will be applied to variable 1. 
c    >>>   'in': indicates that the masking should be done before any
c    >>>     regridding is performed. (See later examples for further
c    >>>     explanation.)
c    >>>   2:    Indicates that variable 2 contains the geography data.
c    >>>   'North America': indicates that this is the region of 
c    >>>         interest and any data outside this region should be
c    >>>         masked.
 
c     -------------------------------------------
c     Extract variable 1 from file and mask data outside region of 
c     interest.

      lons = 0
      lats = 0
      mons = 12
      i4 = 0
      call getdata(1,  nlon,nlat,nmon,n4,  lons,lats,mons,i4,
     &             wtsmask, adata)


c     -------------------------------------------
c     Retrieve longitude and latitude coordinates of variable 1.

      call getcoord(1, 1, alon)
      call getcoord(1, 2, alat)

c     -------------------------------------------
c     Close all files opened by EzGet.

      call closeget

c     -------------------------------------------
c     Write out weights and data retrieved (last month only):

      write(*,'("longitudes = " / (10f8.3))') (alon(i), i=1,lons)
    
      write(*,'(/ "month = ", i3)') mons

      do 100 j=1,lats
        write(*,'(/"latitude = ", f8.3)') alat(j)
        write(*,'("weights:")')
        write(*,'(9f9.6)') (wtsmask(i,j,mons), i=1,lons)
        write(*,'("temperature:")')
        write(*,'(9f9.3)') (adata(i,j,mons), i=1,lons)
  100 continue

      end
