c f77 -e extract.f -L$PCMDI/lib -lezget -lcdms -ldrs -L/usr/local/netcdf-2.3.2/libsrc -lnetcdf -o extract
c f77 -e -g extract.f -L/zooks0/ktaylor/pcmdi/util/ezget -lezgetdebug -L$PCMDI/lib -lcdms -ldrs -L/usr/local/netcdf-2.3.2/libsrc -lnetcdf -o extract
c  IBM
c xlf extract.f  -L/u/williams/ezget -lezget -L/u/williams/devel/cdms/lib -lcdms -L/u/williams/drs/lib -ldrs -L/u/williams/netcdf/libsrc -lnetcdf -o extract
c   HP
c fort77  +U77 extract.f -L/afs/nersc.gov/u/ahs/u4160/ezget -lezget -L/afs/nersc.gov/software/pub/drs/lib -lcdms -ldrs -lm -L/afs/nersc.gov/projects/graphics/netcdf/2.3.2pl4/lib -lnetcdf -o extract
c   sun/SOLARIS
c f77 -e extract.f -L/pcmdi/cirrus/ktaylor/ezget/solaris -lezget -L/pcmdi/cirrus/drs_solaris -lcdms -ldrs -L/pcmdi/cirrus/netcdf-2.3.2/libsrc -lnetcdf -o extract
c   SGI
c f77 extract.f -L/pcmdi/ktaylor/pcmdi/util/ezget/sgi -lezget -L/usr/local/cdms/lib -lcdms -L/usr/local/drs/drslib -ldrs -L/usr/local/lib -lnetcdf -o extract


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
c
c   Concerning the structure of the stored data, we assume the
c   following:
c
c    1.  The field is a function of longitude, latitude, and time, but
c        the dimension order is unknown.
c    2.  The data are stored in a rectangular array.
c    3.  The size of the array retrieved is limited as follows:
c        a. The longitude dimension is no longer than 100 elements.
c        b. The latitude dimension is no longer than 50 elements.
c        c. The time dimension is no longer than 12 elements.
c    4.  The names of the dimensions (as stored in the file) are:
c        'latitude', 'longitude', and 'time' (but not necessarily in
c        that order).
c    5.  The longitude coordinates are evenly spaced.
c    6.  The latitude coordinates are evenly spaced.
c    7.  The units for longitude and latitude are degrees.
c    8.  Only one field in the file has the name, 'tas' (i.e., this 
c        variable name is unique in the file).
c    9.  The number, 1.0e20, is stored in the place of the true data
c        anywhere that data are missing (i.e., this is the "missing 
c        data" value or indicator).
c
c          
c    This program obtains the following information:
c
c    lons = the actual length of the longitude dimension of the data 
c             retrieved from storage.
c    lats = the actual length of the latitude dimension of the data
c             retrieved from storage.
c    mons = the actual number of months of data retrieved from storage.
c    alon(100) = a vector containing the longitude coordinates for the
c                 data.
c    alat(50) = a vector containing the latitude coordinates for the
c                 data.  
c    adata(100,50,12) = the retrieved array of data, which according to
c       our specifications given below, will be put in the following 
c       structure (regardless of how it was originally stored):
c       1. The dimension order will be: longitude, latitude, time 
c          (i.e., in the array, "adata", the first index is associated
c          with longitude, the second with latitude, and the third with 
c          time).
c       2. The longitudes will be ordered from west to east, starting
c          with the longitude nearest to 180 W.
c       3. The latitudes will be ordered north to south.
c       4. The months will be ordered consecutively.
c    wtsmask(100,50,12) = the created "missing data" mask which will
c       be ordered the same as "adata", with elements set proportional
c       to the grid cell area (except for grid cells with missing data.
c       where the elements will be set to 0.0)
c 
c    Note that if lons < nlon and/or lats < nlat, and/or mons < nmon,
c    then the extra elements of the array, "adata" and "wtsmask" will be
c    assigned a value of 0.0 by EzGet.  

      parameter (nlon=100, nlat=50, nmon=12, n4=0)

      real adata(nlon,nlat,nmon), wtsmask(nlon,nlat,nmon),
     &      alon(nlon), alat(nlat)
      integer lons, lats, mons, i4

c     -------------------------------------------
c     Initialize EzGet:

      call initget

c     -------------------------------------------
c     Tell EzGet to consider "missing" any data that have 
c     values (within a small tolerance) of 1.0e20. 

      call defmisc('input missing value', 'real', 1.0e20) 

c     -------------------------------------------
c     Define variable 1 as 'tas' and indicate path/filename where 
c     it is stored.  In subsequent calls to EzGet (e.g., defdim
c     getdata, and getcoord), this field will be referenced by
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

c     where in the first call to defdim above the arguments indicate
c     the following:
c
c      1, 1, 'longitude': 
c         The first 3 arguments indicate that for variable 1 (as 
c         indicated by 1st argument), the data should be retrieved such
c         that the dimension named 'longitude' (as indicated by the 
c         3rd argument) is the first dimension  (as indicated by the 2nd
c         argument).
c      'width': this argument controls the creation of the array of
c           "weights" that will be associated with the data. For each 
c           array element the weight is equal to the product of the 
c           weights defined by each dimension and here 'width' specifies 
c           that the "weights" should be proportional to the 
c           longitudinal width of each grid cell.
c      'nearest', -180.0, 180.0, 360.0:
c           these 4 arguments define the domain that will be retrieved 
c           (and also the order of retrieval).  In this case all data 
c           with a longitude coordinate roughly in the range -180 to 180 
c           will be extracted, starting near -180.  The value 360.0 
c           indicates that this coordinate is cyclical with a period of 
c           360, so that EzGet will recognize equivalences such as 
c           0. = 360. = -360. and -90. = 270.  Note that with 'nearest'
c           specified and with the range covering a complete cycle, 
c           EzGet may shift the domain slightly so as to prevent a grid
c           cell near -180.0 (=180.0) from being split across the 
c           boundary.  If, for example, the center of grid cells
c           are located at -180, -170, -160, ... 170, then EzGet will 
c           shift the requested domain to -185 to 175 so the grid cell
c           at -180 will not be split. 
c
c     and in the second call to defdim above the arguments indicate the
c     following:
c      'cosine': this option for controlling creation of weights
c           specifies that weights be generated equal approximately
c           to the cosine of latitude (i.e.,  abs(sin(bdry1) - 
c           sin(bdry2)), where bdry1 and bdry2 are edges of the latitude
c           grid-cell, assumed to be half-way between grid-cell centers)
c           be assigned to all elements (except the weight will be 0.0 
c           for grid cells with missing data).  Another option 
c           ('gaussian') would be appropriate for spectral models with 
c           gaussian grids (as opposed to the evenly spaced grid 
c           accessed here).
c      'range', 90.0, -90.0, 0.0: 
c           these 4 arguments define the domain that will be retrieved 
c           (and also the order of retrieval).  In this case all data 
c           with latitude coordinates in the range 90 to -90 will
c           be retrieved, starting near 90 (i.e. data will be retrieved
c           from north to south).  The value 0.0 indicates that this
c           coordinate is not cyclical.
c  
c     and in the third call to defdim above the fourth argument 
c     indicates the following:
c      'unit': this option for controlling creation of weights specifies 
c           that unit weight should be given to each element of this 
c           dimension (except the weight will be 0. for grid cells with 
c           missing data). 
c
c     -------------------------------------------
c     Extract variable 1 from file and create missing data mask.
c     The lengths of the longitude, latitude and time dimensions of 
c     the data stored in the file are unknown, so initialize the
c     "expected" dimension lengths to 0.  Return the actual dimensions
c     of the retrieved data.

      lons = 0
      lats = 0
      mons = 12
      i4 = 0
      call getdata(1,  nlon,nlat,nmon,n4,  lons,lats,mons,i4,
     &             wtsmask, adata)

c     where on calling getdata the arguments have been defined as 
c     follows:
c
c      1: The first argument has the value 1 and indicates that data 
c         will be extracted for variable 1 (defined in defvar above).
c      nlon, nlat, nmon, n4:
c         These are the declared dimensions of the arrays, wtsmask and
c         adata.  Note that n4=0 because these arrays have only 3
c         dimensions.
c      lons, lats, mons, i4:
c         In this example, the user does not know how large the actual
c         longitude and latitude dimensions of the data being retrieved 
c         will be, so these are set to 0.  The time dimension is 
c         expected to be 12 (months) as specified in the earlier call to 
c         defdim.  In general, if the user knows the size of the domain 
c         to be retrieved, it is usually prudent to set these arguments 
c         to the expected size of the domain, because then EzGet can 
c         error exit if the actual size differs from what is expected.
c
c     On return from getdata:
c
c       1, nlon, nlat, nmon, n4:
c          The first 5 arguments will be unchanged.
c       lons, lats, mons, i4:
c          returns the actual dimensions of the data retrieved.
c       wtsmask: returns the weights associated with the data.
c          In this example, the weights will either be 0.0 or
c          proportional to grid cell area, depending on whether or not 
c          the data are missing.
c       adata: returns the extracted data.

c     -------------------------------------------
c     Retrieve longitude and latitude coordinates of variable 1 (in 
c     the same order as the retrieved data).

      call getcoord(1, 1, alon)
      call getcoord(1, 2, alat)

c     where in the first call to getcoord above, the arguments indicate
c     the following:
c
c      1: The 1st argument has the value 1 and indicates that we want
c           to obtain the coordinates of data extracted for variable 1 
c           (as defined in the earlier call to defvar).
c      1: The 2nd argument has the value 1 and indicates that we want
c           to obtain the coordinates for the 1st dimension (as defined
c           in the earlier call to defdim).
c      alon: returns the coordinate values for the 1st dimension of
c           of variable 1 (i.e., 'longitude', as defined by defdim).

c     -------------------------------------------
c     Close all files opened by EzGet.

      call closeget

c     -------------------------------------------
c     Write data retrieved for last month.  Note that any missing data 
c     will have been assigned the value 1.e20, so when written with
c     f8.3 format will appear as '********':

      write(*,'("longitudes = " / (10f8.3))') (alon(i), i=1,lons)
    
      write(*,'(/ "month = ", i3)') mons

      do 100 j=1,lats
        write(*,'(/"latitude = ", f8.3)') alat(j)
        write(*,'(10f8.3)') (adata(i,j,mons), i=1,lons)
  100 continue

      end
