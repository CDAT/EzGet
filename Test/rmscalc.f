c f77 -e rmscalc.f -L$PCMDI/lib -lezget -lcdms -ldrs -L/usr/local/netcdf-2.3.2/libsrc -lnetcdf -o rmscalc
c f77 -e -g rmscalc.f -L/zooks0/ktaylor/pcmdi/util/ezget -lezgetdebug -L$PCMDI/lib -lcdms -ldrs -L/usr/local/netcdf-2.3.2/libsrc -lnetcdf -o rmscalc
c  IBM
c xlf rmscalc.f  -L/u/williams/ezget -lezget -L/u/williams/devel/cdms/lib -lcdms -L/u/williams/drs/lib -ldrs -L/u/williams/netcdf/libsrc -lnetcdf -o rmscalc
c   HP
c fort77  +U77 rmscalc.f -L/afs/nersc.gov/u/ahs/u4160/ezget -lezget -L/afs/nersc.gov/software/pub/drs/lib -lcdms -ldrs -lm -L/afs/nersc.gov/projects/graphics/netcdf/2.3.2pl4/lib -lnetcdf -o rmscalc
c   sun/SOLARIS
c f77 -e rmscalc.f -L/pcmdi/cirrus/ktaylor/ezget/solaris -lezget -L/pcmdi/cirrus/drs_solaris -lcdms -ldrs -L/pcmdi/cirrus/netcdf-2.3.2/libsrc -lnetcdf -o rmscalc
c   SGI
c f77 rmscalc.f -L/pcmdi/ktaylor/pcmdi/util/ezget/sgi -lezget -L/usr/local/cdms/lib -lcdms -L/usr/local/drs/drslib -ldrs -L/usr/local/lib -lnetcdf -o rmscalc

      program rmscalc
      
c   This program loops through the AMIP models and computes the area- 
c   weighted RMS difference between the modeled and observed annual mean
c   surface air temperature over North America (north of the Tropic of 
c   Cancer) for the year 1988.  Only grid cells with observed values for
c   every month of this year are included in computing the RMS 
c   difference.  Because some models may have unrealistic 
c   representations of the continental geography, some data from the 
c   model's North American grid may not map onto the observed North 
c   American grid. Any grid cells where either observations or model 
c   data are missing are excluded from the RMS difference computed 
c   here.

      parameter (nlon=35, nlat=20, nmon=12, n4=0, nmods=28)

      real adata(nlon, nlat, nmon), wtsmask(nlon, nlat, nmon),
     &     yrobs(nlon,nlat), yrmodel(nlon,nlat), 
     &     wtobs(nlon,nlat), wtmodel(nlon,nlat)
      double precision asum, rmsdiff, wtsum
      integer lons, lats, mons, i4, i, j, m, n, mm
      character*3 modlname(nmods)
      data modlname /   'bmr', 'ccc', 'col', 'cnr', 'csi', 'csu', 'der',
     &    'dnm', 'ecm', 'gfd', 'gis', 'gla', 'gsf', 'iap', 'jma', 'lmd',
     &    'mgo', 'mpi', 'mri', 'nca', 'nmc', 'nrl', 'sng', 'sun', 'ucl', 
     &    'uiu', 'ukm', 'yon' /

c      Note:  the 'rpn' and 'uga' models were left out of the above list
c      because the variable 'tas' is not available for these 
c      models.

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
     &                                             -190.0, -40.0, 360.0)
      call defdim(1, 2, 'latitude', 'cosine', 'range', 23.5, 90.0, 0.0)
      call defdim(1, 3, 'time',    'unit', 'nearest', 109., 120., 0.0)

c     -------------------------------------------
c     Define variable 2 as the geography data on the observational grid,
c     and specify that North American data should be selected:

      call defvar(2, 'sftbyrgn', '/amipsp/drs/sftbyrgn/sftbyrgn_gla') 
      call defgeog(1, 'in', 2, 'North America')

c     -------------------------------------------
c     Extract variable 1 (observed field). 

      lons = 0
      lats = 0
      mons = nmon
      i4 = 0
      call getdata(1,  nlon,nlat,nmon,n4,  lons,lats,mons,i4,
     &             wtsmask, adata)

c     -------------------------------------------
c     Compute annual mean at each grid cell, but only if all 12 months
c     of observed data are available.

      do 300 j=1,nlat
        do 200 i=1,nlon

          asum = 0.0
          wtsum = 0.0
          mm = 0

          do 100 m=1,12

            if (wtsmask(i,j,m) .gt. 0.) then
              mm = mm + 1
              wtsum = wtsum + wtsmask(i,j,m)
              asum = asum + adata(i,j,m)
            endif

  100     continue

          if (mm .eq. 12) then
            yrobs(i,j) = asum/12.
            wtobs(i,j) = wtsum/12.
          else
            yrobs(i,j) = 0.0
            wtobs(i,j) = 0.0
          endif

  200   continue
  300 continue
    
c     -------------------------------------------
c     Loop through AMIP models and compute RMS differences between 
c     modeled and observed fields.

      do 1000 n=1,nmods

c       -----------------------------------------
c       Define variable 3 as the model simulated surface air temperature 
c       and specify the type of grid this model has: 

        call defvar(3, 'tas', '/amipsp/drs/tas/tas_'//modlname(n))

        call defdim(3, 1, 'longitude', modlname(n), 'nearest',
     &                                             0.0, 0.0, 360.0)
        call defdim(3, 2, 'latitude', modlname(n), 'nearest',
     &                                               0.0, 0.0, 0.0)
        call defdim(3, 3, 'time',   'unit', 'nearest', 109., 120., 0.0)

c       >>> Note that the domains specified above for latitude and 
c       >>> longitude are ignored by EzGet because the data for this
c       >>> variable will be regridded to a target grid, which 
c       >>> determines the domain.

c       >>> Note that the name of the model can be used to specify the 
c       >>> type of longitude and latitude weights that will be 
c       >>> generated by EzGet.  EzGet contains a table that allows
c       >>> it to generate the proper weights proportional to grid cell
c       >>> area, if the model name is a standard AMIP or PMIP name. 
 
c       >>> Note that for AMIP data, months 109 through 120 are the 
c       >>> months of the calendar year 1988.
  
c       -----------------------------------------
c       Define variable 4 as the geography data on the model grid, and
c       specify that North American data should be selected for
c       variable 3.

        call defvar(4, 'sftbyrgn', 
     &              '/amipsp/drs/sftbyrgn/sftbyrgn_'//modlname(n)) 
        call defgeog(3, 'in', 4, 'North America')

c       -----------------------------------------
c       Instruct EzGet to map modeled field to observational grid upon
c       retrieving data. 

        call defregrd(3, 'to', 1,  'area-weighted', 
     &                                0, 0.0, 0.0,  0, 0.0, 0.0)

c       -----------------------------------------
c       Extract variable 3 and map to observed grid. 

        lons = 0
        lats = 0
        mons = nmon
        i4 = 0
        call getdata(3,  nlon,nlat,nmon,n4,  lons,lats,mons,i4,
     &             wtsmask, adata)

c       -----------------------------------------
c       Compute annual mean at each grid cell, but only if all 12 months
c       of data are available.

        do 600 j=1,nlat
          do 500 i=1,nlon

            asum = 0.0
            wtsum = 0.0
            mm = 0

            do 400 m=1,12

              if (wtsmask(i,j,m) .gt. 0.) then
                mm = mm + 1
                wtsum = wtsum + wtsmask(i,j,m)
                asum = asum + adata(i,j,m)
              endif

  400       continue

            if (mm .eq. 12) then
              yrmodel(i,j) = asum/12.
              wtmodel(i,j) = wtsum/12.
            else
              yrmodel(i,j) = 0.0
              wtmodel(i,j) = 0.0
            endif

  500     continue
  600   continue

c       -----------------------------------------
c       Compute area-weighted RMS difference between model and observed
c       fields.

        wtsum = 0.0
        rmsdiff = 0.0

          do 800 j=1,nlat
          do 700 i=1,nlon

            if (wtmodel(i,j) .gt. 0.0) then
              wtsum = wtsum + wtobs(i,j)
              rmsdiff = rmsdiff + 
     &                      wtobs(i,j)*((yrobs(i,j)-yrmodel(i,j))**2)
            endif
  
  700     continue
  800   continue

c       -----------------------------------------
c       Write area-weighted RMS difference between model and observed
c       fields.

        if (wtsum .gt. 0.0) then 

          rmsdiff = dsqrt(rmsdiff/wtsum)
          write(*,'(1x, a8, f12.7, f12.3)') modlname(n), wtsum, rmsdiff

        else

          write(*,'(1x, a8, " data missing")') modlname(n)

        endif

 1000 continue

c     -------------------------------------------
c     Close all files opened by EzGet.

      call closeget

      end
