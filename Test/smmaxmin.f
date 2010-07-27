c f77 -e smmaxmin.f -L$PCMDI/lib -lezget -lcdms -ldrs -L/usr/local/netcdf-2.3.2/libsrc -lnetcdf -o smmaxmin
c f77 -e -g smmaxmin.f -L/zooks0/ktaylor/pcmdi/util/ezget -lezgetdebug -L$PCMDI/lib -lcdms -ldrs -L/usr/local/netcdf-2.3.2/libsrc -lnetcdf -o smmaxmin
c  IBM
c xlf smmaxmin.f  -L/u/williams/ezget -lezget -L/u/williams/devel/cdms/lib -lcdms -L/u/williams/drs/lib -ldrs -L/u/williams/netcdf/libsrc -lnetcdf -o smmaxmin
c   HP
c fort77  +U77 smmaxmin.f -L/afs/nersc.gov/u/ahs/u4160/ezget -lezget -L/afs/nersc.gov/software/pub/drs/lib -lcdms -ldrs -lm -L/afs/nersc.gov/projects/graphics/netcdf/2.3.2pl4/lib -lnetcdf -o smmaxmin
c   sun/SOLARIS
c f77 -e smmaxmin.f -L/pcmdi/cirrus/ktaylor/ezget/solaris -lezget -L/pcmdi/cirrus/drs_solaris -lcdms -ldrs -L/pcmdi/cirrus/netcdf-2.3.2/libsrc -lnetcdf -o smmaxmin
c   SGI
c f77 smmaxmin.f -L/pcmdi/ktaylor/pcmdi/util/ezget/sgi -lezget -L/usr/local/cdms/lib -lcdms -L/usr/local/drs/drslib -ldrs -L/usr/local/lib -lnetcdf -o smmaxmin


      program smmaxmin

c   This program
c     ** obtains the structure of a variable (the number and length
c        of each of its dimensions).
c     ** finds the maximum and minimum value stored for the variable
c        (but skips 'missing' data).
c     ** prints out the information retrieved.
c
c   This program dynamically allocates enough memory to accomodate 
c     the retrieved data. It also avoids creating a data mask.

      pointer (ptadata, adata)
      real adata(*)
      real begdom(4), enddom(4), rmax, rmin, err
      integer ldim(4), isize, ndim, n
      character*16 dimnames(4)
      
c     -------------------------------------------
c     Initialize EzGet:

      call initget

c     -------------------------------------------
c     Define "missing" data value.  (Neither of these calls is 
c     actually necessary, since the default input and output missing
c     value is 1.0e20.

      call defmisc('input missing value', 'real', 1.0e20) 
      call defmisc('output missing value', 'real', 1.0e20) 

c     -------------------------------------------
c     Define variable 1 as 'tas' and indicate path/filename where 
c     it is stored.  In subsequent calls to EzGet (e.g., defdim
c     getdata, and getcoord), this field will be referred to by
c     the index assigned here (i.e., 1).

      call defvar(1, 'tas', '/scratch/staff/lisa/amip_obs/nasa-amip_t') 

c     -------------------------------------------
c     Obtain from the file the number of dimensions, the names of the
c     dimensions and the full domain of each dimension of the variable 
c     (as originally stored).

      call domain(1, ndim, dimnames, begdom, enddom)

c     -------------------------------------------
c     Define domain for variable 1:  retrieve all data stored

      do 100 n=1,ndim
         call defdim(1, n, dimnames(n), 'unit', 'as saved', 0., 0., 0.)
  100 continue

c     -------------------------------------------
c     Retrieve the length of each dimension of the variable and the 
c     total length of the vector of data that will be extracted:

      call shape(1, ldim(1), ldim(2), ldim(3), ldim(4), isize)

c     where
c     1:    The first argument has the value 1 and specifies that 
c           the dimension lengths that will be returned for variable
c           1 are being requested.
c     ldim(1), ldim(2), etc.:  return the lengths of dimensions 1, 2, 3,
c           and 4, respectively.
c     isize:   returns the size of the array needed to accomodate the
c           retrieved data (which may be different from the number
c           of elements stored in the file for this variable).

c     -------------------------------------------
c     Allocate memory for array. 
c     Note: some platforms may have slightly different functions
c         for allocating memory dynamically.

      ptadata = malloc(isize*4)

c     where
c     ptadata:  is a pointer to array adata as declared at the beginning
c           of this program.
c     isize:    specifies how many words to allocate for 
c           array adata.

c     -------------------------------------------
c     Inform EzGet of array size. 

      call defmisc('data size', 'integer', isize) 

c     -------------------------------------------
c     Extract variable 1 from file.

      call getfield(1, adata)

c     where
c     1:    The first argument specifes that data should be retrieved
c           from variable 1.
c     adata:  returns the extracted data.

c     -------------------------------------------
c     Close all files opened by EzGet.

      call closeget

c     -------------------------------------------
c     Find the maximum and minimum values retrieved (skipping "missing"
c     values).

      rmax = -1.e40
      rmin =  1.e40

      do 200 n=1,isize
         if (abs(adata(n)-1.0e20) .gt. 1.e15) then 
            rmax = amax1(rmax, adata(n))
            rmin = amin1(rmin, adata(n))
         endif
  200 continue

c     -------------------------------------------
c     Release memory allocated for adata.
c     Note: some platforms may have slightly different functions
c         for releasing memory.

      err = free(ptadata)

c     where
c     ptadata:  is the pointer to array adata, as declared at the 
c     beginning of this program. 

c     -------------------------------------------
c     Report the structure of the data and the maximum and minimum
c     values stored.

      write(*,'("         Data Structure" / 
     &          " Dimension Name       Length " / /
     &      4(a16, i10 /))') 
     &      (dimnames(n), ldim(n), n=1,ndim)

      if (rmax .ge. rmin) then
        write(*,'(" Maximum value found: ", 1pe14.5)') rmax
        write(*,'(" Minimum value found: ", 1pe14.5)') rmin
      else
        write(*,'(" Data set includes only ''missing'' data.")')
      endif

      end
