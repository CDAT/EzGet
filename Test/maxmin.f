c f77 -e maxmin.f -L$PCMDI/lib -lezget -lcdms -ldrs -L/usr/local/netcdf-2.3.2/libsrc -lnetcdf -o maxmin
c f77 -e -g maxmin.f -L/zooks0/ktaylor/pcmdi/util/ezget -lezgetdebug -L$PCMDI/lib -lcdms -ldrs -L/usr/local/netcdf-2.3.2/libsrc -lnetcdf -o maxmin
c  IBM
c xlf maxmin.f  -L/u/williams/ezget -lezget -L/u/williams/devel/cdms/lib -lcdms -L/u/williams/drs/lib -ldrs -L/u/williams/netcdf/libsrc -lnetcdf -o maxmin
c   HP
c fort77  +U77 maxmin.f -L/afs/nersc.gov/u/ahs/u4160/ezget -lezget -L/afs/nersc.gov/software/pub/drs/lib -lcdms -ldrs -lm -L/afs/nersc.gov/projects/graphics/netcdf/2.3.2pl4/lib -lnetcdf -o maxmin
c   sun/SOLARIS
c f77 -e maxmin.f -L/pcmdi/cirrus/ktaylor/ezget/solaris -lezget -L/pcmdi/cirrus/drs_solaris -lcdms -ldrs -L/pcmdi/cirrus/netcdf-2.3.2/libsrc -lnetcdf -o maxmin
c   SGI
c f77 maxmin.f -L/pcmdi/ktaylor/pcmdi/util/ezget/sgi -lezget -L/usr/local/cdms/lib -lcdms -L/usr/local/drs/drslib -ldrs -L/usr/local/lib -lnetcdf -o maxmin

      program maxmin

c   This program
c     ** obtains the structure of a variable (the number and length
c        of each of its dimensions).
c     ** finds the maximum and minimum value stored for the variable
c        (but skips 'missing' data).
c     ** prints out the information retrieved.
c
c   This program will error exit if a data set's size exceeds maxsize 
c     (declared in the parameter statement below) 

      parameter (maxsize=500000)

      real adata(maxsize), wtsmask(maxsize)
      real begdom(4), enddom(4), rmax, rmin
      integer ldim(4), isize, ndim, n
      character*16 dimnames(4)
      
c     -------------------------------------------
c     Initialize EzGet:

      call initget

c     -------------------------------------------
c     Define "missing" data value. 

      call defmisc('input missing value', 'real', 1.0e20) 

c     -------------------------------------------
c     Inform EzGet of maximum array size treatable by this program. 

      call defmisc('data size', 'integer', maxsize) 

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

c     where
c      1:    The first argument has the value 1 and indicates that the
c            you want to retrieve dimension information for variable 1.
c      ndim: This argument returns the number of dimensions there are 
c            for variable 1.
c      dimnames: This vector argument returns the dimension names.
c      begdom:  This vector argument returns the first coordinate value
c            stored for each dimension.
c      enddom:  This vector argument returns the last coordinate value
c            stored for each dimension.

c     -------------------------------------------
c     Define domain for variable 1:  retrieve all data stored

      do 100 n=1,ndim
         call defdim(1, n, dimnames(n), 'unit', 'as saved', 0., 0., 0.)
  100 continue

c     -------------------------------------------
c     Extract variable 1 from file and create missing data mask.

      call getvdata(1, wtsmask, adata)

c     where
c     1:     The first argument has the value 1 and indicates that the
c            you want to retrieve data from variable 1.
c     wtsmask: returns the weights associated with the data.
c     adata:   returns the extracted data.

c     -------------------------------------------
c     Retrieve the length of each dimension of the variable and the 
c     total length of the vector of data retrieved:

      call lendims(1, ldim(1), ldim(2), ldim(3), ldim(4), isize)

c     where
c     1:    The first argument has the value 1 and specifies that the
c           the dimension lengths for variable 1 should be obtained.
c     ldim(1), ldim(2), ldim(3), ldim(4):  return the lengths of
c           dimensions 1, 2, 3, and 4, respectively.
c     isize:  isize = ldim(1)*ldim(2)*ldim(3)*ldim(4), but if any
c     dimension is length 0, a 1 is substituted in the above formula.

c     -------------------------------------------
c     Close all files opened by EzGet.

      call closeget

c     -------------------------------------------
c     Find the maximum and minimum values retrieved (skipping "missing"
c     values).

      rmax = -1.e40
      rmin =  1.e40

      do 200 n=1,isize
         if (wtsmask(n) .gt. 0.0) then
            rmax = amax1(rmax, adata(n))
            rmin = amin1(rmin, adata(n))
         endif
  200 continue

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

