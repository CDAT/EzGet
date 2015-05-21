c
c    16 May 1997
c
c    Karl E. Taylor
c
c     compile these subroutines
c
c   Preprocessor "includes":
c      These subroutines rely on the preprocessor to include cliches 
c         containing dimension and common statements (hence the .F 
c         suffix on ezget.F).  Be sure to properly set your path 
c         in the Makefile (or change the path in all the following
c          include statements) so that these cliches can be found.
c       All but 2 include files should be stored in a subdirectory
c          called "cliche".  The files should be accessible through
c          statements of the form:
c
c     #include "cdimtbl"
c     #include "cvartbl"
c     #include "cdomain"
c     #include "ciotbl"
c     #include "cregddim"
c     #include "cscalars"
c
c     In addition to the above, the following 2 cliches must be 
c       accessible statements of the form:
c

c     #include "usersget"
c     #include "drsdef.h"
c
c
c   ***  Procedure for Sun
c
c f77 -c -g  ezget.F -I$PCMDI/include/ezget  -o ezget.o
c f77 -c -g  ezget.F -o ezget_dbg.o
c f77 -c ezget1.F -o ezget.o
c
c   create library
c
c     ar rcv libezget.a ezget.o
c
c    create table of contents
c
c    ranlib libezget.a


c
c    compile and load program that calls libezget.a
c
c    f77 -g -e sampledrs.F -L$PCMDI -lezget -ldrs -o sampledrs
c
c   *** Procedure for CRAY-2
c
c cf77  ezget.F  -o ezget.o
c cf77 -g ezget.F  -o ezget.o
c cf77 -Zv ezget.F  -o ezget.o
c
c    Create a makefile: (name it Makefile)
c
c   FFLAGS = -Zv (optimize) OR -g (debug)
c   AR = bld
c
c   .PRECIOUS: libezget.a
c
c   libezget.a: libezget.a(ezget.o)
c
c     then  to run the Makefile, enter
c   make libezget.a
c
c   obtain from cfs:
c 
c      /017374/drs/unicos/dev/libdrs.a
c      /017374/drs/unicos/dev/drsdef.h
c   
c    compile and load program that calls libezget.a
c
c    cf77 -g sampledrs.F -lezget -ldrs -o sampledrs
c    cf77 -Zv sampledrs.F -lezget -ldrs -o sampledrs
c

c 3456789012345678901234567890123456789012345678901234567890123456789012
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     integer function applywts
c     -------
c
c     Description:
c     -----------
c    This function intersects the input mask with "areas" determined
c      by the dimension weights.
c
c
c     Usage:
c     ------
c
c      ierr = applywts(ivar, amask)
c
c
c     Exit Status:
c     -----------
c     0 if success
c     1 if an error occurred
c
c     ------
c     Date: 1/19/92
c     ----
c
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

      integer function applywts(ivar, amask)

      implicit none
#include "usersget"
#include "cvartbl"
#include "cdimtbl"
     
      integer ivar
      real amask(*)

      real wt2, wt3, wt4
      integer ipt1, ipt2, ipt3, ipt4, len1, len2, len3, len4, i, j, k,
     + l, n, ion1, ion2, ion3, ion4

      applywts = 0

c     get info. from variable table


      if (iptdimrg(1,ivar) .eq. 0) then
          len1 = 1
          ipt1 = lencoord+1
          ion1 = 0
      else
          ipt1 = iptdimrg(1,ivar)
          len1 = idimlen(ipt1)
          ipt1 = iptwts(ipt1)
          ion1 = 1
          if (ipt1 .lt. 0) then
            ipt1 = lencoord+1
            ion1 = 0
          endif
      endif

      if (iptdimrg(2,ivar) .eq. 0) then
          len2 = 1
          ipt2 = lencoord+1
          ion2 = 0
      else
          ipt2 = iptdimrg(2,ivar)
          len2 = idimlen(ipt2)
          ipt2 = iptwts(ipt2)
          ion2 = 1
          if (ipt2 .lt. 0) then
            ipt2 = lencoord+1
            ion2 = 0
          endif
      endif

      if (iptdimrg(3,ivar) .eq. 0) then
          len3 = 1
          ipt3 = lencoord+1
          ion3 = 0
      else
          ipt3 = iptdimrg(3,ivar)
          len3 = idimlen(ipt3)
          ipt3 = iptwts(ipt3)
          ion3 = 1
          if (ipt3 .lt. 0) then
            ipt3 = lencoord+1
            ion3 = 0
          endif
      endif

      if (iptdimrg(4,ivar) .eq. 0) then
          len4 = 1
          ipt4 = lencoord+1
          ion4 = 0
      else
          ipt4 = iptdimrg(4,ivar)
          len4 = idimlen(ipt4)
          ipt4 = iptwts(ipt4)
          ion4 = 1
          if (ipt4 .lt. 0) then
            ipt4 = lencoord+1
            ion4 = 0
          endif
      endif

      n = 0
      do 400 l=1,len4
         wt4 = coordvls(ipt4-1+l*ion4)
         do 300 k=1,len3
            wt3 = wt4*coordvls(ipt3-1+k*ion3)
            do 200 j=1,len2
               wt2 = wt3*coordvls(ipt2-1+j*ion2)
               do 100 i=1,len1
                  n = n + 1
                  amask(n) = wt2*coordvls(ipt1-1+i*ion1)*amask(n)
  100           continue
  200        continue
  300    continue
  400 continue

      return
      end




c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     logical function caseindp
c     -------
c
c     Description:
c     -----------
c    This function tests for character string equivalence
c
c     Usage:
c     ------
c
c      if (caseindp(a, b)) then
c
c     ------
c     Date: 9/17/94
c     ----
c
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      logical function caseindp(a, b)

c    This function tests for character string equivalence
c
c    compare what follows any leading blanks
c    check for equivalence of strings ignoring case.
c    only check through the last character of the shorter string 
c    (where "shorter" is computed after removing leading blanks)
c    The only exception to this is if 1 string is null or completely
c    blank; then equivalence requires both strings to be null or
c    completely blank.   

      implicit none
      character*(*) a, b
      character*1 c, d
      integer i, j, k, m, n, j1, k1, ii, j2, k2

c     look for leading blanks and neglect

      j = len(a)
      k = len(b)

c     check whether both strings are length zero

      if ((j .eq. 0) .or. (k .eq. 0)) then
        if (j .eq. k) then
          caseindp = .true.
        else
          caseindp = .false.
        endif
        return
      endif


      j1 = 0
      do 10 i=1,j
        if (a(i:i) .ne. ' ') go to 20
        j1 = j1 + 1
   10 continue
   20 continue

      k1 = 0
      do 30 i=1,k
        if (b(i:i) .ne. ' ') go to 40
        k1 = k1 + 1
   30 continue
   40 continue

c   check whether either string is completely blank

      if ((j1 .eq. j) .or. (k1 .eq. k)) then
        if ((j1 .eq. j) .and. (k1 .eq. k)) then
          caseindp = .true.
        else
          caseindp = .false.
        endif
        return
      endif


c   look for trailing blanks and neglect
      
      j2 = j
      do 50 i=1,j
        if (a(j-i+1:j-i+1) .ne. ' ') go to 60
        j2 = j2 - 1
   50 continue
   60 continue

      k2 = k
      do 70 i=1,k
        if (b(k-i+1:k-i+1) .ne. ' ') go to 80
        k2 = k2 - 1
   70 continue
   80 continue

      j = j2 - j1
      k = k2 - k1
      ii = min0(j, k)
 
      do 100 i=1,ii
        m = i + j1
        n = i + k1
        c = a(m:m)
        if (lge(c,'A') .and. lle(c,'Z')) c = char(ichar(c)+32)
        d = b(n:n)
        if (lge(d,'A') .and. lle(d,'Z')) d = char(ichar(d)+32)
        if (c .ne. d) then
           caseindp = .false.
           return
        endif
  100 continue

      caseindp = .true. 
      return
      end
          



c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     subroutine closeget 
c     -------
c                
c     Description:
c     -----------
c   This subroutine closes any drs-files opened (and not already closed)
c       by ezget.
c
c     ------
c     Date: 1/19/93
c     ----
c
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine closeget

c      common variables

      implicit none
#include "drsdef.h"
#include "fcddrs.h"


#include "ciotbl"
#include "cscalars"

      integer ierr, i, cllun
      
      do 100 i = 1,maxiodev
         if (iofiluse(i) .ge. 0) ierr = cllun(90+i)
         iofilnam(i) = 'not defined'
         iofiluse(i) = i - 1 - maxiodev
  100 continue

      if (warncnt .gt. 0) then 
        print*, ' '
        if (warncnt .eq. 1) then
           print*,  warncnt, ' warning was encountered in ezget,'
        else
         print*,  warncnt, ' warnings were encountered in ezget,'
        endif
        print*, 'but execution was not halted.'
      endif

      if (errcount .gt. 0) then
        print*, ' '
        if (errcount .eq. 1) then
           print*,  errcount, ' error was encountered in ezget,'
        else
         print*,  errcount, ' errors were encountered in ezget,'
        endif
        print*, 'but execution was not halted.'
      endif

      return
      end

  

c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      subroutine clrtable
c     -------
c
c     Description:
c     -----------
c    This function clears the dimension table.  It erases all 
c        information about each dimension to make room for new
c        dimensions.  It only should be called if you no longer
c        are planning to obtain dimension information about any
c        of the currently defined variables.
c
c
c     Usage:
c     ------
c
c      call clrtable
c
c     ------
c     Date: 1/19/93
c     ----
c
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine clrtable

      implicit none
#include "usersget"
#include "cdimtbl"
#include "cvartbl"

      integer i, j

      do 100 i = 1,maxdims
         dimnam(i) = 'not defined'
         idimlen(i) = 0
         ioriglen(i) = 0
         idimtype(i) = 0
         iptcoord(i) = 0
         iptedge(i) = 0
         iptwts(i) = 0
         dimdoman(i) = 'not def.'
         dimcycle(i) = 0.0
         dimwtflg(i) = 'not def.'
         begbdy(i) = 0.0
         endbdy(i) = 0.0
         felretr(i) = 0.0
         lelretr(i) = 0.0
  100 continue
  

      icount = 0
      ifreecoo = 1
      coordvls(lencoord) = 1.

      do 200 i=1,maxvars
         do 150  j=1,4
            iptdim(j,i) = 0
            iptdimrg(j,i) = 0
  150    continue
  200 continue

      return
      end



c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     subroutine defdim
c     -------
c
c     Description:
c     -----------
c     this subroutine sets controls for dimensions
c
c     the dimensions should be defined in ascending order:
c       i.e., define the second dimension after the first dimension,
c         the third after the second and the fourth after the third.   
c  
c     Usage:
c     ------
c
c      call defdim(ivar, idim, name, weight, domain, bdry1, bdry2,
c                  dcycle)
c
c     if ivar .ne. 0,  then assign the dimension characteristics
c             specified here to variable ivar only.
c     otherwise, assign these dimension characteristics to all
c             variables already defined or subsequently defined.  These
c             characteristics can be subsequently overridden for a 
c             specific variable by a call to defdim with ivar set 
c             appropriately, or for all variables with ivar=0. 
c
c     ivar = a defined variable for which the dimension will be define
c            (if 0, all variables will be assigned this dimension)
c     idim = dimension number (as you wish to assign it to your array
c             not necessarily as it was stored)
c            if 0, reset all dimensions for ivar to 'not defined'
c            if -1, -2, -3, or -4, set the appropriate dimension to
c               not defined (dimension 1 if -1, dimension 2 if -2, etc.)
c     name = dimension name (a character string which will be truncated
c              here to 16 characters if necessary) 
c     weight = a string indicating the type of weighting associated with
c                this dimension:
c              'width' = weight by the width of the grid cell
c              'unit' = weight by 1
c              'equal' = weight by 1/ngrid, where ngrid is the number of
c                         grid points in this dimension
c              'gaussian' = weight by gaussian weights
c              'cosine' = weight by cosine of latitude
c              'month' = weight by the number of days in the month
c              'leapyr' =weight by the number of days in the months of a
c                          leap year.
c              a standard AMIP model name (e.g., 'csu') = weight by 
c                       proper amount for model indicated
c     domain = controls type of domain retrieved.  
c              'range' = retrieve domain specified by bdry1 and bdry2
c              'nearest' = retrieve domain specified by bdry1 and bdry2
c              'assaved' = retrieve full domain stored
c     bdry1 =first boundary of domain to be retrieved (if domain='range' 
c              or domain = 'nearest')
c            if domain .eq. 'assaved', this parameter is ignored.     
c     bdry2=other boundary of domain to be retrieved (if domain='range' 
c              or domain = 'nearest')
c            if domain .eq. 'assaved', this parameter is ignored.     
c     dcycle=the domain extent for cyclical domains such as longitude or
c                   months. 
c
c     ------
c     Date: 9/17/94
c     ----
c
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine defdim(ivar, idim, name, weight, domain, bdry1,
     +         bdry2, dcycle) 

c ???  check the following include statements
      implicit none
#include "usersget"
#include "cvartbl"
#include "cdimtbl"
#include "cscalars"

      integer ivar, idim
      logical caseindp
      character*(*) name, weight, domain
      real bdry1, bdry2, dcycle
                   
      integer  i, m, m1, m2, univunit, i1, i2, ii
      real bdry1l, bdry2l

      character domaincl*8, cwtl*8
                 
      if (ivar .gt. maxvars) then
 
        print*, 'You are attempting to assign a dimension to'
        print*, 'a nonexistent variable, which you identify'
        print*, 'by the integer ', ivar, ' when calling subroutine'
        print*, 'defdim.  Note that no more than, ', maxvars
        print*, 'may defined unless you redefine maxvars in'
        print*, 'common usersget' 
        call errcheck('subroutine defdim', ' ',1)
        return
 
      elseif (ivar .lt. 0) then

        print*, 'You have improperly passed ivar = ', ivar
        print*, 'as the argument that identifies the variable'
        print*, 'for which you want to define a dimension in '
        print*, 'subroutine defdim.  This integer should be .ge. 0'  
        call errcheck('subroutine defdim', ' ',1)
        return

      endif

      if ((idim .gt. 4) .or. (idim .lt. -4)) then
        print*, 'You have improperly passed idim = ', idim
        print*, 'as the argument that identifies the dimension'
        print*, 'you want defined in subroutine defdim.  This integer'
        print*, 'should be in the following range: -5 < idim < 5'
        call errcheck('subroutine defdim', ' ',1)
        return
      endif   

      if (ivar .eq. 0) then
        m1 = 1
        m2 = maxvars
      else
        m1 = ivar
        m2 = ivar
      endif                                                         
                     
      do 80 m=m1,m2
        do 70 i=1,4
          iptdim(i,m) = 0
          iptdimrg(i,m) = 0
   70   continue
   80 continue
                     
      if (idim .gt. 1) then

        do 100 m=m1,m2              

          do 90 i=1,idim-1
          
          if (caseindp(name, vdimname(i,m))) then
            print*, 'You are attempting to name dimension ', idim
            print*, 'to be the same as dimension ', i
            print*, 'Make sure you are defining your dimensions'
            print*, 'in ascending order.  (e.g., define dimension 3'
            print*, 'after defining dimension 1 and 2 and before '
            print*, 'defining dimension 4)'
            print*, ' '
            if (ivar .gt. 0) then
              print*, 'You were defining the dimension improperly for '
              print*, 'variable ', ivar 
              print*, 'which is currently defined as ', varnam(ivar)
              print*, 'in file: ', filnam(ivar)  
            else
              print*, 'It is possible that this same error will be'
              print*, 'encountered for other variables too.'
            endif 
            call errcheck('subroutine defdim', ' ', 1)
            return
          endif

   90     continue

  100   continue

      endif 
             
      if (idim .le. 0) then
        if (idim .eq. 0) then
          i1 = 1
          i2 = 4
        else
          i1 = -idim
          i2 = -idim
        endif
        do 140 m=m1,m2
          do 130 i=i1,i2
            vdimname(i,m) = 'not defined'
            vdimwt(i,m) = 'not def.'
            vdimdom(i,m) = 'not def.'
            vdimbdry(1,i,m) = 0.   
            vdimbdry(2,i,m) = 0.   
            ivdimpln(1,i,m) = 0.   
            ivdimpln(2,i,m) = 0.   
            vdimcycl(i,m) = 0.
  130     continue
  140   continue
        return
      endif

      if (len(weight) .lt. 3) then
        print*, weight, ' not recognized as a valid weighting option'
        print*, 'for dimension ', idim, ' named ', name
        print*, 'You were defining the dimension improperly for '
        if (ivar .gt. 0) then    
          print*, 'variable ', ivar 
          print*, 'which is currently defined as ', varnam(ivar)
          print*, 'in file: ', filnam(ivar)  
        else
          print*, 'all variables.'
        endif  
        call errcheck('subroutine defdim', ' ', 1)
        return
      endif                     

      if (caseindp(weight, 'wid')) then
        cwtl = 'width'
      elseif (caseindp(weight, 'uni')) then    
        cwtl = 'unit'
      elseif (caseindp(weight, 'equ')) then    
        cwtl = 'equal'
      elseif (caseindp(weight, 'gau')) then    
        cwtl = 'gaussian'
      elseif (caseindp(weight, 'cos')) then    
        cwtl = 'cosine'
      elseif (caseindp(weight, 'mon')) then    
        cwtl = 'month'  
      elseif (caseindp(weight, 'lea')) then    
        cwtl = 'leapyr'
      elseif  (caseindp(weight, 'csu') .or. caseindp(weight, 'ucl')
     +  .or. caseindp(weight, 'lmd') .or. caseindp(weight, 'lmc')) then
        if (univunit(name) .eq. intlong) then
          cwtl = 'width'
        elseif (univunit(name) .eq. intlat) then
          cwtl = weight(1:3)
        else
          print*, weight, ' not recognized as a valid weighting option'
          print*, 'for dimension ', idim, ' named ', name
          print*, 'you may only specify this weighting option for'
          print*, 'latitude or longitude dimensions.'
          print*, 'You were defining the dimension improperly for '
          if (ivar .gt. 0) then    
            print*, 'variable ', ivar 
            print*, 'which is currently defined as ', varnam(ivar)
            print*, 'in file: ', filnam(ivar)  
          else
            print*, 'all variables.'
          endif  
          call errcheck('subroutine defdim', ' ', 1)
          return
        endif
      elseif (caseindp(weight, 'bmr') .or. caseindp(weight, 'ccc') 
     +   .or. caseindp(weight, 'cnr') .or. caseindp(weight, 'col') 
     +   .or. caseindp(weight, 'csi') .or. caseindp(weight, 'der') 
     +   .or. caseindp(weight, 'ecm') .or. caseindp(weight, 'gfd') 
     +   .or. caseindp(weight, 'mgo') .or. caseindp(weight, 'mpi')
     +   .or. caseindp(weight, 'nca') .or. caseindp(weight, 'nmc') 
     +   .or. caseindp(weight, 'rpn') .or. caseindp(weight, 'sun') 
     +   .or. caseindp(weight, 'uga') .or. caseindp(weight, 'sng')
     +   .or. caseindp(weight, 'gen') .or. caseindp(weight, 'cer')
     +   .or. caseindp(weight, 'ccm') .or. caseindp(weight, 'ccs')
     +   .or. caseindp(weight, 'ech') .or. caseindp(weight, 'nce')) then
        if (univunit(name) .eq. intlong) then
          cwtl = 'width'
        elseif (univunit(name) .eq. intlat) then
          cwtl = 'gaussian'
        else
          print*, weight, ' not recognized as a valid weighting option'
          print*, 'for dimension ', idim, ' named ', name
          print*, 'you may only specify this weighting option for'
          print*, 'latitude or longitude dimensions.'
          print*, 'You were defining the dimension improperly for '
          if (ivar .gt. 0) then    
            print*, 'variable ', ivar 
            print*, 'which is currently defined as ', varnam(ivar)
            print*, 'in file: ', filnam(ivar)  
          else
            print*, 'all variables.'
          endif  
          call errcheck('subroutine defdim', ' ', 1)
          return
        endif
      elseif (caseindp(weight, 'dnm') .or. caseindp(weight, 'gis') 
     +   .or. caseindp(weight, 'gla') .or. caseindp(weight, 'gsf') 
     +   .or. caseindp(weight, 'iap') .or. caseindp(weight, 'mri') 
     +   .or. caseindp(weight, 'uiu') 
     +   .or. caseindp(weight, 'jma') .or. caseindp(weight, 'nrl') 
     +   .or. caseindp(weight, 'ukm') .or. caseindp(weight, 'yon')) 
     +      then
        if (univunit(name) .eq. intlong) then
          cwtl = 'width'
        elseif (univunit(name) .eq. intlat) then
          cwtl = 'cosine'
        else
          print*, weight, ' not recognized as a valid weighting option'
          print*, 'for dimension ', idim, ' named ', name
          print*, 'you may only specify this weighting option for'
          print*, 'latitude or longitude dimensions.'
          print*, 'You were defining the dimension improperly for '
          if (ivar .gt. 0) then    
            print*, 'variable ', ivar 
            print*, 'which is currently defined as ', varnam(ivar)
            print*, 'in file: ', filnam(ivar)  
          else
            print*, 'all variables.'
          endif  
          call errcheck('subroutine defdim', ' ', 1)
          return
        endif
      else
        print*, weight, ' not recognized as a valid weighting option'
        print*, 'for dimension ', idim, ' named ', name
        print*, 'You were defining the dimension improperly for '
        if (ivar .gt. 0) then    
          print*, 'variable ', ivar 
          print*, 'which is currently defined as ', varnam(ivar)
          print*, 'in file: ', filnam(ivar)  
        else
          print*, 'all variables.'
        endif  
        call errcheck('subroutine defdim', ' ', 1)
        return
      endif 

      if (len(domain) .lt. 3) then 
        print*, domain, ' not recognized as a valid domain option'
        print*, 'for dimension ', idim, ' named ', name
        print*, 'You were defining the dimension improperly for '
        if (ivar .gt. 0) then    
          print*, 'variable ', ivar 
          print*, 'which is currently defined as ', varnam(ivar)
          print*, 'in file: ', filnam(ivar)  
        else
          print*, 'all variables.'
        endif  
        call errcheck('subroutine defdim', ' ', 1)
        return
      endif

      bdry1l = bdry1
      bdry2l = bdry2
      if (caseindp(domain, 'ran')) then
        domaincl = 'range'
      elseif (caseindp(domain, 'nea')) then
        domaincl = 'nearest'          
      elseif (caseindp(domain, 'ass') .or. caseindp(domain, 'as s'))
     +    then
        domaincl = 'assaved'
        bdry1l = 0.0
        bdry2l = 0.0
      else
        print*, domain, ' not recognized as a valid domain option'
        print*, 'for dimension ', idim, ' named ', name
        call errcheck('subroutine defdim', ' ', 1)
        return
      endif    

      ii = min0(16, len(name))                           
      do 200 m=m1,m2

          vdimname(idim,m) = ' '
          vdimname(idim,m) = name(1:ii)
          vdimwt(idim,m) = cwtl
          vdimdom(idim,m) = domaincl
          vdimbdry(1,idim,m) = bdry1l 
          vdimbdry(2,idim,m) = bdry2l
          vdimcycl(idim,m) = dcycle
          ivdimpln(1,idim,m) = 0      
          ivdimpln(2,idim,m) = 0      
            
          ivarndim(m) = 0
          ivarlen(m) = 0
          iptdim(idim,m) = 0
          iptdimrg(idim,m) = 0

  200 continue

      return
      end       



c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     subroutine defdimi
c     -------
c
c     Description:
c     -----------
c     this subroutine sets controls for dimensions
c
c     the dimensions should be defined in ascending order:
c       i.e., define the second dimension after the first dimension,
c         the third after the second and the fourth after the third.   
c  
c     Usage:
c     ------
c
c      call defdimi(ivar, idim, name, weight, indx1, indx2)
c
c     if ivar .ne. 0,  then assign the dimension characteristics
c             specified here to variable ivar only.
c     otherwise, assign these dimension characteristics to all
c             variables already defined or subsequently defined.  These
c             characteristics can be subsequently overridden for a 
c             specific variable by a call to defdimi with ivar set 
c             appropriately, or for all variables with ivar=0. 
c
c     ivar = a defined variable for which the dimension will be define
c            (if 0, all variables will be assigned this dimension)
c     idim = dimension number (as you wish to assign it to your array
c             not necessarily as it was stored)
c            if 0, reset all dimensions for ivar to 'not defined'
c            if -1, -2, -3, or -4, set the appropriate dimension to
c               not defined (dimension 1 if -1, dimension 2 if -2, etc.)
c     name = dimension name (a character string which will be truncated
c              here to 16 characters if necessary) 
c     weight = a string indicating the type of weighting associated with
c                this dimension:
c              'width' = weight by the width of the grid cell
c              'unit' = weight by 1
c              'equal' = weight by 1/ngrid, where ngrid is the number of
c                         grid points in this dimension
c              'gaussian' = weight by gaussian weights
c              'cosine' = weight by cosine of latitude
c              'month' = weight by the number of days in the month
c              'leapyr' =weight by the number of days in the months of a
c                          leap year.
c              a standard AMIP model name (e.g., 'csu') = weight by 
c                       proper amount for model indicated
c     indx1 = index of first plane to be retrieved
c     indx2 = index of last plane to be retrieved
c
c     ------
c     Date: 9/17/94
c     ----
c
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine defdimi(ivar, idim, name, weight, indx1, indx2) 

c ???  check the following include statements
      implicit none
#include "usersget"
#include "cvartbl"
#include "cdimtbl"
#include "cscalars"

      integer ivar, idim, indx1, indx2
      logical caseindp
      character*(*) name, weight

      integer i, m, m1, m2, univunit, i1, i2, ii     
      character cwtl*8, domaincl*8
                 
      if (ivar .gt. maxvars) then
 
        print*, 'You are attempting to assign a dimension to'
        print*, 'a nonexistent variable, which you identify'
        print*, 'by the integer ', ivar, ' when calling subroutine'
        print*, 'defdimi.  Note that no more than, ', maxvars
        print*, 'may defined unless you redefine maxvars in'
        print*, 'common usersget' 
        call errcheck('subroutine defdimi', ' ',1)
        return
 
      elseif (ivar .lt. 0) then

        print*, 'You have improperly passed ivar = ', ivar
        print*, 'as the argument that identifies the variable'
        print*, 'for which you want to define a dimension in '
        print*, 'subroutine defdimi.  This integer should be .ge. 0'  
        call errcheck('subroutine defdimi', ' ',1)
        return

      endif

      if ((idim .gt. 4) .or. (idim .lt. -4)) then
        print*, 'You have improperly passed idim = ', idim
        print*, 'as the argument that identifies the dimension'
        print*, 'you want defined in subroutine defdimi.  This integer'
        print*, 'should be in the following range: -5 < idim < 5'
        call errcheck('subroutine defdimi', ' ',1)
        return
      endif   

      if (ivar .eq. 0) then
        m1 = 1
        m2 = maxvars
      else
        m1 = ivar
        m2 = ivar
      endif   

      do 80 m=m1,m2
        do 70 i=1,4
          iptdim(i,m) = 0
          iptdimrg(i,m) = 0
   70   continue
   80 continue
                     
      if (idim .gt. 1) then

        do 100 m=m1,m2              

          do 90 i=1,idim-1
          
          if (caseindp(name, vdimname(i,m))) then
            print*, 'You are attempting to name dimension ', idim
            print*, 'to be the same as dimension ', i
            print*, 'Make sure you are defining your dimensions'
            print*, 'in ascending order.  (e.g., define dimension 3'
            print*, 'after defining dimension 1 and 2 and before '
            print*, 'defining dimension 4)'
            print*, ' '
            if (ivar .gt. 0) then
              print*, 'You were defining the dimension improperly for '
              print*, 'variable ', ivar 
              print*, 'which is currently defined as ', varnam(ivar)
              print*, 'in file: ', filnam(ivar)  
            else
              print*, 'It is possible that this same error will be'
              print*, 'encountered for other variables too.'
            endif 
            call errcheck('subroutine defdimi', ' ', 1)
            return
          endif

   90     continue

  100   continue

      endif

      if (idim .le. 0) then
        if (idim .eq. 0) then
          i1 = 1
          i2 = 4
        else
          i1 = -idim
          i2 = -idim
        endif
        do 140 m=m1,m2
          do 130 i=i1,i2
            vdimname(i,m) = 'not defined'
            vdimwt(i,m) = 'not def.'
            vdimdom(i,m) = 'not def.'
            vdimbdry(1,i,m) = 0.   
            vdimbdry(2,i,m) = 0.   
            ivdimpln(1,i,m) = 0.   
            ivdimpln(2,i,m) = 0.   
            vdimcycl(i,m) = 0.
  130     continue
  140   continue
        return
      endif

      if (len(weight) .lt. 3) then
        print*, weight, ' not recognized as a valid weighting option'
        print*, 'for dimension ', idim, ' named ', name
        print*, 'You were defining the dimension improperly for '
        if (ivar .gt. 0) then    
          print*, 'variable ', ivar 
          print*, 'which is currently defined as ', varnam(ivar)
          print*, 'in file: ', filnam(ivar)  
        else
          print*, 'all variables.'
        endif  
        call errcheck('subroutine defdimi', ' ', 1)
        return
      endif                     

      if (caseindp(weight, 'wid')) then
        cwtl = 'width'
      elseif (caseindp(weight, 'uni')) then    
        cwtl = 'unit'
      elseif (caseindp(weight, 'equ')) then    
        cwtl = 'equal'
      elseif (caseindp(weight, 'gau')) then    
        cwtl = 'gaussian'
      elseif (caseindp(weight, 'cos')) then    
        cwtl = 'cosine'
      elseif (caseindp(weight, 'mon')) then    
        cwtl = 'month'  
      elseif (caseindp(weight, 'lea')) then    
        cwtl = 'leapyr'
      elseif  (caseindp(weight, 'csu') .or. caseindp(weight, 'ucl')
     +  .or. caseindp(weight, 'lmd') .or. caseindp(weight, 'lmc')) then
        if (univunit(name) .eq. intlong) then
          cwtl = 'width'
        elseif (univunit(name) .eq. intlat) then
          cwtl = weight(1:3)
        else
          print*, weight, ' not recognized as a valid weighting option'
          print*, 'for dimension ', idim, ' named ', name
          print*, 'you may only specify this weighting option for'
          print*, 'latitude or longitude dimensions.'
          print*, 'You were defining the dimension improperly for '
          if (ivar .gt. 0) then    
            print*, 'variable ', ivar 
            print*, 'which is currently defined as ', varnam(ivar)
            print*, 'in file: ', filnam(ivar)  
          else
            print*, 'all variables.'
          endif  
          call errcheck('subroutine defdimi', ' ', 1)
          return
        endif
      elseif (caseindp(weight, 'bmr') .or. caseindp(weight, 'ccc') 
     +   .or. caseindp(weight, 'cnr') .or. caseindp(weight, 'col') 
     +   .or. caseindp(weight, 'csi') .or. caseindp(weight, 'der') 
     +   .or. caseindp(weight, 'ecm') .or. caseindp(weight, 'gfd') 
     +   .or. caseindp(weight, 'mgo') .or. caseindp(weight, 'mpi')
     +   .or. caseindp(weight, 'nca') .or. caseindp(weight, 'nmc') 
     +   .or. caseindp(weight, 'rpn') .or. caseindp(weight, 'sun') 
     +   .or. caseindp(weight, 'uga') .or. caseindp(weight, 'sng') 
     +   .or. caseindp(weight, 'gen') .or. caseindp(weight, 'cer')
     +   .or. caseindp(weight, 'ccm') .or. caseindp(weight, 'ccs')
     +   .or. caseindp(weight, 'ech') .or. caseindp(weight, 'nce')) then
        if (univunit(name) .eq. intlong) then
          cwtl = 'width'
        elseif (univunit(name) .eq. intlat) then
          cwtl = 'gaussian'
        else
          print*, weight, ' not recognized as a valid weighting option'
          print*, 'for dimension ', idim, ' named ', name
          print*, 'you may only specify this weighting option for'
          print*, 'latitude or longitude dimensions.'
          print*, 'You were defining the dimension improperly for '
          if (ivar .gt. 0) then    
            print*, 'variable ', ivar 
            print*, 'which is currently defined as ', varnam(ivar)
            print*, 'in file: ', filnam(ivar)  
          else
            print*, 'all variables.'
          endif  
          call errcheck('subroutine defdimi', ' ', 1)
          return
        endif
      elseif (caseindp(weight, 'dnm') .or. caseindp(weight, 'gis') 
     +   .or. caseindp(weight, 'gla') .or. caseindp(weight, 'gsf') 
     +   .or. caseindp(weight, 'iap') .or. caseindp(weight, 'mri') 
     +   .or. caseindp(weight, 'uiu') 
     +   .or. caseindp(weight, 'jma') .or. caseindp(weight, 'nrl') 
     +   .or. caseindp(weight, 'ukm') .or. caseindp(weight, 'yon')) 
     +      then
        if (univunit(name) .eq. intlong) then
          cwtl = 'width'
        elseif (univunit(name) .eq. intlat) then
          cwtl = 'cosine'
        else
          print*, weight, ' not recognized as a valid weighting option'
          print*, 'for dimension ', idim, ' named ', name
          print*, 'you may only specify this weighting option for'
          print*, 'latitude or longitude dimensions.'
          print*, 'You were defining the dimension improperly for '
          if (ivar .gt. 0) then    
            print*, 'variable ', ivar 
            print*, 'which is currently defined as ', varnam(ivar)
            print*, 'in file: ', filnam(ivar)  
          else
            print*, 'all variables.'
          endif  
          call errcheck('subroutine defdimi', ' ', 1)
          return
        endif
      else
        print*, weight, ' not recognized as a valid weighting option'
        print*, 'for dimension ', idim, ' named ', name
        print*, 'You were defining the dimension improperly for '
        if (ivar .gt. 0) then    
          print*, 'variable ', ivar 
          print*, 'which is currently defined as ', varnam(ivar)
          print*, 'in file: ', filnam(ivar)  
        else
          print*, 'all variables.'
        endif  
        call errcheck('subroutine defdimi', ' ', 1)
        return
      endif 

      domaincl = 'byindex'   

      if (indx1*indx2 .le. 0) then
        print*, 'You have improperly passed ', indx1, ' and ', indx2
        print*, 'as the indices of the planes you are requesting by '
        print*, 'calling subroutine defdimi.  These indices must both'
        print*, 'be positive or must both be negative.' 
        call errcheck('subroutine defdimi', ' ',1)
        return
      endif

      ii = min0(16, len(name))                           
      do 200 m=m1,m2
                                         
          vdimname(idim,m) = ' '
          vdimname(idim,m) = name(1:ii)
          vdimwt(idim,m) = cwtl
          vdimdom(idim,m) = domaincl
          ivdimpln(1,idim,m) = indx1
          ivdimpln(2,idim,m) = indx2
          vdimbdry(1,idim,m) = 0.0 
          vdimbdry(2,idim,m) = 0.0 
          vdimcycl(idim,m) = 0.0

          ivarndim(m) = 0
          ivarlen(m) = 0
          iptdim(idim,m) = 0
          iptdimrg(idim,m) = 0

  200 continue

      return
      end       




c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     subroutine defgeog      
c     -------
c
c     Description:
c     -----------
c     this subroutine sets controls for selecting only certain
c       geographical regions.
c                  
c     Usage:
c     ------
c
c      call defgeog(ivar, inorout, imaskvar, select) 
c
c    ivar = defined variable that will be masked
c    inorout = 'in' if field is to be masked before any regridding 
c            = 'out' if field is to be masked after regridding
c    imaskvar = defined variable containing geography
c             = 0 if no masking of inorout
c               (note set select = 'all' if no masking before regridding
c                 AND no masking after regridding)
c    select = string containing the list of regions to be selected  
c                   (i.e., not masked)
c             If 'all', then no geography masking will be done.
c             the available list of regions is given below and should
c                  be separated by commas 
c                 (e.g., select = 'N. Amer., S. Amer.,Greenland')
c                 Note that blanks following commas are optional.

c           
c         REGION                    ACCEPTABLE ABBREVIATIONS
c
c                                     (where '*' represents any
c                                     string of characters, but no
c                                     blanks except in the group: ' & ')
c
c         'land'                     'lan*'
c         'ocean'                    'oce*'
c         'seaice'                   'seai*' or 'sea-i*' or 'sea i*'
c         'North America'            'n* ame*'
c         'South America'            's* ame*'
c         'Africa'                   'afr*'
c         'Greenland'                'gre*'
c         'Europe & Asia'            'eur*asia*'
c         'Antarctica'               'ant*'
c         'Australia'                'aus*'
c         'Indo-Pacific Islands'     'ind* i*'
c         'North Pacific Ocean'      'n* pac*' (which optionally may be 
c                                          followed by 'ocea*' or 'sea')
c         'South Pacific Ocean'      's* pac*' (which optionally may be 
c                                          followed by 'ocea*' or 'sea')
c         'South Pacific 1'          's* pac* 1*' 
c         'South Pacific 2'          's* pac* 2*'
c         'South Pacific 3'          's* pac* 3*'
c         'North Atlantic Ocean'     'n* atl*' (which optionally may be 
c                                          followed by 'ocea*' or 'sea')
c         'South Atlantic Ocean'     's* atl*' (which optionally may be 
c                                          followed by 'ocea*' or 'sea')
c         'South Atlantic 1'         's* atl* 1*'
c         'South Atlantic 2'         's* atl* 2*'
c         'Indian Ocean'             'ind*' (which optionally may be 
c                                          followed by 'ocea*' or 'sea')
c         'Indian 1'                 'ind* 1*'
c         'Indian 2'                 'ind* 2*'
c         'Arctic Ocean'             'arc' (which optionally may be 
c                                          followed by 'ocea*' or 'sea')
c         'Baffin & Hudson_Bays      'baf*
c         'Mediterranean Sea'        'med*' (which optionally may be 
c                                          followed by 'ocea*' or 'sea')
c         'American Lakes'           'ame* l*'
c         'Asian & African Lakes'    'asi* l*'
c
c
c
c     Notes:
c        1. ezget is case insensitive
c        2. blanks (i.e., ' ') is not allowed as part of a 'wild string'
c             except as part of the string ' & '
c        3. '&' may not be substituted for ','  
c            for example, 'Africa & Australia' will select only Africa.
c            to select both 'Africa' and 'Australia', specify:
c             'Africa, Australia'
c        4. 'South Atlantic' = 'South Atlantic 1, South Atlantic 2'
c        5. 'South Pacific' = 
c                   'South Pacific 1, South Pacific 2, South Pacific 3'
c        6. 'Indian Ocean' = 'Indian 1, Indian 2'
c
c
c     ------
c     Date: 9/17/94
c     ----
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine defgeog(ivar, inorout, imaskvar, select)
         
c ??? check following includes
c    need: all, maxvars, ishowin, ishowout, igeogin, igeogout, varnam
      implicit none
#include "usersget"
#include "cscalars"
#include "cvartbl"

      integer ivar, imaskvar, iloc(28), iflg(28)
      logical caseindp
      character*(*) inorout, select
      integer m1, m2, m, j2, j1, j, nreg, kk, k, klen, i2, ia, ib, i3,
     &   mreg, ish, i, jj
      character cloc(28)*5, copy*512, reg*32, c*1, cccc*4, 
     +      abbrev*512 
      data cloc /  'xocey', 'xlany', 'xseay', 'xgrey', 'samey', 'namey',
     +    'xafry', 'xeury', 'xanty', 'xausy', 'xindi', 'npacy', 'spacy',
     +    'xindy', 'xmedy', 'xasil', 'natly', 'xarcy', 'xamel', 'xbafy',
     +    'satly', 'spac1', 'spac2', 'spac3', 'xind1', 'xind2', 'satl1',
     +    'satl2' / 

c    special cases:  eurasia, baffin, indo-pacific islands

      iloc(1) = 2**30
      iloc(2) = 2**29
      iloc(3) = 2**28
      iloc(4) = 2**8
      iloc(5) = 2**12
      iloc(6) = 2**13
      iloc(7) = 2**11
      iloc(8) = 2**14
      iloc(9) = 2**9 
      iloc(10) = 2**10
      iloc(11) = 2**7
      iloc(12) = 2**21
      iloc(13) = 2**24 + 2**25 + 2**26
      iloc(14) = 2**18 + 2**19
      iloc(15) = 2**16
      iloc(16) = 2**6
      iloc(17) = 2**20
      iloc(18) = 2**17
      iloc(19) = 2**5
      iloc(20) = 2**15
      iloc(21) = 2**23 + 2**22
      iloc(22) = 2**26
      iloc(23) = 2**25
      iloc(24) = 2**24
      iloc(25) = 2**19
      iloc(26) = 2**18
      iloc(27) = 2**23
      iloc(28) = 2**22


      if (ivar .gt. maxvars) then
 
        print*, 'You are attempting to define a geography mask for'
        print*, 'a nonexistent variable, which you identify'
        print*, 'by the integer ', ivar, ' when calling subroutine'
        print*, 'defgeog.  Note that the value of of this integer'
        print*, 'parameter must be no larger than ', maxvars 
        print*, 'unless you redefine maxvars in common usersget'
        call errcheck('subroutine defgeog', ' ',1)
        return
 
      elseif (ivar .lt. 0) then

        print*, 'You have improperly passed ivar = ', ivar
        print*, 'as the argument that identifies the variable'
        print*, 'for which you want to define a geography mask. '
        print*, 'This integer should be .ge. 0'  
        call errcheck('subroutine defgeog', ' ',1)
        return                    

      endif                              

      if (imaskvar .gt. maxvars) then
 
        print*, 'You are attempting to define a geography mask'
        print*, 'for the variable, ', ivar, ' but you have passed'
        print*, 'imaskvar = ', imaskvar, ' This is the argument that'
        print*, 'identifies the mask and is too large.  Note that'
        print*, 'imaskvar (an integer) must be no larger than, ',
     +              maxvars 
        print*, 'unless you redefine maxvars in common usersget'
        call errcheck('subroutine defgeog', ' ',1)
        return
 
      elseif (imaskvar .lt. 0) then
                                                      
        print*, 'You are attempting to define a geography mask'
        print*, 'for the variable, ', ivar, ' but you have passed'
        print*, 'imaskvar = ', imaskvar 
        print*, ' This integer must be .gt. 0'
        call errcheck('subroutine defgeog', ' ',1)
        return

      endif
        
      if (ivar .eq. 0) then
        m1 = 1
        m2 = maxvars 
      else
        m1 = ivar
        m2 = ivar
      endif  
                         
      if (caseindp(inorout, 'in')) then
        do 100 m=m1,m2
            igeogin(m) = imaskvar
  100   continue
      elseif (caseindp(inorout, 'out')) then
        do 200 m=m1,m2
            igeogout(m) = imaskvar 
  200   continue
      else
        print*, 'Error ... You have passed the following inappropriate'
        print*, 'string to defgeog:'
        print*, inorout
        print*, 'This string should either be "in" or "out".'
        go to 2000   
      endif

      if (caseindp(select, 'all') .or. (imaskvar .eq. 0)) then
        if (caseindp(inorout, 'in')) then
          do 80 m=m1,m2
            ishowin(m) = all
   80     continue
        elseif (caseindp(inorout, 'out')) then
          do 90 m=m1,m2
            ishowout(m) = all
   90     continue
        endif
        return
      endif

      if (imaskvar .eq. 0) return

c   consider string of the form: "  n Pacif , S. Pac ocean, Ind. Ocean"

c    eliminate blanks immediately preceeding or following ","
c    or immediately preceeding or following "+"
c    or immediately preceeding or following "&" and also
c    reduce multiple consecutive blanks to a single blank 

      j2 = len(select)

c    find first non-blank character       
      j1 = 1
      do 400 j=1,j2
        if (select(j:j) .ne. ' ') go to 450
        j1 = j1 + 1
  400 continue
  450 continue

      if (j1 .ge. j2) then
        if (j1 .gt. j2) then
          print*, 'Error ... You have passed a completely blank string'
          print*, 'to defgeog for selecting geographical regions.' 
          go to 2000
        else
          print*, 'Error ... You have passed a string with only 1 '
          print*, 'non-blank character (the last character)' 
          print*, 'for selecting geographical regions.' 
          go to 2000
        endif
      endif
      
c   copy select but omitting leading blanks, consecutive blanks, and
c      blanks next to commas.  Also put a blank at beginning of string 
c      and a comma at end (if not already present). 
      m2 = 1
      copy = ' ' 
      do 500 j=j1,j2-1
        if (select(j:j) .eq. ' ') then
          c = select(j+1:j+1)
          if ((c .eq. ' ') .or. (c .eq. ',') .or. (copy(m2:m2) .eq. ',')
     +        .or. (c .eq. '&') .or. (copy(m2:m2) .eq. '&')) go to 500
        endif
        m2 = m2+1
        if (m2 .gt. 510) then
          print*, 'Error ...  The string you have passed to defgeog'
          print*, 'has too many characters: '
          print*, select
          go to 2000  
        endif   
        copy(m2:m2) = select(j:j)
  500 continue

      if (select(j2:j2) .ne. ' ') then
        m2 = m2+1
        if ((m2 .gt. 510) .and. (select(j2:j2) .ne. ',')) then
          print*, 'Error ...  The string you have passed to defgeog'
          print*, 'has too many characters: '
          print*, select
          go to 2000  
        endif   
        copy(m2:m2) = select(j2:j2)
      endif
  
      if (copy(m2:m2) .ne. ',') then
        m2 = m2+1
        copy(m2:m2) = ','                         
      endif

      nreg = 0
c      skip leading blank in string: copy
      kk=1
         
c   extract substring (up to  ",")
      m = 1
c    note:  abbrev will have a leading blank, no matter what because
c                             m=1 
      abbrev = ' '

      m2 = m2 + 1 

  600 continue 
        k = kk+1         
        klen = index(copy(k:m2), ',')     
c         jump to 1000 when no more regions selected  
        if (klen .eq. 0) go to 1000
        kk = klen + k - 1 

        if ((klen .lt. 3) .or. (i2 .eq. 2)) then
          print*, 'Error ...  You passed the following string to ',
     +         'defgeog ' 
          print*, select
          print*, 'which contains the following substring '
          print*, copy(k:m2)
          print*, 'which is unrecognizable'    
          go to 2000  
        endif

        if (klen .gt. 32) then
          print*, 'Error ...  You passed the following string to ',
     +         'defgeog ' 
          print*, select
          print*, 'which contains the following substring '
          print*, copy(k:m2)
          print*, 'which is too long.  Perhaps you have forgotten a ',
     +         'comma.'    
          go to 2000  
        endif
        if (klen .eq. 1) then
          print*, 'Error ...  You passed the following string to ',
     +         'defgeog ' 
          print*, select
          print*, 'which has a comma in the wrong place.'
          go to 2000  
        endif

        nreg = nreg + 1 
        if (nreg .gt. 28) then
          print*, 'Error ...  You passed the following string to ',
     +         'defgeog ' 
          print*, select
         print*, 'in which you select more than the number of regions',
     +        ' available.'
          go to 2000  
        endif

        reg = ' '

        do 650 j=1,klen
          c = copy(k+j-1:k+j-1)
          if (lge(c,'A') .and. lle(c,'Z')) c = char(ichar(c)+32)
          reg(j:j) = c
  650   continue

c   extract first character         
        
        c = reg(1:1)

        if (c .eq. 'n') then
          m = m+1
          abbrev(m:m) = 'n'
        elseif (c .eq. 's') then
          c = reg(2:2)
          if (c .ne. 'e') then
            m = m+1
            abbrev(m:m) = 's'
          else        
            m = m+1
            abbrev(m:m) = 'x'
            i2 = 1
            go to 700             
          endif
        else                                
          m = m+1
          abbrev(m:m) = 'x'
          i2 = 1
          go to 700
        endif

c      find first character following first blank

        i2 = 2 + index(reg(2:klen), ' ')

        if (i2 .gt. klen-3) then
          print*, 'Error ...  You passed the following string to ',
     +         'defgeog ' 
          print*, select
          print*, 'which contains the following substring '
          print*, reg
          print*, 'which is inappropriate.'
          go to 2000  
        endif        

c       extract 3 characters                            
                                
  700   c = reg(i2:i2)
        m = m+1
        abbrev(m:m) = c

        c = reg(i2+1:i2+1)
        if (c .eq. ' ') then 
          print*, 'Error ...  You passed the following string to ',
     +         'defgeog '
          print*, select
          print*, 'which contains the following substring '
          print*, reg
          print*, 'which is inappropriate.'
          go to 2000  
        endif  
        m = m+1
        abbrev(m:m) = c
                 
        c = reg(i2+2:i2+2)
        if (c .eq. ' ') then 
          print*, 'Error ...  You passed the following string to ',
     +         'defgeog '
          print*, select
          print*, 'which contains the following substring '
          print*, reg
          print*, 'which is inappropriate.'
          go to 2000  
        endif  
        m = m+1
        abbrev(m:m) = c

c       check whether sea ice requested

        if (abbrev(m-2:m) .eq. 'sea') then
          if (i2 .gt. klen-4) then
c             'sea' is stand-alone without 'ice'
              print*, 'Error ...  You passed the following string to ',
     +              'defgeog '
               print*, select
               print*, 'which contains the following substring '
               print*, reg
               print*, 'which is inappropriate.'
               go to 2000  
          endif
          c = reg(i2+3:i2+3)
          if ((c .eq. ' ') .or. (c .eq. '-')) then
            if (i2 .gt. klen-5) then
              print*, 'Error ...  You passed the following string to ',
     +              'defgeog '
               print*, select
               print*, 'which contains the following substring '
               print*, reg
               print*, 'which is inappropriate.'
               go to 2000  
            endif
            c = reg(i2+4:i2+4)
            if (c .eq. ' ') i2 = i2+4
          endif 
          if (c .ne. 'i') then
              print*, 'Error ...  You passed the following string to ',
     +              'defgeog '
               print*, select
               print*, 'which contains the following substring '
               print*, reg
               print*, 'which is inappropriate. Do you want sea ice?'
               go to 2000  
          endif
        endif

        if (abbrev(m-2:m) .eq. 'eur') then
c         check whether asia also appears in string
          ia = index(reg, 'asia')
          if (ia .eq. 0) then 
            print*, 'Error ...  You passed the following string to ',
     +            'defgeog '
             print*, select
             print*, 'which contains the following substring '
             print*, reg
             print*, 'which is inappropriate. Do you want euro-asia?'
             go to 2000  
          endif
c          check whether more than 1 blank between Europe and Asia
          ib = index(reg(i2:ia), ' ')
          if (ib .ne. 0) then
            i2 = i2 + ib
            ib = index(reg(i2:ia), ' ')
            if (ib .ne. 0) then 
             print*, 'Error ...  You passed the following string to ',
     +             'defgeog '
              print*, select
              print*, 'which contains the following substring '
              print*, reg
              print*, 'which is inappropriate.'
              print*, 'Do you want europe & asia?'
              go to 2000  
            endif
          endif
          i2 = ia + 4       
        endif

c      find first character following blank

        i3 = i2 + index(reg(i2:klen), ' ')

c       check whether we have reached end of string

        if (i3 .eq. i2) then
          m = m+2
          abbrev(m-1:m) = 'y+'
          go to 600
        endif

c       extract 1 character following next blank
                                
        c = reg(i3:i3)
        if (c .eq. 'i') then
          if (abbrev(m-2:m) .eq. 'ind') then
            m = m+2
            abbrev(m-1:m) = 'i+'
          else
              print*, 'Error ...  You passed the following string to ',
     +              'defgeog '
               print*, select
               print*, 'which contains the following substring '
               print*, reg
               print*, 'which is inappropriate.'
               print*, 'Do you want Indo Pacific Islands?'
               go to 2000  
          endif
        elseif (c .eq. 'l') then
          m = m+2
          abbrev(m-1:m) = c//'+'
          abbrev(m-5:m-5) = 'x'
        elseif ((c .eq. '1') .or. (c .eq. '2') .or. (c .eq. '3')) then 
          m = m+2
          abbrev(m-1:m) = c//'+'
        else
c         extract 3 more characters
          if (i3 .gt. 29) then
            print*, 'Error ...  You passed the following string to ',
     +           'defgeog '
            print*, select
            print*, 'which contains the following substring '
            print*, reg
            print*, 'which is inappropriate.'
            go to 2000  
          endif  
          cccc = c
          c = reg(i3+1:i3+1)
          cccc(2:2) = c               
          c = reg(i3+2:i3+2)
          cccc(3:3) = c               
          c = reg(i3+3:i3+3)
          cccc(4:4) = c               
          if ((cccc .eq. 'sea,') .or. (cccc .eq. 'ocea')) then
            m = m+2
            abbrev(m-1:m) = 'y+'
          else 
            print*, 'Error ...  You passed the following string to ',
     +           'defgeog '
            print*, select
            print*, 'which contains the following substring '
            print*, reg
            print*, 'which is inappropriate.'
            go to 2000  
          endif  

        endif
                                      
c      find first character following blank

        i2 = i3 + index(reg(i3:klen), ' ')
        if (i3 .ne. i2) then
            print*, 'Error ...  You passed the following string to ',
     +           'defgeog '
            print*, select
            print*, 'which contains the following substring '
            print*, reg
            print*, 'which is inappropriate.'
            go to 2000  
        endif 

        go to 600      

 1000 continue                           

      do 1050 j=1,28
        iflg(j) = 0
 1050 continue  
        
      mreg = 0
      ish = 0
      do 1100 j=1,28
        i = index(abbrev, cloc(j))
        if (i .gt. 0) then
          mreg = mreg + 1
          iflg(j) = iflg(j) + 1 
          if (iflg(j) .gt. 1) then
            print*, 'Error ... You are attempting to define a'
            print*, 'geography mask, but have included duplicate' 
            print*, 'regions, which is not allowed.'
            print*, 'You passed the following string to defgeog '
            print*, select
            print*, 'You should eliminate one of the duplicate regions.'
            go to 2000  
          endif   
          ish = ish + iloc(j)
        endif
 1100 continue                
     
      if (mreg .gt. nreg) then
         print*, 'Error ... Something wrong in the following string'
         print*, 'which you passed to defgeog: '
         print*, select
         go to 2000  
      endif

      if (mreg .lt. nreg) then
        print*, 'Error ... You have included an unknown region in'
        print*, 'the string you have passed to defgeog.'
        print*, 'You passed the following string to defgeog: '
        print*, select
        go to 2000  
      endif

      if ((iflg(13) .gt. 0) .and.
     +    ((iflg(22)+iflg(23)+iflg(24)) .gt. 0)) then
        print*, 'Error ... in selecting geographical regions, you'
        print*,'may not simultaneously select "S. Pacific" and also any'
        print*, 'of the following: "S. Pacific 1", "S. Pacific 2", '
        print*, '"S. Pacific 3"'
        print*, 'You passed the following string to defgeog: '
        print*, select
        go to 2000
      endif

      if ((iflg(14) .gt. 0) .and. ((iflg(25)+iflg(26)) .gt. 0)) then
        print*, 'Error ... in selecting geographical regions, you'
        print*,'may not simultaneously select "Indian" and also any'
        print*, 'of the following: "Indian 1", "Indian 2" '
        print*, 'You passed the following string to defgeog: '
        print*, select
        go to 2000
      endif

      if ((iflg(21) .gt. 0) .and. ((iflg(27)+iflg(28)) .gt. 0)) then
        print*, 'Error ... in selecting geographical regions, you'
       print*,'may not simultaneously select "S. Atlantic" and also any'
        print*, 'of the following: "S. Atlantic 1", "S. Atlantic 2" '
        print*, 'You passed the following string to defgeog: '
        print*, select
        go to 2000
      endif

      if ((iflg(1)+iflg(3)) .gt. 0) then
        jj=0
        do 1200 j=12,28
          jj = jj + iflg(j)
 1200   continue
        if (jj .gt. 0) then
c        can''t specify ocean+sea ice and also individual ocean basins
          print*, 'Error ... in selecting geographical regions, you'
          print*,'may not simultaneously select "ocean, sea-ice"' 
          print*, 'and also individual ocean basins. '
          print*, 'You passed the following string to defgeog: '
          print*, select
          go to 2000
        endif
        if (iflg(1)+iflg(3) .ne. 2) then
          jj=0
          do 1300 j=4,11
            jj = jj + iflg(j)
 1300     continue
          if (jj .ne. 0) then
            print*, 'Error ... in selecting geographical regions, you'
            print*,'may not select "ocean" separately from "sea-ice"'
            print*, 'while also selecting individual land regions.'
            print*, 'You passed the following string to defgeog: '
            print*, select
            go to 2000
          endif
        endif
      endif

      if (iflg(2) .gt. 0) then
        jj=0
        do 1400 j=4,11
          jj = jj + iflg(j)
 1400   continue
        if (jj .ne. 0) then
          print*, 'Error ... in selecting geographical regions, you'
          print*,'may not select "land" and also individual land' 
          print*, 'regions. '
          print*, 'You passed the following string to defgeog: '
          print*, select
          go to 2000
        endif
      endif            

      if (caseindp(inorout, 'in')) then
        if (ivar .eq. 0) then
          do 1500 m=1,maxvars
            ishowin(m) = ish         
 1500     continue
        else 
          ishowin(ivar) = ish
        endif
      elseif (caseindp(inorout, 'out')) then
        if (ivar .eq. 0) then
          do 1550 m=1,maxvars
            ishowout(m) = ish         
 1550     continue
        else 
          ishowout(ivar) = ish
        endif
      endif

      return

 2000 print*, 'You passed the string:'
      print*, select
      print*, 'You were defining a geography mask for variable ', 
     +     varnam(ivar)
      print*, 'in file ', filnam(ivar)
      call errcheck('subroutine defgeog', ' ',1) 
      return

      end




c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     subroutine defmisc
c     -------
c
c     Description:
c     -----------
c     this subroutine sets controls for the following miscellaneous
c     parameters:
c
c       lendata
c       ittrunc
c       errcntl
c       amissing
c       omit
c       namelong
c       namelat
c       masktypi, masktypo
c       
c       
c     Usage:
c     ------
c
c      call defmisc(param, type, value)
c
c      param = string that describes parameter being set:
c          'data size' = the size of the data vector being passed
c                         to getvdata
c          'truncation' = the spherical harmonic truncation for the
c                    grid (needed if less than global domain data 
c                     are stored.
c          'error control' = an integer controlling number and 
c                            verbosity of error messages.
c                         .le. 0  absolute value gives number of errors 
c                                    allowed before error exiting.
c                                 no warning messages will be displayed
c                         .eq. 0  no limit on number of errors allowed
c                                 no warnings or error messages will be
c                                  displayed
c                         .gt.  0 value gives number of errors allowed 
c                                     before error exiting.
c                                 warnings and error messages will be 
c                                    displayed
c          'input missing value' = value you expect missing data to
c                         have on drs files you are reading.
c          'output missing value' = value you request missing data to 
c                         be assigned after retrieval from drs.
c          'longitude name' = string with name of longitude coordinate
c          'latitude name' = string with name of latitude coordinate
c          'maximum dimension kept' = integer giving maximum size of 
c                   dimension that will be saved and will be retrievable
c                   through calls to getcoord, getedges, and getdimwt.
c                   Longer dimensions must be assigned "unit" weighting.
c         'longest dimension' = integer giving longest dimension
c                for which weights (other than "unit" weights) can be 
c                generated.  The dimension length here refers to the
c                full length as it appears in the file from which the
c                data will be retrieved.   In some cases for cyclic data
c                the length specified here may need to be even larger.
c         'mask type'  integer identifying the type of geography file
c                          that will be accessed on input and output.
c                     = 0 for default: i.e., geography files containing 
c                           0, 1, 2 or 0., 1., 2. or detailed geography 
c                           i.d., or if the geography variable name is
c                           'lndfrc', % land, or if the geography
c                           variable name is 'sifrc', % seaice.
c                     = 1 for geography files containing land fraction
c                            (expressed as a percent)
c                     = 2 for geography files containing sea ice 
c                            fraction (expressed as a percent).
c         'mask type in'  as above except for input mask only
c         'mask type out' as above except for output mask only
c
c     type = 'integer' or 'real' describing type of variable
c                passed as 3rd dummy variable (i.e., value).
c    
c      value = value to assign the parameter
c
c     ------
c     Date: 9/17/94
c     ----
c
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine defmisc(param, type, value)

      implicit none
#include "usersget"
#include "cscalars"

      integer ival
      real value, val
      logical caseindp
      character*(*) param, type
      equivalence(val, ival)
                   
      if (caseindp(param, 'dat')) then

        val = value
        if (caseindp(type, 'in')) then
          if (ival .lt. 0) then
            print*, 'error ... The third parameter in your call to'
            print*, 'defmisc should be a non-negative integer.'
            print*, 'You have incorrectly passed ', ival, ' as your'
            print*, 'third parameter'
            print*, 'You were attempting to define ', param
            call errcheck('subroutine defmisc', ' ',1) 
            return
          else
            lendata = ival
          endif
        else
          print*, 'error ... The second parameter in your call to'
          print*, 'defmisc should be a string containing, "integer"'
          print*, 'You have incorrectly passed ', type, ' as your'
          print*, 'second parameter'
          print*, 'You were attempting to define ', param
          call errcheck('subroutine defmisc', ' ',1) 
          return
        endif

      elseif (caseindp(param, 'trun')) then

        val = value
        if (caseindp(type, 'in')) then
          ittrunc = ival
        else
          print*, 'error ... The second parameter in your call to'
          print*, 'defmisc should be a string containing, "integer"'
          print*, 'You have incorrectly passed ', type, ' as your'
          print*, 'second parameter'
          print*, 'You were attempting to define ', param
          call errcheck('subroutine defmisc', ' ',1) 
          return
        endif

      elseif (caseindp(param, 'err')) then

        val = value
        if (caseindp(type, 'in')) then
          errcntl = ival
        else
          print*, 'error ... The second parameter in your call to'
          print*, 'defmisc should be a string containing, "integer"'
          print*, 'You have incorrectly passed ', type, ' as your'
          print*, 'second parameter'
          print*, 'You were attempting to define ', param
          call errcheck('subroutine defmisc', ' ',1) 
          return
        endif

      elseif (caseindp(param, 'maximum dimension kept')) then

        val = value
        if (caseindp(type, 'in')) then
          if (ival .ge. lencoord/3) then
            mxdmkeep = lencoord/3
            print*, 'warning ... The parameter "maximum dimension kept"'
            print*, 'must be smaller than 1/3 lencoord, which in this'
            print*, 'version of the EZGET library is ', lencoord
            print*, 'EZGET will therefore set this parameter to ',
     &             mxdmkeep
            call errcheck('subroutine defmisc', ' ',-1) 
          else
            mxdmkeep = ival
          endif
        else
          print*, 'error ... The second parameter in your call to'
          print*, 'defmisc should be a string containing, "integer"'
          print*, 'You have incorrectly passed ', type, ' as your'
          print*, 'second parameter'
          print*, 'You were attempting to define ', param
          call errcheck('subroutine defmisc', ' ',1) 
          return
        endif

      elseif (caseindp(param, 'longest dimension')) then

        val = value
        if (caseindp(type, 'in')) then
          mxdimelm = ival
        else
          print*, 'error ... The second parameter in your call to'
          print*, 'defmisc should be a string containing, "integer"'
          print*, 'You have incorrectly passed ', type, ' as your'
          print*, 'second parameter'
          print*, 'You were attempting to define ', param
          call errcheck('subroutine defmisc', ' ',1) 
          return
        endif


      elseif (caseindp(param, 'input missing value')) then

        if (caseindp(type, 're')) then
          amissing = value
          misstype = 1
        elseif (caseindp(type, 'in')) then
          amissing = value
          misstype = 0
        else
          print*, 'error ... The second parameter in your call to'
          print*, 'defmisc should be a string containing, "real"'
          print*, 'or "integer".'
          print*, 'You have incorrectly passed ', type, ' as your'
          print*, 'second parameter'
          print*, 'You were attempting to define ', param
          call errcheck('subroutine defmisc', ' ',1) 
          return
        endif

      elseif (caseindp(param, 'output missing value')) then

        if (caseindp(type, 're')) then
          omit = value
          misstype = 1
        elseif (caseindp(type, 'in')) then
          omit = value
          misstype = 0
        else
          print*, 'error ... The second parameter in your call to'
          print*, 'defmisc should be a string containing, "real"'
          print*, 'or "integer".'
          print*, 'You have incorrectly passed ', type, ' as your'
          print*, 'second parameter'
          print*, 'You were attempting to define ', param
          call errcheck('subroutine defmisc', ' ',1) 
          return
        endif

      elseif (caseindp(param, 'mask type')) then

        val = value
        if (caseindp(type, 'in')) then

          if ((ival .lt. 0) .or. (ival .gt. 4)) then

            print*, 'error ... The third parameter in your call to'
            print*, 'defmisc should be an integer with one of the '
            print*, 'following values:  0 (default mask type), '
            print*, '1 (mask with % land fraction), or 2 (mask with '
            print*, "% sea ice fraction), or 3 (mask with 0's, 1's," 
            print*, "and 2's or 0.0's, 1.0's, and 2.0's for ocean, "
            print*, 'land, and seaice, respectively.'
            print*, 'You were attempting to define ', param
            print*, ' '
            call errcheck('subroutine defmisc', ' ',1) 
            return

          endif

          if (param .eq. 'mask type') then
            masktypi = ival
            masktypo = ival 
          elseif (caseindp(param, 'mask type in')) then
            masktypi = ival
          elseif (caseindp(param, 'mask type out')) then
            masktypo = ival
          else
           print*, ' '
           print*, 'error ...  You have incorrectly passed the'
           print*, 'following string as the first parameter in your'
           print*, 'call to defmisc:', param
           print*, 'Perhaps you meant to select:'
           print*, '"mask type", mask type in", or "mask type out"?'
           print*, ' '
           call errcheck('subroutine defmisc', ' ',1) 
           return
          endif

        else

          print*, 'error ... The second parameter in your call to'
          print*, 'defmisc should be a string containing, "integer"'
          print*, 'You have incorrectly passed ', type, ' as your'
          print*, 'second parameter'
          print*, 'You were attempting to define ', param
          call errcheck('subroutine defmisc', ' ',1) 
          return

        endif

      elseif (caseindp(param, 'longitude name')) then

        if (caseindp(type, 'string')) then
          print*, 'error: In calling defmisc you planned to define the'
          print*, '"longitude name", but you incorrectly specified '
          print*, '"string" as your SECOND parameter (var-type).'
          print*, 'Contrary to the April 1996 EzGet documentation, '
          print*, 'you should pass the name of the longitude dimension'
          print*, 'as the second parameter. '
          print*, 'The last argument (param-value) is ignored in this'
          print*, 'case.'
          print*, ' '
        else
          namelong = type
        endif
              
      elseif (caseindp(param, 'latitude name')) then 

        if (caseindp(type, 'string')) then
          print*, 'error: In calling defmisc you planned to define the'
          print*, '"latitude name", but you incorrectly specified '
          print*, '"string" as your second parameter (var-type).'
          print*, 'Contrary to the April 1996 EzGet documentation, '
          print*, 'you should pass the name of the latitude dimension'
          print*, 'as the SECOND parameter. '
          print*, 'The last argument (param-value) is ignored in this'
          print*, 'case.'
          print*, ' '
        else
          namelat = type
        endif

      elseif (caseindp(param, 'version')) then

        print*, ' '
        print*, 'Executing EZGET version 1.2 '
        print*, 'Released 11 February 2000'
        print*, ' '

      else
 
        print*, 'error ...  You have incorrectly passed the following'
        print*, 'string as the first parameter in your call to defmisc:'
        print*, param
        print*, 'Only the following strings are acceptible choices:'
        print*, '"data length", "truncation", "error control", '
        print*, '"input missing value", "output missing value",'
        print*, '"maximum dimension kept", "longest dimension",'
        print*, '"mask type", "mask type in", "mask type out",'
        print*, '"longitude name", "latitude name", "version"'
        call errcheck('subroutine defmisc', ' ',1) 
        return

      endif     

      return
      end         




c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     subroutine defregrd      
c     -------
c
c     Description:
c     -----------
c
c    This subroutine allows the user to set parameters controlling
c      regridding. 
c
c
c     Usage:
c     ------
c
c      call defregrd(ivar, targtype, itarg, method, nlat, alat, dellat, 
c    +         nlon, alon, dellon)
c
c    ivar = defined variable which will be regridded
c          if = 0, then assign to all variables
c    targtype = type of target grid:
c             'none' don''t regrid
c             'to'  regrid to grid of defined variable itarg
c                 (if itarg=0, then no regridding will be done)
c             'uniform'  regrid to a user-defined lat-long grid
c             'gaussian' regrid to a gaussian grid
c    itarg = defined variable that defines the target grid.
c                 (if itarg=0, then no regridding will be done)
c             if type = 'to' this parameter defines the target grid
c             if type = 'uniform' or 'gaussian', this parameter is
c                        ignored
c    method = type of regridding
c             'area-weighted'  
c             'bilinear'  (not yet available) 
c    nlat = defines the number of target latitudes
c             if type = 'to' this parameter is ignored     
c             if type = 'gaussian' this defines the number of grid
c                        cells for the number of wave numbers
c                        (e.g., for t42, set nlat=64)
c    alat = southern most grid cell if dellat>0
c               northern most grid cell if dellat<0 
c             if type = 'to' this parameter is ignored     
c             if type = 'gaussian', this parameter is ignored
c    dellat = increment between latitudes (set to negative number to
c             define a grid from north to south.)
c             if type = 'to' this parameter is ignored     
c             if type = 'gaussian', this parameter is ignored
c    nlon = defines the number of target longitudes
c             if type = 'to' this parameter is ignored     
c    alon = western most grid cell if dellon>0
c               eastern most grid cell if dellon<0 
c             if type = 'to' this parameter is ignored     
c    dellon = increment between latitudes (set to negative number to
c             define a grid from east to west.)
c             if type = 'to' this parameter is ignored
c     
c     ------
c     Date: 1/1/95
c     ----
c
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine defregrd(ivar, targtype, itarg, method, nlat, alat, 
     +         dellat, nlon, alon, dellon)

      implicit none
#include "usersget"
#include "cvartbl"

      integer ivar, itarg, nlat, nlon, n1, n2, n3, n4, isize, m1, m2, m
      real alat, dellat, alon, dellon
      character*(*) targtype, method
      logical caseindp
                               

      if (.not. caseindp(method, 'area')) then
        print*, method, ' is not an available method of regridding.'
        print*, 'You should currently specify "area-weighted" as '
        print*, 'your method.'
        call errcheck('subroutine defregrd', ' ',1)
        return
      endif
                                
      if (ivar .gt. maxvars) then
 
        print*,'You are attempting to define regridding specifications'
        print*, 'for a nonexistent variable, which you identify'
        print*, 'by the integer ', ivar, ' when calling subroutine'
        print*, 'defregrd.  Note that the value of of this integer'
        print*, 'parameter must be no larger than ', maxvars 
        print*, 'unless you redefine maxvars in common usersget'
        call errcheck('subroutine defregrd', ' ',1)
        return
 
      elseif (ivar .lt. 0) then

        print*, 'You have improperly passed ivar = ', ivar
        print*, 'as the argument that identifies the variable'
        print*, 'for which you want to define a geography mask. '
        print*, 'This integer should be .ge. 0'  
        call errcheck('subroutine defregrd', ' ',1)
        return                    

      endif                              

      if (ivar .eq. 0) then
        m1 = 1
        m2 = maxvars 
      else
        m1 = ivar
        m2 = ivar
      endif  

      if (caseindp(targtype, 'none')) then
         do 80 m=m1,m2
            iregrid(m) = 0
   80    continue
         return
      endif
         
      if (targtype(1:2) .eq.  'to') then
 
        if (itarg .gt. maxvars) then
          print*, 'You are attempting to regrid variable ', ivar
          print*, 'to a grid defined by a nonexistent variable.'
          print*, 'the nonexistent variable is identified by'
          print*, 'index ', itarg, ' which is larger than maxvars.'
          print*, 'In calling defregrd itarg should be less than ',
     +          maxvars     
          print*, 'unless you redefine maxvars in common usersget.'
          call errcheck('subroutine defregrd', ' ',1)
          return                          

        elseif (itarg .lt. 0) then
                                                      
          print*, 'You are attempting to regrid'
          print*, 'variable, ', ivar, ' but you have passed'
          print*, 'itarg = ', itarg 
          print*, ' This integer must be .ge. 0'
          call errcheck('subroutine defgeog', ' ',1)
          return

        endif                             
                  
        do 100 m=m1,m2
            iregrid(m) = itarg
  100   continue

        if (itarg .ne. 0) call shape(itarg, n1, n2, n3, n4, isize)

      elseif (caseindp(targtype, 'gau')) then

        if (dellon .eq. 0.0) then
          print*, 'Error in specifying regridding options.'
          print*, 'The longitude increment should be non-zero'
          call errcheck('subroutine defregrd', ' ', 1)
          return
        endif
         

        do 200 m=m1,m2
            iregrid(m) = -1
            nlatrgdv(m) = nlat
            dlatrgdv(m) = 0.0
            alatrgdv(m) = 0.0
            nlonrgdv(m) = nlon 
            dlonrgdv(m) = dellon
            alonrgdv(m) = alon     
  200   continue

      elseif (caseindp(targtype, 'uni')) then

        if (dellon .eq. 0.0) then
          print*, 'Error in specifying regridding options.'
          print*, 'The longitude increment should be non-zero'
          call errcheck('subroutine defregrd', ' ', 1)
          return
        endif
        if (dellat .eq. 0.0) then
          print*, 'Error in specifying regridding options.'
          print*, 'The latitude increment should be non-zero'
          call errcheck('subroutine defregrd', ' ', 1)
          return
        endif

        do 300 m=m1,m2
            iregrid(m) = -1   
            nlatrgdv(m) = nlat
            dlatrgdv(m) = dellat
            alatrgdv(m) = alat
            nlonrgdv(m) = nlon
            dlonrgdv(m) = dellon
            alonrgdv(m) = alon     
  300   continue
         
      else

        print*, 'Error in specifying regridding options'
        print*, targtype, ' is not a recognized option'
        call errcheck('subroutine defregrd', ' ', 1)
        return

      endif
            
      return
      end        


      

c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     subroutine defvar
c     -------
c
c     Description:
c     -----------
c     This subroutine defines a variable that will subsequently be
c        referred to by the integer, i.  
c
c     Usage:
c     ------
c      call defvar(integer-i.d., variable-name, file-name) 
c
c     ------
c     Date: 1/19/92
c     ----
c
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine defvar(ii, varname, filename)

      implicit none

c      common variables

#include "usersget"
#include "cscalars"
#include "cvartbl"

      character*(*) varname, filename
      integer i, j, ii, opendrs, idicun, ierr
      logical caseindp                                       

         
      i = iabs(ii)

      if (i .eq. 0) then
             print*, 'You have improperly passed ivar = ', i
             print*, 'as the argument that identifies the variable'
             print*, 'you are defining through a call to defvar or'
             print*, 'defvarex.  This integer should not equal 0'  
             call errcheck('subroutine defvar', ' ',1)
             return
      endif

      if (i .gt. maxvars) then
             
             print*, 'You have improperly passed ivar = ', i
             print*, 'as the argument that identifies the variable'
             print*, 'you are defining through a call to defvar or'
             print*, 'defvarex.  The absolute value of this integer'
             print*, 'should be no greater than ', maxvars   
             print*, 'unless you redefine maxvars in common usersget'
             call errcheck('subroutine defvar', ' ',1)
             return
                  
       endif     

         if (len(varname) .gt. 64) then
           varnam(i) = varname(1:64)
           if (errcntl .gt. 0) then
             print*, 'Warning: You are attempting to define a variable'
             print*, 'named ', varname, ' which is longer than 64 '
             print*, 'characters.  This name will be truncated to' 
             print*, varnam(i)
           endif
           call errcheck('subroutine defvar', ' ',-1)
         else
           varnam(i) = varname
         endif              
         
         if (filename .eq. ' ') then
            do 100 j=1,maxvars
              if (.not. (caseindp(filnam(maxvars+1-j), 'not def'))) then
                    filnam(i) = filnam(maxvars+1-j)
                    go to 200
              endif
  100       continue

            print*,'Error .... You must define a filename'
            print*,'in first call to defvar or defvarex'
            call errcheck('subroutine defvar', ' ',1)
            return

  200       continue
      
         else
           if (len(filename) .gt. 120) then
             filnam(i) = filename(1:120)
             if (errcntl .gt. 0) then
               print*, 'Warning: You are attempting to define the '
               print*, 'following directory path and file-name: '
               print*,  filename
               print*, 'which is longer than 120 characters. '
               print*, 'This string will be truncated to: '
               print*, filnam(i) 
             endif
             call errcheck('subroutine defvar', ' ',-1)
           else                            
             filnam(i) = filename
           endif
         endif

      vartitle(i) = 'not defined'
      varunits(i) = 'not defined'
      vartype(i) = 'not def.'
      varsourc(i) = 'not defined'
      vardate(i) = 'not def.'
      vartime(i) = 'not def.'

      ivarndim(i) = 0
      ivarlen(i) = 0
         
      do 250  j=1,4
         iptdim(j,i) = 0
         iptdimrg(j,i) = 0
  250 continue

      if (ii .gt. 0) then
        ierr = opendrs(i, 1, idicun)
        if (ierr .ne. 0) call errcheck(' ', 'subroutine defvar',1)
      else
        ierr = opendrs(i, 0, idicun)
        if (ierr .eq. 0) then
          ii = i
        else
          ii = ierr
          if ((ii .ne. -1000) .and. (ii .ne. -1001)) ii = -1002
        endif
      endif

      return
      end




c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     subroutine defvarex
c     -------
c
c     Description:
c     -----------
c     This subroutine defines a variable that will subsequently be
c        referred to by the integer, i.  It should be used if the
c        the variable name is not unique in the file.
c
c     Usage:
c     ------
c      call defvar(integer-i.d., variable-name, variable-title,
c            source, file-name) 
c
c     Date: 9/20/94
c     ----
c
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine defvarex(ivar, varname, vtitle, vsource, filename)

      implicit none
#include "usersget"
#include "cvartbl"

      character*(*) varname, vtitle, vsource, filename
      integer ivar

      call defvar(ivar, varname, filename)

      if (len(vtitle) .gt. 80) then
        print*, 'Error  ... '
        print*, 'You are attempting to define the following title:'
        print*, vtitle
        print*, 'which is longer than the allowed 80 characters.' 
        print*, 'You were defining variable ', varname, ' in file'
        print*, filename 
        call errcheck('subroutine defvarex', ' ',1)
        return
      else  
        vartitle(ivar) = vtitle
      endif                        

      if (len(vsource) .gt. 120) then  
        print*, 'Error  ... '
        print*, 'You are attempting to define the following source:'
        print*, vsource
        print*, 'which is longer than the allowed 120 characters.' 
        print*, 'You were defining variable ', varname, ' in file'
        print*, filename 
        call errcheck('subroutine defvarex', ' ',1)
        return
      else  
        varsourc(ivar) = vsource
      endif

      return
      end      


c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     subroutine diminfo
c     -------
c
c     Description:
c     -----------
c     this subroutine returns information about a defined dimension
c
c     Usage:
c     ------
c
c      call diminfo
c
c     ------
c     Date: 9/29/98
c     ----
c
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine diminfo(ivar, idim, param, info)

      implicit none
#include "usersget"
#include "cvartbl"
#include "cdimtbl"

      integer ivar, idim
      character*(*) param, info
      integer i, n, iflag

      character dname*16, dsource*120, dtitle*80, dunits*40
      integer idtype, idlen, ndim, ierr, idicun
      real abeg, aend
cjfp was      integer cluvdb, getedim, opendrs
      integer fcw_cluvdb, getedim, opendrs
      logical caseindp


      if (ivar .le. 0) then
             print*, 'You have improperly passed ivar = ', ivar
             print*, 'as the argument that identifies the variable'
             print*, 'you are requesting information from by calling'
             print*, 'subroutine diminfo.  This integer should be' 
             print*, '.gt. 0'  
             call errcheck('subroutine diminfo', ' ',1)
             return
      endif

      if (ivar .gt. maxvars) then
             
             print*, 'You have improperly passed ivar = ', ivar
             print*, 'as the argument that identifies the variable'
             print*, 'you are requesting information from by calling'
             print*, 'subroutine diminfo.  This integer should be' 
             print*, 'no larger than ', maxvars
             print*, 'unless you redefine maxvars in common usersget'
             call errcheck('subroutine diminfo', ' ',1)
             return
                  
       endif     

      if (idim .le. 0) then
             print*, 'You have improperly passed idim = ', idim
             print*, 'as the argument that identifies the dimension'
             print*, 'you are requesting information from by calling'
             print*, 'subroutine diminfo.  This integer should be' 
             print*, '.gt. 0'  
             call errcheck('subroutine diminfo', ' ',1)
             return
      endif

      if (idim .gt. 4) then
             
             print*, 'You have improperly passed idim = ', idim
             print*, 'as the argument that identifies the dimension'
             print*, 'you are requesting information from by calling'
             print*, 'subroutine diminfo.  This integer should be' 
             print*, 'no larger than ', 4
             call errcheck('subroutine diminfo', ' ',1)
             return
                  
       endif     

      i = iptdimrg(idim,ivar)

      if (i .eq. 0) then
       
c       check that variable and file have been defined

        if (caseindp(varnam(ivar), 'not def') .or. 
     +    caseindp(filnam(ivar), 'not def')) then
         print*,'you must call defvar to define variable name and file'
         print*, 'for the variable identified by index: ', ivar
         call errcheck('subroutine diminfo', ' ',1)
         return
        endif

cjfp was       ierr = cluvdb()
        ierr = fcw_cluvdb()
        if (ierr .ne. 0) then
          call errcheck('DRS function cluvdb', 'subroutine diminfo',1)
          return
        endif

        ierr = opendrs(ivar, 1, idicun)
        if (ierr .ne. 0) then
          call errcheck(' ', 'subroutine diminfo',1)
          return
        endif

c         get info on each dimension
        ndim = ivarndim(ivar)

        do  100 i=1,ndim
          ierr = getedim(ivar,dsource,dname,dtitle,dunits,
     $        idtype, idlen, abeg, aend)
          if (ierr .ne. 0) then
            call errcheck('DRS function getedim','subroutine diminfo',1)
            return
          endif
         
  100   continue

      endif

      n = len(info)
      iflag = 0

      if (caseindp(param, 'unit')) then

        if (n .lt. 40) then
          iflag = 40
        else
          if (i .gt. 0) then
            info = dimunits(i)
          else
            info = dunits
          endif
        endif

      elseif (caseindp(param, 'sourc')) then

        if (n .lt. 120) then
          iflag = 120
        else
          if (i .gt. 0) then
            info = dimsourc(i)
          else
            info = dsource
          endif
        endif

      elseif (caseindp(param, 'titl')) then

        if (n .lt. 80) then
          iflag = 80
        else
          if (i .gt. 0) then
            info = dimtitle(i)
          else
            info = dtitle
          endif
        endif

      elseif (caseindp(param, 'domain')) then

        if (n .lt. 8) then
          iflag = 8
        else
          if (i .gt. 0) then
            info = dimdoman(i)
          else
            info = 'unknown'
          endif
        endif

      elseif (caseindp(param, 'weight')) then

        if (n .lt. 8) then
          iflag = 8
        else
          if (i .gt. 0) then
            info = dimwtflg(i)
          else
            info = 'unknown'
          endif
        endif

      elseif (caseindp(param, 'name')) then

        if (n .lt. 16) then
          iflag = 16
        else
          if (i .gt. 0) then
            info = dimnam(i)
          else
            info = dname
          endif
        endif

      else

        print*, 'error ...  You have incorrectly passed the following'
        print*, 'string as the first parameter in your call to diminfo:'
        print*, param
        print*, 'Only the following strings are acceptible choices:'
        print*, '"units", "source", "title", "domain", "weight", "name"'
        call errcheck('subroutine diminfo', ' ',1) 
        return

      endif

      if (iflag .ne. 0) then
        print*, 'Error in retrieving variable information.  The length'
        print*, 'of the character string you pass to diminfo when'
        print*, 'retrieving the ', param, ' should be no shorter than ',
     &           iflag
        print*, 'Your string was only ', n, ' characters long.'
        call errcheck('subroutine diminfo', ' ',1) 
        return
      endif


      return
      end


c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     subroutine domain
c     -------
c
c     Description:
c     -----------
c     This subroutine retrieves the number of dimensions, the names of
c        the dimensions and the domain limits of each dimension as they
c        appear in the drs file.  
c        It may be called after defining a variable.
c
c    
c     Usage:
c     ------
c     call domain(ivar, ndim, dnames, beg, end)
c
c
c     where
c       Input:
c        ivar is the index pointing to the defined variable 
c       Output:
c        ndim is the number of dimensions
c        dnames is a vector of length ndim (or longer) that will contain
c           the names of the dimensions.  Each element of dnames should
c           be large enough to hold 16 characters.
c        beg is a vector of length ndim (or longer) that will contain
c           the values of the first coordinates stored 
c        end is a vector of length ndim (or longer) that will contain
c           the values of the last coordinates stored
c
c                         
c
c     ------
c     Date: 3/15/95
c     ----
c
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine domain(ivar, ndim, dnames, beg, end)

      implicit none
#include "usersget"
#include "cvartbl"

      character dsource*120, dtitle*80, dunits*40
      character*(*) dnames(*)
      integer ivar, idtype, idlen, i, ndim, ierr, idicun
      real beg(*), end(*)
cjfp was      integer cluvdb, getedim, opendrs
      integer fcw_cluvdb, getedim, opendrs
      logical caseindp

       
c       check that variable and file have been defined

      if (caseindp(varnam(ivar), 'not def') .or. 
     +    caseindp(filnam(ivar), 'not def')) then
         print*,'you must call defvar to define variable name and file'
         print*, 'for the variable identified by index: ', ivar
         call errcheck('subroutine domain', ' ',1)
         return
      endif

      ierr = fcw_cluvdb()
      if (ierr .ne. 0) then
        call errcheck('DRS function cluvdb', 'subroutine domain',1)
        return
      endif

      ierr = opendrs(ivar, 1, idicun)
      if (ierr .ne. 0) then
        call errcheck(' ', 'subroutine domain',1)
        return
      endif

c         get info on each dimension
      ndim = ivarndim(ivar)

      do  100 i=1,ndim
        ierr = getedim(i,dsource,dnames(i),dtitle,dunits,
     $        idtype, idlen, beg(i), end(i))
        if (ierr .ne. 0) then
          call errcheck('DRS function getedim','subroutine domain',1)
          return
        endif
         
  100 continue

      return
      end

c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     subroutine domlimit
c     -------
c
c     Description:
c     -----------
c     This subroutine retrieves the domain limits of a dimension and  
c        may be called after retrieving the data or calling shape.
c
c    
c     Usage:
c     ------
c     call domlimit(ivar, dim, beg, end)
c
c
c     where
c       Input:
c        ivar is the index pointing to the defined variable that will
c                be extracted
c        dim is the name of the dimension
c       Output:
c        beg is the beginning edge of the domain
c        end is the ending edge of the domain
c
c                         
c
c     ------
c     Date: 3/1/95
c     ----
c
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine domlimit(ivar, dim, beg, end)

      implicit none
#include "usersget"
#include "cvartbl"
#include "cdimtbl"

      character*(*) dim
      integer ivar, univunit, j, i
      real beg, end
      logical caseindp

      do 50 j=1,4
        i = iptdimrg(j,ivar)
        if (i .ne. 0) then 
          if (caseindp(dimnam(i), dim) .or. 
     &       ((univunit(dim) .gt. 0) .and.
     &        (univunit(dim) .eq. univunit(dimnam(i))))) then
            beg = begbdy(i)
            end = endbdy(i)
            return
          endif
        endif
   50 continue

      print*, 'You are attempting to obtain the domain retrieved '
      print*, 'for the variable, ', varnam(ivar)
      print*, 'stored in file, ', filnam(ivar)
      if (iptdimrg(1,ivar) .eq. 0) then
        print*, 'but you have not yet retrieved it (or called shape)'
        print*, 'so the domain is not available.'
      else
        print*, 'but this dimension does not seem to exist.'
      endif
      call errcheck('function domlimit', ' ', 1)

      return
      end

      



c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     function doregrid 
c     -------
c
c     Description:
c     -----------
c
c           if iregrid > 0  it will point to a variable in the
c                  table that has the target grid
c                   or if > 100 user-specified regridding options
c                 (either given in common or read in or set as default)
c           if iregrid < 0 it will represent a spherical harmonic
c                   truncation (t21=-21, t42 = -42, etc.)
c
c
c     Usage:
c     ------
c     ierr = doregrid(ivar, amaskin, arrayin, amask, array)
c
c     ------
c     Date: 1/19/92
c     ----
c
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer function doregrid(ivar, amaskin, arrayin, amask, array)

          
      implicit none
#include "usersget"
#include "cscalars"
#include "cvartbl"
#include "cdimtbl"
#include "cdomain"
#include "cregddim"
      
      pointer (p1, londx), (p2, lonpt), (p3, latdx), (p4, latpt)
      pointer (p5, wtlon), (p6, wtlat)
      integer ivar
      real amaskin(*), arrayin(*), amask(*), array(*)
      
      character dname*16, dsource*128, dtitle*80, dunits*40
      character cpln(2)*8, cwtloc(2)*8
c      integer idname(2)
      integer idrow(2), idlen(2), idtype(2), 
     +        ilonlat(2), lendim(2), ibegbloc(2), iendbloc(2) 
      integer i, ii, iii, ilong, ilat, k, jj, irloni, irlati,irlono,
     +   irlato, nloni, nlati, nlono, nlato, lonlon, latlat, ibwin, 
     +   ibein, ibsin, ibnin, ibwout, ibeout, ibnout, ibsout, ik,
     &   itim1, itim2, ijk, j1, j2, ierr1, ierr2, ierr3, ierr4, ierr5,
     &   ierr6, ierr, ll, istore, irtim1, irtim2, ntim1, ntim2, 
     &   lentot, maxdimn, idum
   
      pointer (ptv, vecdim), (ptw, twt), (pte, tedge)
      real dum1, dum2, dum3, begbloc(2), cyclel(2), endbloc(2), 
     &     vecdim(*)
      double precision tedge(*), twt(*)
      
      integer mkdimtbl, getdimen, univunit, miallocw, mifree
         
      real    wtlon(*), wtlat(*)
      integer londx(*), lonpt(*), latdx(*), latpt(*)

      doregrid = 0
    
c     look at dimensions for variable ivar, extract longitude and
c       latitude

      ilat = 0
      ilong = 0
      itim1 = 0
      itim2 = 0
      do 100 i = 1,4
         ii = iptdim(i,ivar)
         iii = 0
         if (ii .ne. 0) iii = univunit(dimnam(ii))
         if (iii .eq. intlong)  then
            ilong = i
         elseif (iii .eq. intlat) then
            ilat = i
         elseif (itim1 .eq. 0) then
            itim1 = i
         else
            itim2 = i
         endif

  100 continue
  
      if (ilat .eq. 0) then
        print*, 'Error encountered in attempting to regrid data:'
        print*, varnam(ivar), ' which is stored in '
        print*, filnam(ivar), ' is not a function of latitude.'
        doregrid = 1
        call errcheck('function doregrid', ' ', 1)
        return
      endif 

      if (ilong .eq. 0) then
        print*, 'Error encountered in attempting to regrid data:'
        print*, varnam(ivar), ' which is stored in '
        print*, filnam(ivar), ' is not a function of longitude.'
        doregrid = 1
        call errcheck('function doregrid', ' ', 1)
        return
      endif 

      i = iptdim(ilat,ivar)
  
      if ((dimwtflg(i) .ne. 'gaussian') .and.(dimwtflg(i) .ne. 'cosine')
     +   .and. (dimwtflg(i) .ne. 'csu') .and. (dimwtflg(i) .ne. 'ucl')
     +   .and. (dimwtflg(i) .ne. 'lmd') .and. (dimwtflg(i) .ne. 'lmc')) 
     +     then
         if (errcntl .gt. 0) then
         print*, 'Warning: In attempting to regrid data,'
         print*, 'function doregrid assumes that the weights for the'
         print*, 'input grid represent grid-cell area.  Contrary'
         print*, 'to your current prescription, normally you should'
         print*, 'define the latitude weights to be either'
         print*, 'gaussian or cosine' 
         print*, 'You are regridding a variable named ', varnam(ivar)
         print*, 'which is stored in file ', filnam(ivar)
         endif 
         call errcheck('function doregrid', ' ',-1)
      endif


      if (iregrid(ivar) .ge. 0) then
             
c        look at dimensions for output variable, extract longitude and
c         latitude

         ijk = iregrid(ivar)
         irlono = 0
         irlato = 0
         do 500 i = 1,4
         
            ii = iptdim(i,ijk)
            iii = 0
            if (ii .ne. 0)  iii = univunit(dimnam(ii))
            if (iii .eq. intlong) irlono = ii
            if (iii .eq. intlat)  irlato = ii
                    
 500     continue

         if ((irlono .eq. 0) .or. (irlato .eq. 0)) then
           print*, 'You are attempting to obtain a target grid from a'
           print*, 'variable that is not a function of longitude and '
           print*, 'latitude.'
           doregrid = 1
           call errcheck('function doregrid', ' ',1)
           return
         endif
             
         if ((dimwtflg(irlato) .ne. 'gaussian') .and.
     +        (dimwtflg(irlato) .ne. 'cosine') .and.
     +        (dimwtflg(irlato) .ne. 'csu') .and. 
     +        (dimwtflg(irlato) .ne. 'ucl') .and. 
     +        (dimwtflg(irlato) .ne. 'lmd') .and.  
     +        (dimwtflg(irlato) .ne. 'lmc'))  
     +          then
           if (errcntl .gt. 0) then
           print*, 'Warning: In attempting to regrid data,'
           print*, 'function doregrid assumes that the weights for the'
           print*, 'target grid represent grid-cell area (contrary'
           print*, 'to your current prescription).  Normally you should'
           print*, 'define the latitude weights to be either'
           print*, 'gaussian or cosine'                  
           print*, 'You are regridding a variable named ', varnam(ivar) 
           print*, 'which is stored in file ', filnam(ivar)
           print*, 'The target grid is defined by a variable named ', 
     +                varnam(ijk) 
           print*, 'which is stored in file ', filnam(ijk)
           endif
           call errcheck('function doregrid', ' ',-1)
         endif 
      
      else
         
          ilonlat(1) = -1
          ilonlat(2) = -2
c   ???  check next two assignments
          idrow(1) = iptdim(ilong,ivar)
          idrow(2) = iptdim(ilat,ivar)
c          idname(1) = intlong
c          idname(2) = intlat
c          cpln(1) = 'range'
c          cpln(2) = 'range'
          idtype(1) = 1
c          cwtloc(1) = 'width'

          if (dlatrgdv(ivar) .eq. 0.) then
             idtype(2) = 2
c             cwtloc(2) = 'gaussian'
          else
             idtype(2) = 1
c             cwtloc(2) = 'cosine'
          endif

         nlonrgd = nlonrgdv(ivar)
         nlatrgd = nlatrgdv(ivar)
         dlonrgd = dlonrgdv(ivar)          
         dlatrgd = dlatrgdv(ivar)
         alonrgd = alonrgdv(ivar)
         alatrgd = alatrgdv(ivar)  

c ?/?  check that following works with regrdagain in compall
         begbloc(1) = bdry(1,ilong)
         endbloc(1) = bdry(2,ilong)
         cyclel(1) = cycle(ilong)
         if (cyclel(1) .eq. 0.0) cyclel(1) = 360.0
         if (nlonrgd .gt. 0) then
           idlen(1) = nlonrgd
         else
           j1 = nint((begbloc(1)-alonrgd)/abs(dlonrgd))
           j2 = nint((endbloc(1)-alonrgd)/abs(dlonrgd))
           idlen(1) = iabs(j2-j1) + 1
         endif
         ibegbloc(1) = 0
         iendbloc(1) = 0

c  fuzzy = abs(dlonrgd)*0.0001

c  if (   (amin1(begbloc(1), endbloc(1))+fuzzy  .lt. 
c &    amin1((alonrgd-dlonrgd/2.), (alonrgd+(nlonrgd-0.5)*dlonrgd)))
c &     .or. (amax1(begbloc(1), endbloc(1))-fuzzy  .gt. 
c &    amax1((alonrgd-dlonrgd/2.), (alonrgd+(nlonrgd-0.5)*dlonrgd))))
c &     then

c print*, 'There is a problem with your specifications for'
c print*, 'generating a new grid.  The longitude domain you'
c print*, 'requested is larger than the grid you are creating.'
c print*, 'Either reduce the domain in your call to defdim'
c print*, 'or modify the last 3 arguments in your call to' 
c print*, 'defregrd.' 
c print*, 'You are regridding a variable named ', varnam(ivar) 
c print*, 'which is stored in file ', filnam(ivar)
c call errcheck('function doregrid', ' ',1)

c endif
c ?/?  check that following works with regrdagain in compall
         begbloc(2) = bdry(1,ilat)
         endbloc(2) = bdry(2,ilat)
         cyclel(2) = 0.0
         if (nlatrgd .gt. 0) then
           idlen(2) = nlatrgd
         else
           j1 = nint((begbloc(2)-alatrgd)/abs(dlatrgd))
           j2 = nint((endbloc(2)-alatrgd)/abs(dlatrgd))
           idlen(2) = iabs(j2-j1) + 1
         endif
         ibegbloc(2) = 0
         iendbloc(2) = 0

c         fuzzy = abs(dlatrgd)*0.0001

c         if ((cwtloc(2)(1:6) .eq. 'cosine') .and. 
c     &          (   (amin1(begbloc(2), endbloc(2))+fuzzy  .lt. 
c     &    amin1((alatrgd-dlatrgd/2.), (alatrgd+(nlatrgd-0.5)*dlatrgd)))
c     &     .or. (amax1(begbloc(2), endbloc(2))-fuzzy  .gt. 
c     &   amax1((alatrgd-dlatrgd/2.), (alatrgd+(nlatrgd-0.5)*dlatrgd)))))
c     &     then

c           print*, 'There is a problem with your specifications for'
c           print*, 'generating a new grid.  The latitude domain you'
c           print*, 'requested is larger than the grid you are creating.'
c           print*, 'Either reduce the domain in your call to defdim'
c           print*, 'or modify the last 3 arguments in your call to' 
c           print*, 'defregrd' 
c           print*, 'You are regridding a variable named ', varnam(ivar) 
c           print*, 'which is stored in file ', filnam(ivar)
c           print*, 'The target grid is defined by a variable named ', 
c     +                varnam(ijk) 
c           print*, 'which is stored in file ', filnam(ijk)
c           call errcheck('function doregrid', ' ',1)

c         endif

         do 300 k = 1,2
             i = idrow(k)
          
c              retrieve dimension

            lendim(k) = idlen(k)
            ik = ilonlat(k)
            ierr1 = miallocw(lendim(k), ptv)
            ll = 2*max0(ittrunc,lendim(k))
            ierr2 = miallocw(ll, ptw)
            ll = 2*max0(ittrunc,(lendim(k)+1))
            ierr3 = miallocw(ll, pte)
            istore = 1

            if ((ierr1 .ne. 0) .or. (ierr2 .ne. 0) .or. 
     &          (ierr3 .ne. 0)) then
           print*, 'Error in attempting to allocate memory dynamically.'
           print*, 'You may need more memory.'
           print*, 'Try running on cirrus or retrieving smaller '
           print*, 'portions of the dimension.  Error encountered '
           print*, 'while processing', filnam(ivar), varnam(ivar)
           doregrid = 1
           call errcheck('function miallocw','function doregrid',1)
           return
            endif

            dum1=0.0
            dum2=0.0
            dum3=0.0
            ierr = getdimen(ik, cpln(k), cyclel(k), ibegbloc(k),
     +          iendbloc(k), begbloc(k), endbloc(k), dum1, dum2,
     +          dum3, cwtloc(k), lendim(k), dsource, dname, 
     &          dtitle, dunits, idum, istore, vecdim, tedge, twt)

            if (ierr .ne. 0) then
               print*, 'while processing', filnam(ivar), varnam(ivar)
               doregrid = ierr
               call errcheck('function getdimen','function doregrid',1)
               return
            endif

c           check whether this dimension appears in dimension table
c           if not, add it and also store vecdim, tedge, twt in coordvls
            jj = lendim(k)

            ierr= mkdimtbl(mxdmkeep, istore, dimnam(i), lendim(k),  
     +          idlen(k), dsource, dtitle, dunits, idtype(k),
     &          cpln(k), cyclel(k), begbloc(k), 
     +          endbloc(k), vecdim(1), vecdim(jj), idrow(k), cwtloc(k),
     +          vecdim, tedge, twt)

            if (ierr .ne. 0) then
               doregrid = ierr
               call errcheck(' ', 'function doregrid',0)
               return
            endif

            ierr = mifree(ptv)
            ierr = mifree(ptw)
            ierr = mifree(pte)
            
  300    continue
  
         irlono = idrow(1)
         irlato = idrow(2)

      endif

  
      irloni = iptdim(ilong,ivar)
      irlati = iptdim(ilat,ivar)
      irtim1 = iptdim(itim1,ivar)
      irtim2 = iptdim(itim2,ivar)
      if (irtim1 .gt. 0) then
        ntim1 = idimlen(irtim1)
      else 
        ntim1 = 0
      endif
      if (irtim2 .gt. 0) then
        ntim2 = idimlen(irtim2)
      else 
        ntim2 = 0
      endif
      nloni = idimlen(irloni)
      nlati = idimlen(irlati)
      nlono = idimlen(irlono)
      nlato = idimlen(irlato)
      lonlon = nloni + nlono
      latlat = nlati + nlato
     
c       fill variable tables
      iptdimrg(ilong,ivar) = irlono
      iptdimrg(ilat,ivar) = irlato
      
      lentot = 1
      do 400 i=1,4
          k = iptdimrg(i,ivar)
          if (k .ne. 0) then
               maxdimn=idimlen(k)
          else 
               maxdimn = 0
          endif
          lentot = lentot*max0(maxdimn,1)
  400 continue
  
      lendatan = lentot
      if (lendata .lt. 0) return 
      if ((lendatan .gt. lendata) .and. (lendata .ne. 0)) then
          print*, 'You must dimension your array and mask to be larger.'
          print*, 'You were about to regrid the data when the error was'
          print*, 'detected.'
          print*, 'The minimum dimension of your arrays should be ',
     &      lendatan
          print*, 'Note that if the wrap-around or '
          print*, 'cycling option is in effect for longitude and'
          print*, 'the domain option for the dimension is "range",' 
          print*, 'then one grid cell might be split into 2, so'
          print*, 'you might need to increase that dimension by 1.'
          doregrid = 1
          call errcheck('function doregrid', ' ',1)
          return                                          
      endif
      
      i = iptedge(irloni)
      if (i .lt. 0) then
        print*, 'Dimension vector was not saved but is needed.'
        print*, 
     &    'User should allow longer dimensions to be saved by calling'
        print*, "call defmisc('max dimension kept', 'integer', val)"
        print*, "where val .gt. ", ioriglen(irloni)
        doregrid = 1
        call errcheck('function doregrid', ' ',1)
        return
      endif
 
      if (coordvls(i+1) .gt. coordvls(i)) then
          ibwin = i
          ibein = i+1
      else
          ibwin = i+1
          ibein = i
      endif
      
      i = iptedge(irlati)
      if (i .lt. 0) then
        print*, 'Dimension vector was not saved but is needed.'
        print*, 
     &    'User should allow longer dimensions to be saved by calling'
        print*, "call defmisc('max dimension kept', 'integer', val)"
        print*, "where val .gt. ", ioriglen(irlati)
        doregrid = 1
        call errcheck('function doregrid', ' ',1)
        return
      endif

      if (coordvls(i+1) .gt. coordvls(i)) then
          ibsin = i
          ibnin = i+1
      else
          ibsin = i+1
          ibnin = i
      endif
      
      i = iptedge(irlono)
      
      if (i .lt. 0) then
        print*, 'Dimension vector was not saved but is needed.'
        print*, 
     &    'User should allow longer dimensions to be saved by calling'
        print*, "call defmisc('max dimension kept', 'integer', val)"
        print*, "where val .gt. ", ioriglen(irlono)
        doregrid = 1
        call errcheck('function doregrid', ' ',1)
        return
      endif

      if (coordvls(i+1) .gt. coordvls(i)) then
          ibwout = i
          ibeout = i+1
      else
          ibwout = i+1
          ibeout = i
      endif
      
      i = iptedge(irlato)

      if (i .lt. 0) then
        print*, 'Dimension vector was not saved but is needed.'
        print*, 
     &    'User should allow longer dimensions to be saved by calling'
        print*, "call defmisc('max dimension kept', 'integer', val)"
        print*, "where val .gt. ", ioriglen(irlato)
        doregrid = 1
        call errcheck('function doregrid', ' ',1)
        return
      endif

      if (coordvls(i+1) .gt. coordvls(i)) then
          ibsout = i
          ibnout = i+1
      else
          ibsout = i+1
          ibnout = i
      endif
      
c      pointer (p1, londx), (p2, lonpt), (p3, latdx), (p4, latpt)
c      pointer (p5, wtlon), (p6, wtlat)
     
      ierr1 = miallocw(nlono,p1)
      ierr2 = miallocw(lonlon,p2)
      ierr3 = miallocw(nlato,p3)
      ierr4 = miallocw(latlat,p4)
      ierr5 = miallocw(lonlon,p5)
      ierr6 = miallocw(latlat,p6)

      if ((ierr1 .ne. 0) .or. (ierr2 .ne. 0) .or. (ierr3 .ne. 0) .or.
     +    (ierr4 .ne. 0) .or. (ierr5 .ne. 0) .or. (ierr6 .ne. 0)) then
           print*, 'Error in attempting to allocate memory dynamically.'
           print*, 'You may need more memory.'
           print*, 'Try running on cirrus or retrieving smaller '
           print*, 'portions of the array.  Error encountered '
           print*, 'while processing', filnam(ivar), varnam(ivar)
           doregrid = 1
           call errcheck('function miallocw','function doregrid',1)
           return
      endif

      call maparea(nloni, nlono, nlati, nlato,
     +    coordvls(ibnin), coordvls(ibnout), coordvls(ibsin), 
     +    coordvls(ibsout), coordvls(ibein), coordvls(ibeout),
     +    coordvls(ibwin), coordvls(ibwout), 
     +    londx, lonpt, wtlon, latdx, latpt, wtlat) 

      call rgdarea(ilong, ilat, itim1, itim2,
     +   ntim1, ntim2, nloni, nlono, nlati, nlato, 
     +   londx, lonpt, wtlon, latdx, latpt, wtlat, amissing,
     +   amaskin, arrayin, amask, array)
 
      ierr = mifree(p1)
      ierr = mifree(p2)
      ierr = mifree(p3)
      ierr = mifree(p4)
      ierr = mifree(p5)
      ierr = mifree(p6)
       
     
c       fill variable tables

      ivarlen(ivar) = lentot      

      return
      end
     


c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     subroutine errcheck
c     -------
c
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      subroutine errcheck(sub, from, iadd)

      implicit none
#include "cscalars"

      integer iadd
      character*(*) sub, from
c      character*1 cdum

      if (iadd .lt. 0) then

        warncnt = warncnt + 1
        if (errcntl .le. 0)  return
        print*, ' '
        if (sub(1:1) .ne. ' ')  print*, 'warning detected in ', sub
        if (from(1:1) .ne. ' ') print*, 'which was called by ', from

      else

        errcount = errcount + iadd
        if (errcntl .eq. 0) return
        print*, ' '
        if (sub(1:1) .ne. ' ')  print*, 'error detected in ', sub
        if (from(1:1) .ne. ' ') print*, 'which was called by ', from

        if ((iadd .gt. 0) .and. (errcount .ge. iabs(errcntl))) then
            print*, ' '
          if (errcount .eq. 1) then
            print*, errcount, ' error has been encountered by ezget.'
          else
            print*, 
     +      'So far ',errcount,' errors have been encountered by ezget.'
          endif
          print*, 'Your specified error limit has been exceeded.'
          print*, ' (You may increase the limit by assigning a larger'
          print*, ' integer to errcntl or assign 0 for no limit.)'
          call exit(1)
        endif
        print*, ' '

      endif

      return
      end





c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      function genwts
c     -------
c
c     Description:
c     -----------
c      This function generates a vector of weights and a vector
c      of grid-cell edges based on the value of cwt, obtaining
c      information from the variable table and the information
c      table.  All dimensions are treated with one call.
c
c     Usage:
c     ------
c
c      genwts(dname, istore, fuzzy, cwt, cdomloc, jj, vecdim,
c     +         cycle, jbeg, jend, begb, endb, tedge, twt)
c
c
c     Exit Status:
c     -----------
c     0 if success
c     1 if an error occurred
c
c     ------
c     Date: 1/19/92
c     ----
c
c
c
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



      integer function genwts(dname, istore, fuzzy, cwtloc, cdomloc, jj,
     +     tc, dcycle, jbeg, jend, ibegb, iendb, begb, endb, tedge, twt)


c      common variables

      implicit none
#include "cscalars"
#include "cdomain"
      
      character dname*16, cdomloc*8, cwtloc*8
      logical caseindp                              
      integer univunit, ibegb, iendb
      integer j, jj, jbeg, jend, istore, ij, jjj, jb, je
      real dcycle, begb, endb, fuzzy, tnorm, abc, tc(*)
      double precision pi, eb, bb, ds, dsp, tedge(*), twt(*)
      save daypermo
      real daypermo(12)
      data daypermo/31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,
     +      30.,31./

      genwts = 0
       
      if (istore .eq. 0) then

        if (cdomloc(1:7) .eq. 'assaved') then

          jbeg = 1
          jend = jj

        elseif (cdomloc(1:7) .eq. 'byindex') then

          jbeg = ibegb
          jend = iendb

        elseif (cdomloc(1:7) .eq. 'nearest') then

          jbeg = jj
          jend = 1
          do 40 j=1,jj
            if (tc(jj+1-j) .gt. begb-fuzzy) jbeg = jj+1-j
            if (tc(j) .lt. endb+fuzzy)  jend = j       
   40     continue    
            
          if (jend .lt. jbeg) then
            if (abs(tc(jbeg)-begb) .gt. abs(tc(jend)-endb)) jbeg = jend
            jend = jbeg
          endif

        else

          jbeg = jj
          jend = 1
          do 50 j=1,jj
            if ((tc(j) .ge. begb-fuzzy) .and. (tc(j) .le. endb+fuzzy))
     &          jend = j
            if ((tc(jj+1-j) .ge. begb-fuzzy) .and. 
     &          (tc(jj+1-j) .le. endb-fuzzy))  jbeg = jj+1-j
   50     continue

          if (jend .lt. jbeg) then

            print*, 'All weights for dimension, ', dname, 'are 0.0'
            print*, 'You may have specified range for this dimension'
            print*, 'with boundaries completely outside the domain of'
            print*, 'the data.  You might correct this error by'
            print*, 'specifying "nearest" for the domain option of ' 
            print*, 'this dimension.'
            genwts = 1
            call errcheck('function genwts', ' ',1)
            return

          endif

        endif

        if (jbeg .eq. 1) then
          begb = 1.5*tc(1) - 0.5*tc(2)
        else
          begb = 0.5*(tc(jbeg) + tc(jbeg-1))
        endif

        if (jend .eq. jj) then
          endb = 1.5*tc(jend) - 0.5*tc(jend-1)
        else
          endb = 0.5*(tc(jend) + tc(jend+1))
        endif

        return

      endif

      if ((cwtloc(1:5) .eq. 'equal') .or.
     &    (cwtloc(1:4) .eq. 'unit')) then

c           weight each grid-cell equally except for those outside
c              begbdy and endbdy (set those to 0.)
c              set edges half-way between grid points

            tedge(1) = 1.5*tc(1) - .5*tc(2)
            tedge(jj+1) = 1.5*tc(jj) - .5*tc(jj-1)

            do 100 j=1,jj-1
               tedge(j+1) = (tc(j)+tc(j+1))/2.
               if ((cdomloc(1:7) .eq. 'nearest') .or. 
     +          ((tc(j) .ge. begb-fuzzy) .and. (tc(j) .le.endb+fuzzy)))
     +          then
                  twt(j) = 1.
               else
                  twt(j) = 0.
               endif
  100       continue
            if ((cdomloc(1:7) .eq. 'nearest') .or. 
     +       ((tc(jj) .ge. begb-fuzzy) .and. (tc(jj) .le. endb+fuzzy)))
     +        then
               twt(jj) = 1.
            else
               twt(jj) = 0.
            endif
            tnorm = jj
            if (cwtloc(1:4) .eq. 'unit') tnorm = 1.


      elseif ((cwtloc(1:5) .eq. 'month') .or. 
     +           (cwtloc(1:6) .eq. 'leapyr')) then

c           weight each grid-cell by number of days in month

            if (cwtloc(1:5) .eq. 'month') then
              daypermo(2) = 28.
            else
              daypermo(2) = 29.
            endif    

            tedge(jj+1) = tc(jj) + .5
            do 120 j=1,jj
               tedge(j) = tc(j)-.5
               twt(j) = 0.
             if ((tc(j) .ge. begb-fuzzy) .and. (tc(j) .le. endb+fuzzy))
     +                  then
                 if (tc(j) .gt. 0.1) then
                    ij = (tc(j)-.1)/12.
                    ij = tc(j) - ij*12 +.1
                 else
                    ij = (tc(j)-.1)/12.
                    ij = tc(j) - ij*12 +12.1
                 endif
                 twt(j) = daypermo(ij)
            endif 
  120       continue
  
            if (endb .lt. tedge(1)) endb = tedge(1)
            if (begb .gt. tedge(jj+1)) begb = tedge(jj+1)  
  
            if (cwtloc(1:5) .eq. 'month') then
              tnorm = 365.
            else  
              tnorm = 366.
            endif
                                                       


         elseif (cwtloc(1:5) .eq. 'width') then
c            assume grid-cell boundaries are half-way between
c              grid-coordinates.  Weight by width of cells.

            tedge(1) = 1.5*tc(1) - .5*tc(2)
            tedge(jj+1) = 1.5*tc(jj) - .5*tc(jj-1)

            if ((caseindp(dname, 'pres') .or. caseindp(dname, 'sig'))
     +           .and. (tedge(1) .lt. 0.)) tedge(1) = 0.  
            if (begb .ge. endb) endb = begb + 10.*fuzzy

            if (cdomloc(1:7) .eq. 'nearest') then
               bb = tedge(1) - fuzzy
               eb = tedge(jj+1) + fuzzy
            else
               bb = begb
               eb = endb
            endif
            
            do 200 j=1,jj-1
               tedge(j+1) = (tc(j)+tc(j+1))/2.
               twt(j) = dmin1(eb,tedge(j+1)) - dmax1(bb, tedge(j))
               if (twt(j) .lt. fuzzy) twt(j) = 0.
  200       continue
            twt(jj) = dmin1(eb,tedge(jj+1)) - dmax1(bb, tedge(jj))
            if (twt(jj) .lt. fuzzy) twt(jj) = 0.

            tnorm = tedge(jj+1) - tedge(1)
            if (univunit(dname) .eq. intlong)   tnorm = 360.

         elseif (cwtloc(1:8) .eq. 'gaussian') then
c            use gaussian weighting

c            if ittrunc = 0,
c            assume that number of elements in the drs file for this
c            dimension equals the number of latitudes

c            the argument, tedge, is used as a dummy
c            the coordinate values returned are for the center
c            of the grid-cells, not the edges.  gauaw is called
c            only to retrieve the weights.

            jjj = jj
            if (ittrunc .gt. 0) jjj=ittrunc
            call gauaw(tedge, twt, jjj)

            pi = 4.*datan(1.d0)
            tedge(1) = -1.
            tedge(jjj+1) = 1.
            begb = amin1(amax1(begb, -90.), 90.)
            endb = amin1(amax1(endb, -90.), 90.)
            if (endb .le. begb) then
               endb = amin1((begb + .05), 90.)
               begb = amax1((begb - .05), -90.)
            endif

            do 300 j=1,jjj/2
               tedge(j+1) = twt(j) + tedge(j)
               tedge(jjj+1-j) = -tedge(j+1)
  300       continue

c           note:  tedge(jjj/2 + 1) should be very nearly zero
            if (abs(tedge(jjj/2 + 1)) .gt. .001) then
                print*,'warning calculation of gaussian latitude edges'
                print*, 'and weights may be incorrect .... program'
                print*, 'expects an even number of latitudes ',
     +                    tedge(jjj/2+1)
                call errcheck('function genwts', ' ',-1)
            endif
            
            if (cdomloc(1:7) .eq. 'nearest') then
               bb = -1.
               eb = 1.
            else
               bb = dsin(begb*pi/180.)
               eb = dsin(endb*pi/180.)
            endif

            do 400 j=1,jjj
               twt(j) = dmin1(eb,tedge(j+1)) - dmax1(bb, tedge(j))
               if (twt(j) .lt. 1.e-7) twt(j) = 0.
               tedge(j) = (dasin(tedge(j)))*180./pi
  400       continue
            tedge(jjj+1) = (dasin(tedge(jjj+1)))*180./pi
            tnorm = 2.
            if (ittrunc .gt. 0) then
               do 420 j=1,jjj
                  if ((tedge(j).lt.tc(1)) .and. (tedge(j+1).gt.tc(1)))
     +                        jb = j
                 if ((tedge(j).lt.tc(jj)) .and. (tedge(j+1).gt.tc(jj)))
     +                        je = j
  420          continue
  
               if ((je-jb) .ne. (jj-1)) then 
                  print*, 'error in creating gaussian weights'
                  print*, 'something wrong with the domain'
                  genwts = 1
                  call errcheck('function genwts', ' ',1)
                  return
               endif
               
               do 430 j=1,jj
                  twt(j) = twt(jb-1+j)
                  tedge(j) = tedge(jb-1+j)
  430          continue
               tedge(jj+1) = tedge(jb+jj)
               
            endif
                  

         elseif ((cwtloc(1:6) .eq. 'cosine') .or.
     &           (cwtloc .eq. 'csu') .or.
     +           (cwtloc .eq. 'ucl')) then
c            weight by cosine of latitude, assume boundaries half way
c                between grid-cell centers 

c            assume grid-cell boundaries are half-way between
c              grid-coordinates.  Weight by width of cells.

            tedge(1) = 1.5*tc(1) - .5*tc(2)
            tedge(jj+1) = 1.5*tc(jj) - .5*tc(jj-1)

            if ((cwtloc .eq. 'csu') .or. (cwtloc .eq. 'ucl')) then
               if (tedge(1) .lt. -87.9) tedge(1) = -90.
               if (tedge(jj+1) .gt. 87.9) tedge(jj+1) = 90.
            endif
            

            tedge(1) = dmin1(dmax1(tedge(1),-90.d0), 90.d0)
            tedge(jj+1) = dmin1(dmax1(tedge(jj+1),-90.d0), 90.d0)

            begb = amin1(amax1(begb, -90.), 90.)
            endb = amin1(amax1(endb, -90.), 90.)
            if (endb .le. begb) then
               endb = amin1((begb + 10.*fuzzy), 90.)
               begb = amax1((begb - 10.*fuzzy), -90.)
            endif

            pi = 4.*datan(1.d0)
            
            if (cdomloc(1:7) .eq. 'nearest') then
               bb = dsin(tedge(1)*pi/180.)
               eb = dsin(tedge(jj+1)*pi/180.)
            else
               bb = dsin(begb*pi/180.)
               eb = dsin(endb*pi/180.)
            endif

            dsp = dsin(tedge(1)*pi/180.)

            do 450 j=2,jj
               tedge(j) =(tc(j-1)+tc(j))/2.
               ds = dsp
               dsp = dsin(tedge(j)*pi/180.)
               twt(j-1) = dmin1(eb, dsp) -  dmax1(bb, ds)
               if (twt(j-1) .lt. 1.e-7)  twt(j-1) = 0.
  450       continue
            twt(jj) = dmin1(eb, dsin(tedge(jj+1)*pi/180.)) -
     +                    dmax1(bb, dsp)
            if (twt(jj) .lt. 1.e-7)  twt(jj) = 0.

            tnorm = 2.
            
         elseif ((cwtloc .eq. 'lmd') .or. (cwtloc .eq. 'lmc')) then
c            create special weights for lmd/lmc models

            if (abs(float(jj/2)-(jj/2)) .gt. 0.1) then
              print*, 'Error: You have selected latitude weighting '
              print*, 'appropriate for the lmd model'
              print*, 'but this model must have an even number of '
              print*, 'latitudes, not ', jj, ' latitudes.'
              genwts = 1
              call errcheck('function genwts', ' ',1)
c             return
            endif

            pi = 4.*datan(1.d0)

            tedge(1) = -90.0
            tedge(jj/2+1) = 0.0
            tedge(jj+1) = 90.0

            abc = 2.0d0/jj
            do 460 j=2,jj/2
              tedge(j) = dasin(-1.d0+(j-1)*abc)*180./pi
              tedge(jj+2-j) = dasin(1.d0-(j-1)*abc)*180./pi
  460       continue

            begb = amin1(amax1(begb, -90.), 90.)
            endb = amin1(amax1(endb, -90.), 90.)
            if (endb .le. begb) then
               endb = amin1((begb + 10.*fuzzy), 90.)
               begb = amax1((begb - 10.*fuzzy), -90.)
            endif

            if (cdomloc(1:7) .eq. 'nearest') then
               bb = dsin(tedge(1)*pi/180.)
               eb = dsin(tedge(jj+1)*pi/180.)
            else
               bb = dsin(begb*pi/180.)
               eb = dsin(endb*pi/180.)
            endif

            dsp = dsin(tedge(1)*pi/180.)

            do 470 j=2,jj
               ds = dsp
               dsp = dsin(tedge(j)*pi/180.)
               twt(j-1) = dmin1(eb, dsp) -  dmax1(bb, ds)
               if (twt(j-1) .lt. 1.e-7)  twt(j-1) = 0.
  470       continue
            twt(jj) = dmin1(eb, dsin(tedge(jj+1)*pi/180.)) -
     +                    dmax1(bb, dsp)
            if (twt(jj) .lt. 1.e-7)  twt(jj) = 0.

            tnorm = 2.
            
         else

            print*, 'parameter determining the method of weighting'
            print*, 'has an unrecognized value of ', cwtloc
            genwts = 1
            call errcheck('function genwts', ' ',1)
            return

         endif

         if (cdomloc(1:7) .eq. 'byindex') then
         
            if (iendb .gt. jj) then
              print*, 'Requested data plane, ',iendb ,' not available'
              print*, 'Maximum dimension of ', dname, ' = ', jj
              genwts = 1
              call errcheck('function genwts', ' ',1)
              return
            endif

            if (ibegb .gt. 1) then
              do 500 j=1,ibegb-1
                 twt(j) = 0.
  500         continue
            endif

            if (iendb .lt. jj) then
               do 600 j=iendb+1,jj
                  twt(j) = 0.
  600          continue
            endif

            begb = tedge(ibegb)
            endb = tedge(iendb+1)
         endif

         if (cdomloc(1:7) .eq. 'nearest') then
        
      
            jbeg = jj
            jend = 1
            do 700 j=1,jj
               if (tc(jj+1-j) .gt. begb-fuzzy) jbeg = jj+1-j
               if (tc(j) .lt. endb+fuzzy)  jend = j       
  700       continue    
            
            if (jend .lt. jbeg) then
               if (abs(tc(jbeg)-begb) .gt. abs(tc(jend)-endb)) 
     +              jbeg = jend
               jend = jbeg
            endif

            if ((dcycle .ne. 0.) .and.
     +           (abs(endb-begb-dcycle) .lt. 20.*fuzzy)) 
     +            jend = float(jbeg) - 0.5 + dcycle/(tedge(2)-tedge(1)) 
               
            begb = tedge(jbeg)
            endb = tedge(jend+1)
            
            if (jbeg .gt. 1) then
               do 725 j=1,jbeg-1
                  twt(j) = 0.
  725          continue
            endif
            
            if (jend .lt. jj) then
               do 775 j=jend+1,jj
                  twt(j) = 0.
  775          continue
            endif
            
         endif

c     normalize weights
      jend = 0
      do 800 j=1,jj
         twt(j) = twt(j)/tnorm
         if (twt(j) .ne. 0.) jend=j
         if (twt(jj+1-j) .ne. 0.) jbeg = jj+1-j
  800 continue
  
      if (jend .eq. 0) then 
      
            print*, 'All weights for dimension, ', dname, 'are 0.0'
            print*, 'You may have specified range for this dimension'
            print*, 'with boundaries completely outside the domain of'
            print*, 'the data.  You might correct this error by'
            print*, 'specifying "nearest" for the domain option '
            print*, 'for this dimension.'
            genwts = 1
            call errcheck('function genwts', ' ',1)
            return
            
      endif

      return
      end
      
      

      subroutine getcoord(ivar, idim, coords)
              
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Description:
c     -----------
c     This subroutine retrieves a vector of coordinate values or a 
c        vector of grid-cell edge values for 
c        dimension idim of variable ivar (after the data have been
c        retrieved by a call to getdata or getvdata or after calling
c        shape.
c
c     Usage:
c     ------
c       call getcoord(ivar, idim, coords)
c
c        ivar = variable i.d.
c        idim = dimension number (i.e., 1, 2, 3, or 4, if you want the 
c                 coordinate values associated with the first, second,
c                third or fourth dimension retrieved; 
c                 -1, -2, -3, -4, if you want the coordinate values 
c                associated with the grid cell edges or the first, 
c                second, third or fourth dimension retrieved)
c        coords = vector long enough to receive all the coordinate
c                  values for the requested dimension (should be 1 
c                  greater than the number of grid cells if you
c                  are retrieving grid-cell edges).
c     
c     ------
c     Date: 6/15/94
c     ----
c     
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
#include "usersget"
#include "cvartbl"
#include "cdimtbl"   

      integer ivar, idim
      real coords(*)

      integer i, j, i1, ii

      j = iabs(idim)
      j = iptdimrg(j,ivar)
      
      if (j .eq. 0) then
        print*, 'variable ', ivar, ' doesn''t have a dimension ', idim
        print*, 'or you have cleared the dimension table, so that'
        print*, 'information on the dimensions is no longer available.'
        call errcheck('subroutine getcoord', ' ',1)
        return
      endif
            
      if ((idim .gt. 0) .and. (idim .lt. 5)) then
                  
        i1 = iptcoord(j) - 1
        ii = idimlen(j)

      if (i1 .lt. -1) then
        print*, 'Dimension vector was not saved but is needed.'
        print*, 
     &    'User should allow longer dimensions to be saved by calling'
        print*, "call defmisc('max dimension kept', 'integer', val)"
        print*, "where val .gt. ", ioriglen(j)
        call errcheck('subroutine getcoord', ' ',1)
        return
      endif
      
        do 100 i = 1,ii
          coords(i) = coordvls(i1+i)
  100   continue

      elseif ((idim .lt. 0) .and. (idim .gt. -5)) then

        i1 = iptedge(j) - 1
        ii = idimlen(j) + 1

       if (i1 .lt. 0) then
        print*, 'Dimension vector was not saved but is needed.'
        print*, 
     &    'User should allow longer dimensions to be saved by calling'
        print*, "call defmisc('max dimension kept', 'integer', val)"
        print*, "where val .gt. ", ioriglen(j)
        call errcheck('subroutine getcoord', ' ',1)
        return
      endif
      
        do 200 i = 1,ii
          coords(i) = coordvls(i1+i)
  200   continue

      else

        print*, 'The dimension number passed to getcoord must be in'
        print*, 'the range -5 < i < 0 or 0 < i < 5'
        call errcheck('subroutine getcoord', ' ',1)
        return

      endif

      return
      end


   

      
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     subroutine getdata
c     -------
c
c     Description:
c     -----------
c     This subroutine reads a variable from a DRS file.
c     It checks whether the dimensions expected are correct
c     or at least that the arrays are dimensioned large enough to 
c     accomodate the retrieved data.
c
c     
c     The variable may have up to four dimensions. The only restriction
c     is that the longitude dimension have, approximately, a constant
c     increment between successive values and all dimensions be
c     stored monotonically
c
c     Usage:
c     ------
c     call getdata(ivar, ndim1, ndim2, ndim3, ndim4,
c              len1, len2, len3, len4, amask, array)
c
c
c     where
c        ivar is the index pointing to the defined variable that will
c                be extracted
c        ndim1 is the length of the 1st dimension of amask and array
c               as declared in the calling subroutine
c        ndim2, ndim3, ndim4 are defined as ndim1 but for 2nd, 3rd, and
c                4th dimensions
c
c        len1  on entry:
c                the expected 1st dimension of the data that
c                is being retrieved. (If user specifies len1 to be 0,
c                then he will not be warned when the actual 1st 
c                dimension of the data being retreived is different
c                from what he expects.)
c               on return:
c                the actual 1st dimension of the data that is being
c                    retrieved.
c        len2, len3, len4 are defined as len1, but for 2nd, 3rd, and
c                4th dimensions
c        amask contains weights associated with each element in array.
c        array contains retrieved data.
c                         
c
c     ------
c     Date: 1/14/93
c     ----
c
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine getdata(ivar, ndim1, ndim2, ndim3, ndim4, 
     +                        len1,  len2,  len3,  len4, amask, array)
     
      implicit none
#include "usersget"
#include "cvartbl"  
#include "cdomain"  
#include "cscalars"
     
c      real amask(ndim1,ndim2,ndim3,ndim4),
c     +     array(ndim1,ndim2,ndim3,ndim4)
      integer ivar, ndim1, ndim2, ndim3, ndim4, len1, len2, len3, len4
      real amask(*), array(*)
      integer univunit, iii(4), ll(4), n(4), ln1, ln2, ln3, ln4, 
     &     lendatal, n1, n2, n3, n4, ierr, j, lenarray, kkk, m, l, lll,
     &     k, jjj, i
            
      lendatal = lendata
      n1 = max0(1,ndim1)
      n2 = max0(1,ndim2)
      n3 = max0(1,ndim3)
      n4 = max0(1,ndim4)
      n(1) = ndim1
      n(2) = ndim2
      n(3) = ndim3
      n(4) = ndim4
      ll(1) = len1
      ll(2) = len2
      ll(3) = len3
      ll(4) = len4
      
      lendata = n1*n2*n3*n4
      ierr = 0
      
      do 50 j=1,4
      
         if (ll(j) .gt. n(j)) then 
            ierr = 1
            print*, 'In calling getdata you have specified dimension ',
     +          j, ' to be smaller than'
            print*, 'the length you expect it to be'
            print*, 'You indicate that this dimension has a length of: '
     +          , n(j) 
            print*, 'but you expect the length to be: ', ll(j)
         endif
         
   50 continue                                                

      if (ierr .gt. 0) then
        print*, 'error occurred in subroutine getdata'
        print*, 'while processing', filnam(ivar), varnam(ivar)
        call errcheck(' ', ' ',1)
  
        print*, ' '
        print*, ' '
      endif      
      
      call getvdata(ivar, amask, array)
      
      call lendims(ivar, iii(1), iii(2), iii(3), iii(4), lenarray)
             
      
      ierr = 0
      
      do 100 j=1,4
      
         if (iii(j) .gt. n(j)) then
            ierr = 1
            print*, 'Error: Dimension ', j, ' of the array and mask '
            print*, 'passed to getdata is too small to accomodate'
            print*, 'the data you have requested.'
            print*, 'Either change the domain requested or increase'
            print*, 'this dimension of the arrays to ', iii(j) 
            if ((iii(j) .eq. n(j)+1) .and. (cycle(j) .ne. 0.)) then
              kkk = univunit(namedim(j)) 
              if (kkk .eq. intlong) then 
                print*, 'You have invoked the wrap-around or cycle'
                print*, 'option, which may split 1 longitude'
                print*, 'cell into 2 cells, so you probably need to'
                print*, 'increase your longitude dimension by 1.'
                print*, 'Alternatively, for this dimension specify'
                print*, '"nearest" as the domain option and the range'
                print*, 'will be shifted so as not to split the cell.'
              endif      
            endif
         endif
         
         if ((ll(j) .gt. 0) .and. (iii(j) .ne. ll(j))) then
            
            if (errcntl .gt. 0) then
            print*,'Warning: In calling getdata you expected dimension '
     +           , j
            print*, 'to be ', ll(j), ' but in fact it turns out to be '
     +           , iii(j) 
            print*, 'To correct this problem, set the expected'
            print*, 'length of this dimension to ', iii(j)
            print*, 'or to supress this warning, set expected length'
            print*, 'of this dimension equal to 0.'
            endif
            call errcheck('subroutine getdata', ' ',-1)
         endif

  100 continue

      if (ierr .gt. 0) then
        print*, 'error occurred in subroutine getdata'
        print*, 'while processing', filnam(ivar), varnam(ivar)
        call errcheck(' ', ' ',1)
  
        print*, ' '
        print*, ' '
      
      else  

        ln1 = max0(1,iii(1))
        ln2 = max0(1,iii(2))
        ln3 = max0(1,iii(3))
        ln4 = max0(1,iii(4))
  
        m = lenarray + 1
        do 500 l=ln4,1,-1
           lll = n3*(l-1)
           do 400 k=ln3,1,-1
              kkk = n2*(lll+k-1)
              do 300 j=ln2,1,-1
                 jjj = n1*(kkk+j-1)
                 do 200 i=ln1,1,-1
                    m = m - 1
                    amask(jjj+i) = amask(m)
                    array(jjj+i) = array(m)
  200            continue
  300         continue
  400      continue
  500   continue
  
       if (ln4 .lt. n4) then
        do 640 l=ln4+1,n4
           lll = n3*(l-1)
           do 630 k=1,n3
              kkk = n2*(lll+k-1)
              do 620 j=1,n2
                 jjj = n1*(kkk+j-1)
                 do 610 i=1,n1
                    amask(jjj+i) = 0.0
                    array(jjj+i) = 0.0
  610            continue
  620         continue
  630      continue
  640   continue
       endif

       if (ln3 .lt. n3) then
        do 740 k=ln3+1,n3
           kkk = n1*n2*(k-1)
           do 730 l=1,ln4
              lll = kkk + n1*n2*n3*(l-1)
              do 720 j=1,n2
                 jjj = lll + n1*(j-1)
                 do 710 i=1,n1
                    amask(jjj+i) = 0.0
                    array(jjj+i) = 0.0
  710            continue
  720         continue
  730      continue
  740   continue
       endif
            
       if (ln2 .lt. n2) then   
        do 840 j=ln2+1,n2
           jjj = n1*(j-1)
           do 830 l=1,ln4
              lll = jjj + n1*n2*n3*(l-1)
              do 820 k=1,ln3
                 kkk = lll + n1*n2*(k-1)
                 do 810 i=1,n1
                    amask(kkk+i) = 0.0
                    array(kkk+i) = 0.0
  810            continue
  820         continue
  830      continue
  840   continue
       endif
               
       if (ln1 .lt. n1) then
        do 940 i=ln1+1,n1
           do 930 l=1,ln4
              lll = n3*(l-1)
              do 920 k=1,ln3
                 kkk = n2*(lll+k-1)
                 do 910 j=1,ln2
                    jjj = n1*(kkk+j-1)
                    amask(jjj+i) = 0.0
                    array(jjj+i) = 0.0
  910            continue
  920         continue
  930      continue
  940   continue
       endif
    
      endif
        
      len1 = iii(1)
      len2 = iii(2)
      len3 = iii(3)
      len4 = iii(4)
      
      lendata = lendatal

      return
      end
      

      
    
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     integer function getdimen
c     -------
c
c     Description:
c     -----------
c     This function retrieves information about a single dimension of
c       a single variable and stores the information in vartabl and
c       dimtabl
c
c
c     Usage:
c     ------
c
c     ierr =  getdimen(idim, domainc, cycle, ibegbdy, iendbdy
c    +      begbdy, endbdy, felretr, lelretr, delta, cwtl, lendim, 
c    +      istore, vecdim, tedge, twt)
c
c
c     Exit Status:
c     -----------
c     0 if success
c     1 if an error occurred
c
c     ------
c     Date: 1/19/92
c     ----
c
c     Output:
c           felretr = value of first element to be retrieved for this
c                          dimension
c           lelretr = value of last element to be retrieved for this
c                          dimension
c           lendim = number of elements to be retrieved for this
c                          dimension
c
c
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      integer function getdimen(idim, domainc, cycle, ibegbdy, iendbdy,
     +      begbdy, endbdy, felretr, lelretr, delta, cwtl, lendim, 
     +      dsource, dname, dtitle, dunits, idtype, istore, vecdim, 
     &      tedge, twt)


c      common variables
          
      implicit none
#include "drsdef.h"
#include "fcddrs.h"

#include "usersget"
#include "cscalars"
#include "cregddim"


      integer idim, ibegbdy, iendbdy, lendim, istore, idtype
      character domainc*8, cwtl*8, dname*16, dsource*120, dtitle*80,
     &      dunits*40
      real cycle, begbdy, endbdy, felretr, lelretr, delta
      double precision pi, tedge(*), twt(*)
      real vecdim(*)

cjfp was      integer getcdim
      integer fcw_getcdim
      integer ierr, genwts, jj, nc, jbeg, jend,
     +       ncb, nce, j, n, lbegbdy, lendbdy, j1, j2, ij, i, je,
     &       itemp
      real fuzzy, begb, endb, vect, dx, yrlen
      
      save daypermo
      real daypermo(12)
      data daypermo/31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,
     +      30.,31./
     
      getdimen = 0

      if (idim .gt. 0) then
      
         ierr = fcw_getcdim(idim, dsource, dname, dtitle, dunits, 
     +                     idtype, lendim, vecdim, jj)
         if (ierr .ne. 0) then
            getdimen = 1
            call errcheck('DRS function fcw_getcdim',
     1           'function getdimen',1)
            return
         endif

         if (domainc(1:7) .eq. 'assaved') then

            begbdy = vecdim(1)
            endbdy = vecdim(jj)
  
         endif

         if (domainc(1:7) .eq. 'byindex') then

           if ((ibegbdy .gt. jj) .or. (iendbdy .gt. jj)) then
            print*, 'you have requested a plane from dimension ', dname
            print*, 'that is too large.  This dimension extends from 1'
            print*, 'to ', jj, '  Check your call to defdimi.'
            getdimen = 1
            call errcheck('function getdimen', ' ',1)
            return
           endif

           if (ibegbdy .le. iendbdy) then
             begbdy = vecdim(1)
             endbdy = vecdim(jj)
           else
             begbdy = vecdim(jj)
             endbdy = vecdim(1)
           endif

         endif
        
      elseif (idim .eq. -1) then
c        create longitude grid

         dname='Longitude'
         domainc = 'range'
         cwtl = 'width'
         ibegbdy = 0
         iendbdy = 0
         jj = lendim

         if (nlonrgd .eq. 0) then
c ?/?
c          create longitude grid that spans domain needed
           j1 = nint((begbdy-alonrgd)/abs(dlonrgd))
           j2 = nint((endbdy-alonrgd)/abs(dlonrgd))
           do 22 j=1,jj
             vecdim(j) = alonrgd + (min0(j1,j2)+j-1)*abs(dlonrgd)
   22      continue

         else

           do 25 j=1,jj 
              vecdim(j) = alonrgd + (j-1)*dlonrgd
   25      continue

         endif
   
      elseif (idim .eq. -2) then
c        create latitude grid

         dname='Latitude'
         domainc = 'range'
         ibegbdy = 0
         iendbdy = 0
         pi = 4.*datan(1.d0)
         jj = lendim

         if (dlatrgd .eq. 0.) then

           cwtl = 'gaussian'
c           create gaussian latitudes
           call gauaw(tedge, twt, jj)
           do 65 j=1,jj 
              vecdim(j) = (dasin(tedge(j)))*180./pi
   65      continue

         else

           cwtl = 'cosine'
           if (nlatrgd .eq. 0) then

c            create latitude grid that spans domain needed
             j1 = nint((begbdy-alatrgd)/abs(dlatrgd))
             j2 = nint((endbdy-alatrgd)/abs(dlatrgd))
             do 27 j=1,jj
               vecdim(j) = alatrgd + (min0(j1,j2)+j-1)*abs(dlatrgd)
   27        continue
    
           else

             jj = nlatrgd
             do 75 j=1,jj 
                vecdim(j) = alatrgd + (j-1)*dlatrgd
   75       continue

           endif

         endif
     
      endif
   
      if (begbdy .le. endbdy) then
         begb = begbdy
         endb = endbdy
      else
         begb = endbdy
         endb = begbdy
      endif

c           assure ascending order for dimension and 
c           find maximum (absolute) of vector and define
c               a fuzzy such that fuzzy tests of equivalence
c               can be performed on floating point numbers.


      if (vecdim(1) .gt. vecdim(jj)) then
        je = jj/2
        do 100 j=1,je
            vect = vecdim(j)
            vecdim(j) = vecdim(jj+1-j) 
            vecdim(jj+1-j) = vect
100     continue
      endif

      fuzzy = 0.0
      do 200 j=1,jj
        fuzzy = amax1(fuzzy, abs(vecdim(j)))
  200 continue
      fuzzy = fuzzy* 1.e-5
         
      nc = 0

      if (domainc(1:7) .eq. 'byindex') then

        if (ibegbdy .gt. 0) then 
          lbegbdy = ibegbdy
          lendbdy = iendbdy
        else
          lbegbdy = jj+1+ibegbdy
          lendbdy = jj+1+iendbdy
        endif

        if (iendbdy .ge. ibegbdy) then

          if (endbdy .ge. begbdy) then
            lbegbdy = lbegbdy
            lendbdy = lendbdy
          else
            itemp = lbegbdy
            lbegbdy = jj+1-lendbdy
            lendbdy = jj+1-itemp
          endif

        else

          if (endbdy .ge. begbdy) then
            lbegbdy = jj+1-lbegbdy
            lendbdy = jj+1-lendbdy
          else
            itemp = lbegbdy
            lbegbdy = lendbdy
            lendbdy = itemp
          endif

        endif
      endif

      if ((domainc(1:7) .eq. 'assaved') .or.
     +                 (domainc(1:7) .eq. 'byindex')) then

c        temporarily define begb and endb to be well beyond first
c          and last elements in vecdim
         begb = vecdim(1) - 10.*(vecdim(jj) - vecdim(1)) - 10000.
         endb = vecdim(jj) + 10.*(vecdim(jj) - vecdim(1)) + 10000.

      elseif ((cycle .ne. 0.) .and. (idim .gt. 0)) then
c        possibly change begb and endb to fall in range of vecdim
c ***    assume that grid points are equally spaced

        dx = 0.
        if (jj .ne. 1) then
          dx = (vecdim(jj) - vecdim(1))/(jj-1) 
          if (abs(dx*jj-cycle) .gt. fuzzy) then
         print*, 'cycle prescribed is not consistent with evenly spaced'
           print*, 'grid-cells for dimension ', dname
           print*, 'first and last values are: ', vecdim(1), vecdim(jj)
           print*, 'number of elements is: ', jj
           getdimen = 1
           call errcheck('function getdimen', ' ',1)
           return
          endif
        endif
        ncb = 0

  300   if (begb .lt. vecdim(1)-dx/2.-fuzzy) then
           begb = begb + cycle
           ncb = ncb + 1
           go to 300
        endif

  325   if (begb .gt. vecdim(jj)+dx/2.+fuzzy) then
           begb = begb - cycle
           ncb = ncb - 1
           go to 325
        endif

        nce = 0
  350   if (endb .lt. vecdim(1)-dx/2.-fuzzy) then
           endb = endb + cycle
           nce = nce + 1
           go to 350
        endif

  375   if (endb .gt. vecdim(jj)+dx/2.+fuzzy) then
           endb = endb - cycle
           nce = nce - 1
           go to 375
        endif

        nc = ncb

        if ((ncb .ne. 0) .or. (nce .ne. 0)) wrap = .true.

        if (ncb .lt. nce) then
          print*, 'coding error in function getdimen'
          print*, 'check logic'
          getdimen=1
          call errcheck(' ', ' ',1)
          return
        endif

        if (ncb .gt. nce) then
c          extend vecdim

          if (lendim .lt. jj*(1+ncb-nce)) then
            lendim = (1+ncb-nce)*jj
            getdimen = -1
            return
          endif

          endb = endb+ cycle*(ncb-nce)
          do 390 n=1,ncb-nce
            do 380 j=1,jj
               vecdim(jj*n+j) = vecdim(j) + n*cycle
  380       continue
  390     continue

          jj = (1+ncb-nce)*jj

        endif

      endif

c      generate weights and grid-cell edges for all coordinates
c

      if (domainc(1:7) .eq. 'nearest') then
          if (begb .gt. vecdim(jj)) begb = vecdim(jj)
          if (endb .lt. vecdim(1)) then
            if (cycle .eq. 0.0) then
              endb = vecdim(1)
            else
              if ((vecdim(1)-endb) .lt. (begb+cycle-vecdim(jj))) then
                endb = vecdim(1)
              else
                begb = vecdim(jj)
                endb = begb + dx*(jj-1) 
                nc = nc+1
              endif
           endif
        endif
      endif

      if (jj .eq. 1) then

           if ((domainc(1:7) .eq. 'byindex') .and. 
     +         ((lbegbdy .ne. 1) .or. (lendbdy .ne. 1))) then
              print*, 'Error ... you have requested a domain "byindex"'
              print*, 'for the dimension, ', dname
              print*, 'but only a single plane is stored for this'
              print*, 'dimension.'          
              getdimen = 1
              call errcheck('function getdimen', ' ',1)
              return
            endif

           if ((domainc(1:5) .eq. 'range') .and.
     +         ((begb .gt. vecdim(1)) .or. (endb .lt. vecdim(1)))) then
              print*, 'Error ... no data falls within the range you'  
              print*, 'have requested for dimension, ', dname
              print*, 'range requested: ', begb, endb
              print*, 'Only a single plane exists for this dimension'
              print*, 'and the coordinate value is: ', vecdim(1)
              getdimen = 1
              call errcheck('function getdimen', ' ',1)
              return
           endif
      
           if ((cwtl(1:5) .eq. 'month') .or. (cwtl(1:6) .eq. 'leapyr')) 
     +           then

              if (cwtl(1:6) .eq. 'leapyr') then
                 daypermo(2) = 29.
                 yrlen = 366.
              else  
                 daypermo(2) = 28.
                 yrlen = 365. 
              endif
   
              tedge(1) = vecdim(1) - 0.5
              tedge(2) = vecdim(1) + 0.5
              ij = (vecdim(1)-.99)/12.
              ij = vecdim(1) - ij*12 +.01
              twt(1) = daypermo(ij)/yrlen

           elseif (idim .eq. -1) then
 
              tedge(1) = amax1(begb, vecdim(1) - abs(dlonrgd)/2.)
              tedge(2) = amin1(endb, vecdim(1) + abs(dlonrgd)/2.)
              twt(1) = amin1((endb-begb), abs(dlonrgd))/360.

           elseif (idim .eq. -2) then
       
              tedge(1) =
     &            amin1(amax1(begb, vecdim(1) - abs(dlatrgd)/2.), 90.)
              tedge(2) =
     &            amax1(amin1(endb, vecdim(1) + abs(dlatrgd)/2.), -90.)
           twt(1) = (dsin(tedge(2)*pi/180.) - dsin(tedge(1)*pi/180.))/2.

           else

              tedge(1) = vecdim(1)-fuzzy
              tedge(2) = vecdim(1)+fuzzy
              twt(1) = 1.

           endif
           
           begb = tedge(1)
           endb = tedge(2)
           jbeg = 1
           jend = 1

      else
          ierr = genwts(dname, istore, fuzzy, cwtl, domainc, jj, vecdim,
     +     cycle, jbeg, jend, lbegbdy, lendbdy, begb, endb, tedge, twt)
          if (ierr .ne. 0) then
              getdimen=1
              call errcheck(' ', 'function getdimen',0)
              return
          endif

      endif

      lendim = jend - jbeg + 1

      if (begbdy .le. endbdy) then

         begbdy = begb - nc*cycle
         endbdy = endb - nc*cycle
         delta = (vecdim(jend)-vecdim(jbeg))/(10000.*lendim)
         felretr = vecdim(jbeg) - nc*cycle - delta
         lelretr = vecdim(jend) - nc*cycle + delta

         if (istore .gt. 0) then

           do 500 j=jbeg,jend
             vecdim(j-jbeg+1) = vecdim(j) - nc*cycle
             tedge(j-jbeg+1) = tedge(j) - nc*cycle
             twt(j-jbeg+1) = twt(j)
  500      continue
           tedge(jend-jbeg+2) = tedge(jend+1) - nc*cycle
           if (begbdy .lt. tedge(1)) begbdy = tedge(1)
           if (tedge(1) .lt. begbdy) tedge(1) = begbdy
           if (endbdy .gt. tedge(lendim+1)) endbdy = tedge(lendim+1)
           if (tedge(lendim+1) .gt. endbdy) tedge(lendim+1) = endbdy

         else

           do 550 j=jbeg,jend
             vecdim(j-jbeg+1) = vecdim(j) - nc*cycle
550        continue

         endif

      else

         begbdy = endb - nc*cycle
         endbdy = begb - nc*cycle
         delta = -(vecdim(jend)-vecdim(jbeg))/(10000.*lendim)
         felretr = vecdim(jend) - nc*cycle - delta
         lelretr = vecdim(jbeg) - nc*cycle + delta

         if (istore .gt. 0) then
           i = 0
           do 600 j=jbeg,jend
             i = i+1
             vecdim(i) = vecdim(j) - nc*cycle
             twt(i) = twt(j)
             tedge(i) = tedge(j) - nc*cycle
  600      continue
           tedge(i+1) = tedge(jend+1) - nc*cycle

           je = lendim/2
           do 625 j=1,je

             vect = vecdim(j)
             vecdim(j) = vecdim(lendim+1-j) 
             vecdim(lendim+1-j) = vect

             vect = twt(j)
             twt(j) = twt(lendim+1-j)
             twt(lendim+1-j) = vect

             vect = tedge(j)
             tedge(j) = tedge(lendim+2-j)
             tedge(lendim+2-j) = vect

  625      continue

           if ((lendim+1)/2 .gt. je) then

             je = (lendim+1)/2
             vect = tedge(je)
             tedge(je) = tedge(lendim+2-je)
             tedge(lendim+2-je) = vect

           endif

           if (begbdy .gt. tedge(1)) begbdy = tedge(1)
           if (tedge(1) .gt. begbdy) tedge(1) = begbdy
           if (endbdy .lt. tedge(lendim+1)) endbdy = tedge(lendim+1)
           if (tedge(lendim+1) .lt. endbdy) tedge(lendim+1) = endbdy

         else

           i = 0
           do 650 j=jbeg,jend
             i = i+1
             vecdim(i) = vecdim(j) - nc*cycle
650        continue
           tedge(i+1) = tedge(jend+1) - nc*cycle

           je = lendim/2
           do 675 j=1,je

             vect = vecdim(j)
             vecdim(j) = vecdim(lendim+1-j) 
             vecdim(lendim+1-j) = vect

  675      continue

         endif
      endif

      return
      end






      subroutine getdimwt(ivar, idim, wts)
              
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Description:
c     -----------
c     This subroutine retrieves a vector of weights associated with
c        each of the grid cells of dimension idim of variable ivar
c        (after the data have been retrieved by a call to getdata
c        or getvdata, or after calling shape.
c
c     Usage:
c     ------
c       call getdimwt(ivar, idim, wts)
c
c        ivar = variable i.d.
c        idim = dimension number (i.e., 1, 2, 3, or 4, if you want the 
c                 weights associated with the first, second,
c                third or fourth dimension retrieved; 
c        wts = vector long enough to receive all the weights
c                  for the requested dimension.
c     
c     ------
c     Date: 1/16/96
c     ----
c     
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
#include "usersget"
#include "cvartbl"
#include "cdimtbl"   

      integer ivar, idim
      real wts(*)
      integer i, j, ii, i1

      j = iabs(idim)
      j = iptdimrg(j,ivar)
      
      if (j .eq. 0) then
        print*, 'variable ', ivar, ' doesn''t have a dimension ', idim
        print*, 'or you have cleared the dimension table, so that'
        print*, 'information on the dimensions is no longer available.'
        call errcheck('subroutine getdimwt', ' ',1)
        return
      endif
            
      if ((idim .gt. 0) .and. (idim .lt. 5)) then
                  
        i1 = iptwts(j) - 1
        ii = idimlen(j)
 
       if (i1 .lt. 0) then
        print*, 'Dimension vector was not saved but is needed.'
        print*, 
     &    'User should allow longer dimensions to be saved by calling'
        print*, "call defmisc('max dimension kept', 'integer', val)"
        print*, "where val .gt. ", ioriglen(j)
        call errcheck('subroutine getdimwt', ' ',1)
        return
      endif
     
        do 100 i = 1,ii
          wts(i) = coordvls(i1+i)
  100   continue

      else

        print*, 'The dimension number passed to getdimwt must be in'
        print*, 'the range 0 < i < 5'
        call errcheck('subroutine getdimwt', ' ',1)
        return

      endif

      return
      end




      subroutine getedges(ivar, idim, edges)
              
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Description:
c     -----------
c     This subroutine retrieves a vector of the coordinate positions
c        of grid-cell edges for dimension idim of variable ivar
c        (after the data have been retrieved by a call to getdata or
c        getvdata or after calling shape.
c
c     Usage:
c     ------
c       call getedges(ivar, idim, edges)
c
c        ivar = variable i.d.
c        idim = dimension number (i.e., 1, 2, 3, or 4, if you want the 
c                 edges associated with the first, second,
c                third or fourth dimension retrieved; 
c        edges = vector long enough to receive all the coordinate
c                  values for the requested dimension (should be 1 
c                  greater than the number of grid cells).
c     
c     ------
c     Date: 1/15/96
c     ----
c     
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
#include "usersget"
#include "cvartbl"
#include "cdimtbl"   

      integer ivar, idim
      real edges(*)
      integer i, j, i1, ii

      j = iabs(idim)
      j = iptdimrg(j,ivar)
      
      if (j .eq. 0) then
        print*, 'variable ', ivar, ' doesn''t have a dimension ', idim
        print*, 'or you have cleared the dimension table, so that'
        print*, 'information on the dimensions is no longer available.'
        call errcheck('subroutine getedges', ' ',1)
        return
      endif
            
      if ((idim .gt. 0) .and. (idim .lt. 5)) then
                  
        i1 = iptedge(j) - 1
        ii = idimlen(j) + 1

       if (i1 .lt. 0) then
        print*, 'Dimension vector was not saved but is needed.'
        print*, 
     &    'User should allow longer dimensions to be saved by calling'
        print*, "call defmisc('max dimension kept', 'integer', val)"
        print*, "where val .gt. ", ioriglen(j)
        call errcheck('subroutine getedges', ' ',1)
        return
      endif
      
        do 100 i = 1,ii
          edges(i) = coordvls(i1+i)
  100   continue

      else

        print*, 'The dimension number passed to getcoord must be in'
        print*, 'the range -5 < i < 0 or 0 < i < 5'
        call errcheck('subroutine getedges', ' ',1)
        return

      endif

      return
      end


   
 
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     subroutine getfield
c     -------
c
c     Description:
c     -----------
c     This subroutine reads a variable from a DRS file.
c
c     The variable may have up to four dimensions. The only restriction
c     is that the longitude dimension have, approximately, a constant
c     increment between successive values and all dimensions be
c     stored monotonically
c
c     Usage:
c     ------
c     call getfield(ivar, array)
c
c     where
c        ivar is the index pointing to the defined variable that will
c                be extracted
c        array contains retrieved data.
c
c     ------
c     Date: 6/5/95
c     ----
c
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine getfield(ivar, array)


c      common variables

      implicit none
#include "usersget"
#include "cvartbl"
#include "cdomain"
#include "cscalars"
      
c     passed variables

      real array(*)
      integer ivar 


c      local variables:

      real a1
      integer ierr, i, j, jj, lendatal
      integer  getfld

      do 40 i = 1,4
            namedim(i) = vdimname(i,ivar)
            domainc(i) = vdimdom(i,ivar)
            bdry(1,i) = vdimbdry(1,i,ivar) 
            bdry(2,i) = vdimbdry(2,i,ivar) 
            ibdry(1,i) = ivdimpln(1,i,ivar) 
            ibdry(2,i) = ivdimpln(2,i,ivar) 
            cycle(i) = vdimcycl(i,ivar)
            cwt(i) = vdimwt(i,ivar) 
   40 continue

      lendatal = lendata

c      if (lendata .eq. 0) then
c        ierr = getfld(ivar, array)
c        if (ierr .ne. 0) then
c          call errcheck(' ', 'subroutine getfield', 0)
c          return 
c        endif
c        lendata = lendatan
c      endif

      ierr = getfld(ivar, array)
      if (ierr .ne. 0) then
        call errcheck(' ', 'subroutine getfield', 0)
        return 
      endif

      lendata = lendatal
      jj = ivarlen(ivar)

      if (misstype .eq. 1) then
        a1 = abs(1.0e-5*amissing)
        do 100 j=1,jj
          if (abs(array(j)-amissing) .le. a1) array(j) = omit
 100    continue  
      else
        do 150 j=1,jj
          if (abs(array(j)) .eq. amissing) array(j) = omit
 150    continue
      endif
         
      return
      end             




c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     integer function getfld
c     -------
c
c     Description:
c     -----------
c     This function reads a variable from a DRS file.
c
c     The variable may have up to four dimensions. The only restriction
c     is that the longitude dimension have, approximately, a constant
c     increment between successive values and all dimensions be
c     stored monotonically
c
c     Usage:
c     ------
c  ierr = getfld(ivar, array)
c
c     Exit Status:
c     -----------
c     0 if success
c     1 if an error occurred
c
c     ------
c     Date: 1/19/92
c     ----
c
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer function getfld(ivar, array)


c      common variables

      implicit none
#include "drsdef.h"
#include "fcddrs.h"

#include "usersget"
#include "cscalars"
#include "cvartbl"
#include "cdimtbl"
#include "cdomain"

c     passed variables

      real array(*)
      integer ivar


c      local variables:


      character dname(4)*16, dsource(4)*120, dtitle(4)*80, dunits(4)*40
      logical caseindp


      integer ierr, ndim, i, jj, idicun, idtype(4), irow, lvarlen,
     +         iorders(4), lendim(4), idlen(4), icnt, jcnt, ll, j, 
     &         nm, idrs, ierr1, ierr2, ierr3, istore, k, isize
cjfp was      integer cluvdb, fcw_getslab, getdat
      integer fcw_cluvdb, fcw_getslab, fcw_getdat
cjfp was     integer setname, getedim, getcdim, setdim
      integer fcw_setname, fcw_getedim, fcw_getcdim, fcw_setdim
      integer fcw_error
      integer
     +     opendrs, getdimen, mkdimtbl,
     +     univunit, miallocw, mifree

      real fes, les, felrloc(4), lelrloc(4), delta(4), cycleloc(4)
      pointer (ptv, vecdim), (ptw, twt), (pte, tedge)
      real vecdim(*)
      double precision tedge(*), twt(*)

      integer drstest

      getfld = 0
      
      do 10 i=1,4
        lendim(i) = 0
        idlen(i) = 0
        iorders(i) = 0
        felrloc(i) = 0.0
        lelrloc(i) = 0.0
        cycleloc(i) = 0.0
        delta(i) = 0.0
c       These strings should have values when passed into other routines.
c       Note that Fortran pads with blanks.
        dname(i) = " "
        dsource(i) = " "
        dtitle(i) = " "
        dunits(i) = " "
   10 continue
      
c      check that variable and file have been defined
      if (caseindp(varnam(ivar), 'not def') .or. 
     +    caseindp(filnam(ivar), 'not def')) then
         print*,'you must call defvar to define variable name and file'
         print*, 'for the variable identified by index: ', ivar
         getfld = 1
         call errcheck('function getfld', ' ',1)
         return
      endif

      ierr = fcw_cluvdb()
      if (ierr .ne. 0) then
         getfld = ierr
         call errcheck('DRS function fcw_cluvdb', 'function getfld',1)
         return
      endif

      ierr = opendrs(ivar, 1, idicun)
      if (ierr .ne. 0) then
         getfld = ierr
         call errcheck(' ', 'function getfld',1)
         return
      endif

c         get info on each dimension
      ndim = ivarndim(ivar)
      do  300 i=1,ndim
         ierr = fcw_getedim(i,dsource(i),dname(i),dtitle(i),dunits(i),
     $        idtype(i) , idlen(i), fes, les)
         if (ierr .ne. 0) then
            getfld = ierr
            call errcheck('DRS function fcw_getedim',
     1           'function getfld',1)
            return
         endif

  300 continue

c     Check for dimensions that are defined by user but not in 
c          file

      icnt = 0
      jcnt = 0
      do 410 i=1,4
        if (namedim(i) .ne. 'not defined') then

          icnt = icnt + 1

          do 404 j=1,ndim
            if ((caseindp(namedim(i), dname(j)) .or. 
     &         ((univunit(namedim(i)) .gt. 0) .and.
     &          (univunit(namedim(i)) .eq. univunit(dname(j)))))) then
              go to 410
            endif
  404     continue

          jcnt = jcnt + 1

          if (errcntl .gt. 0) then
            print*, 'Warning --'
            print*, 'you have defined a dimension named ', namedim(i)
            print*, 'but there is no dimension by this name for the' 
            print*, 'variable named ', varnam(ivar), ' in file ',
     +                filnam(ivar) 
            print*,
     &          'This dimension will be ignored, at least temporarily.'
          endif
          call errcheck(' ', ' ',-1)

        endif
  410 continue

c   set pointers to dimensions requested

      do 400 i=1,ndim

        do 350 j=1,4                               

          if (caseindp(namedim(j), dname(i)) .or. 
     &           ((univunit(dname(i)) .gt. 0) .and. 
     &           (univunit(namedim(j)) .eq. univunit(dname(i))))) then
            iorders(j) = i
            go to 400
          endif

  350   continue

        if ((idlen(i) .eq. 1) .and. (icnt .lt. 4)) then

          if (errcntl .gt. 0) then
            print*, 'Warning --'
            print*,
     &         'A dimension has been found that you have not defined,'
            print*, 'but since the dimension is length 1 (and therefore'
            print*, 'probably a "dummy" dimension), EzGet will proceed'
            print*, 'as if you had defined it.'
            print*, 'Dimensions in file: ', (dname(nm), nm=1,ndim)
            print*, 'You are in the process of extracting a'
            print*, 'variable named ', varnam(ivar), ' in file ',
     +                filnam(ivar) 
          endif
          call errcheck(' ', ' ',-1)

          icnt = icnt + 1
          namedim(icnt) = dname(i)
          domainc(icnt) = 'assaved'
          cwt(icnt) = 'unit'
          iorders(icnt) = i
          cycle(icnt) = 0.0

        else

          print*, 'You must call defdim or defdimi to define the '
          print*, 'dimension named ', dname(i), 'for variable ',
     +            varnam(ivar) 
          print*, 'in file,', filnam(ivar)
          print*, 'Dimensions in file: ', (dname(nm), nm=1,ndim)
          print*, 'Dimensions defined by user: ', (namedim(nm), nm=1,4)
          if (icnt .ge. 4) then
            print*, 'Also check that all the dimensions you defined '
            print*, 'are actually needed for this variable.'
          endif
          getfld = 1
          call errcheck('function getfld', ' ',1)
          return

        endif

  400 continue

c    cycle through each dimension
c     set index information for retrieving data

      wrap = .false.
      do 1000 i=1,4

       cycleloc(i) = cycle(i)

       if (iorders(i) .ne. 0) then                                            
        idrs = iorders(i)

c              retrieve dimension
c           allocate space for vecdim, tedge, twt
        lendim(i) = idlen(idrs)
        ierr1 = miallocw(lendim(i), ptv)

        if (lendim(i) .le. mxdimelm) then
          ll = 2*max0(ittrunc,lendim(i))
          ierr2 = miallocw(ll, ptw)
          ll = 2*max0(ittrunc,(lendim(i)+1))
          ierr3 = miallocw(ll, pte)
          istore = 1
        elseif (.not. caseindp(cwt(i),'uni')) then

          print*, 'You are attempting to retrieve a dimension that is'
          print*, 'longer than maximum allowed which is ', mxdimelm
          print*, 'It is possible to retrieve this dimension if '
          print*, 'you specify "unit" for the weighting option, or'
          print*, 'you may increase the limit on dimension size by'
          print*, 'calling defmisc with the "longest dimension"'
          print*, 'option and its value set larger than ', lendim(i)
          print*, ' '
          call errcheck('function getfld',' ',2)
          getfld = 1
          return

        else
          ierr2 = miallocw(2, ptw)
          ierr3 = miallocw(4, pte)
          istore = 0
        endif

        if ((ierr1 .ne. 0) .or. (ierr2 .ne. 0) .or. 
     &          (ierr3 .ne. 0)) then
           print*, 'Error in attempting to allocate memory dynamically.'
           print*, 'You may need more memory.'
           print*, 'Try running on cirrus or retrieving smaller '
           print*, 'portions of the dimension.  Error encountered '
           print*, 'while processing', filnam(ivar), varnam(ivar)
           getfld = 1
           call errcheck('function miallocw','function getfld',1)
           return
        endif

        ierr = getdimen(idrs, domainc(i), cycle(i), ibdry(1,i), 
     +          ibdry(2,i), bdry(1,i), bdry(2,i), felrloc(i),
     +          lelrloc(i), delta(i), cwt(i), lendim(i),
     +          dsource(i), dname(i), dtitle(i), dunits(i), idtype(i),
     &          istore, vecdim, tedge, twt)

        if (ierr .lt. 0) then
          ierr = mifree(ptv)
          ierr = mifree(ptw)
          ierr = mifree(pte)
          ierr1 = miallocw(lendim(i), ptv)
          if (lendim(i) .le. mxdimelm) then
            ll = 2*max0(ittrunc,lendim(i))
            ierr2 = miallocw(ll, ptw)
            ll = 2*max0(ittrunc,(lendim(i)+1))
            ierr3 = miallocw(ll, pte)
            istore = 1
          elseif (.not. caseindp(cwt(i),'uni')) then

            print*, 'You are attempting to retrieve a dimension that is'
            print*, 'longer than maximum allowed which is ', mxdimelm
            print*, 'It is possible to retrieve this dimension if '
            print*, 'you specify "unit" for the weighting option, or'
            print*, 'you may increase the limit on dimension size by'
            print*, 'calling defmisc with the "longest dimension"'
            print*, 'option and its value set larger than ', lendim(i)
            print*, 'Because you have invoked the wrap-around option,'
            print*, '"longest dimension" might need to be larger than'
            print*, 'you expected.'
            print*, ' '
            call errcheck('function getfld',' ',2)
            getfld = 1
            return

          else
            ierr2 = miallocw(2, ptw)
            ierr3 = miallocw(4, pte)
            istore = 0
          endif

          if ((ierr1 .ne. 0) .or. (ierr2 .ne. 0) .or. 
     &          (ierr3 .ne. 0)) then
           print*, 'Error in attempting to allocate memory dynamically.'
           print*, 'You may need more memory.'
           print*, 'Try running on cirrus or retrieving smaller '
           print*, 'portions of the dimension.  Error encountered '
           print*, 'while processing', filnam(ivar), varnam(ivar)
           getfld = 1
           call errcheck('function miallocw','function getfld',1)
           return
          endif

          ierr = getdimen(idrs, domainc(i), cycle(i), ibdry(1,i), 
     +          ibdry(2,i), bdry(1,i), bdry(2,i), felrloc(i),
     +          lelrloc(i), delta(i), cwt(i), lendim(i),
     +          dsource(i), dname(i), dtitle(i), dunits(i), idtype(i),
     +          istore, vecdim, tedge, twt)
          if (ierr .lt. 0) then
            print*, 'Error: some problem in allocating array to '
            print*, 'retrieve dimension.'
          endif

        endif
          

        if (ierr .ne. 0) then
           print*, 'while processing', filnam(ivar), varnam(ivar)
           getfld = ierr
           call errcheck('function getdimen', 'function getfld',0)
           return
        endif

c         check whether this dimension appears in dimension table
c         if not, add it and also store vecdim, tedge, twt in coordvls
         jj = lendim(i)
         idrs = iorders(i)
        ierr=mkdimtbl(mxdmkeep, istore, dname(i), lendim(i),  
     +          idlen(idrs), dsource(i), dtitle(i), dunits(i),      
     &          idtype(i), domainc(i), cycle(i), 
     +          bdry(1,i), bdry(2,i),vecdim(1), vecdim(jj), irow, 
     +          cwt(i), vecdim, tedge, twt)
         if (ierr .ne. 0) then
             getfld = ierr
             call errcheck(' ', 'function getfld',0)
             return
         endif

         ierr = mifree(ptv)
         ierr = mifree(ptw)
         ierr = mifree(pte)
         
         iptdim(i,ivar) = irow
         iptdimrg(i,ivar) = irow

       else

         iptdim(i,ivar) = 0
         iptdimrg(i,ivar) = 0
      
       endif    
 1000 continue

      i = 1
      do 1100 k = 1,4

c          remove dimensions that were missing 
           if ((iorders(i) .eq. 0) .and. (i .lt. 4)) then
             do 1090 j=i,3
               iorders(j) = iorders(j+1)
               lendim(j) = lendim(j+1)
               felrloc(j) = felrloc(j+1)
               lelrloc(j) = lelrloc(j+1)
               delta(j) = delta(j+1)
               cycleloc(j) = cycleloc(j+1)
 1090        continue
             iorders(4) = 0
             lendim(4) = 0
             felrloc(4) = 0.0
             lelrloc(4) = 0.0
             cycleloc(4) = 0.0
             delta(4) = 0.0
           else
             i = i + 1
           endif 

 1100 continue

      lvarlen = 1
      do 1110 i=1,4   
                                     
           lvarlen=lvarlen*max0(lendim(i),1)

 1110 continue

       lendatan = lvarlen
       if (lendata .lt. 0) then
          return
       elseif ((lvarlen .gt. lendata) .and. (lendata .ne. 0)) then
          print*, 'data array was dimensioned to be too small while'
          print*, 'preparing to retrieve ', filnam(ivar), varnam(ivar)
          print*, 'dimensioned: ', lendim
          print*, 'If calling getdata, the 4 dimensions should equal or'
          print*, 'exceed the above values'
          print*, 'If calling getvdata, define scalar lendata'
          print*, 'to be at least', lvarlen
          print*, 'Also dimension the array and mask to be at least' 
          print*, 'that size.  Note that if the wrap-around or '
          print*, 'cycling option is in effect and if'
          print*, 'domainc = range for the dimension, then one'  
          print*, 'grid cell or time sample might be split into 2, so'
          print*, 'you might need to increase that dimension by 1.'
          print*, 'Alternatively, to prevent splitting, specify'
          print*, 'domainc = nearest.'  
          getfld = 1
          call errcheck('function getfld', ' ',1)
          return
       endif


      ierr = fcw_cluvdb()
      if (ierr .ne. 0) then
         getfld = 1
         call errcheck('DRS function fcw_cluvdb', 'function getfld',1)
         return
      endif
         
      ierr = fcw_setname(varsourc(ivar),varnam(ivar),vartitle(ivar),    
     +      varunits(ivar), ' ')
      if (ierr .ne. 0) then
         getfld = 1
         call errcheck('DRS function setname', 'function getfld',1)
         return
      endif

      if (wrap) then
        if (drstest( fcw_getslab(idicun, ndim, iorders, felrloc,
     +       lelrloc, cycleloc, array, lendim) ) .ne. 0 ) then
          print *,'DRS error reading data for ', filnam(ivar),
     +            varnam(ivar) 
          getfld = 1
          call errcheck('DRS function getslab', 'function getfld' ,1)
          return
        endif
      else
        do 1200 j=1,ndim
          i = iorders(j)
c          print*, dname(i)
c          print*, lendim(j)
c          print*, felrloc(j)
c          print*, lelrloc(j)
          ierr = fcw_setdim(j, dname(j), ' ', lendim(j), felrloc(j),
     +             lelrloc(j))                                 
          if (ierr .ne. 0) then
            print*, 'DRS error while setting dimension, ', dname(i)
            print*, 'preparing to read data from ', filnam(ivar),
     +                varnam(ivar) 
            getfld = 1
            call errcheck('DRS function fcw_setdim','function getfld',1)
            return
          endif
 1200   continue
        If (vartype(ivar) .eq. 'R*8') then
c           IDRS_BYTESPERWORD is a parameter in drsdef.h which is
c              4 on Suns and 8 on Crays.
c   ????
          isize = lvarlen*8
        else
           isize = lvarlen*4
        endif
c ???   check if isize is correct argument
        ierr = fcw_getdat(idicun, array, isize)
        if (ierr .ne. 0) then 
           print *, 'DRS error reading data from ', filnam(ivar),
     +              varnam(ivar) 
          getfld = 1
          call errcheck('DRS function fcw_getdat','function getfld',1)
          return
        endif
      endif

c Get the dimension vectors for the data AS RETRIEVED. Note that this
c must be done after the call to getslab and before any calls to fcw_cluvdb.
c Also put info. into vartable and dimtable


      ivarlen(ivar) = 1
      do 1300 i=1,ndim
         ierr = miallocw(lendim(i), ptv)

         if (ierr .ne. 0) then
           print*, 'Error in attempting to allocate memory dynamically.'
           print*, 'You may need more memory.'
           print*, 'Try running on cirrus or retrieving smaller '
           print*, 'portions of the dimension.  Error encountered '
           print*, 'while processing', filnam(ivar), varnam(ivar)
           getfld = 1
           call errcheck('function miallocw','function getfld',1)
           return
          endif



         ierr = fcw_getcdim(i, dsource(i), dname(i), dtitle(i), 
     &          dunits(i), idtype(i), lendim(i), vecdim, jj)
         ivarlen(ivar) = ivarlen(ivar)*jj
         if (jj .ne. lendim(i)) then
           if (errcntl .gt. 0) then
           print*, 'Warning: Possible error in program getfld'
           print*, 'Expected number of elements in vecdim vector = ',
     +           lendim(i)
           print*, 'Actual number of elements in vecdim vector = ', jj
           lendim(i) = jj
           endif
           call errcheck('function getfld',' ',-1)
         endif

         if (abs(felrloc(i)+delta(i)-vecdim(1)) .gt.
     +              1.e-5*(abs(lelrloc(i)-felrloc(i)) + 
     +                     abs(vecdim(jj)-vecdim(1)))) then
           if (errcntl .gt. 0) then
           print*, 'Warning: Possible error in program getfld'
           print*, 'Expected first element of vecdim vector = ',
     +           felrloc(i)
          print*, 'Actual first element of vecdim vector = ', vecdim(1)
           endif
           call errcheck('function getfld',' ',-1)
         endif

         if (abs(lelrloc(i)-delta(i)-vecdim(jj)) .gt.
     +              1.e-5*(abs(lelrloc(i)-felrloc(i)) + 
     +                     abs(vecdim(jj)-vecdim(1)))) then
           if (errcntl .gt. 0) then
           print*, 'Warning: Possible error in program getfld'
           print*, 'Expected last element of vecdim vector = ',
     +           lelrloc(i)
          print*, 'Actual last element of vecdim vector = ', vecdim(jj)
           endif
           call errcheck('function getfld',' ',-1)
         endif
         if (ierr .ne. 0) then
            getfld = 1
            call errcheck('DRS function fcw_getcdim',
     1           'function getfld',1)
            return
         endif

         ierr = mifree(ptv)
         
 1300 continue

      return
      end




c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     subroutine getgeog
c     -------
c
c     Description:
c     -----------
c     This subroutine reads a geography mask from a DRS file.
c     It checks whether the dimensions expected are correct
c     or at least that the arrays are dimensioned large enough to 
c     accomodate the retrieved data.
c
c     
c     The variable may have up to four dimensions. The only restriction
c     is that the longitude dimension have, approximately, a constant
c     increment between successive values and all dimensions be
c     stored monotonically
c
c     Usage:
c     ------
c     call getgeog(ivar, ndim1, ndim2, ndim3, ndim4,
c              len1, len2, len3, len4, amask, select)
c
c
c     where
c        ivar is the index pointing to the defined variable that will
c                be extracted
c        ndim1 is the length of the 1st dimension of amask and array
c               as declared in the calling subroutine
c        ndim2, ndim3, ndim4 are defined as ndim1 but for 2nd, 3rd, and
c                4th dimensions
c
c    Note:   ndim1*ndim2*ndim3*ndim4  (where values equal to 0 are
c               reset to 1) should be large enough to
c               accomodate the extracted data array even if it
c               is subsequently regridded to a smaller size. 
c
c        len1  on entry:
c                the expected 1st dimension of the data that
c                is being retrieved. (If user specifies len1 to be 0,
c                then he will not be warned when the actual 1st 
c                dimension of the data being retreived is different
c                from what he expects.)
c               on return:
c                the actual 1st dimension of the data that is being
c                    retrieved.
c        len2, len3, len4 are defined as len1, but for 2nd, 3rd, and
c                4th dimensions
c        amask contains 0. or 1. depending on whether on not the grid
c               cell should be masked or not.
c        select is a string containing a list of the regions selected 
c                (i.e., left unmasked
c                         
c
c     ------
c     Date: 1/5/95
c     ----
c
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine getgeog(ivar, ndim1, ndim2, ndim3, ndim4, 
     +                        len1, len2, len3, len4, amask, select)
     
      implicit none
#include "usersget"
#include "cvartbl"  
#include "cdomain"  
#include "cscalars"

      integer ivar, ndim1, ndim2, ndim3, ndim4, len1, len2, len3, len4
      character*(*) select
     
c      real amask(ndim1,ndim2,ndim3,ndim4),
c     +     array(ndim1,ndim2,ndim3,ndim4)
      real amask(*)
      integer univunit, iii(4), ll(4), n(4), getfld, xregion, ishowl,
     &     igeoginl, lendatal, n1, n2, n3, n4, j, ierr, i, lenarray,
     &     kkk, m, l, lll, k, jjj

      ishowl = ishowin(ivar)
      igeoginl = igeogin(ivar)
            
      lendatal = lendata
      n1 = max0(1,ndim1)
      n2 = max0(1,ndim2)
      n3 = max0(1,ndim3)
      n4 = max0(1,ndim4)
      n(1) = ndim1
      n(2) = ndim2
      n(3) = ndim3
      n(4) = ndim4
      ll(1) = len1
      ll(2) = len2
      ll(3) = len3
      ll(4) = len4
      
      lendata = n1*n2*n3*n4
      
      do 50 j=1,4
            
         if (ll(j) .gt. n(j)) then 
            ierr = 1
            print*, 'In calling getgeog you have specified dimension ',
     +          j, ' to be smaller than'
            print*, 'the length you expect it to be'
            print*, 'You indicate that this dimension has a length of: '
     +          , n(j) 
            print*, 'but you expect the length to be: ', ll(j)
         endif
         
   50 continue                                                

      if (ierr .gt. 0) then
        print*, 'error occurred in subroutine getgeog'
        print*, 'while processing', filnam(ivar), varnam(ivar)
        call errcheck(' ', ' ',1)
  
        print*, ' '
        print*, ' '
      endif      

      do 70 i = 1,4
            namedim(i) = vdimname(i,ivar)
            domainc(i) = vdimdom(i,ivar)
            bdry(1,i) = vdimbdry(1,i,ivar) 
            bdry(2,i) = vdimbdry(2,i,ivar) 
            ibdry(1,i) = ivdimpln(1,i,ivar) 
            ibdry(2,i) = ivdimpln(2,i,ivar) 
            cycle(i) = vdimcycl(i,ivar)
            cwt(i) = vdimwt(i,ivar) 
   70 continue

      
      ierr = getfld(ivar, amask)

      call defgeog(ivar, 'in', ivar, select)

      masktype = 0
      if (masktypi .ne. 0) then
        masktype = masktypi
      elseif (masktypo .ne. 0) then
        masktype = masktypo
      endif

      ierr = xregion(ishowin(ivar), vartype(ivar), varnam(ivar),
     &          ivarlen(ivar), amask)
      
      call lendims(ivar, iii(1), iii(2), iii(3), iii(4), lenarray)
             
      
      ierr = 0
      
      do 100 j=1,4
      
         if (iii(j) .gt. n(j)) then
            ierr = 1
            print*, 'Error: Dimension ', j, ' of the array and mask '
            print*, 'passed to getgeog is too small to accomodate'
            print*, 'the data you have requested.'
            print*, 'Either change the domain requested or increase'
            print*, 'this dimension of the arrays to ', iii(j) 
            if ((iii(j) .eq. n(j)+1) .and. (cycle(j) .ne. 0.)) then
              kkk = univunit(namedim(j)) 
              if (kkk .eq. intlong) then 
                print*, 'You have invoked the wrap-around or cycle'
                print*, 'option, which may split 1 longitude'
                print*, 'cell into 2 cells, so you probably need to'
                print*, 'increase your longitude dimension by 1.'
                print*, 'Alternatively, for this dimension set'
                print*, 'domainc(',j,')= nearest and the range will'  
                print*, 'be shifted so as not to split the cell.'
              endif      
            endif
         endif
         
         if ((ll(j) .gt. 0) .and. (iii(j) .ne. ll(j))) then
            
           if (errcntl .gt. 0) then
            print*,'Warning: In calling getgeog you expected dimension '
     +           , j
            print*, 'to be ', ll(j), ' but in fact it turns out to be '
     +           , iii(j) 
            print*, 'To correct this problem, set the expected'
            print*, 'length of this dimension to ', iii(j)
            print*, 'or to supress this warning, set expected length'
            print*, 'of this dimension equal to 0.'
            endif
            call errcheck('subroutine getgeog', ' ',-1)
         endif

  100 continue

      if (ierr .gt. 0) then
        print*, 'error occurred in subroutine getgeog'
        print*, 'while processing', filnam(ivar), varnam(ivar)
        call errcheck(' ', ' ',1)
  
        print*, ' '
        print*, ' '
      
      else
                 
      len1 = max0(1,iii(1))
      len2 = max0(1,iii(2))
      len3 = max0(1,iii(3))
      len4 = max0(1,iii(4))
      
        m = lenarray + 1
        do 500 l=len4,1,-1
           lll = n3*(l-1)
           do 400 k=len3,1,-1
              kkk = n2*(lll+k-1)
              do 300 j=len2,1,-1
                 jjj = n1*(kkk+j-1)
                 do 200 i=len1,1,-1
                    m = m - 1
                    amask(jjj+i) = amask(m)
  200            continue
  300         continue
  400      continue
  500   continue
  
       if (len4 .lt. n4) then
        do 640 l=len4+1,n4
           lll = n3*(l-1)
           do 630 k=1,n3
              kkk = n2*(lll+k-1)
              do 620 j=1,n2
                 jjj = n1*(kkk+j-1)
                 do 610 i=1,n1
                    amask(jjj+i) = 0.0
  610            continue
  620         continue
  630      continue
  640   continue
       endif

       if (len3 .lt. n3) then
        do 740 k=len3+1,n3
           kkk = n1*n2*(k-1)
           do 730 l=1,len4
              lll = kkk + n1*n2*n3*(l-1)
              do 720 j=1,n2
                 jjj = lll + n1*(j-1)
                 do 710 i=1,n1
                    amask(jjj+i) = 0.0
  710            continue
  720         continue
  730      continue
  740   continue
       endif
            
       if (len2 .lt. n2) then   
        do 840 j=len2+1,n2
           jjj = n1*(j-1)
           do 830 l=1,len4
              lll = jjj + n1*n2*n3*(l-1)
              do 820 k=1,len3
                 kkk = lll + n1*n2*(k-1)
                 do 810 i=1,n1
                    amask(kkk+i) = 0.0
  810            continue
  820         continue
  830      continue
  840   continue
       endif
               
       if (len1 .lt. n1) then
        do 940 i=len1+1,n1
           do 930 l=1,len4
              lll = n3*(l-1)
              do 920 k=1,len3
                 kkk = n2*(lll+k-1)
                 do 910 j=1,len2
                    jjj = n1*(kkk+j-1)
                    amask(jjj+i) = 0.0
  910            continue
  920         continue
  930      continue
  940   continue
       endif
    
      endif
                 
      len1 = iii(1)
      len2 = iii(2)
      len3 = iii(3)
      len4 = iii(4)
      
      lendata = lendatal
      ishowin(ivar) = ishowl
      igeogin(ivar) = igeoginl

      return
      end
      


      
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     subroutine getnogap
c     -------
c
c     Description:
c     -----------
c     This subroutine reads a variable from a file and places it in
c        contiguous memory.
c
c     The variable may have up to four dimensions. The only restriction
c     is that the longitude dimension have, approximately, a constant
c     increment between successive values and all dimensions be
c     stored monotonically
c
c     Usage:
c     ------
c     call getnogap(ivar, len, amask, array)
c
c     where
c        ivar is the index pointing to the defined variable that will
c                be extracted
c        len is the dimension of amask, and array
c        amask contains weights associated with each element in array.
c        array contains retrieved data.
c
c     ------
c     Date: 2/28/97
c     ----
c
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine getnogap(ivar, len, amask, array)

      implicit none
#include "cscalars"
     
      integer ivar, len, lendatal
      real amask(len), array(len)

      lendatal = lendata
      lendata = len

      call getvdata(ivar, amask, array)

      lendata = lendatal

      return
      end




c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     subroutine getvdata
c     -------
c
c     Description:
c     -----------
c     This subroutine reads a variable from a DRS file.
c
c     The variable may have up to four dimensions. The only restriction
c     is that the longitude dimension have, approximately, a constant
c     increment between successive values and all dimensions be
c     stored monotonically
c
c     Usage:
c     ------
c     call getvdata(ivar, amask, array)
c
c     where
c        ivar is the index pointing to the defined variable that will
c                be extracted
c        amask contains weights associated with each element in array.
c        array contains retrieved data.
c
c     ------
c     Date: 1/19/92
c     ----
c
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine getvdata(ivar, amask, array)


c      common variables

      implicit none
#include "usersget"
#include "cscalars"
#include "cvartbl"
#include "cdimtbl"
#include "cdomain"
      
c     passed variables

      real array(*), amask(*)
      integer ivar 


c      local variables:

      pointer (p1, amaskin), (p2, arrayin)           
      real amaskin(*), arrayin(*)
      integer ierr, j, jj, i, iloni, ilati, iii, ijk, irlono, irlato,
     &    ii, lendatal, ilono, ilato, ierr1, ierr2, len
      integer  xgeog, applywts, getfld, doregrid, univunit, miallocw,
     &    mifree
      real adum(1), a1
      logical caseindp

      if (caseindp(vdimname(1,ivar),  'not def')) then

         print*,'Error: You must define dimensions for the variable'
         print*, 'being extracted: ',  varnam(ivar)
         print*, 'stored in file: ', filnam(ivar)
         call errcheck('subroutine getvdata', ' ', 1)

      else

         do 40 i = 1,4
            namedim(i) = vdimname(i,ivar)
            domainc(i) = vdimdom(i,ivar)
            bdry(1,i) = vdimbdry(1,i,ivar) 
            bdry(2,i) = vdimbdry(2,i,ivar) 
            ibdry(1,i) = ivdimpln(1,i,ivar) 
            ibdry(2,i) = ivdimpln(2,i,ivar) 
            cycle(i) = vdimcycl(i,ivar)
            cwt(i) = vdimwt(i,ivar) 
   40    continue

      endif

      if (iregrid(ivar) .ne. 0) then

         iloni = 0
         ilati = 0
         do 45 i = 1,4
            iii = univunit(namedim(i))
            if ((iloni .eq. 0) .and. (iii .eq. intlong)) iloni = i
            if ((ilati .eq. 0) .and. (iii .eq. intlat))  ilati = i    
  45     continue

         if (cycle(iloni) .eq. 0.) then 
c           don't override cycle prescription because if only a 
c           limited domain is stored, then program won't work.
c           cycle(iloni) = 360.
           if (errcntl .gt. 0) then
            print*, 'Warning:  You are regridding the data, so'
            print*, 'instead of setting cycle = 0.0 for longitude,'
            print*, 'it is usually prudent to set the longitude'
            print*, 'cycle to 360. or 2*pi, or whatever is appropriate.'
            print*, 'If, however, input data are not available for a'
            print*, 'complete cycle, you should not change the current'
            print*, 'specification (i.e., cycle = 0.0)'    
c           print*, 'ezget has set cycle to 360. before proceeding.'
           endif
           call errcheck('subroutine getvdata', ' ',-1)
         endif

         if (iregrid(ivar) .lt. 0) then
           if ((dlatrgdv(ivar) .ne. 0.0) .and. (nlatrgdv(ivar) .gt. 0))
     &           then
              domainc(ilati) = 'range'
              bdry(1,ilati) = alatrgdv(ivar) - 0.5*dlatrgdv(ivar)
              bdry(2,ilati) = alatrgdv(ivar) +
     &                        (nlatrgdv(ivar)-0.5)*dlatrgdv(ivar)
           endif
           if (nlonrgdv(ivar) .gt. 0) then
              domainc(iloni) = 'range'
              bdry(1,iloni) = alonrgdv(ivar) - 0.5*dlonrgdv(ivar)
              bdry(2,iloni) = alonrgdv(ivar) + 
     &                        (nlonrgdv(ivar)-0.5)*dlonrgdv(ivar)
           endif
         endif

      endif

      if ((iregrid(ivar) .gt. 0) .and.(iregrid(ivar) .le. maxvars)) then
                         
c          get longitude-latitude range from target grid
                 
         ijk = iregrid(ivar)

         if (iregrid(ijk) .ne. 0) then
           print*, 'Error ... You have specified a target grid by '
           print*, 'variable index, ', ijk, ' but you have also '
           print*, 'specified that this variable be regridded.'
           print*, 'EZGET will not permit target grid variables to be'
           print*, 'mapped to a new grid.' 
           print*, 'Your target grid was specified as: ', varnam(ijk)
           print*, 'stored in file: ', filnam(ijk)
           call errcheck('subroutine getvdata', ' ', 1)
         endif

         irlono = 0
         irlato = 0
         do 48 i = 1,4
            ii = iptdim(i,ijk)
            iii = 0
            if (ii .ne. 0)  iii = univunit(dimnam(ii))
            if (iii .eq. intlong) irlono = ii
            if (iii .eq. intlat)  irlato = ii               
  48     continue

         if ((irlono .eq. 0) .or. (irlato .eq. 0)) then

c          we have not retrieved target grid yet so:


           if (caseindp(vdimname(1,ijk),  'not def')) then

             print*,'Error: You must define dimensions for target grid.'
             print*, 'Your target grid was specified as: ', varnam(ijk)
             print*, 'stored in file: ', filnam(ijk)
             call errcheck('subroutine getvdata', ' ', 1)

           else

c             retrieve target grid according to user specifications
c               for domain of target-grid variable ijk.

             do 50 i = 1,4
                namedim(i) = vdimname(i,ijk)
                domainc(i) = vdimdom(i,ijk)
                bdry(1,i) = vdimbdry(1,i,ijk) 
                bdry(2,i) = vdimbdry(2,i,ijk) 
                ibdry(1,i) = ivdimpln(1,i,ijk) 
                ibdry(2,i) = ivdimpln(2,i,ijk) 
                cycle(i) = vdimcycl(i,ijk)
                cwt(i) = vdimwt(i,ijk) 
   50        continue

           endif 
      
           lendatal = lendata
           lendata = -1  
           ierr = getfld(ijk, adum)

           if (ierr .ne. 0) then
             call errcheck(' ', 'subroutine getvdata', 0)
             return
           endif
    
           lendata = lendatal 

                        
c        look at dimensions for target grid, extract longitude
c           and latitude 

           irlono = 0
           irlato = 0
           ijk = iregrid(ivar)
           do 60 i = 1,4
              ii = iptdim(i,ijk)
              iii = 0
              if (ii .ne. 0)  iii = univunit(dimnam(ii))
              if (iii .eq. intlong) irlono = ii
              if (iii .eq. intlat)  irlato = ii               
  60       continue

           if (irlono .eq. 0) then
             print*, 'You are attempting to obtain a target grid from a'
             print*, 'variable that is not a function of longitude' 
             call errcheck('subroutine getvdata', ' ',1)
             return
           endif

           if (irlato .eq. 0) then
             print*, 'You are attempting to obtain a target grid from a'
             print*, 'variable that is not a function of latitude' 
             call errcheck('subroutine getvdata', ' ',1)
             return
           endif

           do 70 i = 1,4
              namedim(i) = vdimname(i,ivar)
              domainc(i) = vdimdom(i,ivar)
              bdry(1,i) = vdimbdry(1,i,ivar) 
              bdry(2,i) = vdimbdry(2,i,ivar) 
              ibdry(1,i) = ivdimpln(1,i,ivar) 
              ibdry(2,i) = ivdimpln(2,i,ivar) 
              cycle(i) = vdimcycl(i,ivar)
              cwt(i) = vdimwt(i,ivar) 
   70      continue

           if (cycle(iloni) .eq. 0.) cycle(iloni) = 360.
      
         endif 

         domainc(iloni) = 'range'
         bdry(1,iloni) = begbdy(irlono)
         bdry(2,iloni) = endbdy(irlono)

         domainc(ilati) = 'range'
         bdry(1,ilati) = begbdy(irlato)
         bdry(2,ilati) = endbdy(irlato)

      endif


      if (lendata .lt. 0) then
                 
c       find out how much space is needed but don''t retrieve data

        ierr = getfld(ivar, array)

        if (ierr .ne. 0) then
          call errcheck(' ', 'subroutine getvdata', 0)
          return
        endif

        if (iregrid(ivar) .ne. 0) then

          if (iregrid(ivar) .lt. 0) then

c           get longitude-latitude range from source grid
                 
            ilono = 0
            ilato = 0
            do 78 i = 1,4
              ii = iptdim(i,ivar)
              iii = 0
              if (ii .ne. 0)  iii = univunit(dimnam(ii))
              if (iii .eq. intlong) ilono = ii
              if (iii .eq. intlat)  ilato = ii               
   78       continue
c ?/?
            domainc(iloni) = 'range'
            if (nlonrgdv(ivar) .eq. 0) then
              bdry(1,iloni) = begbdy(ilono)
              bdry(2,iloni) = endbdy(ilono)
            else
              bdry(1,iloni) = alonrgdv(ivar) - 0.5*dlonrgdv(ivar)
              bdry(2,iloni) = alonrgdv(ivar) + 
     &                        (nlonrgdv(ivar)-0.5)*dlonrgdv(ivar)
            endif

            domainc(ilati) = 'range'
            if ((nlatrgdv(ivar) .gt. 0) .and. 
     &          (dlatrgdv(ivar) .ne. 0.0)) then
              bdry(1,ilati) = alatrgdv(ivar) - 0.5*dlatrgdv(ivar)
              bdry(2,ilati) = alatrgdv(ivar) +
     &                        (nlatrgdv(ivar)-0.5)*dlatrgdv(ivar)
            else
              bdry(1,ilati) = begbdy(ilato)
              bdry(2,ilati) = endbdy(ilato)
            endif

          endif

          ierr = doregrid(ivar, amaskin, arrayin, amask, array)

          if (ierr .ne. 0) then
            call errcheck(' ', 'subroutine getvdata', 0)
          endif

        endif

        return

      endif     

      if ((ishowin(ivar) .ne. all) .and. (igeogin(ivar) .eq. 0)) then
        if (errcntl .gt. 0) then
            print*, 'Warning:  If you want to mask geography before '
            print*, 'regridding, you must call defgeog.'
            call errcheck('subroutine getvdata', ' ',-1)
         endif
      endif
      if ((ishowout(ivar) .ne. all) .and. (igeogout(ivar) .eq. 0)) then
        if (errcntl .gt. 0) then
            print*, 'Warning:  If you want to mask geography after ' 
            print*, 'regridding, you must call defgeog.'
            call errcheck('subroutine getvdata', ' ',-1)
         endif
      endif


c   if regridding data, find out how much space is needed for original
c        grid

      if (iregrid(ivar) .ne. 0) then

        lendatal = lendata
        lendata = -1         
        ierr = getfld(ivar, array)
        if (ierr .ne. 0) then 
          call errcheck(' ', 'subroutine getvdata', 0)
          return
        endif

c        allocate space for data and mask on original grid

        ierr1 = miallocw(lendatan, p1)
        ierr2 = miallocw(lendatan, p2)
        if ((ierr1 .ne. 0) .or. (ierr2 .ne. 0)) then
           print*, 'Error in attempting to allocate memory dynamically.'
           print*, 'You may need more memory.'
           print*, 'Try running on cirrus or retrieving smaller '
           print*, 'portions of the array.  Error encountered '
           print*, 'while processing', filnam(ivar), varnam(ivar)
           call errcheck('function miallocw','subroutine getvdata',1)
           return
        endif

        lendata = lendatan
        ierr = getfld(ivar, arrayin)
        if (ierr .ne. 0) then 
          call errcheck(' ', 'subroutine getvdata', 0)
          return
        endif

        do 80 i=1,lendata
          amaskin(i) = 1.
   80   continue 
             
        if (igeogin(ivar) .ne. 0) then
      
          masktype = masktypi
          ierr = xgeog(ivar, igeogin(ivar), ishowin(ivar), amaskin)
          if (ierr .ne. 0) then
             call errcheck(' ', 'subroutine getvdata',0)
             return
          endif
         
        endif
       
c        mask areas with missing data and set array to 0.

        if (misstype .eq. 1) then
          a1 = abs(1.0e-5*amissing)                                    
          do 90 i=1,lendata
            if (abs(arrayin(i)-amissing) .le. a1) amaskin(i) = 0. 
            if(amaskin(i) .eq. 0.) arrayin(i) = amissing                  
   90     continue 
        else
          do 95 i=1,lendata
            if (arrayin(i) .eq. amissing) amaskin(i) = 0.0
            if (amaskin(i) .eq. 0.0) arrayin(i) = amissing
   95     continue
        endif

c         prepare to regrid field:  set values in common  
c           if iregrid > 0  it will point to a variable in the
c                  table that has the target grid
c           if iregrid < 0 it will represent a spherical harmonic
c                   truncation (t21=-21, t42 = -42, etc.)
c                   or user-specified regridding options
c                  (either given in common or read in or set as default)
c
        if (iregrid(ivar) .lt. 0) then

c           get longitude-latitude range from source grid
                 
          ilono = 0
          ilato = 0
          do 98 i = 1,4
            ii = iptdim(i,ivar)
            iii = 0
            if (ii .ne. 0)  iii = univunit(dimnam(ii))
            if (iii .eq. intlong) ilono = ii
            if (iii .eq. intlat)  ilato = ii               
   98     continue

c ?/?
          domainc(iloni) = 'range'
          if (nlonrgdv(ivar) .eq. 0) then
            bdry(1,iloni) = begbdy(ilono)
            bdry(2,iloni) = endbdy(ilono)
          else
            bdry(1,iloni) = alonrgdv(ivar) - 0.5*dlonrgdv(ivar)
            bdry(2,iloni) = alonrgdv(ivar) + 
     &                      (nlonrgdv(ivar)-0.5)*dlonrgdv(ivar)
          endif

          domainc(ilati) = 'range'
          if ((nlatrgdv(ivar) .gt. 0) .and.
     &        (dlatrgdv(ivar) .ne. 0.0)) then
            bdry(1,ilati) = alatrgdv(ivar) - 0.5*dlatrgdv(ivar)
            bdry(2,ilati) = alatrgdv(ivar) +
     &                      (nlatrgdv(ivar)-0.5)*dlatrgdv(ivar)
          else
            bdry(1,ilati) = begbdy(ilato)
            bdry(2,ilati) = endbdy(ilato)
          endif

        endif

        lendata = lendatal                           
        ierr = doregrid(ivar, amaskin, arrayin, amask, array)             

        if (ierr .ne. 0) then
          call errcheck(' ', 'subroutine getvdata',0)
          return
        endif
          
        ierr = mifree(p1)
        ierr = mifree(p2)

      else
                
        ierr = getfld(ivar, array)
      
        if (ierr .ne. 0) then
          call errcheck(' ', 'subroutine getvdata',0)
          return
        endif
               
        len = ivarlen(ivar)
        do 100 i=1,len
            amask(i) = 1.
  100   continue

        if (igeogin(ivar) .ne. 0) then
      
          masktype = masktypi
          ierr = xgeog(ivar, igeogin(ivar), ishowin(ivar), amask)
          if (ierr .ne. 0) then
             call errcheck(' ', 'subroutine getvdata',0)
             return
          endif
         
        endif
       
c        mask areas with missing data and set array to 0.

        if (misstype .eq. 1) then
          a1 = abs(1.0e-5*amissing)
          do 250 i=1,len
             if (abs(array(i)-amissing) .le. a1) amask(i) = 0.
             if (amask(i) .eq. 0.) array(i) = amissing
  250     continue
        else
          do 255 i=1,len
            if (array(i) .eq. amissing) amask(i) = 0.0
            if (amask(i) .eq. 0.0) array(i) = amissing
  255     continue
        endif  

      endif

      if (igeogout(ivar) .ne. 0) then

c           allocate space for temporary mask that will be fraction of 
c            grid cell not masked

            jj = ivarlen(ivar)
            ierr1 = miallocw(jj, p1)
            do 260 j=1,jj
              amaskin(j) = 1.0
  260       continue

            masktype = masktypo
            ierr = xgeog(ivar, igeogout(ivar), ishowout(ivar), amaskin)
            if (ierr .ne. 0) then
                call errcheck(' ', 'subroutine getvdata',0)
                return
            endif

            do 270 j=1,jj
              amask(j) = amin1(amask(j), amaskin(j))
  270       continue

            ierr = mifree(p1)

      endif
             
      
c        apply grid-cell weights to current mask

         ierr = applywts(ivar,amask) 

c        mask areas with missing data and set array to omit

         jj = ivarlen(ivar)

         if (misstype .eq. 1) then
           a1 = abs(1.0e-5*amissing)
           do 1300 j=1,jj
             if (abs(array(j)-amissing) .le. a1) amask(j) = 0.
             if (amask(j) .eq. 0.)  array(j) = omit
 1300      continue
         else   
           do 1305 j=1,jj
             if (abs(array(j)) .eq. amissing) amask(j) = 0.0
             if (amask(j) .eq. 0.0) array(j) = omit
 1305      continue
         endif
 
      return
      end




c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     subroutine initget
c     -------
c
c     Description:
c     -----------
c     this function initializes common variables with default values
c
c     Usage:
c     ------
c
c      call initget
c
c     ------
c     Date: 1/19/92
c     ----
c
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine initget()
          
      implicit none
#include "drsdef.h"
#include "fcddrs.h"

#include "usersget"
#include "cscalars"
#include "cvartbl"
#include "cdimtbl"
#include "ciotbl"
#include "cdomain"
      integer cllun, seterr
      integer ierr, i, ii, j

c    intialize error controls
      warncnt = 0
      errcount = 0
      errcntl = 2
      
c    initialize various variables

      masktypi = 0
      masktypo = 0
      masktype = 0

      lendata = 0
      lendatan = 0

      if (mxsfctyp .gt. 26) then
          print*, 'mxsfctyp must be less than 27 on a 32-bit machine'
          call errcheck('subroutine initget', ' ',1)
          return
      endif
      mxdimelm = 20000
      mxdmkeep = 2000
      ifreecoo = 1
      icount = 0
      coordvls(lencoord) = 1.
      ierr = seterr(6,IDRS_NOREPORT)
      do 100 i = 1,maxiodev
        ii = 90 + i
        ierr = cllun(ii)
        iofilnam(i) = 'not defined'
        iofiluse(i) = i - 1 - maxiodev
 100  continue
      ierr = seterr(6,IDRS_FATAL)

      do 200 i = 1, maxvars
         varnam(i) = 'not defined'
         filnam(i) = 'not defined'
         vartitle(i) = 'not defined'
         varunits(i) = 'not defined'
         vartype(i) = 'not def.'
         varsourc(i) = 'not defined'
         vardate(i) = 'not def.'
         vartime(i) = 'not def.'

         ivarndim(i) = 0
         ivarlen(i) = 0
         igeogin(i) = 0
         igeogout(i) = 0
         ishowin(i) = all
         ishowout(i) = all
         iregrid(i) = 0

         nlatrgdv(i) = 0
         dlatrgdv(i) = 0.0
         alatrgdv(i) = 0.0
         nlonrgdv(i) = 0
         dlonrgdv(i) = 0.0
         alonrgdv(i) = 0.0  
         
         do 250  j=1,4
            iptdim(j,i) = 0
            iptdimrg(j,i) = 0
            vdimname(j,i) = 'not defined'
            vdimwt(j,i) = 'not def.'
            vdimdom(j,i) = 'not def.'
            vdimbdry(1,j,i) = 0.   
            vdimbdry(2,j,i) = 0.   
            ivdimpln(1,j,i) = 0.   
            ivdimpln(2,j,i) = 0.   
            vdimcycl(j,i) = 0.
  250    continue

  200 continue

      misstype = 1
      amissing = 1.e20
      omit = 1.e20
      ittrunc = 0


      do 300 i = 1,maxdims
         dimnam(i) = 'not defined'
         idimlen(i) = 0
         ioriglen(i) = 0
         idimtype(i) = 0
         iptcoord(i) = 0
         iptedge(i) = 0
         iptwts(i) = 0
         dimdoman(i) = 'not def.'
         dimcycle(i) = 0.0
         dimwtflg(i) = 'not def.'
         begbdy(i) = 0.0
         endbdy(i) = 0.0
         felretr(i) = 0.0
         lelretr(i) = 0.0
  300 continue

      namelong = 'Longitude'
      namelat = 'Latitude'
      
      do 400 i=1,4
        namedim(i) = 'not defined'
        bdry(1,i) = 0.0
        bdry(2,i) = 0.0
        ibdry(1,i) = 0
        ibdry(1,i) = 0
        domainc(i) = 'not def.'
        cycle(i) = 0.0
        cwt(i) = 'not def.'
  400 continue   

      return
      end

 
 
      
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     subroutine lendims
c     -------
c
c     Description:
c     -----------
c     this function returns the length of each dimension of variable 
c            ivar  (after the variable has been extracted by getdata
c            or getvdata or after calling shape.)
c
c     Usage:
c     ------
c
c      call lendims(ivar, i1, i2, i3, i4, ii)
c
c         ivar = variable i.d.
c         i1, i2, i3, i4 = length of dimension 1 through 4, respectively
c         isize = i1*i2*i3*i4 (after replacing any 0's with 1's) 
c     ------
c     Date: 1/19/92
c     ----
c
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine lendims(ivar, i1, i2, i3, i4, ii)
           
      implicit none
#include "usersget"
#include "cvartbl"
#include "cdimtbl"
#include "cscalars"

      integer ivar, i1, i2, i3, i4, ii, i, n, j, nn(4), univunit
      logical caseindp 

      j = 1                                            
      do 100 n=1,4
        i = iptdimrg(j,ivar)
        if (i .eq. 0) then
          nn(n) = 0
          j = j + 1 
        elseif (caseindp(vdimname(n,ivar), dimnam(i)) .or.
     &      ((univunit(dimnam(i)) .gt. 0) .and.
     &      (univunit(dimnam(i)) .eq. univunit(vdimname(n,ivar))))) then
          nn(n) =  idimlen(i)
          j = j + 1
        else
          nn(n) = 0
        endif
  100 continue

      i1 = nn(1)
      i2 = nn(2)
      i3 = nn(3)
      i4 = nn(4)
      
      ii = max0(1,i1)*max0(1,i2)*max0(1,i3)*max0(1,i4)
      
      if ((i1+i2+i3+i4) .eq. 0) then
         if (errcntl .gt. 0) then
         print*, 'Warning:  the length of variable ', ivar 
         print*, 'appears to be 0.  It is possible that you have'
         print*, 'cleared the dimension table, so that information'
         print*, 'on the dimensions is no longer available.'
         endif
         call errcheck('subroutine lendims', ' ',-1)
      endif
      
      return
      end
      
  
      subroutine maparea(nloni, nlono, nlati, nlato,
     &    bnin, bnout, bsin, bsout, bein, beout, bwin, bwout,
     &    londx, lonpt, wtlon, latdx, latpt, wtlat)
     
     
c  *********************************************************************
c
c  This subroutine provides the information needed to interpolate
c  from the input grid to the output grid (i.e., it calculates
c  the weight of the contribution from each of the input cells
c  to each of the output cells.
c
c  The input and output grid-cell boundaries must be specified by the
c  user(bnin, bsin, bein, bwin, bnout, bsout, beout, and bwout), and
c  the ordering must be monotonic, but need not be continuous (input or 
c  output cells may be missing). The units are degrees latitude and 
c  longitude.  
c  
c   The input and output domains do not need to be identical, and
c   the longitude coordinate will be "wrapped around" if necessary.
c   The coordinate ordering may differ between input and output grids
c   (for example N to S on input, S to N on output.
c
c
c   Input:
c
c    nloni = number of input grid cell longitudes.
c    nlono = number of output grid cells longitudes.
c    nlati = number of input grid cell latitudes.
c    nlato = number of output grid cells latitudes.
c    bnin(nlati) = northern boundary of each grid cell of input field.
c    bsin(nlati) = southern boundary of each grid cell of input field.
c    bein(nloni) = eastern boundary of each grid cell of input field.
c    bwin(nloni) = western boundary of each grid cell of input field.
c    bnout(nlato) =northern boundary of each grid cell of output field.
c    bsout(nlato) =southern boundary of each grid cell of output field.
c    beout(nlono) = eastern boundary of each grid cell of output field.
c    bwout(nlono) = western boundary of each grid cell of output field.
c     
c   Output
c
c    londx(nlono) = index of the last element in lonpt and wtlon that
c                   apply to each of the the output longitudes
c    lonpt(nloni+nlono) = array of indices pointing to input grid cells
c                         that contribute to each output grid cell.
c    wtlon(nloni+nlono) = array of weights indicating how much each
c                         input grid cell contributes to each 
c                         output grid cell.
c
c    latdx(nlato) = index of the last element in latpt and wtlat that
c                   apply to each of the the output latitudes
c    latpt(nlati+nlato) = array of indices pointing to input grid cells
c                         that contribute to each output grid cell.
c    wtlat(nlati+nlato) = array of weights indicating how much each
c                         input grid cell contributes to each 
c                         output grid cell.
c
c  *********************************************************************

     
      implicit none
      integer nloni, nlono, nlati, nlato, londx(nlono), lonpt(*), 
     &       latdx(nlato), latpt(*) 
      real  bnin(nlati), bnout(nlato), bsin(nlati), bsout(nlato)
      real   bein(nloni), beout(nlono), bwin(nloni), bwout(nlono),
     &       wtlon(*), wtlat(*)
      real pi, bnorth, bsouth
      integer istradle, j1, j2, j3, last, jo, iflag, ji, ii, i, ip,
     &     isize, ierr1, ierr2, ierr
      integer miallocw, mifree
      pointer (p1, bwinl), (p2, beinl)
      real westout, eastout, bwinl(*), beinl(*)
      
      pi = 4.*atan(1.0)
      
c    cycle through output zones

      if ((bsin(nlati)-bsin(1)) * (bsout(nlato)-bsout(1)) .gt. 0.) then
        j1=1
        j2=nlati
        j3=1
      else
        j1=nlati
        j2=1
        j3=-1
      endif

      last = 0
      do 200 jo = 1, nlato
        
        iflag = 0
c       find index of input latitude zones that will contribute.
        do 100 ji = j1,j2,j3
          if ((bsin(ji) .lt. bnout(jo)) .and.
     &        (bnin(ji) .gt. bsout(jo))) then 
            iflag = 1
            last = last + 1
            latpt(last) = ji
            bnorth = amin1(bnout(jo),bnin(ji))
            bsouth = amax1(bsout(jo),bsin(ji))
            wtlat(last) = sin(bnorth*pi/180.) - sin(bsouth*pi/180.)
          elseif (iflag .eq. 1) then
             goto 150
          endif
  100   continue
  
  150   latdx(jo) = last
        if (last .ne. 0) j1 = latpt(last)
  
  200 continue
     
            

      westout = amin1(bwout(1), bwout(nlono))
      eastout = amax1(beout(1), beout(nlono))

c    find input grid cell that straddles the output domain boundary.

      if (bwin(1) .lt. bwin(nloni)) then
        ii = nint((westout-bwin(nloni))/360.+0.5) + 1
      else
        ii = nint((westout-bwin(nloni))/360.+0.5) - 1
      endif
      do 225 i=1,nloni
        ip = ii
        ii = nint((westout-bwin(i))/360.+0.5)
        if (ii .ne. ip) then
             if (bwin(1) .lt. bwin(nloni)) then
               if (westout .eq. bwin(i)) then
                 istradle=i
               else
                 istradle=i-1
                 if (istradle .lt. 1) istradle=nloni
               endif
             else
               istradle=i
             endif
             go to 250
        endif
  225 continue
      if (bwin(1) .lt. bwin(nloni)) then
        istradle=nloni
      else
        istradle=1
      endif
  250 continue

c     allocate enough space for one more than the number of longitude
c       grid cells and make a copy of longitudes that coincides
c       with output domain.

      isize = nloni+1
      ierr1 = miallocw(isize,p1)
      ierr2 = miallocw(isize,p2)
      if ((ierr1 .ne. 0) .or. (ierr2 .ne. 0)) then
           print*, 'Error in attempting to allocate memory dynamically.'
           print*, 'You may need more memory.'
           print*, 'Try running on cirrus or retrieving smaller '
           print*, 'portions of the array.'
           call errcheck('function miallocw','subroutine maparea',1)
           return
      endif

c      isize = (nloni+1)*4
c      p1 = malloc(isize)
c      p2 = malloc(isize)


      do 300 ii=1,nloni+1

         i = istradle + ii - 1
         if (i .gt. nloni) i = i - nloni
     
c        make western bdry of input grid cell fall within output 
c           domain.
         bwinl(ii) = bwin(i) + 360.*nint((westout-bwin(i))/360.+0.5) 
         beinl(ii) = bein(i) + 360.*nint((westout-bwin(i))/360.+0.5)

  300 continue

      if (bwin(nloni) .gt. bwin(1)) then

        if (bwinl(1) .gt. westout) then
           bwinl(1) = bwinl(1) - 360.0
           beinl(1) = beinl(1) - 360.0
        endif

      else

        if (bwinl(nloni+1) .gt. westout) then
           bwinl(nloni+1) = bwinl(nloni+1) - 360.0
           beinl(nloni+1) = beinl(nloni+1) - 360.0
        endif

      endif

c    cycle through output longitudes

      if ((bwinl(nloni+1)-bwinl(1)) * (bwout(nlono)-bwout(1)) .gt. 0.) 
     +       then
        j1=1
        j2=nloni+1
        j3=1
      else
        j1=nloni+1
        j2=1
        j3=-1
      endif

      last = 0
      do 500 jo = 1, nlono
        
        iflag = 0
c       find index of input longitude zones that will contribute.
        do 400 ji = j1,j2,j3
          if ((bwinl(ji) .lt. beout(jo)) .and.
     &        (beinl(ji) .gt. bwout(jo))) then 
            iflag = 1          
            last = last + 1
            lonpt(last) = ji + istradle - 1
            if (lonpt(last) .gt. nloni) lonpt(last) = lonpt(last)-nloni
            wtlon(last) = amin1(beout(jo),beinl(ji)) - 
     +                    amax1(bwout(jo),bwinl(ji))
          elseif (iflag .eq. 1) then
             goto 450
          endif
  400   continue
  
  450   londx(jo) = last
        if (last .ne. 0) j1 = lonpt(last) - istradle + 1
        if (j1 .lt. 1) j1 = j1 + nloni 
  
  500 continue
     
      ierr = mifree(p1)
      ierr = mifree(p2)
c      call free(p1)
c      call free(p2)
      
      return
      end
  



                       
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     integer function mkdimtbl
c     -------
c
c     ------
c     Date: 1/19/92
c     ----
c
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
      integer function mkdimtbl (mxdmkeep, istore, dname, lendim, idlen,
     +      dsource, dtitle, dunits, idtype, ddomainc, cycle, bdry1,
     &      bdry2, felrloc, lelrloc, irow, cwt, vecdim, tedge, twt)

c      common variables
          
      implicit none
#include "usersget"
#include "cdimtbl"

      integer mxdmkeep, istore
      character*(*) dname, ddomainc, cwt, dsource, dtitle, dunits 
      character*8 domaincl
      integer lendim, idlen, idtype, irow, idum, i, j, k, m, j3
      real   cycle, bdry1, bdry2, felrloc, lelrloc, vecdim(*)
      real v, t, e
      double precision tedge(*), twt(*)
      logical caseindp

      mkdimtbl = 0

c     check whether dimension already exists in table
      if (ddomainc(1:7) .eq. 'assaved') ddomainc = 'range'
      do 100 i = 1, icount
        domaincl = dimdoman(i)
        if ((domaincl(1:7) .eq. 'assaved')) domaincl = 'range'
        if (caseindp(dname, dimnam(i)) .and.
     +   (lendim .eq. idimlen(i)) .and. (idlen .eq. ioriglen(i)) .and.
     +   (ddomainc .eq. domaincl) .and. (cycle .eq. dimcycle(i))
     +   .and. (bdry1 .eq. begbdy(i)) .and. (bdry2 .eq. endbdy(i)) .and.
     +  (felrloc .eq. felretr(i)) .and. (lelrloc .eq. lelretr(i)) .and.
     +  (cwt .eq. dimwtflg(i)) .and. (idtype .eq. idimtype(i))) then 

          if ((lendim .le. mxdmkeep) .and. (istore .gt. 0)) then
            j = iptcoord(i)
            k = iptwts(i)
            m = iptedge(i)
            if ((j .gt. 0) .and. (k .gt. 0) .and. (m .gt. 0)) then
              j3 = 1+(lendim/3)
              v = vecdim(j3)
              t = twt(j3)
              e = tedge(j3+1)
              if ((abs(v-coordvls(j+j3-1)) .le. abs(v)*1.e-5) .and. 
     &            (abs(t-coordvls(k+j3-1)) .le. abs(t)*1.e-5) .and.
     &            (abs(e-coordvls(m+j3)) .le. abs(e)*1.e-5)) then
                irow = i
                return
              endif
            endif
          else  
            irow = i
            return
          endif
        endif
  100 continue

c     check whether array will fit in coordvls
      if ((istore .gt. 0) .and. (lendim .le. mxdmkeep) .and. 
     &    (3*lendim+1+ifreecoo-1 .ge. lencoord)) then
          idum = 3*lendim+1+ifreecoo-1
          print*, 'array coordvls is too small.'
          print*, 'When creating the EZGET library, you should set'
          print*, 'lencoord in common usersget'
          print*, 'to a value no smaller than ', idum  
          print*, 'Alternatively, you might call clrtable at a'
          print*, 'convenient point to make room for new dimensions'
          print*, 'if the old ones are no longer needed.'
          mkdimtbl = 1
          call errcheck('function mkdimtbl', ' ',1)
          return
      endif

      icount = icount + 1
      irow = icount

c     check whether another row will fit in dimension table
      if (irow .gt. maxdims) then
         print*, 'parameter maxdims is too small'
         print*, 'In common usedsdrs, you should set '
         print*, 'maxdims to a value no smaller than ', irow
         print*, 'Alternatively, you may clear the dimension table by'
         print*, 'calling subroutine clrtable, but you should only do'
         print*, 'this if you do not plan to retrieve dimension'
         print*, 'information from currently defined variables (unless'
         print*, 'you again obtain data from the defined variable' 
         print*, 'which will redefine its dimensions).'
         mkdimtbl = 1
         call errcheck('function mkdimtbl', ' ',1)
         return
      endif

      dimnam(irow) = dname
      idimlen(irow) = lendim
      ioriglen(irow) = idlen
      idimtype(irow) = idtype
      dimsourc(irow) = dsource
      dimtitle(irow) = dtitle
      dimunits(irow) = dunits
      dimdoman(irow) = ddomainc
      dimcycle(irow) = cycle
      begbdy(irow) = bdry1
      endbdy(irow) = bdry2
      felretr(irow) = felrloc
      lelretr(irow) = lelrloc
      dimwtflg(irow) = cwt

      if ((istore .eq. 0) .or. (lendim .gt. mxdmkeep)) then

        if (caseindp(cwt, 'unit')) then
          iptcoord(irow) = -1
          iptwts(irow) = -1
          iptedge(irow) = -1
        else
          print*, 'A dimension named ', dname, 'is too long to save' 
          print*, 'under the default limits.'
          print*, 'If: 1) you will not be actually retrieving '
          print*, 'the dimension by calling getcoord, getedges, or '
          print*, 'getdimwt, and 2) this dimension is not a longitude'
          print*, 'or latitude dimension required when mapping data'
          print*, 'to a new grid or for geography masks,  and '
          print*, '3) unit weights can be assigned to this dimension,'
          print*, 'then you may successfully execute by assigning '
          print*, 'unit weighting to this dimension.'
          print*, 'Alternatively, you may increase the limit on '
          print*, 'the dimension size by calling subroutine defmisc'
          print*, 'with the "maximum dimension kept" option and its '
          print*, 'value set larger than ', lendim
          print*, ' '
          call errcheck('function mkdimtbl',' ',2)
          mkdimtbl = 1
          return
        endif

      else

        iptcoord(irow) = ifreecoo
        iptwts(irow)  = ifreecoo+lendim
        iptedge(irow) = ifreecoo+2*lendim

        do 200 i = 1,lendim
          coordvls(ifreecoo-1+i) = vecdim(i)
          coordvls(ifreecoo+lendim-1+i) = twt(i)
          coordvls(ifreecoo+2*lendim-1+i) = tedge(i)
  200   continue
        coordvls(ifreecoo+3*lendim) = tedge(lendim+1)

        ifreecoo = ifreecoo + 3*lendim + 1

      endif
      
      return
      end



c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     integer function opendrs
c     -------
c
c     Description:
c     -----------
c     This function opens a DRS file, checks whether the variable sought
c       is actually in the file, and puts information into the i.o. unit
c       table.
c
c
c     Usage:
c     ------
c     ierr = opendrs(variable#, error-cntrl, i.o.-unit#)
c
c
c
c
c     Exit Status:
c     -----------
c     0 if success
c     nonzero if an error occurred
c
c     ------
c     Date: 1/19/92
c     ----
c
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer function opendrs(ivar, ierrcntl, idicun)


      implicit none
#include "drsdef.h"
#include "fcddrs.h"

c      common variables

#include "usersget"
#include "ciotbl"
#include "cvartbl"

c     local variables

      character dumnam*16, vtitle*80, vsourc*120, vunits*40, vtype*8
      data dumnam /"                "/

      integer idicun, idatun, iuse, lastuse, nextunit, ierr, ivar,
     &     ierrcntl, i, icount
cjfp was      integer inqdict, aslun, setname, seterr,
cjfp was     +      cluvdb, cllun, getname, drstest
      integer fcw_inqdict, fcw_aslun, fcw_setname, fcw_seterr
      integer fcw_cluvdb, fcw_cllun, fcw_getname
      logical fcw_drstest
      logical  caseindp

      opendrs = 0
c
c    if ivar .le. 0 then close all open units and return
      
      if (ivar .le. 0) then
         do 30 i = 1,maxiodev
            if (iofiluse(i) .ge. 0) ierr = fcw_cllun(90+i)
            iofiluse(i) = i - 1 - maxiodev
   30    continue
         return
      endif

c    search file table to see if file already opened
c    if so, return unit number, if not, open file and fill file table
c   if table is full, close oldest file and replace it with current file

c      check whether user remembered to define file path
      if (caseindp(filnam(ivar), 'not def')) then
         print*, 'You must define filnam for variable, ', varnam(ivar)
         opendrs = 1
         call errcheck('function opendrs', ' ',1)
         return
      endif
      
      lastuse = 0
      idicun = 0
      do 40 i = 1,maxiodev
           lastuse = max0(lastuse,iofiluse(i))
           if (filnam(ivar) .eq. iofilnam(i)) then
              nextunit = i
              idicun = 90 + i
              idatun = 95 + i
           endif
   40 continue
   
      if (idicun .ne. 0) then
      
         iofiluse(nextunit) = lastuse + 1
         
      else
      
         iuse = 100000

         do 50 i = 1,maxiodev
           if (iofiluse(i) .lt. iuse) then
              nextunit = i
              iuse = iofiluse(i)
           endif
   50    continue

         idicun = nextunit + 90
         idatun = nextunit + 95
         iofilnam(nextunit) = filnam(ivar)

         if (iofiluse(nextunit) .gt. 0) ierr = fcw_cllun(idicun)

      iofiluse(nextunit) = lastuse + 1

         icount = 0
c        Open the input DRS file. Exit if an error occurs.
c   60    if (cw_drstest(cw_aslun(idicun, filnam(ivar), idatun, ' ',
c     +        IDRS_READ)) .ne. 0 ) then
 60      continue
         ierr=fcw_aslun(idicun, filnam(ivar), idatun, ' ',IDRS_READ)
         if (fcw_drstest(ierr)) then
            icount = icount + 1
            if (icount .lt. 3) go to 60
            opendrs = -1000
            if (ierrcntl .gt. 0) then
            print*, 'error opening file --', filnam(ivar)
            print*,'check that you have properly defined the path/name'
            call errcheck('DRS function aslun', 'function opendrs',1)
            endif
            return
         endif

      endif

c     Determine whether the variable is in the file. If not, return.
c           else put info. in vartable

      if (caseindp(vartitle(ivar), 'not def')) then 
          vtitle = ' '
      else
          vtitle = vartitle(ivar)
      endif

      if (caseindp(varsourc(ivar), 'not def')) then 
          vsourc = ' '
      else
          vsourc = varsourc(ivar)
      endif

      if (caseindp(varunits(ivar), 'not def')) then 
          vunits = ' '
      else
          vunits = varunits(ivar)
      endif

      vtype = ' '

      ierr = fcw_cluvdb()
      if (ierr .ne. 0) then
         opendrs = 1
         call errcheck('DRS function fcw_cluvdb', 'function opendrs',1)
         return
      endif
      ierr = fcw_setname(vsourc,varnam(ivar),vtitle,vunits,vtype)
      if (ierr .ne. 0) then
         opendrs = 1
         call errcheck('DRS function fcw_setname', 'function opendrs',1)
         return
      endif
      ierr = fcw_seterr(6,IDRS_FATAL)
cjfp was      ierr = inqdict(idicun,IDRS_GETFIRSTVAR)
      ierr = fcw_inqdict(idicun,IDRS_GETFIRSTVAR)

      if(ierr.ne.IDRS_SUCCESS) then
         opendrs = -1001
         print *,'Variable ',varnam(ivar),' was not found in file ',
     +                     filnam(ivar)
         if (.not. caseindp(vartitle(ivar), 'not def'))  
     +      print*, '  Var. Title: ', vartitle(ivar)
         if (.not. caseindp(varsourc(ivar), 'not def')) 
     +      print*, '  Var. Source: ', varsourc(ivar) 
         if (.not. caseindp(varunits(ivar), 'not def'))
     +      print*, '  Var. Units: ', varunits(ivar)  
         if (ierrcntl .gt. 0) then
         call errcheck('function opendrs', ' ',1)
         endif
         return
      endif
      ierr = fcw_seterr(6,IDRS_WARNING)

c         retrieve number of dimensions (and other info)
      ierr = fcw_getname(
     +     varsourc(ivar), dumnam, vartitle(ivar),
     +     varunits(ivar), vardate(ivar), vartime(ivar),
     +     vartype(ivar), ivarndim(ivar))
      if (ierr .ne. 0) then
         opendrs = 1
         call errcheck('DRS function fcw_getname', 'function opendrs',1)
         return
      endif

      return
      end
       
      subroutine rgdarea(ilong, ilat, itim1, itim2, 
     &    ntim1, ntim2, nloni, nlono, nlati,nlato,
     &    londx, lonpt, wtlon, latdx, latpt, wtlat, omit,
     &    amskin, ain, amskout, aout)
     
     
c  *********************************************************************
c
c  This subroutine takes an input field (fldin) of grid cells and fills
c  an output field (fldout) of grid cells with area-weighted values,
c  thereby achieving a transformation from one grid to another.
c
c  The program uses pointers and weights generated outside this
c   routine to determine how much each input grid cell contributes
c   to each output grid cell
c
c   Besides the automatic area weighting, an additional user-supplied
c   weighting factor (amskin) must be passed to this subroutine, which
c   multiplies the area weighting.  If simple area-weighting is needed,
c   supply a vector of 1.''s in the input.  However, if for example only
c   land values should contribute to the new gridded field, then 1.''s
c   should be assigned for land points, and 0.''s for ocean points.
c
c
c   Input:
c
c    ilong = dimension position for longitude in input and output arrays
c          = 1 if first dimension
c          = 2 if second dimension
c          = 3 if third dimension
c          = 4 if fourth dimension
c    ilat = dimension position for latitude in input and output arrays
c    itim1 = dimension position for time or level in input and output
c              arrays (if no itim1 dimension, then this parameter is 
c              ignored) 
c    itim2 = dimension position for time or level in input and output 
c              arrays (if no itim2 dimension, then this parameter is
c              ignored) 
c    nloni = number of input grid cell longitudes.
c    nlono = number of output grid cells longitudes.
c    nlati = number of input grid cell latitudes.
c    nlato = number of output grid cells latitudes.
c
c    ntim1 = length of itim1 dimension (if this is a dummy dimension, 
c              you should set ntim1 to 0)
c    ntim2 = length of itim2 dimension (if this is a dummy dimension, 
c              you should set ntim2 to 0)
c
c    londx(nlono) = index of the last element in lonpt and wtlon that
c                   apply to each of the the output longitudes
c    lonpt(nloni+nlono) = array of indices pointing to input grid cells
c                         that contribute to each output grid cell.
c    wtlon(nloni+nlono) = array of weights indicating how much each
c                         to input grid cells that contribute to each 
c                         output grid cell.
c
c    latdx(nlato) = index of the last element in latpt and wtlat that
c                   apply to each of the the output latitudes
c    latpt(nlati+nlato) = array of indices pointing to input grid cells
c                         that contribute to each output grid cell.
c    wtlat(nlati+nlato) = array of weights indicating how much each
c                         to input grid cells that contribute to each 
c                         output grid cell.
c    omit = value that will be assigned grid points that have either
c       been masked out or cannot be computed for lack of input data.
c
c    amskin(*) = input weighting factors (for masking purposes).
c    ain(*) = input field
c
c    Input and Output:
c
c    amskout(*) = output weighting factors
c    aout(*) = output (regridded) field
c     
c
c  *********************************************************************
c

      implicit none
      pointer (p1, accum), (p2, wtmsk)
      integer  nloni, nlono, nlati, nlato
      integer londx(nlono), lonpt(*), latdx(nlato), latpt(*)
      real wtlon(*), wtlat(*), amskin(*), ain(*), amskout(*), aout(*)
      real omit
      
      integer ilat, ilong, ii, itim1, itim2, nn, ntim1, ntim2,
     +          je, jo, jb, ie, io, ib, k, ji, k1, k2, kk,
     +          mtim1, mtim2, jtim1, jtim2, i, j, imiss1, imiss2, 
     &          isize, ierr1, ierr2, ierr
      integer leni(4), leno(4), iai(4), iao(4), id(4), nd(4)
      integer miallocw, mifree
            
      double precision accum(*), wtmsk(*), wt
      

      id(1) = ilong
      id(2) = ilat
      id(3) = itim1
      id(4) = itim2

      nd(1) = nloni
      nd(2) = nlati
      nd(3) = ntim1
      nd(4) = ntim2

      do 80 i=3,4
        if (nd(i) .eq. 0) id(i) = 0
   80 continue 

      imiss1=0
      imiss2=0

      do 105 k=1,4       
        do 100 i=1,4   
          if (id(i) .eq. k) then
            if (k .lt. 4) then
c             look for duplicates
              do 90 j=i+1,4
                if (id(j) .eq. k) then
           print*, 'Error in specifying data structure while attempting'
           print*, 'to regrid (horizontally interpolate) data.  You '
           print*, 'have set:'
           print*, 'ilong = ', ilong
           print*, 'ilat =  ',  ilat
           print*, 'itim1 = ', itim1
           print*, 'itim2 = ', itim2
           print*, 'One of these should = 1, another = 2, another = 3, '
           print*, 'and another = 4'
           print*, 'If itim1 and/or itim2 are dummy dimensions, then'
           print*, 'be sure that ntim1 and/or ntim2 are 0, so that'
           print*, 'itim1 and/or itim2 will be ignored.' 
           call errcheck('subroutine rgdarea', ' ', 1)
           return
                endif
   90         continue
            endif
            go to 105
          endif 
  100   continue
        if (imiss1 .eq. 0) then
          imiss1 = k 
        else
          imiss2 = k 
        endif
  105 continue


      if (ntim1 .gt. 0) then
        jtim1 = itim1
        mtim1 = ntim1                                       
      else
        mtim1 = 1
        jtim1 = imiss1
      endif

      if (ntim2 .gt. 0) then
        jtim2 = itim2
        mtim2 = ntim2                                       
      elseif (imiss2 .gt. 0) then
        jtim2 = imiss2
        mtim2 = 1
      else 
        jtim2 = imiss1
        mtim2 = 1
      endif

      if (ilong+ilat+jtim1+jtim2 .ne. 10) then
        print*, 'Error in specifying data structure while attempting'
        print*, 'to regrid (horizontally interpolate) data.  You '
        print*, 'have set:'
        print*, 'ilong = ', ilong
        print*, 'ilat =  ',  ilat
        print*, 'itim1 = ', itim1
        print*, 'itim2 = ', itim2
        print*, 'One of these should = 1, another = 2, another = 3, '
        print*, 'and another = 4'
        print*, 'If itim1 and/or itim2 are dummy dimensions, then'
        print*, 'be sure that ntim1 and/or ntim2 are 0' 
        call errcheck('subroutine rgdarea', ' ', 1)
        return
      endif

      leni(ilong) = nloni
      leni(ilat) = nlati
      leni(jtim1) = mtim1 
      leni(jtim2) = mtim2 
    
      leno(ilong) = nlono
      leno(ilat) = nlato
      leno(jtim1) = mtim1 
      leno(jtim2) = mtim2 
                 
c     allocate enough space for double precision variables
      isize = mtim1*mtim2*2
      ierr1 = miallocw(isize,p1)
      ierr2 = miallocw(isize,p2)
      if ((ierr1 .ne. 0) .or. (ierr2 .ne. 0)) then
           print*, 'Error in attempting to allocate memory dynamically.'
           print*, 'You may need more memory.'
           print*, 'Try running on cirrus or retrieving smaller '
           print*, 'portions of the array.'
           call errcheck('function miallocw','subroutine rgdarea',1)
           return
      endif

c    cycle through output zones and longitudes

      je = 0
      
      do 2000 jo = 1, nlato
      
        iao(ilat) = jo
        jb = je + 1
        je = latdx(jo)
        ie = 0

        do 1000 io = 1,nlono
        
          iao(ilong) = io
          ib = ie + 1
          ie = londx(io)
          
          do 150 k=1,mtim1*mtim2
              wtmsk(k) = 0.
              accum(k) = 0.
  150     continue
   
          if ((je .ge. jb) .and. (ie .ge. ib)) then
          
            wt = 0.
          
c           cycle through input zones and longitudes
            do 500 ji = jb,je
            
              iai(ilat) = latpt(ji)
              
              do 400 ii = ib,ie
              
                iai(ilong) = lonpt(ii)
                wt = wt +wtlat(ji)*wtlon(ii)
                
                do 300 k2=1,mtim2
                
                   iai(jtim2) = k2
                   
                   do 200 k1=1,mtim1
                   
                      iai(jtim1) = k1
                   
                      kk = (k2-1)*mtim1 + k1
                      nn = iai(1) + leni(1) * (iai(2)-1+leni(2) * 
     +                               (iai(3)-1+leni(3) * (iai(4)-1)))  
     
                      
                      wtmsk(kk) = wtmsk(kk) + wtlat(ji)*wtlon(ii) *
     +                                                       amskin(nn)
                      accum(kk) = accum(kk) + wtlat(ji)*wtlon(ii) *
     +                                            amskin(nn)* ain(nn)
  200              continue
  300           continue

  400         continue
  500       continue
  
            do 700 k2=1,mtim2
                
               iao(jtim2) = k2
                   
               do 600 k1=1,mtim1
                   
                  iao(jtim1) = k1
                   
                  kk = (k2-1)*mtim1 + k1
                  nn = iao(1) + leno(1) * (iao(2)-1+leno(2) * 
     +                               (iao(3)-1+leno(3) * (iao(4)-1)))  
          
            if (wtmsk(kk) .gt. 0.) then
              aout(nn) = accum(kk)/wtmsk(kk)
              amskout(nn) = wtmsk(kk)/wt
            else
              aout(nn) = omit
              amskout(nn) = 0.
            endif

  600              continue
  700           continue

  
          else
c           no input grid cells will contribute
            do 900 k2=1,mtim2
                
               iao(jtim2) = k2
                   
               do 800 k1=1,mtim1
                   
                  iao(jtim1) = k1
                   
                  nn = iao(1) + leno(1) * (iao(2)-1+leno(2) * 
     +                               (iao(3)-1+leno(3) * (iao(4)-1)))  
                  aout(nn) = omit
                  amskout(nn) = 0.
                  
  800          continue
  900       continue
            
          endif
  
 1000   continue
 2000 continue
 
      ierr = mifree(p1)
      ierr = mifree(p2)

      return
      end

      
      

c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     subroutine  shape
c     -------
c
c     Description:
c     -----------
c    This subroutine gets the size of each of the dimensions for
c        the requested domain of the variable ivar, but it 
c        doesn''t obtain the data or create a mask
c
c        If the dimensions for the variable have not been defined
c           by calls to defdim or defdimi, then the shape of the 
c           data is returned as it was saved.
c
c
c     Usage:
c     ------
c
c      call shape(ivar, idim1, idim2, idim3, idim4, lentot)
c                                            
c         ivar = variable i.d.
c         i1, i2, i3, i4 = length of dimension 1 through 4, respectively
c         isize = i1*i2*i3*i4 (after replacing any 0's with 1's) 
c
c     ------
c     Date: 1/24/93
c     ----
c
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

      subroutine shape(ivar, idim1, idim2, idim3, idim4, lentot)

      implicit none
#include "cscalars" 

      integer ivar, idim1, idim2, idim3, idim4, lentot, lendatal, 
     &    lenarray  
      real array(1), amask(1)
             
      lendatal = lendata
      lendata = -1
      
      call getvdata(ivar, amask, array)
      
      call lendims(ivar, idim1, idim2, idim3, idim4, lenarray)

      lentot = lendatan
      
      lendata = lendatal 
 
      return
      end




c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     integer functions univunit
c     -------
c
c     Description:
c     -----------
c      These functions are called when testing whether a dimension
c        is latitude, longitude, pressure or time so that the
c         test will be insensitive to case (upper or lower) and
c         the exact name of the variable
c
c
c     Exit Status:
c     -----------
c     0 if success
c     1 if an error occurred
c
c     ------
c     Date: 1/19/92
c     ----
c
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer function univunit(namedim)

      implicit none
#include "cscalars"

      character*16 namedim
      logical caseindp

      univunit = 0

      if (caseindp(namedim, namelong) .or.
     +    caseindp(namedim, 'lon')) then
                                      univunit = intlong

      elseif (caseindp(namedim, namelat) .or.
     +       caseindp(namedim, 'lat')) then
                                      univunit = intlat

      endif

      return
      end



c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     subroutine varinfo
c     -------
c
c     Description:
c     -----------
c     this subroutine returns information about a defined variable
c
c     Usage:
c     ------
c
c      call varinfo
c
c     ------
c     Date: 1/17/96
c     ----
c
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine varinfo(ivar, param, info)

      implicit none
#include "usersget"
#include "cvartbl"
#include "cdimtbl"

      integer ivar
      character*(*) param, info(*)
      integer i, j, n, iflag
      logical caseindp

      if (ivar .le. 0) then
             print*, 'You have improperly passed ivar = ', ivar
             print*, 'as the argument that identifies the variable'
             print*, 'you are requesting information from by calling'
             print*, 'subroutine varinfo.  This integer should be' 
             print*, '.gt. 0'  
             call errcheck('subroutine varinfo', ' ',1)
             return
      endif

      if (ivar .gt. maxvars) then
             
             print*, 'You have improperly passed ivar = ', ivar
             print*, 'as the argument that identifies the variable'
             print*, 'you are requesting information from by calling'
             print*, 'subroutine varinfo.  This integer should be' 
             print*, 'no larger than ', maxvars
             print*, 'unless you redefine maxvars in common usersget'
             call errcheck('subroutine varinfo', ' ',1)
             return
                  
       endif     

      n = len(info(1))
      iflag = 0

      if (caseindp(param, 'unit')) then

        if (n .lt. 40) then
          iflag = 40
        else
          info(1) = varunits(ivar)
        endif

      elseif (caseindp(param, 'sourc')) then

        if (n .lt. 120) then
          iflag = 120
        else
          info(1) = varsourc(ivar)
        endif

      elseif (caseindp(param, 'titl')) then

        if (n .lt. 80) then
          iflag = 80
        else
          info(1) = vartitle(ivar)
        endif

      elseif (caseindp(param, 'date')) then

        if (n .lt. 8) then
          iflag = 8
        else
          info(1) = vardate(ivar)
        endif

      elseif (caseindp(param, 'time')) then

        if (n .lt. 8) then
          iflag = 8
        else
          info(1) = vartime(ivar)
        endif

      elseif (caseindp(param, 'type')) then

        if (n .lt. 8) then
          iflag = 8
        else
          info(1) = vartype(ivar)
        endif

      elseif (caseindp(param, 'weight')) then

        if (n .lt. 8) then
          iflag = 8
        else
          do 100 j=1,4
            i = iptdimrg(j,ivar)
            if (i .ne. 0) info(j) = dimwtflg(i)
  100     continue
        endif

      else

        print*, 'error ...  You have incorrectly passed the following'
        print*, 'string as the first parameter in your call to varinfo:'
        print*, param
        print*, 'Only the following strings are acceptible choices:'
        print*, '"units", "source", "title", "date", "time", "type"'
        call errcheck('subroutine varinfo', ' ',1) 
        return

      endif

      if (iflag .ne. 0) then
        print*, 'Error in retrieving variable information.  The length'
        print*, 'of the character string you pass to varinfo when'
        print*, 'retrieving the ', param, ' should be no shorter than ',
     &           iflag
        print*, 'Your string was only ', n, ' characters long.'
        call errcheck('subroutine varinfo', ' ',1) 
        return
      endif


      return
      end


c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      function xgeog
c     -------
c
c     Description:
c     -----------
c
c      This function obtains a geography mask for a variable
c      filling gmask with 0. and 1. depending on whether the
c      point should be included or not.
c
c     Usage:
c     ------
c
c      xgeog (ivar, igeog, ishow, gmask)
c
c      name of variable is typically 'sfctype'
c
c     Exit Status:
c     -----------
c     0 if success
c     1 if an error occurred
c
c     ------
c     Date: 1/19/92
c     ----
c
c
c
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



      integer function xgeog(ivar, igeog, ishowd, gmask)


c      common variables

      implicit none
#include "drsdef.h"
#include "fcddrs.h"

#include "usersget"
#include "cscalars"
#include "cvartbl"
#include "cdimtbl"
#include "cdomain"

      integer ivar, ishowd
      real gmask(*)

cjfp was      integer getedim, cluvdb, opendrs
      integer getedim, fcw_cluvdb, opendrs
      integer univunit, getfld, xregion
      logical caseindp
      
      integer i, iii, ii, i1, i2, j, jj, idlen(4), idtype, igeog, nvdim,
     +         ia(4), j1, j2, j3, j4, k2, k3, k4, ngdim, 
     +         ja(4), i3, i4, ivlat, ivlon, iglat, iglon, iv, ig, jv, 
     &         jg, ierr, idicun, nlen, lendatal, k, ivp, igp, iim,
     &         jj1, jj2, jj3, jj4

      integer miallocw, mifree
      pointer (pt, tspace)
      real asum, fes, les, tspace(*)

      character dsource*120, dtitle*80, dunits*40, dname(4)*16
            
      xgeog = 0
      
c    how many dimensions?
      nvdim = ivarndim(ivar)
      
      if (ishowd .eq. all) return
         
cjfp was      ierr = cluvdb()
      ierr = fcw_cluvdb()
      if (ierr .ne. 0) then
          xgeog = 1
          call errcheck('DRS function cluvdb', 'function xgeog',1)
          return
      endif
      ierr = opendrs(igeog, 1, idicun)
      if (ierr .ne. 0) then
          xgeog = 1
          call errcheck(' ', 'function xgeog',0)
          return
      endif

         
      ngdim = ivarndim(igeog) 

      do 300 i=1,4
        namedim(i) = 'not defined'
  300 continue
            
      i1 = 1
      ivlat = 0
      ivlon = 0
      iglat = 0
      iglon = 0

      do 400 i=1,ngdim
         ierr = getedim(i,dsource,dname(i),dtitle,dunits,
     +         idtype, idlen(i), fes, les)
         if (ierr .ne. 0) then
            xgeog = 1
            call errcheck('DRS function getedim','function xgeog',1)
            return
         endif
         iii = univunit(dname(i))
         if (iii .eq. intlong) iglon = i
         if (iii .eq. intlat) iglat = i
  400 continue

      i1 = 1
      ivlat = 0
      ivlon = 0
      nlen = 1

      do 600 j=1,4

         jj = iptdimrg(j,ivar)
      
         if (jj .ne. 0) then

          do 500 i=1,ngdim
            if (caseindp(dimnam(jj), dname(i))  .or. 
     &           ((univunit(dname(i)) .gt. 0) .and. 
     &           (univunit(dimnam(jj)) .eq. univunit(dname(i)))))  then
               if (iglon .eq. i) ivlon = j
               if (iglat .eq. i) ivlat = j
               namedim(i1) = dname(i)
               if (idlen(i) .eq. 1) then
                  domainc(i1) = 'assaved'
                  cwt(i1) = 'unit' 
               else
                 namedim(i1) = dname(i)
                 domainc(i1) = 'range'
                 bdry(1,i1) = begbdy(jj)
                 bdry(2,i1) = endbdy(jj)
                 cycle(i1) = dimcycle(jj)
                 cwt(i1) = dimwtflg(jj) 
                 nlen = nlen*idimlen(jj)
                 if (ioriglen(jj) .ne. idlen(i)) then
                   if (errcntl .gt. 0) then
                   print*, 'Warning -- extracted field and geography'
                   print*, '   mask might be stored on different grids.'
                   print*, '   Dimension flagged:  ', dname(i)
                   print*, '   You are extracting ', varnam(ivar)
                   print*, '   found in file: ', filnam(ivar)
                   print*,'   You are also extracting: ', varnam(igeog)
                   print*, '   found in file: ', filnam(igeog)
                   endif
                   call errcheck('function xgeog', ' ', -1)
                 endif 
               endif
               i1 = i1+1
               go to 550
            endif
  500    continue

c        dimension in variable not found in geography

  550    continue
        endif

  600 continue

c   see if any dimension in geography field is not in variable
c   if so, determine which plane to extract

      do 615 i=1,ngdim
        do 605 j=1,4
          jj=iptdimrg(j,ivar)
          if (jj .ne. 0) then
            if (caseindp(dimnam(jj), dname(i)) .or. 
     &           ((univunit(dname(i)) .gt. 0) .and. 
     &           (univunit(dimnam(jj)) .eq. univunit(dname(i))))) 
     &           go to 615
          endif
  605   continue

        if (idlen(i) .eq. 1) then
          namedim(i1) = dname(i)
          domainc(i1) = 'assaved'
          cwt(i1) = 'unit'
          i1 = i1+1
        else
c         find name in the vdimname table
          do 610 k=1,4
            if (caseindp(dname(i), vdimname(k,igeog)) .or. 
     &           ((univunit(dname(i)) .gt. 0) .and. 
     &           (univunit(vdimname(k,igeog)) .eq. univunit(dname(i)))))
     &           then
                  
              if (((vdimdom(k,igeog)(1:7) .eq. 'byindex') .and.
     +           (ivdimpln(1,k,igeog) .eq. ivdimpln(2,k,igeog))) .or.
     +           ((vdimdom(k,igeog)(1:7) .eq. 'nearest') .and.
     +           (vdimbdry(1,k,igeog) .eq. vdimbdry(2,k,igeog)))) then              
                namedim(i1) = vdimname(k,igeog)
                domainc(i1) = vdimdom(k,igeog)
                bdry(1,i1) = vdimbdry(1,k,igeog) 
                bdry(2,i1) = vdimbdry(2,k,igeog) 
                ibdry(1,i1) = ivdimpln(1,k,igeog) 
                ibdry(2,i1) = ivdimpln(2,k,igeog) 
                cycle(i1) = vdimcycl(k,igeog)
                cwt(i1) = vdimwt(k,igeog) 
                 
                i1 = i1+1
                go to 615

              else

                print*, 'Error in obtaining geography grid:'
                print*, 'You should select a single plane for the'
                print*, 'dimension, ', dname(i), ', of variable, '
                print*, varnam(igeog), ', which is stored in file:'
                print*, filnam(igeog)
                print*, 'Do this by calling defdimi or defdim for this'
                print*,'dimension with "byindex" or "nearest" specified'
                print*, 'and the domain boundaries set to the same'
                print*, 'values.'
                xgeog = 1
                call errcheck('function xgeog', ' ',1)
                return
         
              endif

            endif

  610     continue

          print*, 'Error in obtaining geography grid:'
          print*, 'You have failed to specify dimension information for'
          print*, 'dimension, ', dname(i), ', of variable, '
          print*, varnam(igeog), ', which is stored in file:'
          print*, filnam(igeog)
          print*, 'You should select a single plane for this dimension'
          print*, 'by calling defdimi or defdim with "byindex" or' 
          print*, '"nearest" specified and the domain boundaries '   
          print*, 'set to the same values.'
          xgeog = 1
          call errcheck('function xgeog', ' ',1)
          return

        endif

  615 continue

                
      if (iglon .eq. 0) then
      
         print*, 'longitude dimension missing from geography field:'
         print*, 'stored in ', filnam(igeog)
         xgeog = 1
         call errcheck('function xgeog', ' ',1)
         return
            
      endif
      
      if (iglat .eq. 0) then
      
         print*, 'latitude dimension missing from geography field'
         print*, 'stored in ', filnam(igeog)
         xgeog = 1
         call errcheck('function xgeog', ' ',1)
         return
            
      endif
      
      if (ivlon .eq. 0) then
      
         print*, 'longitude dimension missing from variable field'
         print*, 'stored in ', filnam(ivar)
         print*,  'that you are attempting to mask with geography.'
         xgeog = 1
         call errcheck('function xgeog', ' ',1)
         return
            
      endif
      
      if (ivlat .eq. 0) then
      
         print*, 'latitude dimension missing from variable field'
         print*, 'stored in ', filnam(ivar)
         print*,  'that you are attempting to mask with geography.'
         xgeog = 1
         call errcheck('function xgeog', ' ',1)
         return
            
      endif
      
      lendatal = lendata
      lendata = nlen
      ierr = miallocw(nlen, pt)
      if (ierr .ne. 0) then
           print*, 'Error in attempting to allocate memory dynamically.'
           print*, 'You may need more memory.'
           print*, 'Try running on stargate or retrieving smaller '
           print*, 'portions of the array.  Error encountered '
           print*, 'while processing', filnam(igeog), varnam(igeog)
           xgeog = 1
           call errcheck('function miallocw','function xgeog',1)
           return
      endif

      ierr = getfld(igeog, tspace)

      lendata = lendatal

c    check whether geography latitude-longitude grid is identical to
c             grid of variable it will mask

c    check latitude:
      iv = iptdimrg(ivlat,ivar)
      ig = iptdimrg(iglat,igeog)
      jv = idimlen(iv)
      jg = idimlen(ig)

      if (jv .eq. jg) then

        ivp = iptcoord(iv)-1
        igp = iptcoord(ig)-1

        if ((ivp .lt. -1) .or. (igp .lt. -1)) then
        iim = max0(ioriglen(iv), ioriglen(ig))
        print*, 'Dimension vector was not saved but is needed.'
        print*, 
     &    'User should allow longer dimensions to be saved by calling'
        print*, "call defmisc('max dimension kept', 'integer', val)"
        print*, "where val .gt. ", iim
        xgeog = 1
        call errcheck('function xgeog', ' ',1)
        return
        endif

        asum = 0.

        do 625 j=1,jv
           asum = asum + abs(coordvls(ivp+j) - coordvls(igp+j))
  625   continue
        asum = asum/jv

        if (asum .gt. 0.1) then
           print*, 'latitude coordinates for extracted field and'
           print*, 'geography mask disagree on average by more than' 
           print*, '0.1 degree latitude'
           xgeog = 1
           call errcheck('function xgeog', ' ',1)
           return
        elseif (asum .gt. 0.0004) then
           if (errcntl .gt. 0) then
           print*, 'Warning --'
           print*, 'latitude coordinates for extracted field and'
           print*, 'geography mask disagree on average by more than'
           print*, '0.0004 degree latitude'
           endif
           call errcheck('function xgeog', ' ',-1)
        endif

      else
        print*, 'extracted latitude dimension for extracted field'
        print*, 'and geography mask disagree'
        print*, 'geography dimension = ', jg
        print*, 'extracted field dimension = ', jv
        xgeog = 1
        call errcheck('function xgeog', ' ',1)
        return
      endif
          
c    check longitude:
      iv = iptdimrg(ivlon,ivar)
      ig = iptdimrg(iglon,igeog)
      jv = idimlen(iv)
      jg = idimlen(ig)

      if (jv .eq. jg) then

         ivp = iptcoord(iv)
         igp = iptcoord(ig)

        if ((ivp .lt. -1) .or. (igp .lt. -1)) then
        iim = max0(ioriglen(iv), ioriglen(ig))
        print*, 'Dimension vector was not saved but is needed.'
        print*, 
     &    'User should allow longer dimensions to be saved by calling'
        print*, "call defmisc('max dimension kept', 'integer', val)"
        print*, "where val .gt. ", iim
        xgeog = 1
        call errcheck('function xgeog', ' ',1)
        return
        endif

         asum = 0.

         do 675 j=1,jv
           asum = asum + abs(coordvls(ivp-1+j) - coordvls(igp-1+j))
  675    continue
         asum = asum/jv

         if (asum .gt. 0.1) then
           print*, 'longitude coordinates for extracted field and'
           print*, 'geography mask disagree on average by more than' 
           print*, '0.1 degree longitude'
           xgeog = 1
           call errcheck('function xgeog', ' ',1)
           return
         elseif (asum .gt. 0.0004) then
           if (errcntl .gt. 0) then
           print*, 'Warning --'
           print*, 'longitude coordinates for extracted field and'
           print*, 'geography mask disagree on average by more than'
           print*, '0.0004 degree longitude'
           endif
           call errcheck('function xgeog', ' ',-1)
         endif
      else
         print*, 'extracted longitude dimension for extracted field'
         print*, 'and geography mask disagree'
         print*, 'geography dimension = ', jg
         print*, 'extracted field dimension = ', jv
         xgeog = 1
         call errcheck('function xgeog', ' ',1)
         return
      endif
         
      if (ierr .ne. 0) then
         xgeog = 1
         call errcheck(' ', 'function xgeog',0)
         return
      endif
         
      if (ishowd .le. 0) then

         print*, 'unrecognized value of ishow: ',ishowd
         print*, 'must specify correct type of geography masking.'
         xgeog = 1
         call errcheck('function xgeog', ' ',1)
         return

      endif

      ierr = xregion(ishowd, vartype(igeog), varnam(igeog),
     +       ivarlen(igeog), tspace)
         
      ii = 1
      do 1000 j=1,4
          i = iptdimrg(j,ivar)
          if (i .eq. 0) then
              ia(j) = 1
              ja(j) = 1
          else
              ja(j) = idimlen(i)
              iii = iptdim(ii,igeog)
              ia(j) = 1
              if (iii .gt. 0) then
                if (caseindp(dimnam(iii), dimnam(i)) .or.
     &            (univunit(dimnam(iii)) .gt. 0) .and.
     &            (univunit(dimnam(iii)) .eq. univunit(dimnam(i)))) then
                   if (idimlen(iii) .gt. 1) ia(j) = ja(j)
                   ii = ii + 1
                endif
              endif
          endif
1000  continue


      jj1 = ja(1)
      jj2 = ja(2)
      jj3 = ja(3)
      jj4 = ja(4)


c    copy tspace to other planes of field

      do 1400 j4 = 1,jj4
         k4 = ja(3)*(j4-1)
         i4 = ia(3)*(min0(j4,ia(4))-1)
         do 1300 j3 = 1,jj3
            k3 = ja(2)*(j3-1+k4)
            i3 = ia(2)*(min0(j3,ia(3))-1+i4)
            do 1200 j2 = 1,jj2
               k2 = ja(1)*(j2-1+k3)
               i2 = ia(1)*(min0(j2,ia(2))-1+i3)
               do 1100 j1 = 1,jj1
                  i1 = i2 + min0(j1,ia(1))
                  gmask(k2+j1) = tspace(i1)*gmask(k2+j1)
 1100          continue
 1200       continue
 1300    continue
 1400 continue

      ierr = mifree(pt)

      return
      end





c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      function xregion
c     -------
c
c     Description:
c     -----------
c
c      This function creates a geography mask (geogmask)
c      filling it with 0. and 1. depending on whether the
c      point should be included or not.
c
c     Usage:
c     ------
c
c      xregion (ishowd, vartyped, varnaml, jj, geogmask)
c
c
c     Exit Status:
c     -----------
c     0 if success
c     1 if an error occurred
c
c     ------
c     Date: 1/4/95
c     ----
c
c
c
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



      integer function xregion (ishowd, vartyped, varnaml, jj, geogmask)


c      common variables

      implicit none
#include "cscalars"

      integer ishowd, jj
      real geogmask(*)
      character*(*) vartyped, varnaml

      integer ka(50), ma(mxsfctyp), imsk, iii, j, ia(4), kk, k, ii,
     +    jmsk
      real rmsk 

      equivalence (imsk, rmsk)

      xregion = 0

      ii = 2**30         

      if ((masktype .eq. 1) .or. ((masktype .eq. 0) .and. 
     &    ((varnaml .eq. 'sftl') .or. (varnaml .eq. 'sftland')))) then

c      the geography is expressed as the percent land in each grid cell.

        if (ishowd .lt. seaice) then
         print*, ' '
         print*, 'Accessed geography file contains land fraction '
         print*, '   (expressed as a percent).  You therefore must not'
         print*, '   ask for individual geographical regions.'
         print*, 'your only options are: '
         print*, '        "land" and/or "ocean,  sea ice".'
         print*, ' '
         xregion = 1
         call errcheck('function xregion', ' ',1)
         return
        endif

        iii = ishowd
        ia(1) = 0
        ia(2) = 0
        ia(3) = 0
        if (iii .ge. ocean) then
          ia(1) = 1
          iii = iii - ocean
        endif
        if (iii .ge. land) then
          ia(2) = 1
          iii = iii - land
        endif
        if (iii .ge. seaice) then 
          ia(3) = 1
          iii = iii - seaice
        endif

        if (iii .ne. 0) then
         print*, ' '
         print*, 'Accessed geography file contains land fraction '
         print*, '   (expressed as a percent).  You therefore must not'
         print*, '   ask for individual geographical regions.'
         print*, 'your only options are: '
         print*, '        "land" and/or "ocean, sea ice".'
         print*, 'Program will proceed and ignore individual '
         print*, '       regions specified.'
         print*, ' '
        endif

        if ((ia(1) + ia(3)) .eq. 1) then
         print*, ' '
         print*, 'Accessed geography file contains land fraction '
         print*, '   (expressed as a percent).  You therefore must not'
         print*, '   select "ocean" without also selecting "sea ice".'
         print*, 'Program will proceed and select both "ocean" and'
         print*, '      "sea ice".' 
         print*, ' '
         ia(1) = 1
         ia(3) = 1
        endif

        if (ishowd .ge. ocean+seaice+land) then

          do 100 j=1,jj
            geogmask(j) = 1.0
  100     continue
          
        elseif (ia(1) .eq. 1) then
       
          if ((vartyped .eq. 'R*4').or.(vartyped .eq. 'R*8')) then

            do 200 j=1,jj
              if (geogmask(j) .gt. 99.999) then
                geogmask(j) = 0.0
              else
                geogmask(j) = 1.0 - geogmask(j)/100.
              endif
  200       continue

          else

            do 250 j=1,jj
              rmsk = geogmask(j)
              if (imsk .lt. 100) then
                geogmask(j) = 1.0 - imsk/100.
              else
                geogmask(j) = 0.0
              endif
250         continue
   
          endif

        else

          if ((vartyped .eq. 'R*4') .or. (vartyped .eq. 'R*8')) then

            do 300 j=1,jj
              geogmask(j) = geogmask(j)/100.
  300       continue

          else

            do 350 j=1,jj
              rmsk = geogmask(j)
              geogmask(j) = imsk/100.
  350       continue
   
          endif

        endif

        return
      endif

      if ((masktype .eq. 2) .or. ((masktype .eq. 0) .and. 
     &    ((varnaml .eq. 'sic') .or. (varnaml .eq. 'seaice')))) then

c    the geography is expressed as the percent seaice in each grid cell.

        if (ishowd .lt. seaice) then
         print*, ' '
         print*, 'Accessed geography file contains sea ice fraction '
         print*, '   (expressed as a percent).  You therefore must not'
         print*, '   ask for individual geographical regions.'
         print*, 'your only options are: '
         print*, '        "sea ice" and/or "land, ocean" '
         print*, ' '
         xregion = 1
         call errcheck('function xregion', ' ',1)
         return
        endif

        iii = ishowd
        ia(1) = 0
        ia(2) = 0
        ia(3) = 0
        if (iii .ge. ocean) then
          ia(1) = 1
          iii = iii - ocean
        endif
        if (iii .ge. land) then
          ia(2) = 1
          iii = iii - land
        endif
        if (iii .ge. seaice) then
          ia(3) = 1
          iii = iii - seaice
        endif

        if (iii .ne. 0) then
         print*, ' '
         print*, 'Accessed geography file contains sea ice fraction '
         print*, '   (expressed as a percent).  You therefore must not'
         print*, '   ask for individual geographical regions.'
         print*, 'your only options are: '
         print*, '        "sea ice" and/or "ocean, land".'
         print*, 'Program will proceed and ignore individual '
         print*, '       regions specified.'
         print*, ' '
        endif

        if ((ia(1) + ia(2)) .eq. 1) then
         print*, ' '
         print*, 'Accessed geography file contains sea ice fraction '
         print*, '   (expressed as a percent).  You therefore must not'
         print*, '   select "ocean" without also selecting "land".'
         print*, 'Program will proceed and select both "ocean" and'
         print*, '      "land".' 
         print*, ' '
         ia(1) = 1
         ia(3) = 1
        endif

        if (ishowd .ge. ocean+seaice+land) then

          do 400 j=1,jj
            geogmask(j) = 1.0
  400     continue
          
        elseif (ia(3) .eq. 1) then
       
          if ((vartyped .eq. 'R*4').or.(vartyped .eq. 'R*8')) then

            do 500 j=1,jj
              geogmask(j) = geogmask(j)/100.
  500       continue

          else

            do 550 j=1,jj
              rmsk = geogmask(j)
              geogmask(j) = imsk/100.
550       continue
   
          endif

        else

          if ((vartyped .eq. 'R*4').or.(vartyped .eq. 'R*8')) then

            do 600 j=1,jj
              if (geogmask(j) .gt. 99.999) then
                geogmask(j) = 0.0
              else
                geogmask(j) = 1.0 - geogmask(j)/100.
              endif
  600       continue

          else

            do 650 j=1,jj
              rmsk = geogmask(j)
              if (imsk .lt. 100) then
                geogmask(j) = 1.0 - imsk/100.
              else
                geogmask(j) = 0.0
              endif
  650       continue
   
          endif

        endif

        return
      endif


      if ((vartyped .eq. 'R*4') .or. 
     +            (vartyped .eq. 'R*8')) then 
           imsk = geogmask(1) + 0.01
      else
           rmsk = geogmask(1)
      endif

      if ((masktype .eq. 3) .or. ((masktype .eq. 0) .and.
     &              ((imsk .lt. 3) .and. (imsk .ge. 0)))) then
         
c          assume we have read a geography map with 3 surface types.
         
c          determine which surface types were requested:
c          assume imsk contains 0=ocean, 1=land, 2=sea-ice
                   
         iii = ishowd
         do 800 j = 1,3
            ia(j) = -1
           
            if (iii .ge. ii) then
               ia(j) = j - 1
               iii = iii - ii
            endif
      
            ii = ii/2
  800    continue

         if ((ishowd .lt. 2**28) .or. (iii .ne. 0)) then 

           print*, ' '
           print*, "Accessed geography file contains 0's, 1's, and"
           print*, "possibly 2's (or 0.0's, 1.0's, and 2.0's), "
           print*, 'indicating presence of ocean, land and sea ice'
           print*, '  You therefore must not ask for individual'
           print*, '  geographical regions. Your only options are: '
           print*, '        "ocean", and/or "land" and/or "sea ice".'
           print*, ' '
           xregion = 1
           call errcheck('function xregion', ' ',1)
           return
        
         endif

              
         if ((vartyped .eq. 'R*4') .or. (vartyped .eq. 'R*8')) then 
           do 805 j=1,jj
               imsk = geogmask(j) + 0.01
               if  ((imsk .eq. ia(1)) .or. (imsk .eq. ia(2))
     +             .or. (imsk .eq. ia(3))) then
                  geogmask(j)=1.
               else
                  geogmask(j)=0.
               endif
  805      continue
         else
           do 810 j=1,jj
               rmsk = geogmask(j)
               if  ((imsk .eq. ia(1)) .or. (imsk .eq. ia(2))
     +             .or. (imsk .eq. ia(3))) then
                  geogmask(j)=1.
               else
                  geogmask(j)=0.
               endif
  810      continue
         endif

         return

      endif
            
      if ((masktype .eq. 4) .or. (masktype .eq. 0)) then   

c          assume we have read a geography map with regions designated
c          by numbers > 100.
                
c          ocean and seaice:
          ma(1) = 202
          ma(2) = 203
          ma(3) = 204
      
          ma(4) = 213
          ma(5) = 214
      
          ma(6) = 201
          ma(7) = 209
      
          ma(8) = 205
          ma(9) = 206
      
          ma(10) = 210
          ma(11) = 207
          ma(12) = 212
      
c          land:
      
          ma(13) = 219
          ma(14) = 217
          ma(15) = 216
          ma(16) = 218
          ma(17) = 221
          ma(18) = 220
          ma(19) = 215
      
c          islands:
          ma(20) = 222
     
c          lakes:
          ma(21) = 208
          ma(22) = 211
      
c      determine which surface types were requested:
      
          kk = 0
          iii = ishowd
          j = 1
         
          if (((iii .ge. ocean) .and. (iii-ocean .ge. land)) .or.
     +       ((iii .lt. ocean) .and. (iii .ge. land))) then
c     show all land points ( and islands but not inland seas and lakes)

             do 820 k=1,8
                ka(kk+k) = ma(12+k)
  820        continue
  
          kk=kk+8
          iii = iii - land
         endif
                  
         if (iii .ge. seaice) then
            if (iii .lt. seaice+ocean) then
         
              print*, ' '
              print*, 'Accessed geography file contains detailed '
              print*, '  designation of geographical regions, but no'
              print*, '  indication of sea ice.  You'
              print*, '  therefore must not select "ocean" without also' 
              print*,'  selecting "sea ice".  Program will proceed and ' 
              print*, '  select both "ocean" and "sea ice".' 
              print*, ' '

               if (iii .ge. ocean)  then
                   iii = iii + seaice
               else
                   iii = iii + ocean
               endif
            endif
            
c            show all ocean and seaice points (and seas and lakes)

            do 840 k=1,12
               ka(kk+k) = ma(k)
  840       continue
  
            ka(kk+13) = ma(21)
            ka(kk+14) = ma(22)
           
            kk=kk+14
            iii = iii - ocean - seaice
         endif


         ii = 2**(30-4)
         do 850 k=1, mxsfctyp
            
            if (iii .ge. ii) then
               kk=kk+1
               ka(kk) = ma(k)
               iii = iii-ii
            endif
            if (iii .eq. 0) goto 875
            ii = ii/2
  850    continue
            
  875    continue
            
         jmsk  = -1    
         ii = 1
         if ((vartyped .eq. 'R*4') .or. (vartyped .eq. 'R*8')) then 

            do 975 j=1,jj
               imsk = geogmask(j) + 0.01
               if (imsk .eq. ka(ii)) then
                  geogmask(j) = 1.
                  go to 975
               elseif (imsk .eq. jmsk) then
                  geogmask(j) = 0.
                  go to 975
               endif
               do 950 k=1,kk
                  if (imsk .eq. ka(k)) then
                     geogmask(j) = 1.
                     ii=k
                     go to 975
                  endif
  950          continue
  
               jmsk = imsk
               geogmask(j) = 0.
  975       continue
  
         else

            do 985 j=1,jj
               rmsk = geogmask(j)
               if (imsk .eq. ka(ii)) then
                  geogmask(j) = 1.
                  go to 985
               elseif (imsk .eq. jmsk) then
                  geogmask(j) = 0.
                  go to 985
               endif
               do 980 k=1,kk
                  if (imsk .eq. ka(k)) then
                     geogmask(j) = 1.
                     ii=k
                     go to 985
                  endif
  980          continue
  
               jmsk = imsk
               geogmask(j) = 0.
  985       continue
  
         endif

         return

      endif

      print*, ' '
      print*, 'EzGet was unable to mask geographical regions'
      print*, 'because the file containing land fraction or'
      print*, 'or the sea ice fraction or other geography data'
      print*, 'could not be interpreted.  Please call subroutine'
      print*, 'defmisc and specify "mask type", or "mask type in"'
      print*, 'or "mask type out", as appropriate.'
      print*, ' '
      xregion = 1
      call errcheck('function xregion', ' ',1)

      return
      end




      SUBROUTINE gauaw(PA,PW,K)
C
C*** *GAUAW* - COMPUTE ABSCISSAS AND WEIGHTS FOR *GAUSSIAN INTEGRATION.
C
C     PURPOSE.
C     --------
C
C          *GAUAW* IS CALLED TO COMPUTE THE ABSCISSAS AND WEIGHTS REQUIR
C     TO PERFORM *GAUSSIAN INTEGRATION.
C
C**   INTERFACE.
C     ----------
C
C          *CALL* *GAUAW(PA,PW,K)*
C
C               *PA*     - ARRAY, LENGTH AT LEAST *K,* TO RECEIVE ABSCIS
C                          ABSCISSAS.
C               *PW*     - ARRAY, LENGTH AT LEAST *K,* TO RECEIVE
C                          WEIGHTS.
C
C     METHOD.
C     -------
C
C          THE ZEROS OF THE *BESSEL FUNCTIONS ARE USED AS STARTING
C     APPROXIMATIONS FOR THE ABSCISSAS. NEWTON ITERATION IS USED TO
C     IMPROVE THE VALUES TO WITHIN A TOLLERENCE OF *EPS.*
C
C     EXTERNAL.
C     ---------
C
C          *BSSLZR* - ROUTINE TO OBTAIN ZEROS OF *BESSEL FUNCTIONS.
C
C     REFERENCE.
C     ----------
C
c      implicit double precision (a-h, o-z)
      implicit none
      integer K, IS, ITER, N, L, KK
      double precision EPS, PI, C, FK, XZ, PKM2, PKM1, FN, PK, PKMRK,
     &    SP, AVSP, PA(K), PW(K)
      DATA EPS/1.d-13/
C
C     ------------------------------------------------------------------
C
C*         1.     SET CONSTANTS AND FIND ZEROS OF BESSEL FUNCTION.
C                 --- --------- --- ---- ----- -- ------ ---------
C
      PI = 4.*datan(1.d0)
  100 CONTINUE
      C=(1.-(2./PI)**2)*0.25
      FK=K
      KK=K/2
      CALL BSSLZR(PA,KK)
C
      DO 290 IS=1,KK
      XZ=dcos(PA(IS)/dsqrt((FK+0.5)**2+C))
C*                 GIVING THE FIRST APPROXIMATION FOR *XZ.*
      ITER=0
C
C     ------------------------------------------------------------------
C
C*         2.     COMPUTE ABSCISSAS AND WEIGHTS.
C                 ------- --------- --- -------
C
  200 CONTINUE
C
C*         2.1     SET VALUES FOR NEXT ITERATION.
  210 CONTINUE
      PKM2=1.
      PKM1=XZ
      ITER=ITER+1
      IF(ITER.GT.10) GO TO 300
C
C*         2.2     COMPUTATION OF THE *LEGENDRE POLYNOMIAL.
  220 CONTINUE
C
      DO 222 N=2,K
      FN=N
      PK=((2.*FN-1.)*XZ*PKM1-(FN-1.)*PKM2)/FN
      PKM2=PKM1
      PKM1=PK
  222 CONTINUE
C
      PKM1=PKM2
      PKMRK=(FK*(PKM1-XZ*PK))/(1.-XZ**2)
      SP=PK/PKMRK
      XZ=XZ-SP
      AVSP=dabs(SP)
      IF(AVSP.GT.EPS) GO TO 210
C
C*         2.3     ABSCISSAS AND WEIGHTS.
  230 CONTINUE
      PA(IS)=XZ
      PW(IS)=(2.*(1.-XZ**2))/(FK*PKM1)**2
C
C*         2.4     ODD *K* COMPUTATION OF WEIGHT AT THE EQUATOR.
  240 CONTINUE
      IF (K.NE.KK*2) THEN
         PA(KK+1)=0.
         PK=2./FK**2
C
         DO 242 N=2,K,2
         FN=N
         PK=PK*FN**2/(FN-1.)**2
  242    CONTINUE
C
         PW(KK+1)=PK
      ELSE
C
C*         2.5     USE SYMMETRY TO OBTAIN REMAINING VALUES.
C
  250    CONTINUE
C
         DO 252 N=1,KK
         L=K+1-N
         PA(L)=-PA(N)
         PW(L)=PW(N)
  252    CONTINUE
C
      ENDIF
  290 CONTINUE
C      write(6,*) iter
C
      RETURN
C
C     ------------------------------------------------------------------
C
C*         3.     ERROR PROCESSING.
C                 ----- -----------
C
  300 CONTINUE
      print*, '  GAUAW FAILED TO CONVERGE AFTER 10 ITERATIONS.'
c      WRITE (NOUT,9901)
c 9901 FORMAT(//,'  GAUAW FAILED TO CONVERGE AFTER 10 ITERATIONS.')
      STOP
C
C     ------------------------------------------------------------------
C
      END




      SUBROUTINE BSSLZR(PBES,KNUM)
C
C**** *BSSLZR* - ROUTINE TO RETURN ZEROS OF THE J0 *BESSEL FUNCTION.
C
C     PURPOSE.
C     --------
C
C          *BSSLZR* RETURNS *KNUM* ZEROS, OR IF *KNUM>50,* *KNUM*
C     APPROXIMATE ZEROS OF THE *BESSEL FUNCTION J0.
C
C**   INTERFACE.
C     ----------
C
C          *CALL* *NSSLZR(PBES,KNUM)*
C
C               *PBES*   - ARRAY, DIMENSIONED *KNUM,* TO RECEIVE THE
C                          VALUES.
C               *KNUM*   - NUMBER OF ZEROS REQUESTED.
C
C     METHOD.
C     -------
C
C          THE FIRST 50 VALUES ARE OBTAINED FROM A LOOK UP TABLE. ANY
C     ADDITIONAL VALUES REQUESTED ARE INTERPOLATED.
C
C     EXTERNALS.
C     ----------
C
C          NONE.
C
C     REFERENCE.
C     ----------
C
      implicit none
      integer KNUM, INUM, J 
      double precision PBES(KNUM), ZBES(50), API
      DATA ZBES        / 2.4048255577,   5.5200781103,
     X    8.6537279129,  11.7915344391,  14.9309177086,  18.0710639679,
     X   21.2116366299,  24.3524715308,  27.4934791320,  30.6346064684,
     X   33.7758202136,  36.9170983537,  40.0584257646,  43.1997917132,
     X   46.3411883717,  49.4826098974,  52.6240518411,  55.7655107550,
     X   58.9069839261,  62.0484691902,  65.1899648002,  68.3314693299,
     X   71.4729816036,  74.6145006437,  77.7560256304,  80.8975558711,
     X   84.0390907769,  87.1806298436,  90.3221726372,  93.4637187819,
     X   96.6052679510,  99.7468198587, 102.8883742542, 106.0299309165,
     X  109.1714896498, 112.3130502805, 115.4546126537, 118.5961766309,
     X  121.7377420880, 124.8793089132, 128.0208770059, 131.1624462752,
     X  134.3040166383, 137.4455880203, 140.5871603528, 143.7287335737,
     X  146.8703076258, 150.0118824570, 153.1534580192, 156.2950342685/
C
C     ------------------------------------------------------------------
C
C*         1.     EXTRACT VALUES FROM LOOK UP TABLE.
C                 ------- ------ ---- ---- -- ------
C
C                 SET API
C
      API=4.*datan(1.d0)
  100 CONTINUE
      INUM=MIN0(KNUM,50)
C
      DO 110 J=1,INUM
      PBES(J)=ZBES(J)
  110 CONTINUE
C
C     ------------------------------------------------------------------
C
C*         2.     INTERPOLATE REMAINING VALUES.
C                 ----------- --------- -------
C
  200 CONTINUE
C
      DO 210 J=51,KNUM
      PBES(J)=PBES(J-1)+API
  210 CONTINUE
C
C     ------------------------------------------------------------------
C
      RETURN
      end
