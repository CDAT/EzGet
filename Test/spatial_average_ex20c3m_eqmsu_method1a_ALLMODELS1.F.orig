c-------------------------------------------------------------------------------
c-------------------------------------------------------------------------------
c     SPATIAL_AVERAGE_EX20C3M_EQMSU_METHOD1a_ALLMODELS1.F       
c
c     Date:    07-JAN-2011.
c     Author:  Ben Santer, LLNL, Livermore.
c     Version: # 1.0 / Single precision.
c
c
c
c     Input data:    Monthly-mean tam4/tam2/tam6 data from "extended" 20c3m
c                    runs. The extension involves appending data from the
c                    relevant ex20c3m run to the end of the 20c3m run.
c
c     Purpose:       Calculates spatial averages of upper-air temperature data
c                    over user-defined regions.
c
c                    Upper-air data are either tam4, tam2, or tam6, where:
c
c                      tam4 = MSU channel 4
c                      tam2 = MSU channel 2
c                      tam6 = MSU 2LT retrieval (LT = "Lower Troposphere")
c
c                    After calculation of spatial averages, temperature data 
c                    are expressed as anomalies with respect to user-selected 
c                    reference period. Processes model data only. Processes
c                    only one variable and one region at a time.
c
c     HISTORY:
c     04-MAY-2010    Initial version.         
c     04-JUN-2010    Added new regions (NHEX and SHEX).
c     14-JUN-2010    Corrected error affecting masking of land grid-points.
c     24-JUN-2010    Added calculation and output of multi-model average time
c                     series (average of ensemble means).
c     25-JUN-2010    Added new regions (NHM1, NHM2, NHM3, NHM4, NHM5, NHM6,  
c                     SHM1, SHM2, SHM3, SHM4, SHM5, and SHM6).
c     22-JUL-2010    Added new region (SHM7).
c     05-AUG-2010    Added new regions (GLB2 and GLB3).
c     06-JAN-2011    Added new model (csiro_mk3_5).  
c-------------------------------------------------------------------------------
c-------------------------------------------------------------------------------
c
c
c
      program amplify1
      implicit real(o-z, a-h)
c
c
c
c-------------------------------------------------------------------------------
c     ** parameters **
c
      parameter(nmod = 20)
      parameter(nr1   = 1)     ! bccr_bcm_2_0        no. of realizations
      parameter(nr2   = 5)     ! cccma_cgcm3_1       no. of realizations
      parameter(nr3   = 1)     ! cccma_cgcm3_1_t63   no. of realizations
      parameter(nr4   = 1)     ! csiro_mk3_5         no. of realizations
      parameter(nr5   = 1)     ! gfdl_cm2_0          no. of realizations
      parameter(nr6   = 1)     ! gfdl_cm2_1          no. of realizations
      parameter(nr7   = 2)     ! giss_aom            no. of realizations
      parameter(nr8   = 3)     ! giss_model_e_h      no. of realizations
      parameter(nr9   = 5)     ! giss_model_e_r      no. of realizations
      parameter(nr10  = 3)     ! iap_fgoals_1_0_g    no. of realizations
      parameter(nr11  = 1)     ! ingv_echam4         no. of realizations
      parameter(nr12  = 1)     ! inmcm3_0            no. of realizations
      parameter(nr13  = 1)     ! miroc3_2_hires      no. of realizations
      parameter(nr14  = 3)     ! miroc3_2_medres     no. of realizations
      parameter(nr15  = 4)     ! mpi_echam5          no. of realizations
      parameter(nr16  = 5)     ! mri_cgcm2_3_2a      no. of realizations
      parameter(nr17  = 7)     ! ncar_ccsm3_0        no. of realizations
      parameter(nr18  = 4)     ! ncar_pcm1           no. of realizations
      parameter(nr19  = 1)     ! ukmo_hadcm3         no. of realizations
      parameter(nr20  = 1)     ! ukmo_hadgem1        no. of realizations
c
      parameter(nrmax = 7)
      parameter(ntot  = nrmax * nmod)
c
      parameter(nreg =   33)
      parameter(nx   =  361)
      parameter(ny   =  160)
      parameter(nm   =   12)
      parameter(nt   = 1001)
      parameter(nmt  =  nm * nt)
      parameter(lag  =   10)
c
c     nmod      Maximum number of models processed
c     nr1-nr20  Maximum number of model realizations processed
c     nrmax     Maximum number of realizations per model      
c     ntot      Maximum number of input files processed       
c     nreg      Maximum number of study areas                          
c     nx        Maximum input spatial dimension (longitude; MIROC T106)
c     ny        Maximum input spatial dimension (latitude; MIROC T106)
c     nt        Maximum number of years processed
c     nmt       Maximum compressed time dimension (months)
c     lag       Maximum number of lags computed            
c
      include '/usr/local/include/drsdef.h'
c
c
c
c-------------------------------------------------------------------------------
c     ** User options **
c
c     character*16 pv1          1. Name of input variable
c
c     character*16 creg1        2. Selected region
c                                   TO15 = 15N-15S,    180W-180E (Ocean only)
c                                   TO20 = 20N-20S,    180W-180E (Ocean only)
c                                   TO30 = 30N-30S,    180W-180E (Ocean only)
c                                   TL15 = 15N-15S,    180W-180E (Land only)
c                                   TL20 = 20N-20S,    180W-180E (Land only)
c                                   TL30 = 30N-30S,    180W-180E (Land only)
c                                   TLO1 = 15N-15S,    180W-180E (Land and ocean)
c                                   TLO2 = 20N-20S,    180W-180E (Land and ocean)
c                                   TLO3 = 30N-30S,    180W-180E (Land and ocean)
c                                   NI34 =  5N-5S,     170W-120W (Nino 3.4 region)
c                                   OC50 = 50N-50S,    180W-180E (Ocean only)
c                                   SHEM = 90S- 0 ,    180W-180E (Land and ocean)
c                                   NHEM =  0 -90N,    180W-180E (Land and ocean)
c                                   GLBL = 90S-90N,    180W-180E (Land and ocean)
c                                   GLAN = 90S-90N,    180W-180E (Land only)
c                                   GOCN = 90S-90N,    180W-180E (Ocean only)
c                                   SHEX = 90S-20S,    180W-180E (Land and ocean)
c                                   NHEX = 20N-90N,    180W-180E (Land and ocean)
c                                   NHM1 =  0 -10N,    180W-180E (Land and ocean)
c                                   NHM2 = 10N-30N,    180W-180E (Land and ocean)
c                                   NHM3 = 30N-50N,    180W-180E (Land and ocean)
c                                   NHM4 = 50N-70N,    180W-180E (Land and ocean)
c                                   NHM5 = 70N-82.5N,  180W-180E (Land and ocean)
c                                   NHM6 =  0 -82.5N,  180W-180E (Land and ocean)
c                                   SHM1 =  0 -10S,    180W-180E (Land and ocean)
c                                   SHM2 = 10S-30S,    180W-180E (Land and ocean)
c                                   SHM3 = 30S-50S,    180W-180E (Land and ocean)
c                                   SHM4 = 50S-70S,    180W-180E (Land and ocean)
c                                   SHM5 = 70S-82.5S,  180W-180E (Land and ocean)
c                                   SHM6 =  0 -70S,    180W-180E (Land and ocean)
c                                   SHM7 =  0 -82.5S   180W-180E (Land and ocean)
c                                   GLB2 = 82.5S-82.5N 180W-180E (Land and ocean)
c                                   GLB3 = 70S-82.5N   180W-180E (Land and ocean)
c
c     character*16 ctype        3. Determines type of output data
c                                   AVE = Spatial average of data 
c
c     integer iyr1, iyr2        4. First, last year for time-mean computation
c     integer iyr3, iyr4        5. First, last year for time series statistics           
c                                   If either iyr1 or iyr2 = 0, no time-mean
c                                   is computed or subtracted
c
c     integer itail             6. 1=one-tailed test / 2=two-tailed test
c
c     integer ires              7. 0=AR1 from raw data / 1=AR1 from residuals
c
c     real savcrit              8. Threshold for computation of spatial mean
c
      character pv1*16, pv2*16, creg1*16, creg2*16, ctype*16
      integer iyr1, iyr2, iyr3, iyr4, itail, ires
      real savcrit
c
c-------------------------------------------------------------------------------
c     ** Data and mask variables **
c
c     integer iyear1(nrmax, nmod)     First complete year of ex20c3m run
c     integer iyear2(nrmax, nmod)     Last complete year of ex20c3m run
c     dimension c0(nx, ny)            Current input data
c     dimension c2(nmt, nrmax, nmod)  Stored spatial averages
c     dimension fmask(nx, ny)         Fractional land/sea mask
c     dimension bmask1(nx, ny)        Input data mask
c     dimension sm1(nm, nrmax, nmod)  Climatological monthly means
c
      integer iyear1(nrmax, nmod), iyear2(nrmax, nmod)
      dimension c0(nx, ny)
      dimension c2(nmt, nrmax, nmod)
      dimension fmask(nx, ny), bmask1(nx, ny)
      dimension sm1(nm, nrmax, nmod)
c
c-------------------------------------------------------------------------------
c     ** Other data variables **
c
c     integer ileng1(nrmax, nmod)     Total number of months of data
c     real x1(nreg), x2(nreg)         Longitude boundaries of study area
c     real y1(nreg), y2(nreg)         Latitude boundaries of study area
c     dimension ts1(nmt, nrmax, nmod) Spatially-averaged temperature (all times)
c     dimension ts2(nmt, nrmax, nmod) Spatially-averaged temperature (iyr3-iyr4)
c     dimension ts3(nmt)              Multi-model average time series (ts2)
c     dimension ts4(nmt)              Inter-model variability (ts2)   
c     dimension ts5(nmt)              No. of valid model results (ts2)   
c     dimension em1(nmt, nmod)        Ensemble mean time series
c     real xm                         Climatological mean
c
      integer ileng1(nrmax, nmod)
      real x1(nreg), x2(nreg), y1(nreg), y2(nreg)
      dimension ts1(nmt, nrmax, nmod), ts2(nmt, nrmax, nmod)
      dimension ts3(nmt), ts4(nmt), ts5(nmt)
      dimension em1(nmt, nmod)
      real xm
c
c-------------------------------------------------------------------------------
c     ** Workspace **
c
c     dimension qt(nmod)              Workspace
c     dimension xt(nt)                Workspace
c     dimension yt(nmt)               Workspace
c     dimension zt(nmt)               Workspace
c     dimension work(nx, ny)          Workspace
c     dimension yy(16384)             Storage
c     dimension zz(16384)             Storage
c     dimension xf(16384)             Storage
c     dimension wk3(4, nmod)          Workspace
c     dimension wk4(nmod)             Workspace
c     dimension wk5in(nmod)           Workspace
c     dimension wk5out(nmod)          Workspace
c     dimension wk6(5, 4)             Workspace
c
      dimension qt(nmod), xt(nt), yt(nmt), zt(nmt), work(nx, ny)
      dimension yy(16384), zz(16384), xf(16384)
      dimension wk3(4, nmod), wk4(nmod)
      dimension wk5in(nmod), wk5out(nmod), wk6(5, 4)
c
c-------------------------------------------------------------------------------
c     ** Time series statistics **
c
c     integer ngptmax(nmt, nrmax, nmod)   Number of grid-points in domain
c     integer ngptact(nmt, nrmax, nmod)   Number of grid-points with valid data
c     integer miss1(nmt, nrmax, nmod)     Number of months with missing data   
c     integer iqcrit(nrmax, nmod)         Number of unusual time series values
c     dimension pqcrit(nrmax, nmod)       Percent unusual time series values
c     dimension fcover(nmt, nrmax, nmod)  Fractional coverage            
c     dimension qnorm1(nmt, nrmax, nmod)  Normalized first differences
c     dimension qnorm2(nrmax, nmod)       Ave. normalized first differences
c     dimension cxx(0:lag)                Autocovariance
c     dimension rxx(0:lag)                Autocorrelation
c     dimension result(31)                Output from REG_AR0
c     dimension res(nmt)                  Regression residuals 
c     dimension ybar(nrmax, nmod)         Mean
c     dimension sdy(nrmax, nmod)          Standard deviation
c     dimension r1(nrmax, nmod)           Lag-1 autocorrelation
c     dimension trend(nrmax, nmod)        Least-squares linear trend
c     dimension tlc(nrmax, nmod)          Total linear change
c     dimension serr1(nrmax, nmod)        Standard error (unadjusted)
c     dimension serr2(nrmax, nmod)        Standard error (adjusted)
c     dimension ptt(nrmax, nmod)          p-value for t-value (AdjSE+AdjDF)
c     dimension rmiss(nrmax, nmod)        Missing data
c     dimension rneff(nrmax, nmod)        Effective sample size
c
      integer ngptmax(nmt, nrmax, nmod), ngptact(nmt, nrmax, nmod)
      integer miss1(nrmax, nmod), iqcrit(nrmax, nmod)
      dimension pqcrit(nrmax, nmod), fcover(nmt, nrmax, nmod)
      dimension qnorm1(nmt, nrmax, nmod), qnorm2(nrmax, nmod)
      dimension cxx(0:lag), rxx(0:lag), result(31), res(nmt)
      dimension ybar(nrmax, nmod), sdy(nrmax, nmod), r1(nrmax, nmod)
      dimension trend(nrmax, nmod), tlc(nrmax, nmod)
      dimension serr1(nrmax, nmod), serr2(nrmax, nmod)
      dimension ptt(nrmax, nmod)
      dimension rmiss(nrmax, nmod), rneff(nrmax, nmod)
c
c-------------------------------------------------------------------------------
c     ** Coordinate information **
c
c     dimension rlon1(nx)                 Longitude coordinates
c     dimension rlat1(ny)                 Latitude coordinates
c     dimension rtime1(nmt)               Time coordinates (current model)
c     dimension rtime2(nmt, nrmax, nmod)  Time coordinates (all models)
c     dimension rtime3(nmt)               Time coordinates (analysis period)
c     dimension ctime(1)                  Current time coordinate
c
      dimension rlon1(nx), rlat1(ny)
      dimension rtime1(nmt), rtime2(nmt, nrmax, nmod), rtime3(nmt)
      dimension ctime(1)
c
c-------------------------------------------------------------------------------
c     ** Character variables **
c
c     character *42 p1(nmod)           First part of path to input data
c     character *42 p1a, p1b           Path to index data
c     character* 80 f0                 First part of path to data
c     character*120 f1(ntot)           Full path to dataset
c     character*120 f2(nrmax, nmod)    Full path to dataset
c     character*120 file3              Full path to land/sea mask
c     character *31 p10(3)             Special mask cases
c     character  *8 yweight1(nmod)     Weights
c
      character p1(nmod)*42
      character p1a*42, p1b*42
      character f0*80, f1(ntot)*120, f2(nrmax, nmod)*120
      character file3*120
      character p10(3)*31
      character yweight1(nmod)*8
c
c-------------------------------------------------------------------------------
c     ** Other character variables **
c
c     character *19 p2a(nr1)           Path to index data
c     character *19 p2b(nr2)           Path to index data
c     character *19 p2c(nr3)           Path to index data
c     character *19 p2d(nr4)           Path to index data
c     character *19 p2e(nr5)           Path to index data
c     character *19 p2f(nr6)           Path to index data
c     character *19 p2g(nr7)           Path to index data
c     character *19 p2h(nr8)           Path to index data
c     character *19 p2i(nr9)           Path to index data
c     character *19 p2j(nr10)          Path to index data
c     character *19 p2k(nr11)          Path to index data
c     character *19 p2l(nr12)          Path to index data
c     character *19 p2m(nr13)          Path to index data
c     character *19 p2n(nr14)          Path to index data
c     character *19 p2o(nr15)          Path to index data
c     character *19 p2p(nr16)          Path to index data
c     character *19 p2q(nr17)          Path to index data
c     character *19 p2r(nr18)          Path to index data
c     character *19 p2s(nr19)          Path to index data
c     character *19 p2t(nr20)          Path to index data
c     character  *4 p3                 Path to index data
c     character *23 p4a                Path to index data
c     character *23 v                  Path to index data
c     character *19 p5                 Path to index data
c     character *10 p6a                Path to index data
c     character *19 p7(nrmax, nmod)    Path to index data
c     character *19 p8                 Path to index data
c     character *85 fo1                String
c     character *85 fo2                String
c     character  *4 cyr1, cyr2         String
c     character  *4 cyr3, cyr4         String
c     character *22 cyr                String
c     character *14 s1                 String
c     character *12 s2a                String
c     character *12 s2b                String
c     character *35 s2                 String
c     character  *7 p9                 String
c     character *40 units1             Units of variable
c
      character p2a(nr1)*19, p2b(nr2)*19, p2c(nr3)*19
      character p2d(nr4)*19, p2e(nr5)*19, p2f(nr6)*19
      character p2g(nr7)*19, p2h(nr8)*19, p2i(nr9)*19
      character p2j(nr10)*19, p2k(nr11)*19, p2l(nr12)*19
      character p2m(nr13)*19, p2n(nr14)*19, p2o(nr15)*19
      character p2p(nr16)*19, p2q(nr17)*19, p2r(nr18)*19
      character p2s(nr19)*19, p2t(nr19)*19
      character p3*4, p4a*23, v*23, p5*19, p6a*10
      character p7(nrmax, nmod)*19, p8*19
      character fo1*85, fo2*85
      character cyr1*4, cyr2*4, cyr3*4, cyr4*4, cyr*22
      character s1*14, s2a*12, s2b*12, s2*35
      character p9*7
      character units1*40
c
c-------------------------------------------------------------------------------
c     ** Data statements **
c
      data s1 / 'ex20c3m_20mod_' /
c
      data x1(1),  x2(1),  y1(1),  y2(1)  / -180.,  180., -15., 15.  / ! TO15         
      data x1(2),  x2(2),  y1(2),  y2(2)  / -180.,  180., -20., 20.  / ! TO20         
      data x1(3),  x2(3),  y1(3),  y2(3)  / -180.,  180., -30., 30.  / ! TO30         
c
      data x1(4),  x2(4),  y1(4),  y2(4)  / -180.,  180., -15., 15.  / ! TL15         
      data x1(5),  x2(5),  y1(5),  y2(5)  / -180.,  180., -20., 20.  / ! TL20         
      data x1(6),  x2(6),  y1(6),  y2(6)  / -180.,  180., -30., 30.  / ! TL30         
c
      data x1(7),  x2(7),  y1(7),  y2(7)  / -180.,  180., -15., 15.  / ! TL01         
      data x1(8),  x2(8),  y1(8),  y2(8)  / -180.,  180., -20., 20.  / ! TL02         
      data x1(9),  x2(9),  y1(9),  y2(9)  / -180.,  180., -30., 30.  / ! TL03         
c
      data x1(10), x2(10), y1(10), y2(10) / -170., -120.,  -5.,  5.  / ! NI34
      data x1(11), x2(11), y1(11), y2(11) / -180.,  180., -50., 50.  / ! OC50
c
      data x1(12), x2(12), y1(12), y2(12) / -180.,  180., -90.,  0.  / ! SHEM
      data x1(13), x2(13), y1(13), y2(13) / -180.,  180.,   0., 90.  / ! NHEM
      data x1(14), x2(14), y1(14), y2(14) / -180.,  180., -90., 90.  / ! GLBL
      data x1(15), x2(15), y1(15), y2(15) / -180.,  180., -90., 90.  / ! GLAN
      data x1(16), x2(16), y1(16), y2(16) / -180.,  180., -90., 90.  / ! GOCN
c
      data x1(17), x2(17), y1(17), y2(17) / -180.,  180., -90.,-20.  / ! SHEX
      data x1(18), x2(18), y1(18), y2(18) / -180.,  180.,  20., 90.  / ! NHEX
c
      data x1(19), x2(19), y1(19), y2(19) / -180.,  180.,   0., 10.  / ! NHM1
      data x1(20), x2(20), y1(20), y2(20) / -180.,  180.,  10., 30.  / ! NHM2
      data x1(21), x2(21), y1(21), y2(21) / -180.,  180.,  30., 50.  / ! NHM3
      data x1(22), x2(22), y1(22), y2(22) / -180.,  180.,  50., 70.  / ! NHM4
      data x1(23), x2(23), y1(23), y2(23) / -180.,  180.,  70., 82.5 / ! NHM5
c
      data x1(24), x2(24), y1(24), y2(24) / -180.,  180.,   0.,-10.  / ! SHM1
      data x1(25), x2(25), y1(25), y2(25) / -180.,  180., -10.,-30.  / ! SHM2
      data x1(26), x2(26), y1(26), y2(26) / -180.,  180., -30.,-50.  / ! SHM3
      data x1(27), x2(27), y1(27), y2(27) / -180.,  180., -50.,-70.  / ! SHM4
      data x1(28), x2(28), y1(28), y2(28) / -180.,  180., -70.,-82.5 / ! SHM5
c
      data x1(29), x2(29), y1(29), y2(29) / -180.,  180.,   0., 82.5 / ! NHM6
      data x1(30), x2(30), y1(30), y2(30) / -180.,  180.,   0.,-70.  / ! SHM6
      data x1(31), x2(31), y1(31), y2(31) / -180.,  180.,   0.,-82.5 / ! SHM7
      data x1(32), x2(32), y1(32), y2(32) / -180.,  180., -82.5, 82.5 / ! GLB2
      data x1(33), x2(33), y1(33), y2(33) / -180.,  180., -70., 82.5 / ! GLB3
c
      data bogus / 1.0e20 /
c
      data bogus1 / -99.99999 /
c
      data bogus2 / -99.00000 /
c
      data units1 / 'Degrees C   ' /
c
      data p1 / '/pcmdi/bsanter1/Syn01/bccr_bcm2_0/        ', !  1(a)
     &          '/pcmdi/bsanter1/Syn01/cccma_cgcm3_1/      ', !  2(b)
     &          '/pcmdi/bsanter1/Syn01/cccma_cgcm3_1_t63/  ', !  3(c)
     &          '/pcmdi/bsanter1/Syn01/csiro_mk3_5/        ', !  4(d)
     &          '/pcmdi/bsanter1/Syn01/gfdl_cm2_0/         ', !  5(e)
     &          '/pcmdi/bsanter1/Syn01/gfdl_cm2_1/         ', !  6(f)
     &          '/pcmdi/bsanter1/Syn01/giss_aom/           ', !  7(g)
     &          '/pcmdi/bsanter1/Syn01/giss_model_e_h/     ', !  8(h)
     &          '/pcmdi/bsanter1/Syn01/giss_model_e_r/     ', !  9(i)
     &          '/pcmdi/bsanter1/Syn01/iap_fgoals1_0_g/    ', ! 10(j)
     &          '/pcmdi/bsanter1/Syn01/ingv_echam4/        ', ! 11(k)
     &          '/pcmdi/bsanter1/Syn01/inmcm3_0/           ', ! 12(l)
     &          '/pcmdi/bsanter1/Syn01/miroc3_2_hires/     ', ! 13(m)
     &          '/pcmdi/bsanter1/Syn01/miroc3_2_medres/    ', ! 14(n)
     &          '/pcmdi/bsanter1/Syn01/mpi_echam5/         ', ! 15(o)
     &          '/pcmdi/bsanter1/Syn01/mri_cgcm2_3_2a/     ', ! 16(p)
     &          '/pcmdi/bsanter1/Syn01/ncar_ccsm3_0/       ', ! 17(q)
     &          '/pcmdi/bsanter1/Syn01/ncar_pcm1/          ', ! 18(r)
     &          '/pcmdi/bsanter1/Syn01/ukmo_hadcm3/        ', ! 19(s)
     &          '/pcmdi/bsanter1/Syn01/ukmo_hadgem1/       '/ ! 20(t)
c
      data yweight1 / 'gaussian',     !  1 bccr_bcm2_0
     &                'gaussian',     !  2 cccma_cgcm3_1
     &                'gaussian',     !  3 cccma_cgcm3_1_t63
     &                'gaussian',     !  4 csiro_mk3_5 
     &                'cosine  ',     !  5 gfdl_cm2_0  
     &                'cosine  ',     !  6 gfdl_cm2_1    
     &                'cosine  ',     !  7 giss_aom 
     &                'cosine  ',     !  8 giss_model_e_h
     &                'cosine  ',     !  9 giss_model_e_r
     &                'cosine  ',     ! 10 iap_fgoals1_0_g
     &                'gaussian',     ! 11 ingv_echam4     
     &                'cosine  ',     ! 12 inmcm3_0
     &                'gaussian',     ! 13 miroc3_2_hires
     &                'gaussian',     ! 14 miroc3_2_medres
     &                'gaussian',     ! 15 mpi_echam5    
     &                'gaussian',     ! 16 mri_cgcm2_3_2a
     &                'gaussian',     ! 17 ncar_ccsm_3_0
     &                'gaussian',     ! 18 ncar_pcm1    
     &                'cosine  ',     ! 19 ukmo_hadcm3
     &                'cosine  '/     ! 20 ukmo_hadgem1
c
      data p2a/ 'ex20c3m.run1.mo.nc '/ ! bccr_bcm2_0
c
      data p2b/ 'ex20c3m.run1.mo.nc ', ! cccma_cgcm3_1
     &          'ex20c3m.run2.mo.nc ', ! cccma_cgcm3_1
     &          'ex20c3m.run3.mo.nc ', ! cccma_cgcm3_1
     &          'ex20c3m.run4.mo.nc ', ! cccma_cgcm3_1
     &          'ex20c3m.run5.mo.nc '/ ! cccma_cgcm3_1
c
      data p2c/ 'ex20c3m.run1.mo.nc '/ ! cccma_cgcm3_1_t63
c
      data p2d/ 'ex20c3m.run1.mo.nc '/ ! csiro_mk3_5 
c
      data p2e/ 'ex20c3m.run1.mo.nc '/ ! gfdl_cm2_0  
c
      data p2f/ 'ex20c3m.run1.mo.nc '/ ! gfdl_cm2_1   
c
      data p2g/ 'ex20c3m.run1.mo.nc ', ! giss_aom 
     &          'ex20c3m.run2.mo.nc '/ ! giss_aom 
c
      data p2h/ 'ex20c3m.run1.mo.nc ', ! giss_model_e_h
     &          'ex20c3m.run2.mo.nc ', ! giss_model_e_h
     &          'ex20c3m.run3.mo.nc '/ ! giss_model_e_h
c
      data p2i/ 'ex20c3m.run1.mo.nc ', ! giss_model_e_r
     &          'ex20c3m.run2.mo.nc ', ! giss_model_e_r
     &          'ex20c3m.run3.mo.nc ', ! giss_model_e_r
     &          'ex20c3m.run4.mo.nc ', ! giss_model_e_r
     &          'ex20c3m.run5.mo.nc '/ ! giss_model_e_r
c
      data p2j/ 'ex20c3m.run1.mo.nc ', ! iap_fgoals1_0_g
     &          'ex20c3m.run2.mo.nc ', ! iap_fgoals1_0_g
     &          'ex20c3m.run3.mo.nc '/ ! iap_fgoals1_0_g
c
      data p2k/ 'ex20c3m.run1.mo.nc '/ ! ingv_echam4
c
      data p2l/ 'ex20c3m.run1.mo.nc '/ ! inmcm3_0
c
      data p2m/ 'ex20c3m.run1.mo.nc '/ ! miroc3_2_hires
c
      data p2n/ 'ex20c3m.run1.mo.nc ', ! miroc3_2_medres
     &          'ex20c3m.run2.mo.nc ', ! miroc3_2_medres
     &          'ex20c3m.run3.mo.nc '/ ! miroc3_2_medres
c
      data p2o/ 'ex20c3m.run1.mo.nc ', ! mpi_echam5  
     &          'ex20c3m.run2.mo.nc ', ! mpi_echam5  
     &          'ex20c3m.run3.mo.nc ', ! mpi_echam5  
     &          'ex20c3m.run4.mo.nc '/ ! mpi_echam5  
c  
      data p2p/ 'ex20c3m.run1.mo.nc ', ! mri_cgcm2_3_2a
     &          'ex20c3m.run2.mo.nc ', ! mri_cgcm2_3_2a
     &          'ex20c3m.run3.mo.nc ', ! mri_cgcm2_3_2a
     &          'ex20c3m.run4.mo.nc ', ! mri_cgcm2_3_2a
     &          'ex20c3m.run5.mo.nc '/ ! mri_cgcm2_3_2a
c
      data p2q/ 'ex20c3m.run1.mo.nc ', ! ncar_ccsm_3_0
     &          'ex20c3m.run2.mo.nc ', ! ncar_ccsm_3_0
     &          'ex20c3m.run3.mo.nc ', ! ncar_ccsm_3_0
     &          'ex20c3m.run5.mo.nc ', ! ncar_ccsm_3_0
     &          'ex20c3m.run6.mo.nc ', ! ncar_ccsm_3_0
     &          'ex20c3m.run7.mo.nc ', ! ncar_ccsm_3_0
     &          'ex20c3m.run9.mo.nc '/ ! ncar_ccsm_3_0
c
      data p2r/ 'ex20c3m.run1.mo.nc ', ! ncar_pcm1
     &          'ex20c3m.run2.mo.nc ', ! ncar_pcm1
     &          'ex20c3m.run3.mo.nc ', ! ncar_pcm1
     &          'ex20c3m.run4.mo.nc '/ ! ncar_pcm1
c
      data p2s/ 'ex20c3m.run1.mo.nc '/ ! ukmo_hadcm3
c
      data p2t/ 'ex20c3m.run1.mo.nc '/ ! ukmo_hadgem1
c
      data p3 / '/Xy/' /           
c
      data p4a/ '_mm_xy_wf_r0000_0000.nc' / ! String for model file names.
c
      data p10 / 'sftlf_giss_model_e_h_atmgrid.nc',
     &           'sftlf_giss_model_e_r_atmgrid.nc',
     &           'sftlf_ukmo_hadgem1_atmgrid.nc' /
c
c
c
c
c
c
c-------------------------------------------------------------------------------
c     ** 1.0 **
c
c     INITIAL PROCESSING.
c
c     First read batch inputs.
c-------------------------------------------------------------------------------
c
      ic = 0
c
      open(unit=1, file='BOIL', form='formatted', status='old')
      read(1, 4) creg1
      read(1, 4) pv1
      read(1, *) iyr1
      read(1, *) iyr2
      read(1, *) iyr3, iyr4
      read(1, 4) ctype  
      read(1, *) itail     
      read(1, *) ires      
      read(1, *) savcrit
      close(unit=1)
    4 format(a16)
c
c-------------------------------------------------------------------------------
c     ** 1.1 **
c
c     Write output strings for selected processing periods.
c
c     Note: In the case of the ex20c3m experiment, the user selects the 
c     processing period for calculating the reference period (iyr1, iyr2) and 
c     all time series statistics (iyr3, iyr4). Spatial averages are output for 
c     the entire length of the ex20c3m experiment. The first and last years of 
c     the ex20c3m experiment are stored in iyear1 and iyear2 (respectively).
c-------------------------------------------------------------------------------
c
      if (iyr1 .eq. 0) then
        cyr1(1:4) = '0000'
      else
        write(cyr1, '(i4)') iyr1
      end if
c
      if (iyr2 .eq. 0) then
        cyr2(1:4) = '0000'
      else
        write(cyr2, '(i4)') iyr2
      end if
c
      write(cyr3, '(i4)') iyr3
      write(cyr4, '(i4)') iyr4
c
      if (pv1(1:3) .eq. 'tls') pv2 = 'eqmsu_tls'
      if (pv1(1:3) .eq. 'tmt') pv2 = 'eqmsu_tmt'
      if (pv1(1:3) .eq. 'tlt') pv2 = 'eqmsu_tlt'
c
c-------------------------------------------------------------------------------
c     ** 1.2 **
c
c     Loop over the total number of models processed.
c-------------------------------------------------------------------------------
c
      do 3000 ll = 1, nmod ! Loop over number of models 
c
        if (ll .eq.  1) nreal = nr1
        if (ll .eq.  2) nreal = nr2
        if (ll .eq.  3) nreal = nr3
        if (ll .eq.  4) nreal = nr4
        if (ll .eq.  5) nreal = nr5
        if (ll .eq.  6) nreal = nr6
        if (ll .eq.  7) nreal = nr7
        if (ll .eq.  8) nreal = nr8
        if (ll .eq.  9) nreal = nr9
        if (ll .eq. 10) nreal = nr10
        if (ll .eq. 11) nreal = nr11
        if (ll .eq. 12) nreal = nr12
        if (ll .eq. 13) nreal = nr13
        if (ll .eq. 14) nreal = nr14
        if (ll .eq. 15) nreal = nr15
        if (ll .eq. 16) nreal = nr16
        if (ll .eq. 17) nreal = nr17
        if (ll .eq. 18) nreal = nr18
        if (ll .eq. 19) nreal = nr19
        if (ll .eq. 20) nreal = nr20
c
c-------------------------------------------------------------------------------
c     ** 1.3 **
c
c     Set path to current model dataset.               
c-------------------------------------------------------------------------------
c
        p1a = p1(ll)   
        p1b = p1a(23:42)
        print *, 'p1a = ', p1a
        print *, 'p1b = ', p1b
        ida = index(p1a, ' ') - 1
        idb = index(p1b, ' ') - 2
        idc = ida-1                
c
c-------------------------------------------------------------------------------
c     ** 1.4 **
c
c     Loop over the number of realizations per model.
c-------------------------------------------------------------------------------
c
        do 2000 kk = 1, nreal
c
          if (ll .eq.  1) p5 = p2a(kk)
          if (ll .eq.  2) p5 = p2b(kk)
          if (ll .eq.  3) p5 = p2c(kk)
          if (ll .eq.  4) p5 = p2d(kk)
          if (ll .eq.  5) p5 = p2e(kk)
          if (ll .eq.  6) p5 = p2f(kk)
          if (ll .eq.  7) p5 = p2g(kk)
          if (ll .eq.  8) p5 = p2h(kk)
          if (ll .eq.  9) p5 = p2i(kk)
          if (ll .eq. 10) p5 = p2j(kk)
          if (ll .eq. 11) p5 = p2k(kk)
          if (ll .eq. 12) p5 = p2l(kk)
          if (ll .eq. 13) p5 = p2m(kk)
          if (ll .eq. 14) p5 = p2n(kk)
          if (ll .eq. 15) p5 = p2o(kk)
          if (ll .eq. 16) p5 = p2p(kk)
          if (ll .eq. 17) p5 = p2q(kk)
          if (ll .eq. 18) p5 = p2r(kk)
          if (ll .eq. 19) p5 = p2s(kk)
          if (ll .eq. 20) p5 = p2t(kk)
c
          p7(kk, ll) = p5 ! Store run number.
c
c-------------------------------------------------------------------------------
c     ** 1.5 **
c
c     Calculate various character string lengths.    
c-------------------------------------------------------------------------------
c
          ic   = ic + 1
          id1  = index(p5, ' ') - 1
          p6a  = pv1
          id2a = index(p6a, ' ') - 1
          print *, 'p5 = ', p5
          print *, 'p6a = ', p6a
c
c-------------------------------------------------------------------------------
c     ** 1.6 **
c
c     Determine which geographical region has been selected.      
c-------------------------------------------------------------------------------
c
          if (creg1(1:4) .eq. 'TO15') then ! Process variable.
            idd = 1
          else if (creg1(1:4) .eq. 'TO20') then
            idd = 2
          else if (creg1(1:4) .eq. 'TO30') then
            idd = 3
          else if (creg1(1:4) .eq. 'TL15') then
            idd = 4
          else if (creg1(1:4) .eq. 'TL20') then
            idd = 5
          else if (creg1(1:4) .eq. 'TL30') then
            idd = 6
          else if (creg1(1:4) .eq. 'TLO1') then
            idd = 7
          else if (creg1(1:4) .eq. 'TLO2') then
            idd = 8
          else if (creg1(1:4) .eq. 'TLO3') then
            idd = 9
          else if (creg1(1:4) .eq. 'NI34') then
            idd = 10
          else if (creg1(1:4) .eq. 'OC50') then
            idd = 11
          else if (creg1(1:4) .eq. 'SHEM') then
            idd = 12
          else if (creg1(1:4) .eq. 'NHEM') then
            idd = 13
          else if (creg1(1:4) .eq. 'GLBL') then
            idd = 14
          else if (creg1(1:4) .eq. 'GLAN') then
            idd = 15
          else if (creg1(1:4) .eq. 'GOCN') then
            idd = 16
          else if (creg1(1:4) .eq. 'SHEX') then
            idd = 17
          else if (creg1(1:4) .eq. 'NHEX') then
            idd = 18
          else if (creg1(1:4) .eq. 'NHM1') then
            idd = 19
          else if (creg1(1:4) .eq. 'NHM2') then
            idd = 20
          else if (creg1(1:4) .eq. 'NHM3') then
            idd = 21
          else if (creg1(1:4) .eq. 'NHM4') then
            idd = 22
          else if (creg1(1:4) .eq. 'NHM5') then
            idd = 23
          else if (creg1(1:4) .eq. 'SHM1') then
            idd = 24
          else if (creg1(1:4) .eq. 'SHM2') then
            idd = 25
          else if (creg1(1:4) .eq. 'SHM3') then
            idd = 26
          else if (creg1(1:4) .eq. 'SHM4') then
            idd = 27
          else if (creg1(1:4) .eq. 'SHM5') then
            idd = 28
          else if (creg1(1:4) .eq. 'NHM6') then
            idd = 29
          else if (creg1(1:4) .eq. 'SHM6') then
            idd = 30
          else if (creg1(1:4) .eq. 'SHM7') then
            idd = 31
          else if (creg1(1:4) .eq. 'GLB2') then
            idd = 32
          else if (creg1(1:4) .eq. 'GLB3') then
            idd = 33
          else
            print *, 'Error in creg1!'
            print *, 'creg = ', creg1(1:4)
            stop
          end if
c
c-------------------------------------------------------------------------------
c     ** 1.7 **
c
c     Set path names to model data files. Note that the giss_model_e_h, 
c     giss_model_e_r, and ukmo_hadgem models are a special case, since 
c     atmospheric variables (like ta) are on a different grid than surface 
c     variables (like ts).
c-------------------------------------------------------------------------------
c
          f0  = p1a(1:ida)//'ex20c3m/'//p6a(1:id2a)//'_both_'
          idf = index(f0, ' ') - 1
c
          f1(ic)     = f0(1:idf)//p1b(1:idb)//'.'//p5(1:id1)
          f2(kk, ll) = f0(1:idf)//p1b(1:idb)//'.'//p5(1:id1)
c
          if (pv2(1:5) .eq. 'eqmsu') then
            if (p1a(23:idc) .eq. 'giss_model_e_h') then
              file3 = p1a(1:ida)//'Masks/'//p10(1)
            else if (p1a(23:idc) .eq. 'giss_model_e_r') then
              file3 = p1a(1:ida)//'Masks/'//p10(2)
            else if (p1a(23:idc) .eq. 'ukmo_hadgem1') then
              file3 = p1a(1:ida)//'Masks/'//p10(3)
            else
              file3 = p1a(1:ida)//'Masks/sftlf_'//p1a(23:idc)//'.nc'
            end if
          else
            file3 = p1a(1:ida)//'Masks/sftlf_'//p1a(23:idc)//'.nc'
          end if
c
          print *, 'THE FOLLOWING INPUT DATASETS WILL BE USED:'
          print *, 'DATASET #: ', ic
          print *, 'f1 = ', f1(ic)
          print *, 'f2 = ', f2(kk, ll)
          print *, 'file3 = ', file3
c
c-------------------------------------------------------------------------------
c     ** 1.8 **
c
c     Retrieve land/sea mask.
c-------------------------------------------------------------------------------
c
          call rdmask(nx, ny, nx1, ny1, x1(idd), x2(idd),
     &    y1(idd), y2(idd), file3, work, fmask, yweight1(ll))
c
          print *, 'Now in Section 1.8!'
          print *, 'LON/LAT dimensions, mask dataset: ', nx1, ny1
c
c-------------------------------------------------------------------------------
c     ** 1.9 **
c
c     Determine the time dimension of the data. 
c-------------------------------------------------------------------------------
c
          print *, 'Now in Section 1.9!'
          print *, 'Calling SHAPE!'
          call initget
          call defdim(1, 1, 'longitude', 'width', 'nearest', -180.0,
     &    180.0, 360.0)
c
          call defdim(1, 2, 'latitude', yweight1(ll), 'range', -90.0,
     &    90.0, 0.0)
c
          call defdim(1, 3, 'time', 'unit', 'assaved', 1.0, 1.0, 0.0)
          call defmisc('maximum dimension kept', 'integer', 12012)
c
          call defvar(1, pv2, f1(ic))
          call shape(1, ldim1, ldim2, ldim3, ldim4, isize)
          call getcoord(1, 3, rtime1)
          call closeget
c
          nmt1 = ldim3
          print *, 'After SHAPE!'
          print *, 'Dimension 1 of retrieved array: ', ldim1
          print *, 'Dimension 2 of retrieved array: ', ldim2
          print *, 'Dimension 3 of retrieved array: ', ldim3
          print *, 'Dimension 4 of retrieved array: ', ldim4
          print *, 'nmt, nmt1 = ', nmt, nmt1
c
          ileng1(kk, ll) = nmt1 ! Store retrieved length of data series.
c
          do 190 mm = 1, nmt1
            rtime2(mm, kk, ll) = rtime1(mm) ! Store retrieved time.
  190     continue
c
c
c
c
c
c
c-------------------------------------------------------------------------------
c     ** 2.0 **
c
c     DATA RETRIEVAL.
c
c     Retrieve temperature data, looping over total number of months processed 
c     (nmt1). 
c
c     Note: isumm1 is (for each model and ex20c3m realization) the total number 
c     of months with missing data.
c-------------------------------------------------------------------------------
c
          print *, 'Now in Section 2.0!'
c
          isumm1 = 0 ! Counter for total number of  months with missing data.
          do 1000 mm = 1, nmt1 
c
            call rd2d_xy(nx, ny, nx1, ny1, mm, x1(idd), x2(idd), 
     &      y1(idd), y2(idd), f1(ic), rlon1, rlat1, ctime, bmask1, 
     &      c0, pv2, yweight1(ll))
c
c-------------------------------------------------------------------------------
c     ** 2.1 **
c
c     Convert model atmospheric temperature data from Kelvin to degrees Celsius.
c-------------------------------------------------------------------------------
c
            do 201 j = 1, ny1
              do 200 i = 1, nx1
                if (c0(i ,j) .ne. bogus) then
                  c0(i, j) = c0(i, j) - 273.15
                else
                  c0(i, j) = bogus
                end if
  200         continue
  201       continue
c
            if ((ll .eq. nmod) .and. (mm .eq. 1)) then
              print *, 'DIAGNOSTICS!'
              do 203 j = 1, ny1
                do 202 i = 1, nx1
c                 print *, 'i, j, c0 = ', i, j, c0(i, j) 
  202           continue
  203         continue
            end if
c
            call clrtable
c
c-------------------------------------------------------------------------------
c     ** 2.2 **
c
c     Compute spatial means. 
c
c     For "land only" spatial averages, mask out all ocean points in region.
c     Note: If data are missing, bmask is zero.
c-------------------------------------------------------------------------------
c
            sum1 = 0.0
            sum2 = 0.0
            ixy1 = 0
            do 222 j = 1, ny1
              do 221 i = 1, nx1
                term1 = c0(i, j)
                term2 = bmask1(i, j)        ! Area weights for data
                term3 = fmask(i, j) / 100.0 ! Land percentage / 100.0
c
                if (term1 .ne. bogus) then
                  ixy1 = ixy1 + 1
c
c-------------------------------------------------------------------------------
c     ** 2.3 **
c
c     CASE 1: "Ocean only" spatial averages, mask out all land points in region.
c-------------------------------------------------------------------------------
c
                  if (creg1(1:4) .eq. 'TO15') then
                    sum1 = sum1 + (term1 * term2 * (1.0 - term3))
                    sum2 = sum2 + (term2 * (1.0 - term3))
                  else if (creg1(1:4) .eq. 'TO20') then
                    sum1 = sum1 + (term1 * term2 * (1.0 - term3))
                    sum2 = sum2 + (term2 * (1.0 - term3))
                  else if (creg1(1:4) .eq. 'TO30') then
                    sum1 = sum1 + (term1 * term2 * (1.0 - term3))
                    sum2 = sum2 + (term2 * (1.0 - term3))
                  else if (creg1(1:4) .eq. 'NI34') then
                    sum1 = sum1 + (term1 * term2 * (1.0 - term3))
                    sum2 = sum2 + (term2 * (1.0 - term3))
                  else if (creg1(1:4) .eq. 'OC50') then
                    sum1 = sum1 + (term1 * term2 * (1.0 - term3))
                    sum2 = sum2 + (term2 * (1.0 - term3))
                  else if (creg1(1:4) .eq. 'GOCN') then
                    sum1 = sum1 + (term1 * term2 * (1.0 - term3))
                    sum2 = sum2 + (term2 * (1.0 - term3))
c
c-------------------------------------------------------------------------------
c     ** 2.4 **
c
c     CASE 2: "Land only" spatial averages, mask out all ocean points in region.
c-------------------------------------------------------------------------------
c
                  else if (creg1(1:4) .eq. 'TL15') then
                    sum1 = sum1 + (term1 * term2 * term3)
                    sum2 = sum2 + (term2 * term3)
                  else if (creg1(1:4) .eq. 'TL20') then
                    sum1 = sum1 + (term1 * term2 * term3)
                    sum2 = sum2 + (term2 * term3)
                  else if (creg1(1:4) .eq. 'TL30') then
                    sum1 = sum1 + (term1 * term2 * term3)
                    sum2 = sum2 + (term2 * term3)
                  else if (creg1(1:4) .eq. 'GLAN') then
                    sum1 = sum1 + (term1 * term2 * term3)
                    sum2 = sum2 + (term2 * term3)
c
c-------------------------------------------------------------------------------
c     ** 2.5 **
c
c     CASE 3: "Land + ocean" spatial averages, no masking.
c-------------------------------------------------------------------------------
c
                  else if (creg1(1:4) .eq. 'TLO1') then
                    sum1 = sum1 + (term1 * term2)
                    sum2 = sum2 + term2
                  else if (creg1(1:4) .eq. 'TLO2') then
                    sum1 = sum1 + (term1 * term2)
                    sum2 = sum2 + term2
                  else if (creg1(1:4) .eq. 'TLO3') then
                    sum1 = sum1 + (term1 * term2)
                    sum2 = sum2 + term2
                  else if (creg1(1:4) .eq. 'NHEM') then
                    sum1 = sum1 + (term1 * term2)
                    sum2 = sum2 + term2
                  else if (creg1(1:4) .eq. 'SHEM') then
                    sum1 = sum1 + (term1 * term2)
                    sum2 = sum2 + term2
                  else if (creg1(1:4) .eq. 'GLBL') then
                    sum1 = sum1 + (term1 * term2)
                    sum2 = sum2 + term2
                  else if (creg1(1:4) .eq. 'NHEX') then
                    sum1 = sum1 + (term1 * term2)
                    sum2 = sum2 + term2
                  else if (creg1(1:4) .eq. 'SHEX') then
                    sum1 = sum1 + (term1 * term2)
                    sum2 = sum2 + term2
                  else if (creg1(1:4) .eq. 'NHM1') then
                    sum1 = sum1 + (term1 * term2)
                    sum2 = sum2 + term2
                  else if (creg1(1:4) .eq. 'NHM2') then
                    sum1 = sum1 + (term1 * term2)
                    sum2 = sum2 + term2
                  else if (creg1(1:4) .eq. 'NHM3') then
                    sum1 = sum1 + (term1 * term2)
                    sum2 = sum2 + term2
                  else if (creg1(1:4) .eq. 'NHM4') then
                    sum1 = sum1 + (term1 * term2)
                    sum2 = sum2 + term2
                  else if (creg1(1:4) .eq. 'NHM5') then
                    sum1 = sum1 + (term1 * term2)
                    sum2 = sum2 + term2
                  else if (creg1(1:4) .eq. 'NHM6') then
                    sum1 = sum1 + (term1 * term2)
                    sum2 = sum2 + term2
                  else if (creg1(1:4) .eq. 'SHM1') then
                    sum1 = sum1 + (term1 * term2)
                    sum2 = sum2 + term2
                  else if (creg1(1:4) .eq. 'SHM2') then
                    sum1 = sum1 + (term1 * term2)
                    sum2 = sum2 + term2
                  else if (creg1(1:4) .eq. 'SHM3') then
                    sum1 = sum1 + (term1 * term2)
                    sum2 = sum2 + term2
                  else if (creg1(1:4) .eq. 'SHM4') then
                    sum1 = sum1 + (term1 * term2)
                    sum2 = sum2 + term2
                  else if (creg1(1:4) .eq. 'SHM5') then
                    sum1 = sum1 + (term1 * term2)
                    sum2 = sum2 + term2
                  else if (creg1(1:4) .eq. 'SHM6') then
                    sum1 = sum1 + (term1 * term2)
                    sum2 = sum2 + term2
                  else if (creg1(1:4) .eq. 'SHM7') then
                    sum1 = sum1 + (term1 * term2)
                    sum2 = sum2 + term2
                  else if (creg1(1:4) .eq. 'GLB2') then
                    sum1 = sum1 + (term1 * term2)
                    sum2 = sum2 + term2
                  else if (creg1(1:4) .eq. 'GLB3') then
                    sum1 = sum1 + (term1 * term2)
                    sum2 = sum2 + term2
                  else
                    print *, 'Error in creg1!'
                    print *, 'creg = ', creg1(1:4)
                    stop
                  end if
                end if
  221         continue
  222       continue
c
            trm1   = ((nx1 * 1.0) * (ny1 * 1.0)) * savcrit
            jcrit1 = nint(trm1)
            if (mm .eq. 1) print *, 'trm1, jcrit1 = ', trm1, jcrit1
c
            if (ixy1 .gt. jcrit1) then ! Test for adequate no. of valid points.
              c2(mm, kk, ll) = sum1 / sum2 ! Calculate spatial average.
            else
              c2(mm, kk, ll) = bogus1 ! Set spatial average to missing.
              isumm1 = isumm1 + 1
            end if
c
            iqq1 = ny1 * nx1
            iqq2 = ixy1
            ngptmax(mm, kk, ll) = ny1 * nx1  
            ngptact(mm, kk, ll) = ixy1       
            fcover(mm, kk, ll)  = (iqq2 * 1.0) / (iqq1 * 1.0)
 1000     continue ! End loop over time.
c
          print *, 'c2(nmt1, kk, ll) = ', c2(nmt1, kk, ll)
          miss1(kk, ll) = isumm1
c
c
c
c
c
c
c-------------------------------------------------------------------------------
c     ** 3.0 **
c
c     BASIC QUALITY CONTROL.
c
c     First calculate mean and temporal standard deviation of current model
c     time series of first differences (FDs). In calculating these quantities,
c     it is important to exclude missing data.
c
c     Note: The -1 in the calculation of tcrit1 arises because the 300 loop for
c     calculation of first differences is over nmt1 - 1, not over nmt1.
c-------------------------------------------------------------------------------
c
          iq1   = 0   
          sumq1 = 0.0
          do 300 mm = 2, nmt1 ! Loop over number of months.
            ttt1 = c2(mm, kk, ll)
            ttt2 = c2(mm - 1, kk, ll)
            if ((ttt1 .gt. bogus2) .and. (ttt2 .gt. bogus2)) then
              termq1 = abs(ttt1 - ttt2) ! Absolute value of FDs.
              iq1    = iq1 + 1
              sumq1  = sumq1 + termq1
            end if
  300     continue ! End loop over time.
c
          tcrit1 = (nmt1 - 1.0) * 0.5 ! Half length of input time series (-1).
          tcrit2 = iq1 * 1.0
          if (tcrit2 .gt. tcrit1) then
            aveq1 = sumq1 / (iq1 * 1.0) ! Average of FD absolute values.
          else
            aveq1 = bogus1
          end if
c
          iq2   = 0   
          sumq2 = 0.0
          do 301 mm = 2, nmt1 ! Loop over number of months.
            ttt1 = c2(mm, kk, ll)
            ttt2 = c2(mm - 1, kk, ll)
            if ((ttt1 .gt. bogus2) .and. (ttt2 .gt. bogus2)) then
              termq2 = abs(ttt1 - ttt2) ! Absolute value of FDs.
              if (aveq1 .gt. bogus2) then
                iq2   = iq2 + 1
                sumq2 = sumq2 + ((termq2 - aveq1) ** 2) ! Sum of squares.
              end if
            end if
  301     continue ! End loop over time.
c
          tcrit3 = iq2 * 1.0
          if (tcrit3 .gt. tcrit1) then
            sdevq1 = sqrt(sumq2 / (iq2 - 1.0)) ! Standard deviation.
          else
            sdevq1 = bogus1                    
          end if
c
c-------------------------------------------------------------------------------
c     ** 3.1 **
c
c     Calculate first differences for current time series. Then normalize by the
c     temporal standard deviation of the first differences. This is qnorm1. The
c     time average of the qnorm1 values is qnorm2.
c
c     Note: In calculating qnorm2, the number of first differences is nmt1 - 1.
c-------------------------------------------------------------------------------
c
          iq3   = 0   
          iq4   = 0   
          sumq3 = 0.0
          qnorm1(1, kk, ll) = 0.0 ! Initialize qnorm1.
c
          do 310 mm = 2, nmt1 ! Loop over number of months.
            ttt1 = c2(mm, kk, ll)
            ttt2 = c2(mm - 1, kk, ll)
            if ((ttt1 .gt. bogus2) .and. (ttt2 .gt. bogus2)) then
              termq3 = abs(ttt1 - ttt2) ! Absolute value of FDs.
              if (sdevq1 .gt. bogus2) then
                iq3                = iq3 + 1
                qnorm1(mm, kk, ll) = termq3 / sdevq1
                sumq3              = sumq3 + qnorm1(mm, kk, ll)
                qq1                = qnorm1(mm, kk, ll)
                if (qq1 .gt. 5.0) then
                  iq4 = iq4 + 1
                end if
              end if
            end if
  310     continue ! End loop over time.
c
          qnorm2(kk, ll)  = sumq3 / (nmt1 - 1.0)
          iqcrit(kk, ll)  = iq4                     
          pqcrit(kk, ll)  = 100.0 * ((iq4 * 1.0) / (nmt1 - 1.0))
c
          print *, 'qnorm2(kk, ll) = ', kk, ll, qnorm2(kk, ll)
c
c
c
c
c
c
c-------------------------------------------------------------------------------
c     ** 4.0 **
c
c     COMPUTE CLIMATOLOGICAL MONTHLY MEANS AND ANOMALIES.
c
c     First determine and store, for each model and ex20c3m realization:
c
c     iyear1    Time coordinate of first year of complete data   
c     iyear2    Time coordinate of last year of complete data   
c     nyear1    Number of complete years of data
c
c     Note: jyr1 and jyr2 are (respectively) the indices of the first year and
c     the last complete year of data on the input file for the current model and
c     current ex20c3m realization. The variables iyy1 and iyy2 are the values of 
c     the first year and the last complete year of data on the input file for 
c     the current model.
c-------------------------------------------------------------------------------
c
          print *, 'Now in Section 4.0!'
          refcrit = 0.5
c
          call timecoord3(nmt, nmt1, rtime1, nyear1,
     &    jyr1, jyr2, iyy1, iyy2, refcrit, mcrit)
c
          iyear1(kk, ll) = iyy1
          iyear2(kk, ll) = iyy2
c
c-------------------------------------------------------------------------------
c     ** 4.1 **
c
c     Next calculate:                                                   
c
c     kyr1      Time coordinate of first year of reference period
c     kyr2      Time coordinate of last year of reference period
c
c     Do this only if a valid reference period (iyr1, iyr2 .ne. 0) has been
c     selected.
c
c     Note: The reference period is fixed. It does not vary from model to model
c     or from realization to realization.
c-------------------------------------------------------------------------------
c
          if ((iyr1 .ne. 0) .and. (iyr2 .ne. 0)) then
            call timecoord1(nmt, nmt1, rtime1, nyear1, iyr1, iyr2,
     &      kyr1, kyr2, refcrit, mcrit)
          end if
c
c-------------------------------------------------------------------------------
c     ** 4.2 **
c
c     Calculate anomalies.                        
c-------------------------------------------------------------------------------
c
          print *, 'Now in Section 4.2!'
          print *, 'kyr1, kyr2 = ', kyr1, kyr2
          print *, 'Compute climatological monthly means!'
          print *, 'nt, nyear1 = ', nt, nyear1
c
          imm = 0 ! Counter over total number of months processed.
          do 422 m = 1, nm
            do 420 n = 1, nyear1
              inn   = ((n - 1) * 12) + m
              xt(n) = c2(inn, kk, ll)
  420       continue
c
            call climanom1(nt, nyear1, kyr1, kyr2, mcrit, bogus1,
     &      xt, xm)
c
            sm1(m, kk, ll) = xm
c
            do 421 n = 1, nyear1
              imm             = imm + 1
              im              = ((n - 1) * 12) + m
              ts1(im, kk, ll) = xt(n)
  421       continue
  422     continue
c
          print *, 'ANOM! c2(1, kk, ll) = ', c2(1, kk, ll)
c
c-------------------------------------------------------------------------------
c     ** 4.3 **
c
c     Terminate loop over number of models and realizations.
c-------------------------------------------------------------------------------
c
 2000   continue ! End loop over number of ex20c3m realizations.
 3000 continue ! End loop over number of models.               
c
      ntot1 = ic
      print *, 'Total number of realizations processed = ', ntot1
      print *, 'Total number of time points processed  = ', imm  
c
c
c
c
c
c
c-------------------------------------------------------------------------------
c     ** 5.0 **
c
c     CALCULATION OF TIME SERIES STATISTICS FOR EACH MODEL AND EACH EX20C3M
c     REALIZATION.
c
c     First loop over the total number of models.
c-------------------------------------------------------------------------------
c
      do 5000 ll = 1, nmod
        if (ll .eq.  1) nreal = nr1
        if (ll .eq.  2) nreal = nr2
        if (ll .eq.  3) nreal = nr3
        if (ll .eq.  4) nreal = nr4
        if (ll .eq.  5) nreal = nr5
        if (ll .eq.  6) nreal = nr6
        if (ll .eq.  7) nreal = nr7
        if (ll .eq.  8) nreal = nr8
        if (ll .eq.  9) nreal = nr9
        if (ll .eq. 10) nreal = nr10
        if (ll .eq. 11) nreal = nr11
        if (ll .eq. 12) nreal = nr12
        if (ll .eq. 13) nreal = nr13
        if (ll .eq. 14) nreal = nr14
        if (ll .eq. 15) nreal = nr15
        if (ll .eq. 16) nreal = nr16
        if (ll .eq. 17) nreal = nr17
        if (ll .eq. 18) nreal = nr18
        if (ll .eq. 19) nreal = nr19
        if (ll .eq. 20) nreal = nr20
c
        sum1 = 0.0
        sum2 = 0.0
        sum3 = 0.0
        sum4 = 0.0
c
c-------------------------------------------------------------------------------
c     ** 5.1 **
c
c     Now calculate:
c     
c     myr1      Time coordinate, first month of first year of reference period
c     myr2      Time coordinate, last month of last year of reference period
c
c     Note: The reference period is fixed. It does not vary from model to model
c     or from realization to realization.
c-------------------------------------------------------------------------------
c
        do 4000 kk = 1, nreal ! Loop over number of ex20c3m realizations.
c
          rtt0   = rtime2(1, kk, ll) ! Set initial time for current dataset.
          rbase1 = ((rtt0 / 12.0) + 1800.0) - 1.0
          ibase1 = ifix(rbase1)
c
          myr1   = ((iyr3 - ibase1) * 12) - 11
          myr2   = (iyr4 - ibase1) * 12
          kcrit  = ((myr2 - myr1) + 1 ) / 2
c
          print *, 'rtime2     = ', rtime2(1, kk, ll)
          print *, 'rtt0       = ', rtt0   
          print *, 'rbase1     = ', rbase1
          print *, 'ibase1     = ', ibase1
          print *, 'kk, ll, myr1, myr2 = ', kk, ll, myr1, myr2
c
          kcc = 0 ! Initial counter for number of months with valid data.
c
c-------------------------------------------------------------------------------
c     ** 5.2 **
c
c     Copy time series for current model and ex20c3m realization into yt. Also
c     store this data (for the selected analysis period) in ts2. Note that
c     rtime3 contains time coordinate information for the selected analysis
c     period.
c-------------------------------------------------------------------------------
c
          do 520 mm = myr1, myr2 ! Loop over selected analysis period.
            term1 = ts1(mm, kk, ll)
            if (term1 .gt. bogus2) then
              kcc              = kcc + 1
              yt(kcc)          = term1
              zt(kcc)          = kcc * 1.0
              ts2(kcc, kk, ll) = term1
              rtime3(kcc)      = rtime2(mm, kk, ll)
            end if
  520     continue
c
c-------------------------------------------------------------------------------
c     ** 5.3 **
c
c     Calculate time series statistics for each model and realization.
c-------------------------------------------------------------------------------
c
          print *, 'kcc, kcrit = ', kcc, kcrit
c
          if (kcc .gt. kcrit) then
c
            call reg_ar0(nmt, kcc, lag, rnff, zt, yt, result, res,
     &      cxx, rxx, ires, bogus1)
c
            ybar(kk, ll)  = result(2)
            sdy(kk, ll)   = result(4)
            r1(kk, ll)    = rxx(1)
            trend(kk, ll) = result(6) * 120.0
            tlc(kk, ll)   = result(6) * (kcc * 1.0)
c
            if (result(8) .gt. bogus2) then
              serr1(kk, ll) = result(8) * 120.0
            else
              serr1(kk, ll) = bogus1             
            end if
c
            if (result(9) .gt. bogus2) then
              serr2(kk, ll) = result(9) * 120.0
            else
              serr2(kk, ll) = bogus1             
            end if
c
            if (result(28) .gt. bogus2) then
              if (itail .eq. 1) ptt(kk, ll) = result(28) / 2.0
              if (itail .eq. 2) ptt(kk, ll) = result(28)
            else
              ptt(kk, ll) = bogus1    
            end if 
c
            rmiss(kk, ll) = kcc * 1.0
            rneff(kk, ll) = rnff
          else
            ybar(kk, ll)  = bogus1   
            sdy(kk, ll)   = bogus1   
            r1(kk, ll)    = bogus1   
            trend(kk, ll) = bogus1
            tlc(kk, ll)   = bogus1
            serr2(kk, ll) = bogus1
            ptt(kk, ll)   = bogus1
            rmiss(kk, ll) = bogus1   
            rneff(kk, ll) = bogus1   
          end if
c
c-------------------------------------------------------------------------------
c     ** 5.4 **
c
c     For each model, accumulate sum over number of realizations of the values
c     of each of the saved time series statistics. Then calculate each
c     individual model's ensemble-mean value of the time series statistic.
c-------------------------------------------------------------------------------
c
          sum1 = sum1 + result(2)         ! STAT1: Mean
          sum2 = sum2 + result(4)         ! STAT2: Standard deviation
          sum3 = sum3 + rxx(1)            ! STAT3: Lag-1 autocorrelation
          sum4 = sum4 + result(6) * 120.0 ! STAT4: Trend
 4000   continue ! End loop over realizations.
c
        aaa = sum1 / (nreal * 1.0)
        bbb = sum2 / (nreal * 1.0)
        ccc = sum3 / (nreal * 1.0)
        ddd = sum4 / (nreal * 1.0)
c
c-------------------------------------------------------------------------------
c     ** 5.5 **
c
c     Store ensemble-mean values of time series statistics.          
c-------------------------------------------------------------------------------
c
        wk3(1, ll) = aaa ! STAT1: Mean
        wk3(2, ll) = bbb ! STAT2: Standard deviation
        wk3(3, ll) = ccc ! STAT3: Lag-1 autocorrelation
        wk3(4, ll) = ddd ! STAT4: Trend
        print *, 'wk3(1, ll) = ', wk3(1, ll), ll
c
 5000 continue ! End loop over models.
c
c-------------------------------------------------------------------------------
c     ** 5.6 **
c
c     Calculate the ensemble-mean time series for each individual model for
c     which multiple realizations of extended 20c3m runs are available.
c-------------------------------------------------------------------------------
c
      do 562 ll = 1, nmod ! Loop over number of models.
        if (ll .eq.  1) nreal = nr1
        if (ll .eq.  2) nreal = nr2
        if (ll .eq.  3) nreal = nr3
        if (ll .eq.  4) nreal = nr4
        if (ll .eq.  5) nreal = nr5
        if (ll .eq.  6) nreal = nr6
        if (ll .eq.  7) nreal = nr7
        if (ll .eq.  8) nreal = nr8
        if (ll .eq.  9) nreal = nr9
        if (ll .eq. 10) nreal = nr10
        if (ll .eq. 11) nreal = nr11
        if (ll .eq. 12) nreal = nr12
        if (ll .eq. 13) nreal = nr13
        if (ll .eq. 14) nreal = nr14
        if (ll .eq. 15) nreal = nr15
        if (ll .eq. 16) nreal = nr16
        if (ll .eq. 17) nreal = nr17
        if (ll .eq. 18) nreal = nr18
        if (ll .eq. 19) nreal = nr19
        if (ll .eq. 20) nreal = nr20
c
        do 561 mm = 1, kcc ! Loop over total number of months for selected 
c                          !  reference period.
          sum1 = 0.0
          if (nreal .gt. 1) then
            icnt = 0
            do 560 kk = 1, nreal ! Loop over number of realizations.
              term1 = ts2(mm, kk, ll)
              if (term1 .gt. bogus2) then
                icnt = icnt + 1
                sum1 = sum1 + term1
              end if
  560       continue
c
            if (icnt .ge. 2) then
              aaa = sum1 / (icnt * 1.0)
              em1(mm, ll) = aaa
            else
              em1(mm, ll) = bogus1
            end if
          else
            em1(mm, ll) = ts2(mm, 1, ll)
          end if
  561   continue
  562 continue
c
c-------------------------------------------------------------------------------
c     ** 5.7 **
c
c     Calculate multi-model average time series. This is the average of the
c     ensemble-mean time series of individual models.
c-------------------------------------------------------------------------------
c
      jcrit = nmod / 2
c
      do 571 mm = 1, kcc ! Loop over total number of months for selected
c                        !  reference period.
        iq = 0
        do 570 ll = 1, nmod ! Loop over number of models.
          term1 = em1(mm, ll)
          if (term1 .gt. bogus2) then
            iq     = iq + 1
            qt(iq) = term1
          end if
  570   continue
c
        print *, 'iq = ', iq
c
        if (iq .gt. jcrit) then
          call moment(qt, nmod, iq, avv1, adv1, sdd1, var1, sk, rkt)
c
          ts3(mm) = avv1
          ts4(mm) = sdd1
        else
          ts3(mm) = bogus1
          ts4(mm) = bogus1
        end if
        ts5(mm) = iq
  571 continue
c
c
c
c
c
c
c-------------------------------------------------------------------------------
c     ** 6.0 **
c
c     CALCULATION OF SIMPLE STATISTICAL PROPERTIES OF THE MULTI-MODEL
c     DISTRIBUTION OF TIME SERIES STATISTICS.
c
c     First loop over the number of time series statistics processed (typically
c     four). For each time series statistic, calculate the mean and standard 
c     deviation of the multi-model distribution of statistic values. 
c-------------------------------------------------------------------------------
c
      do 601 mm = 1, 4 ! Loop over number of time series statistics.
        do 600 ll = 1, nmod ! Loop over number of models.
          wk4(ll) = wk3(mm, ll)
  600   continue
c
        call moment(wk4, nmod, nmod, ave1, adev1, sdev1, var1,
     &  skew1, curt1)
c
c-------------------------------------------------------------------------------
c     ** 6.1 **
c
c     For each time series statistic, sort multi-model distribution and 
c     calculate median.
c-------------------------------------------------------------------------------
c
        call sort(nmod, wk4)
c
        n2 = nmod / 2
        if (2 * n2 .eq. nmod) then
          xmed1 = 0.5 * (wk4(n2) + wk4(n2 + 1))
        else
          xmed1 = wk4(n2 + 1)
        end if
c
        print *, 'Min  = ', wk4(1)
        print *, 'Max  = ', wk4(nmod)
        print *, 'Med  = ', xmed1
c
        wk6(1, mm) = wk4(1)     ! Distribution minimum.
        wk6(2, mm) = wk4(nmod)  ! Distribution maximum.
        wk6(3, mm) = ave1       ! Distribution average.
        wk6(4, mm) = xmed1      ! Distribution median. 
        wk6(5, mm) = sdev1      ! Distribution standard deviation.
  601 continue
c
c
c
c
c
c
c-------------------------------------------------------------------------------
c     ** 9.0 **
c
c     OUTPUT RESULTS.       
c-------------------------------------------------------------------------------
c
      do 9000 ll = 1, nmod ! Loop over number of models
c
        if (ll .eq.  1) nreal = nr1
        if (ll .eq.  2) nreal = nr2
        if (ll .eq.  3) nreal = nr3
        if (ll .eq.  4) nreal = nr4
        if (ll .eq.  5) nreal = nr5
        if (ll .eq.  6) nreal = nr6
        if (ll .eq.  7) nreal = nr7
        if (ll .eq.  8) nreal = nr8
        if (ll .eq.  9) nreal = nr9
        if (ll .eq. 10) nreal = nr10
        if (ll .eq. 11) nreal = nr11
        if (ll .eq. 12) nreal = nr12
        if (ll .eq. 13) nreal = nr13
        if (ll .eq. 14) nreal = nr14
        if (ll .eq. 15) nreal = nr15
        if (ll .eq. 16) nreal = nr16
        if (ll .eq. 17) nreal = nr17
        if (ll .eq. 18) nreal = nr18
        if (ll .eq. 19) nreal = nr19
        if (ll .eq. 20) nreal = nr20
c
c
        p1a  = p1(ll)
        p1b  = p1a(23:42)
        ida  = index(p1a, ' ') - 1
        idb  = index(p1b, ' ') - 2
        p9   = '_nofilt'
c
        do 8000 kk = 1, nreal ! Loop over number of realizations.
          imm2 = ileng1(kk, ll) ! Set number of months processed for each model.
          p8   = p7(kk, ll)     ! Set run number.
c
c-------------------------------------------------------------------------------
c     ** 9.1 **
c
c     Write full text strings for all output files (individual models).
c-------------------------------------------------------------------------------
c
          print *, 'DATA OUTPUT'
c
          cyr(1:22) = '_r'//cyr1//'_'//cyr2//'_s'//cyr3//'_'//cyr4
          s2a = p6a(1:id2a)//'-'//creg1(1:4)//'_'
          isa = index(s2a, ' ') - 1
c
          fo1 = s1//s2a(1:isa)//p1b(1:idb)//'_'//p8(9:12)//cyr//p9//'.d'
c
          print *, 'fo1 = ', fo1
c
c-------------------------------------------------------------------------------
c     ** 9.2 **
c
c     Open output files (one per model and realization).
c-------------------------------------------------------------------------------
c
          open(unit=1, file=fo1, form='formatted', status='new')
c
c-------------------------------------------------------------------------------
c     ** 9.3 **
c
c     Write header information regarding domain selection.
c-------------------------------------------------------------------------------
c
          write(1, 900)
          if (creg1(1:4) .eq. 'TO15') write(1, 7011) 
          if (creg1(1:4) .eq. 'TO20') write(1, 7012)
          if (creg1(1:4) .eq. 'TO30') write(1, 7013)
          if (creg1(1:4) .eq. 'TL15') write(1, 7014)
          if (creg1(1:4) .eq. 'TL20') write(1, 7015)
          if (creg1(1:4) .eq. 'TL30') write(1, 7016)
          if (creg1(1:4) .eq. 'TLO1') write(1, 7017)
          if (creg1(1:4) .eq. 'TLO2') write(1, 7018)
          if (creg1(1:4) .eq. 'TLO3') write(1, 7019)
          if (creg1(1:4) .eq. 'NI34') write(1, 7020)
          if (creg1(1:4) .eq. 'OC50') write(1, 7021)
          if (creg1(1:4) .eq. 'SHEM') write(1, 7022)
          if (creg1(1:4) .eq. 'NHEM') write(1, 7023)
          if (creg1(1:4) .eq. 'GLBL') write(1, 7024)
          if (creg1(1:4) .eq. 'GLAN') write(1, 7025)
          if (creg1(1:4) .eq. 'GOCN') write(1, 7026)
          if (creg1(1:4) .eq. 'SHEX') write(1, 7027)
          if (creg1(1:4) .eq. 'NHEX') write(1, 7028)
          if (creg1(1:4) .eq. 'NHM1') write(1, 7029)
          if (creg1(1:4) .eq. 'NHM2') write(1, 7030)
          if (creg1(1:4) .eq. 'NHM3') write(1, 7031)
          if (creg1(1:4) .eq. 'NHM4') write(1, 7032)
          if (creg1(1:4) .eq. 'NHM5') write(1, 7033)
          if (creg1(1:4) .eq. 'SHM1') write(1, 7034)
          if (creg1(1:4) .eq. 'SHM2') write(1, 7035)
          if (creg1(1:4) .eq. 'SHM3') write(1, 7036)
          if (creg1(1:4) .eq. 'SHM4') write(1, 7037)
          if (creg1(1:4) .eq. 'SHM5') write(1, 7038)
          if (creg1(1:4) .eq. 'NHM6') write(1, 7039)
          if (creg1(1:4) .eq. 'SHM6') write(1, 7040)
          if (creg1(1:4) .eq. 'SHM7') write(1, 7041)
          if (creg1(1:4) .eq. 'GLB2') write(1, 7042)
          if (creg1(1:4) .eq. 'GLB3') write(1, 7043)
c
c-------------------------------------------------------------------------------
c     ** 9.4 **
c
c     Write header information regarding selected variable.
c-------------------------------------------------------------------------------
c
          if (p6a(1:3) .eq. 'tls') then
            write(1, 9031) p6a(1:id2a)
          else if (p6a(1:3) .eq. 'tmt') then
            write(1, 9032) p6a(1:id2a)
          else if (p6a(1:3) .eq. 'tlt') then
            write(1, 9033) p6a(1:id2a)
          end if
c
          write(1, 9034)
          write(1, 9035)
          write(1, 9036)
          write(1, 9037)
c
c-------------------------------------------------------------------------------
c     ** 9.5 **
c
c     Write header information regarding other user selections. 
c-------------------------------------------------------------------------------
c
          write(1, 904) p6a(1:id2a)
          write(1, 905)
          write(1, 906)
          write(1, 907) x1(idd), x2(idd), y1(idd), y2(idd)
          write(1, 908)
          write(1, 909) iyr1, iyr2                
          write(1, 910) iyr3, iyr4                 
          write(1, 911)
          write(1, 912)
          write(1, 913) ileng1(kk, ll)
c
c-------------------------------------------------------------------------------
c     ** 9.6 **
c
c     Write information from basic quality control.      
c-------------------------------------------------------------------------------
c
          write(1, 920)
          write(1, 921)
          write(1, 922)
          write(1, 923)
          write(1, 924)
          write(1, 925)
          write(1, 926)
c
          write(1, 9260)
          write(1, 9261) p8(12:12)
c
          write(1, 9270) miss1(kk, ll)
          write(1, 9271) qnorm2(kk, ll)
          write(1, 9272) iqcrit(kk, ll)
          write(1, 9273) pqcrit(kk, ll)
c
c-------------------------------------------------------------------------------
c     ** 9.7 **
c
c     Write full path to input model data.
c-------------------------------------------------------------------------------
c
          write(1, 9100) f2(kk, ll)
c
c-------------------------------------------------------------------------------
c     ** 9.8 **
c
c     Write column headers for output model data (individual models).
c-------------------------------------------------------------------------------
c
          write(1, 941) p8(12:12)
c
c-------------------------------------------------------------------------------
c     ** 9.9 **
c
c     Write output data for individual models. The data are output for the
c     entire length of the ex20c3m experiment - not for the analysis period 
c     only.
c-------------------------------------------------------------------------------
c
          do 990 mm = 1, imm2 ! Loop over number of months 
            t1 = (rtime2(mm, kk, ll) / 12.0) + 1800.0
c           print *, 'mm, t1 = ', mm, t1
c
            write(1, 951) mm, t1, ts1(mm, kk, ll), 
     &      ngptact(mm, kk, ll), fcover(mm, kk, ll),
     &      qnorm1(mm, kk, ll)
  990     continue
          close(unit=1)
c
c-------------------------------------------------------------------------------
c     ** 9.10 **
c
c     End loop over number of realizations and number of individual models.
c-------------------------------------------------------------------------------
c
 8000   continue ! End loop over number of realizations.
 9000 continue ! End loop over number of models.
c
c-------------------------------------------------------------------------------
c     ** 9.11 **
c
c     Generate file name for file containing time series statistics for 
c     individual model ex20c3m realizations.
c-------------------------------------------------------------------------------
c
      p1b = 'STATISTICS '
      idb = index(p1b, ' ') - 1
c
      print *, 'DATA OUTPUT'
c
      fo1 = s1//s2a(1:isa)//p1b(1:idb)//cyr//p9//'.d'
c
      print *, 'fo1 = ', fo1
c
c-------------------------------------------------------------------------------
c     ** 9.12 **
c
c     Open file for writing time series statistics for individual model ex20c3m 
c     realizations.
c-------------------------------------------------------------------------------
c
      open(unit=1, file=fo1, form='formatted', status='new')
c
c-------------------------------------------------------------------------------
c     ** 9.13 **
c
c     Write header information for file containing time series statistics for 
c     individual model ex20c3m realizations.
c-------------------------------------------------------------------------------
c
      write(1, 900)
      if (creg1(1:4) .eq. 'TO15') write(1, 7011)
      if (creg1(1:4) .eq. 'TO20') write(1, 7012)
      if (creg1(1:4) .eq. 'TO30') write(1, 7013)
      if (creg1(1:4) .eq. 'TL15') write(1, 7014)
      if (creg1(1:4) .eq. 'TL20') write(1, 7015)
      if (creg1(1:4) .eq. 'TL30') write(1, 7016)
      if (creg1(1:4) .eq. 'TLO1') write(1, 7017)
      if (creg1(1:4) .eq. 'TLO2') write(1, 7018)
      if (creg1(1:4) .eq. 'TLO3') write(1, 7019)
      if (creg1(1:4) .eq. 'NI34') write(1, 7020)
      if (creg1(1:4) .eq. 'OC50') write(1, 7021)
      if (creg1(1:4) .eq. 'SHEM') write(1, 7022)
      if (creg1(1:4) .eq. 'NHEM') write(1, 7023)
      if (creg1(1:4) .eq. 'GLBL') write(1, 7024)
      if (creg1(1:4) .eq. 'GLAN') write(1, 7025)
      if (creg1(1:4) .eq. 'GOCN') write(1, 7026)
      if (creg1(1:4) .eq. 'SHEX') write(1, 7027)
      if (creg1(1:4) .eq. 'NHEX') write(1, 7028)
      if (creg1(1:4) .eq. 'NHM1') write(1, 7029)
      if (creg1(1:4) .eq. 'NHM2') write(1, 7030)
      if (creg1(1:4) .eq. 'NHM3') write(1, 7031)
      if (creg1(1:4) .eq. 'NHM4') write(1, 7032)
      if (creg1(1:4) .eq. 'NHM5') write(1, 7033)
      if (creg1(1:4) .eq. 'SHM1') write(1, 7034)
      if (creg1(1:4) .eq. 'SHM2') write(1, 7035)
      if (creg1(1:4) .eq. 'SHM3') write(1, 7036)
      if (creg1(1:4) .eq. 'SHM4') write(1, 7037)
      if (creg1(1:4) .eq. 'SHM5') write(1, 7038)
      if (creg1(1:4) .eq. 'NHM6') write(1, 7039)
      if (creg1(1:4) .eq. 'SHM6') write(1, 7040)
      if (creg1(1:4) .eq. 'SHM7') write(1, 7041)
      if (creg1(1:4) .eq. 'GLB2') write(1, 7042)
      if (creg1(1:4) .eq. 'GLB3') write(1, 7043)
c
c-------------------------------------------------------------------------------
c     ** 9.14 **
c
c     Write header information regarding selected variable.
c-------------------------------------------------------------------------------
c
        if (p6a(1:3) .eq. 'tls') then
          write(1, 9031) p6a(1:id2a)
        else if (p6a(1:3) .eq. 'tmt') then
          write(1, 9032) p6a(1:id2a)
        else if (p6a(1:3) .eq. 'tlt') then
          write(1, 9033) p6a(1:id2a)
        end if
c
        write(1, 9034)
        write(1, 9035)
        write(1, 9036)
        write(1, 9037)
c
c-------------------------------------------------------------------------------
c     ** 9.15 **
c
c     Write header information regarding other user selections. 
c-------------------------------------------------------------------------------
c
        write(1, 904) p6a(1:id2a)
        write(1, 905)
        write(1, 906)
        write(1, 907) x1(idd), x2(idd), y1(idd), y2(idd)
        write(1, 908)
        write(1, 909) iyr1, iyr2
        write(1, 910) iyr3, iyr4
        write(1, 911)
        write(1, 912)
c
c-------------------------------------------------------------------------------
c     ** 9.16 **
c
c     Write header information regarding output statistics.     
c-------------------------------------------------------------------------------
c
      write(1, 8410)
      write(1, 8411)
c     write(1, 8411) p6a(1:id2a)
      write(1, 8412)
      write(1, 8413)
      write(1, 8414)
      write(1, 8415)
      write(1, 8416)
      write(1, 8417)
      write(1, 8418)
      write(1, 8419)
      write(1, 8420)
      write(1, 8421)
      write(1, 8422)
c
      write(1, 8430)
c
c-------------------------------------------------------------------------------
c     ** 9.17 **
c
c     Write time series statistics for individual model ex20c3m realizations.
c-------------------------------------------------------------------------------
c
      ipp = 0
      do 9171 ll = 1, nmod ! Loop over the number of models.
        if (ll .eq.  1) nreal = nr1
        if (ll .eq.  2) nreal = nr2
        if (ll .eq.  3) nreal = nr3
        if (ll .eq.  4) nreal = nr4
        if (ll .eq.  5) nreal = nr5
        if (ll .eq.  6) nreal = nr6
        if (ll .eq.  7) nreal = nr7
        if (ll .eq.  8) nreal = nr8
        if (ll .eq.  9) nreal = nr9
        if (ll .eq. 10) nreal = nr10
        if (ll .eq. 11) nreal = nr11
        if (ll .eq. 12) nreal = nr12
        if (ll .eq. 13) nreal = nr13
        if (ll .eq. 14) nreal = nr14
        if (ll .eq. 15) nreal = nr15
        if (ll .eq. 16) nreal = nr16
        if (ll .eq. 17) nreal = nr17
        if (ll .eq. 18) nreal = nr18
        if (ll .eq. 19) nreal = nr19
        if (ll .eq. 20) nreal = nr20
c
        p1a = p1(ll)
        p1b = p1a(23:42)
        ida = index(p1a, ' ') - 1
        idb = index(p1b, ' ') - 2
c
        do 9170 kk = 1, nreal
          p8  = p7(kk, ll)
          ipp = ipp + 1
          write(1, 8431) ipp, p1b(1:idb), p8(12:12),
     &    miss1(kk, ll),    sdy(kk, ll),    r1(kk, ll),
     &    trend(kk, ll),    tlc(kk, ll), serr2(kk, ll),
     &      ptt(kk, ll),  rmiss(kk, ll), rneff(kk, ll),
     &   qnorm2(kk, ll), iqcrit(kk, ll), pqcrit(kk, ll)
 9170   continue
 9171 continue
      write(1, 805)
c
c-------------------------------------------------------------------------------
c     ** 9.18 **
c
c     Write ensemble-mean values of time series statistics.
c-------------------------------------------------------------------------------
c
      write(1, 8510)
      write(1, 8511)
      write(1, 8512)
      write(1, 8513)
c
      do 9200 ll = 1, nmod
        p1a = p1(ll)
        p1b = p1a(23:42)
        ida = index(p1a, ' ') - 1
        idb = index(p1b, ' ') - 2
c
        write(1, 8514) ll, p1b(1:idb),
     &  (wk3(mm, ll), mm = 1, 4)
 9200 continue
c
      write(1, 8515) (wk6(1, mm), mm = 1, 4)
      write(1, 8516) (wk6(2, mm), mm = 1, 4)
      write(1, 8517) (wk6(3, mm), mm = 1, 4)
      write(1, 8518) (wk6(4, mm), mm = 1, 4)
      write(1, 8519) (wk6(5, mm), mm = 1, 4)
      write(1, 805)
c
      close(unit=1)
c
c-------------------------------------------------------------------------------
c     ** 9.19 **
c
c     Write multi-model average time series. First generate file name for 
c     output file.
c-------------------------------------------------------------------------------
c
      p1b = 'MEAN_MODEL '
      idb = index(p1b, ' ') - 1
c
      print *, 'DATA OUTPUT'
c
      fo2 = s1//s2a(1:isa)//p1b(1:idb)//cyr//p9//'.d'
c
      print *, 'fo2 = ', fo2
c
c-------------------------------------------------------------------------------
c     ** 9.20 **
c
c     Open file for writing multi-model average time series.
c-------------------------------------------------------------------------------
c
      open(unit=1, file=fo2, form='formatted', status='new')
c
c-------------------------------------------------------------------------------
c     ** 9.21 **
c
c     Write header information for file containing multi-model average time 
c     series.
c-------------------------------------------------------------------------------
c
      write(1, 900)
      if (creg1(1:4) .eq. 'TO15') write(1, 7011)
      if (creg1(1:4) .eq. 'TO20') write(1, 7012)
      if (creg1(1:4) .eq. 'TO30') write(1, 7013) 
      if (creg1(1:4) .eq. 'TL15') write(1, 7014)
      if (creg1(1:4) .eq. 'TL20') write(1, 7015)
      if (creg1(1:4) .eq. 'TL30') write(1, 7016)
      if (creg1(1:4) .eq. 'TLO1') write(1, 7017)
      if (creg1(1:4) .eq. 'TLO2') write(1, 7018)
      if (creg1(1:4) .eq. 'TLO3') write(1, 7019)
      if (creg1(1:4) .eq. 'NI34') write(1, 7020)
      if (creg1(1:4) .eq. 'OC50') write(1, 7021)
      if (creg1(1:4) .eq. 'SHEM') write(1, 7022)
      if (creg1(1:4) .eq. 'NHEM') write(1, 7023)
      if (creg1(1:4) .eq. 'GLBL') write(1, 7024)
      if (creg1(1:4) .eq. 'GLAN') write(1, 7025)
      if (creg1(1:4) .eq. 'GOCN') write(1, 7026)
      if (creg1(1:4) .eq. 'SHEX') write(1, 7027)
      if (creg1(1:4) .eq. 'NHEX') write(1, 7028)
      if (creg1(1:4) .eq. 'NHM1') write(1, 7029)
      if (creg1(1:4) .eq. 'NHM2') write(1, 7030)
      if (creg1(1:4) .eq. 'NHM3') write(1, 7031)
      if (creg1(1:4) .eq. 'NHM4') write(1, 7032)
      if (creg1(1:4) .eq. 'NHM5') write(1, 7033)
      if (creg1(1:4) .eq. 'SHM1') write(1, 7034)
      if (creg1(1:4) .eq. 'SHM2') write(1, 7035)
      if (creg1(1:4) .eq. 'SHM3') write(1, 7036)
      if (creg1(1:4) .eq. 'SHM4') write(1, 7037)
      if (creg1(1:4) .eq. 'SHM5') write(1, 7038)
      if (creg1(1:4) .eq. 'NHM6') write(1, 7039)
      if (creg1(1:4) .eq. 'SHM6') write(1, 7040)
      if (creg1(1:4) .eq. 'SHM7') write(1, 7041)
      if (creg1(1:4) .eq. 'GLB2') write(1, 7042)
      if (creg1(1:4) .eq. 'GLB3') write(1, 7043)
c
c-------------------------------------------------------------------------------
c     ** 9.22 **
c
c     Write header information regarding selected variable.
c-------------------------------------------------------------------------------
c
      if (p6a(1:3) .eq. 'tls') then
        write(1, 9031) p6a(1:id2a)
      else if (p6a(1:3) .eq. 'tmt') then
        write(1, 9032) p6a(1:id2a)
      else if (p6a(1:3) .eq. 'tlt') then
        write(1, 9033) p6a(1:id2a)
      end if
c
      write(1, 9034)
      write(1, 9035)
      write(1, 9036)
      write(1, 9037)
c
c-------------------------------------------------------------------------------
c     ** 9.23 **
c
c     Write header information regarding other user selections. 
c-------------------------------------------------------------------------------
c
      write(1, 904) p6a(1:id2a)
      write(1, 905)
      write(1, 906)
      write(1, 907) x1(idd), x2(idd), y1(idd), y2(idd)
      write(1, 908)
      write(1, 909) iyr1, iyr2
      write(1, 910) iyr3, iyr4
      write(1, 911)
      write(1, 912)
      write(1, 913) kcc            
c
c-------------------------------------------------------------------------------
c     ** 9.24 **
c
c     Write multi-model average time series.
c-------------------------------------------------------------------------------
c
      write(1, 960)
      do 6000 mm = 1, kcc ! Loop over number of months
        t1 = (rtime3(mm) / 12.0) + 1800.0
        write(1, 961) mm, t1, ts3(mm), ts4(mm), ts5(mm)
 6000 continue
      close(unit=1)
c
c
c
c
c
c
c-------------------------------------------------------------------------------
c     ** 10.0 **
c
c     FORMAT STATEMENTS.
c-------------------------------------------------------------------------------
c
  803 format(i4, a16, 4i10, f10.4, i10)
  804 format(i4, a16, 5f10.4, i10)
  805 format(/)
c
c-------------------------------------------------------------------------------
c     ** 10.1 **
c
c     Format statements for output files containing model time series. Program
c     name and "dictionary" information regarding selected domain.
c-------------------------------------------------------------------------------
c
  900 format(t2, 'RESULTS FROM PROGRAM ',
     &'SPATIAL_AVERAGE_EX20C3M_EQMSU_METHOD1a_ALLMODELS1.F')       
c
 7011 format(t2, 'Domain: Tropical Ocean ',
     &'Region 1 (TO15; 180W-180E; 15N-15S)')
 7012 format(t2, 'Domain: Tropical Ocean ',
     &'Region 2 (TO20; 180W-180E; 20N-20S)')
 7013 format(t2, 'Domain: Tropical Ocean ',
     &'Region 3 (TO30; 180W-180E; 30N-30S)')
c
 7014 format(t2, 'Domain: Tropical Land ',
     &'Region 1 (TL15; 180W-180E; 15N-15S)')
 7015 format(t2, 'Domain: Tropical Land ',
     &'Region 2 (TL20; 180W-180E; 20N-20S)')
 7016 format(t2, 'Domain: Tropical Land ',
     &'Region 3 (TL30; 180W-180E; 30N-30S)')
c
 7017 format(t2, 'Domain: Tropical Land ',
     &'and Ocean Region 1 (TLO1; 180W-180E; 15N-15S)')
 7018 format(t2, 'Domain: Tropical Land ',
     &'and Ocean Region 2 (TLO2; 180W-180E; 20N-20S)')
 7019 format(t2, 'Domain: Tropical Land ',
     &'and Ocean Region 3 (TLO3; 180W-180E; 30N-30S)')
c
 7020 format(t2, 'Domain: Nino 3.4 ',
     &'Region (NI34; 170W-120W; 5N-5S)')
 7021 format(t2, 'Domain: Oceans equatorward ',
     &'of 50 degrees latitude (OC50; 180W-180E; 50N-50S)')
 7022 format(t2, 'Domain: Southern Hemisphere ',
     &'land and ocean (SHEM; 180W-180E; 90S-0)')
 7023 format(t2, 'Domain: Northern Hemisphere ',
     &'land and ocean (NHEM; 180W-180E; 0-90N)')
 7024 format(t2, 'Domain: Global ',
     &'land and ocean (GLBL; 180W-180E; 90S-90N)')
 7025 format(t2, 'Domain: Global ',
     &'land (GLAN; 180W-180E; 90S-90N)')
 7026 format(t2, 'Domain: Global ',
     &'ocean (GOCN; 180W-180E; 90S-90N)')
 7027 format(t2, 'Domain: Southern Hemisphere extratropics ',
     &'land and ocean (SHEX; 180W-180E; 90S-20S)')
 7028 format(t2, 'Domain: Northern Hemisphere extratropics ',
     &'land and ocean (NHEX; 180W-180E; 20N-90N)')
c
 7029 format(t2, 'Domain: Northern Hemisphere ',
     &'Land and ocean Region 1 (NHM1; 180W-180E; 0-10N)')
 7030 format(t2, 'Domain: Northern Hemisphere ',
     &'Land and ocean Region 2 (NHM2; 180W-180E; 10N-30N)')
 7031 format(t2, 'Domain: Northern Hemisphere ',
     &'Land and ocean Region 3 (NHM3; 180W-180E; 30N-50N)')
 7032 format(t2, 'Domain: Northern Hemisphere ',
     &'Land and ocean Region 4 (NHM4; 180W-180E; 50N-70N)')
 7033 format(t2, 'Domain: Northern Hemisphere ',
     &'Land and ocean Region 5 (NHM5; 180W-180E; 70N-82.5N)')
c
 7034 format(t2, 'Domain: Southern Hemisphere ',
     &'Land and ocean Region 1 (SHM1; 180W-180E; 0-10S)')
 7035 format(t2, 'Domain: Southern Hemisphere ',
     &'Land and ocean Region 2 (SHM2; 180W-180E; 10S-30S)')
 7036 format(t2, 'Domain: Southern Hemisphere ',
     &'Land and ocean Region 3 (SHM3; 180W-180E; 30S-50S)')
 7037 format(t2, 'Domain: Southern Hemisphere ',
     &'Land and ocean Region 4 (SHM4; 180W-180E; 50S-70S)')
 7038 format(t2, 'Domain: Southern Hemisphere ',
     &'Land and ocean Region 5 (SHM5; 180W-180E; 70S-82.5S)')
c
 7039 format(t2, 'Domain: Northern Hemisphere ',
     &'Land and ocean Region 6 (NHM6; 180W-180E; 0-82.5N)')
 7040 format(t2, 'Domain: Southern Hemisphere ',
     &'Land and ocean Region 6 (SHM6; 180W-180E; 0-70S)')
 7041 format(t2, 'Domain: Southern Hemisphere ',
     &'Land and ocean Region 7 (SHM7; 180W-180E; 0-82.5S)')
 7042 format(t2, 'Domain: Near-global ',
     &'land and ocean (GLB2; 180W-180E; 82.5S-82.5N)')
 7043 format(t2, 'Domain: Near-global ',
     &'land and ocean (GLB3; 180W-180E; 70S-82.5N)')
c
c-------------------------------------------------------------------------------
c     ** 10.2 **
c
c     Header information regarding selected variable.
c-------------------------------------------------------------------------------
c
 9031 format(t2, 'Selected variable: ', a3, ' (synthetic ',
     &'MSU lower stratospheric temperature)')
 9032 format(t2, 'Selected variable: ', a3, ' (synthetic ',
     &'MSU mid- to upper-tropospheric temperature)')
 9033 format(t2, 'Selected variable: ', a3, ' (synthetic ',
     &'MSU lower tropospheric temperature)')
 9034 format(t2, 'Synthetic MSU temperature calculated using ',
     &'Method 1a developed by Dr. Carl Mears (RSS)')
 9035 format(t2, 'Method 1a uses model surface temperature (ts), ',
     &'atmospheric temperature (ta), surface pressure (ps), and ',
     &'land fraction (sftlf)')
 9036 format(t2, 'Model sea-ice coverage data are not used in ',
     &'Method 1a')
 9037 format(t2, 'Synthetic MSU temperatures calculated with ',
     &'Method 1a are designated as "both" in the filename')
c
c-------------------------------------------------------------------------------
c     ** 10.3 **
c
c     Header information regarding other user selections.
c-------------------------------------------------------------------------------
c
  904 format(t2, 'Program calculates spatially-averaged ', a3, 
     &' changes for selected domain')
  905 format(t2, 'Input data are monthly means from IPCC ',
     &'EX20C3M century simulations')
  906 format(t2, 'Analysis performed on native model grid')
  907 format(t2, 'Selected longitude/latitude ',
     &'boundaries of domain (x1, x2, y1, y2): ', 4f7.1)
  908 format(t2, 'Output data: Spatial mean of ',
     &'domain')
  909 format(t2, 'Reference period: ', i4, '-', i4)
  910 format(t2, 'Analysis period:  ', i4, '-', i4)
c9090 format(t2, 'Reference period is entire model control run')
c9091 format(t2, 'Analysis period is entire model control run')
  911 format(t2, 'No filtering of output time series')
  912 format(t2, 'No detrending selected')
  913 format(t2, 'Total number of months processed: ', i5 /)       
c
c-------------------------------------------------------------------------------
c     ** 10.4 **
c
c     Quality control information.
c-------------------------------------------------------------------------------
c
  920 format(t2, 'Basic quality control information')
  921 format(t2, 'NGPT:   Actual number of model grid-points ',
     &'(with valid data) in selected domain')
  922 format(t2, 'COVER:  Fractional data coverage ',
     &'in selected domain')
  923 format(t2, 'MISS1:  Total number of months for which a ',
     &'valid spatial average could not be computed')
  924 format(t2, 'QNORM1: Absolute value of ',
     &'normalized first difference')
  925 format(t2, 'QNORM2: Time-average of absolute value of ',
     &'normalized first differences')
  926 format(t2, 'IQCRIT: Number of normalized first differences ',
     &'(FD) > 5 std. deviations')
 9260 format(t2, 'PQCRIT: 100 x [IQCRIT / (total number of first ',
     &'differences in time series)]'/)
c
 9261 format(t2, 'TYPE', 9x, '   RUN', a1)
c
 9270 format(t2, 'MISS1 ', 4x, i10)
 9271 format(t2, 'QNORM2', 4x, f10.5)
 9272 format(t2, 'IQCRIT', 4x, i10)
 9273 format(t2, 'PQCRIT', 4x, f10.5/)
c
 9280 format(t2, 'MISS1 ', 4x, 2i10)
 9281 format(t2, 'QNORM2', 4x, 2f10.5)
 9282 format(t2, 'IQCRIT', 4x, 2i10)
 9283 format(t2, 'PQCRIT', 4x, 2f10.5/)
c
 9290 format(t2, 'MISS1 ', 4x, 3i10)
 9291 format(t2, 'QNORM2', 4x, 3f10.5)
 9292 format(t2, 'IQCRIT', 4x, 3i10)
 9293 format(t2, 'PQCRIT', 4x, 3f10.5/)
c
c-------------------------------------------------------------------------------
c     ** 10.5 **
c
c     Format statements for output files containing model time series. Path to
c     input climate model data.
c-------------------------------------------------------------------------------
c
 9100 format(t2, 'Input dataset #1: ', a120)                     
c
c-------------------------------------------------------------------------------
c     ** 10.6 **
c
c     Format statements for output files containing model time series. Column
c     headers for output data.
c-------------------------------------------------------------------------------
c
  941 format(1x, 'No.', 6x, 'Time', 6x, 'RUN', a1, 9x, 'NGPT', 3x,
     &'COVER', 5x, 'QNORM1')
c
c-------------------------------------------------------------------------------
c     ** 10.7 **
c
c     Format statements for output files containing model time series. Data   
c     output format.
c-------------------------------------------------------------------------------
c
  951 format(i5, f10.3, 1x, f10.5, i10, 2f10.5)
c
c-------------------------------------------------------------------------------
c     ** 10.8 **
c
c     Format statements for output file containing time series statistics for
c     individual model ex20c3m realizations. Brief explanation of output 
c     statistics.
c-------------------------------------------------------------------------------
c
 8410 format(t2, 'BASIC STATISTICS FOR MODEL TIME SERIES')       
c8411 format(t2, 'MEAN:   Spatial average of ', a3)
 8411 format(t2, 'MISS1:  Number of months in entire time series ',
     &'for which a valid spatial average could not be computed')
 8412 format(t2, 'SDEV:   Temporal standard deviation of ',
     &'spatial average (degrees C)')
 8413 format(t2, 'AR-1:   Lag-1 temporal autocorrelation of ',
     &'spatial average')
 8414 format(t2, 'TREND:  Least-squares linear trend (degrees C ',
     &'per decade)')        
 8415 format(t2, 'TOTLC:  Total linear change over analysis period ',
     &'(degrees C)')
 8416 format(t2, 'SERR2:  Standard error of linear trend ',
     &'(adjusted for temporal autocorrelation; degrees C ',
     &'per decade)')
 8417 format(t2, 'TPROB3: p-value for test of trend significance',
     &' (adjusted for temporal autocorrelation in DOF and SE)')
 8418 format(t2, 'NMISS:  Number of non-missing monthly values ',
     &'during analysis period')
 8419 format(t2, 'RNEFF:  Effective sample size', 
     &' (adjusted for temporal autocorrelation)')
 8420 format(t2, 'QNORM2: Time-average of normalized first ',
     &'differences')
 8421 format(t2, 'IQCRIT: Number of normalized first differences ',
     &'(FD) > 5 std. deviations')
 8422 format(t2, 'PQCRIT: 100 x [IQCRIT / (total number of first ',
     &'differences in time series)]'/)
c
 8430 format('NO.', 6x, 'MODEL', 8x, 'R#', 3x, 'MISS1', 6x, 
     &'SDEV', 6x, 'AR-1', 6x, 'TREND', 5x, 'TOTLC', 5x, 'SERR2', 5x,
     &'TPROB3', 5x, 'NMISS', 5x, 'RNEFF', 4x, 'QNORM2', 4x, 
     &'IQCRIT', 3x, 'PQCRIT')
 8431 format(i3, 1x, a17, 2x, a1, i8, 2x, 6f10.4, 2f10.1, f10.4, 
     &i9, f10.4)
c
c-------------------------------------------------------------------------------
c     ** 10.9 **
c
c     Format statements for output file containing "ensemble mean" time series 
c     statistics for model ex20c3m realizations. Brief explanation of output 
c     statistics.
c-------------------------------------------------------------------------------
c
 8510 format(t2, 'BASIC STATISTICS FOR MODEL ',
     &'TIME SERIES: ENSEMBLE-MEAN VALUES')
 8511 format(t2, 'Results are for realization #1 if only one ',
     &'ensemble member is available')     
 8512 format(t2, 'Results are ensemble-means of time series ',
     &'statistics, not statistics of ensemble-mean time series')
 8513 format('NO.', 10x, 'MODEL', 9x, 'MEAN', 6x,
     &'SDEV', 6x, 'AR-1', 5x, 'TREND')                               
 8514 format(i3, 1x, a17, 4f10.4)
 8515 format(t5, 'Minimum', t22, 4f10.4)
 8516 format(t5, 'Maximum', t22, 4f10.4)
 8517 format(t5, 'Mean', t22, 4f10.4)
 8518 format(t5, 'Median', t22, 4f10.4)
 8519 format(t5, 'Std. Dev.', t22, 4f10.4)
c
c-------------------------------------------------------------------------------
c     ** 10.10 **
c
c     Format statements for output file containing multi-model average time   
c     series.
c-------------------------------------------------------------------------------
c
  960 format(1x, 'No.', 5x, 'Time', 7x, 'MEAN', 6x, 'SDEV', 5x,
     &'NMODELS')
  961 format(i4, f10.3, 2x, 3f10.4)
c
c
      stop
      end
c-------------------------------------------------------------------------------
c-------------------------------------------------------------------------------
c     SUBROUTINE RD2D_XY
c
c     Date:    21-MAY-2001.
c     Author:  Ben Santer, LLNL, Livermore.
c     Version: # 1.0 / Single precision.
c
c
c
c     This subroutine reads an input two-dimensional longitude-latitude
c     field.
c
c
c
      subroutine rd2d_xy(nx, ny, nxi, nyi, mm, bx1, bx2, by1, by2,
     &infil, rlon, rlat, rtime, amask, field1, varnam, yweights)
c
c
c
      implicit real (o-z, a-h)
c
c     ** INPUT **
c     integer nx                     Input spatial dimension (longitude)
c     integer ny                     Input spatial dimension (latitude)
c     integer mm                     Current time (months)
c     real bx1, bx2                  Longitude bounds of current region
c     real by1, by2                  Latitude bounds of current region
c     character*(*) infil            Input file name
c     character*16 varnam            Variable name
c     character* 8 yweights          Weight type
c
c     ** OUTPUT **
c     integer nxi                    Output spatial dimension (longitude)
c     integer nyi                    Output spatial dimension (latitude)
c     dimension field1(nx, ny)       Output data
c     dimension amask(nx, ny)        Area weights for data
c     dimension rlon(nx)             Longitude coordinates
c     dimension rlat(ny)             Latitude coordinates
c     dimension rtime(1)             Time coordinates
c
c
c
      integer nx, ny, nxi, nyi, mm
      integer i1, i2, i3, i4
      real bx1, bx2, by1, by2
      dimension field1(nx, ny)
      dimension amask(nx, ny)
      dimension rlon(nx), rlat(ny), rtime(1)
      character infil*(*), varnam*16, yweights*8
c
c
c
c-------------------------------------------------------------------------------
c     ** 1.0 **
c
c     Define domain and method of computing weights.
c-------------------------------------------------------------------------------
c
c     print *, 'Weighting = ', yweights
c
      call defdim(1, 1, 'longitude', 'width', 'range', bx1,
     &bx2, 360.0)
c
      call defdim(1, 2, 'latitude', yweights, 'range', by1,
     &by2, 0.0)
c
      call defdimi(1, 3, 'time', 'unit', mm, mm)
c
c-------------------------------------------------------------------------------
c     ** 1.1 **
c
c     Define expected dimensions of data.
c-------------------------------------------------------------------------------
c
      i1 = 0
      i2 = 0
      i3 = 1
      i4 = 0
c
c-------------------------------------------------------------------------------
c     ** 1.2 **
c
c     Define variable 1.
c-------------------------------------------------------------------------------
c
      call defvar(1, varnam, infil)
c
c-------------------------------------------------------------------------------
c     ** 1.3 **
c
c       Extract data.
c-------------------------------------------------------------------------------
c
      call getdata(1, nx, ny, 1, i4, i1, i2, i3, i4, amask,
     &field1)
      call getcoord(1, 1, rlon)
      call getcoord(1, 2, rlat)
      call getcoord(1, 3, rtime)
      nxi = i1
      nyi = i2
c
      return
      end
c-------------------------------------------------------------------------------
c-------------------------------------------------------------------------------
c     SUBROUTINE RDMASK.
c
c     Date:    05-DEC-2005.
c     Author:  Ben Santer and Karl Taylor, LLNL, Livermore.
c     Version: # 1.0 / Single precision.
c
c
c
c     This subroutine reads the fractional land coverage mask.
c
c
c
      subroutine rdmask(nx, ny, nxi, nyi, bx1, bx2, by1, by2,
     &infil, amask, fmask, yweights)
c
c
c
      implicit real (o-z, a-h)
c
c     ** INPUT **
c     integer nx                  Spatial dimension (longitude)
c     integer ny                  Spatial dimension (latitude)
c     real bx1, bx2               Longitude bounds of current region
c     real by1, by2               Latitude bounds of current region
c     character*(*)infil          Input file name
c     character* 8 yweights       Weight type
c
c     ** OUTPUT **
c     integer nxi                 Output spatial dimension (longitude)
c     integer nyi                 Output spatial dimension (latitude)
c     dimension amask(nx, ny)     Area weights
c     dimension fmask(nx, ny)     Fractional data coverage mask
c
c
c
      integer nx, ny, nxi, nyi
      real bx1, bx2, by1, by2
      dimension fmask(nx, ny), amask(nx, ny)
      character infil*(*), yweights*8
c
c
c
c-------------------------------------------------------------------------------
c     ** 1.0 **
c
c     Initialize EZGET.
c-------------------------------------------------------------------------------
c
      call initget
c
c-------------------------------------------------------------------------------
c     ** 1.1 **
c
c     Define domain and method of computing weights.
c-------------------------------------------------------------------------------
c
c     print *, 'Weighting = ', yweights
c
      call defdim(1, 1, 'longitude', 'width', 'range', bx1,
     &bx2, 360.0)
c
      call defdim(1, 2, 'latitude', yweights, 'range', by1,
     &by2, 0.0)
c
c-------------------------------------------------------------------------------
c     ** 1.2 **
c
c     Define expected dimensions of data.
c-------------------------------------------------------------------------------
c
      i1 = 0
      i2 = 0
      i3 = 0
      i4 = 0
      n3 = 0
      n4 = 0
c
c-------------------------------------------------------------------------------
c     ** 1.3 **
c
c     Define variable 2.
c-------------------------------------------------------------------------------
c
      call defvar(1, 'sftlf', infil)
c
c-------------------------------------------------------------------------------
c     ** 1.4 **
c
c       Extract data.
c-------------------------------------------------------------------------------
c
      call getdata(1, nx, ny, n3, n4, i1, i2, i3, i4, amask,
     &fmask)
      nxi = i1
      nyi = i2
c
c-------------------------------------------------------------------------------
c     ** 1.5 **
c
c     Release memory.
c-------------------------------------------------------------------------------
c
      call closeget
c
      return
      end
