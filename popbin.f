***
      PROGRAM popbin
***
*
* Evolves a population of binaries using input parameters 
* read from input file binaries.in (M1, M2, P, e, Z, Tmax). 
*
***
      implicit none
*
      INCLUDE 'const_bse.h'
*
      integer i,j,k,jj,nm1
      integer kw,kw2,kwx,kwx2,kstar(2)
      integer i1,i2,kdum
*
      real*8 m1,m2,tmax
      real*8 mass0(2),mass(2),z,zpars(20)
      real*8 epoch(2),tms(2),tphys,tphysf,dtp
      real*8 rad(2),lum(2),ospin(2)
      real*8 massc(2),radc(2),menv(2),renv(2)
      real*8 sep0,tb0,tb,ecc0,ecc,aursun,yeardy,yearsc,tol
      PARAMETER(aursun=214.95d0,yeardy=365.25d0,yearsc=3.1557d+07)
      PARAMETER(tol=1.d-07)
      real*8 t1,t2,mx,mx2,tbx,eccx



*
************************************************************************
* BSE parameters:
*
* neta is the Reimers mass-loss coefficent (neta*4x10^-13: 0.5 normally). 
* bwind is the binary enhanced mass loss parameter (inactive for single).
* hewind is a helium star mass loss factor (1.0 normally).
* alpha1 is the common-envelope efficiency parameter (1.0).  
* lambda is the binding energy factor for common envelope evolution (0.5).
*
* ceflag > 0 activates spin-energy correction in common-envelope (0). 
* tflag > 0 activates tidal circularisation (1).
* ifflag > 0 uses WD IFMR of HPE, 1995, MNRAS, 272, 800 (0). 
* wdflag > 0 uses modified-Mestel cooling for WDs (0). 
* bhflag > 0 allows velocity kick at BH formation (0). 
* nsflag > 0 takes NS/BH mass from Belczynski et al. 2002, ApJ, 572, 407 (1). 
* mxns is the maximum NS mass (1.8, nsflag=0; 3.0, nsflag=1). 
* idum is the random number seed used by the kick routine. 
*
* Next come the parameters that determine the timesteps chosen in each
* evolution phase:
*                 pts1 - MS                  (0.05) 
*                 pts2 - GB, CHeB, AGB, HeGB (0.01)
*                 pts3 - HG, HeMS            (0.02)
* as decimal fractions of the time taken in that phase.
*
* sigma is the dispersion in the Maxwellian for the SN kick speed (190 km/s). 
* beta is wind velocity factor: proportional to vwind**2 (1/8). 
* xi is the wind accretion efficiency factor (1.0). 
* acc2 is the Bondi-Hoyle wind accretion factor (3/2). 
* epsnov is the fraction of accreted matter retained in nova eruption (0.001). 
* eddfac is Eddington limit factor for mass transfer (1.0).
* gamma is the angular momentum factor for mass lost during Roche (-1.0). 
*
      neta = 0.5
      bwind = 0.0
      hewind = 1.0
      alpha1 = 1.0
      lambda = 0.5 
      ceflag = 0
      tflag = 1
      ifflag = 0 
      wdflag = 1 
      bhflag = 1
      nsflag = 1
      mxns = 2.5
      pts1 = 0.05
      pts2 = 0.01
      pts3 = 0.02
      sigma = 265.0
      beta = 0.125
      xi = 1.0 
      acc2 = 1.5
      epsnov = 0.001
      eddfac = 10.0
      gamma = -1.0
*
* Set the seed for the random number generator. 
*
      idum = 1764
      if(idum.gt.0) idum = -idum
*
* Set the collision matrix.
*
      CALL instar
*
* Open the input file - list of binary initial parameters. 
*
      OPEN(10,file='binariesBH1e6.in',status='unknown')
      READ(10,*)nm1
*
* Open the output files. 
*
      OPEN(11,file='binariesBH1e6.out',status='unknown')
*
      do i = 1,nm1
*
* Read in parameters and set coefficients which depend on metallicity. 
*
         READ(10,*)m1,m2,tb,ecc,z,tmax
         CALL zcnsts(z,zpars)
*        Resetting the log other arrays  
         bpp(:,:)=0
*


         M1presn = 0.0
         M2presn = 0.0
         Apresn = 0.0
         epresn = 0.0
         M1postsn = 0.0
         M2postsn = 0.0
         Apostsn = 0.0
         epostsn = 0.0
         Vkick = 0.0
         Vsystem = 0.0
         CE = 0


         ecc0 = ecc
*
* Initialize the binary. 
*
         kstar(1) = 1
         mass0(1) = m1
         mass(1) = m1
         massc(1) = 0.0
         ospin(1) = 0.0
         epoch(1) = 0.0
*
         kstar(2) = 1
         mass0(2) = m2
         mass(2) = m2
         massc(2) = 0.0
         ospin(2) = 0.0
         epoch(2) = 0.0
*

         tb0 = tb/yeardy
         sep0 = aursun*(tb0*tb0*(mass(1) + mass(2)))**(1.d0/3.d0)
         tb0 = tb



         tphys = 0.0
         tphysf = tmax
         dtp = 0.0
*
* Evolve the binary. 
*        
         CALL evolv2(kstar,mass0,mass,rad,lum,massc,radc,
     &               menv,renv,ospin,epoch,tms,
     &               tphys,tphysf,dtp,z,zpars,tb,ecc)
*
* Search the BCM array for the formation of binaries of 
* interest (data on unit 12 if detected) and also output 
* the final state of the binary (unit 11). 
*
* In this example we will search for XRB progenitors. 
*



C       bpp contains: '     TIME M1 M2 K1 K2 SEP ECC R1/ROL1 R2/ROL2  TYPE'
     
C       label(1) = 'INITIAL '
C       label(2) = 'KW CHNGE'
C       label(3) = 'BEG RCHE'
C       label(4) = 'END RCHE'
C       label(5) = 'CONTACT '
C       label(6) = 'COELESCE'
C       label(7) = 'COMENV  '
C       label(8) = 'GNTAGE  '
C       label(9) = 'NO REMNT'
C       label(10) = 'MAX TIME'
C       label(11) = 'DISRUPT '
C       label(12) = 'BEG SYMB'
C       label(13) = 'END SYMB'
C       label(14) = 'BEG BSS'


          do j=2,size(bpp,1)
            if (bpp(j,4).eq.14 .and. bpp(j-1,4).ne.14) then
                Tsn = bpp(j,1)               
                exit
            endif
          enddo


          do j=1,size(bpp,1)
            if (bpp(j,4).eq.14 .and. bpp(j,10).eq.3) then
                write(11,*) i,m1,m2,tb0,sep0, ecc0,
     &           M1presn, M2presn, Apresn, epresn,
     &           M1postsn, M2postsn, Apostsn, epostsn,
     &           Tsn, Vkick, Vsystem, CE, bpp(j,:)
                exit
            endif
          enddo





*
      enddo
*
      CLOSE(10)
      CLOSE(11)
*
************************************************************************
*
      STOP
      END
***
