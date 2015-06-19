***
      PROGRAM bse
***
*
* Evolves a binary by calling evolv2.f 
* (see header of subroutine for algorithm description). 
*
* Required input is described below. 
***
* See Tout et al., MNRAS, 1997, 291, 732 for a description of many of the
* processes in this code as well as the relevant references mentioned
* within the code.
* Updated reference is:
*           Hurley J.R., Tout C.A., & Pols O.R., 2002, MNRAS, 329, 897
* (please use this one).
***
* For single star evolution see Hurley, Pols & Tout, 2000, MNRAS, 315, 543.
* or Hurley, 2000, PhD Thesis, University of Cambridge (Chapter 2).
* The binary evolution algorithm is described in Chapter 3 of the thesis.
***
*
*           B I N A R Y
*           ***********
*
*       Roche lobe overflow.
*       --------------------
*
*       Developed by Jarrod Hurley, IOA, Cambridge.
*       .........................................................
*
*       Advice by Christopher Tout, Onno Pols & Sverre Aarseth.
*       ++++++++++++++++++++++++++++++++++++++++++++++++++
***
      implicit none
*
      INCLUDE 'const_bse.h'
*
      integer kw,kw2,kstar(2),j,k,time
*
      real*8 mass0(2),mass(2),z,zpars(20)
      real*8 epoch(2),tms(2),tphys,tphysf,dtp,aj
      real*8 rad(2),lum(2),ospin(2)
      real*8 massc(2),radc(2),menv(2),renv(2)
      real*8 tb,ecc,yearsc
      real*8 M1i, M2i, M3i, a1i, a2i, e1i, e2i
      real*8 M1preSN, M2preSN, M3preSN, TpreSN
      real*8 a1preSN, a2preSN, e1preSN, e2preSN
      real*8 M1postSN, M2postSN, M3postSN, TpostSN
      real*8 a1postSN, a2postSN, e1postSN, e2postSN
      CHARACTER*32 ArgName
      integer narg, i, flag


      PARAMETER(yearsc=3.1557d+07)
      CHARACTER*8 label(14)
*
************************************************************************
* Input:
*
* mass is in solar units.
* tphysf is the maximum evolution time in Myr.
* tb is the orbital period in days.
* kstar is the stellar type: 0 or 1 on the ZAMS - unless in evolved state. 
* z is metallicity in the range 0.0001 -> 0.03 where 0.02 is Population I.
* eccentricity can be anywhere in the range 0.0 -> 1.0.
*
* neta is the Reimers mass-loss coefficent (neta*4x10^-13: 0.5 normally). 
* bwind is the binary enhanced mass loss parameter (inactive for single).
* hewind is a helium star mass loss factor (1.0 normally).
* alpha1 is the common-envelope efficiency parameter (1.0).  
* lambda is the binding energy factor for common envelope evolution (0.5).
*
* ceflag > 0 activates spin-energy correction in common-envelope (0). #defunct#
* ceflag = 3 activates de Kool common-envelope model (0). 
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
* If you enter a negative kstar then parameters for an evolved star are
* required in the order of:
* current age, initial mass and spin rate, 
* otherwise the star will start on the ZAMS.
*
      OPEN(22,file='binary.in', status='old')
      READ(22,*)mass0(1),mass0(2),tphysf,tb,kstar(1),kstar(2),z,ecc
      READ(22,*)neta,bwind,hewind,alpha1,lambda
      READ(22,*)ceflag,tflag,ifflag,wdflag,bhflag,nsflag,mxns,idum
      READ(22,*)pts1,pts2,pts3
      READ(22,*)sigma,beta,xi,acc2,epsnov,eddfac,gamma
      if(kstar(1).lt.0.or.kstar(2).lt.0)then
         READ(22,*)tphys
         READ(22,*)aj,mass(1),ospin(1)
         epoch(1) = tphys - aj
         kstar(1) = ABS(kstar(1))
         READ(22,*)aj,mass(2),ospin(2)
         epoch(2) = tphys - aj
         kstar(2) = ABS(kstar(2))
      else

C       !Check if any arguments are found
      narg=command_argument_count()
C       !Loop over the arguments
      if(narg .ne. 7) stop "Wrong number of arguments"

      call get_command_argument(1,ArgName)
      read(ArgName,*) M1i
      mass0(1) = M1i
      call get_command_argument(2,ArgName)
      read(ArgName,*) M2i
      mass0(2) = M2i
      call get_command_argument(3,ArgName)
      read(ArgName,*) M3i
      call get_command_argument(4,ArgName)
      read(ArgName,*) a1i
      tb = sqrt(a1i**3)/(M1i+M2i) * 365.0d0
      call get_command_argument(5,ArgName)
      read(ArgName,*) a2i
      call get_command_argument(6,ArgName)
      read(ArgName,*) e1i
      ecc = e1i
      call get_command_argument(7,ArgName)
      read(ArgName,*) e2i



*
* Initialize the parameters.
* Set the initial spin of the stars. If ospin is zero (actually < 0.001)
* at time zero then evolv2 will set an appropriate ZAMS spin. If 
* ospin is greater than zero then it will start with that spin regardless
* of the time. If you want to start at time zero with negligible spin 
* then I suggest using a negligible value (but greater than 0.001).
* If ospin is negative then the stars will be in co-rotation with the orbit.
*
         tphys = 0.d0
         mass(1) = mass0(1)
         epoch(1) = 0.d0
         ospin(1) = 0.d0
         mass(2) = mass0(2)
         epoch(2) = 0.d0
         ospin(2) = 0.d0
      endif
      if(idum.gt.0) idum = -idum
      CLOSE(22)
*
* Note that this routine can be used to evolve a single star if you 
* simply set mass(2) = 0.0 or tb = 0.0 (setting both is advised as  
* well as some dummy value for ecc). 
*
************************************************************************
*
* Set parameters which depend on the metallicity 
*
      CALL zcnsts(z,zpars)
*
* Set the collision matrix.
*
      CALL instar
*
      label(1) = 'INITIAL '
      label(2) = 'KW CHNGE'
      label(3) = 'BEG RCHE'
      label(4) = 'END RCHE'
      label(5) = 'CONTACT '
      label(6) = 'COELESCE'
      label(7) = 'COMENV  '
      label(8) = 'GNTAGE  '
      label(9) = 'NO REMNT'
      label(10) = 'MAX TIME'
      label(11) = 'DISRUPT '
      label(12) = 'BEG SYMB'
      label(13) = 'END SYMB'
      label(14) = 'BEG BSS'
*
* Set the data-save parameter. If dtp is zero then the parameters of the 
* star will be stored in the bcm array at each timestep otherwise they 
* will be stored at intervals of dtp. Setting dtp equal to tphysf will 
* store data only at the start and end while a value of dtp greater than 
* tphysf will mean that no data is stored.
*
      dtp = 0.d0
*
* Evolve the binary.
* 
      CALL evolv2(kstar,mass0,mass,rad,lum,massc,radc,
     &            menv,renv,ospin,epoch,tms,
     &            tphys,tphysf,dtp,z,zpars,tb,ecc)


      flag = 0
      DO i=2, 80
         if (bpp(i,4) .eq. 14 .and. bpp(i-1,4) .ne. 14 .and. all(bpp(1:i,10) .ne. 3)) then
            TpostSN = bpp(i,1)
            M1postSN = bpp(i,2)
            M2postSN = bpp(i,3)
            a1postSN = bpp(i,6)/215.09539480759717d0
            e1postSN = bpp(i,7)
            flag = 1
            DO j = 2,50000
               if (bcm(j,2) .eq. 14 .and. bcm(j-1,2) .ne. 14) then
                  TpreSN = bcm(j-1,1)
                  M1preSN = bcm(j-1,4)
                  M2preSN = bcm(j-1,18)
                  a1preSN = bcm(j-1,31)/215.09539480759717d0
                  e1preSN = bcm(j-1,32)
                  exit
               endif
            ENDDO
            exit
         endif
      ENDDO


      if (flag .eq. 1) then
         call triple_orbit(M1i, M2i, M3i, a1i, a2i, e1i, e2i,
     &   M1preSN, M2preSN, M3preSN, a1preSN, a2preSN, e1preSN, e2preSN,
     &   M1postSN,M2postSN,M3postSN,a1postSN,a2postSN,e1postSN,e2postSN)

         if (a2postSN.lt.0.0d0) then
            TpreSN = -2.0d0
            TpostSN= -2.0d0
         endif

         write(*,*) TpreSN, M1preSN, M2preSN, M3preSN, 
     &               a1preSN, a2preSN, e1preSN, e2preSN
         write(*,*) TpostSN, M1postSN, M2postSN, M3postSN, 
     &               a1postSN, a2postSN, e1postSN, e2postSN
      else
         write(*,*) -1, -1, -1, -1, -1, -1, -1, -1
         write(*,*) -1, -1, -1, -1, -1, -1, -1, -1
      endif

      STOP
      END
***
