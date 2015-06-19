***
      PROGRAM sse
c-------------------------------------------------------------c
c
c     Evolves a single star.
c     Mass loss is an option.
c     The timestep is not constant but determined by certain criteria.
c
c     Written by Jarrod Hurley 26/08/97 at the Institute of
c     Astronomy, Cambridge.
c
c-------------------------------------------------------------c
c
c     STELLAR TYPES - KW
c
c        0 - deeply or fully convective low mass MS star
c        1 - Main Sequence star
c        2 - Hertzsprung Gap
c        3 - First Giant Branch
c        4 - Core Helium Burning
c        5 - First Asymptotic Giant Branch
c        6 - Second Asymptotic Giant Branch
c        7 - Main Sequence Naked Helium star
c        8 - Hertzsprung Gap Naked Helium star
c        9 - Giant Branch Naked Helium star
c       10 - Helium White Dwarf
c       11 - Carbon/Oxygen White Dwarf
c       12 - Oxygen/Neon White Dwarf
c       13 - Neutron Star
c       14 - Black Hole
c       15 - Massless Supernova
c
c-------------------------------------------------------------c
      implicit none
*
      INCLUDE 'const_bse.h'
*
      integer kw,j,k
*
      real*8 mass,mt,z,zpars(20)
      real*8 epoch,tms,tphys,tphysf,dtp
      real*8 r,lum,ospin
      real*8 mc,rc,menv,renv
      character*50 text1,text2,text3
      character*30 label(16)
      CHARACTER*32 ArgName
      integer narg, nbreaks
      data label /' Low Mass MS Star ',' Main sequence Star ',
     &            ' Hertzsprung Gap ',' Giant Branch ',
     &            ' Core Helium Burning ',
     &            ' First AGB ',' Second AGB ',
     &            ' Naked Helium MS ',' Naked Helium HG ',
     &            ' Naked Helium GB ',' Helium WD ',
     &            ' Carbon/Oxygen WD ',' Oxygen/Neon WD ',
     &            ' Neutron Star ',' Black Hole ',
     &            ' Massless Supernova '/
*
************************************************************************
* Input:
*
* mass is in solar units.
* z is metallicity in the range 0.0001 -> 0.03 where 0.02 is Population I.
* tphysf is the maximum evolution time in Myr.
*
* neta is the Reimers mass-loss coefficent (neta*4x10^-13; 0.5 normally).
* bwind is the binary enhanced mass loss parameter (inactive for single).
* hewind is a helium star mass loss factor (1.0 normally). 
* sigma is the dispersion in the Maxwellian for the SN kick speed (190 km/s). 
*
* ifflag > 0 uses WD IFMR of HPE, 1995, MNRAS, 272, 800 (0). 
* wdflag > 0 uses modified-Mestel cooling for WDs (0). 
* bhflag > 0 allows velocity kick at BH formation (0). 
* nsflag > 0 takes NS/BH mass from Belczynski et al. 2002, ApJ, 572, 407 (1). 
* mxns is the maximum NS mass (1.8, nsflag=0; 3.0, nsflag=1). 
* idum is the random number seed used in the kick routine. 
*
* Next come the parameters that determine the timesteps chosen in each
* evolution phase:
*                 pts1 - MS                  (0.05)
*                 pts2 - GB, CHeB, AGB, HeGB (0.01)
*                 pts3 - HG, HeMS            (0.02)
* as decimal fractions of the time taken in that phase.
*
* If you enter a negative mass then parameters for an evolved star are
* required in the order of:
* initial mass, current mass, type, current time & epoch,
* otherwise the star will start on the ZAMS.
*


C       !Check if any arguments are found
      narg=command_argument_count()
C       !Loop over the arguments
      if(narg .ne. 1) stop "Wrong number of arguments"

      call get_command_argument(1,ArgName)
      read(ArgName,*) mass



C       OPEN(22,file='evolve.in',status='old')
C       READ(22,*)mass,z,tphysf
C       READ(22,*)neta,bwind,hewind,sigma
C       READ(22,*)ifflag,wdflag,bhflag,nsflag,mxns
C       READ(22,*)pts1,pts2,pts3
	   z = 0.008
	   tphysf = 13800
	   neta = 0.5
	   bwind = 0.0
	   hewind = 0.5
	   sigma = 190.0
	   ifflag = 0
	   wdflag = 1
	   bhflag = 0
	   nsflag = 1
	   mxns = 3.0
	   pts1 = 0.05
	   pts2 = 0.01
	   pts3 = 0.02
	   idum = -3453


*
************************************************************************
*
* Set parameters which depend on the metallicity 
*
      CALL zcnsts(z,zpars)
      if(idum.gt.0) idum = -idum
*
      if(mass.gt.0.0)then
*
* Initialize the star to begin on the ZAMS.
*
         mt = mass
         kw = 1
         tphys = 0.d0
         epoch = 0.d0
      else
         stop "Mass must be positive"
      endif
*
* Set the initial spin of the star. If ospin is less than or equal to 
* zero at time zero then evolv1 will set an appropriate ZAMS spin. If 
* ospin is greater than zero then it will start with that spin regardless
* of the time. If you want to start at time zero with negligible spin 
* then I suggest using a negligible value (but greater than 0.001).
*
      ospin = 0.d0
*
* Set the data-save parameter. If dtp is zero then the parameters of the 
* star will be stored in the scm array at each timestep otherwise they 
* will be stored at intervals of dtp. Setting dtp equal to tphysf will 
* store data only at the start and end while a value of dtp greater than 
* tphysf will mean that no data is stored.
*
      dtp = 0.d0
* 
      CALL evolv1(kw,mass,mt,r,lum,mc,rc,menv,renv,ospin,
     &            epoch,tms,tphys,tphysf,dtp,z,zpars)
*
************************************************************************
* Output:
*
		j=2
		nbreaks=0
		do while (scm(j,1).ge.0.0)
		  if (scm(j,2) .ne. scm(j-1,2))  nbreaks=nbreaks+1
		  j=j+1
		enddo
		write(*,*) nbreaks

		j=2

		do while (scm(j,1).ge.0.0)
		  if (scm(j,2) .ne. scm(j-1,2))  write(*,*) scm(j,1)
		  j=j+1
		enddo







      j = 0
*
* The scm array stores the stellar parameters at the specified output 
* times. The parameters are (in order of storage):
*    
*    Time, stellar type, initial mass, current mass, log10(L), log10(r),
*    log10(Teff), core mass, epoch and spin.
*
      text1 = ' Tev(Myr)    type      Mo        Mt      log10(L) '
      text2 = ' log10(R) log10(Teff)  Mc        Menv     ' 
      text3 = ' epoch      spin' 
      WRITE(*,'(a,a,a)')text1,text2,text3
 30   j = j + 1
      if(scm(j,1).lt.0.0)then
         scm(j-1,1) = scm(j,1)
         j = j - 1
      endif
      if(scm(j,1).ge.0.0) then
	    WRITE(*,99)(scm(j,k),k=1,8),scm(j,10),scm(j,12),scm(j,13)
        goto 30
      endif
 99   FORMAT(8f10.4,1p,e12.4,0p,f12.4,1p,e12.4)
*
      STOP
      END
***
