		SUBROUTINE triple_orbit(M1i, M2i, M3i, a1i, a2i, e1i, e2i,
     & M1preSN, M2preSN, M3preSN, a1preSN, a2preSN, e1preSN, e2preSN,
     & M1postSN,M2postSN,M3postSN,a1postSN,a2postSN,e1postSN,e2postSN)

		implicit none

		real*8 M1i, M2i, M3i, a1i, a2i, e1i, e2i
      real*8 M1preSN, M2preSN, M3preSN 
      real*8 a1preSN, a2preSN, e1preSN, e2preSN
      real*8 M1postSN, M2postSN, M3postSN
      real*8 a1postSN, a2postSN, e1postSN, e2postSN
      real*8 func, rtbis, Manom, Eanom, pi
      EXTERNAL func, rtbis

      real*8 r,DM, M0, Eanom_cr
      integer flag_disrupt



      flag_disrupt = 0
c     We assume that M2 AND M3 REMAIN CONSTANT DURING THE EVOLUTION.
      M2preSN = M2i
      M2postSN = M2i
      M3preSN = M3i
      M3postSN = M3i
c     We assume that tides are not important and that eccentricity in the
c     outer orbit remains constant until the SN
      e2preSN = e2i 

      a2preSN = a2i * (M1preSN + M2preSN + M3preSN)/(M1i + M2i + M3i) * (M1i/M1preSN)**2

c     Choose randomly a mean anomaly between 0 and 2pi
      CALL init_random_seed()         ! see example of RANDOM_SEED
      CALL RANDOM_NUMBER(Manom)
      pi = acos(-1.0d0)
      Manom = Manom*2.0d0*pi
c     Find the eccentric anomaly by finding the root of Manom-Eanom+ecc*sin(Eanom)     
      Eanom = rtbis(func,0.0d0,2.0d0*pi,1d-4,Manom,e2preSN)
      if (Eanom .gt. 2.0d0*pi) Eanom = Eanom - 2.0d0*pi
c     THe mass of the system before the explosion
      M0 = M1preSN+M2preSN+M3preSN
c     mass lost from the system during the explosion
      DM = M1preSN - M1postSN
c      Equation (3) from Hills (1983)

      Eanom_cr = acos((1.0d0-2.0d0*DM/M0)/e2preSN)
      if (Eanom .le. (2.0d0*pi-Eanom_cr) .or. Eanom .ge. Eanom_cr) flag_disrupt = 1
      
      r = a2preSN*(1.0d0-e2preSN*cos(Eanom))
c      Equation (1) from Hills (1983)
      a2postSN = a2preSN/2.0d0 * (M0-DM)/(M0/2.0d0-(a2preSN/r)*DM)
c      Equation (6) from Hills (1983)
      e2postSN = 1.0d0-(2.0d0*a2preSN/r)*(DM/M0)
      e2postSN = e2postSN/(1.0d0-DM/M0)**2
      e2postSN = 1.0d0 - (1.0d0-e2preSN**2)*e2postSN
      e2postSN = sqrt(e2postSN)

      if (flag_disrupt .eq. 1) then
         M1preSN = -2.0d0
         M2preSN = -2.0d0
         M3preSN = -2.0d0
         M1postSN= -2.0d0
         M2postSN= -2.0d0
         M3postSN= -2.0d0
         a1preSN = -2.0d0
         a2preSN = -2.0d0
         a1postSN= -2.0d0
         a2postSN= -2.0d0
         e1preSN = -2.0d0
         e2preSN = -2.0d0
         e1postSN= -2.0d0
         e2postSN= -2.0d0
      endif
	   END SUBROUTINE



! 		subroutine init_random_seed()
! 		use iso_fortran_env, only: int64
! 		implicit none
! 		integer, allocatable :: seed(:)
! 		integer :: i, n, un, istat, dt(8), pid
! 		integer(int64) :: t

! 		call random_seed(size = n)
! 		allocate(seed(n))
! 		! First try if the OS provides a random number generator
! 		open(newunit=un, file="/dev/urandom", access="stream", form="unformatted", action="read", status="old", iostat=istat)
! 		if (istat == 0) then
! 		   read(un) seed
! 		   close(un)
! 		else
! 		   ! Fallback to XOR:ing the current time and pid. The PID is
! 		   ! useful in case one launches multiple instances of the same
! 		   ! program in parallel.
! 		   call system_clock(t)
! 		   if (t == 0) then
! 		      call date_and_time(values=dt)
! 		      t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 
!             t = t + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 
!             t = t + dt(3) * 24_int64 * 60 * 60 * 1000 
! 	         t = t + dt(5) * 60 * 60 * 1000 
! 	         t = t + dt(6) * 60 * 1000 + dt(7) * 1000 
! 	         t = t + dt(8)
! 		   end if
! 		   pid = getpid()
! 		   t = ieor(t, int(pid, kind(t)))
! 		   do i = 1, n
! 		      seed(i) = lcg(t)
! 		   end do
! 		end if
! 		call random_seed(put=seed)
! 		contains
! 		! This simple PRNG might not be good enough for real work, but is
! 		! sufficient for seeding a better PRNG.
! 		function lcg(s)
! 		  integer :: lcg
! 		  integer(int64) :: s
! 		  if (s == 0) then
! 		     s = 104729
! 		  else
! 		     s = mod(s, 4294967296_int64)
! 		  end if
! 		  s = mod(s * 279470273_int64, 4294967291_int64)
! 		  lcg = int(mod(s, int(huge(0), int64)), kind(0))
! 		end function lcg
! 		end subroutine init_random_seed

    subroutine init_random_seed()
    implicit none
    integer, allocatable :: seed(:)
    integer :: i, n, un, istat, dt(8), pid
    integer(kind=16) :: t

    call random_seed(size = n)
    allocate(seed(n))
    ! First try if the OS provides a random number generator
    open(newunit=un, file="/dev/urandom", access="stream", form="unformatted", action="read", status="old", iostat=istat)
    if (istat == 0) then
       read(un) seed
       close(un)
    else
       ! Fallback to XOR:ing the current time and pid. The PID is
       ! useful in case one launches multiple instances of the same
       ! program in parallel.
       call system_clock(t)
       if (t == 0) then
          call date_and_time(values=dt)
          t = (dt(1) - 1970) * 365L * 24 * 60 * 60 * 1000 
            t = t + dt(2) * 31L * 24 * 60 * 60 * 1000 
            t = t + dt(3) * 24L * 60 * 60 * 1000 
           t = t + dt(5) * 60 * 60 * 1000 
           t = t + dt(6) * 60 * 1000 + dt(7) * 1000 
           t = t + dt(8)
       end if
       pid = getpid()
       t = ieor(t, int(pid, kind(t)))
       do i = 1, n
          seed(i) = lcg(t)
       end do
    end if
    call random_seed(put=seed)
    contains
    ! This simple PRNG might not be good enough for real work, but is
    ! sufficient for seeding a better PRNG.
    function lcg(s)
      integer :: lcg
      integer(int64) :: s
      if (s == 0) then
         s = 104729
      else
         s = mod(s, 4294967296_int64)
      end if
      s = mod(s * 279470273_int64, 4294967291_int64)
      lcg = int(mod(s, int(huge(0), int64)), kind(0))
    end function lcg
    end subroutine init_random_seed


      FUNCTION rtbis(func,x1,x2,xacc,Manom,ecc)
      INTEGER JMAX
      REAL*8 rtbis,x1,x2,xacc,func, Manom, ecc
      EXTERNAL func
      PARAMETER (JMAX=40)
      INTEGER j
      REAL*8 dx,f,fmid,xmid
      fmid=func(x2,Manom,ecc)
      f=func(x1,Manom,ecc)
      if(f*fmid.ge.0.) stop 'root must be bracketed in rtbis'
      if(f.lt.0.)then
        rtbis=x1
        dx=x2-x1
      else
        rtbis=x2
        dx=x1-x2
      endif
      do 11 j=1,JMAX
        dx=dx*.5
        xmid=rtbis+dx
        fmid=func(xmid,Manom,ecc)
        if(fmid.le.0.)rtbis=xmid
        if(abs(dx).lt.xacc .or. fmid.eq.0.) return
11    continue
      stop 'too many bisections in rtbis'
      END

      FUNCTION func(x,Manom,ecc)
         real*8 Manom,x,func,ecc
         func = Manom-x+ecc*sin(x)
      end
