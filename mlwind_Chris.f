***
      real*8 FUNCTION mlwind(kw,lum,r,mt,mc,rl,z)
      implicit none
      integer kw
      real*8 lum,r,mt,mc,rl,z
      real*8 dml,dms,dmt,p0,x,mew,lum0,kap,neta,bwind,hewind,mxns
      real*8 effTemp, whichWR_wind, fLBV, dms_LBV, dms_VKL, dms_WRl
      parameter(lum0=7.0d+04,kap=-0.5d0)
      common /value1/ neta,bwind,hewind,mxns

      effTemp =  1000.d0*((1130.d0*lum/(r**2.d0))**(1.d0/4.d0))

*	Wind prescription from Belczynski et al. (2010)

*	FOllowingHamman & Koesterke (1998)==1, following Nugis and Lamers(2000) ==2
      whichWR_wind = 1.d0      

	  fLBV  = 7.d0
	  if (whichWR_wind.eq.1.d0) fLBV  = 1.5d0

* Apply mass loss of Vink JS, de Koter A, Lamers, 2001
* For H-rich O/B stars

		dms = 0.d0
		dms_VKL = 0.d0
	   if(kw.ge.1.d0.and.kw.le.6.d0) then
	       if(effTemp.gt.25000.d0.and.effTemp.le.50000.d0)then
		     dms_VKL = -6.697d0+2.194d0*log10(lum/1.d5)
	    	 dms_VKL = dms_VKL-1.313d0*log10(mt/30.d0)
		     dms_VKL = dms_VKL-1.226d0*log10(2.6d0/2.d0)     
		     dms_VKL = dms_VKL+0.933d0*log10(effTemp/4.d4)
	    	 dms_VKL = dms_VKL-10.92d0*log10((effTemp/4.d4))**2d0     
	         dms_VKL = dms_VKL+0.85d0*log10(z/0.019d0)
		     dms_VKL = 10.d0**dms_VKL
	       else if(effTemp.ge.12500.0.and.effTemp.le.25000.0)then
		     dms_VKL = -6.688d0+2.210d0*log10(lum/1.d5)
	    	 dms_VKL = dms_VKL-1.339d0*log10(mt/30.d0)
		     dms_VKL = dms_VKL-1.601d0*log10(1.3d0/2.d0)     
		     dms_VKL = dms_VKL+1.07d0*log10(effTemp/2.d4)     
        	 dms_VKL = dms_VKL+0.85d0*log10(z/0.019d0)
		     dms_VKL = 10.d0**dms_VKL
		  endif
		dms = dms_VKL
		endif

	   ! LBV-like mass loss
	  dms_LBV = 0.D0
      x = 1.0d-5*r*sqrt(lum)
      if(lum.gt.6.0d+05.and.x.gt.1.0)then
         dms_LBV = fLBV*1.0d-4
         dms = dms_LBV
      endif


      dms_WRl = 0.d0  
      if (kw.ge.7.and.kw.le.9) then
		  if (whichWR_wind.eq.1.d0) then
    	     dms_WRl = 1.0d-13*lum**1.5*(z/0.019d0)**0.86
      	  else
	         dms_WRl = -5.73 + 0.88*log10(mt)
 	         dms_WRl = 10.d0**dms_WRl
      	endif
	    dms = dms_WRL      	  
	  endif
	
	  if(dms.le.0.d0) then
* If wind is null use old one Hurley et al (2000)

* Apply mass loss of Nieuwenhuijzen & de Jager, A&A, 1990, 231, 134,
* for massive stars over the entire HRD.
      dms = 0.d0
      if(lum.gt.4000.d0)then
         x = MIN(1.d0,(lum-4000.d0)/500.d0)
         dms = 9.6d-15*x*(r**0.81d0)*(lum**1.24d0)*(mt**0.16d0)
         dms = dms*(z/0.02d0)**(1.d0/2.d0)
      endif
      if(kw.ge.2.and.kw.le.9)then
* 'Reimers' mass loss
         dml = neta*4.0d-13*r*lum/mt
         if(rl.gt.0.d0) dml = dml*(1.d0 + bwind*(MIN(0.5d0,(r/rl)))**6)
* Apply mass loss of Vassiliadis & Wood, ApJ, 1993, 413, 641, 
* for high pulsation periods on AGB.
         if(kw.eq.5.or.kw.eq.6)then
            p0 = -2.07d0 - 0.9d0*log10(mt) + 1.94d0*log10(r)
            p0 = 10.d0**p0
            p0 = MIN(p0,2000.d0)
            dmt = -11.4d0+0.0125d0*(p0-100.d0*MAX(mt-2.5d0,0.d0))
            dmt = 10.d0**dmt
            dmt = 1.d0*MIN(dmt,1.36d-09*lum)
            dml = MAX(dml,dmt)
         endif
         if(kw.gt.6)then
            dms = MAX(dml,1.0d-13*hewind*lum**(3.d0/2.d0))
         else
            dms = MAX(dml,dms)
            mew = ((mt-mc)/mt)*MIN(5.d0,MAX(1.2d0,(lum/lum0)**kap))
* reduced WR-like mass loss for small H-envelope mass
            if(mew.lt.1.d0)then
               dml = 1.0d-13*lum**(3.d0/2.d0)*(1.d0 - mew)
               dms = MAX(dml,dms)
            end if
* LBV-like mass loss beyond the Humphreys-Davidson limit.
            x = 1.0d-5*r*sqrt(lum)
            if(lum.gt.6.0d+05.and.x.gt.1.d0)then
               dml = 0.1d0*(x-1.d0)**3*(lum/6.0d+05-1.d0)
               dms = dms + dml
            endif
         endif
      endif
	  endif

*
      mlwind = dms
*
      write (*,*) dms_VKL, dms_WRl, dms_LBV, kw, dms, fLBV

      return
      end
***
