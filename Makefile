IGNORE:
CMPLR = gfortran
FFLAGS = -O2 -ffixed-line-length-132 
LFLAGS = const_bse.h zdata.h 

.f.o:
	$(CMPLR) -c $(FFLAGS) $<

SRCE1 = \
sse.f deltat.f evolv1.f hrdiag.f kick.f mlwind.f mrenv.f \
ran3.f star.f zcnsts.f zfuncs.f

OBJT1 = $(SRCE1:.f=.o)

sse:    $(OBJT1) $(LFLAGS)
	$(CMPLR) $(FFLAGS) $(OBJT1) -o sse 

SRCE2 = \
bse.f comenv.f corerd.f deltat.f dgcore.f evolv2.f gntage.f \
hrdiag.f instar.f kick.f mix.f mlwind.f mrenv.f ran3.f rl.f \
star.f zcnsts.f zfuncs.f
 
OBJT2 = $(SRCE2:.f=.o)

bse:    $(OBJT2) $(LFLAGS)
	$(CMPLR) $(FFLAGS) $(OBJT2) -o bse 

SRCE3 = \
popbin.f comenv.f corerd.f deltat.f dgcore.f evolv2.f gntage.f \
hrdiag.f instar.f kick.f mix.f mlwind.f mrenv.f ran3.f rl.f \
star.f zcnsts.f zfuncs.f
 
OBJT3 = $(SRCE3:.f=.o)

popbin: $(OBJT3) $(LFLAGS)
	$(CMPLR) $(FFLAGS) $(OBJT3) -o popbin


SRCE4 = \
popbin_cl.f comenv.f corerd.f deltat.f dgcore.f evolv2.f gntage.f \
hrdiag.f instar.f kick.f mix.f mlwind.f mrenv.f ran3.f rl.f \
star.f zcnsts.f zfuncs.f
 
OBJT4 = $(SRCE4:.f=.o)



popbin_cl: $(OBJT4) $(LFLAGS)
	$(CMPLR) $(FFLAGS) $(OBJT4) -o popbin_cl


SRCE5 = \
sse_cl.f deltat.f evolv1.f hrdiag.f kick.f mlwind.f mrenv.f \
ran3.f star.f zcnsts.f zfuncs.f

OBJT5 = $(SRCE5:.f=.o)

sse_cl:    $(OBJT5) $(LFLAGS)
	$(CMPLR) $(FFLAGS) $(OBJT5) -o sse_cl



SRCE6 = \
bse_cl.f comenv.f corerd.f deltat.f dgcore.f evolv2.f gntage.f \
hrdiag.f instar.f kick.f mix.f mlwind.f mrenv.f ran3.f rl.f \
star.f zcnsts.f zfuncs.f triple_orbit.f
 
OBJT6 = $(SRCE6:.f=.o)

bse_cl:    $(OBJT6) $(LFLAGS)
	$(CMPLR) $(FFLAGS) $(OBJT6) -o bse_cl 


