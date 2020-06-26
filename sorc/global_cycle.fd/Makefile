SHELL=	/bin/ksh
FC      =$(FCMP)
LIBS    =$(LIBSM) -L/scratch2/BMC/gsienkf/Tseganeh.Gichamo/libs -llapack
OBJS	=num_parthds.o machine.o cycle.o sfcsub_oida.o read_write_data.o snow_da_OI.o snow_analysis.o

CMD     =global_cycle
$(CMD):	$(OBJS)
	${FC} ${DEBUG} $(FFLAGS) $(LDFLG) $(OBJS) $(LIBS) -o $(CMD)

num_parthds.o:	num_parthds.f90
	$(FC) $(FFLAGS) $(OMPFLAG) -c num_parthds.f90

machine.o:	machine.f90
	${FC} $(FFLAGS) -c machine.f90

read_write_data.o:	read_write_data.f90
	${FC} $(FFLAGS) -c read_write_data.f90

snow_da_OI.o: snow_da_OI.f90
	$(FC)  $(FFLAGS) $(OMPFLAG) -c snow_da_OI.f90

snow_analysis.o: snow_da_OI.o snow_analysis.f90
	$(FC)  $(FFLAGS) $(OMPFLAG) -c snow_analysis.f90

sfcsub_oida.o: machine.o sfcsub_oida.F
	$(FC)  $(FFLAGS) $(OMPFLAG) -c sfcsub_oida.F

cycle.o: read_write_data.o snow_analysis.o cycle.f90
	$(FC)  $(FFLAGS) -c cycle.f90

install:
	-cp -p $(CMD) ../../exec/.
clean:
	-rm -f $(OBJS) *.mod $(CMD)
