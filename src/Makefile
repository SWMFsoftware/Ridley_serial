#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission
#  For more information, see http://csem.engin.umich.edu/tools/swmf
SHELL =/bin/sh

include ../Makefile.def
include ../Makefile.conf
-include Makefile.DEPEND

#
#	Menu of make options
#
help:
	@echo ' '
	@echo '  You can "make" the following:'
	@echo ' '
	@echo '    LIB          (libIONO.a - ionosphere library)'
	@echo '    IONO         (Stand-Alone Ionosphere Model for Earth)'
	@echo '    PostIONO.exe (convert IDL files to TecPlot format)'
	@echo ' '
	@echo '    help         (Makefile option list)'
	@echo '    clean        (rm -f *~ *.o *.kmo *.mod *.T *.lst core)'
	@echo '    distclean    (make clean; rm ModSize.f90)'
	@echo '    backup       (tar -cf IONO.`date +%m_%d_%y_%H%M.tar` ...)'
	@echo ' '


MODULES = \
	ModSize.o\
	IE_ModMain.o\
	IE_ModIo.o\
	ModIeRlm.o\
	ModIonosphere.o\
	ModProcIE.o\
	ModIE_Interface.o\
	ModFiles.o\
	ModMagnit.o\
	ModConductance.o

#
#	Object Lists
#

OBJECTS = \
	IE_solve.o\
	IE_wrapper.o\
	IE_Initialize.o\
	IE_Library.o\
	IE_gather.o\
	create_auroral_oval.o\
	iono_currents.o\
	iono_joule.o\
	iono_solver.o\
	Create_Region2_Currents.o\
	ionosphere.o\
	facs_to_fluxes.o

#
#	Libraries
#

SEARCHDIR = -I${SHAREDIR} -I${LIBRARYDIR} -I${COUPLERDIR} -I${EMPIRICALIEDIR}

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCHDIR} ${MODULES} ${OBJECTS}

MY_LIB = ${LIBDIR}/libIE.a

LIB:	DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}:	${MF} ${MODULES} ${OBJECTS}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${MODULES} ${OBJECTS}

#
#	Stand-alone executables
#

PostIONO.exe:	PostIONO.o
	$(LINK.f90) -o ${BINDIR}/PostIONO.exe $(Lflag) PostIONO.o

IONO: ${EXE2}
	@echo Program IONO has been brought up to date.

backup:
	tar -cf IONO.`date +%m_%d_%y_%H%M.tar` Make* *.f90 *.h *.sub *.script *.llscript.* *.in mpi* perf_mjp.a BATio BATioPP* IOcat*

clean: cleanfiles

distclean: clean
	rm -f ModSize.f90

