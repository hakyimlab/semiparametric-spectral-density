# %Z%%M% %I% %E%
# This makefile is invoked by the command line
#	make -f makefile
# it will compile the fortran modules indicated by SOURCES into the object
# modules and produce an executable
#
# To remove all the objects but leave the executables use the command line
#	make -f makefile clean
#

.DEFAULT:
DBX=
#DBX= -g
# specify DBX=-g for dbx version

DESTDIR  = 
PROGRAM  =  prog
INCLUDES =
LDLIBS = 
#LDLIBS   = -L$(M3HOME)/lib/ioapi_22/Linux2_x86pg -lioapi \
	   -L$(M3HOME)/lib/netCDF/Linux -lnetcdf
	   

SOURCES  = xex3.f bsplpp.f bsplvb.f ppvalu.f interv.f

OBJECTS  = $(SOURCES:.f=.o)

FC       = g77 
FFLAGS   = $(DBX) 
#FFLAGS   = $(DBX) -e -fast -O4
LDFLAGS  = $(DBX)
#LDFLAGS  = $(DBX) -e -fast -O4
#CPPFLAGS = -DSUBST_GRID_ID='' \
           -DSUBST_IOPARMS=\'$(M3HOME)/lib/ioapi_22/PARMS3.EXT\' \
           -DSUBST_IOFDESC=\'$(M3HOME)/lib/ioapi_22/FDESC3.EXT\' \
           -DSUBST_IODECL=\'$(M3HOME)/lib/ioapi_22/IODECL3.EXT\' \
           -DSUBST_FILES_ID=\'FILES_CTM.EXT\' 
.KEEP_STATE:

##### beginning of dependency lines #####

all: $(INCLUDES) $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(FC) $(FFLAGS) -o $@ $(OBJECTS) $(LDFLAGS) $(LDLIBS)

clean: 
	-rm -f $(PROGRAM) $(OBJECTS) core






















