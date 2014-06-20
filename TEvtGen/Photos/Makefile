include make.inc

LIB_VER = 1.0.0

#Name of libraries to create
LIB_PHOTOS_CXX_INT_SO = libPhotosCxxInterface.so
LIB_PHOTOS_CXX_INT_A  = libPhotosCxxInterface.a
LIB_PHOTOS_FORTRAN_SO = libPhotosFortran.so
LIB_PHOTOS_FORTRAN_A  = libPhotosFortran.a

PHOTOS_CXX_INT_OBJECTS = src/$(EVENT_RECORD_INTERFACE_DIR)/*.o \
                         src/$(C_PHOTOS_INTERFACE_DIR)/*.o \
                         src/$(UTILITIES_DIR)/*.o \
                         src/$(FORTRAN_PHOTOS_INTERFACE_DIR)/*.o

PHOTOS_FORTRAN_OBJECTS = src/$(PHOTOS_FORTRAN_DIR)/*.o

#directories containing source code
EVENT_RECORD_INTERFACE_DIR   = eventRecordInterfaces
FORTRAN_PHOTOS_INTERFACE_DIR = photosFortranInterfaces
C_PHOTOS_INTERFACE_DIR       = photosCInterfaces
PHOTOS_FORTRAN_DIR           = photos-fortran
UTILITIES_DIR                = utilities

##### Link objects to make library ######
all: include_dir $(EVENT_RECORD_INTERFACE_DIR) $(FORTRAN_PHOTOS_INTERFACE_DIR) $(C_PHOTOS_INTERFACE_DIR) $(PHOTOS_FORTRAN_DIR) $(UTILITIES_DIR)
	ar cr lib/$(LIB_PHOTOS_CXX_INT_A) $(PHOTOS_CXX_INT_OBJECTS)
	$(LD) $(LDFLAGS) $(SOFLAGS) -o lib/$(LIB_PHOTOS_CXX_INT_SO).$(LIB_VER) $(PHOTOS_CXX_INT_OBJECTS)
	ar cr lib/$(LIB_PHOTOS_FORTRAN_A) $(PHOTOS_FORTRAN_OBJECTS)
	$(LD) $(LDFLAGS) $(SOFLAGS) -o lib/$(LIB_PHOTOS_FORTRAN_SO).$(LIB_VER) $(PHOTOS_FORTRAN_OBJECTS)
	ln -sf $(LIB_PHOTOS_CXX_INT_SO).$(LIB_VER) lib/$(LIB_PHOTOS_CXX_INT_SO)
	ln -sf $(LIB_PHOTOS_FORTRAN_SO).$(LIB_VER) lib/$(LIB_PHOTOS_FORTRAN_SO)
	@echo "##################################################################"	
	@echo " Photos C++ Interface library created and moved to lib/ directory "
	@echo "##################################################################"
	@echo ""
	@echo "##################################################################"	
	@echo " To run examples, cd examples/ directory and there './configure'  "
	@echo " and 'make' again. Examples require Pythia8, ROOT and MC-Tester   "
	@echo "  installed. For details see examples/README.                     "
	@echo "##################################################################"

include_dir:
	mkdir -p include/Photos

####### Make object files ########
$(FORTRAN_PHOTOS_INTERFACE_DIR):
	make -C src/$(FORTRAN_PHOTOS_INTERFACE_DIR)
	cp src/$(FORTRAN_PHOTOS_INTERFACE_DIR)/*.h include/Photos

$(EVENT_RECORD_INTERFACE_DIR):
	make -C src/$(EVENT_RECORD_INTERFACE_DIR)
	cp src/$(EVENT_RECORD_INTERFACE_DIR)/*.h include/Photos

$(C_PHOTOS_INTERFACE_DIR):
	make -C src/$(C_PHOTOS_INTERFACE_DIR)
	cp src/$(C_PHOTOS_INTERFACE_DIR)/*.h include/Photos

$(UTILITIES_DIR):
	make -C src/$(UTILITIES_DIR)
	cp src/$(UTILITIES_DIR)/*.h include/Photos

$(PHOTOS_FORTRAN_DIR):
	make -C src/$(PHOTOS_FORTRAN_DIR)

install:
	mkdir -p $(PREFIX)/include/Photos
	cp include/Photos/* $(PREFIX)/include/Photos/.
	mkdir -p $(PREFIX)/lib
	cp lib/* $(PREFIX)/lib/.

clean:
	make clean -C src/$(EVENT_RECORD_INTERFACE_DIR)
	make clean -C src/$(FORTRAN_PHOTOS_INTERFACE_DIR)
	make clean -C src/$(C_PHOTOS_INTERFACE_DIR)
	make clean -C src/$(PHOTOS_FORTRAN_DIR)
	make clean -C src/$(UTILITIES_DIR)
	rm -f *~

Clean: clean
	rm -f lib/* include/Photos/*
	rm -rf config.log config.status autom4te.cache 
	rm -rf configure.paths.sh configure.paths.csh
	rm -f platform/make.inc make.inc
	rm -f examples/make.inc

make.inc:
	@echo ""
	@echo "Please execute ./configure first!"
	@echo "(Consider using 'source platform/afs.paths.sh' [or .csh] first)"
	@echo ""
	@false

always:
	@true
