CXX =		g++
CC =		$(CXX)
CXXFLAGS =      -g -O -s -frounding-math
CXXFLAGS +=     -Wextra -Wall -Werror -fPIC
CXXFLAGS +=	-Wno-unused-variable -Wno-unused-parameter
CXXFLAGS +=	-Wno-unused-function
ifeq ($(findstring g++,$(CXX)),g++)
    CXXFLAGS += -no-canonical-prefixes -Wno-builtin-macro-redefined \
                -D__DATE__="redacted" -D__TIMESTAMP__="redacted" \
                -D__TIME__="redacted" -U_FORTIFY_SOURCE \
                -D_FORTIFY_SOURCE=1 -fstack-protector
endif
CXXFLAGS +=     -std=c++11
CXXFLAGS +=	$(shell root-config --cflags)
LDFLAGS =	$(shell root-config --ldflags)
LDLIBS =	$(shell root-config --libs)

CXXFLAGS +=	-I$(BOOST_ROOT)/include

ifneq (,$(HDF5_DIR))
    CXXFLAGS +=	-I$(HDF5_DIR)/include
    LDFLAGS +=	-L$(HDF5_DIR)/lib64 -L$(HDF5_DIR)/lib
    LDLIBS +=	-lhdf5_cpp -lhdf5
endif

CXXFLAGS +=	-I../.. -I.

TARGET =	convert_sample convert_sample_full inference \
		update_ntuple eta_data

all:		check-nersc $(TARGET)

check-nersc:
ifeq (,$(ROOT_DIR))
		$(error At NERSC, run \`module load boost cray-hdf5 \
			root\' first)
endif
ifeq (,$(BOOST_DIR))
		$(error At NERSC, run \`module load boost cray-hdf5 \
			root\' first)
endif
ifeq (,$(HDF5_DIR))
		$(error At NERSC, run \`module load boost cray-hdf5 \
			root\' first)
endif

clean:
		/usr/bin/rm -f *~ *.o $(TARGET) ml_out.py

%.cc:		../../emcal.h ../../special_function.h \
		AliVCaloCells.h AliVCluster.h AliVVZERO.h

inference:	inference.cc ../../keras_model.cc

eta_data:	eta_data.cc ../../keras_model.cc

update_ntuple:	update_ntuple.cc ../../keras_model.cc
