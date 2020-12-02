CXX =		g++
CC =		$(CXX)
CXXFLAGS =      -g -O -s -frounding-math
CXXFLAGS +=     -Wextra -Wall -Werror -fPIC
CXXFLAGS +=	-Wno-deprecated-declarations -Wno-unused-variable
CXXFLAGS +=	-Wno-unused-parameter -Wno-unused-function 
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

# We should be completely autarky of any external CGAL

ifeq ($(findstring fastjet-3.,$(FASTJET_ROOT)),fastjet-3.)
    CXXFLAGS +=	-I$(FASTJET_ROOT)/include
    LDLIBS +=	-Wl,-rpath,$(FASTJET_ROOT)/src/.libs \
		-Wl,-rpath,$(FASTJET_ROOT)/tools/.libs -lm \
		-L$(FASTJET_ROOT)/src/.libs \
		-L$(FASTJET_ROOT)/tools/.libs -lfastjettools -lfastjet
else
    CXXFLAGS +=	$(shell $(FASTJET_ROOT)/bin/fastjet-config --cxxflags)
    LDLIBS +=	$(shell $(FASTJET_ROOT)/bin/fastjet-config --libs)
endif
LDLIBS +=	-lboost_thread

CXXFLAGS +=	-I$(BOOST_ROOT)/include

CXXFLAGS +=     -I../..

TARGET =	ue_subtract_jet

all:		$(TARGET)

clean:
		/usr/bin/rm -f *~ *.o $(TARGET)

%.cc:		../../jet.h
