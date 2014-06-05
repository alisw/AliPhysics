//--------------------------------------------------------------------------
#ifndef HEPMC_STREAM_HELPERS_H
#define HEPMC_STREAM_HELPERS_H

//////////////////////////////////////////////////////////////////////////
// garren@fnal.gov, March 2009
//
// This header contains helper functions used by streaming IO
//////////////////////////////////////////////////////////////////////////

#include <ostream>
#include <istream>

#include "HepMC/GenEvent.h"
#include "HepMC/TempParticleMap.h"

namespace HepMC {

namespace detail {

/// used by IO_GenEvent constructor
std::ostream & establish_output_stream_info( std::ostream & );
/// used by IO_GenEvent constructor
std::istream & establish_input_stream_info( std::istream & );

/// get a GenVertex from ASCII input
/// TempParticleMap is used to track the associations of particles with vertices
std::istream & read_vertex( std::istream &, TempParticleMap &, GenVertex * );

/// get a GenParticle from ASCII input
/// TempParticleMap is used to track the associations of particles with vertices
std::istream & read_particle( std::istream&, TempParticleMap &, GenParticle * );

/// write a double - for internal use by streaming IO
inline std::ostream & output( std::ostream & os, const double& d ) {
    if( os  ) {
	if ( d == 0. ) {
	    os << ' ' << (int)0;
	} else {
	    os << ' ' << d;
	}
    }
    return os;
}

/// write a float - for internal use by streaming IO
inline std::ostream & output( std::ostream & os, const float& d ) {
    if( os  ) {
	if ( d == 0. ) {
	    os << ' ' << (int)0;
	} else {
	    os << ' ' << d;
	}
    }
    return os;
}

/// write an int - for internal use by streaming IO
inline std::ostream & output( std::ostream & os, const int& i ) { 
    if( os  ) {
	if ( i == 0. ) {
	    os << ' ' << (int)0;
	} else {
            os << ' ' << i; 
	}
    }
    return os;
}

/// write a long - for internal use by streaming IO
inline std::ostream & output( std::ostream & os, const long& i ) {
    if( os  ) {
	if ( i == 0. ) {
	    os << ' ' << (int)0;
	} else {
            os << ' ' << i; 
	}
    }
    return os;
}

/// write a single char - for internal use by streaming IO
inline std::ostream & output( std::ostream & os, const char& c ) {
    if( os  ) {
	if ( c ) {
            os << c; 
	} else {
	    os << ' ' ;
	}
    }
    return os;
}

/// used to read to the end of a bad event
std::istream & find_event_end( std::istream & );

} // detail

} // HepMC

#endif  // HEPMC_STREAM_HELPERS_H
//--------------------------------------------------------------------------
