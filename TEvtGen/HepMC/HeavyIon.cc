//--------------------------------------------------------------------------
//
// HeavyIon.cc
// Author:  Lynn Garren
//
// Implement operator >> and operator <<
//
// ----------------------------------------------------------------------

#include <iostream>
#include <ostream>
#include <istream>
#include <sstream>

#include "HepMC/HeavyIon.h"
#include "HepMC/StreamHelpers.h"
#include "HepMC/IO_Exception.h"

namespace HepMC {

/// Write the contents of HeavyIon to an output stream.
/// GenEvent stores a pointer to a HeavyIon.
std::ostream & operator << (std::ostream & os, HeavyIon const * ion)
{
    if ( !os ) {
	std::cerr << "HeavyIon output stream !os, "
		  << " setting badbit" << std::endl;
	os.clear(std::ios::badbit); 
	return os;
    }
    os << 'H';
    // HeavyIon* is set to 0 by default
    if ( !ion  ) {
	detail::output( os, 0 );
	detail::output( os, 0 );
	detail::output( os, 0 );
	detail::output( os, 0 );
	detail::output( os, 0 );
	detail::output( os, 0 );
	detail::output( os, 0 );
	detail::output( os, 0 );
	detail::output( os, 0 );
	detail::output( os, 0. );
	detail::output( os, 0. );
	detail::output( os, 0. );
	detail::output( os, 0. );
	detail::output( os,'\n');
	return os;
    }
    //
    detail::output( os, ion->Ncoll_hard() );
    detail::output( os, ion->Npart_proj() );
    detail::output( os, ion->Npart_targ() );
    detail::output( os, ion->Ncoll() );
    detail::output( os, ion->spectator_neutrons() );
    detail::output( os, ion->spectator_protons() );
    detail::output( os, ion->N_Nwounded_collisions() );
    detail::output( os, ion->Nwounded_N_collisions() );
    detail::output( os, ion->Nwounded_Nwounded_collisions() );
    detail::output( os, ion->impact_parameter() );
    detail::output( os, ion->event_plane_angle() );
    detail::output( os, ion->eccentricity() );
    detail::output( os, ion->sigma_inel_NN() );
    detail::output( os,'\n');

    return os;
}

/// Read the contents of HeavyIon from an input stream.
/// GenEvent stores a pointer to a HeavyIon.
std::istream & operator >> (std::istream & is, HeavyIon * ion)
{
    // make sure the stream is valid
    if ( !is ) { 
      std::cerr << "HeavyIon input stream setting badbit." << std::endl;
      is.clear(std::ios::badbit); 
      return is; 
    }
    // get the HeavyIon line
    std::string line;
    std::getline(is,line);
    std::istringstream iline(line);
    std::string firstc;
    iline >> firstc;
    // test to be sure the next entry is of type "H" 
    if( firstc != "H" ) { 
	std::cerr << "HeavyIon input stream invalid line type: " 
	          << firstc << std::endl;
	// The most likely problem is that we have found a HepMC block line
	throw IO_Exception("HeavyIon input stream encounterd invalid data");
    } 
    // read values into temp variables, then create a new HeavyIon object
    int nh =0, np =0, nt =0, nc =0, 
        neut = 0, prot = 0, nw =0, nwn =0, nwnw =0;
    float impact = 0., plane = 0., xcen = 0., inel = 0.; 
    iline >> nh ;
    if(!iline) throw IO_Exception("HeavyIon input stream encounterd invalid data");
    iline >> np ;
    if(!iline) throw IO_Exception("HeavyIon input stream encounterd invalid data");
    iline >> nt ;
    if(!iline) throw IO_Exception("HeavyIon input stream encounterd invalid data");
    iline >> nc ;
    if(!iline) throw IO_Exception("HeavyIon input stream encounterd invalid data");
    iline >> neut ;
    if(!iline) throw IO_Exception("HeavyIon input stream encounterd invalid data");
    iline >> prot;
    if(!iline) throw IO_Exception("HeavyIon input stream encounterd invalid data");
    iline >> nw ;
    if(!iline) throw IO_Exception("HeavyIon input stream encounterd invalid data");
    iline >> nwn ;
    if(!iline) throw IO_Exception("HeavyIon input stream encounterd invalid data");
    iline >> nwnw ;
    if(!iline) throw IO_Exception("HeavyIon input stream encounterd invalid data");
    iline >> impact ;
    if(!iline) throw IO_Exception("HeavyIon input stream encounterd invalid data");
    iline >> plane ;
    if(!iline) throw IO_Exception("HeavyIon input stream encounterd invalid data");
    iline >> xcen ;
    if(!iline) throw IO_Exception("HeavyIon input stream encounterd invalid data");
    iline >> inel;
    if(!iline) throw IO_Exception("HeavyIon input stream encounterd invalid data");
    if( nh == 0 ) {
        return is;
    }
    
    ion->set_Ncoll_hard(nh);
    ion->set_Npart_proj(np);
    ion->set_Npart_targ(nt);
    ion->set_Ncoll(nc);
    ion->set_spectator_neutrons(neut);
    ion->set_spectator_protons(prot);
    ion->set_N_Nwounded_collisions(nw);
    ion->set_Nwounded_N_collisions(nwn);
    ion->set_Nwounded_Nwounded_collisions(nwnw);
    ion->set_impact_parameter(impact);
    ion->set_event_plane_angle(plane);
    ion->set_eccentricity(xcen);
    ion->set_sigma_inel_NN(inel);

    return is;
}


} // HepMC
