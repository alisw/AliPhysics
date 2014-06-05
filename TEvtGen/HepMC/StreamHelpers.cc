//--------------------------------------------------------------------------
//
// GenEventStreamIO.cc
// Author:  Lynn Garren
//
// helper functions used by streaming IO
//
// ----------------------------------------------------------------------

#include <ostream>
#include <istream>
#include <sstream>

#include "HepMC/GenVertex.h"
#include "HepMC/GenParticle.h"
#include "HepMC/StreamHelpers.h"
#include "HepMC/IO_Exception.h"

namespace HepMC {

namespace detail {

std::istream & read_vertex( std::istream & is, 
                            TempParticleMap & particle_to_end_vertex, 
			    GenVertex * v )
{
    //
    // make sure the stream is valid
    if ( !is ) {
	std::cerr << "StreamHelpers::detail::read_vertex setting badbit." << std::endl;
	is.clear(std::ios::badbit); 
	return is;
    } 
    //
    // get the vertex line
    std::string line;
    std::getline(is,line);
    std::istringstream iline(line);
    std::string firstc;
    iline >> firstc;
    //
    // test to be sure the next entry is of type "V" 
    if ( firstc != "V"  ) {
	std::cerr << "StreamHelpers::detail::read_vertex invalid line type: " 
	          << firstc << std::endl;
	std::cerr << "StreamHelpers::detail::read_vertex setting badbit." << std::endl;
	is.clear(std::ios::badbit); 
	return is;
    } 
    // read values into temp variables, then create a new GenVertex object
    int identifier =0, id =0, num_orphans_in =0, 
        num_particles_out = 0, weights_size = 0;
    double x = 0., y = 0., z = 0., t = 0.; 
    iline >> identifier ;
    if(!iline) { throw IO_Exception("read_vertex input stream encounterd invalid data"); }
    iline >> id ;
    if(!iline) { throw IO_Exception("read_vertex input stream encounterd invalid data"); }
    iline >> x ;
    if(!iline) { throw IO_Exception("read_vertex input stream encounterd invalid data"); }
    iline >> y ;
    if(!iline) { throw IO_Exception("read_vertex input stream encounterd invalid data"); }
    iline >> z ;
    if(!iline) { throw IO_Exception("read_vertex input stream encounterd invalid data"); }
    iline >> t;
    if(!iline) { throw IO_Exception("read_vertex input stream encounterd invalid data"); }
    iline >> num_orphans_in ;
    if(!iline) { throw IO_Exception("read_vertex input stream encounterd invalid data"); }
    iline >> num_particles_out ;
    if(!iline) { throw IO_Exception("read_vertex input stream encounterd invalid data"); }
    iline >> weights_size;
    if(!iline) { throw IO_Exception("read_vertex input stream encounterd invalid data"); }
    WeightContainer weights(weights_size);
    for ( int i1 = 0; i1 < weights_size; ++i1 ) {
        iline >> weights[i1];
        if(!iline) { throw IO_Exception("read_vertex input stream encounterd invalid data"); }
    }
    v->set_position( FourVector(x,y,z,t) );
    v->set_id( id );
    v->weights() = weights;
    v->suggest_barcode( identifier );
    //
    // read and create the associated particles. outgoing particles are
    //  added to their production vertices immediately, while incoming
    //  particles are added to a map and handled later.
    for ( int i2 = 1; i2 <= num_orphans_in; ++i2 ) {
        GenParticle* p1 = new GenParticle( ); 
	detail::read_particle(is,particle_to_end_vertex,p1);
    }
    for ( int i3 = 1; i3 <= num_particles_out; ++i3 ) {
        GenParticle* p2 = new GenParticle( ); 
	detail::read_particle(is,particle_to_end_vertex,p2);
	v->add_particle_out( p2 );
    }

    return is;
}

std::istream & find_event_end( std::istream & is ) {
    // since there is no end of event flag, 
    // read one line at time until we find the next event 
    // or the end of event block
    // don't throw until we find the end of the event
    std::string line, firstc;
    while ( is ) { 
	is >> firstc;
        if( firstc=="E" ) {	// next event
	    is.unget();
            throw IO_Exception("input stream encountered invalid data");
	    return is;
	} else if( firstc.size() > 1 ) { // no more events in this block
            throw IO_Exception("input stream encountered invalid data, now at end of event block");
	    return is;
	}
        std::getline(is,line);
    }
    // the stream is bad 
    throw IO_Exception("input stream encountered invalid data, stream is now corrupt");
    return is;
}

} // detail

} // HepMC
