//////////////////////////////////////////////////////////////////////////
// Author: Lynn Garren
// search vectors for a GenParicle* instance
//////////////////////////////////////////////////////////////////////////

#include "HepMC/SearchVector.h"

namespace HepMC {

  
bool not_in_vector( std::vector<GenParticle*>* v, GenParticle* p )
{
    if( already_in_vector(v,p) == v->end()  ) return true;
    return false;
}

/// returns true if GenParticle is in the vector
std::vector<HepMC::GenParticle*>::iterator already_in_vector( std::vector<GenParticle*>* v, GenParticle* p )
{
    // if the vectors are mostly large, the search should be coded differently
    std::vector<GenParticle*>::iterator it;
    for( it = v->begin(); it != v->end(); ++it ) {
        // bail as soon as we find a match
        if( (*it) == p ) return it;
    }
    return v->end();
}

}	// HepMC
