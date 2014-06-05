#ifndef SearchVector_h
#define SearchVector_h
// ----------------------------------------------------------------------
//
// SearchVector.h
// Author: Lynn Garren
//
// Utilities to search std::vector<GenParticle*> a GenParticle instance
// ----------------------------------------------------------------------

#include "HepMC/GenVertex.h"
#include "HepMC/GenParticle.h"

namespace HepMC {
  
/// returns true if it cannot find GenParticle* in the vector
bool not_in_vector( std::vector<HepMC::GenParticle*>*, GenParticle* );

/// Returns the index of a GenParticle* within a vector.
/// Returns -1 if GenParticle* is not in the vector.
std::vector<HepMC::GenParticle*>::iterator already_in_vector( std::vector<HepMC::GenParticle*>*, GenParticle* );

}	// HepMC

#endif // SearchVector_h
