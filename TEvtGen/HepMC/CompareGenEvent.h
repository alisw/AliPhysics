//--------------------------------------------------------------------------
#ifndef HEPMC_COMPARE_GENEVENT_H
#define HEPMC_COMPARE_GENEVENT_H

/////////////////////////////////////////////////////////////////////////
// CompareGenEvent.h
//
// garren@fnal.gov, January 2008
// Free functions used to compare two copies of GenEvent
//////////////////////////////////////////////////////////////////////////
//

#include <iostream>

#include "HepMC/GenEvent.h"

namespace HepMC {

bool compareGenEvent( GenEvent*, GenEvent* );
bool compareSignalProcessVertex( GenEvent*, GenEvent* );
bool compareBeamParticles( GenEvent*, GenEvent* );
bool compareWeights( GenEvent*, GenEvent* );
bool compareVertices( GenEvent*, GenEvent* );
bool compareParticles( GenEvent*, GenEvent* );
bool compareVertex( GenVertex* v1, GenVertex* v2 );

} // HepMC

#endif  // HEPMC_COMPARE_GENEVENT_H
//--------------------------------------------------------------------------
