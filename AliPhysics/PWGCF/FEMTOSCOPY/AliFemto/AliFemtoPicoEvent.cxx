////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoPicoEvent - internal AliFemto representation of the event.         //
// Created for the sake of minimizing the used size. Multiple instances are   //
// stored in memory when several mixing bins are used.                        //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
/***************************************************************************
 *
 * $Id$
 *
 * Author: Mike Lisa, Ohio State, lisa@mps.ohio-state.edu
 ***************************************************************************
 *
 * Description: part of STAR HBT Framework: AliFemtoMaker package
 *  PicoEvents are last-step ultra-compressed "events" just containing
 *  bare information about the particles of interest.  They have already
 *  gone through Event and Track cuts, so only Pair cuts are left.
 *  PicoEvents are *internal* to the code, and are stored in the
 *  Event-mixing buffers.
 *           
 *
 ***************************************************************************
 *
 * $Log$
 * Revision 1.1.2.1  2007/10/12 14:28:37  akisiel
 * New wave of cleanup and rule conformance
 *
 * Revision 1.1  2007/05/16 10:22:11  akisiel
 * Making the directory structure of AliFemto flat. All files go into one common directory
 *
 * Revision 1.2  2007/05/03 09:42:29  akisiel
 * Fixing Effective C++ warnings
 *
 * Revision 1.1.1.1  2007/04/25 15:38:41  panos
 * Importing the HBT code dir
 *
 * Revision 1.1.1.1  2007/03/07 10:14:49  mchojnacki
 * First version on CVS
 *
 * Revision 1.4  2000/07/16 21:38:23  laue
 * AliFemtoCoulomb.cxx AliFemtoSectoredAnalysis.cxx : updated for standalone version
 * AliFemtoV0.cc AliFemtoV0.h : some cast to prevent compiling warnings
 * AliFemtoParticle.cc AliFemtoParticle.h : pointers mTrack,mV0 initialized to 0
 * AliFemtoIOBinary.cc : some printouts in #ifdef STHBTDEBUG
 * AliFemtoEvent.cc : B-Field set to 0.25Tesla, we have to think about a better
 *                 solution
 *
 * Revision 1.3  2000/06/01 20:40:13  laue
 * AliFemtoIO.cc: updated for new V0s
 * AliFemtoPicoEvent.cc: collections especially cleared
 * franks1DHistoD.h, include changed from  <stdio> to <cstdio>
 * franks1DHistoD.cc, cout statement deleted
 *
 * Revision 1.2  2000/03/17 17:23:05  laue
 * Roberts new three particle correlations implemented.
 *
 * Revision 1.1.1.1  1999/06/29 16:02:57  lisa
 * Installation of AliFemtoMaker
 *
 **************************************************************************/

#include "AliFemtoPicoEvent.h"
#include "AliFemtoParticleCollection.h"

//________________
AliFemtoPicoEvent::AliFemtoPicoEvent() :
  fFirstParticleCollection(0),
  fSecondParticleCollection(0),
  fThirdParticleCollection(0)
{
  // Default constructor
  fFirstParticleCollection = new AliFemtoParticleCollection;
  fSecondParticleCollection = new AliFemtoParticleCollection;
  fThirdParticleCollection = new AliFemtoParticleCollection;
}
//_________________
AliFemtoPicoEvent::AliFemtoPicoEvent(const AliFemtoPicoEvent& aPicoEvent) :
  fFirstParticleCollection(0),
  fSecondParticleCollection(0),
  fThirdParticleCollection(0)
{
  // Copy constructor
  AliFemtoParticleIterator iter;

  fFirstParticleCollection = new AliFemtoParticleCollection;
  if (aPicoEvent.fFirstParticleCollection) {
    for (iter=aPicoEvent.fFirstParticleCollection->begin();iter!=aPicoEvent.fFirstParticleCollection->end();iter++){
      fFirstParticleCollection->push_back(*iter);
    }
  }
  fSecondParticleCollection = new AliFemtoParticleCollection;
  if (aPicoEvent.fSecondParticleCollection) {
    for (iter=aPicoEvent.fSecondParticleCollection->begin();iter!=aPicoEvent.fSecondParticleCollection->end();iter++){
      fSecondParticleCollection->push_back(*iter);
    }
  }
  fThirdParticleCollection = new AliFemtoParticleCollection;
  if (aPicoEvent.fThirdParticleCollection) {
    for (iter=aPicoEvent.fThirdParticleCollection->begin();iter!=aPicoEvent.fThirdParticleCollection->end();iter++){
      fThirdParticleCollection->push_back(*iter);
    }
  }
}
//_________________
AliFemtoPicoEvent::~AliFemtoPicoEvent(){
  // Destructor
  AliFemtoParticleIterator iter;
  
  if (fFirstParticleCollection){
    for (iter=fFirstParticleCollection->begin();iter!=fFirstParticleCollection->end();iter++){
      delete *iter;
    }
    fFirstParticleCollection->clear();
    delete fFirstParticleCollection;
    fFirstParticleCollection = 0;
  }
  
  if (fSecondParticleCollection){
    for (iter=fSecondParticleCollection->begin();iter!=fSecondParticleCollection->end();iter++){
      delete *iter;
    }
    fSecondParticleCollection->clear();
    delete fSecondParticleCollection;
    fSecondParticleCollection = 0;
  }

  if (fThirdParticleCollection){
    if (fThirdParticleCollection->size() != 0 ) {
      for (iter=fThirdParticleCollection->begin();iter!=fThirdParticleCollection->end();iter++){
	delete *iter;
      }
    }
    fThirdParticleCollection->clear();
    delete fThirdParticleCollection;
    fThirdParticleCollection = 0;
  }
}
//_________________
AliFemtoPicoEvent& AliFemtoPicoEvent::operator=(const AliFemtoPicoEvent& aPicoEvent) 
{
  // Assignment operator
  if (this == &aPicoEvent) 
    return *this;

  AliFemtoParticleIterator iter;
   
  if (fFirstParticleCollection){
      for (iter=fFirstParticleCollection->begin();iter!=fFirstParticleCollection->end();iter++){
	delete *iter;
      }
      fFirstParticleCollection->clear();
      delete fFirstParticleCollection;
      fFirstParticleCollection = 0;
  }

  if (fSecondParticleCollection){
    for (iter=fSecondParticleCollection->begin();iter!=fSecondParticleCollection->end();iter++){
      delete *iter;
    }
    fSecondParticleCollection->clear();
    delete fSecondParticleCollection;
    fSecondParticleCollection = 0;
  }

  if (fThirdParticleCollection){
    if (fThirdParticleCollection->size() != 0 ) {
      for (iter=fThirdParticleCollection->begin();iter!=fThirdParticleCollection->end();iter++){
	delete *iter;
      }
    }
    fThirdParticleCollection->clear();
    delete fThirdParticleCollection;
    fThirdParticleCollection = 0;
  }

  fFirstParticleCollection = new AliFemtoParticleCollection;
  if (aPicoEvent.fFirstParticleCollection) {
    for (iter=aPicoEvent.fFirstParticleCollection->begin();iter!=aPicoEvent.fFirstParticleCollection->end();iter++){
      fFirstParticleCollection->push_back(*iter);
    }
  }
  fSecondParticleCollection = new AliFemtoParticleCollection;
  if (aPicoEvent.fSecondParticleCollection) {
    for (iter=aPicoEvent.fSecondParticleCollection->begin();iter!=aPicoEvent.fSecondParticleCollection->end();iter++){
      fSecondParticleCollection->push_back(*iter);
    }
  }
  fThirdParticleCollection = new AliFemtoParticleCollection;
  if (aPicoEvent.fThirdParticleCollection) {
    for (iter=aPicoEvent.fThirdParticleCollection->begin();iter!=aPicoEvent.fThirdParticleCollection->end();iter++){
      fThirdParticleCollection->push_back(*iter);
    }
  }

  return *this;
}

