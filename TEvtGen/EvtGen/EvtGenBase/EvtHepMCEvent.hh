//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package. If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 2011      University of Warwick, UK
//
// Module: EvtHepMCEvent
//
// Description: Create an HepMC::GenEvent for the complete EvtParticle 
//              decay tree.
//
// Modification history:
//
//    John Back       June 2011            Module created
//
//------------------------------------------------------------------------

#ifndef EVTHEPMCEVENT_HH
#define EVTHEPMCEVENT_HH

#include "EvtGenBase/EvtVector4R.hh"

#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"
#include "HepMC/SimpleVector.h"

class EvtParticle;

class EvtHepMCEvent {

public:

  EvtHepMCEvent();
  virtual ~EvtHepMCEvent();

  // Select what frame a given GenParticle is in:
  // its own restframe, the lab frame (first mother), or its mother's frame
  enum HepMCFrame {RESTFRAME = 1, LAB = 2, MOTHER = 3};
  // Select the GenParticle status
  enum HepMCStatus {STABLE = 1, DECAYED = 2, HISTORY = 3};

  void constructEvent(EvtParticle* baseParticle);
  void constructEvent(EvtParticle* baseParticle, EvtVector4R& translation);
  
  HepMC::GenEvent* getEvent() {return _theEvent;}

  // Methods used to create GenParticles and FourVectors of vertices.
  // Make these public so that other classes may call them if they use EvtHepMCEvent.

  // Create a GenParticle using info from the EvtParticle, specifying what frame
  // the 4-momentum is from.
  HepMC::GenParticle* createGenParticle(EvtParticle* theParticle, int frameType);

  // Find out the decay vertex position for the given EvtParticle.
  HepMC::FourVector getVertexCoord(EvtParticle* theParticle);

protected:

private:

  // Delete the event structure (called by destructor)
  void deleteEvent();

  // Add a vertex to the event. This is called by the constructEvent function
  // and is recursive, i.e. it loops through all possible daughter particles and
  // their descendents.
  void addVertex(EvtParticle* inEvtParticle, HepMC::GenParticle* inGenParticle);

  HepMC::GenEvent* _theEvent;
  EvtVector4R _translation;

};

#endif
