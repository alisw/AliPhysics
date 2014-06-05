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

#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtHepMCEvent.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPDL.hh"

#include "HepMC/Units.h"

EvtHepMCEvent::EvtHepMCEvent() : 
  _theEvent(0), 
  _translation(0.0, 0.0, 0.0, 0.0)
{
}

EvtHepMCEvent::~EvtHepMCEvent() {
  this->deleteEvent();
}

void EvtHepMCEvent::deleteEvent() {

  if (_theEvent != 0) {
    _theEvent->clear();
    delete _theEvent; _theEvent = 0;
  }

}

void EvtHepMCEvent::constructEvent(EvtParticle* baseParticle) {

  EvtVector4R origin(0.0, 0.0, 0.0, 0.0);
  this->constructEvent(baseParticle, origin);

}

void EvtHepMCEvent::constructEvent(EvtParticle* baseParticle, EvtVector4R& translation) {

  // This class does not take ownership of the base particle pointer.
  // Rather, it uses the base particle to construct the event.

  this->deleteEvent();
  if (baseParticle == 0) {return;}

  _theEvent = new HepMC::GenEvent(HepMC::Units::GEV, HepMC::Units::MM);
  _translation = translation;

  // Use the recursive function addVertex to add a vertex with incoming/outgoing
  // particles. Adds a new vertex for any EvtParticles with decay daughters.
  // All particles are in the rest frame of the base particle ("lab frame").

  HepMC::GenParticle* hepMCGenParticle = this->createGenParticle(baseParticle, EvtHepMCEvent::LAB);

  this->addVertex(baseParticle, hepMCGenParticle);

}

HepMC::GenParticle* EvtHepMCEvent::createGenParticle(EvtParticle* theParticle, int frameType) {

  // Create an HepMC GenParticle, with the 4-momenta in the frame given by the frameType integer
  HepMC::GenParticle* genParticle = 0;

  if (theParticle != 0) {

    // Set the particle status integer to either stable or decayed
    int status(EvtHepMCEvent::STABLE);
    int nDaug = theParticle->getNDaug();
    if (nDaug > 0) {status = EvtHepMCEvent::DECAYED;}

    // Get the 4-momentum (E, px, py, pz) for the EvtParticle.
    EvtVector4R p4(0.0, 0.0, 0.0, 0.0);

    if (frameType == EvtHepMCEvent::RESTFRAME) {
      p4 = theParticle->getP4Restframe();
    } else if (frameType == EvtHepMCEvent::LAB) {
      p4 = theParticle->getP4Lab();
    } else {
      p4 = theParticle->getP4();
    }
  
    // Convert this to the HepMC 4-momentum
    double E = p4.get(0);
    double px = p4.get(1);
    double py = p4.get(2);
    double pz = p4.get(3);

    HepMC::FourVector hepMC_p4(px, py, pz, E);

    // Get the particle PDG integer id
    int PDGInt = EvtPDL::getStdHep(theParticle->getId());

    genParticle = new HepMC::GenParticle(hepMC_p4, PDGInt, status);

  }

  return genParticle;

}

void EvtHepMCEvent::addVertex(EvtParticle* inEvtParticle, HepMC::GenParticle* inGenParticle) {

  // This is a recursive function that adds GenVertices to the GenEvent for
  // the incoming EvtParticle and its daughters. We use two separate
  // pointers for the EvtParticle and GenParticle information: the former
  // to obtain the PDGId, 4-momenta, daughter and vertex positions, the latter to
  // set the incoming particle to the vertex. Note that the outgoing particle for
  // one vertex might be the incoming particle for another vertex - this needs to
  // be the same GenParticle pointer, hence the reason for using it as a 2nd argument
  // in this function.

  if (_theEvent == 0 || inEvtParticle == 0 || inGenParticle == 0) {return;}

  // Create the decay vertex
  HepMC::FourVector vtxCoord = this->getVertexCoord(inEvtParticle);
  HepMC::GenVertex* theVertex = new HepMC::GenVertex(vtxCoord);

  // Add the vertex to the event
  _theEvent->add_vertex(theVertex);

  // Set the incoming particle
  theVertex->add_particle_in(inGenParticle);

  // Set the outgoing particles (decay products)
  int nDaug = inEvtParticle->getNDaug();
  int iDaug(0);
  // Loop over the daughters
  for (iDaug = 0; iDaug < nDaug; iDaug++) {

    EvtParticle* evtDaughter = inEvtParticle->getDaug(iDaug);
    HepMC::GenParticle* genDaughter = this->createGenParticle(evtDaughter, EvtHepMCEvent::LAB);

    if (genDaughter != 0) {

      // Add a new GenParticle (outgoing) particle daughter to the vertex
      theVertex->add_particle_out(genDaughter);

      // Find out if the daughter also has decay products.
      // If so, recursively run this function again.
      int nDaugProducts = evtDaughter->getNDaug();

      if (nDaugProducts > 0) {
	  
	// Recursively process daughter particles and add their vertices to the event
	this->addVertex(evtDaughter, genDaughter);

      } // Have daughter products

    } // hepMCDaughter != 0

  } // Loop over daughters

}

HepMC::FourVector EvtHepMCEvent::getVertexCoord(EvtParticle* theParticle) {

  HepMC::FourVector vertexCoord(0.0, 0.0, 0.0, 0.0);

  if (theParticle != 0 && theParticle->getNDaug() != 0) {

    // Get the position (t,x,y,z) of the EvtParticle, offset by the translation vector.
    // This position will be the point where the particle decays. So we ask
    // the position of the (1st) daughter particle.
    EvtParticle* daugParticle = theParticle->getDaug(0);
    
    if (daugParticle != 0) {

      EvtVector4R vtxPosition = daugParticle->get4Pos() + _translation;

      // Create the HepMC 4 vector of the position (x,y,z,t)
      vertexCoord.setX(vtxPosition.get(1));
      vertexCoord.setY(vtxPosition.get(2));
      vertexCoord.setZ(vtxPosition.get(3));
      vertexCoord.setT(vtxPosition.get(0));

    }

  }

  return vertexCoord;

}
