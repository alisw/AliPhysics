// EvtGen.h is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.
// Author: Philip Ilten.

// This file contains an EvtGen interface. HepMC and EvtGen must be enabled.

#ifndef Pythia8_EvtGen_H
#define Pythia8_EvtGen_H

#include "Pythia8/Pythia.h"
#include "EvtGen/EvtGen.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtParticleFactory.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtRandomEngine.hh"
#include "EvtGenBase/EvtDecayBase.hh"
#include "EvtGenExternal/EvtExternalGenList.hh"

using namespace Pythia8;

//==========================================================================

// A class to wrap the Pythia random number generator for use by EvtGen.

class EvtGenRandom : public EvtRandomEngine {

public:

  // Constructor.
  EvtGenRandom(Rndm *rndmPtrIn) {rndmPtr = rndmPtrIn;}

  // Return a random number.
  double random() {if (rndmPtr) return rndmPtr->flat(); else return -1.0;}

  // The random number pointer.
  Rndm *rndmPtr;

};

//==========================================================================

// A class to perform decays via the external EvtGen decay program,
// see http://evtgen.warwick.ac.uk/, the program manual provided with
// the EvtGen distribution, and D. J. Lange,
// Nucl. Instrum. Meth. A462, 152 (2001) for details.

// EvtGen performs a series of decays from some initial particle
// decay, rather than just a single decay, and so EvtGen cannot be
// interfaced through the standard external DecayHandler class without
// considerable complication. Consequently, EvtGen is called on the
// complete event record after all steps of Pythia are completed.

class EvtGenDecays {

public:

  // Constructor.
  EvtGenDecays(Pythia *pythiaPtrIn, string decayFile, string particleDataFile,
               EvtAbsRadCorr *isrPtr = 0, int mixing = 1, bool xml = false,
               bool limit = true);

  // Perform all decays.
  void decay();

  // Stop EvtGen decaying a particle and use Pythia unless "pythia = false".
  void exclude(int id, bool pythia = true) {
    if (!pythiaPtr) return; set<int>::iterator itr = ids.find(id);
    if (itr != ids.end()) {ids.erase(itr);
      if (pythia) pythiaPtr->particleData.mayDecay(id, true);}}

  // Update the Pythia particle database from EvtGen.
  void updatePythia();

  // Update the EvtGen particle database from Pythia.
  void updateEvtGen();

  // Read an EvtGen user decay file.
  void readDecayFile(string decayFile, bool xml = false) {
    evtgen.readUDecay(decayFile.c_str(), xml);}

  // List of additional models.
  EvtExternalGenList genlist;
  list<EvtDecayBase*> models;

protected:

  // The pointer to the associated Pythia object.
  Pythia *pythiaPtr;

  // Random number wrapper for EvtGen.
  EvtGenRandom rndm;

  // The EvtGen object.
  EvtGen evtgen;

  // Set of particle IDs to decay with EvtGen.
  set<int> ids;

  // Parameters used to check if a particle should decay (as set via Pythia).
  double tau0Max, tauMax, rMax, xyMax, zMax;
  bool limitTau0, limitTau, limitRadius, limitCylinder, limitDecay;

};

//--------------------------------------------------------------------------

// The constructor.

// The EvtGenDecays object is associated with a single Pythia
// instance. This is to ensure a consistent random number generator
// across the two, as well as any updates to particle data, etc. Note
// that if multiple EvtGenDecays objects exist, that they will modify
// one anothers particle databases due to the design of EvtGen.

// This constructor also sets all particles to be decayed by EvtGen as
// stable within Pythia. The parameters within Pythia used to check if
// a particle should be decayed, as described in the "Particle Decays"
// section of the Pythia manual, are set. Note that if the variable
// "limit" is set to "false", then no check will be made before
// decaying a particle with EvtGen.

// The constructor is designed to have the exact same form as the
// EvtGen constructor except for three differences.
// (1) The first variable is the pointer to the Pythia object.
// (2) The last variable is a flag to limit decays based on the Pythia
//     criteria.
// (3) No random engine pointer is passed, as this is obtained from
//     Pythia.

//   pythiaPtrIn:      the pointer to the associated Pythia object.
//   decayFile:        the name of the decay file to pass to EvtGen.
//   particleDataFile: the name of the particle data file to pass to EvtGen.
//   isrPtr:           the EvtAbsRadCorr pointer to pass to EvtGen.
//   mixing:           the mixing type to pass to EvtGen.
//   xml:              flag to use XML files to pass to EvtGen.
//   limit:            flag to limit particle decays based on Pythia criteria.

EvtGenDecays::EvtGenDecays(Pythia *pythiaPtrIn, string decayFile,
                           string particleDataFile, EvtAbsRadCorr *isrPtr,
                           int mixing, bool xml, bool limit) :
  models(genlist.getListOfModels()), pythiaPtr(pythiaPtrIn),
  rndm(&pythiaPtr->rndm), evtgen(decayFile.c_str(), particleDataFile.c_str(),
                                   &rndm, isrPtr, &models, mixing, xml) {

  if (!pythiaPtr) return;
  for (int iPrt = 0; iPrt < (int)EvtPDL::entries(); ++iPrt) {
    int id = EvtPDL::getStdHep(EvtPDL::getEntry(iPrt));
    ids.insert(id);
    pythiaPtr->particleData.mayDecay(id, false);
  }
  limitTau0     = pythiaPtr->settings.flag("ParticleDecays:limitTau0");
  tau0Max       = pythiaPtr->settings.parm("ParticleDecays:tau0Max");
  limitTau      = pythiaPtr->settings.flag("ParticleDecays:limitTau");
  tauMax        = pythiaPtr->settings.parm("ParticleDecays:tauMax");
  limitRadius   = pythiaPtr->settings.flag("ParticleDecays:limitRadius");
  rMax          = pythiaPtr->settings.parm("ParticleDecays:rMax");
  limitCylinder = pythiaPtr->settings.flag("ParticleDecays:limitCylinder");
  xyMax         = pythiaPtr->settings.parm("ParticleDecays:xyMax");
  zMax          = pythiaPtr->settings.parm("ParticleDecays:zMax");
  limitDecay    = limit && (limitTau0 || limitTau ||
                            limitRadius || limitCylinder);

}

//--------------------------------------------------------------------------

// The decay routine.

// All particles in the event record that can be decayed by EvtGen are
// decayed. The production vertex of each particle (which can also be
// obtained in EvtGen via EvtParticle::get4Pos()) is set by the decay
// vertex of its mother, which in turn is calculated from the mother's
// lifetime. The status code 93 is used to indicate an external decay.

void EvtGenDecays::decay() {

  // Loop over all particles in the Pythia event.
  Event &event = pythiaPtr->event;
  for (int iPro = 0; iPro < event.size(); ++iPro) {

    // Check particle is final and can be decayed by EvtGen.
    Particle* pyPro = &event[iPro];
    if (!pyPro->isFinal()) continue;
    if (ids.find(pyPro->idAbs()) == ids.end()) continue;

    // Perform the decay of the progenitor.
    EvtParticle *egPro = EvtParticleFactory::particleFactory
      (EvtPDL::evtIdFromStdHep(pyPro->id()),
       EvtVector4R(pyPro->e(), pyPro->px(), pyPro->py(), pyPro->pz()));
    egPro->setDiagonalSpinDensity();
    evtgen.generateDecay(egPro);
    if (egPro->getNDaug() == 0) {egPro->deleteTree(); continue;}
    pyPro->tau(egPro->getLifetime());

    // Add the decay tree to the event record.
    vector< pair<EvtParticle*, int> >
      moms(1, pair<EvtParticle*, int>(egPro, iPro));
    while (moms.size() != 0) {

      // Check if particle should decay.
      EvtParticle* egMom = moms.back().first;
      int          iMom  = moms.back().second;
      Particle*    pyMom = &event[iMom];
      moms.pop_back();
      if (limitDecay) {
        if (limitTau0 && pyMom->tau0() > tau0Max)    continue;
        else if (limitTau && pyMom->tau() > tauMax)  continue;
        else if (limitRadius && pow2(pyMom->xDec())
                 + pow2(pyMom->zDec()) > pow2(rMax)) continue;
        else if (limitCylinder && (pow2(pyMom->xDec()) + pow2(pyMom->yDec())
                                   > pow2(xyMax) || abs(pyMom->zDec()) > zMax))
          continue;
      }

      // Set the children of the mother.
      pyMom->daughters(event.size(), event.size() + egMom->getNDaug() - 1);
      pyMom->statusNeg();
      Vec4 vProd = pyMom->vDec();
      for (int iDtr = 0 ; iDtr < (int)egMom->getNDaug(); ++iDtr) {
        EvtParticle *egDtr = egMom->getDaug(iDtr);
        int          id    = egDtr->getPDGId();
        EvtVector4R  p     = egDtr->getP4Lab();
        int idx = event.append(id, 93, iMom, 0, 0, 0, 0, 0, p.get(1),
                               p.get(2), p.get(3), p.get(0), egDtr->mass());
        Particle *pyDtr = &event.back();
        pyDtr->vProd(vProd);
        pyDtr->tau(egDtr->getLifetime());
        if(egDtr->getNDaug() > 0)
          moms.push_back(pair<EvtParticle*, int>(egDtr, idx));
      }
    }
    egPro->deleteTree();
  }
}

//--------------------------------------------------------------------------

// Update the Pythia particle database from EvtGen.

// Note that only the particle spin type, charge type, nominal mass,
// width, minimum mass, maximum mass, and nominal lifetime are
// set. The name string is not set.

void EvtGenDecays::updatePythia() {
  if (!pythiaPtr) return;
  for (int entry = 0; entry < (int)EvtPDL::entries(); ++entry) {
    EvtId egid = EvtPDL::getEntry(entry);
    int   pyid = EvtPDL::getStdHep(egid);
    pythiaPtr->particleData.spinType  (pyid, EvtPDL::getSpinType(egid));
    pythiaPtr->particleData.chargeType(pyid, EvtPDL::chg3(egid));
    pythiaPtr->particleData.m0        (pyid, EvtPDL::getMass(egid));
    pythiaPtr->particleData.mWidth    (pyid, EvtPDL::getWidth(egid));
    pythiaPtr->particleData.mMin      (pyid, EvtPDL::getMinMass(egid));
    pythiaPtr->particleData.mMax      (pyid, EvtPDL::getMaxMass(egid));
    pythiaPtr->particleData.tau0      (pyid, EvtPDL::getctau(egid));
  }
}


//--------------------------------------------------------------------------

// Update the EvtGen particle database from Pythia.

// The particle update database between Pythia and EvtGen is not
// symmetric. Particularly, it is not possible to update the spin
// type, charge, or nominal lifetime in the EvtGen particle database.

void EvtGenDecays::updateEvtGen() {
  if (!pythiaPtr) return;
  int pyid = pythiaPtr->particleData.nextId(0);
  while (pyid != 0) {
    EvtId egid = EvtPDL::evtIdFromStdHep(pyid);
    EvtPDL::reSetMass   (egid, pythiaPtr->particleData.m0(pyid));
    EvtPDL::reSetWidth  (egid, pythiaPtr->particleData.mWidth(pyid));
    EvtPDL::reSetMassMin(egid, pythiaPtr->particleData.mMin(pyid));
    EvtPDL::reSetMassMax(egid, pythiaPtr->particleData.mMax(pyid));
    pyid = pythiaPtr->particleData.nextId(pyid);
  }
}

//==========================================================================

#endif // end Pythia8_EvtGen_H
