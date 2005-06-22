// $Id$

/****************************************************************************
 * AliTkJetTrigger.cxx                                                         *
 * Test version of jet-trigger, includes many analysis functions to tune    *
 * paramters                                                                *
 * Thorsten Kollegger <kollegge@ikf.physik.uni-frankfurt.de>                *
 ***************************************************************************/

#include <Riostream.h>
#include <TFile.h>
#include <TROOT.h>
#include <TParticle.h>
#include <TClonesArray.h>
#include <TTree.h>
#include "AliTkChargedJetFinder.h"
#include "AliTkJetTriggerDecision.h"
#include "AliTkJetTriggerEvent.h"
#include "AliTkJetTrigger.h"

void AliTkJetTrigger::init() {
  // open file and create tree...
  if (getEvOutFilename() != NULL) {
    cout << "Writing events to " << getEvOutFilename() << endl;
    mEvOutFile = new TFile(getEvOutFilename(),"RECREATE");
    mEvOutEvent = new AliTkJetTriggerEvent();
    mEvOutTree = new TTree("jetTrigger","");
    mEvOutTree->Branch("event","AliTkJetTriggerEvent",&mEvOutEvent,32000,99);
  }
}

void AliTkJetTrigger::initEvent() {
  // can not be implemented, we need the particles...
}

void AliTkJetTrigger::initEvent(TClonesArray *particles, Int_t type) {
  setParticleType(type);
  initEvent(particles);
}

void AliTkJetTrigger::initEvent(TClonesArray *particles) {
  const UInt_t mParticlesDefaultSize = getExpectedParticleNumber();
  // make a copy of all accepted particles
  // just to be sure that they are not modified...

  clearParticles();
  if (!particles) {
    return;
  }
  mParticles = new TClonesArray("TParticle",mParticlesDefaultSize);
  TIterator *iter = particles->MakeIterator();
  TParticle *particle;
  UInt_t pos = 0;
  while ((particle = (TParticle *) iter->Next()) != NULL) {
    if (isParticleAccepted(particle)) {
      if(particle->Pt()>fParticleMinPt)
	new ((*mParticles)[pos++]) TParticle(*particle);
    }
  }
  delete iter;
  createSeedPoints();
}

void AliTkJetTrigger::make() {
  // fill the matrix containing trigger decision for all different cone radia,
  // nParticles and PtThresholds...
  // loop over all seedpoints...
  cout << "generated " << seedPoints.size() << " seedpoints" << endl;
  for (list<AliTkEtaPhiVector>::const_iterator seedIter = seedPoints.begin();
       seedIter != seedPoints.end(); ++seedIter) {
    AliTkEtaPhiVector coneCenter(*seedIter);
    for (list<Float_t>::const_iterator radiusIter = mConeRadia.begin();
	 radiusIter != mConeRadia.end(); ++radiusIter) {
      // get particles in cone
      list<TParticle*> *myParticles = getParticlesInCone(coneCenter,
							 *radiusIter);
      // loop over all PtThresholds...
      for (list<Float_t>::const_iterator PtThreshold =mConePtThreshold.begin();
	   PtThreshold != mConePtThreshold.end(); ++PtThreshold) {
	// loop over all nParticles thresholds...
	for (list<UInt_t>::const_iterator nParticles = mConeNParticles.begin();
	     nParticles != mConeNParticles.end(); ++nParticles) {
	  // what's the decision???
	  if (decide(myParticles,*PtThreshold,*nParticles)) {
	    AliTkJetTriggerDecision *d = new AliTkJetTriggerDecision();
	    d->setConeRadius(*radiusIter);
	    d->setPtThreshold(*PtThreshold);
	    d->setNParticles(*nParticles);
	    mEvOutEvent->addDecision(d);
	    // new change, requires testing 
	    delete d;
	  }
	} // end of nParticle threshold loop
      } // end of PtThreshold loop
      myParticles->erase(myParticles->begin(),myParticles->end());
      delete myParticles;
    } // end of cone radia loop
  } // end of seed point loop
}

void AliTkJetTrigger::finishEvent() {
  // let's write out all triggered settings...
  // not finish yet...
  clearParticles();
  if (mEvOutFile) {
    mEvOutTree->Fill();
  }
  mEvOutEvent->clear();
  seedPoints.erase(seedPoints.begin(),seedPoints.end());
}

void AliTkJetTrigger::finish() {
  // write object containing settings
  // close file
  if (mEvOutFile) {
    if (!mEvOutFile->IsOpen()) { cerr << "WHAT?" << endl; }
    mEvOutFile->cd();
    mEvOutTree->Write();
    mEvOutFile->Close();
  }
}

void AliTkJetTrigger::setDefaults() {
  fParticleMinPt=0;
  fParticleMinPtSeed=1;
  setParticleType(3);
  addConeRadius(0.7);
  addConeRadius(0.5);
  addConeRadius(0.3);
  addConeRadius(0.2);
  addConeRadius(0.1);
  addConeNParticles(1);
  addConeNParticles(2);
  addConeNParticles(3);
  addConeNParticles(4);
  addConeNParticles(5);
  addConePtThreshold(1.);
  addConePtThreshold(2.);
  addConePtThreshold(3.);
  addConePtThreshold(4.);
  addConePtThreshold(5.);

  setEvOutFilename("$JF_DATADIR/trigger.evout.root");
  setExpectedParticleNumber(30000);
}

void AliTkJetTrigger::setParticleType(Int_t type) {
  mParticleType = type;
  // current used particles types:
  // 1 == TParticle from TPythia6
  // 2 == TParticle from THijing
  // 3 == TParticle from TPythia6 + ALICE acceptance approximation
  // 4 == TParticle from THijing  + ALICE acceptance approximation
}

void AliTkJetTrigger::addConeRadius(Float_t radius) {
  mConeRadia.push_back(radius);
}

void AliTkJetTrigger::addConeNParticles(UInt_t nParticles) {
  mConeNParticles.push_back(nParticles);
}

void AliTkJetTrigger::addConePtThreshold(Float_t ptThreshold) {
  mConePtThreshold.push_back(ptThreshold);
}

void AliTkJetTrigger::mixParticles(TClonesArray *particles, Int_t type) {
  // mixes particles to the existing event
  if (!mParticles) {
    return;
  }
  if (!particles) {
    return;
  }
  TParticle *particle = NULL;
  UInt_t pos = 0;
  pos = mParticles->GetLast()+1;
  TIterator *iter = particles->MakeIterator();
  while ((particle = (TParticle *) iter->Next()) != NULL) {
    if (isParticleAccepted(particle,type)) {
      new ((*mParticles)[pos++]) TParticle(*particle);
    }
  }
  delete iter;
  // we need to generate new seedPoints...
  seedPoints.erase(seedPoints.begin(),seedPoints.end());
  createSeedPoints();
}

void AliTkJetTrigger::createSeedPoints() {
  TParticle *particle;
  TIterator *iter = mParticles->MakeIterator();
  while ((particle = (TParticle *)iter->Next()) != NULL) {
    // particles are checked for acceptance during initEvent, mixParticles
    //    if (isParticleAccepted(particle)) {
      if (particle->Pt() > fParticleMinPtSeed) {
	AliTkEtaPhiVector v(particle->Eta(),particle->Phi());
	seedPoints.push_back(v);
      }
      //    }
  }
  delete iter;
}

Bool_t AliTkJetTrigger::isParticleAccepted(TParticle *particle) {
  return isParticleAccepted(particle,getParticleType());
}

Bool_t AliTkJetTrigger::isParticleAccepted(TParticle *particle, Int_t type) {
  switch (type) {
  case 1:
    return isParticleAcceptedPythia(particle);
    break;
  case 2:
    return isParticleAcceptedHijing(particle);
    break;
  case 3:
    return isParticleAcceptedPythiaAliceGeo(particle);
    break;
  case 4:
    return isParticleAcceptedHijingAliceGeo(particle);
    break;
  };
  return kFALSE;
}

Bool_t AliTkJetTrigger::isParticleAcceptedPythia(TParticle *particle) {
  // check for stable particle
  if (particle->GetStatusCode() != 1) {
    return kFALSE;
  }
  return kTRUE;
}

Bool_t AliTkJetTrigger::isParticleAcceptedPythiaAliceGeo(TParticle *particle) {
  // check if valid PYTHIA particle
  if (!isParticleAcceptedPythia(particle)) {
    return kFALSE;
  }
  // check if particle is in |eta| < 1
  if (TMath::Abs(particle->Eta()) > 1) {
    return kFALSE;
  }
  // check if charged particle
  TParticlePDG *partPDG = particle->GetPDG();
  if (partPDG->Charge() == 0) {
    return kFALSE;
  }
  return kTRUE;
}

Bool_t AliTkJetTrigger::isParticleAcceptedHijing(TParticle */*particle*/) {
  return kTRUE;
}

Bool_t AliTkJetTrigger::isParticleAcceptedHijingAliceGeo(TParticle *particle) {
  if (!isParticleAcceptedHijing(particle)) {
    return kFALSE;
  }
  // check if particle is in |eta| < 1
  if (TMath::Abs(particle->Eta()) > 1) {
    return kFALSE;
  }
  // check if charged particle
  TParticlePDG *partPDG = particle->GetPDG();
  if (partPDG->Charge() == 0) {
    return kFALSE;
  }
  return kTRUE;
}

void AliTkJetTrigger::clearParticles() {
  // delete my particle copy...
  if (mParticles) {
    mParticles->Delete();
    delete mParticles;
  }
  mParticles = NULL;
}

list<TParticle*> *AliTkJetTrigger::getParticlesInCone(AliTkEtaPhiVector center,
						   Float_t radius) {
  return getParticlesInCone(center,radius,mParticles);
}

list<TParticle*> *AliTkJetTrigger::getParticlesInCone(AliTkEtaPhiVector center,
						   Float_t radius,
						   TClonesArray *particles) {
  // returns a list with all particles which are in a cone with radius around
  // center
  // one can easily speed up this part by saving the last cone+parameters and
  // calculate the new one from this one if there is some overlap
  list<TParticle*> *myParticles = new list<TParticle*>;
  Float_t radiusSq = radius*radius;
  TIterator *iter = particles->MakeIterator();
  TParticle *particle;
  while ((particle = (TParticle *)iter->Next()) != NULL) {
    AliTkEtaPhiVector pos(particle->Eta(),particle->Phi());
    if (center.diffSq(pos) < radiusSq) {
      myParticles->push_back(particle);
    }
  }
  delete iter;
  return myParticles;
}

Bool_t AliTkJetTrigger::decide(list<TParticle*> *particles, Float_t ptThreshold, 
			    UInt_t nParticles) {
  UInt_t n = 0;
  for (list<TParticle*>::const_iterator iter = particles->begin();
       iter != particles->end(); ++iter) {
    if ((*iter)->Pt() > ptThreshold) {
      n++;
    }
  }
  if (n >= nParticles) {
    return kTRUE;
  }
  return kFALSE;
}

void AliTkJetTrigger::setEvOutFilename(Char_t *filename) {
  if (!mEvOutName) {
    mEvOutName = new Char_t[4096];
  }
  strcpy(mEvOutName,filename);
}


ClassImp(AliTkJetTrigger)
