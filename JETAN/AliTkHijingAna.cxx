// $Id$

#include <Riostream.h>
#include <TROOT.h>
#include <TParticle.h>
#include <TClonesArray.h>

#include "AliTkChargedJetFinder.h"
#include "AliTkHijingAna.h"

void AliTkHijingAna::setCenter(Float_t eta, Float_t phi) {
  center.setEta(eta);
  center.setPhi(phi);
  // new center makes particle list invalid
  clearCone();
}

Float_t AliTkHijingAna::getCenterEta() const {
  return center.Eta();
}

Float_t AliTkHijingAna::getCenterPhi() const {
  return center.Phi();
}

Bool_t AliTkHijingAna::isParticleInRadius(TParticle *particle,Float_t radius) {
  AliTkEtaPhiVector p(particle->Eta(),particle->Phi());
  if (center.diff(p) < radius) {
    return kTRUE;
  } else {
    return kFALSE;
  }
}

Bool_t AliTkHijingAna::isParticleAccepted(TParticle *particle,Float_t pt,Bool_t neutral) {

  //check if particle was final  
  UInt_t status = 0;
  status = (particle->GetStatusCode() % 100);
  if (status != 1) {
    return kFALSE;
  }

  // if !neutral accept only charged particles, else accept everything
  if(particle->Pt()<pt) return kFALSE;
  if(!neutral){
    TParticlePDG *partPDG = particle->GetPDG();
    if (partPDG->Charge() == 0) {
      return kFALSE;
    }
  }

  return kTRUE;
}

Bool_t AliTkHijingAna::isParticleAcceptedALICE(TParticle *particle,Float_t /*pt*/,Bool_t /*neutral*/) {
  // fake ALICE acceptance, only particles in |eta|<0.9

  if (/*(isParticleAccepted(particle,pt,neutral)) &&*/ (TMath::Abs(particle->Eta()) < 0.9)) {
    return kTRUE;
  }

  return kFALSE;
}

void AliTkHijingAna::clear() {
  // clear TClonesArrays
  if (oParticles) {
    oParticles->Delete();
    delete oParticles;
    oParticles = NULL;
  }
  if (mParticles) {
    mParticles->Delete();
    delete mParticles;
    mParticles = NULL;
  }
  mRadius = -1;
  updated = kFALSE;
}

void AliTkHijingAna::clearCone() {
  if (mParticles) {
    mParticles->Delete();
    delete mParticles;
    mParticles = NULL;
  }
  mRadius = -1;
  updated = kFALSE;
}

void AliTkHijingAna::setParticles(TClonesArray *particles) {
  // take particles and dont apply any further cuts
  if (!particles) {
    return;
  }

  clear();
  oParticles = new TClonesArray("TParticle",100000);
  TIterator *iter = particles->MakeIterator();
  TParticle *particle;
  UInt_t pos = 0;
  while ((particle = (TParticle *) iter->Next()) != NULL) {
    new ((*oParticles)[pos++]) TParticle(*particle);
  }
  delete iter;

}

void AliTkHijingAna::setParticles(TClonesArray *particles,Float_t pt,Bool_t neutral) {
  // we are only interested in (charged) particles, to reduce computing time
  // create a copy which contains only the charged particles
  if (!particles) {
    return;
  }
  // new particle list makes old one invalid...
  clear();
  oParticles = new TClonesArray("TParticle",100000);
  TIterator *iter = particles->MakeIterator();
  TParticle *particle;
  UInt_t pos = 0;
  while ((particle = (TParticle *) iter->Next()) != NULL) {
    if (isParticleAccepted(particle,pt,neutral)) {
      new ((*oParticles)[pos++]) TParticle(*particle);
    }
  }
  delete iter;
}

TClonesArray *AliTkHijingAna::getParticlesInRadius(Float_t radius) {
  return getParticlesInRadius(NULL,radius);
}

TClonesArray *AliTkHijingAna::getParticlesInRadius(TClonesArray *particles, 
						Float_t radius) {
  // this creates a TClonesArray with all charged particles in radius
  // and ALICE acceptance

  // uses the list of all charged particles (if created)
  // OR
  // used a previous list if radius < mRadius

  if (updated && (mParticles != NULL) && (radius < mRadius)) {
    particles = mParticles;
  } else {
    if (oParticles) {
      particles = oParticles;
    }
  }
  if (!particles) {
    return NULL;
  }

  mRadius = radius;
  TClonesArray *MyParticles = new TClonesArray("TParticle",100000);
  TParticle *particle;

  TIterator *iter = particles->MakeIterator();
  UInt_t pos = 0;
  while ((particle = (TParticle *) iter->Next()) != NULL) {
    if (/*isParticleAcceptedALICE(particle) &&*/
	isParticleInRadius(particle,mRadius)) {
      new ((*MyParticles)[pos++]) TParticle(*particle);
    }
  }

  // clean up
  delete iter;
  if (mParticles) {
    mParticles->Delete();
    delete mParticles;
  }
  mParticles = MyParticles;
  updated = kTRUE;
  
  return mParticles;
}

Int_t AliTkHijingAna::getNParticlesInRadius(Float_t ptCut) {
  if (!updated) {
    // must call getParticlesInRadius before!
    return -1;
  }
  // returns number of particles in cone with pt>ptCut
  TIterator *iter = mParticles->MakeIterator();
  TParticle *particle;
  Int_t nPart = 0;
  while ((particle = (TParticle *) iter->Next()) != NULL) {
    if (particle->Pt() > ptCut) {
      nPart++;
    }
  }
  delete iter;
  return nPart;
}

Float_t AliTkHijingAna::getEtInRadius(Float_t ptCut) {
  if (!updated) {
    // must call getParticlesInRadius before!
    return -1;
  }
  // returns Et in cone using massless particles with pT>ptCut
  TIterator *iter = mParticles->MakeIterator();
  TParticle *particle;
  Float_t Et = 0;
  while((particle = (TParticle*) iter->Next()) != NULL) {
    if (particle->Pt() > ptCut) {
      Et += particle->Pt(); // assume massless particle...
    }
  }
  delete iter;
  return Et;
}  

Float_t AliTkHijingAna::getEtInRadius(TClonesArray *particles,Float_t radius) {
  TIterator *iter = particles->MakeIterator();
  TParticle *particle;
  Float_t Et = 0;
  while((particle = (TParticle*) iter->Next()) != NULL) {
    if (/*isParticleAccepted(particle) &&*/
	isParticleInRadius(particle,radius)) {
      Et += particle->Pt(); // assume massless particle...
    }
  }
  delete iter;
  return Et;
}


ClassImp(AliTkHijingAna)
