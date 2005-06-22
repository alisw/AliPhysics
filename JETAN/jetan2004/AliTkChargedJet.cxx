// $Id$

#include <Riostream.h>
#include <TROOT.h>
#include <TParticle.h>
#include <TClonesArray.h>

#include "AliTkChargedJet.h"

//-------------------------------------------------------------------------
// implemenatation of AliTkChargedJet
//-------------------------------------------------------------------------
ClassImp(AliTkChargedJet)

AliTkChargedJet::AliTkChargedJet() : TObject() 
{
  fEta=-999; 
  fPhi=-999; 
  fPt=-999; 
  fParticles=new TClonesArray("TParticle",1000);
}

AliTkChargedJet::AliTkChargedJet(jet j) : TObject() 
{
  // save center of the jet (eg. of seed particle)
  AliTkEtaPhiVector center = j.getCentroid();
  fEta = center.Eta();
  fPhi = center.Phi();
  fPt = j.getPt();
  fNParticles = j.getNParticles();

  fParticles = j.getParticles();
  if (fParticles->GetEntries() != j.getNParticles()) {
    cerr << "AliTkChargedJet: cannot happen!" << endl;
  }
}

AliTkChargedJet::~AliTkChargedJet() 
{
 if(fParticles) {
   fParticles->Clear(); 
   delete fParticles;
 }
}

Float_t AliTkChargedJet::getPtInRadius(Float_t r) const
{
  Float_t pt = 0.0;
  Float_t rSq = r*r;
  AliTkEtaPhiVector center(this->getEta(), this->getPhi());
  TIterator *iter = this->fParticles->MakeIterator();
  TParticle *particle;
  while ((particle = (TParticle *) iter->Next()) != NULL) {
    AliTkEtaPhiVector v(particle->Eta(),particle->Phi());
    if (center.diffSq(v) < rSq) {
      pt += particle->Pt();
    }
  }
  return pt;
}

Int_t AliTkChargedJet::getParticlesInRadius(Float_t r) const
{
  Int_t part = 0;
  Float_t rSq = r*r;
  AliTkEtaPhiVector center(this->getEta(), this->getPhi());
  TIterator *iter = this->fParticles->MakeIterator();
  TParticle *particle;
  while ((particle = (TParticle *) iter->Next()) != NULL) {
    AliTkEtaPhiVector v(particle->Eta(),particle->Phi());
    if (center.diffSq(v) < rSq) {
      part ++;
    }
  }
  return part;
}

Int_t AliTkChargedJet::getNChargedInRadius(Float_t /*r*/) const
{
  return -1;
}

Int_t AliTkChargedJet::getNNeutralInRadius(Float_t /*r*/) const
{
  return -1;
}
