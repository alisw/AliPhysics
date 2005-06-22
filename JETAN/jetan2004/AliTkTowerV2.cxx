//$Id$

#include <Riostream.h>
#include <TROOT.h>
#include <TClonesArray.h>
#include <TParticle.h>
#include <TMath.h>

#include "AliTkTowerV2.h"

#define Nparticles__ 100

ClassImp(AliTkTowerV2)

AliTkTowerV2::AliTkTowerV2() : TObject(), 
                         fEta(-999), fPhi(-999), fEt(-999), 
			 fNParticles(0),fEtCharged(0), fEtEM(0),
			 fUpdate(kFALSE), fParticles(0)
{
  fParticles = new TClonesArray("TParticle",Nparticles__);
}

AliTkTowerV2::AliTkTowerV2(const AliTkTowerV2 &t) : TObject(), 
                         fEta(-999), fPhi(-999), fEt(-999), 
			 fNParticles(0),fEtCharged(0), fEtEM(0),
			 fUpdate(kFALSE), fParticles(0)
{
  fEta = t.getEta();
  fPhi = t.getPhi();
  fEt = t.getEt();
  fEtCharged = 0;
  fEtEM = 0;
  fNParticles = t.getNParticles();
  if(fNParticles) fUpdate=kTRUE;
  else fUpdate=kFALSE;
  TClonesArray *otherParticles = t.getParticleList();
  TParticle *particle;
  fParticles = new TClonesArray("TParticle",Nparticles__);
  TIterator *iter = otherParticles->MakeIterator();
  Int_t i = 0;
  while ((particle = (TParticle *) iter->Next()) != NULL) {
    new ((*fParticles)[i]) TParticle(*particle);
    i++;
    if(i>fNParticles) 
      cerr << "TKTowerV2: should not happen " << i << " " << fNParticles << endl;
  }
  delete iter;
}

AliTkTowerV2::~AliTkTowerV2() 
{
  delete fParticles;
}

TClonesArray *AliTkTowerV2::getChargedParticleList() const 
{
  TClonesArray *chargedParticles = new TClonesArray("TParticle",Nparticles__);
  Int_t nChargedParticles = 0;
  TIterator *iter = getParticleList()->MakeIterator();
  TParticle *particle = NULL;
  while ((particle = (TParticle *) iter->Next()) != NULL) {
    if (isChargedParticle(particle)) {
      new ((*chargedParticles)[nChargedParticles]) TParticle(*particle);
      nChargedParticles++;
    }
  }
  delete iter;
  return chargedParticles;
}

TClonesArray *AliTkTowerV2::getNeutralParticleList() const 
{
  TClonesArray *neutralParticles = new TClonesArray("TParticle",Nparticles__);
  Int_t nNeutralParticles = 0;
  TIterator *iter = getParticleList()->MakeIterator();
  TParticle *particle = NULL;
  while ((particle = (TParticle *) iter->Next()) != NULL) {
    if (!isChargedParticle(particle)) {
      new ((*neutralParticles)[nNeutralParticles]) TParticle(*particle);
      nNeutralParticles++;
    }
  }
  delete iter;
  return neutralParticles;
}

void AliTkTowerV2::addParticle(const TParticle *particle) 
{
  new ((*fParticles)[fNParticles]) TParticle(*particle);
  fNParticles++;
  Float_t et_=TMath::Sqrt(particle->Pt()*particle->Pt());
  fEt+=et_;
  fUpdate=kTRUE;
}

void AliTkTowerV2::setParticleList(TClonesArray *ptr)
{
  delete fParticles;
  fParticles=ptr;
  fNParticles=ptr->GetEntries();
  fUpdate=kTRUE;
}

void AliTkTowerV2::update()  
{
  if(!fNParticles) return;
  if(!fUpdate) return;

  fEtEM=0;
  fEtCharged=0;

  TIterator *iter = fParticles->MakeIterator();
  TParticle *particle = NULL;
  while ((particle = (TParticle *) iter->Next()) != NULL) {
    Float_t et_=TMath::Sqrt(particle->Pt()*particle->Pt());
    if (isChargedParticle(particle))
      fEtCharged += et_;
    
    if (isEMParticle(particle))
      fEtEM += et_;
  }
  fUpdate=kFALSE;
}

void AliTkTowerV2::Clear(Option_t *option) 
{
  TObject::Clear(option);
  fParticles->Clear("C");
  fEta=-999;
  fPhi=-999; 
  fEt=-999; 
  fNParticles = 0;
  fUpdate=kFALSE;
  fEtCharged = 0;
  fEtEM = 0;
}

Bool_t AliTkTowerV2::isChargedParticle(const TParticle *particle) const 
{
  Bool_t isCharged = kFALSE;
  TParticle *part = new TParticle(*particle);
  TParticlePDG *partPDG = part->GetPDG(0);
  if (partPDG->Charge() != 0) {
    isCharged = kTRUE;
  }
  delete part;
  return isCharged;
}

Bool_t AliTkTowerV2::isEMParticle(const TParticle *particle) const 
{
  Bool_t isEM = kFALSE;
  // define electrons and gammas as EM particles...
  // not so sure if this is right...
  if ((particle->GetPdgCode() == 11) ||
      (particle->GetPdgCode() == -11) ||
      (particle->GetPdgCode() == 22)) {
    isEM = kTRUE;
  }
  return isEM;
}

Int_t AliTkTowerV2::Compare(const TObject *obj) const
{
  Double_t val=((AliTkTowerV2*)obj)->getEt();

  if(fEt>val) return -1; //qsort is ascending
  else if (fEt<val) return 1;
  else return 0;
}
