// $Id$

#include <Riostream.h>

#include <TParticle.h>
#include <TClonesArray.h>
#include <TIterator.h>
#include <TMath.h>

#include "AliJFMCJet.h"

ClassImp(AliJFMCJet)

AliJFMCJet::AliJFMCJet(Int_t n) : AliJFJet(n)
{
}

AliJFMCJet::~AliJFMCJet()
{
  Clean();
}

void AliJFMCJet::AddParticle(TParticle *p)
{
  if(!p) return;

  new(fParticles[fN]) TParticle(*p);
  fN++;
  fIsUpdated=kFALSE;
}

void AliJFMCJet::Clean()
{
  fParticles.Delete();
}

void AliJFMCJet::Update()
{
  if(fIsUpdated) return;

  fIsUpdated=kTRUE;

  fPhi=0;
  fEta=0;
  fY=0;
  fPt=0;
  fPx=0;
  fPy=0;
  fPz=0;
  fE=0;
  fE_=0;

  //  Float_t fMaxParticlePt=0;

  TIterator *iter=fParticles.MakeIterator();
  TParticle *p;

  while((p=(TParticle*)iter->Next()) != NULL){

    Float_t px=p->Px();
    Float_t py=p->Py();
    Float_t pz=p->Pz();
    Float_t E=p->Energy();
    Float_t E_=TMath::Sqrt(px*px+py*py+pz*pz); //massless particles
    //    Float_t pt=p->Pt();

    fPx+=px;
    fPy+=py;
    fPz+=pz;
    fE+=E;
    fE_+=E_;
    /* to do: max particle, charge particle, etc.
    if(pt>fMaxParticlePt){
      fMaxParticlePt=pt;
      fMaxParticle=*p;
    }
    */
  } 

  fPt=TMath::Sqrt(fPx*fPx+fPy*fPy);
  fPhi=TMath::ATan(fPy/fPx);
  fEta=0.5*TMath::Log((fE+fPz)/(fE-fPz));
  fY=0.5*TMath::Log((fE_+fPz)/(fE_-fPz));
}

void AliJFMCJet::Debug()
{
  TIterator *iter=fParticles.MakeIterator();
  TParticle *p;

  Int_t i=0;
  while((p=(TParticle*)iter->Next()) != NULL){
    cout << i++ << ": " << p->Energy() << endl;
  }
}
