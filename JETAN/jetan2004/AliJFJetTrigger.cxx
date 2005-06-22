// $Id$

#include <Riostream.h>

#include <TClonesArray.h>
#include <TIterator.h>
#include <TParticle.h>

#include "AliJFJetTrigger.h"

ClassImp(AliJFJetTrigger)

AliJFJetTrigger::AliJFJetTrigger(Int_t n) : AliJFJetFinder(n)
{
  fParticles=new TClonesArray("TParticle",100000);
}

AliJFJetTrigger::~AliJFJetTrigger()
{
  delete fParticles;
}

Bool_t AliJFJetTrigger::IsAcceptedParticle(TParticle *p)
{
  if(p->GetStatusCode()%100!=1) return kFALSE;

  Int_t pcode=p->GetPdgCode();  

  if((!fEM) && ((pcode==11)||(pcode==-11)||(pcode==22))) return kFALSE;

  TParticlePDG *pdg=p->GetPDG();
  Float_t ch=pdg->Charge(); 
  if((!fCharged)&&(ch)) return kFALSE;
  if((!fNeutral)&&(!ch)) return kFALSE;

  Float_t eta=p->Eta();
  if((eta<fEtaMin)||(eta>fEtaMax)) return kFALSE;

  Float_t phi=p->Phi();
  if((phi<fPhiMin)||(phi>fPhiMax)) return kFALSE;

  Float_t pt=p->Pt();
  if((pt<fPtMin)||(pt>fPtMax)) return kFALSE;

  return kTRUE;
}

