// $Id$

#include <Riostream.h>

#include <TParticle.h>
#include <TClonesArray.h>
#include <TH2F.h>

#include "AliJFParticlesCut.h"


AliJFParticlesCut::AliJFParticlesCut(TClonesArray *p) 
               : fPtMin(0),fPtMax(1000),
		 fEtaMin(-1),fEtaMax(1),
		 fPhiMin(0),fPhiMax(2*TMath::Pi()),
		 fNeutral(kTRUE),fCharged(kTRUE),fEM(kTRUE),
		 fParts(0)
{
  SetParticles(p);
}

Int_t AliJFParticlesCut::Cut()
{
  if(fParts==0) return -1;

  Int_t n=0;
  TParticle *particle;
  TIterator *iter = fParts->MakeIterator();
  while ((particle = (TParticle *) iter->Next()) != NULL) {
    if(IsAcceptedParticle(particle)) n++;
    else fParts->Remove(particle);
  } //end loop particles

  delete iter;
  return n;
}

Int_t AliJFParticlesCut::Cut(TClonesArray *p)
{
  SetParticles(p);
  return Cut();
}

TH2F* AliJFParticlesCut::CreateHistogram(Char_t *title,Char_t *text,Int_t phibins,Int_t etabins)
{
  if(fParts==0) return 0;

  TH2F *h=new TH2F(title,text,etabins,fEtaMin,fEtaMax,phibins,fPhiMin,fPhiMax);
  TParticle *particle;
  TIterator *iter = fParts->MakeIterator();
  while ((particle = (TParticle *) iter->Next()) != NULL) {
    if(IsAcceptedParticle(particle)){
      h->Fill(particle->Eta(),particle->Phi(),particle->Pt());
    }
  } //end loop particles

  //h->GetZaxis()->SetTitle("Pt");
  return h;
}

Bool_t AliJFParticlesCut::IsAcceptedParticle(TParticle *p)
{
#ifndef ALICEINTERFACE
  if(p->GetStatusCode()%100!=1) return kFALSE;
#endif

  Int_t pcode=p->GetPdgCode();  

  if ((pcode==11)||(pcode==-11)||(pcode==22)) {
    if(!fEM) return kFALSE;
  }  else {
    TParticlePDG *pdg=p->GetPDG();
    Float_t ch=pdg->Charge(); 
    if((!fCharged)&&(ch)) return kFALSE;
    if((!fNeutral)&&(!ch)) return kFALSE;
  }

  Float_t eta=p->Eta();
  if((eta<fEtaMin)||(eta>fEtaMax)) return kFALSE;

  Float_t phi=p->Phi();
  if((phi<fPhiMin)||(phi>fPhiMax)) return kFALSE;

  Float_t pt=p->Pt();
  if((pt<fPtMin)||(pt>fPtMax)) return kFALSE;

  return kTRUE;
}

void AliJFParticlesCut::SetPtCut(Float_t ptmin, Float_t ptmax)
{
  fPtMin=ptmin;
  fPtMax=ptmax;
}

void AliJFParticlesCut::SetPhiCut(Float_t phimin, Float_t phimax)
{
  fPhiMin=phimin;
  fPhiMax=phimax;
}

void AliJFParticlesCut::SetEtaCut(Float_t emin, Float_t emax)
{
  fEtaMin=emin;
  fEtaMax=emax;
}
