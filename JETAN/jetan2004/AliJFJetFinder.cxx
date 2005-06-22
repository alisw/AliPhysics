// $Id$

#include <Riostream.h>

#include <TClonesArray.h>
#include <TIterator.h>
#include <TParticle.h>

#include "AliJFJetFinder.h"

ClassImp(AliJFJetFinder)

AliJFJetFinder::AliJFJetFinder(Int_t n) : fNJets(0),fNJetsMax(n),fJets(n),
				    fPtMin(0),fPtMax(1000),
                                    fEtaMin(-1),fEtaMax(1),
				    fPhiMin(0),fPhiMax(6.4),
				    fNeutral(kTRUE),fCharged(kTRUE),fEM(kTRUE)
{
}

AliJFJetFinder::~AliJFJetFinder()
{
}

Bool_t AliJFJetFinder::IsAcceptedParticle(TParticle *p)
{
  if(p->GetStatusCode()%100!=1) return kFALSE;

  return kTRUE;
}

void AliJFJetFinder::SetPtCut(Float_t ptmin, Float_t ptmax)
{
  fPtMin=ptmin;
  fPtMax=ptmax;
}

void AliJFJetFinder::SetPhiCut(Float_t phimin, Float_t phimax)
{
  fPhiMin=phimin;
  fPhiMax=phimax;
}

void AliJFJetFinder::SetEtaCut(Float_t emin, Float_t emax)
{
  fEtaMin=emin;
  fEtaMax=emax;
}
