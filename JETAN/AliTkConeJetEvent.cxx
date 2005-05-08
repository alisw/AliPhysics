//$Id$

#include <Riostream.h>
#include <TROOT.h>
#include <TClonesArray.h>
#include <TParticle.h>
#include <TMath.h>

#include "AliTkTowerV2.h"
#include "AliTkConeJet.h"
#include "AliTkConeJetEvent.h"

#define Njets__ 100

ClassImp(AliTkConeJetEvent)

AliTkConeJetEvent::AliTkConeJetEvent() : TObject(),
				   fNJets(0),fJets(0),
#ifdef ALICEINTERFACE
				   fParticles(0),
#endif
				   fDesc(""),fRadius(0.),fPtCut(0.),fEtCut(0.)
{
  fJets = new TClonesArray("AliTkConeJet",Njets__);
}

AliTkConeJetEvent::~AliTkConeJetEvent() 
{
  delete fJets;
#ifdef ALICEINTERFACE
  if(fParticles) delete fParticles;
#endif
}

void AliTkConeJetEvent::Clear(Option_t *option)
{
  TObject::Clear(option);
  delete fJets;
  fJets = new TClonesArray("AliTkConeJet",Njets__);
  //fJets->Clear("C"); // I really dont understand this!!!
  fNJets=0;
#ifdef ALICEINTERFACE
  if(fParticles) delete fParticles;
  fParticles=0;
#endif
}

void AliTkConeJetEvent::addJet(AliTkConeJet *jet) 
{
  if (!jet) {
    return;
  }
  new ((*fJets)[fNJets]) AliTkConeJet(*jet);
  fNJets++;
}
