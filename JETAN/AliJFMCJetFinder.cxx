// $Id$

#include <Riostream.h>

#include <TCollection.h>
#include <TParticle.h>
#include <TIterator.h>

#include "AliJFMCJet.h"
#include "AliJFMCJetFinder.h"

ClassImp(AliJFMCJetFinder)

AliJFMCJetFinder::AliJFMCJetFinder(Int_t n) : AliJFJetFinder(n), fParticles(NULL), fJet(NULL)
{
}

AliJFMCJetFinder::~AliJFMCJetFinder()
{
  Clean();
}

Int_t AliJFMCJetFinder::Init(TClonesArray *particles)
{
  if(particles==NULL) return 0;

  fParticles=particles;
  return fParticles->GetEntries();
}

Int_t AliJFMCJetFinder::Run()
{
  if(fParticles==NULL) return 0;

  TIterator *iter=fParticles->MakeIterator();
  TParticle *p;

  while((p=(TParticle*)iter->Next()) != NULL){
    if(p->GetPdgCode()==92){
      fJets.AddAt(new AliJFMCJet(),fNJets);
      fJet=(AliJFMCJet*)fJets.At(fNJets);
      fJet->SetNJet(++fNJets);
      FollowDaughters(p->GetFirstDaughter(),p->GetLastDaughter()); 
      //cout << "Add new Jet " << fNJets << "with " << fJet->GetE() << " " << fJet->GetE_() << " " << fJet->GetNPart() << endl;
      fJet->Update();
    }
  }    

  fJets.Sort();
  //fJets.Print();
  
  return fNJets;
}

void AliJFMCJetFinder::Debug()
{
  TIterator *iter=fJets.MakeIterator();
  AliJFMCJet *j;

  while((j=(AliJFMCJet*)iter->Next()) != NULL){
    j->Update();
    //cout << j << ": " << j->GetNPart() << " " << j->GetNJet() << " " << j->GetPtMax() << endl;
    j->Print("");
  }

  /*
  for(Int_t i=0;i<fNJets;i++){
    JFJet *j=(JFMCJet*)fJets->At(i);
    cout << i << ": " << j->GetNPart()<< " " << j->GetNJet() << endl;
  }
  */
}

void AliJFMCJetFinder::FollowDaughters(Int_t first, Int_t last)
{
  if(last<first) return;

  for(Int_t i=first-1;i<=last-1;i++){
    //Fortran (Pythia/Hijing) has different numbering than TParticle

    TParticle *p = (TParticle *)fParticles->At(i);
    if(p->GetStatusCode()!=1){
       FollowDaughters(p->GetFirstDaughter(),p->GetLastDaughter());
    }
    else {
      if(!IsAcceptedParticle(p)) continue;

      fJet->AddParticle(p);
      //p->Print();
    } 
  }
}

/*
  void AliJFMCJetFinder::Clean()
  {
  }
*/
