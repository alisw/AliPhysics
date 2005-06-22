// $Id$

#ifndef ALIJFJETTRIGGERH
#define ALIJFJETTRIGGERH

#include <TObject.h>
#include "AliJFJetFinder.h"

class AliJFJetTrigger : public AliJFJetFinder
{
 public:
  AliJFJetTrigger(Int_t n=50);
  virtual ~AliJFJetTrigger();

  Bool_t IsAcceptedParticle(TParticle *p);

  void Clean(){fParticles->Delete();}

 protected:

  TClonesArray *fParticles;

  ClassDef(AliJFJetTrigger,1) //AliJFJetTrigger class
};

#endif /*ALIJFJETTRIGGERH*/
