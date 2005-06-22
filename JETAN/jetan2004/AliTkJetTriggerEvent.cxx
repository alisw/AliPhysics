// $Id$

#include "AliTkJetTriggerEvent.h"

AliTkJetTriggerEvent::AliTkJetTriggerEvent() : TObject() {
  decisions = new TClonesArray("AliTkJetTriggerDecision",10000);
  counter = 0;
}

AliTkJetTriggerEvent::AliTkJetTriggerEvent(AliTkJetTriggerEvent &t) : TObject() {
  counter = t.counter;
  decisions = new TClonesArray("AliTkJetTriggerDecision",counter);
  TIterator *iter = t.decisions->MakeIterator();
  UInt_t pos = 0;
  AliTkJetTriggerDecision *dec;
  while ((dec = (AliTkJetTriggerDecision *) iter->Next()) != NULL) {
    new ((*decisions)[pos++]) AliTkJetTriggerDecision(*dec);
  }
}

void AliTkJetTriggerEvent::addDecision(AliTkJetTriggerDecision *d) {
  // let's check if we have already a decision for this trigger configuration
  TIterator *iter = decisions->MakeIterator();
  AliTkJetTriggerDecision *myDec;
  while ((myDec = (AliTkJetTriggerDecision *) iter->Next()) != NULL) {
    if ((d->getNParticles() == myDec->getNParticles()) &&
	(d->getConeRadius() == myDec->getConeRadius()) &&
	(d->getPtThreshold() == myDec->getPtThreshold())) {
      delete iter;
      return;
    }
  }
  // we dont have a decision, let's add it...
  new ((*decisions)[counter++]) AliTkJetTriggerDecision(*d);
}

void AliTkJetTriggerEvent::clear() {
  decisions->Delete();
  counter = 0;
}

Bool_t AliTkJetTriggerEvent::isEventTriggered(Float_t coneRadius, 
					   UInt_t nParticles,
					   Float_t ptThreshold) {
  AliTkJetTriggerDecision *dec;
  TIterator *iter = decisions->MakeIterator();
  while ((dec = (AliTkJetTriggerDecision *) iter->Next()) != NULL) {
    if ((dec->getConeRadius() == coneRadius) &&
	(dec->getNParticles() == nParticles) &&
	(dec->getPtThreshold() == ptThreshold)) {
      delete iter;
      return kTRUE;
    }
  }
  delete iter;
  return kFALSE;
}

ClassImp(AliTkJetTriggerEvent)
