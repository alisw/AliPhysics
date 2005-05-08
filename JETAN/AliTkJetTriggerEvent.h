// $Id$

#ifndef ALITKJETTRIGGEREVENT_H
#define ALITKJETTRIGGEREVENT_H

#include <TROOT.h>
#include <TObjArray.h>
#include <TClonesArray.h>
#include "AliTkJetTriggerDecision.h"

class AliTkJetTriggerEvent : public TObject {
 public:
  AliTkJetTriggerEvent();
  AliTkJetTriggerEvent(AliTkJetTriggerEvent &t);

  void addDecision(AliTkJetTriggerDecision *d);
  void clear();
  Bool_t isEventTriggered(Float_t coneRadius, UInt_t nParticles, 
			Float_t ptThreshold);

  TClonesArray *decisions; //->
  UInt_t counter;

 private:

  ClassDef(AliTkJetTriggerEvent,1)
};
#endif
