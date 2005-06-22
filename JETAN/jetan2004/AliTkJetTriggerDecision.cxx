// $Id$

#include "AliTkJetTriggerDecision.h"

AliTkJetTriggerDecision::AliTkJetTriggerDecision(AliTkJetTriggerDecision &d):TObject() {
  setConeRadius(d.getConeRadius());
  setPtThreshold(d.getPtThreshold());
  setNParticles(d.getNParticles());
}

ClassImp(AliTkJetTriggerDecision)
