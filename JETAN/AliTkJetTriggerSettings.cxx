// $Id$

#include "AliTkJetTriggerSettings.h"

AliTkJetTriggerSettings::AliTkJetTriggerSettings() : TObject() {
  mConeRadia = new TArrayF(10);
  mNParticles = new TArrayI(10);
  mPtThresholds = new TArrayF(10);
}

AliTkJetTriggerSettings::AliTkJetTriggerSettings(AliTkJetTriggerSettings &s) :TObject(s) {
  mConeRadia = new TArrayF(*(s.getConeRadia()));
  mNParticles = new TArrayI(*(s.getNParticles()));
  mPtThresholds = new TArrayF(*(s.getPtThresholds()));
}

AliTkJetTriggerSettings::~AliTkJetTriggerSettings() {
  delete mConeRadia;
  delete mNParticles;
  delete mPtThresholds;
}

ClassImp(AliTkJetTriggerSettings)
