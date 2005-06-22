// $Id$

#ifndef ALITKJETTRIGGERSETTINGS_H
#define ALITKJETTRIGGERSETTINGS_H

#include <TROOT.h>
#include <TArrayF.h>
#include <TArrayI.h>

class AliTkJetTriggerSettings : public TObject {
 public:
  AliTkJetTriggerSettings();
  AliTkJetTriggerSettings(AliTkJetTriggerSettings &s);
  ~AliTkJetTriggerSettings();

  TArrayF *getConeRadia() { return mConeRadia; }
  TArrayI *getNParticles() { return mNParticles; }
  TArrayF *getPtThresholds() { return mPtThresholds; }

 private:
  TArrayF *mConeRadia; //->
  TArrayI *mNParticles; //->
  TArrayF *mPtThresholds; //->


  ClassDef(AliTkJetTriggerSettings,1)
};
#endif
