// $Id$

#ifndef ALITKJETTRIGGERDECISION_H
#define ALITKJETTRIGGERDECISION_H

#include <TROOT.h>

class AliTkJetTriggerDecision : public TObject {
 public:
  AliTkJetTriggerDecision() : TObject(), mConeRadius(-999), mPtThreshold(-999),
    mNParticles(0) { }
  AliTkJetTriggerDecision(AliTkJetTriggerDecision &d);

  void setConeRadius(Float_t r) { mConeRadius = r; }
  Float_t getConeRadius() { return mConeRadius; }
  void setPtThreshold(Float_t pt) { mPtThreshold = pt; }
  Float_t getPtThreshold() { return mPtThreshold; }
  void setNParticles(UInt_t n) { mNParticles = n; }
  UInt_t getNParticles() { return mNParticles; }

 private:
  Float_t mConeRadius;
  Float_t mPtThreshold;
  UInt_t mNParticles;

  ClassDef(AliTkJetTriggerDecision,1)
};
#endif
