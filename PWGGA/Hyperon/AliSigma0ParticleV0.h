#ifndef AliSigma0ParticleV0_H
#define AliSigma0ParticleV0_H

#include "AliAODConversionPhoton.h"
#include "AliESDVertex.h"
#include "AliESDv0.h"
#include "AliSigma0ParticleBase.h"
#include "Riostream.h"
#include "TObject.h"

class AliSigma0ParticleV0 : public AliSigma0ParticleBase {
 public:
  AliSigma0ParticleV0();
  virtual ~AliSigma0ParticleV0() {}
  AliSigma0ParticleV0(AliESDv0 *v0, const AliESDtrack *pos,
                      const AliESDtrack *neg, const AliESDVertex *vertex,
                      const int pdg, const float magneticField,
                      AliMCEvent *mcEvent);
  AliSigma0ParticleV0(const AliAODConversionPhoton &gamma,
                      const AliVEvent *inputEvent);
  AliSigma0ParticleV0 &operator=(const AliSigma0ParticleV0 &obj);

  double GetPx() const { return AliSigma0ParticleBase::GetPx(); }
  double GetPy() const { return AliSigma0ParticleBase::GetPy(); }
  double GetPz() const { return AliSigma0ParticleBase::GetPz(); }

  double GetPxMC() const { return AliSigma0ParticleBase::GetPxMC(); }
  double GetPyMC() const { return AliSigma0ParticleBase::GetPyMC(); }
  double GetPzMC() const { return AliSigma0ParticleBase::GetPzMC(); }

  int GetPDGcode() const { return AliSigma0ParticleBase::GetPDGcode(); }
  int GetPDGcodeMother() const {
    return AliSigma0ParticleBase::GetPDGcodeMother();
  }
  double GetMass() const { return AliSigma0ParticleBase::GetMass(); }

  double GetPt() const { return AliSigma0ParticleBase::GetPt(); }
  int GetTrackLabel() const { return AliSigma0ParticleBase::GetTrackLabel(); }
  double GetPhi() const { return AliSigma0ParticleBase::GetPhi(); }
  double GetEta() const { return AliSigma0ParticleBase::GetEta(); }
  bool GetIsUse() const { return AliSigma0ParticleBase::GetIsUse(); }

  int GetTrackLabelPos() const { return fTrackLabelPos; }
  int GetTrackLabelNeg() const { return fTrackLabelNeg; }
  int GetMCLabelPos() const { return fTrackPos.GetMCLabel(); }
  int GetMCLabelNeg() const { return fTrackNeg.GetMCLabel(); }
  double GetCosineAlpha() const { return fCosAlpha; }
  double GetRecMass() const { return fRecMass; }
  double GetPDGMass() const { return fPDGMass; }

  void SetPDGMass(float mass) { fPDGMass = mass; }

  float GetArmenterosAlpha() const;
  float GetArmenterosQt() const;

  AliSigma0ParticleBase GetNegDaughter() const { return fTrackNeg; }
  AliSigma0ParticleBase GetPosDaughter() const { return fTrackPos; }

  int MatchToMC(AliMCEvent *mcEvent, const int PIDmother,
                const std::vector<int> PIDdaughters);

 private:
  int fTrackLabelPos;
  int fTrackLabelNeg;
  int fMCLabelPos;
  int fMCLabelNeg;
  AliSigma0ParticleBase fTrackPos;
  AliSigma0ParticleBase fTrackNeg;
  double fCosAlpha;
  double fRecMass;
  double fPDGMass;

  ClassDef(AliSigma0ParticleV0, 1)
};

#endif
