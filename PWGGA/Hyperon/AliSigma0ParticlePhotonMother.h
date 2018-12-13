#ifndef AliSigma0ParticlePhotonMother_H
#define AliSigma0ParticlePhotonMother_H

#include "AliAODConversionPhoton.h"
#include "AliSigma0ParticleBase.h"
#include "AliSigma0ParticleV0.h"
#include "Riostream.h"
#include "TObject.h"

class AliSigma0ParticlePhotonMother : public AliSigma0ParticleBase {
 public:
  AliSigma0ParticlePhotonMother();
  virtual ~AliSigma0ParticlePhotonMother() {}
  AliSigma0ParticlePhotonMother(const AliSigma0ParticleV0 &lambdaCandidate,
                                const AliSigma0ParticleV0 &photonCandidate,
                                const AliVEvent *inputEvent);
  AliSigma0ParticlePhotonMother &operator=(
      const AliSigma0ParticlePhotonMother &obj);

  double GetPx() const { return AliSigma0ParticleBase::GetPx(); }
  double GetPy() const { return AliSigma0ParticleBase::GetPy(); }
  double GetPz() const { return AliSigma0ParticleBase::GetPz(); }

  double GetPxMC() const { return AliSigma0ParticleBase::GetPxMC(); }
  double GetPyMC() const { return AliSigma0ParticleBase::GetPyMC(); }
  double GetPzMC() const { return AliSigma0ParticleBase::GetPzMC(); }

  int GetPDGcode() const { return AliSigma0ParticleBase::GetPDGcode(); }
  double GetMass() const { return AliSigma0ParticleBase::GetMass(); }

  double GetPt() const { return AliSigma0ParticleBase::GetPt(); }
  double GetP() const { return AliSigma0ParticleBase::GetP(); }
  int GetTrackLabel() const { return AliSigma0ParticleBase::GetTrackLabel(); }
  double GetPhi() const { return AliSigma0ParticleBase::GetPhi(); }
  double GetEta() const { return AliSigma0ParticleBase::GetEta(); }
  bool GetIsUse() const { return AliSigma0ParticleBase::GetIsUse(); }

  double GetRecMassPhoton() const { return fRecMassPhoton; }
  double GetRecMassLambda() const { return fRecMassLambda; }
  double GetRecMass() const { return fRecMass; }
  double GetRapidity() const;

  float GetArmenterosAlpha() const;
  float GetArmenterosQt() const;

  int GetMCLabel() const { return fMCLabel; }
  int GetPDGCode() const { return fPDGCode; }

  AliSigma0ParticleV0 GetV0() const { return fV0; }
  AliSigma0ParticleV0 GetPhoton() const { return fPhoton; }

  int MatchToMC(const AliMCEvent *mcEvent, const int PIDmother,
                const std::vector<int> PIDdaughters, int &pdgLambdaMother,
                int &pdgPhotonMother);

 private:
  short fType;
  double fRecMassPhoton;
  double fRecMassLambda;
  double fRecMass;
  int fMCLabel;
  int fPDGCode;
  AliSigma0ParticleV0 fV0;
  AliSigma0ParticleV0 fPhoton;

  ClassDef(AliSigma0ParticlePhotonMother, 4)
};

#endif
