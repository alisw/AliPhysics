#ifndef ALINANOFILTERPID_H
#define ALINANOFILTERPID_H

#include <Rtypes.h>

#include <AliAnalysisCuts.h>
#include <AliPID.h>

class TObject;
class TList;

class AliESDtrackCuts;

class AliNanoFilterPID : public AliAnalysisCuts {
public:
  AliNanoFilterPID();

  virtual bool IsSelected(TObject *obj);
  virtual bool IsSelected(TList *);

  void TriggerOnSpecies(AliPID::EParticleType sp, AliESDtrackCuts *cuts, ULong_t fb,
                        double nsigmaTPC, double ptRangeTPC[2],
                        double nsigmaTOF, double ptRangeTOF[2]);

  float fMinDeltaM;
  float fMaxDeltaM;

private:
  bool fTriggerOnSpecies[AliPID::kSPECIESC];
  AliESDtrackCuts *fTrackCuts[AliPID::kSPECIESC];
  ULong_t fFilterBits[AliPID::kSPECIESC];
  double fTOFpidTriggerNsigma[AliPID::kSPECIESC];
  double fTOFpidTriggerPtRange[AliPID::kSPECIESC][2];
  double fTPCpidTriggerNsigma[AliPID::kSPECIESC];
  double fTPCpidTriggerPtRange[AliPID::kSPECIESC][2];
  ClassDef(AliNanoFilterPID, 1)
};

#endif