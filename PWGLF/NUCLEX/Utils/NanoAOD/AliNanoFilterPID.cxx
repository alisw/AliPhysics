#include "AliNanoFilterPID.h"

#include <TObject.h>

#include <AliAnalysisManager.h>
#include <AliAODTrack.h>
#include <AliESDtrackCuts.h>
#include <AliInputEventHandler.h>
#include <AliPIDResponse.h>
#include <AliAnalysisTaskNucleiYield.h>

#include <AliVEvent.h>

AliNanoFilterPID::AliNanoFilterPID()
    : AliAnalysisCuts{}, fMinDeltaM{-2.4}, fMaxDeltaM{3.6 }, fTriggerOnSpecies{}, fTrackCuts{}, fFilterBits{}, fTOFpidTriggerNsigma{},
      fTOFpidTriggerPtRange{}, fTPCpidTriggerNsigma{}, fTPCpidTriggerPtRange{} {
  for (int iSpecies{0}; iSpecies < AliPID::kSPECIESC; ++iSpecies) {
    fTriggerOnSpecies[iSpecies] = false;
    fTrackCuts[iSpecies] = nullptr;
    fFilterBits[iSpecies] = 0ul;
    fTOFpidTriggerNsigma[iSpecies] = -1.;
    fTOFpidTriggerPtRange[iSpecies][0] = 0.;
    fTOFpidTriggerPtRange[iSpecies][1] = 1000.;
    fTPCpidTriggerNsigma[iSpecies] = -1.;
    fTPCpidTriggerPtRange[iSpecies][0] = 0.;
    fTPCpidTriggerPtRange[iSpecies][1] = 1000.;
  }
}

void AliNanoFilterPID::TriggerOnSpecies(AliPID::EParticleType sp,
                                        AliESDtrackCuts *cuts, ULong_t fb,
                                        double nsigmaTPC, double ptRangeTPC[2],
                                        double nsigmaTOF, double ptRangeTOF[2]) {
  fTriggerOnSpecies[sp] = true;
  fTrackCuts[sp] = cuts;
  fFilterBits[sp] = fb;
  fTOFpidTriggerNsigma[sp] = nsigmaTOF;
  fTOFpidTriggerPtRange[sp][0] = ptRangeTOF[0];
  fTOFpidTriggerPtRange[sp][1] = ptRangeTOF[1];
  fTPCpidTriggerNsigma[sp] = nsigmaTPC;
  fTPCpidTriggerPtRange[sp][0] = ptRangeTPC[0];
  fTPCpidTriggerPtRange[sp][1] = ptRangeTPC[1];
}

bool AliNanoFilterPID::IsSelected(TObject *obj) {
  fSelected = false;
  fFilterMask = 0u;
  AliVTrack *trk = dynamic_cast<AliVTrack *>(obj);
  if (!trk) {
    return fSelected;
  }
  AliAODTrack* aodTrk = dynamic_cast<AliAODTrack *>(obj);

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *handl =
      (AliInputEventHandler *)mgr->GetInputEventHandler();
  AliPIDResponse *pid = handl->GetPIDResponse();
  
  for (int iSpecies{0}; iSpecies < AliPID::kSPECIESC; ++iSpecies) {
    if (!fTriggerOnSpecies[iSpecies])
      continue;

    if (fTrackCuts[iSpecies])
      if (!fTrackCuts[iSpecies]->AcceptVTrack(trk))
        continue;

    if (aodTrk && fFilterBits[iSpecies])
      if (!aodTrk->TestFilterBit(fFilterBits[iSpecies]))
        continue;

    double pt = trk->Pt() * AliPID::ParticleCharge(iSpecies);
    double nTPCsigma =
        std::abs(pid->NumberOfSigmasTPC(trk, AliPID::EParticleType(iSpecies)));
    double nTOFsigma =
        std::abs(pid->NumberOfSigmasTOF(trk, AliPID::EParticleType(iSpecies)));
    
    float beta = AliAnalysisTaskNucleiYield::HasTOF(aodTrk,pid);
    if (beta > 1. - EPS) beta = -1;
    const float m2 = aodTrk->P() * aodTrk->P() * (1.f / (beta * beta) - 1.f);
    bool goodTOF = m2 > fMinDeltaM && m2 < fMaxDeltaM;

    if ((nTPCsigma < fTPCpidTriggerNsigma[iSpecies] &&
        pt > fTPCpidTriggerPtRange[iSpecies][0] &&
        pt < fTPCpidTriggerPtRange[iSpecies][1]) ||
        (goodTOF &&
        nTOFsigma < fTOFpidTriggerNsigma[iSpecies] &&
        pt > fTOFpidTriggerPtRange[iSpecies][0] &&
        pt < fTOFpidTriggerPtRange[iSpecies][1])) {
      SETBIT(fFilterMask,iSpecies);
      fSelected = true;
    }
  }

  return fSelected;
}

bool AliNanoFilterPID::IsSelected(TList *) {
  Fatal("AliNanoFilterPID::IsSelected","Method not implemented for lists");
  return false;
}
