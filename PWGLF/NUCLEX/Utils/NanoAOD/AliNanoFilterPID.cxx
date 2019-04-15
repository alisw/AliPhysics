#include "AliNanoFilterPID.h"

#include <TObject.h>

#include <AliAnalysisManager.h>
#include <AliESDtrackCuts.h>
#include <AliInputEventHandler.h>
#include <AliPIDResponse.h>

#include <AliVEvent.h>

AliNanoFilterPID::AliNanoFilterPID()
    : AliAnalysisCuts{}, fTrackCuts{}, fTOFpidTriggerNsigma{},
      fTOFpidTriggerPtRange{}, fTPCpidTriggerNsigma{}, fTPCpidTriggerPtRange{} {
  for (int iSpecies{0}; iSpecies < AliPID::kSPECIESC; ++iSpecies) {
    fTrackCuts[iSpecies] = nullptr;
    fTOFpidTriggerNsigma[iSpecies] = -1.;
    fTOFpidTriggerPtRange[iSpecies][0] = 0.;
    fTOFpidTriggerPtRange[iSpecies][1] = 1000.;
    fTPCpidTriggerNsigma[iSpecies] = -1.;
    fTPCpidTriggerPtRange[iSpecies][0] = 0.;
    fTPCpidTriggerPtRange[iSpecies][1] = 1000.;
  }
}

void AliNanoFilterPID::TriggerOnSpecies(AliPID::EParticleType sp,
                                        AliESDtrackCuts *cuts, double nsigmaTPC,
                                        double ptRangeTPC[2], double nsigmaTOF,
                                        double ptRangeTOF[2]) {
  fTrackCuts[sp] = cuts;
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
  AliVEvent *ev = dynamic_cast<AliVEvent *>(obj);
  if (!ev) {
    return fSelected;
  }
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *handl =
      (AliInputEventHandler *)mgr->GetInputEventHandler();
  AliPIDResponse *pid = handl->GetPIDResponse();
  
  for (int iTrack{0}; iTrack < ev->GetNumberOfTracks(); ++iTrack) {
    AliVTrack* trk = (AliVTrack*)ev->GetTrack(iTrack);
    for (int iSpecies{0}; iSpecies < AliPID::kSPECIESC; ++iSpecies) {
      if (!fTrackCuts[iSpecies] && fTPCpidTriggerNsigma[iSpecies] < 0. &&
          fTOFpidTriggerNsigma[iSpecies] < 0.)
        continue;

      if (!fTrackCuts[iSpecies]->AcceptVTrack(trk))
        continue;

      double pt = trk->Pt() * AliPID::ParticleCharge(iSpecies);
      double nTPCsigma =
          std::abs(pid->NumberOfSigmasTPC(trk, AliPID::EParticleType(iSpecies)));
      double nTOFsigma =
          std::abs(pid->NumberOfSigmasTOF(trk, AliPID::EParticleType(iSpecies)));
      if ((nTPCsigma < fTPCpidTriggerNsigma[iSpecies] &&
          pt > fTPCpidTriggerPtRange[iSpecies][0] &&
          pt < fTPCpidTriggerPtRange[iSpecies][1]) ||
          (nTOFsigma < fTOFpidTriggerNsigma[iSpecies] &&
          nTPCsigma < fTPCpidTriggerNsigma[iSpecies] &&
          pt > fTOFpidTriggerPtRange[iSpecies][0] &&
          pt < fTOFpidTriggerPtRange[iSpecies][1])) {
        SETBIT(fFilterMask,iSpecies);
        fSelected = true;
      }
    }
  }

  return fSelected;
}

bool AliNanoFilterPID::IsSelected(TList *) {
  Fatal("AliNanoFilterPID::IsSelected","Method not implemented for lists");
  return false;
}
