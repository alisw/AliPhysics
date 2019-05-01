/**************************************************************************
* Copyright(c) 1998-2017, ALICE Experiment at CERN, All rights reserved.  *
*                                                                         *
* Permission to use, copy, modify and distribute this software and its    *
* documentation strictly for non-commercial purposes is hereby granted    *
* without fee, provided that the above copyright notice appears in all    *
* copies and that both the copyright notice and this permission notice    *
* appear in the supporting documentation. The authors make no claims      *
* about the suitability of this software for any purpose. It is           *
* provided "as is" without express or implied warranty.                   *
**************************************************************************/

#include <TH1F.h>
#include <TH2F.h>
#include "AliVParticle.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliMCEventHandler.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliAODMCParticle.h"
#include "AliAODTrack.h"
#include "AliPIDResponse.h"
#include "AliLog.h"
#include "AliCSTrackMaps.h"
#include "AliCSPIDCuts.h"

/// \file AliCSPIDCuts.cxx
/// \brief Implementation of PID track cuts class within the correlation studies analysis

/// \cond CLASSIMP
ClassImp(AliCSPIDCuts);
/// \endcond

const char *AliCSPIDCuts::fgkCutsNames[AliCSPIDCuts::kNCuts] = {
    "ITS dE/dx n#sigma",
    "TPC dE/dx n#sigma",
    "TOF n#sigma"
};


/// Default constructor for serialization
AliCSPIDCuts::AliCSPIDCuts() :
    AliCSTrackCutsBase(),
    fMinP(0.0),
    fMaxP(9999.0),
//    fITSnSigmaAbove{100.0},
//    fITSnSigmaBelow{-100.0},
//    fTPCnSigmaAbove{100.0},
//    fTPCnSigmaBelow{-100.0},
//    fTOFRequired{kFALSE},
//    fTOFnSigmaAbove{100.0},
//    fTOFnSigmaBelow{-100.0},
    fITSEnabledSpeciesMask(TBits()),
    fTPCEnabledSpeciesMask(TBits()),
    fTOFEnabledSpeciesMask(TBits()),
    fPIDResponse(NULL),
    fTargetSpecies(AliPID::kUnknown),
    fhCutsStatistics(NULL),
    fhCutsCorrelation(NULL)
//    fhITSdEdxSigmaVsP{NULL},
//    fhITSdEdxSignalVsP{NULL},
//    fhTPCdEdxSigmaVsP{NULL},
//    fhTPCdEdxSignalVsP{NULL},
//    fhTOFSigmaVsP{NULL},
//    fhTOFSignalVsP{NULL}
{
}

/// Constructor
/// \param name name of the event cuts
/// \param title title of the event cuts
/// \param target the PID target particle of type AliPID::EParticleType
AliCSPIDCuts::AliCSPIDCuts(const char *name, const char *title, AliPID::EParticleType target) :
    AliCSTrackCutsBase(kNCuts,kNCutsParameters,name,title),
    fMinP(0.0),
    fMaxP(9999.0),
    //    fITSnSigmaAbove{100.0},
    //    fITSnSigmaBelow{-100.0},
    //    fTPCnSigmaAbove{100.0},
    //    fTPCnSigmaBelow{-100.0},
    //    fTOFRequired{kFALSE},
    //    fTOFnSigmaAbove{100.0},
    //    fTOFnSigmaBelow{-100.0},
    fITSEnabledSpeciesMask(TBits(AliPID::kSPECIESC)),
    fTPCEnabledSpeciesMask(TBits(AliPID::kSPECIESC)),
    fTOFEnabledSpeciesMask(TBits(AliPID::kSPECIESC)),
    fPIDResponse(NULL),
    fTargetSpecies(target),
    fhCutsStatistics(NULL),
    fhCutsCorrelation(NULL)
//    fhITSdEdxSigmaVsP{NULL},
//    fhITSdEdxSignalVsP{NULL},
//    fhTPCdEdxSigmaVsP{NULL},
//    fhTPCdEdxSignalVsP{NULL},
//    fhTOFSigmaVsP{NULL},
//    fhTOFSignalVsP{NULL}
{
  for (Int_t spid = AliPID::kElectron; spid < AliPID::kSPECIESC; spid++) {
    fITSnSigmaAbove[spid] = 100.0;
    fITSnSigmaBelow[spid] = -100.0;
    fTPCnSigmaAbove[spid] = 100.0;
    fTPCnSigmaBelow[spid] = -100.0;
    fTOFRequired[spid] = kFALSE;
    fTOFnSigmaAbove[spid] = 100.0;
    fTOFnSigmaBelow[spid] = -100.0;
  }

  for (Int_t i = 0; i < 2; i++) {
    fhITSdEdxSigmaVsP[i] = NULL;
    fhITSdEdxSignalVsP[i] = NULL;
    fhTPCdEdxSigmaVsP[i] = NULL;
    fhTPCdEdxSignalVsP[i] = NULL;
    fhTOFSigmaVsP[i] = NULL;
    fhTOFSignalVsP[i] = NULL;
  }

  fITSEnabledSpeciesMask.ResetAllBits();
  fTPCEnabledSpeciesMask.ResetAllBits();
  fTOFEnabledSpeciesMask.ResetAllBits();
}

/// Destructor
/// We don't own anything, everything we allocate is owned
/// by the output list
AliCSPIDCuts::~AliCSPIDCuts()
{

}

/// Processes a potential change in the run number
///
/// Checks if the current period under analysis has changed and if so
/// updates the needed members
void AliCSPIDCuts::NotifyRun() {

  /* checks the change in the analysis period */
  if (AliCSTrackCutsBase::GetGlobalPeriod() != fDataPeriod) {

    fDataPeriod = AliCSTrackCutsBase::GetGlobalPeriod();

    /* and now we ask for histogram allocation */
    DefineHistograms();
  }
}

/// Processes the start of a new event
///
/// Does nothin for the time being
void AliCSPIDCuts::NotifyEvent() {

}

/// Check whether the passed track is recognized as the target by the different configured PID cuts
/// \param ttrk the track to analyze whether it is recognized as the target or not
/// \return kTRUE if the track is recognized, kFALSE otherwise
///
Bool_t AliCSPIDCuts::IsTrackAccepted(AliVTrack *ttrk) {
  /* just to be sure */
  if (ttrk == NULL) return kFALSE;

  /* if not in the momentum range it is not recognized */
  if (ttrk->P() < fMinP || fMaxP < ttrk->P()) return kFALSE;

  /* for the time being */
  Bool_t accepted = kTRUE;

  /* initialize the mask of activated cuts */
  fCutsActivatedMask.ResetAllBits();

  /* we now need to consider the potential constrained track */
  AliVTrack *trk = ttrk;
  if (ttrk->GetID() < 0)
    /* let's switch to the original one which has the PID information */
    trk = AliCSTrackMaps::GetOriginalTrack(dynamic_cast<AliAODTrack*>(ttrk));


  /* ITS PID cut */
  if (fCutsEnabledMask.TestBitNumber(kITSdEdxSigmaCut)) {
    if( fPIDResponse->NumberOfSigmasITS(trk, fTargetSpecies) < fITSnSigmaBelow[fTargetSpecies] ||
        fITSnSigmaAbove[fTargetSpecies] < fPIDResponse->NumberOfSigmasITS(trk, fTargetSpecies)){
      fCutsActivatedMask.SetBitNumber(kITSdEdxSigmaCut);
      accepted = kFALSE;
    }
    else {
      /* now check the separation if required */
      if (fITSEnabledSpeciesMask.CountBits() > 1) {
        AliPID::EParticleType spid = AliPID::EParticleType(fITSEnabledSpeciesMask.FirstSetBit());
        while (spid < AliPID::kDeuteron) {
          if (spid != fTargetSpecies) {
            if(fITSnSigmaBelow[spid] < fPIDResponse->NumberOfSigmasITS(trk, spid) &&
                fPIDResponse->NumberOfSigmasITS(trk, spid) < fITSnSigmaAbove[spid]){
              fCutsActivatedMask.SetBitNumber(kITSdEdxSigmaCut);
              accepted = kFALSE;
              break;
            }
          }
          spid = AliPID::EParticleType(fITSEnabledSpeciesMask.FirstSetBit(spid+1));
        }
      }
    }
  }

  /* TPC PID cut */
  if (fCutsEnabledMask.TestBitNumber(kTPCdEdxSigmaCut)) {
    if( fPIDResponse->NumberOfSigmasTPC(trk, fTargetSpecies) < fTPCnSigmaBelow[fTargetSpecies] ||
        fTPCnSigmaAbove[fTargetSpecies] < fPIDResponse->NumberOfSigmasTPC(trk, fTargetSpecies)){
      fCutsActivatedMask.SetBitNumber(kTPCdEdxSigmaCut);
      accepted = kFALSE;
    }
    else {
      /* now check the separation if required */
      if (fTPCEnabledSpeciesMask.CountBits() > 1) {
        AliPID::EParticleType spid = AliPID::EParticleType(fTPCEnabledSpeciesMask.FirstSetBit());
        while (spid < AliPID::kDeuteron) {
          if (spid != fTargetSpecies) {
            if(fTPCnSigmaBelow[spid] < fPIDResponse->NumberOfSigmasTPC(trk, spid) &&
                fPIDResponse->NumberOfSigmasTPC(trk, spid) < fTPCnSigmaAbove[spid]){
              fCutsActivatedMask.SetBitNumber(kTPCdEdxSigmaCut);
              accepted = kFALSE;
              break;
            }
          }
          spid = AliPID::EParticleType(fTPCEnabledSpeciesMask.FirstSetBit(spid+1));
        }
      }
    }
  }

  /* TOF PID cut */
  if (fCutsEnabledMask.TestBitNumber(kTOFSigmaCut)) {
    if ((trk->GetStatus() & AliESDtrack::kTOFin) && (!(trk->GetStatus() & AliESDtrack::kTOFmismatch))) {
      if( fPIDResponse->NumberOfSigmasTOF(trk, fTargetSpecies) < fTOFnSigmaBelow[fTargetSpecies] ||
          fTOFnSigmaAbove[fTargetSpecies] < fPIDResponse->NumberOfSigmasTOF(trk, fTargetSpecies)){
        fCutsActivatedMask.SetBitNumber(kTOFSigmaCut);
        accepted = kFALSE;
      }
      else {
        /* now check the separation if required */
        if (fTOFEnabledSpeciesMask.CountBits() > 1) {
          AliPID::EParticleType spid = AliPID::EParticleType(fTOFEnabledSpeciesMask.FirstSetBit());
          while (spid < AliPID::kDeuteron) {
            if (spid != fTargetSpecies) {
              if(fTOFnSigmaBelow[spid] < fPIDResponse->NumberOfSigmasTOF(trk, spid) &&
                  fPIDResponse->NumberOfSigmasTOF(trk, spid) < fTOFnSigmaAbove[spid]){
                fCutsActivatedMask.SetBitNumber(kTOFSigmaCut);
                accepted = kFALSE;
                break;
              }
            }
            spid = AliPID::EParticleType(fTOFEnabledSpeciesMask.FirstSetBit(spid+1));
          }
        }
      }
    }
    else {
      if (fTOFRequired[fTargetSpecies]) {
        fCutsActivatedMask.SetBitNumber(kTOFSigmaCut);
        accepted = kFALSE;
      }
    }
  }

  if (fQALevel > kQALevelNone) {
    /* let's fill the histograms */
    fhCutsStatistics->Fill(fhCutsStatistics->GetBinCenter(fhCutsStatistics->GetXaxis()->FindBin("n tracks")));
    if (!accepted)
      fhCutsStatistics->Fill(fhCutsStatistics->GetBinCenter(fhCutsStatistics->GetXaxis()->FindBin("n cut tracks")));

    for (Int_t i=0; i<kNCuts; i++) {
      if (fhCutsStatistics->GetXaxis()->FindBin(fgkCutsNames[i]) < 1)
        AliFatal(Form("Inconsistency! Cut %d with name %s not found", i, fgkCutsNames[i]));

      if (fCutsActivatedMask.TestBitNumber(i))
        fhCutsStatistics->Fill(fhCutsStatistics->GetBinCenter(fhCutsStatistics->GetXaxis()->FindBin(fgkCutsNames[i])));

      if (fQALevel > kQALevelLight) {
        for (Int_t j=i; j<kNCuts; j++) {
          if (fhCutsStatistics->GetXaxis()->FindBin(fgkCutsNames[j]) < 1)
            AliFatal(Form("Inconsistency! Cut %d with name %s not found", j, fgkCutsNames[j]));

          if (fCutsActivatedMask.TestBitNumber(i) && fCutsActivatedMask.TestBitNumber(j)) {
            Float_t xC = fhCutsCorrelation->GetXaxis()->GetBinCenter(fhCutsCorrelation->GetXaxis()->FindBin(fgkCutsNames[i]));
            Float_t yC = fhCutsCorrelation->GetYaxis()->GetBinCenter(fhCutsCorrelation->GetYaxis()->FindBin(fgkCutsNames[j]));
            fhCutsCorrelation->Fill(xC, yC);
          }
        }
      }
    }

    for (Int_t i = 0; i < 2; i++) {

      /* don't fill before if not requested */
      if (i == 1 || (fQALevel > kQALevelLight)) {
        fhITSdEdxSigmaVsP[i]->Fill(ttrk->P(),fPIDResponse->NumberOfSigmasITS(trk, fTargetSpecies));
        fhITSdEdxSignalVsP[i]->Fill(ttrk->P(),trk->GetITSsignal());
        fhTPCdEdxSigmaVsP[i]->Fill(ttrk->P(),fPIDResponse->NumberOfSigmasTPC(trk, fTargetSpecies));
        fhTPCdEdxSignalVsP[i]->Fill(ttrk->P(),TMath::Abs(trk->GetTPCsignal()));
        if ((trk->GetStatus() & AliESDtrack::kTOFin) && (!(trk->GetStatus() & AliESDtrack::kTOFmismatch))) {
          static const Double_t c_cm_ps = TMath::C() * 1.0e2 * 1.0e-12;
          Double_t tracklen_cm = trk->GetIntegratedLength();
          Double_t toftime_ps = trk->GetTOFsignal() - fPIDResponse->GetTOFResponse().GetStartTime(ttrk->P());
          Double_t beta = tracklen_cm / toftime_ps / c_cm_ps;
          fhTOFSigmaVsP[i]->Fill(ttrk->P(),fPIDResponse->NumberOfSigmasTOF(trk, fTargetSpecies));
          fhTOFSignalVsP[i]->Fill(ttrk->P(),beta);
        }
      }
      /* don't fill after if event not accepted */
      if (!accepted) break;
    }
  }
  return accepted;
}


/// Check whether the true track associated to the passed track is accepted by the PID cuts
/// \param trk the track to analyze whether its associated true track is accepted or not
/// \return kTRUE if the associated true track  is accepted, kFALSE otherwise
///
Bool_t AliCSPIDCuts::IsTrueTrackAccepted(AliVTrack *trk) {

  /* reject ghost tracks */
  if (trk->GetLabel() < 0) return kFALSE;

  return IsTrueTrackAccepted(trk->GetLabel());
}


/// Check whether the passed true track is recognized by the PID cut
/// \param itrk the index of true track to analyze whether it is recognized or not
/// \return kTRUE if the track  is recognized, kFALSE otherwise
///
Bool_t AliCSPIDCuts::IsTrueTrackAccepted(Int_t itrk) {

  Double_t p = 0.0;

  if (fgIsESD) {
    /* get the associated particle */
    AliVParticle *particle = fgMCHandler->MCEvent()->GetTrack(itrk);

    /* just to be sure */
    if (particle == NULL) return kFALSE;

    p = particle->P();

    /* if not in the momentum range it is not recognized */
    if (p < fMinP || fMaxP < p) return kFALSE;

    if (fCutsEnabledMask.TestBitNumber(kITSdEdxSigmaCut) ||
        fCutsEnabledMask.TestBitNumber(kTPCdEdxSigmaCut) ||
        fCutsEnabledMask.TestBitNumber(kTOFSigmaCut)) {
      if (GetTrueSpecies(particle) != fTargetSpecies) {
        return kFALSE;
      }
    }
  }
  else {
    /* get the associated particle */
    AliVParticle *particle = (AliVParticle *) fgMCArray->At(itrk);

    /* just to be sure */
    if (particle == NULL) return kFALSE;

    p = particle->P();

    /* if not in the momentum range it is not recognized */
    if (p < fMinP || fMaxP < p) return kFALSE;

    if (fCutsEnabledMask.TestBitNumber(kITSdEdxSigmaCut) ||
        fCutsEnabledMask.TestBitNumber(kTPCdEdxSigmaCut) ||
        fCutsEnabledMask.TestBitNumber(kTOFSigmaCut)) {
      if (GetTrueSpecies(particle) != fTargetSpecies) {
        return kFALSE;
      }
    }
  }
  return kTRUE;
}

/// Get the true species associated to a reconstructed track
/// \param trk the reconstructed track
/// \return the ID of the associated species
AliPID::EParticleType AliCSPIDCuts::GetTrueSpecies(AliVTrack *trk) {

  if (fgIsESD) {
    AliVParticle *particle = fgMCHandler->MCEvent()->GetTrack(TMath::Abs(trk->GetLabel()));
    return GetTrueSpecies(particle);
  }
  else {
    /* get the associated particle */
    AliVParticle *particle = (AliVParticle *) fgMCArray->At(TMath::Abs(trk->GetLabel()));
    return GetTrueSpecies(particle);
  }
}

/// Get the true species associated to a true particle
/// \param par the true particle
/// \return the ID of the particle species
AliPID::EParticleType AliCSPIDCuts::GetTrueSpecies(AliVParticle *par) {

  switch(par->PdgCode()) {
  case ::kPositron:
  case ::kElectron:
    return AliPID::kElectron;
    break;
  case ::kProton:
  case ::kProtonBar:
    return AliPID::kProton;
    break;
  case ::kMuonPlus:
  case ::kMuonMinus:
    return AliPID::kMuon;
    break;
  case ::kPiPlus:
  case ::kPiMinus:
    return AliPID::kPion;
    break;
  case ::kKPlus:
  case ::kKMinus:
    return AliPID::kKaon;
    break;
  default:
    return AliPID::kUnknown;
  }
}

/// Sets the individual value for the cut parameter ID
/// \param paramID the ID of the cut parameter of interest
/// \param value the value for the cut parameter
/// \return kTRUE if the cut parameter value was accepted
Bool_t AliCSPIDCuts::SetCutAndParams(Int_t paramID, Int_t value) {

  switch (cutsParametersIds(paramID)) {
  case kPRangeCutParam:
    if (SetPRange(value)) {
      fParameters[kPRangeCutParam] = value;
      UpdateCutsString();
      return kTRUE;
    } else return kFALSE;
  case kITSdEdxSigmaCutParam_e:
    if (SetITSdEdxSigmaCut(AliPID::kElectron, value)) {
      fParameters[kITSdEdxSigmaCutParam_e] = value;
      UpdateCutsString();
      return kTRUE;
    } else return kFALSE;
  case kITSdEdxSigmaCutParam_mu:
    if (SetITSdEdxSigmaCut(AliPID::kMuon, value)) {
      fParameters[kITSdEdxSigmaCutParam_mu] = value;
      UpdateCutsString();
      return kTRUE;
    } else return kFALSE;
  case kITSdEdxSigmaCutParam_pi:
    if (SetITSdEdxSigmaCut(AliPID::kPion, value)) {
      fParameters[kITSdEdxSigmaCutParam_pi] = value;
      UpdateCutsString();
      return kTRUE;
    } else return kFALSE;
  case kITSdEdxSigmaCutParam_k:
    if (SetITSdEdxSigmaCut(AliPID::kKaon, value)) {
      fParameters[kITSdEdxSigmaCutParam_k] = value;
      UpdateCutsString();
      return kTRUE;
    } else return kFALSE;
  case kITSdEdxSigmaCutParam_p:
    if (SetITSdEdxSigmaCut(AliPID::kProton, value)) {
      fParameters[kITSdEdxSigmaCutParam_p] = value;
      UpdateCutsString();
      return kTRUE;
    } else return kFALSE;
  case kTPCdEdxSigmaCutParam_e:
    if (SetTPCdEdxSigmaCut(AliPID::kElectron, value)) {
      fParameters[kTPCdEdxSigmaCutParam_e] = value;
      UpdateCutsString();
      return kTRUE;
    } else return kFALSE;
  case kTPCdEdxSigmaCutParam_mu:
    if (SetTPCdEdxSigmaCut(AliPID::kMuon, value)) {
      fParameters[kTPCdEdxSigmaCutParam_mu] = value;
      UpdateCutsString();
      return kTRUE;
    } else return kFALSE;
  case kTPCdEdxSigmaCutParam_pi:
    if (SetTPCdEdxSigmaCut(AliPID::kPion, value)) {
      fParameters[kTPCdEdxSigmaCutParam_pi] = value;
      UpdateCutsString();
      return kTRUE;
    } else return kFALSE;
  case kTPCdEdxSigmaCutParam_k:
    if (SetTPCdEdxSigmaCut(AliPID::kKaon, value)) {
      fParameters[kTPCdEdxSigmaCutParam_k] = value;
      UpdateCutsString();
      return kTRUE;
    } else return kFALSE;
  case kTPCdEdxSigmaCutParam_p:
    if (SetTPCdEdxSigmaCut(AliPID::kProton, value)) {
      fParameters[kTPCdEdxSigmaCutParam_p] = value;
      UpdateCutsString();
      return kTRUE;
    } else return kFALSE;
  case kTOFSigmaCutParam_e:
    if (SetTOFSigmaCut(AliPID::kElectron, value)) {
      fParameters[kTOFSigmaCutParam_e] = value;
      UpdateCutsString();
      return kTRUE;
    } else return kFALSE;
  case kTOFSigmaCutParam_mu:
    if (SetTOFSigmaCut(AliPID::kMuon, value)) {
      fParameters[kTOFSigmaCutParam_mu] = value;
      UpdateCutsString();
      return kTRUE;
    } else return kFALSE;
  case kTOFSigmaCutParam_pi:
    if (SetTOFSigmaCut(AliPID::kPion, value)) {
      fParameters[kTOFSigmaCutParam_pi] = value;
      UpdateCutsString();
      return kTRUE;
    } else return kFALSE;
  case kTOFSigmaCutParam_k:
    if (SetTOFSigmaCut(AliPID::kKaon, value)) {
      fParameters[kTOFSigmaCutParam_k] = value;
      UpdateCutsString();
      return kTRUE;
    } else return kFALSE;
  case kTOFSigmaCutParam_p:
    if (SetTOFSigmaCut(AliPID::kProton, value)) {
      fParameters[kTOFSigmaCutParam_p] = value;
      UpdateCutsString();
      return kTRUE;
    } else return kFALSE;
  default:
    AliError(Form("Cut param id %d out of supported range", paramID));
    return kFALSE;
  }
}


/// Print the whole cut iformatio for the cut ID
/// \param paramID the ID of the cut of interest
void AliCSPIDCuts::PrintCutWithParams(Int_t paramID) const {

  switch (cutsParametersIds(paramID)) {
  case kPRangeCutParam:
    printf("  Cut applicable P range: ");
    printf("%3.1f GeV < P %s", fMinP, (fMaxP < 9990) ? Form("%3.1f GeV\n", fMaxP) : "\n");
    break;
  case kITSdEdxSigmaCutParam_e:
    PrintITSdEdxSigmaCut(AliPID::kElectron);
    break;
  case kITSdEdxSigmaCutParam_mu:
    PrintITSdEdxSigmaCut(AliPID::kMuon);
    break;
  case kITSdEdxSigmaCutParam_pi:
    PrintITSdEdxSigmaCut(AliPID::kPion);
    break;
  case kITSdEdxSigmaCutParam_k:
    PrintITSdEdxSigmaCut(AliPID::kKaon);
    break;
  case kITSdEdxSigmaCutParam_p:
    PrintITSdEdxSigmaCut(AliPID::kProton);
    break;
  case kTPCdEdxSigmaCutParam_e:
    PrintTPCdEdxSigmaCut(AliPID::kElectron);
    break;
  case kTPCdEdxSigmaCutParam_mu:
    PrintTPCdEdxSigmaCut(AliPID::kMuon);
    break;
  case kTPCdEdxSigmaCutParam_pi:
    PrintTPCdEdxSigmaCut(AliPID::kPion);
    break;
  case kTPCdEdxSigmaCutParam_k:
    PrintTPCdEdxSigmaCut(AliPID::kKaon);
    break;
  case kTPCdEdxSigmaCutParam_p:
    PrintTPCdEdxSigmaCut(AliPID::kProton);
    break;
  case kTOFSigmaCutParam_e:
    PrintTOFSigmaCut(AliPID::kElectron);
    break;
  case kTOFSigmaCutParam_mu:
    PrintTOFSigmaCut(AliPID::kMuon);
    break;
  case kTOFSigmaCutParam_pi:
    PrintTOFSigmaCut(AliPID::kPion);
    break;
  case kTOFSigmaCutParam_k:
    PrintTOFSigmaCut(AliPID::kKaon);
    break;
  case kTOFSigmaCutParam_p:
    PrintTOFSigmaCut(AliPID::kProton);
    break;
  default:
    AliError(Form("Cut param id %d out of supported range", paramID));
  }
}

/// Configures the applicable track momentum range for the PID cut
/// \param ptcode the **P** range code
/// | code | minimum **P** (GeV/c) | maximum **P** (GeV/c) |
/// |:--:|:--:|:--:|
/// | 0 | full range | full range |
/// | 1 | 0.2 | 2.0 |
/// | 2 | 0.2 | 3.0 |
/// | 3 | 0.2 | 5.0 |
/// | 4 | 0.2 | 1.4 |
/// | 5 | 0.2 | 1.6 |
/// | 6 | 0.2 | 1.8 |
/// \return kTRUE if proper and supported \f$ p_{T} \f$ cut code
///
Bool_t AliCSPIDCuts::SetPRange(Int_t ptcode)
{
  switch(ptcode){
  case 0:
    fMinP = 0.0;
    fMaxP = 9999.0;
    break;
  case 1:
    fMinP = 0.2;
    fMaxP = 2.0;
    break;
  case 2:
    fMinP = 0.2;
    fMaxP = 3.0;
    break;
  case 3:
    fMinP = 0.2;
    fMaxP = 5.0;
    break;
  case 4:
    fMinP = 0.2;
    fMaxP = 1.4;
    break;
  case 5:
    fMinP = 0.2;
    fMaxP = 1.6;
    break;
  case 6:
    fMinP = 0.2;
    fMaxP = 1.8;
    break;
  default:
    AliError(Form("P range code %d not supported", ptcode));
    return kFALSE;
  }
  return kTRUE;
}

/// Sets the range for the dEdx \f$ n \sigma \f$ cut within the ITS
///
/// The cut establishes an acceptance band around a concrete
/// particle species line within the ITS. For species other than
/// the selected target the band is a separation band.
/// \param id the id of the species for which the band is being configured
/// \param dEdxCode the code for the \f$ n \sigma \f$ cut
/// | code | \f$ n \sigma \f$ below line | \f$ n \sigma \f$ above line | observations |
/// |:--:|:--:|:--:|:--|
/// | 0 | n/a | n/a | passive cut |
/// | 1 | -10 | 10 | |
/// | 2 | -6| 7 | |
/// | 3 | -5 | 5 | |
/// | 4 | -4 | 5 | |
/// | 5 | -3 | 5 | |
/// | 6 | -4 | 4 | |
/// | 7 | -2.5 | 4 | |
/// | 8 | -2 | 3.5 | |
/// \return kTRUE if proper and supported dEdx code

Bool_t AliCSPIDCuts::SetITSdEdxSigmaCut(AliPID::EParticleType id, Int_t dEdxCode){
  switch(dEdxCode){
    case 0:
      fITSEnabledSpeciesMask.ResetBitNumber(id);
      fITSnSigmaBelow[id] = -100.0;
      fITSnSigmaAbove[id] = 100.0;
      break;
    case 1:
      fITSEnabledSpeciesMask.SetBitNumber(id);
      fITSnSigmaBelow[id] = -10.0;
      fITSnSigmaAbove[id] = 10.0;
      break;
    case 2:
      fITSEnabledSpeciesMask.SetBitNumber(id);
      fITSnSigmaBelow[id] = -6.0;
      fITSnSigmaAbove[id] = 7.0;
      break;
    case 3:
      fITSEnabledSpeciesMask.SetBitNumber(id);
      fITSnSigmaBelow[id] = -5.0;
      fITSnSigmaAbove[id] = 5.0;
      break;
    case 4:
      fITSEnabledSpeciesMask.SetBitNumber(id);
      fITSnSigmaBelow[id] = -4.0;
      fITSnSigmaAbove[id] = 5.0;
      break;
    case 5:
      fITSEnabledSpeciesMask.SetBitNumber(id);
      fITSnSigmaBelow[id] = -3.0;
      fITSnSigmaAbove[id] = 5.0;
      break;
    case 6:
      fITSEnabledSpeciesMask.SetBitNumber(id);
      fITSnSigmaBelow[id] = -4.0;
      fITSnSigmaAbove[id] = 4.0;
      break;
    case 7:
      fITSEnabledSpeciesMask.SetBitNumber(id);
      fITSnSigmaBelow[id] = -2.5;
      fITSnSigmaAbove[id] = 4.0;
      break;
    case 8:
      fITSEnabledSpeciesMask.SetBitNumber(id);
      fITSnSigmaBelow[id] = -2.0;
      fITSnSigmaAbove[id] = 3.5;
      break;
    default:
      AliError(Form("ITS dEdx n sigmas cut code %d not supported",dEdxCode));
      return kFALSE;
  }
  if (fITSEnabledSpeciesMask.CountBits() > 0)
    fCutsEnabledMask.SetBitNumber(kITSdEdxSigmaCut);
  else
    fCutsEnabledMask.ResetBitNumber(kITSdEdxSigmaCut);
  return kTRUE;
}


/// Prints the dEdx \f$ n \sigma \f$ cut within the ITS
///
/// The cut establishes an acceptance band around a concrete
/// particle species line within the ITS. For species other than
/// the selected target the band is a separation band.
/// \param id the id of the species for which cut is being printed
void AliCSPIDCuts::PrintITSdEdxSigmaCut(AliPID::EParticleType id) const {
  if (fCutsEnabledMask.TestBitNumber(kITSdEdxSigmaCut)) {
    printf("  ITS PID CUT %s: ", AliPID::ParticleName(fTargetSpecies));
    if (fITSEnabledSpeciesMask.TestBitNumber(id)) {
      if (fTargetSpecies != id) {
        printf("nsigma %s < %3.1f OR %3.1f < nsigma %s\n",
            AliPID::ParticleName(id),
            fITSnSigmaBelow[id],
            fITSnSigmaAbove[id],
            AliPID::ParticleName(id));
      }
      else
        printf("%3.1f < nsigma < %3.1f\n",
            fITSnSigmaBelow[id],
            fITSnSigmaAbove[id]);
    }
    else
      printf("none to %s line\n",AliPID::ParticleName(id));
  }
  else
    printf("  ITS PID CUT %s: none\n", AliPID::ParticleName(fTargetSpecies));
}


/// Sets the range for the dEdx \f$ n \sigma \f$ cut within the TPC
///
/// The cut establishes an acceptance band around a concrete
/// particle species line within the TPC. For species other than
/// the selected target the band is a separation band.
/// \param id the id of the species for which the band is being configured
/// \param dEdxCode the code for the \f$ n \sigma \f$ cut
/// | code | \f$ n \sigma \f$ below line | \f$ n \sigma \f$ above line | observations |
/// |:--:|:--:|:--:|:--|
/// | 0 | n/a | n/a | passive cut |
/// | 1 | -10 | 10 | |
/// | 2 | -6| 7 | |
/// | 3 | -5 | 5 | |
/// | 4 | -4 | 5 | |
/// | 5 | -4 | 4 | |
/// | 6 | -3 | 4 | |
/// | 7 | -3 | 3 | |
/// | 8 | -3 | 5 | |
/// | 9 | -2 | 3 | |
/// \return kTRUE if proper and supported dEdx code

Bool_t AliCSPIDCuts::SetTPCdEdxSigmaCut(AliPID::EParticleType id, Int_t dEdxCode){
  AliInfo(Form("Configuring TPC dEdx cut for %s, with %d code", AliPID::ParticleName(id), dEdxCode));

  switch(dEdxCode){
    case 0:
      fTPCEnabledSpeciesMask.ResetBitNumber(id);
      fTPCnSigmaBelow[id] = -100.0;
      fTPCnSigmaAbove[id] = 100.0;
      break;
    case 1:
      fTPCEnabledSpeciesMask.SetBitNumber(id);
      fTPCnSigmaBelow[id] = -10.0;
      fTPCnSigmaAbove[id] = 10.0;
      break;
    case 2:
      fTPCEnabledSpeciesMask.SetBitNumber(id);
      fTPCnSigmaBelow[id] = -6.0;
      fTPCnSigmaAbove[id] = 7.0;
      break;
    case 3:
      fTPCEnabledSpeciesMask.SetBitNumber(id);
      fTPCnSigmaBelow[id] = -5.0;
      fTPCnSigmaAbove[id] = 5.0;
      break;
    case 4:
      fTPCEnabledSpeciesMask.SetBitNumber(id);
      fTPCnSigmaBelow[id] = -4.0;
      fTPCnSigmaAbove[id] = 5.0;
      break;
    case 5:
      fTPCEnabledSpeciesMask.SetBitNumber(id);
      fTPCnSigmaBelow[id] = -4.0;
      fTPCnSigmaAbove[id] = 4.0;
      break;
    case 6:
      fTPCEnabledSpeciesMask.SetBitNumber(id);
      fTPCnSigmaBelow[id] = -3.0;
      fTPCnSigmaAbove[id] = 4.0;
      break;
    case 7:
      fTPCEnabledSpeciesMask.SetBitNumber(id);
      fTPCnSigmaBelow[id] = -3.0;
      fTPCnSigmaAbove[id] = 3.0;
      break;
    case 8:
      fTPCEnabledSpeciesMask.SetBitNumber(id);
      fTPCnSigmaBelow[id] = -3.0;
      fTPCnSigmaAbove[id] = 5.0;
      break;
    case 9:
      fTPCEnabledSpeciesMask.SetBitNumber(id);
      fTPCnSigmaBelow[id] = -2.0;
      fTPCnSigmaAbove[id] = 3.0;
      break;
    default:
      AliError(Form("TPC dEdx n sigmas cut code %d not supported",dEdxCode));
      return kFALSE;
  }
  if (fTPCEnabledSpeciesMask.CountBits() > 0)
    fCutsEnabledMask.SetBitNumber(kTPCdEdxSigmaCut);
  else
    fCutsEnabledMask.ResetBitNumber(kTPCdEdxSigmaCut);
  return kTRUE;
}


/// Print the dEdx \f$ n \sigma \f$ cut within the TPC
///
/// The cut establishes an acceptance band around a concrete
/// particle species line within the TPC. For species other than
/// the selected target the band is a separation band.
/// \param id the id of the species for which the cut is being printed

void AliCSPIDCuts::PrintTPCdEdxSigmaCut(AliPID::EParticleType id) const {
  if (fCutsEnabledMask.TestBitNumber(kTPCdEdxSigmaCut)) {
    printf("  TPC PID CUT %s: ", AliPID::ParticleName(fTargetSpecies));
    if (fTPCEnabledSpeciesMask.TestBitNumber(id)) {
      if (fTargetSpecies != id) {
        printf("nsigma %s < %3.1f OR %3.1f < nsigma %s\n",
            AliPID::ParticleName(id),
            fTPCnSigmaBelow[id],
            fTPCnSigmaAbove[id],
            AliPID::ParticleName(id));
      }
      else
        printf("%3.1f < nsigma < %3.1f\n",
            fTPCnSigmaBelow[id],
            fTPCnSigmaAbove[id]);
    }
    else
      printf("none to %s line\n",AliPID::ParticleName(id));
  }
  else
    printf("  TPC PID CUT %s: none\n", AliPID::ParticleName(fTargetSpecies));
}


/// Sets the range for the \f$ n \, \sigma \f$ cut within TOF
///
/// The cut establishes an acceptance band around a concrete
/// particle species line within the TOF detector. For species other than
/// the selected target the band is a separation band.
/// \param id the id of the species for which the band is being configured
/// \param tofcode the code for the \f$ n \sigma \f$ cut
/// | code | \f$ n \sigma \f$ below line | \f$ n \sigma \f$ above line | observations |
/// |:--:|:--:|:--:|:--|
/// | 0 | n/a | n/a | passive cut |
/// | 1 | -7 | 7 | |
/// | 2 | -5 | 5 | |
/// | 3 | -3 | 5 | |
/// | 4 | -2 | 3 | |
/// | 5 | -3 | 3 | TOF required |
/// \return kTRUE if proper and supported TOF code

Bool_t AliCSPIDCuts::SetTOFSigmaCut(AliPID::EParticleType id, Int_t tofcode){
  switch(tofcode){
    case 0:
      fTOFEnabledSpeciesMask.ResetBitNumber(id);
      fTOFRequired[id] = kFALSE;
      fTOFnSigmaBelow[id] = -100.0;
      fTOFnSigmaAbove[id] = 100.0;
      break;
    case 1:
      fTOFEnabledSpeciesMask.SetBitNumber(id);
      fTOFRequired[id] = kFALSE;
      fTOFnSigmaBelow[id] = -7.0;
      fTOFnSigmaAbove[id] = 7.0;
      break;
    case 2:
      fTOFEnabledSpeciesMask.SetBitNumber(id);
      fTOFRequired[id] = kFALSE;
      fTOFnSigmaBelow[id] = -5.0;
      fTOFnSigmaAbove[id] = 5.0;
      break;
    case 3:
      fTOFEnabledSpeciesMask.SetBitNumber(id);
      fTOFRequired[id] = kFALSE;
      fTOFnSigmaBelow[id] = -3.0;
      fTOFnSigmaAbove[id] = 5.0;
      break;
    case 4:
      fTOFEnabledSpeciesMask.SetBitNumber(id);
      fTOFRequired[id] = kFALSE;
      fTOFnSigmaBelow[id] = -2.0;
      fTOFnSigmaAbove[id] = 3.0;
      break;
    case 5:
      fTOFEnabledSpeciesMask.SetBitNumber(id);
      fTOFRequired[id] = kTRUE;
      fTOFnSigmaBelow[id] = -3.0;
      fTOFnSigmaAbove[id] = 3.0;
      break;
    default:
      AliError(Form("TOF n sigmas cut code %d not supported", tofcode));
      return kFALSE;
  }
  if (fTOFEnabledSpeciesMask.CountBits() > 0)
    fCutsEnabledMask.SetBitNumber(kTOFSigmaCut);
  else
    fCutsEnabledMask.ResetBitNumber(kTOFSigmaCut);
  return kTRUE;
}

/// Prints the \f$ n \, \sigma \f$ cut within TOF
///
/// The cut establishes an acceptance band around a concrete
/// particle species line within the TOF detector. For species other than
/// the selected target the band is a separation band.
/// \param id the id of the species for which the cut is being printed
void AliCSPIDCuts::PrintTOFSigmaCut(AliPID::EParticleType id) const {
  if (fCutsEnabledMask.TestBitNumber(kTOFSigmaCut)) {
    printf("  TOF (%s) PID CUT %s: ", fTOFRequired ? "REQUIRED" : "NOT required", AliPID::ParticleName(fTargetSpecies));
    if (fTOFEnabledSpeciesMask.TestBitNumber(id)) {
      if (fTargetSpecies != id) {
        printf("nsigma %s < %3.1f OR %3.1f < nsigma %s\n",
            AliPID::ParticleName(id),
            fTOFnSigmaBelow[id],
            fTOFnSigmaAbove[id],
            AliPID::ParticleName(id));
      }
      else
        printf("%3.1f < nsigma < %3.1f\n",
            fTOFnSigmaBelow[id],
            fTOFnSigmaAbove[id]);
    }
    else
      printf("none to %s line\n",AliPID::ParticleName(id));
  }
  else
    printf("  TOF PID CUT %s: none\n", AliPID::ParticleName(fTargetSpecies));
}


/// Initializes the cuts
///
/// Initializes the needed data and allocates the needed histograms list if needed
/// \param name an additional name to precede the cuts string
void AliCSPIDCuts::InitCuts(const char *name){

  if (name == NULL) name = GetName();

  /* let's get the PID response instance */
  AliAnalysisManager *manager = AliAnalysisManager::GetAnalysisManager();
  if(manager != NULL) {
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (manager->GetInputEventHandler());
    fPIDResponse = (AliPIDResponse*) inputHandler->GetPIDResponse();
    /* if we need PID response instance and it is not there we cannot continue */
    if ((fPIDResponse == NULL) &&
        (fCutsEnabledMask.TestBitNumber(kITSdEdxSigmaCut) ||
            fCutsEnabledMask.TestBitNumber(kTPCdEdxSigmaCut) ||
            fCutsEnabledMask.TestBitNumber(kTOFSigmaCut)))
      AliFatal("No PID response instance. ABORTING!!!");
  }
  else {
    AliFatal("No analysis manager instance. ABORTING!!!");
  }

  if (fQALevel > kQALevelNone) {
    Bool_t oldstatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);

    if(fHistogramsList != NULL)
      delete fHistogramsList;

    fHistogramsList =new TList();
    fHistogramsList->SetOwner(kTRUE);
    fHistogramsList->SetName(name);

    TH1::AddDirectory(oldstatus);
  }
}

/// Allocates the different histograms if needed
///
/// It is supposed that the current cuts string is the running one
void AliCSPIDCuts::DefineHistograms(){

  if (fQALevel > kQALevelNone) {
    Bool_t oldstatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);

    /* the original name is used as title for the statistics histogram so, preserve it */
    TString originalTempName = fHistogramsList->GetName();
    fHistogramsList->SetName(Form("%s_%s",fHistogramsList->GetName(),GetCutsString()));

    fhCutsStatistics = new TH1F(Form("CutsStatistics_%s",GetCutsString()),Form("%s tracks cuts statistics",originalTempName.Data()),kNCuts+4,-0.5,kNCuts+3.5);
    fhCutsStatistics->GetXaxis()->SetBinLabel(1,"n tracks");
    fhCutsStatistics->GetXaxis()->SetBinLabel(2,"n cut tracks");
    for (Int_t i = 0; i < kNCuts; i++)
      fhCutsStatistics->GetXaxis()->SetBinLabel(i+4, fgkCutsNames[i]);
    fHistogramsList->Add(fhCutsStatistics);

    if(fQALevel == kQALevelHeavy){
      fhCutsCorrelation = new TH2F(Form("CutCorrelation_%s",GetCutsString()),"Cuts correlation",kNCuts+2,-0.5,kNCuts+1.5,kNCuts+2,-0.5,kNCuts+1.5);
      for (Int_t i=0; i<kNCuts; i++) {
        fhCutsCorrelation->GetXaxis()->SetBinLabel(i+2,fgkCutsNames[i]);
        fhCutsCorrelation->GetYaxis()->SetBinLabel(i+2,fgkCutsNames[i]);
      }
      fHistogramsList->Add(fhCutsCorrelation);
    }

    /* build the P bins */
    const Int_t nPbins = 150;
    Double_t minP = 0.05;
    Double_t maxP = 20.0;
    Double_t *edges = new Double_t [nPbins + 1];
    Double_t factor = TMath::Power(maxP / minP, 1. / nPbins);
    edges[0] = minP; for (Int_t bin = 0; bin < nPbins; bin++) edges[bin+1] = factor * edges[bin];

    fhITSdEdxSigmaVsP[0] = new TH2F(Form("ITSdEdxSigmaB_%s",GetCutsString()),"ITS dE/dx n#sigma before;P (GeV/c);n#sigma",nPbins,edges, 400, -10, 10);
    fhITSdEdxSigmaVsP[1] = new TH2F(Form("ITSdEdxSigmaA_%s",GetCutsString()),"ITS dE/dx n#sigma;P (GeV/c);n#sigma",nPbins,edges, 400, -10, 10);
    fHistogramsList->Add(fhITSdEdxSigmaVsP[0]);
    fHistogramsList->Add(fhITSdEdxSigmaVsP[1]);

    fhITSdEdxSignalVsP[0] = new TH2F(Form("ITSdEdxSignalB_%s",GetCutsString()),"ITS dE/dx signal before;P (GeV/c);#frac{dE}{dx} (au)",nPbins,edges, 800, 0.0, 200);
    fhITSdEdxSignalVsP[1] = new TH2F(Form("ITSdEdxSignalA_%s",GetCutsString()),"ITS dE/dx signal;P (GeV/c);#frac{dE}{dx} (au)",nPbins,edges, 800, 0.0, 200);
    fHistogramsList->Add(fhITSdEdxSignalVsP[0]);
    fHistogramsList->Add(fhITSdEdxSignalVsP[1]);

    fhTPCdEdxSigmaVsP[0] = new TH2F(Form("TPCdEdxSigmaB_%s",GetCutsString()),"TPC dE/dx n#sigma before;P (GeV/c);n#sigma",nPbins,edges, 400, -10, 10);
    fhTPCdEdxSigmaVsP[1] = new TH2F(Form("TPCdEdxSigmaA_%s",GetCutsString()),"TPC dE/dx n#sigma;P (GeV/c);n#sigma",nPbins,edges, 400, -10, 10);
    fHistogramsList->Add(fhTPCdEdxSigmaVsP[0]);
    fHistogramsList->Add(fhTPCdEdxSigmaVsP[1]);

    fhTPCdEdxSignalVsP[0] = new TH2F(Form("TPCdEdxSignalB_%s",GetCutsString()),"TPC dE/dx signal before;P (GeV/c);#frac{dE}{dx} (au)",nPbins,edges, 800, 0.0, 200.0);
    fhTPCdEdxSignalVsP[1] = new TH2F(Form("TPCdEdxSignalA_%s",GetCutsString()),"TPC dE/dx signal;P (GeV/c);#frac{dE}{dx} (au)",nPbins,edges, 800, 0.0, 200.0);
    fHistogramsList->Add(fhTPCdEdxSignalVsP[0]);
    fHistogramsList->Add(fhTPCdEdxSignalVsP[1]);

    fhTOFSigmaVsP[0] = new TH2F(Form("TOFSigmaB_%s",GetCutsString()),"TOF n#sigma before;P (GeV/c);n#sigma",nPbins,edges, 400, -10, 10);
    fhTOFSigmaVsP[1] = new TH2F(Form("TOFSigmaA_%s",GetCutsString()),"TOF n#sigma;P (GeV/c);n#sigma",nPbins,edges, 400, -10, 10);
    fHistogramsList->Add(fhTOFSigmaVsP[0]);
    fHistogramsList->Add(fhTOFSigmaVsP[1]);

    fhTOFSignalVsP[0] = new TH2F(Form("TOFSignalB_%s",GetCutsString()),"TOF signal before;P (GeV/c);#beta",nPbins,edges, 400, 0.0, 1.1);
    fhTOFSignalVsP[1] = new TH2F(Form("TOFSignalA_%s",GetCutsString()),"TOF signal;P (GeV/c);#beta",nPbins,edges, 400, 0.0, 1.1);
    fHistogramsList->Add(fhTOFSignalVsP[0]);
    fHistogramsList->Add(fhTOFSignalVsP[1]);

    TH1::AddDirectory(oldstatus);
  }
}

