/**************************************************************************
* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved.  *
*                                                                         *
* Permission to use, copy, modify and distribute this software and its    *
* documentation strictly for non-commercial purposes is hereby granted    *
* without fee, provided that the above copyright notice appears in all    *
* copies and that both the copyright notice and this permission notice    *
* appear in the supporting documentation. The authors make no claims      *
* about the suitability of this software for any purpose. It is           *
* provided "as is" without express or implied warranty.                   *
**************************************************************************/

#include <TFormula.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliCSEventCuts.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliAODMCParticle.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliCentrality.h"
#include "AliMultSelection.h"
#include "AliESDtrackCuts.h"
#include "AliLog.h"

/// \file AliCSEventCuts.cxx
/// \brief Implementation of event cuts class within the correlation studies analysis

const char *AliCSEventCuts::fgkCutsNames[AliCSEventCuts::kNCuts] = {
    "MC without proper MC data",
    "DAQ incomplete",
    "No tracks",
    "Offline trigger",
    "Vertex contributors",
    "Vertex quality",
    "SPD and tracks vertex distance",
    "z_{vtx}",
    "Pile up",
    "2015 Pile up",
    "SPD clusters vs tracklets",
    "Centrality"
};

Float_t AliCSEventCuts::fgkVertexResolutionThreshold = 0.25;
Float_t AliCSEventCuts::fgkVertexDispersionThreshold = 0.03;
Float_t AliCSEventCuts::fgkSPDTracksVtxDistanceThreshold = 0.2;
Float_t AliCSEventCuts::fgkSPDTracksVtxDistanceSigmas = 10.0;
Float_t AliCSEventCuts::fgkTrackVertexSigmas = 20.0;

/// Default constructor for serialization
AliCSEventCuts::AliCSEventCuts() :
    AliCSAnalysisCutsBase(),
    fSystem(kNoSystem),
    fVertexZ(999.0),
    fCentrality(0.0),
    fCentralityDetector(0),
    fCentralityModifier(0),
    fCentralityMin(0.0),
    fCentralityMax(1e6),
    fOfflineTriggerMask(AliVEvent::kAny),
    fMaxVertexZ(1000.0),
    fUseSPDTracksVtxDist(kFALSE),
    fVertexResolutionTh(fgkVertexResolutionThreshold),
    fVertexDispersionTh(fgkVertexDispersionThreshold),
    fUseNewMultFramework(kFALSE),
    f2015V0MtoTrkTPCout(NULL),
    fCentOutLowCut(NULL),
    fCentOutHighCut(NULL),
    fTOFMultOutLowCut(NULL),
    fTOFMultOutHighCut(NULL),
    fMultCentOutLowCut(NULL),
    fV0MCentrality(0),
    fV0ACentrality(0),
    fV0CCentrality(0),
    fCL0Centrality(0),
    fCL1Centrality(0),
    fReferenceMultiplicity(0),
    fV0Multiplicity(0),
    fNoOfAODTracks(0),
    fNoOfESDTracks(0),
    fNoOfFB32Tracks(0),
    fNoOfFB128Tracks(0),
    fNoOfFB32AccTracks(0),
    fNoOfFB32TOFTracks(0),
    fNoOfTPCoutTracks(0),
    fNoOfInitialTPCoutTracks(0),
    fAnalysisUtils(),
    fESDFB32(NULL),
    fESDFB128(NULL),
    fhCutsStatistics(NULL),
    fhUniqueCutsStatistics(NULL),
    fhCutsCorrelation(NULL),
    fhCentrality{NULL},
    fhVertexZ{NULL},
    fhSPDClustersVsTracklets{NULL},
    fhV0MvsTracksTPCout{NULL},
    fhV0MvsTracksInitialTPCout{NULL},
    fhCentralityAltVsSel{NULL},
    fhCL0vsV0MCentrality{NULL},
    fhESDvsTPConlyMultiplicity{NULL},
    fhTOFvsGlobalMultiplicity{NULL},
    fhAccTrkvsV0MCentrality{NULL},
    fhTriggerClass{NULL}
{
}

/// Constructor
/// \param name name of the event cuts
/// \param title title of the event cuts
AliCSEventCuts::AliCSEventCuts(const char *name, const char *title) :
    AliCSAnalysisCutsBase(kNCuts,kNCutsParameters,name,title),
    fSystem(kNoSystem),
    fVertexZ(999.0),
    fCentrality(0.0),
    fCentralityDetector(0),
    fCentralityModifier(0),
    fCentralityMin(0.0),
    fCentralityMax(1e6),
    fOfflineTriggerMask(AliVEvent::kAny),
    fMaxVertexZ(1000.0),
    fUseSPDTracksVtxDist(kFALSE),
    fVertexResolutionTh(fgkVertexResolutionThreshold),
    fVertexDispersionTh(fgkVertexDispersionThreshold),
    fUseNewMultFramework(kFALSE),
    f2015V0MtoTrkTPCout(NULL),
    fCentOutLowCut(NULL),
    fCentOutHighCut(NULL),
    fTOFMultOutLowCut(NULL),
    fTOFMultOutHighCut(NULL),
    fMultCentOutLowCut(NULL),
    fV0MCentrality(0),
    fV0ACentrality(0),
    fV0CCentrality(0),
    fCL0Centrality(0),
    fCL1Centrality(0),
    fReferenceMultiplicity(0),
    fV0Multiplicity(0),
    fNoOfAODTracks(0),
    fNoOfESDTracks(0),
    fNoOfFB32Tracks(0),
    fNoOfFB128Tracks(0),
    fNoOfFB32AccTracks(0),
    fNoOfFB32TOFTracks(0),
    fNoOfTPCoutTracks(0),
    fNoOfInitialTPCoutTracks(0),
    fAnalysisUtils(),
    fESDFB32(NULL),
    fESDFB128(NULL),
    fhCutsStatistics(NULL),
    fhUniqueCutsStatistics(NULL),
    fhCutsCorrelation(NULL),
    fhCentrality{NULL},
    fhVertexZ{NULL},
    fhSPDClustersVsTracklets{NULL},
    fhV0MvsTracksTPCout{NULL},
    fhV0MvsTracksInitialTPCout{NULL},
    fhCentralityAltVsSel{NULL},
    fhCL0vsV0MCentrality{NULL},
    fhESDvsTPConlyMultiplicity{NULL},
    fhTOFvsGlobalMultiplicity{NULL},
    fhAccTrkvsV0MCentrality{NULL},
    fhTriggerClass{NULL}
{
}

/// Destructor
AliCSEventCuts::~AliCSEventCuts()
{
  if (f2015V0MtoTrkTPCout != NULL)
    delete f2015V0MtoTrkTPCout;
  if (fCentOutLowCut != NULL)
    delete fCentOutLowCut;
  if (fCentOutHighCut != NULL)
    delete fCentOutHighCut;
  if (fTOFMultOutLowCut != NULL)
    delete fTOFMultOutLowCut;
  if (fTOFMultOutHighCut != NULL)
    delete fTOFMultOutHighCut;
  if (fMultCentOutLowCut != NULL)
    delete fMultCentOutLowCut;
  if (fESDFB32 != NULL)
    delete fESDFB32;
  if (fESDFB128 != NULL)
    delete fESDFB128;
}

/// Processes a potential change in the run number
///
/// Checks if the current period under analysis has changed and if so
/// updates the needed members. Finally, once the running cuts string is
/// definitely known, asks for histograms allocation.
void AliCSEventCuts::NotifyRun() {

  /* checks the change in the analysis period */
  if (AliCSAnalysisCutsBase::GetGlobalPeriod() != fDataPeriod) {

    fDataPeriod = AliCSAnalysisCutsBase::GetGlobalPeriod();

    /* set the system type according to the data period */
    SetActualSystemType();

    /* set the trigger according to the data period */
    SetActualActiveTrigger();

    /* set the 2015 pileup rejection according to the data period */
    SetActual2015PileUpRemoval();

    /* we adapt the different cuts accordingly with the data period */
    fUseNewMultFramework = UseNewMultiplicityFramework();

    /* set the counting tracks filter cuts according to data period */
    SetActualFilterTracksCuts();

    /* and now we ask for histogram allocation */
    DefineHistograms();
  }
}

/// A new event is starting to be analyzed
/// Store MC needed data in case of AOD format
void AliCSEventCuts::NotifyEvent() {
  /* let's produce some feedback about MC dataset configuration */
  if (fgIsMC) {
    AliInfo(Form("Handling a MC event with %s format", (fgIsESD ? "ESD" : "AOD")));
    AliInfo(Form("========= MC handler is %s", ((fgMCHandler != NULL) ? "NOT null" : "NULL")));
  }

  if (fgIsMC && !fgIsESD) {
    fgMCArray = dynamic_cast<TClonesArray*>(((AliAODEvent *)fgInputHandler->GetEvent())->FindListObject(AliAODMCParticle::StdBranchName()));
  }
}

/// Check whether the passed event is accepted by the different cuts
/// \param fInputEvent the current event to evaluate
/// \return kTRUE if the event is accepted, kFALSE otherwise
Bool_t AliCSEventCuts::IsEventAccepted(AliVEvent *fInputEvent) {
  Bool_t accepted = kTRUE;

  /* initialize the mask of activated cuts */
  fCutsActivatedMask.ResetAllBits();

  /* check for MC event and its quality */
  if (fgIsMC) {
    if (fgIsESD) {
      if (!fgMCHandler->InitOk() || !fgMCHandler->TreeK() || !fgMCHandler->TreeTR()){
        fCutsActivatedMask.SetBitNumber(kMCdataQuality);
        accepted = kFALSE;
      }
    }
    else {
      TClonesArray *arrayMC = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
      if (arrayMC == NULL) {
        fCutsActivatedMask.SetBitNumber(kMCdataQuality);
        accepted = kFALSE;
      }
    }
  }

  /* is the event information complete */
  if (fInputEvent->IsIncompleteDAQ()) {
    fCutsActivatedMask.SetBitNumber(kDAQIncompleteCut);
    accepted = kFALSE;
  }

  StoreEventCentralities(fInputEvent);
  StoreEventMultiplicities(fInputEvent);

  /* is the event having any track */
  if (fReferenceMultiplicity < 0){
    fCutsActivatedMask.SetBitNumber(kNoTracks);
    accepted = kFALSE;
  }

  /* trigger cut */
  UInt_t offlineTriggerMask = fOfflineTriggerMask;
  if (fgIsMC && fSystem != kpp)
    offlineTriggerMask = AliVEvent::kAny;
  if (!(offlineTriggerMask & fgInputHandler->IsEventSelected())) {
    /* cut activated, event rejected */
    fCutsActivatedMask.SetBitNumber(kOfflineTriggerCut);
    accepted = kFALSE;
  }

  /* vertex cut */
  fVertexZ = GetVertexZ(fInputEvent);
  if (fCutsEnabledMask.TestBitNumber(kVertexCut)) {
    if (GetNumberOfVertexContributors(fInputEvent) < 1) {
      fCutsActivatedMask.SetBitNumber(kVertexContributorsCut);
      accepted = kFALSE;
    }
    if (!PassVertexResolutionAndDispersionTh(fInputEvent)) {
      fCutsActivatedMask.SetBitNumber(kVertexQualityCut);
      accepted = kFALSE;
    }
    if (fUseSPDTracksVtxDist) {
      if (!AcceptSPDTracksVtxDist(fInputEvent)) {
        fCutsActivatedMask.SetBitNumber(kSPDTrackVtxDistance);
        accepted = kFALSE;
      }
    }
    if (fMaxVertexZ < TMath::Abs(fVertexZ)) {
      fCutsActivatedMask.SetBitNumber(kVertexCut);
      accepted = kFALSE;
    }
  }

  /* pile up cut */
  if (fSystem == kpPb) {
    if (fAnalysisUtils.IsFirstEventInChunk(fInputEvent)) {
      fCutsActivatedMask.SetBitNumber(kPileUpCut);
      accepted = kFALSE;
    }
  }
  if (fCutsEnabledMask.TestBitNumber(kPileUpCut)) {
    if(fAnalysisUtils.IsPileUpEvent(fInputEvent)){
      fCutsActivatedMask.SetBitNumber(kPileUpCut);
      accepted = kFALSE;
    }
    if (fAnalysisUtils.IsSPDClusterVsTrackletBG(fInputEvent)){
      fCutsActivatedMask.SetBitNumber(kSPDClsVsTrkaletsCut);
      accepted = kFALSE;
    }
  }

  /* centrality cut */
  fCentrality = GetEventCentrality(fInputEvent);
  AliInfo(Form("Event centrality: %f", Float_t(fCentrality)));
  if (fCutsEnabledMask.TestBitNumber(kCentralityCut)) {
    if (fCentrality < fCentralityMin || fCentralityMax <= fCentrality ) {
      fCutsActivatedMask.SetBitNumber(kCentralityCut);
      accepted = kFALSE;
    }
  }

  /* 2015 additional pile up cut also applicable to 2010h*/
  /* check the additional pile up rejection if required */
  if (fCutsEnabledMask.TestBitNumber(k2015PileUpCut)) {
    if (Is2015PileUpEvent()) {
      fCutsActivatedMask.SetBitNumber(k2015PileUpCut);
      accepted = kFALSE;
    }
  }

  if (fQALevel > kQALevelNone) {
    /* let's fill the histograms */
    fhCutsStatistics->Fill(fhCutsStatistics->GetBinCenter(fhCutsStatistics->GetXaxis()->FindBin("n events")));
    fhUniqueCutsStatistics->Fill(fhUniqueCutsStatistics->GetBinCenter(fhUniqueCutsStatistics->GetXaxis()->FindBin("n events")));
    if (!accepted)
      fhCutsStatistics->Fill(fhCutsStatistics->GetBinCenter(fhCutsStatistics->GetXaxis()->FindBin("n cut events")));
    else
      fhUniqueCutsStatistics->Fill(fhUniqueCutsStatistics->GetBinCenter(fhUniqueCutsStatistics->GetXaxis()->FindBin("n passed events")));

    for (Int_t i=0; i<kNCuts; i++) {
      if (fCutsActivatedMask.TestBitNumber(i)) {
        if ((fhCutsStatistics->GetXaxis()->FindBin(fgkCutsNames[i]) < 1) || (fhUniqueCutsStatistics->GetXaxis()->FindBin(fgkCutsNames[i]) < 1))
          AliFatal(Form("Inconsistency! Cut %d with name %s not found", i, fgkCutsNames[i]));

        fhCutsStatistics->Fill(fhCutsStatistics->GetBinCenter(fhCutsStatistics->GetXaxis()->FindBin(fgkCutsNames[i])));
        if (fCutsActivatedMask.CountBits() == 1)
          fhUniqueCutsStatistics->Fill(fhUniqueCutsStatistics->GetBinCenter(fhUniqueCutsStatistics->GetXaxis()->FindBin(fgkCutsNames[i])));
      }

      if (fQALevel > kQALevelLight) {
        for (Int_t j=i; j<kNCuts; j++) {
          if (fCutsActivatedMask.TestBitNumber(i) && fCutsActivatedMask.TestBitNumber(j)) {
            if (fhCutsStatistics->GetXaxis()->FindBin(fgkCutsNames[j]) < 1)
              AliFatal(Form("Inconsistency! Cut %d with name %s not found", j, fgkCutsNames[j]));

            Float_t xC = fhCutsCorrelation->GetXaxis()->GetBinCenter(fhCutsCorrelation->GetXaxis()->FindBin(fgkCutsNames[i]));
            Float_t yC = fhCutsCorrelation->GetYaxis()->GetBinCenter(fhCutsCorrelation->GetYaxis()->FindBin(fgkCutsNames[j]));
            fhCutsCorrelation->Fill(xC, yC);
          }
        }
      }
    }

    /* prepare additional data in case they were needed */
    Int_t nClustersLayer0 = fInputEvent->GetNumberOfITSClusters(0);
    Int_t nClustersLayer1 = fInputEvent->GetNumberOfITSClusters(1);
    Int_t nTracklets      = fInputEvent->GetMultiplicity()->GetNumberOfTracklets();
    Float_t altcentrality = this->GetEventAltCentrality(fInputEvent);

    for (Int_t i = 0; i < 2; i++) {

      fhCentrality[i]->Fill(fCentrality);
      fhVertexZ[i]->Fill(fVertexZ);

      if (fQALevel > kQALevelLight) {
        fhSPDClustersVsTracklets[i]->Fill(nTracklets, nClustersLayer0+nClustersLayer1);
        fhV0MvsTracksTPCout[i]->Fill(fNoOfTPCoutTracks, fV0Multiplicity);
        fhV0MvsTracksInitialTPCout[i]->Fill(fNoOfInitialTPCoutTracks, fV0Multiplicity);
        fhCentralityAltVsSel[i]->Fill(fCentrality,altcentrality);
        fhCL0vsV0MCentrality[i]->Fill(fV0MCentrality,fCL0Centrality);
        fhESDvsTPConlyMultiplicity[i]->Fill(fNoOfFB128Tracks,fNoOfESDTracks);
        fhTOFvsGlobalMultiplicity[i]->Fill(fNoOfFB32Tracks,fNoOfFB32TOFTracks);
        fhAccTrkvsV0MCentrality[i]->Fill(fV0MCentrality,fNoOfFB32AccTracks);
      }

      /* the trigger histogram */
      UInt_t selectedTrigger = fgInputHandler->IsEventSelected();
      if (selectedTrigger & AliVEvent::kMB)fhTriggerClass[i]->Fill(0);
      if (selectedTrigger & AliVEvent::kINT7)fhTriggerClass[i]->Fill(1);
      if (selectedTrigger & AliVEvent::kMUON)fhTriggerClass[i]->Fill(2);
      if (selectedTrigger & AliVEvent::kHighMult)fhTriggerClass[i]->Fill(3);
      if (selectedTrigger & AliVEvent::kEMC1)fhTriggerClass[i]->Fill(4);
      if (selectedTrigger & AliVEvent::kCINT5)fhTriggerClass[i]->Fill(5);
      if (selectedTrigger & AliVEvent::kCMUS5)fhTriggerClass[i]->Fill(6);
      if (selectedTrigger & AliVEvent::kMUSH7)fhTriggerClass[i]->Fill(7);
      if (selectedTrigger & AliVEvent::kMUL7)fhTriggerClass[i]->Fill(8);
      if (selectedTrigger & AliVEvent::kMUU7)fhTriggerClass[i]->Fill(9);
      if (selectedTrigger & AliVEvent::kEMC7)fhTriggerClass[i]->Fill(10);
      if (selectedTrigger & AliVEvent::kMUS7)fhTriggerClass[i]->Fill(11);
      if (selectedTrigger & AliVEvent::kPHI1)fhTriggerClass[i]->Fill(12);
      if (selectedTrigger & AliVEvent::kPHI7)fhTriggerClass[i]->Fill(13);
      if (selectedTrigger & AliVEvent::kEMCEJE)fhTriggerClass[i]->Fill(14);
      if (selectedTrigger & AliVEvent::kEMCEGA)fhTriggerClass[i]->Fill(15);
      if (selectedTrigger & AliVEvent::kCentral)fhTriggerClass[i]->Fill(16);
      if (selectedTrigger & AliVEvent::kSemiCentral)fhTriggerClass[i]->Fill(17);
      if (selectedTrigger & AliVEvent::kDG5)fhTriggerClass[i]->Fill(18);
      if (selectedTrigger & AliVEvent::kZED)fhTriggerClass[i]->Fill(19);
      if (selectedTrigger & AliVEvent::kSPI7)fhTriggerClass[i]->Fill(20);
      if (selectedTrigger & AliVEvent::kINT8)fhTriggerClass[i]->Fill(21);
      if (selectedTrigger & AliVEvent::kMuonSingleLowPt8)fhTriggerClass[i]->Fill(22);
      if (selectedTrigger & AliVEvent::kMuonSingleHighPt8)fhTriggerClass[i]->Fill(23);
      if (selectedTrigger & AliVEvent::kMuonLikeLowPt8)fhTriggerClass[i]->Fill(24);
      if (selectedTrigger & AliVEvent::kMuonUnlikeLowPt8)fhTriggerClass[i]->Fill(25);
      if (selectedTrigger & AliVEvent::kMuonUnlikeLowPt0)fhTriggerClass[i]->Fill(26);
      if (selectedTrigger & AliVEvent::kUserDefined)fhTriggerClass[i]->Fill(27);
      if (selectedTrigger & AliVEvent::kTRD)fhTriggerClass[i]->Fill(28);
      if (selectedTrigger & AliVEvent::kFastOnly)fhTriggerClass[i]->Fill(29);
      if (selectedTrigger & AliVEvent::kAnyINT)fhTriggerClass[i]->Fill(30);
      if (selectedTrigger & AliVEvent::kAny)fhTriggerClass[i]->Fill(31);
      if (!(selectedTrigger & AliVEvent::kFastOnly))fhTriggerClass[i]->Fill(32);
      if (selectedTrigger == 0x0)fhTriggerClass[i]->Fill(33);

      /* don't fill after if event not accepted */
      if (!accepted) break;
    }
  }
  return accepted;
}



/// Sets the individual value for the cut parameter ID
/// \param paramID the ID of the cut parameter of interest
/// \param value the value for the cut parameter
/// \return kTRUE if the cut parameter value was accepted
Bool_t AliCSEventCuts::SetCutAndParams(Int_t paramID, Int_t value) {

  switch (cutsParametersIds(paramID)) {
  case kSystem:
    if( SetSystemType(SystemType(value))) {
      fParameters[kSystem] = value;
      UpdateCutsString();
      return kTRUE;
    } else return kFALSE;
  case kCentralityType:
    if( SetCentralityType(value)) {
      fParameters[kCentralityType] = value;
      UpdateCutsString();
      return kTRUE;
    } else return kFALSE;
  case kCentralityMin:
    if( SetCentralityMin(value)) {
      fParameters[kCentralityMin] = value;
      UpdateCutsString();
      return kTRUE;
    } else return kFALSE;
  case kCentralityMax:
    if( SetCentralityMax(value)) {
      fParameters[kCentralityMax] = value;
      UpdateCutsString();
      return kTRUE;
    } else return kFALSE;
  case kActiveTrigger:
    if( SetSelectActiveTrigger(value)) {
      fParameters[kActiveTrigger] = value;
      UpdateCutsString();
      return kTRUE;
    } else return kFALSE;
  case kActiveSubTrigger:
    if( SetSelectActiveSubTrigger(value)) {
      fParameters[kActiveSubTrigger] = value;
      UpdateCutsString();
      return kTRUE;
    } else return kFALSE;
  case kRemovePileUp:
    if( SetRemovePileUp(value)) {
      fParameters[kRemovePileUp] = value;
      UpdateCutsString();
      return kTRUE;
    } else return kFALSE;
  case kVertex:
    if( SetVertexCut(value)) {
      fParameters[kVertex] = value;
      UpdateCutsString();
      return kTRUE;
    } else return kFALSE;
  case kRemove2015PileUp:
    if( SetRemove2015PileUp(value)) {
      fParameters[kRemove2015PileUp] = value;
      UpdateCutsString();
      return kTRUE;
    } else return kFALSE;
  default:
    AliError(Form("Cut id %d out of range", paramID));
    return kFALSE;
  }
}


/// Print the whole cut information for the cut parameter ID
/// \param paramID the ID of the cut of interest
void AliCSEventCuts::PrintCutWithParams(Int_t paramID) const {

  switch (cutsParametersIds(paramID)) {
  case kSystem:
    switch(fSystem) {
    case kNoSystem:
      printf("SYSTEM: none\n");
      break;
    case kpp:
      printf("SYSTEM: p-p\n");
      break;
    case kpPb:
      printf("SYSTEM: p-Pb\n");
      break;
    case kPbPb:
      printf("SYSTEM: Pb-Pb\n");
      break;
    default:
      printf("SYSTEM: no proper system configured\n");
    }
    break;
  case kCentralityType:
    if (fCutsEnabledMask.TestBitNumber(kCentralityCut)) {
      TString szEstimator;
      if (fCentralityDetector == 1)
        szEstimator = "CL1";
      else
        if (fSystem == kpPb)
          szEstimator = "V0A";
        else
          szEstimator = "V0M";

      printf("  Centrality cut: using %s. ", szEstimator.Data());
    }
    else {
      printf("  Centrality cut: not cutting in centrality\n");
    }
    break;
  case kCentralityMin:
    if (fCutsEnabledMask.TestBitNumber(kCentralityCut)) {
      printf("Range: %.1f-", fCentralityMin);
    }
    break;
  case kCentralityMax:
    if (fCutsEnabledMask.TestBitNumber(kCentralityCut)) {
      printf("%.1f\n", fCentralityMax);
    }
    break;
  case kActiveTrigger: {
      /* this is a first try to print trigger information for the current */
      /* reduced set of triggers supported. Smarter version should come */
      printf("  Offline trigger: ");
      UInt_t triggerPrintedMask = 0x0;
      PrintTrigger(triggerPrintedMask, AliVEvent::kAny, "kAny");
      PrintTrigger(triggerPrintedMask, AliVEvent::kMB, "kMB");
      PrintTrigger(triggerPrintedMask, AliVEvent::kINT7, "kINT7");
      PrintTrigger(triggerPrintedMask, AliVEvent::kDG5, "kDG5");
    }
    break;
  case kActiveSubTrigger: {
      /* no sub-trigger supported yet */
      printf("    sub-trigger: none\n");
    }
    break;
  case kRemovePileUp:
    /* this is a bit tricky because AnalysisUtils does not allow to get its configuration */
    printf("  Pile up removal:\n");
    switch(fParameters[kRemovePileUp]) {
    case 0:
      printf("    none\n");
      break;
    case 1: /* use SPD multiple vertices */
      printf("    using SPD multiple vertices\n");
      printf("    background rejection with cluster-vs-tracklet default cut\n");
      break;
    case 2: /* use SPD multiple vertices  & out-of-bunch pileup rejection */
      printf("    using SPD multiple vertices\n");
      printf("    background rejection with cluster-vs-tracklet default cut\n");
      printf("    using out-of-bunch pileup rejection\n");
      break;
    case 3: /* use track multivertexer */
      printf("    using track multivertexer\n");
      printf("    background rejection with cluster-vs-tracklet default cut\n");
      break;
    case 4: /* use track multivertexer & out-of-bunch pileup rejection */
      printf("    using track multivertexer\n");
      printf("    background rejection with cluster-vs-tracklet default cut\n");
      printf("    using out-of-bunch pileup rejection\n");
      break;
    default:
      printf("    Pile up removal procedure %d not supported", fParameters[kRemovePileUp]);
    }
    break;
  case kVertex:
    printf("  Z vertex cut: ");
    if (fCutsEnabledMask.TestBitNumber(kVertexCut)) {
      printf("|Zvtx| < %.1f\n", fMaxVertexZ);
      if (fUseSPDTracksVtxDist) {
        printf("    Cutting on distance between SPD and track vertices\n");
      }
    }
    else
      printf("none\n");
    break;
  case kRemove2015PileUp:
    /* this is a bit tricky because AnalysisUtils does not allow to get its configuration */
    printf("  2015 additional pile up removal:\n");
    switch(fParameters[kRemove2015PileUp]) {
    case 0:
      printf("    none\n");
      break;
    case 1: /* use J/psi 2015 pile up rejection initial, faulty, method */
      printf("    using J/psi 2015 pile up rejection, initial (faulty) method of track counting\n");
      printf("    actual cut will depend on data period\n");
      break;
    case 2: /* Centrality and multiplicity correlations for 2015*/
      printf("    using centrality and multiplicity correlations for 2015 pile up rejection\n");
      printf("    actual cut will depend on data period\n");
      break;
    case 3: /* use J/psi 2015 pile up rejection initial, corrected, method*/
      printf("    using J/psi 2015 pile up rejection , initial (corrected) method of track counting\n");
      printf("    actual cut will depend on data period\n");
      break;
    default:
      printf("    2015 additional pile up removal procedure %d not supported\n", fParameters[kRemove2015PileUp]);
    }
    break;
  default:
    AliError(Form("Cut id %d out of range", paramID));
  }
}

/// Prints trigger name
///
/// Trigger name is only printed if active.
/// A check is made to print the continuation after printing the trigger name.
/// The printed mask is updated with the current trigger if active and printed.
/// \param printed the already printed triggers
/// \param trigger the mask associated with the trigger to print
/// \param name the name of the interested trigger
void AliCSEventCuts::PrintTrigger(UInt_t &printed, UInt_t trigger, const char *name) const {
  /* check if done */
  if (fOfflineTriggerMask == printed) return;
  /* check if active */
  if ((fOfflineTriggerMask & trigger) == trigger) {
    printf("%s", name);
    printed = printed | trigger;
    if ((fOfflineTriggerMask & printed) == fOfflineTriggerMask)
      printf("\n"); /* last trigger */
      if (trigger == AliVEvent::kAny) {
        printf("    could depend on data period\n");
      }
    else
      printf("+"); /* still more triggers to print */
  }
}

/// Sets the system type
/// \param system the system type
///    |code| centrality estimator, range modifier |
///    |:---|--------|
///    |  AliCSEventCuts::kNoSystem | no system dependent cut |
///    |  AliCSEventCuts::kpp  | **p-p** system |
///    |  AliCSEventCuts::kpPb | **p-Pb** system |
///    |  AliCSEventCuts::kPbPb  | **Pb-Pb** system |
///    |  AliCSEventCuts::kXeXe  | **Xe-Xe** system |
/// \return kTRUE if a proper and supported system
Bool_t AliCSEventCuts::SetSystemType(SystemType system) {
  switch(system){
  case kNoSystem:
    fSystem = kNoSystem;
    break;
  case kpp:
    fSystem = kpp;
    break;
  case kpPb:
    fSystem = kpPb;
    break;
  case kPbPb:
    fSystem = kPbPb;
    break;
  case kXeXe:
    fSystem = kXeXe;
    break;
  default:
    AliError(Form("SetSystemType not defined %d",Int_t(system)));
    return kFALSE;
  }
  /* re-evaluate centralitiy ranges in case called individualy */
  return SetCentralityType(fParameters[kCentralityType]);
}

/// Sets the acutal system type according to the data period
void AliCSEventCuts::SetActualSystemType() {
  SystemType system;
  switch (GetGlobalAnchorPeriod()) {
  case kLHC10bg:
  case kLHC11a:
  case kLHC11b:
  case kLHC11cg:
  case kLHC12:
  case kLHC13g:
  case kLHC15fm:
  case kLHC15n:
  case kLHC16k:
  case kLHC16l:
    system = kpp;
    AliInfo("SYSTEM: p-p");
    break;
  case kLHC13bc:
  case kLHC13de:
  case kLHC13f:
    system = kpPb;
    AliInfo("SYSTEM: p-Pb");
    break;
  case kLHC10h:
  case kLHC11h:
  case kLHC15oLIR:
  case kLHC15oHIR:
    system = kPbPb;
    AliInfo("SYSTEM: Pb-Pb");
    break;
  case kLHC17n:
    system = kXeXe;
    AliInfo("SYSTEM: Xe-Xe");
    break;
  default:
    system = kNoSystem;
    AliError("SYSTEM: no system set from data files");
  }

  /* re-state system type according to data period */
  SetSystemType(system);
  fParameters[kSystem] = system;
  UpdateCutsString();
}

/// Sets the type of centrality to handle
/// \param ctype the centrality type
///    |code| detector for centrality estimate, range modifier |
///    |:--:|--------|
///    |  0 | no centrality cut |
///    |  1 | default detector for the concerned system, cut in the range 0-100% in steps of 10% |
///    |  2 | alternative detector for the concerned system, cut in the range 0-100% in steps of 10% |
///    |  3 | default detector for the concerned system, cut in the range 0-50% in steps of 5% |
///    |  4 | default detector for the concerned system, cut in the range 50-100% in steps of 5% |
///    |  5 | default detector for the concerned system, cut in the range 0-10% in steps of 1% |
///    |  6 | default detector for the concerned system, cut in the range 10-20% in steps of 1% |
/// \return kTRUE if proper and supported centrality type
///
/// The default and alternative detector for centrality estimation in the different systems
///   | System | default detector | alternative detector |
///   |:--|:--:|:--:|
///   | **p-p** | **V0M** | **CL1** |
///   | **p-pB** | **V0A** | **CL1** |
///   | **pB-pB** | **V0M** | **CL1** |
///
Bool_t AliCSEventCuts::SetCentralityType(Int_t ctype)
{   // Set Cut
  switch(ctype){
  case 0:
    fCutsEnabledMask.ResetBitNumber(kCentralityCut);
    fCentralityDetector=0;
    fCentralityModifier=0;
    break;
  case 1:
    /* default centrality detector for the concerned system */
    /* centrality cut in the range 0-100% */
    fCutsEnabledMask.SetBitNumber(kCentralityCut);
    fCentralityDetector=0;
    fCentralityModifier=0;
    break;
  case 2:
    /* alternative centrality detectorfor the concerned system */
    /* centrality cut in the range 0-100% */
    fCutsEnabledMask.SetBitNumber(kCentralityCut);
    fCentralityDetector=1;
    fCentralityModifier=0;
    break;
  case 3:
    /* default centrality detector for the concerned system */
    /* centrality cut in the range 0-50% in steps of 5% */
    fCutsEnabledMask.SetBitNumber(kCentralityCut);
    fCentralityDetector=0;
    fCentralityModifier=1;
    break;
  case 4:
    /* default centrality detector for the concerned system */
    /* centrality cut in the range 50-100% in steps of 5% */
    fCutsEnabledMask.SetBitNumber(kCentralityCut);
    fCentralityDetector=0;
    fCentralityModifier=2;
    break;
  case 5:
    /* default centrality detector for the concerned system */
    /* centrality cut in the range 0-10% in steps of 1% */
    fCutsEnabledMask.SetBitNumber(kCentralityCut);
    fCentralityDetector=0;
    fCentralityModifier=3;
    break;
  case 6:
    /* default centrality detector for the concerned system */
    /* centrality cut in the range 10-20% in steps of 1% */
    fCutsEnabledMask.SetBitNumber(kCentralityCut);
    fCentralityDetector=0;
    fCentralityModifier=4;
    break;
  default:
    AliError(Form("Centrality type %d not supported",ctype));
    return kFALSE;
  }
  /* re-evaluate centrality ranges in case called individually */
  return SetCentralityMax(fParameters[kCentralityMax]);
}

/// Sets the min value for the centrality cut
/// The action is delayed until getting the max value
/// \param min the min value for the centrality cut
/// \return kTRUE always
/// If \f$max < min\f$ or \f$min = max \neq 0\f$ any positive centrality value is accepted.
///
/// For **p-p** systems \f$min\f$ and \f$max\f$ are indexes of the array
/// ~~~~{.cpp}
/// static const Float_t primaryTracksFor_pp [10] = { 0, 2, 5, 10, 15, 30, 50, 100, 500, 1000};
/// ~~~~
Bool_t AliCSEventCuts::SetCentralityMin(Int_t min)
{
  /* re-evaluate centrality ranges in case called individually */
  return SetCentralityMax(fParameters[kCentralityMax]);
}

/// Sets the max value for the centrality cut
///
/// min and max values are analyzed together with system type
/// and centrality type to extract the actual min and max
/// values for the cut.
/// \param max the max value for the centrality cut
/// \return kTRUE if the min and max values are consistent
/// If \f$max < min\f$ or \f$min = max \neq 0\f$ any positive centrality value is accepted.
///
/// For **p-p** systems \f$min\f$ and \f$max\f$ are indexes of the array
/// ~~~~{.cpp}
/// static const Float_t primaryTracksFor_pp [10] = { 0, 2, 5, 10, 15, 30, 50, 100, 500, 1000};
/// ~~~~
/// If \f$max = 0\f$ then \f$max = 10\f$
Bool_t AliCSEventCuts::SetCentralityMax(Int_t max)
{
  static const Float_t primaryTracksFor_pp [10] = { 0, 2, 5, 10, 15, 30, 50, 100, 500, 1000};

  /* first check if the cut is active */
  if (fCutsEnabledMask.TestBitNumber(kCentralityCut)) {
    /* we rescue the min value */
    Int_t min = fParameters[kCentralityMin];

    if ((max < min) ||
        ((min == max) && (min != 0))) {
      /* accept everything */
      fCentralityMin = 0.0;
      fCentralityMax = 1e6;
    }

    if(fSystem == kpp){
      fCentralityMin = primaryTracksFor_pp[min];
      fCentralityMax = primaryTracksFor_pp[max];
    }
    else {
      /* full range */
      if (max == 0) max = 10;
      switch (fCentralityModifier) {
      case 0:
        fCentralityMin = min * 10;
        fCentralityMax = max * 10;
        break;
      case 1:
        fCentralityMin = min * 5;
        fCentralityMax = max * 5;
        break;
      case 2:
        fCentralityMin = 50 + min * 5;
        fCentralityMax = 50 + max * 5;
        break;
      case 3:
        fCentralityMin = min;
        fCentralityMax = max;
        break;
      case 4:
        fCentralityMin = 10 + min;
        fCentralityMax = 10 + max;
        break;
      default:
        AliError("Inconsistent centrality modifier");
        return kFALSE;
      }
    }
    return kTRUE;
  }
  else
    /* cut not active, accept anything */
    return kTRUE;
}

/// Check for the use of the new multiplicity framework
/// \return kTRUE if the new multiplicity framework has to be use, kFALSE otherwise
Bool_t AliCSEventCuts::UseNewMultiplicityFramework() const{

  switch (GetGlobalAnchorPeriod()) {
  case kLHC15oLIR:
  case kLHC15oHIR:
  case kLHC17n:
    AliInfo("Using NEW mulitplicity framework");
    return kTRUE;
  default:
    AliInfo("Using TRADITIONAL centrality framework");
    return kFALSE;
  }
}

/// Gets the event centrality
///
/// A check on the period under analysis is done to select the way to estimate centrality.
/// Additionally the configured default or alternative detector is used.
/// \param event the current event to analyze
/// \return the event centrality

Float_t AliCSEventCuts::GetEventCentrality(AliVEvent *event) const
{
  AliESDEvent *esdEvent=dynamic_cast<AliESDEvent*>(event);
  AliAODEvent *aodEvent=dynamic_cast<AliAODEvent*>(event);

  /* for p-p systems just return FB32 accepted multiplicity */
  if (fSystem == kpp)
    return this->fNoOfFB32AccTracks;

  if (esdEvent != NULL) {
    AliCentrality *Centrality = event->GetCentrality();
    if (fUseNewMultFramework) {
      AliMultSelection *MultSelection = (AliMultSelection*) event->FindListObject("MultSelection");
      if (MultSelection == NULL) {
        AliError("No MultSelection object instance");
        return -1.0;
      }
      switch(fCentralityDetector) {
      case 0:
        /* default centrality detector */
        if (fSystem == kpPb)
          return Centrality->GetCentralityPercentile("V0A");
        else
          /* we don't leave Multiplicity task cuts precede our own */
          return MultSelection->GetMultiplicityPercentile("V0M", kFALSE);
      case 1:
        /* alternative centrality detector */
        /* we don't leave Multiplicity task cuts precede our own */
        return MultSelection->GetMultiplicityPercentile("CL1", kFALSE);
      default:
        AliError("Wrong stored centrality detector");
        return -1.0;
      }
    }
    else {
      switch(fCentralityDetector) {
      case 0:
        /* default centrality detector */
        if (fSystem == kpPb)
          return Centrality->GetCentralityPercentile("V0A");
        else
          return Centrality->GetCentralityPercentile("V0M");
      case 1:
        /* alternative centrality detector */
        return Centrality->GetCentralityPercentile("CL1");
      default:
        AliError("Wrong stored centrality detector");
        return -1.0;
      }
    }
  }

  if(aodEvent){
    if(fUseNewMultFramework){
      AliMultSelection *MultSelection = (AliMultSelection *)event->FindListObject("MultSelection");
      if (MultSelection == NULL) {
        AliError("No MultSelection object instance");
        return -1.0;
      }
      switch(fCentralityDetector) {
      case 0:
        /* default centrality detector */
        /* we don't leave Multiplicity task cuts precede our own */
        return MultSelection->GetMultiplicityPercentile("V0M",kFALSE);
      case 1:
        /* alternative centrality detector */
        /* we don't leave Multiplicity task cuts precede our own */
        return MultSelection->GetMultiplicityPercentile("CL1",kFALSE);
      default:
        AliError("Wrong stored centrality detector");
        return -1.0;
      }
    }else{
      if(aodEvent->GetHeader()) {
        AliCentrality *aodCentrality = ((AliVAODHeader*)aodEvent->GetHeader())->GetCentralityP();
        if (aodCentrality != NULL) {
          switch(fCentralityDetector) {
          case 0:
            /* default centrality detector */
            if (fSystem == kpPb)
              return aodCentrality->GetCentralityPercentile("V0A");
            else
              return aodCentrality->GetCentralityPercentile("V0M");
          case 1:
            /* alternative centrality detector */
            return aodCentrality->GetCentralityPercentile("CL1");
          default:
            AliError("Wrong stored centrality detector");
            return -1.0;
          }
        }
        else {
          AliError("No AliCentrality attached to AOD header");
          return -1;
        }
      }
      else {
        AliError("Not a standard AOD");
        return -1;
      }
    }
  }
  return -1;
}


/// Gets the event centrality from the alternative detector
///
/// A check on the period under analysis is done to select the way to estimate centrality.
/// Additionally the alternative to the configured detector is used.
/// \param event the current event to analyze
/// \return the event centrality from the alternative detector

Float_t AliCSEventCuts::GetEventAltCentrality(AliVEvent *event) const
{
  AliESDEvent *esdEvent=dynamic_cast<AliESDEvent*>(event);
  AliAODEvent *aodEvent=dynamic_cast<AliAODEvent*>(event);

  if (esdEvent != NULL) {
    AliCentrality *Centrality = event->GetCentrality();
    if (fUseNewMultFramework) {
      AliMultSelection *MultSelection = (AliMultSelection*) event->FindListObject("MultSelection");
      if (MultSelection == NULL) {
        AliError("No MultSelection object instance");
        return -1.0;
      }
      switch(fCentralityDetector) {
      case 0:
        /* default centrality detector  selected so we use the alternative one*/
        /* we don't leave Multiplicity task cuts precede our own */
        return MultSelection->GetMultiplicityPercentile("CL1", kFALSE);

      case 1:
        /* alternative centrality detector so we use the default one */
        if (fSystem == kpPb)
          return Centrality->GetCentralityPercentile("V0A");
        else
          /* we don't leave Multiplicity task cuts precede our own */
          return MultSelection->GetMultiplicityPercentile("V0M", kFALSE);
      default:
        AliError("Wrong stored centrality detector");
        return -1.0;
      }
    }
    else {
      switch(fCentralityDetector) {
      case 0:
        /* default centrality detector so we use the alternative one */
        return Centrality->GetCentralityPercentile("CL1");
      case 1:
        /* alternative centrality detector so we use the default one */
        if (fSystem == kpPb)
          return Centrality->GetCentralityPercentile("V0A");
        else
          return Centrality->GetCentralityPercentile("V0M");
      default:
        AliError("Wrong stored centrality detector");
        return -1.0;
      }
    }
  }

  if(aodEvent){
    if(fUseNewMultFramework){
      AliMultSelection *MultSelection = (AliMultSelection *)event->FindListObject("MultSelection");
      if (MultSelection == NULL) {
        AliError("No MultSelection object instance");
        return -1.0;
      }
      switch(fCentralityDetector) {
      case 0:
        /* default centrality detector so we use the alternative one */
        /* we don't leave Multiplicity task cuts precede our own */
        return MultSelection->GetMultiplicityPercentile("CL1",kFALSE);
      case 1:
        /* alternative centrality detector so we use the default one */
        /* we don't leave Multiplicity task cuts precede our own */
        return MultSelection->GetMultiplicityPercentile("V0M",kFALSE);
      default:
        AliError("Wrong stored centrality detector");
        return -1.0;
      }
    }else{
      if(aodEvent->GetHeader()) {
        AliCentrality *aodCentrality = ((AliVAODHeader*)aodEvent->GetHeader())->GetCentralityP();
        if (aodCentrality != NULL) {
          switch(fCentralityDetector) {
          case 0:
            /* default centrality detector so we use the alternative one */
            return aodCentrality->GetCentralityPercentile("CL1");
          case 1:
            /* alternative centrality detector so we use the default one */
            if (fSystem == kpPb)
              return aodCentrality->GetCentralityPercentile("V0A");
            else
              return aodCentrality->GetCentralityPercentile("V0M");
          default:
            AliError("Wrong stored centrality detector");
            return -1.0;
          }
        }
        else {
          AliError("No AliCentrality attached to AOD header");
          return -1;
        }
      }
      else {
        AliError("Not a standard AOD");
        return -1;
      }
    }
  }
  return -1;
}



/// Gets the names of the centrality estimator detectors
///
/// A check on the period under analysis is done to select the way to estimate centrality.
/// The configured and of the alternative detector is returned.
/// \param sel the name of the selected detector
/// \param alt the name of the alternative selector

void AliCSEventCuts::GetCentralityEstimatorNames(const char *&sel, const char *&alt) const
{
  /* TODO: when ESD - AOD is supported by the base clas this method needs to be adapted */
  /* so far only considers the ESD case */
  static const char V0A[] = "V0A";
  static const char V0M[] = "V0M";
  static const char CL1[] = "CL1";

  sel = NULL;
  alt = NULL;

  switch(fCentralityDetector) {
  case 0:
    /* default centrality detector */
    if (fSystem == kpPb) {
      sel = V0A;
      alt = CL1;
    }
    else {
      sel = V0M;
      alt = CL1;
    }
    break;
  case 1:
    /* alternative centrality detector */
    if (fSystem == kpPb) {
      sel = CL1;
      alt = V0A;
    }
    else {
      sel = CL1;
      alt = V0M;
    }
    break;
  default:
    AliError("Wrong stored centrality detector");
  }
}



/// Manual trigger selection
/// \param trigger the selected trigger combination to activate
///    |code| trigger | Observations |
///    |:--:|--------|
///    |  0 | AliVEvent::kAny | |
///    |  1 | AliVEvent::kMB or AliVEvent::kINT7 | depends on period |
///    |  5 | AliVEvent::kDG or AliVEvent::kDG5 | double gap diffractive |
/// \return kTRUE for proper and supported active trigger selection
Bool_t AliCSEventCuts::SetSelectActiveTrigger(Int_t trigger)
{
  switch(trigger){
  case 0: /* Any */
    fCutsEnabledMask.SetBitNumber(kOfflineTriggerCut);
    break;
  case 1: /* MB */
    fCutsEnabledMask.SetBitNumber(kOfflineTriggerCut);
    break;
  case 5: /* Double gap diffractive */
    fCutsEnabledMask.SetBitNumber(kOfflineTriggerCut);
    fOfflineTriggerMask = AliVEvent::kDG5;
    break;
  default:
    AliError(Form("Active trigger code %d not supported", trigger));
    return kFALSE;
  }
  return kTRUE;
}

/// Manual trigger setting
///
/// Sets the actual active trigger according to the configurated value
/// and the current data period under analysis.
void AliCSEventCuts::SetActualActiveTrigger()
{
  switch(fParameters[kActiveTrigger]){
  case 0: /* Any */
    fOfflineTriggerMask = AliVEvent::kAny;
    break;
  case 1: /* MB */
    switch (GetGlobalAnchorPeriod()) {
    case kLHC15oLIR:
    case kLHC15oHIR:
    case kLHC16k:
    case kLHC16l:
    case kLHC17n:
      fOfflineTriggerMask = AliVEvent::kINT7;
      AliInfo("Using AliVEvent::kINT7 as MB trigger");
      break;
    default:
      fOfflineTriggerMask = AliVEvent::kMB;
      AliInfo("Using AliVEvent::kMB as MB trigger");
    }
    break;
  case 5: /* Double gap diffractive, already set */
      break;
  default:
    AliError(Form("Active trigger code %d not supported", fParameters[kActiveTrigger]));
  }
}

/// Manual sub-trigger selection
/// \param trigger the selected sub-trigger combination to activate
///    |code| sub-trigger |
///    |:--:|--------|
///    |  0 | none |
/// \return kTRUE for proper and supported active sub-trigger selection
Bool_t AliCSEventCuts::SetSelectActiveSubTrigger(Int_t trigger)
{
  switch(trigger){
  case 0: /* none */
    /* the only value so far supported: do nothing */
    break;
  default:
    AliError(Form("Active sub-trigger code %d not supported", trigger));
    return kFALSE;
  }
  return kTRUE;
}

/// Sets and configures the procedure to remove event pile up
/// \param pupcode the pileup removal cut code
///    |code| method |
///    |:--:|--------|
///    |  0 | no pileup rejection |
///    |  1 | PileUpSPD & background rejection with cluster-vs-tracklet default cut |
///    |  2 | PileUpSPD & background rejection with cluster-vs-tracklet default cut & out-of-bunch pileup rejection using trigger information |
///    |  3 | PileUpMV & background rejection with cluster-vs-tracklet default cut |
///    |  4 | PileUpMV & background rejection with cluster-vs-tracklet default cut & out-of-bunch pileup rejection using trigger information |
/// \return kTRUE for proper and supported pile up removal procedures
Bool_t AliCSEventCuts::SetRemovePileUp(Int_t pupcode)
{
  switch(pupcode) {
  case 0:
    fCutsEnabledMask.ResetBitNumber(kPileUpCut);
    break;
  case 1: /* use SPD multiple vertices */
    fCutsEnabledMask.SetBitNumber(kPileUpCut);
    break;
  case 2: /* use SPD multiple vertices  & out-of-bunch pileup rejection */
    fCutsEnabledMask.SetBitNumber(kPileUpCut);
    fAnalysisUtils.SetUseOutOfBunchPileUp(kTRUE);
    break;
  case 3: /* use track multivertexer */
    fCutsEnabledMask.SetBitNumber(kPileUpCut);
    fAnalysisUtils.SetUseMVPlpSelection(kTRUE);
    break;
  case 4: /* use track multivertexer & out-of-bunch pileup rejection */
    fCutsEnabledMask.SetBitNumber(kPileUpCut);
    fAnalysisUtils.SetUseMVPlpSelection(kTRUE);
    fAnalysisUtils.SetUseOutOfBunchPileUp(kTRUE);
    break;
  default:
    AliError(Form("PilEup removal procedure %d not supported", pupcode));
    return kFALSE;
  }
  return kTRUE;
}

/// Sets the z vertex cut
/// \param vtxcut the z vertex cut code
///    |code| cut                        | additional functionality                                                            |
///    |:--:|----------------------------|-------------------------------------------------------------------------------------|
///    |  0 | no z vertex cut            |                                                                                     |
///    |  1 | \f$ |z_{vtx}|\f$ < 12.0 cm |                                                                                     |
///    |  2 | \f$ |z_{vtx}|\f$ < 10.0 cm |                                                                                     |
///    |  3 | \f$ |z_{vtx}|\f$ < 7.0 cm  |                                                                                     |
///    |  4 | \f$ |z_{vtx}|\f$ < 5.0 cm  |                                                                                     |
///    |  5 | \f$ |z_{vtx}|\f$ < 12.0 cm | distance between SPD and tracks vertex less than 0.5 (PbPb 2011) or 0.2 cm (PbPb 2015) |
///    |  6 | \f$ |z_{vtx}|\f$ < 10.0 cm | distance between SPD and tracks vertex less than 0.5 (PbPb 2011) or 0.2 cm (PbPb 2015) |
///    |  7 | \f$ |z_{vtx}|\f$ < 7.0 cm  | distance between SPD and tracks vertex less than 0.5 (PbPb 2011) or 0.2 cm (PbPb 2015) |
///    |  8 | \f$ |z_{vtx}|\f$ < 5.0 cm  | distance between SPD and tracks vertex less than 0.5 (PbPb 2011) or 0.2 cm (PbPb 2015) |
/// \return kTRUE for proper and supported vertex cuts

Bool_t AliCSEventCuts::SetVertexCut(Int_t vtxcut) {

  switch(vtxcut){
  case 0:
    fCutsEnabledMask.ResetBitNumber(kVertexCut);
    fMaxVertexZ     = 1000;
    break;
  case 1:
    fCutsEnabledMask.SetBitNumber(kVertexCut);
    fMaxVertexZ     = 12.0;
    break;
  case 2:
    fCutsEnabledMask.SetBitNumber(kVertexCut);
    fMaxVertexZ     = 10.0;
    break;
  case 3:
    fCutsEnabledMask.SetBitNumber(kVertexCut);
    fMaxVertexZ     = 7.0;
    break;
  case 4:
    fCutsEnabledMask.SetBitNumber(kVertexCut);
    fMaxVertexZ     = 5.;
    break;
  case 5:
    fCutsEnabledMask.SetBitNumber(kVertexCut);
    fMaxVertexZ     = 12.0;
    fUseSPDTracksVtxDist = kTRUE;
    break;
  case 6:
    fCutsEnabledMask.SetBitNumber(kVertexCut);
    fMaxVertexZ     = 10.0;
    fUseSPDTracksVtxDist = kTRUE;
    break;
  case 7:
    fCutsEnabledMask.SetBitNumber(kVertexCut);
    fMaxVertexZ     = 7.0;
    fUseSPDTracksVtxDist = kTRUE;
    break;
  case 8:
    fCutsEnabledMask.SetBitNumber(kVertexCut);
    fMaxVertexZ     = 5.;
    fUseSPDTracksVtxDist = kTRUE;
    break;
  default:
    AliError(Form("Vertex Cut %d not supported",vtxcut));
    return kFALSE;
  }
  return kTRUE;
}

/// Gets the number of contributors for the vertex estimation
/// \param event the current event to analyze
/// \return the number of vertex contributors
Int_t AliCSEventCuts::GetNumberOfVertexContributors(AliVEvent *event) const{

  const AliVVertex *vtx = event->GetPrimaryVertex();

  if (vtx != NULL){
    if(vtx->GetNContributors()>1) {
      return vtx->GetNContributors();
    }
  }

  const AliVVertex *vtxSPD = event->GetPrimaryVertexSPD();

  if(vtxSPD !=NULL)
    return vtxSPD->GetNContributors();
  return 0;
}

/// Compare the vertex resolution and dispersion with the stored limits for both of them
/// \param event the current event to analyze
/// \return kTRUE if values within limits, kFALSE otherwise
Bool_t AliCSEventCuts::PassVertexResolutionAndDispersionTh(AliVEvent *event) const {

  const AliVVertex *vtxSPD = event->GetPrimaryVertexSPD();
  if (vtxSPD != NULL) {
    Double_t cov[6] = {0};
    vtxSPD->GetCovarianceMatrix(cov);
    Double_t zRes = TMath::Sqrt(cov[5]);

    /* check vertex resolution */
    if (vtxSPD->IsFromVertexerZ() && (fVertexResolutionTh < zRes))
      /* bad resolution vertex */
      return kFALSE;

    /* if ESD check vertex dispersion */
    AliESDEvent *esd = dynamic_cast<AliESDEvent*>(event);
    if(esd != NULL){
      const AliESDVertex *esdVtxSPD = dynamic_cast<const AliESDVertex *>(vtxSPD);
      if (esdVtxSPD->IsFromVertexerZ() && (fVertexDispersionTh < esdVtxSPD->GetDispersion()))
        /* poor resolution vertex */
        return kFALSE;
    }
    return kTRUE;
  }
  else
    return kFALSE;
}

/// Check if the distance between SPD and tracks vertex are acceptable by the standar cut
/// \param event the current event to analyze
/// \return kTRUE if the vertex distance is accepted kFALSE otherwise
Bool_t AliCSEventCuts::AcceptSPDTracksVtxDist(AliVEvent *event) const {

  const AliVVertex *vtx = event->GetPrimaryVertex();

  if (vtx != NULL){
    if(vtx->GetNContributors() < 2) {
      return kFALSE;
    }
  }
  else {
    return kFALSE;
  }

  const AliVVertex *vtxSPD = event->GetPrimaryVertexSPD();

  if(vtxSPD !=NULL) {
    if (vtxSPD->GetNContributors() < 1) {
      return kFALSE;
    }
  }
  else {
    return kFALSE;
  }

  /* so, proper contributors and proper vertex */
  Double_t covTrck[6], covSPD[6];
  vtx->GetCovarianceMatrix(covTrck);
  vtxSPD->GetCovarianceMatrix(covSPD);

  Double_t dz = vtx->GetZ() - vtxSPD->GetZ();

  Double_t errTot = TMath::Sqrt(covTrck[5]+covSPD[5]);
  Double_t errTrck = TMath::Sqrt(covTrck[5]);
  Double_t nsigTot = dz/errTot;
  Double_t nsigTrck = dz/errTrck;

  if (TMath::Abs(dz) > fgkSPDTracksVtxDistanceThreshold || TMath::Abs(nsigTot) > fgkSPDTracksVtxDistanceSigmas || TMath::Abs(nsigTrck) > fgkTrackVertexSigmas)
     return kFALSE;

  return kTRUE;
}

/// Gets the the event z vertex coordinate
/// \param event the current event to analyze
/// \return the event z vertex coordinate
Double_t AliCSEventCuts::GetVertexZ(AliVEvent *event) const {
  Double_t vertexZ = 999.0;
  const AliVVertex *vtx = event->GetPrimaryVertex();
  const AliVVertex *vtxSPD = event->GetPrimaryVertexSPD();

  if (vtx != NULL){
    vertexZ = vtx->GetZ();
  }
  else {
    if (vtxSPD != NULL) {
      vertexZ = vtxSPD->GetZ();
    }
  }
  return vertexZ;
}

/// Sets and configures the procedure to remove 2015 additional event pileup
/// \param pupcode the 2015 additional pileup removal cut code
///    |code| method |
///    |:--:|--------|
///    |  0 | no 2015 additional pileup rejection |
///    |  1 | J/psi analysis pileup removal, initial (faulty) track counting method |
///    |  2 | Centrality and multiplicity correlations for 2015 |
///    |  3 | J/psi analysis pileup removal, initial (corrected) track counting method |
/// \return kTRUE for proper and supported 2015 additional pileup removal procedures
Bool_t AliCSEventCuts::SetRemove2015PileUp(Int_t pupcode)
{
  switch(pupcode) {
  case 0:
    fCutsEnabledMask.ResetBitNumber(k2015PileUpCut);
    break;
  case 1: /* J/psi analysis pileup removal method */
    fCutsEnabledMask.SetBitNumber(k2015PileUpCut);
    break;
  case 2: /* Centrality and multiplicity correlations for 2015*/
    fCutsEnabledMask.SetBitNumber(k2015PileUpCut);
    break;
  case 3: /* J/psi analysis pileup removal initial method */
    fCutsEnabledMask.SetBitNumber(k2015PileUpCut);
    break;
  default:
    AliError(Form("2015 additional pileup removal procedure %d not supported", pupcode));
    return kFALSE;
  }
  return kTRUE;
}

/// Sets the actual 2015 pileup removal according to the configured value
/// and the current data period under analysis.
void AliCSEventCuts::SetActual2015PileUpRemoval()
{
  switch(fParameters[kRemove2015PileUp]){
  case 0: /* no additional pileup rejection */
    break;
  case 1: /* J/psi analysis pileup removal method */
    if(f2015V0MtoTrkTPCout){
      delete f2015V0MtoTrkTPCout;
    }
    switch (GetGlobalAnchorPeriod()) {
    case kLHC10h:
      f2015V0MtoTrkTPCout = new TFormula(Form("f2015V0MtoTrkTPCout_%s",GetCutsString()),"-1000+2.8*x"); /* pass2 with the initial, faulty, method for track count */
      break;
    case kLHC15oLIR:
      /* f2015V0MtoTrkTPCout = new TFormula(Form("f2015V0MtoTrkTPCout_%s",GetCutsString()),"-4000+3.8*x"); pass2 */
      f2015V0MtoTrkTPCout = new TFormula(Form("f2015V0MtoTrkTPCout_%s",GetCutsString()),"-800+2.93*x"); /* pass3 */
      break;
    case kLHC15oHIR:
      f2015V0MtoTrkTPCout = new TFormula(Form("f2015V0MtoTrkTPCout_%s",GetCutsString()),"-2000+3.0*x");
      break;
    default:
      f2015V0MtoTrkTPCout = new TFormula(Form("f2015V0MtoTrkTPCout_%s",GetCutsString()),"-450+10.5*x");
      break;
    }
    AliInfo(Form("2015 pileup removal: V0 mult < %s\n", TString(f2015V0MtoTrkTPCout->GetTitle()).ReplaceAll("x","trkTPCout").Data()));
    break;
  case 2: /* Centrality and multiplicity correlations for 2015*/
    if (fCentOutLowCut != NULL)
      delete fCentOutLowCut;
    if (fCentOutHighCut != NULL)
      delete fCentOutHighCut;
    if (fTOFMultOutLowCut != NULL)
      delete fTOFMultOutLowCut;
    if (fTOFMultOutHighCut != NULL)
      delete fTOFMultOutHighCut;
    if (fMultCentOutLowCut != NULL)
      delete fMultCentOutLowCut;
    switch (GetGlobalAnchorPeriod()) {
    case kLHC10h:
      /* inhibited, TODO */
      fCutsEnabledMask.ResetBitNumber(k2015PileUpCut);
      AliWarning("2015 pileup removal: inhibited for LHC10h anchored datasets");
      break;
    case kLHC15oLIR:
      /* inhibited, TODO */
      fCutsEnabledMask.ResetBitNumber(k2015PileUpCut);
      AliWarning("2015 pileup removal: inhibited for LHC15oLIR anchored datasets");
      break;
    case kLHC15oHIR:
      f2015V0MtoTrkTPCout = new TFormula(Form("f2015V0MtoTrkTPCout_%s",GetCutsString()),"-2000+3.0*x");
      /* let's initialize the expressions for 2015 pile up rejection */
      fCentOutLowCut = new TF1("fCentOutLowCut", "[0]+[1]*x - 5.*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
      fCentOutLowCut->SetParameters(0.0157497, 0.973488, 0.673612, 0.0290718, -0.000546728, 5.82749e-06);
      fCentOutHighCut = new TF1("fCentOutHighCut", "[0]+[1]*x + 5.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
      fCentOutHighCut->SetParameters(0.0157497, 0.973488, 0.673612, 0.0290718, -0.000546728, 5.82749e-06);
      fTOFMultOutLowCut = new TF1("fTOFMultOutLowCut", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 4.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x+[9]*x*x*x*x*x)", 0, 10000);
      fTOFMultOutLowCut->SetParameters(-1.0178, 0.333132, 9.10282e-05, -1.61861e-08, 1.47848, 0.0385923, -5.06153e-05, 4.37641e-08, -1.69082e-11, 2.35085e-15);
      fTOFMultOutHighCut = new TF1("fTOFMultOutHighCut", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 4.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x+[9]*x*x*x*x*x)", 0, 10000);
      fTOFMultOutHighCut->SetParameters(-1.0178, 0.333132, 9.10282e-05, -1.61861e-08, 1.47848, 0.0385923, -5.06153e-05, 4.37641e-08, -1.69082e-11, 2.35085e-15);
      fMultCentOutLowCut = new TF1("fMultCentOutLowCut", "[0]+[1]*x+[2]*exp([3]-[4]*x) - 5.*([5]+[6]*exp([7]-[8]*x))", 0, 100);
      fMultCentOutLowCut->SetParameters(-6.15980e+02, 4.89828e+00, 4.84776e+03, -5.22988e-01, 3.04363e-02, -1.21144e+01, 2.95321e+02, -9.20062e-01, 2.17372e-02);

      /* TODO user feedback */
      // AliInfo(Form("2015 pileup removal: V0 mult < %s\n", TString(f2015V0MtoTrkTPCout->GetTitle()).ReplaceAll("x","trkTPCout").Data()));
      break;
    default:
      /* inhibited, TODO */
      fCutsEnabledMask.ResetBitNumber(k2015PileUpCut);
      AliWarning("2015 pileup removal: inhibited for unknown anchored datasets");
      break;
    }
    break;
  case 3: /* J/psi analysis pileup removal initial, corrected, method */
    if(f2015V0MtoTrkTPCout){
      delete f2015V0MtoTrkTPCout;
    }
    switch (GetGlobalAnchorPeriod()) {
    case kLHC10h:
      f2015V0MtoTrkTPCout = new TFormula(Form("f2015V0MtoTrkTPCout_%s",GetCutsString()),"-1000+3.1*x");
      break;
    case kLHC15oLIR:
      /* f2015V0MtoTrkTPCout = new TFormula(Form("f2015V0MtoTrkTPCout_%s",GetCutsString()),"-4000+3.8*x"); pass2 */
      f2015V0MtoTrkTPCout = new TFormula(Form("f2015V0MtoTrkTPCout_%s",GetCutsString()),"-800+2.93*x"); /* pass3 */
      break;
    case kLHC15oHIR:
      f2015V0MtoTrkTPCout = new TFormula(Form("f2015V0MtoTrkTPCout_%s",GetCutsString()),"-2500+5.0*x");
      break;
    case kLHC17n:
      f2015V0MtoTrkTPCout = new TFormula(Form("f2015V0MtoTrkTPCout_%s",GetCutsString()),"-900+6.0*x");
      break;
    default:
      f2015V0MtoTrkTPCout = new TFormula(Form("f2015V0MtoTrkTPCout_%s",GetCutsString()),"-1000+2.8*x");
      break;
    }
    AliInfo(Form("2015 pileup removal (initial method): V0 mult < %s\n", TString(f2015V0MtoTrkTPCout->GetTitle()).ReplaceAll("x","trkTPCout").Data()));
    break;
  default:
    AliError(Form("2015 additional pileup removal code %d not supported", fParameters[kRemove2015PileUp]));
  }
}

/// Checks whether the V0 multiplicity and the number of TPCout on tracks mark the event as pileup
/// \return kTRUE if the event is a pileup event kFALSE otherwise
Bool_t AliCSEventCuts::Is2015PileUpEvent() const {

  switch(fParameters[kRemove2015PileUp]){
  case 1: /* J/psi analysis pileup removal initial, faulty, method */
    if (fV0Multiplicity  < f2015V0MtoTrkTPCout->Eval(fNoOfInitialTPCoutTracks))
      return kTRUE;
    return kFALSE;
    break;
  case 2: /* Centrality and multiplicity correlations for 2015*/ {
      Float_t multDiffESDTPC = fNoOfESDTracks - fNoOfFB128Tracks*3.8;

      if (Float_t(fCL0Centrality) < fCentOutLowCut->Eval(fV0MCentrality))
        return kTRUE;

      if (Float_t(fCL0Centrality) > fCentOutHighCut->Eval(fV0MCentrality))
        return kTRUE;

      if (multDiffESDTPC > 15000) //for stronger cuts use 700
        return kTRUE;

      if (Float_t(fNoOfFB32TOFTracks) < fTOFMultOutLowCut->Eval(Float_t(fNoOfFB32Tracks)))
        return kTRUE;

      if (Float_t(fNoOfFB32TOFTracks) > fTOFMultOutHighCut->Eval(Float_t(fNoOfFB32Tracks)))
        return kTRUE;

      if (Float_t(fNoOfFB32AccTracks) < fMultCentOutLowCut->Eval(fV0MCentrality))
        return kTRUE;

      return kFALSE;
    }
    break;
  case 3: /* J/psi analysis pileup removal initial, corrected, method */
    if (fV0Multiplicity  < f2015V0MtoTrkTPCout->Eval(fNoOfInitialTPCoutTracks))
      return kTRUE;
    return kFALSE;
    break;
  default:
    AliFatal(Form("Inconsistent parameter value %d for removal 2015 pileup", fParameters[kRemove2015PileUp]));
    return kFALSE;
  }
}

/// Stores the event different centralities
/// \param event the current event to handle
/// \return kTRUE if the process proceeded properly kFALSE otherwise
Bool_t AliCSEventCuts::StoreEventCentralities(AliVEvent *event) {

  fV0MCentrality = -1;
  fV0ACentrality = -1;
  fV0CCentrality = -1;
  fCL0Centrality = -1;
  fCL1Centrality = -1;

  AliESDEvent *esdEvent=dynamic_cast<AliESDEvent*>(event);
  AliAODEvent *aodEvent=dynamic_cast<AliAODEvent*>(event);

  if (esdEvent != NULL) {
    if (fUseNewMultFramework) {
      AliMultSelection *MultSelection = (AliMultSelection*) event->FindListObject("MultSelection");
      if (MultSelection != NULL) {
        fV0ACentrality = MultSelection->GetMultiplicityPercentile("V0A");
        fV0CCentrality = MultSelection->GetMultiplicityPercentile("V0M");
        fV0MCentrality = MultSelection->GetMultiplicityPercentile("V0M");
        fCL0Centrality = MultSelection->GetMultiplicityPercentile("CL1");
        fCL1Centrality = MultSelection->GetMultiplicityPercentile("CL1");
      }
      else {
        AliError("No MultSelection object instance");
        return kFALSE;
      }
    }
    else {
      AliCentrality *Centrality = event->GetCentrality();
      if (Centrality != NULL) {
        fV0ACentrality = Centrality->GetCentralityPercentile("V0A");
        fV0CCentrality = Centrality->GetCentralityPercentile("V0C");
        fV0MCentrality = Centrality->GetCentralityPercentile("V0M");
        fCL0Centrality = Centrality->GetCentralityPercentile("CL0");
        fCL1Centrality = Centrality->GetCentralityPercentile("CL1");
      }
      else {
        AliError("No Centrality object instance");
        return kFALSE;
      }
    }
  }

  if(aodEvent){
    if(fUseNewMultFramework){
      AliMultSelection *MultSelection = (AliMultSelection *)event->FindListObject("MultSelection");
      if (MultSelection == NULL) {
        AliError("No MultSelection object instance");
        return kFALSE;
      }
      fV0ACentrality = MultSelection->GetMultiplicityPercentile("V0A");
      fV0CCentrality = MultSelection->GetMultiplicityPercentile("V0M");
      fV0MCentrality = MultSelection->GetMultiplicityPercentile("V0M");
      fCL0Centrality = MultSelection->GetMultiplicityPercentile("CL1");
      fCL1Centrality = MultSelection->GetMultiplicityPercentile("CL1");
    }
    else{
      if(aodEvent->GetHeader()) {
        AliCentrality *aodCentrality = ((AliVAODHeader*)aodEvent->GetHeader())->GetCentralityP();
        if (aodCentrality != NULL) {
          fV0ACentrality = aodCentrality->GetCentralityPercentile("V0A");
          fV0CCentrality = aodCentrality->GetCentralityPercentile("V0C");
          fV0MCentrality = aodCentrality->GetCentralityPercentile("V0M");
          fCL0Centrality = aodCentrality->GetCentralityPercentile("CL0");
          fCL1Centrality = aodCentrality->GetCentralityPercentile("CL1");
        }
        else {
          AliError("No AliCentrality attached to AOD header");
          return kFALSE;
        }
      }
      else {
        AliError("Not a standard AOD");
        return kFALSE;
      }
    }
  }

  AliInfo(Form("Event centralities: V0A: %.2f, V0C: %.2f, V0M: %.2f, CL0: %.2f, CL1: %.2f",
      Float_t(fV0MCentrality),
      Float_t(fV0ACentrality),
      Float_t(fV0CCentrality),
      Float_t(fCL0Centrality),
      Float_t(fCL1Centrality)));

  return kTRUE;
}

/// Set the standard filtering track cuts set according to the data period
void AliCSEventCuts::SetActualFilterTracksCuts() {

  if (fESDFB32 != NULL)
    delete fESDFB32;
  if (fESDFB128 != NULL)
    delete fESDFB128;

  TString system = "";
  TString period = "";
  TString basename = "";
  baseSystems baseSystem = kUnknownBase;


  switch (GetGlobalAnchorPeriod()) {
  case kLHC10bg:
    baseSystem = k2010based;
    basename = "2010";
    system = "p-p";
    period = "2010bg";
    break;
  case kLHC11a:
  case kLHC11b:
  case kLHC11cg:
    baseSystem = k2011based;
    basename = "2011";
    system = "p-p";
    period = "2011ag";
    break;
  case kLHC12:
  case kLHC13g:
  case kLHC15fm:
  case kLHC15n:
  case kLHC16k:
  case kLHC16l:
    baseSystem = k2011based;
    basename = "2011";
    system = "p-p";
    period = "various";
    break;
  case kLHC13bc:
  case kLHC13de:
  case kLHC13f:
    baseSystem = k2011based;
    basename = "2011";
    system = "p-Pb";
    period = "2013bf";
    break;
  case kLHC10h:
    baseSystem = k2010based;
    basename = "2010";
    system = "Pb-Pb";
    period = "2010h";
    break;
  case kLHC11h:
    baseSystem = k2011based;
    basename = "2011";
    system = "Pb-Pb";
    period = "2011h";
    break;
  case kLHC15oLIR:
  case kLHC15oHIR:
    baseSystem = k2011based;
    basename = "2011";
    system = "Pb-Pb";
    period = "2015o";
    break;
  case kLHC17n:
    baseSystem = k2011based;
    basename = "2011";
    system = "Xe-Xe";
    period = "2017n";
    break;
  default:
    fESDFB32 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010();
    fESDFB128 = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    fESDFB32->SetNameTitle(Form("FB32TrackCount_%s",GetCutsString()),"Filter tracks 2010: globals with tight DCA");
    fESDFB128->SetNameTitle(Form("FB128TrackCount_%s",GetCutsString()),"Filter tracks 2010: Standard TPC only");
    AliError("SYSTEM: no system set from data files. Standard 2010 track cuts");
    return;
  }

  switch (baseSystem) {
  case k2010based:
    fESDFB32 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010();
    fESDFB128 = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    fESDFB32->SetNameTitle(Form("FB32TrackCount_%s",GetCutsString()),"Filter tracks 2010: globals with tight DCA");
    fESDFB128->SetNameTitle(Form("FB128TrackCount_%s",GetCutsString()),"Filter tracks 2010: Standard TPC only");
    break;
  case k2011based:
    fESDFB32 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
    fESDFB128 = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    fESDFB32->SetNameTitle(Form("FB32TrackCount_%s",GetCutsString()),"Filter tracks 2011: globals with tight DCA");
    fESDFB128->SetNameTitle(Form("FB128TrackCount_%s",GetCutsString()),"Filter tracks 2011: Standard TPC only");
    break;
  default:
    AliFatal("Internal inconsistency between data period and base period for track counting filter cuts selection. ABORTING!!!");
    return;
  }

  /* report the selected cut */
  AliInfo("=============== Track count filter cut ===========================");
  AliInfo(Form("SYSTEM: %s; PERIOD: %s", system.Data(), period.Data()));
  AliInfo(Form("BASED ON: %s standards cuts", basename.Data()));
  AliInfo("=============== Track count filter cut end =======================");
}


/// Stores the event different multiplicities
/// \param event the current event to handle
/// \return kTRUE if the process proceeded properly kFALSE otherwise
Bool_t AliCSEventCuts::StoreEventMultiplicities(AliVEvent *event) {

  AliESDEvent *esdEvent=dynamic_cast<AliESDEvent*>(event);
  AliAODEvent *aodEvent=dynamic_cast<AliAODEvent*>(event);

  fNoOfAODTracks = 0;
  fNoOfESDTracks = 0;
  fNoOfFB32Tracks = 0;
  fNoOfFB128Tracks = 0;
  fNoOfFB32AccTracks = 0;
  fNoOfFB32TOFTracks = 0;
  fNoOfTPCoutTracks = 0;
  fNoOfInitialTPCoutTracks = 0;
  fReferenceMultiplicity = -1;

  Int_t nTracks = 0;

  fV0Multiplicity = event->GetVZEROData()->GetMTotV0A()+event->GetVZEROData()->GetMTotV0C();

  AliInfo(Form("Event V0M multiplicity: %d", fV0Multiplicity));
  if (aodEvent == NULL && esdEvent == NULL) {
    AliError("Not a proper event");
    return kFALSE;
  }

  if (aodEvent != NULL) {
    if(aodEvent->GetHeader()) {
      fReferenceMultiplicity = ((AliAODHeader*)aodEvent->GetHeader())->GetRefMultiplicityComb08();
      fNoOfAODTracks = aodEvent->GetNumberOfTracks();
      fNoOfESDTracks = ((AliVAODHeader*)aodEvent->GetHeader())->GetNumberOfESDTracks();
      nTracks = fNoOfAODTracks;
    }
    else {
      AliError("Not a proper AOD event. No header!");
    }
  }

  if (esdEvent != NULL) {
    AliESDtrackCuts::MultEstTrackType estType = esdEvent->GetPrimaryVertexTracks()->GetStatus() ? AliESDtrackCuts::kTrackletsITSTPC : AliESDtrackCuts::kTracklets;
    fReferenceMultiplicity = AliESDtrackCuts::GetReferenceMultiplicity(esdEvent,estType,0.8);
    fNoOfESDTracks = esdEvent->GetNumberOfTracks();
    nTracks = fNoOfESDTracks;
  }

  for (Int_t itrk = 0; itrk < nTracks; itrk++) {
    AliVTrack *trk = dynamic_cast<AliVTrack*>(event->GetTrack(itrk));

    if (trk != NULL) {

      /* the initial method of counting TPC out tracks in its two versions */
      if (!(trk->Pt() < 0.15) && (TMath::Abs(trk->Eta()) < 0.8)) {
        if (fParameters[kRemove2015PileUp] == 1) {
          /* the initial method of counting TPC out tracks, faulty */
          /* we have to force the parenthesis for silencing compiler warning */
          /* but this is how it really looks like, that's why is faulty */
          if (!(trk->GetStatus() & (AliVTrack::kTPCout != AliVTrack::kTPCout))) {
            fNoOfInitialTPCoutTracks++;
          }
        }
        else {
          /* the initial method of counting TPC out tracks, corrected */
          if ((trk->GetStatus() & AliVTrack::kTPCout) == AliVTrack::kTPCout) {
            fNoOfInitialTPCoutTracks++;
          }
        }
      }

      if (trk->IsA() == AliESDtrack::Class()) {
        if (fESDFB32->AcceptTrack(dynamic_cast<AliESDtrack *>(trk))) {
          fNoOfFB32Tracks++;

          if (TMath::Abs(trk->GetTOFsignalDz()) <= 10. && trk->GetTOFsignal() >= 12000. && trk->GetTOFsignal() <= 25000.)
            fNoOfFB32TOFTracks++;

          if ((TMath::Abs(trk->Eta()) < 0.8) && (trk->GetTPCNcls() >= 70) && (trk->Pt() >= 0.2) && (trk->Pt() < 50.)) {
            fNoOfFB32AccTracks++;
            if ((trk->GetStatus() & AliVTrack::kTPCout) == AliVTrack::kTPCout )
              fNoOfTPCoutTracks++;
          }
        }
        if (fESDFB128->AcceptTrack(dynamic_cast<AliESDtrack *>(trk))) {
          /* TODO we need to enrich this to really match AOD FB128 */
          /* there are TPC only tracks which will not become AOD FB128 tracks */
          fNoOfFB128Tracks++;
        }
      }
      else if (trk->IsA() == AliAODTrack::Class()) {
        AliAODTrack *aodt = dynamic_cast<AliAODTrack*>(trk);

        if (aodt->TestFilterBit(32)) {
          fNoOfFB32Tracks++;

          if (TMath::Abs(trk->GetTOFsignalDz()) <= 10. && trk->GetTOFsignal() >= 12000. && trk->GetTOFsignal() <= 25000.)
            fNoOfFB32TOFTracks++;

          if ((TMath::Abs(trk->Eta()) < 0.8) && (trk->GetTPCNcls() >= 70) && (trk->Pt() >= 0.2) && (trk->Pt() < 50.)) {
            fNoOfFB32AccTracks++;
            if ((trk->GetStatus() & AliVTrack::kTPCout) == AliVTrack::kTPCout )
              fNoOfTPCoutTracks++;
          }
        }
        if (aodt->TestFilterBit(128)) {
          fNoOfFB128Tracks++;
        }
      }
      else
        continue;
    }
  }

  AliInfo(Form("Event multiplicities: AOD: %d, ESD: %d, FB32: %d, FB128: %d, FB32 acc: %d", fNoOfAODTracks, fNoOfESDTracks, fNoOfFB32Tracks, fNoOfFB128Tracks, fNoOfFB32AccTracks));
  AliInfo(Form("Event multiplicities: FB32 TOF: %d, TPC out: %d, TPC out(initial): %d, Ref: %d", fNoOfFB32TOFTracks, fNoOfTPCoutTracks, fNoOfInitialTPCoutTracks, fReferenceMultiplicity));

  return kTRUE;
}

///
/// Initializes the cuts
///
/// Initializes the needed data and allocates the needed histograms list if needed
/// \param name an additional name to precede the cuts string
void AliCSEventCuts::InitCuts(const char *name){

  if (name == NULL) name = GetName();

  fAnalysisUtils.SetUseSPDCutInMultBins(kTRUE);

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
void AliCSEventCuts::DefineHistograms(){

  if (fQALevel > kQALevelNone) {
    Bool_t oldstatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);

    fHistogramsList->SetName(Form("%s_%s",fHistogramsList->GetName(),GetCutsString()));

    fhCutsStatistics = new TH1F(Form("CutsStatistics_%s",GetCutsString()),"Event cuts statistics",kNCuts+4,-0.5,kNCuts+3.5);
    fhUniqueCutsStatistics = new TH1F(Form("UniqueCutsStatistics_%s",GetCutsString()),"Unique event cuts statistics",kNCuts+4,-0.5,kNCuts+3.5);
    fhCutsStatistics->GetXaxis()->SetBinLabel(1,"n events");
    fhUniqueCutsStatistics->GetXaxis()->SetBinLabel(1,"n events");
    fhCutsStatistics->GetXaxis()->SetBinLabel(2,"n cut events");
    fhUniqueCutsStatistics->GetXaxis()->SetBinLabel(2,"n passed events");
    for (Int_t i = 0; i < kNCuts; i++) {
      if (fgIsMC) {
        fhCutsStatistics->GetXaxis()->SetBinLabel(i+4, fgkCutsNames[i]);
        fhUniqueCutsStatistics->GetXaxis()->SetBinLabel(i+4, fgkCutsNames[i]);
      }
      else {
        if (i != kMCdataQuality) {
          fhCutsStatistics->GetXaxis()->SetBinLabel(i+4, fgkCutsNames[i]);
          fhUniqueCutsStatistics->GetXaxis()->SetBinLabel(i+4, fgkCutsNames[i]);
        }
        else {
          fhCutsStatistics->GetXaxis()->SetBinLabel(i+4, "n/a");
          fhUniqueCutsStatistics->GetXaxis()->SetBinLabel(i+4, "n/a");
        }
      }
    }
    fHistogramsList->Add(fhCutsStatistics);
    fHistogramsList->Add(fhUniqueCutsStatistics);

    if(fQALevel == kQALevelHeavy){
      fhCutsCorrelation = new TH2F(Form("CutCorrelation_%s",GetCutsString()),"Cuts correlation",kNCuts+2,-0.5,kNCuts+1.5,kNCuts+2,-0.5,kNCuts+1.5);
      for (Int_t i=0; i<kNCuts; i++) {
        fhCutsCorrelation->GetXaxis()->SetBinLabel(i+2,fgkCutsNames[i]);
        fhCutsCorrelation->GetYaxis()->SetBinLabel(i+2,fgkCutsNames[i]);
      }
      fHistogramsList->Add(fhCutsCorrelation);
    }

    if(fSystem  > kpp){
      fhCentrality[0] = new TH1F(Form("CentralityB_ %s",GetCutsString()),"Centrality before cut; centrality",400,0,100);
      fhCentrality[1] = new TH1F(Form("CentralityA_ %s",GetCutsString()),"Centrality; centrality",400,0,100);
      fHistogramsList->Add(fhCentrality[0]);
      fHistogramsList->Add(fhCentrality[1]);
    }
    else {
      /* for pp systems use multiplicity instead */
      fhCentrality[0] = new TH1F(Form("MultiplicityB_ %s",GetCutsString()),"Multiplicity before cut; multiplicity",400,0,400);
      fhCentrality[1] = new TH1F(Form("MultiplicityA_ %s",GetCutsString()),"Multiplicity; multiplicity",400,0,400);
      fHistogramsList->Add(fhCentrality[0]);
      fHistogramsList->Add(fhCentrality[1]);
    }

    fhVertexZ[0] = new TH1F(Form("VertexZB_%s",GetCutsString()),"Vertex Z; z_{vtx}",1000,-50,50);
    fhVertexZ[1] = new TH1F(Form("VertexZA_%s",GetCutsString()),"Vertex Z; z_{vtx}",1000,-50,50);
    fHistogramsList->Add(fhVertexZ[0]);
    fHistogramsList->Add(fhVertexZ[1]);

    fhTriggerClass[0] = new TH1F(Form("OfflineTriggerB_%s",GetCutsString()),"OfflineTrigger before cut",34,-0.5,33.5);
    fhTriggerClass[1] = new TH1F(Form("OfflineTriggerA_%s",GetCutsString()),"OfflineTrigger",34,-0.5,33.5);
    for (Int_t i = 0; i < 2; i++) {
      fhTriggerClass[i]->GetXaxis()->SetBinLabel( 1,"kMB");
      fhTriggerClass[i]->GetXaxis()->SetBinLabel( 2,"kINT7");
      fhTriggerClass[i]->GetXaxis()->SetBinLabel( 3,"kMUON");
      fhTriggerClass[i]->GetXaxis()->SetBinLabel( 4,"kHighMult");
      fhTriggerClass[i]->GetXaxis()->SetBinLabel( 5,"kEMC1");
      fhTriggerClass[i]->GetXaxis()->SetBinLabel( 6,"kCINT5");
      fhTriggerClass[i]->GetXaxis()->SetBinLabel( 7,"kCMUS5/kMUSPB");
      fhTriggerClass[i]->GetXaxis()->SetBinLabel( 8,"kMUSH7/kMUSHPB");
      fhTriggerClass[i]->GetXaxis()->SetBinLabel( 9,"kMUL7/kMuonLikePB");
      fhTriggerClass[i]->GetXaxis()->SetBinLabel(10,"kMUU7/kMuonUnlikePB");
      fhTriggerClass[i]->GetXaxis()->SetBinLabel(11,"kEMC7/kEMC8");
      fhTriggerClass[i]->GetXaxis()->SetBinLabel(12,"kMUS7");
      fhTriggerClass[i]->GetXaxis()->SetBinLabel(13,"kPHI1");
      fhTriggerClass[i]->GetXaxis()->SetBinLabel(14,"kPHI7/kPHI8/kPHOSPb");
      fhTriggerClass[i]->GetXaxis()->SetBinLabel(15,"kEMCEJE");
      fhTriggerClass[i]->GetXaxis()->SetBinLabel(16,"kEMCEGA");
      fhTriggerClass[i]->GetXaxis()->SetBinLabel(17,"kCentral");
      fhTriggerClass[i]->GetXaxis()->SetBinLabel(18,"kSemiCentral");
      fhTriggerClass[i]->GetXaxis()->SetBinLabel(19,"kDG5");
      fhTriggerClass[i]->GetXaxis()->SetBinLabel(20,"kZED");
      fhTriggerClass[i]->GetXaxis()->SetBinLabel(21,"kSPI7/kSPI");
      fhTriggerClass[i]->GetXaxis()->SetBinLabel(22,"kINT8");
      fhTriggerClass[i]->GetXaxis()->SetBinLabel(23,"kMuonSingleLowPt8");
      fhTriggerClass[i]->GetXaxis()->SetBinLabel(24,"kMuonSingleHighPt8");
      fhTriggerClass[i]->GetXaxis()->SetBinLabel(25,"kMuonLikeLowPt8");
      fhTriggerClass[i]->GetXaxis()->SetBinLabel(26,"kMuonUnlikeLowPt8");
      fhTriggerClass[i]->GetXaxis()->SetBinLabel(27,"kMuonUnlikeLowPt0");
      fhTriggerClass[i]->GetXaxis()->SetBinLabel(28,"kUserDefined");
      fhTriggerClass[i]->GetXaxis()->SetBinLabel(29,"kTRD");
      fhTriggerClass[i]->GetXaxis()->SetBinLabel(30,"kFastOnly");
      fhTriggerClass[i]->GetXaxis()->SetBinLabel(31,"kAnyINT");
      fhTriggerClass[i]->GetXaxis()->SetBinLabel(32,"kAny");
      fhTriggerClass[i]->GetXaxis()->SetBinLabel(33,"NOT kFastOnly");
      fhTriggerClass[i]->GetXaxis()->SetBinLabel(34,"Failed Physics Selection");
    }
    fHistogramsList->Add(fhTriggerClass[0]);
    fHistogramsList->Add(fhTriggerClass[1]);

    if(fQALevel == kQALevelHeavy){
      fhSPDClustersVsTracklets[0] = new TH2F(Form("SPDClustersVsTrackletsB_%s",GetCutsString()),"SPD clusters vs tracklets before cut;tracklets;SPD clusters",200,0,3000,500,0,10000);
      fhSPDClustersVsTracklets[1] = new TH2F(Form("SPDClustersVsTrackletsA_%s",GetCutsString()),"SPD clusters vs tracklets;tracklets;SPD clusters",200,0,3000,500,0,10000);
      fHistogramsList->Add(fhSPDClustersVsTracklets[0]);
      fHistogramsList->Add(fhSPDClustersVsTracklets[1]);

      fhV0MvsTracksTPCout[0] =
          new TH2F(Form("V0MvsTracksTPCoutB_%s", GetCutsString()),"V0 multiplicity vs tracks with kTPCout on before cut;# tracks with kTPCout on;V0 multiplicity",300,0,5000,300,0,40000);
      fhV0MvsTracksTPCout[1] =
          new TH2F(Form("V0MvsTracksTPCoutA_%s", GetCutsString()),"V0 multiplicity vs tracks with kTPCout on;# tracks with kTPCout on;V0 multiplicity",300,0,5000,300,0,40000);
      fHistogramsList->Add(fhV0MvsTracksTPCout[0]);
      fHistogramsList->Add(fhV0MvsTracksTPCout[1]);

      fhV0MvsTracksInitialTPCout[0] =
          new TH2F(Form("V0MvsTracksInitialTPCoutB_%s", GetCutsString()),"V0 multiplicity vs tracks with kTPCout on before cut;# tracks with kTPCout on (initial method);V0 multiplicity",300,0,30000,300,0,40000);
      fhV0MvsTracksInitialTPCout[1] =
          new TH2F(Form("V0MvsTracksInitialTPCoutA_%s", GetCutsString()),"V0 multiplicity vs tracks with kTPCout on;# tracks with kTPCout on (initial method);V0 multiplicity",300,0,30000,300,0,40000);
      fHistogramsList->Add(fhV0MvsTracksInitialTPCout[0]);
      fHistogramsList->Add(fhV0MvsTracksInitialTPCout[1]);

      const char *sel;
      const char *alt;
      this->GetCentralityEstimatorNames(sel,alt);
      if (sel != NULL && alt != NULL) {
        fhCentralityAltVsSel[0] =
            new TH2F(Form("CentralityAltVsSelB_%s", GetCutsString()),
                Form("Centrality, alternative vs selected before cut;centrality percentile (%s);centrality percentile (%s)", sel, alt),
                100,0,100,100,0,100);
        fhCentralityAltVsSel[1] =
            new TH2F(Form("CentralityAltVsSelA_%s", GetCutsString()),
                Form("Centrality, alternative vs selected;centrality percentile (%s);centrality percentile (%s)", sel, alt),
                100,0,100,100,0,100);
        fHistogramsList->Add(fhCentralityAltVsSel[0]);
        fHistogramsList->Add(fhCentralityAltVsSel[1]);
      }
      fhCL0vsV0MCentrality[0] =
          new TH2F(Form("CL0vsV0MCentralityB_%s",GetCutsString()),
              "Centrality, CL0 vs V0M, before cuts;centrality percentile (V0M);centrality percentile (CL0)",
              100,0,100,100,0,100);
      fhCL0vsV0MCentrality[1] =
          new TH2F(Form("CL0vsV0MCentralityA_%s",GetCutsString()),
              "Centrality, CL0 vs V0M;centrality percentile (V0M);centrality percentile (CL0)",
              100,0,100,100,0,100);
      fHistogramsList->Add(fhCL0vsV0MCentrality[0]);
      fHistogramsList->Add(fhCL0vsV0MCentrality[1]);

      fhESDvsTPConlyMultiplicity[0] =
          new TH2F(Form("ESDvsTPConlyMultiplicityB_%s",GetCutsString()),
              "Multiplicity, ESD vs TPC only, before cuts;multiplicity (TPC only tracks);multiplicity (ESD tracks)",
              400,0,7000,400,0,50000);
      fhESDvsTPConlyMultiplicity[1] =
          new TH2F(Form("ESDvsTPConlyMultiplicityA_%s",GetCutsString()),
              "Multiplicity, ESD vs TPC only;multiplicity (TPC only tracks);multiplicity (ESD tracks)",
              400,0,7000,400,0,50000);
      fHistogramsList->Add(fhESDvsTPConlyMultiplicity[0]);
      fHistogramsList->Add(fhESDvsTPConlyMultiplicity[1]);

      fhTOFvsGlobalMultiplicity[0] =
          new TH2F(Form("TOFvsGlobalMultiplicityB_%s",GetCutsString()),
              "global tracks in TOF vs global tracks, before cuts;multiplicity (global, FB32, tracks);multiplicity (global, FB32, TOF tracks)",
              400,0,4000,400,0,2000);
      fhTOFvsGlobalMultiplicity[1] =
          new TH2F(Form("TOFvsGlobalMultiplicityA_%s",GetCutsString()),
              "global tracks in TOF vs global tracks;multiplicity (global, FB32, tracks);multiplicity (global, FB32, TOF tracks)",
              400,0,4000,400,0,2000);
      fHistogramsList->Add(fhTOFvsGlobalMultiplicity[0]);
      fHistogramsList->Add(fhTOFvsGlobalMultiplicity[1]);

      fhAccTrkvsV0MCentrality[0] =
          new TH2F(Form("AccTrkvsV0MCentralityB_%s",GetCutsString()),
              "accepted global tracks vs V0M centrality, before cuts;centrality percentile (V0M);multiplicity (accepted global, FB32, tracks)",
              100,0,100,400,0,3000);
      fhAccTrkvsV0MCentrality[1] =
          new TH2F(Form("AccTrkvsV0MCentralityA_%s",GetCutsString()),
              "accepted global tracks vs V0M centrality;centrality percentile (V0M);multiplicity (accepted global, FB32, tracks)",
              100,0,100,400,0,3000);
      fHistogramsList->Add(fhAccTrkvsV0MCentrality[0]);
      fHistogramsList->Add(fhAccTrkvsV0MCentrality[1]);
    }

    TH1::AddDirectory(oldstatus);
  }
}



/// \cond CLASSIMP
ClassImp(AliCSEventCuts);
/// \endcond
