/**************************************************************************
 * Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id: AliRDHFJetsCuts.cxx 60376 2013-01-18 14:11:18Z fprino $ */

/////////////////////////////////////////////////////////////
//
// Base class for cuts on AOD reconstructed heavy-flavour Jets
//
// Authors: Andrea Rossi (andrea.rossi@cern.ch), Elena Bruna (elena.bruna@to.infn.it)
/////////////////////////////////////////////////////////////
#include <Riostream.h>

#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliVVertex.h"
#include "AliESDVertex.h"
#include "AliLog.h"
#include "AliAODVertex.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliESDtrackCuts.h"
#include "AliCentrality.h"
#include "AliAODRecoDecayHF.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliVertexerTracks.h"
#include "AliRDHFJetsCuts.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliPIDResponse.h"
#include "TRandom.h"
#include "AliEmcalJet.h"

#include "AliESDCaloCluster.h"
#include "AliEMCALRecoUtils.h"
#include "AliEMCALGeometry.h"

#include "AliPID.h"
#include "AliPIDResponse.h"

using std::cout;
using std::endl;

ClassImp(AliRDHFJetsCuts)

//--------------------------------------------------------------------------
AliRDHFJetsCuts::AliRDHFJetsCuts(const Char_t* name, const Char_t* title) :
  AliAnalysisCuts(name, title),
  fMinVtxType(3),
  fMinVtxContr(1),
  fMaxVtxRedChi2(1e6),
  fMaxVtxZ(10.),
  fMinSPDMultiplicity(0),
  fTriggerMask(AliVEvent::kAnyINT),
  fUseOnlyOneTrigger(kFALSE),
  fTrackCuts(0),
  fUseAOD049(kFALSE),
  fWhyRejection(0),
  fEvRejectionBits(0),
  fRemoveDaughtersFromPrimary(kFALSE),
  fRecomputePrimVertex(kFALSE),
  fUseMCVertex(kFALSE),
  fUsePhysicsSelection(kTRUE),
  fOptPileup(0),
  fMinContrPileup(3),
  fMinDzPileup(0.6),
  fUseCentrality(0),
  fMinCentrality(0.),
  fMaxCentrality(100.),
  fFixRefs(kFALSE),
  fMinPtJet(-1.),
  fMaxPtJet(100000.),
  fMaxEtaJet(1.),
  fJetRadius(1.),
  fKeepSignalMC(kFALSE),
  fApplySPDDeadPbPb2011(kFALSE),
  fApplySPDMisalignedPP2012(kFALSE),
  fRemoveTrackletOutliers(kFALSE),
  fCutOnzVertexSPD(0),
  fKinkReject(kFALSE),
  fUseTrackSelectionWithFilterBits(kTRUE),
  fUseCentrFlatteningInMC(kFALSE),
  fHistCentrDistr(0x0)
{
  //
  // Default Constructor
  //
  fTriggerClass[0] = "CINT1"; fTriggerClass[1] = "";
}
//--------------------------------------------------------------------------
AliRDHFJetsCuts::AliRDHFJetsCuts(const AliRDHFJetsCuts& source) :
  AliAnalysisCuts(source),
  fMinVtxType(source.fMinVtxType),
  fMinVtxContr(source.fMinVtxContr),
  fMaxVtxRedChi2(source.fMaxVtxRedChi2),
  fMaxVtxZ(source.fMaxVtxZ),
  fMinSPDMultiplicity(source.fMinSPDMultiplicity),
  fTriggerMask(source.fTriggerMask),
  fUseOnlyOneTrigger(source.fUseOnlyOneTrigger),
  fTriggerClass(),
  fTrackCuts(0),
  fUseAOD049(source.fUseAOD049),
  fWhyRejection(source.fWhyRejection),
  fEvRejectionBits(source.fEvRejectionBits),
  fRemoveDaughtersFromPrimary(source.fRemoveDaughtersFromPrimary),
  fRecomputePrimVertex(source.fRecomputePrimVertex),
  fUseMCVertex(source.fUseMCVertex),
  fUsePhysicsSelection(source.fUsePhysicsSelection),
  fOptPileup(source.fOptPileup),
  fMinContrPileup(source.fMinContrPileup),
  fMinDzPileup(source.fMinDzPileup),
  fUseCentrality(source.fUseCentrality),
  fMinCentrality(source.fMinCentrality),
  fMaxCentrality(source.fMaxCentrality),
  fFixRefs(source.fFixRefs),
  fMinPtJet(source.fMinPtJet),
  fMaxPtJet(source.fMaxPtJet),
  fMaxEtaJet(source.fMaxEtaJet),
  fJetRadius(source.fJetRadius),
  fKeepSignalMC(source.fKeepSignalMC),


  fApplySPDDeadPbPb2011(source.fApplySPDDeadPbPb2011),
  fApplySPDMisalignedPP2012(source.fApplySPDMisalignedPP2012),
  fRemoveTrackletOutliers(source.fRemoveTrackletOutliers),
  fCutOnzVertexSPD(source.fCutOnzVertexSPD),
  fKinkReject(source.fKinkReject),
  fUseTrackSelectionWithFilterBits(source.fUseTrackSelectionWithFilterBits),
  fUseCentrFlatteningInMC(source.fUseCentrFlatteningInMC),
  fHistCentrDistr(0x0)
{
  //
  // Copy constructor
  //
  fTriggerClass[0] = source.fTriggerClass[0];
  fTriggerClass[1] = source.fTriggerClass[1];
  if (source.GetTrackCuts()) AddTrackCuts(source.GetTrackCuts());




  if (source.fHistCentrDistr) fHistCentrDistr = (TH1F*)(source.fHistCentrDistr->Clone());
  PrintAll();

}
//--------------------------------------------------------------------------
AliRDHFJetsCuts& AliRDHFJetsCuts::operator=(const AliRDHFJetsCuts& source)
{
  //
  // assignment operator
  //
  if (&source == this) return *this;

  AliAnalysisCuts::operator=(source);

  fMinVtxType = source.fMinVtxType;
  fMinVtxContr = source.fMinVtxContr;
  fMaxVtxRedChi2 = source.fMaxVtxRedChi2;
  fMaxVtxZ = source.fMaxVtxZ;
  fMinSPDMultiplicity = source.fMinSPDMultiplicity;
  fTriggerMask = source.fTriggerMask;
  fUseOnlyOneTrigger = source.fUseOnlyOneTrigger;
  fTriggerClass[0] = source.fTriggerClass[0];
  fTriggerClass[1] = source.fTriggerClass[1];
  fUseAOD049 = source.fUseAOD049;
  fWhyRejection = source.fWhyRejection;
  fEvRejectionBits = source.fEvRejectionBits;
  fRemoveDaughtersFromPrimary = source.fRemoveDaughtersFromPrimary;
  fRecomputePrimVertex = source.fRecomputePrimVertex;
  fUseMCVertex = source.fUseMCVertex;
  fUsePhysicsSelection = source.fUsePhysicsSelection;
  fOptPileup = source.fOptPileup;
  fMinContrPileup = source.fMinContrPileup;
  fMinDzPileup = source.fMinDzPileup;
  fUseCentrality = source.fUseCentrality;
  fMinCentrality = source.fMinCentrality;
  fMaxCentrality = source.fMaxCentrality;
  fFixRefs = source.fFixRefs;
  fMinPtJet = source.fMinPtJet;
  fMaxPtJet = source.fMaxPtJet;
  fMaxEtaJet = source.fMaxEtaJet;
  fJetRadius = source.fJetRadius;
  fKeepSignalMC = source.fKeepSignalMC;
  fApplySPDDeadPbPb2011 = source.fApplySPDDeadPbPb2011;
  fApplySPDMisalignedPP2012 = source.fApplySPDMisalignedPP2012;
  fRemoveTrackletOutliers = source.fRemoveTrackletOutliers;
  fCutOnzVertexSPD = source.fCutOnzVertexSPD;
  fKinkReject = source.fKinkReject;
  fUseTrackSelectionWithFilterBits = source.fUseTrackSelectionWithFilterBits;
  if (fHistCentrDistr) delete fHistCentrDistr;
  fUseCentrFlatteningInMC = source.fUseCentrFlatteningInMC;
  if (source.fHistCentrDistr)fHistCentrDistr = (TH1F*)(source.fHistCentrDistr->Clone());

  if (source.GetTrackCuts()) {delete fTrackCuts; fTrackCuts = new AliESDtrackCuts(*(source.GetTrackCuts()));}


  PrintAll();

  return *this;
}
//--------------------------------------------------------------------------
AliRDHFJetsCuts::~AliRDHFJetsCuts()
{
  //
  // Default Destructor
  //

  if (fHistCentrDistr)delete fHistCentrDistr;
  if (fTrackCuts) delete fTrackCuts;
}
//---------------------------------------------------------------------------
Int_t AliRDHFJetsCuts::IsEventSelectedInCentrality(AliVEvent* event)
{
  //
  // Centrality selection
  //
  if (fUseCentrality < kCentOff || fUseCentrality >= kCentInvalid) {
    AliWarning("Centrality estimator not valid");
    return 3;
  } else {
    Float_t centvalue = GetCentrality((AliAODEvent*)event);
    if (centvalue < -998.) { //-999 if no centralityP
      return 0;
    } else {
      if (centvalue < fMinCentrality || centvalue > fMaxCentrality) {
        return 2;
      }
      // selection to flatten the centrality distribution
      if (fHistCentrDistr) {
        if (!IsEventSelectedForCentrFlattening(centvalue))return 4;
      }
    }
  }
  return 0;
}


//-------------------------------------------------
void AliRDHFJetsCuts::SetHistoForCentralityFlattening(TH1F* h, Double_t minCentr, Double_t maxCentr, Double_t centrRef, Int_t switchTRand)
{
  // set the histo for centrality flattening
  // the centrality is flatten in the range minCentr,maxCentr
  // if centrRef is zero, the minimum in h within (minCentr,maxCentr) defines the reference
  //                positive, the value of h(centrRef) defines the reference (-> the centrality distribution might be not flat in the whole desired range)
  //                negative, h(bin with max in range)*centrRef is used to define the reference (-> defines the maximum loss of events, also in this case the distribution might be not flat)
  // switchTRand is used to set the unerflow bin of the histo: if it is < -1 in the analysis the random event selection will be done on using TRandom

  if (maxCentr < minCentr) {
    AliWarning("AliRDHFJetsCuts::Wrong centralities values while setting the histogram for centrality flattening");
  }

  if (fHistCentrDistr)delete fHistCentrDistr;
  fHistCentrDistr = (TH1F*)h->Clone("hCentralityFlat");
  fHistCentrDistr->SetTitle("Reference histo for centrality flattening");
  Int_t minbin = fHistCentrDistr->FindBin(minCentr * 1.00001); // fast if fix bin width
  Int_t maxbin = fHistCentrDistr->FindBin(maxCentr * 0.9999);
  fHistCentrDistr->GetXaxis()->SetRange(minbin, maxbin);
  Double_t ref = 0., bincont = 0., binrefwidth = 1.;
  Int_t binref = 0;
  if (TMath::Abs(centrRef) < 0.0001) {
    binref = fHistCentrDistr->GetMinimumBin();
    binrefwidth = fHistCentrDistr->GetBinWidth(binref);
    ref = fHistCentrDistr->GetBinContent(binref) / binrefwidth;
  } else if (centrRef > 0.) {
    binref = h->FindBin(centrRef);
    if (binref < 1 || binref > h->GetNbinsX()) {
      AliWarning("AliRDHFJetsCuts::Wrong centrality reference value while setting the histogram for centrality flattening");
    }
    binrefwidth = fHistCentrDistr->GetBinWidth(binref);
    ref = fHistCentrDistr->GetBinContent(binref) / binrefwidth;
  } else {
    if (centrRef < -1) AliWarning("AliRDHFJetsCuts: with this centrality reference no flattening will be applied");
    binref = fHistCentrDistr->GetMaximumBin();
    binrefwidth = fHistCentrDistr->GetBinWidth(binref);
    ref = fHistCentrDistr->GetMaximum() * TMath::Abs(centrRef) / binrefwidth;
  }

  for (Int_t j = 1; j <= h->GetNbinsX(); j++) { // Now set the "probabilities"
    if (h->GetBinLowEdge(j) * 1.0001 >= minCentr && h->GetBinLowEdge(j + 1) * 0.9999 <= maxCentr) {
      bincont = h->GetBinContent(j);
      fHistCentrDistr->SetBinContent(j, ref / bincont * h->GetBinWidth(j));
      fHistCentrDistr->SetBinError(j, h->GetBinError(j)*ref / bincont);
    } else {
      h->SetBinContent(j, 1.1); // prob > 1 to assure that events will not be rejected
    }
  }

  fHistCentrDistr->SetBinContent(0, switchTRand);
  return;

}

//-------------------------------------------------
Bool_t AliRDHFJetsCuts::IsEventSelectedForCentrFlattening(Float_t centvalue)
{
  //
  //  Random event selection, based on fHistCentrDistr, to flatten the centrality distribution
  //  Can be faster if it was required that fHistCentrDistr covers
  //  exactly the desired centrality range (e.g. part of the lines below should be done during the
  // setting of the histo) and TH1::SetMinimum called
  //

  if (!fHistCentrDistr) return kTRUE;
  // Int_t maxbin=fHistCentrDistr->FindBin(fMaxCentrality*0.9999);
  //   if(maxbin>fHistCentrDistr->GetNbinsX()){
  //     AliWarning("AliRDHFJetsCuts: The maximum centrality exceeds the x-axis limit of the histogram for centrality flattening");
  //   }

  Int_t bin = fHistCentrDistr->FindBin(centvalue); // Fast if the histo has a fix bin
  Double_t bincont = fHistCentrDistr->GetBinContent(bin);
  Double_t centDigits = centvalue - (Int_t)(centvalue * 100.) / 100.; // this is to extract a random number between 0 and 0.01

  if (fHistCentrDistr->GetBinContent(0) < -0.9999) {
    if (gRandom->Uniform(1.) < bincont)return kTRUE;
    return kFALSE;
  }

  if (centDigits * 100. < bincont)return kTRUE;
  return kFALSE;

}
// //---------------------------------------------------------------------------
// void AliRDHFJetsCuts::SetupPID(AliVEvent *event) {
//   // Set the PID response object in the AliAODPidHF
//   // in case of old PID sets the TPC dE/dx BB parameterization

//   if(fPidHF){
//     if(fPidHF->GetPidResponse()==0x0){
//       AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
//       AliInputEventHandler *inputHandler=(AliInputEventHandler*)mgr->GetInputEventHandler();
//       AliPIDResponse *pidResp=inputHandler->GetPIDResponse();
//       fPidHF->SetPidResponse(pidResp);
//     }
//     if(fPidHF->GetUseCombined()) fPidHF->SetUpCombinedPID();
//     if(fPidHF->GetOldPid()) {

//       Bool_t isMC=kFALSE;
//       TClonesArray *mcArray = (TClonesArray*)((AliAODEvent*)event)->GetList()->FindObject(AliAODMCParticle::StdBranchName());
//       if(mcArray) {isMC=kTRUE;fUseAOD049=kFALSE;}

//       // pp, from LHC10d onwards
//       if((event->GetRunNumber()>121693 && event->GetRunNumber()<136851) ||
//   event->GetRunNumber()>139517) fPidHF->SetOnePad(kTRUE);
//       // pp, 2011 low energy run
//       if((event->GetRunNumber()>=146686 && event->GetRunNumber()<=146860)){
//  fPidHF->SetppLowEn2011(kTRUE);
//  fPidHF->SetOnePad(kFALSE);
//       }
//       // PbPb LHC10h
//       if(event->GetRunNumber()>=136851 && event->GetRunNumber()<=139517) fPidHF->SetPbPb(kTRUE);
//       // MC
//       if(isMC) fPidHF->SetMC(kTRUE);
//       if(isMC && (event->GetRunNumber()>=146686 && event->GetRunNumber()<=146860))
//  fPidHF->SetMClowenpp2011(kTRUE);
//       fPidHF->SetBetheBloch();
//     }else{
//       // check that AliPIDResponse object was properly set in case of using OADB
//       if(fPidHF->GetPidResponse()==0x0) AliFatal("AliPIDResponse object not set");
//     }
//   }
// }
//---------------------------------------------------------------------------
Bool_t AliRDHFJetsCuts::IsEventSelected(AliVEvent* event)
{
  //
  // Event selection
  //

  fWhyRejection = 0;
  fEvRejectionBits = 0;
  Bool_t accept = kTRUE;

  if (fRecomputePrimVertex) {
    Bool_t vertOK = RecomputePrimaryVertex((AliAODEvent*)event);
    if (!vertOK) {
      fWhyRejection = 6;
      return kFALSE;
    }
  }

  // check if it's MC
  Bool_t isMC = kFALSE;
  TClonesArray* mcArray = (TClonesArray*)((AliAODEvent*)event)->GetList()->FindObject(AliAODMCParticle::StdBranchName());
  if (mcArray) {isMC = kTRUE; fUseAOD049 = kFALSE;}


  // trigger class
  TString firedTriggerClasses = ((AliAODEvent*)event)->GetFiredTriggerClasses();

  // don't do for MC and for PbPb 2010 data
  // The pPb MC productions do have fired trigger info so this must be altered
  if (!isMC && (event->GetRunNumber() < 136851 || event->GetRunNumber() > 139517)) {
    // if(event->GetRunNumber()<136851 || event->GetRunNumber()>139517) {
    if (!firedTriggerClasses.Contains(fTriggerClass[0].Data()) &&
        (fTriggerClass[1].CompareTo("") == 0 || !firedTriggerClasses.Contains(fTriggerClass[1].Data())) ) {
      fWhyRejection = 5;
      fEvRejectionBits += 1 << kNotSelTrigger;
      accept = kFALSE;
    }
  }

  // TEMPORARY FIX FOR GetEvent
  Int_t nTracks = ((AliAODEvent*)event)->GetNumberOfTracks();
  for (Int_t itr = 0; itr < nTracks; itr++) {
    AliAODTrack* tr = (AliAODTrack*)((AliAODEvent*)event)->GetTrack(itr);
    tr->SetAODEvent((AliAODEvent*)event);
  }

  // physics selection requirements
  if (fUsePhysicsSelection) {
    Bool_t isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & fTriggerMask);
    if (!isSelected) {
      if (accept) fWhyRejection = 7;
      fEvRejectionBits += 1 << kPhysicsSelection;
      accept = kFALSE;
    } else {
      if (fUseOnlyOneTrigger) {
        if (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() != fTriggerMask) {
          if (accept) fWhyRejection = 7;
          fEvRejectionBits += 1 << kPhysicsSelection;
          accept = kFALSE;
        }
      }
    }
  }

  // vertex requirements

  const AliVVertex* vertex = event->GetPrimaryVertex();

  if (!vertex) {
    accept = kFALSE;
    fEvRejectionBits += 1 << kNoVertex;
  } else {
    TString title = vertex->GetTitle();
    if (title.Contains("Z") && fMinVtxType > 1) {
      accept = kFALSE;
      fEvRejectionBits += 1 << kNoVertex;
    } else if (title.Contains("3D") && fMinVtxType > 2) {
      accept = kFALSE;
      fEvRejectionBits += 1 << kNoVertex;
    }
    if (vertex->GetNContributors() < fMinVtxContr) {
      accept = kFALSE;
      fEvRejectionBits += 1 << kTooFewVtxContrib;
    }
    if (TMath::Abs(vertex->GetZ()) > fMaxVtxZ) {
      fEvRejectionBits += 1 << kZVtxOutFid;
      if (accept) fWhyRejection = 6;
      accept = kFALSE;
    }
  }

  if (fCutOnzVertexSPD > 0) {
    const AliVVertex* vSPD = ((AliAODEvent*)event)->GetPrimaryVertexSPD();
    if (!vSPD || (vSPD && vSPD->GetNContributors() < fMinVtxContr)) {
      accept = kFALSE;
      fEvRejectionBits += 1 << kBadSPDVertex;
    } else {
      if (fCutOnzVertexSPD == 1 && TMath::Abs(vSPD->GetZ()) > 12.) {
        fEvRejectionBits += 1 << kZVtxSPDOutFid;
        if (accept) fWhyRejection = 6;
        accept = kFALSE;
      }
      if (fCutOnzVertexSPD == 2 && vertex) {
        if (TMath::Abs(vSPD->GetZ() - vertex->GetZ()) > 0.5) {
          fEvRejectionBits += 1 << kZVtxSPDOutFid;
          if (accept) fWhyRejection = 6;
          accept = kFALSE;
        }
      }
    }
  }

  // pile-up rejection
  if (fOptPileup == kRejectPileupEvent) {
    Int_t cutc = (Int_t)fMinContrPileup;
    Double_t cutz = (Double_t)fMinDzPileup;
    if (event->IsPileupFromSPD(cutc, cutz, 3., 2., 10.)) {
      if (accept) fWhyRejection = 1;
      fEvRejectionBits += 1 << kPileupSPD;
      accept = kFALSE;
    }
  }

  // centrality selection
  if (fUseCentrality != kCentOff) {
    Int_t rejection = IsEventSelectedInCentrality(event);
    Bool_t okCent = kFALSE;
    if (rejection == 0) okCent = kTRUE;
    if (isMC && rejection == 4 && !fUseCentrFlatteningInMC) okCent = kTRUE;
    if (!okCent) {
      if (!isMC) {
        if (accept) fWhyRejection = rejection;
        if (fWhyRejection == 4)fEvRejectionBits += 1 << kCentralityFlattening;
        else fEvRejectionBits += 1 << kOutsideCentrality;
        accept = kFALSE;
      } else {
        if (rejection != 4) {
          if (accept) fWhyRejection = rejection;
          fEvRejectionBits += 1 << kOutsideCentrality;
          accept = kFALSE;
        }
      }
    }

  }

  // PbPb2011 outliers in tracklets vs. VZERO
  if (event->GetRunNumber() >= 167693 && event->GetRunNumber() <= 170593) {
    if (fRemoveTrackletOutliers) {
      Double_t v0cent = GetCentrality((AliAODEvent*)event, kCentV0M);
      Double_t ntracklets = ((AliAODEvent*)event)->GetTracklets()->GetNumberOfTracklets();
      Double_t cutval = 60. - 0.08 * ntracklets + 1. / 50000.*ntracklets * ntracklets;
      if (ntracklets < 1000. && v0cent < cutval) {
        if (accept) fWhyRejection = 2;
        fEvRejectionBits += 1 << kOutsideCentrality;
        accept = kFALSE;
      }
    }

    if (fApplySPDMisalignedPP2012 && !(event->GetRunNumber() >= 195681 && event->GetRunNumber() <= 197388)) fApplySPDMisalignedPP2012 = false;
  }

  return accept;
}
//---------------------------------------------------------------------------
Bool_t AliRDHFJetsCuts::IsJetSelected(const AliEmcalJet* jet) const
{

  Bool_t accepted = kFALSE;
  if (jet->Pt() > fMinPtJet && jet->Pt() < fMaxPtJet && TMath::Abs(jet->Eta()) < fMaxEtaJet) accepted = kTRUE;

  return accepted;
}

//---------------------------------------------------------------------------
Bool_t AliRDHFJetsCuts::AreDaughtersSelected(AliAODRecoDecayHF* d) const
{
  //
  // Daughter tracks selection
  //
  if (!fTrackCuts) return kTRUE;

  Int_t ndaughters = d->GetNDaughters();
  AliAODVertex* vAOD = d->GetPrimaryVtx();
  Double_t pos[3], cov[6];
  vAOD->GetXYZ(pos);
  vAOD->GetCovarianceMatrix(cov);
  const AliESDVertex vESD(pos, cov, 100., 100);

  Bool_t retval = kTRUE;

  for (Int_t idg = 0; idg < ndaughters; idg++) {
    AliAODTrack* dgTrack = (AliAODTrack*)d->GetDaughter(idg);
    if (!dgTrack) {retval = kFALSE; continue;}
    //printf("charge %d\n",dgTrack->Charge());
    if (dgTrack->Charge() == 0) continue; // it's not a track, but a V0


    if (!IsDaughterSelected(dgTrack, &vESD, fTrackCuts)) retval = kFALSE;
  }

  return retval;
}
//---------------------------------------------------------------------------
void AliRDHFJetsCuts::IsDaugElectron(AliAODEvent* aod, const AliAODJet* jet, Bool_t& fFlagElec)
{
  //
  // Convert to ESDtrack, get emcal info and id
  //
  Bool_t flagElec = kFALSE;

  TRefArray* reftracks = (TRefArray*)jet->GetRefTracks();
  Int_t ntrks = reftracks->GetEntriesFast();

  for (Int_t ktracks = 0; ktracks < ntrks; ktracks++) {

    AliAODTrack* aodtrack = (AliAODTrack*)(jet->GetRefTracks()->At(ktracks));
    AliESDtrack* esdTrack = new AliESDtrack(aodtrack);

    //Double_t fClsE = -999, p = -999, fEovp=-999, pt = -999, fdEdx=-999, fTPCnSigma=0, charge=-99;
    //Double_t fClsE = -999, p = -999, fEovp=-999, pt = -999, fdEdx=-999, charge=-99;
    //Double_t fClsE = -999, p = -999, fEovp=-999, fdEdx=-999, charge=-99;
    Double_t fClsE = -999, p = -999, fEovp = -999, fdEdx = -999;
    // Track extrapolation to EMCAL
    Int_t fClsId = esdTrack->GetEMCALcluster();
    p = esdTrack->P();
    //pt = esdTrack->Pt();
    //charge = esdTrack->Charge();

    if (fClsId > 0) {
      AliVCluster* cluster = aod->GetCaloCluster(fClsId);
      if (TMath::Abs(cluster->GetTrackDx()) < 0.05 && TMath::Abs(cluster->GetTrackDz()) < 0.05) {

        fClsE = cluster->E();
        fEovp = fClsE / p;
        fdEdx = esdTrack->GetTPCsignal();

        if ((fEovp > 0.8 && fEovp < 1.1) && (fdEdx > 75 && fdEdx < 105)  && !flagElec) flagElec = kTRUE;
      }
    }
  }
  fFlagElec = flagElec;
}
//---------------------------------------------------------------------------
Bool_t AliRDHFJetsCuts::IsDaughterSelected(AliAODTrack* track, const AliESDVertex* primary, AliESDtrackCuts* cuts) const
{
  //
  // Convert to ESDtrack, relate to vertex and check cuts
  //
  if (!cuts) return kTRUE;

  if (cuts->GetFlagCutTOFdistance()) cuts->SetFlagCutTOFdistance(kFALSE);

  Bool_t retval = kTRUE;

  // convert to ESD track here
  AliESDtrack esdTrack(track);
  // set the TPC cluster info
  esdTrack.SetTPCClusterMap(track->GetTPCClusterMap());
  esdTrack.SetTPCSharedMap(track->GetTPCSharedMap());
  esdTrack.SetTPCPointsF(track->GetTPCNclsF());
  // needed to calculate the impact parameters
  esdTrack.RelateToVertex(primary, 0., 3.);

  if (!cuts->IsSelected(&esdTrack)) {
    AliDebug(3, "IsDaughterSelected Track rejected by IsSelected");
    retval = kFALSE;
  }


  if (fKinkReject) {
    AliAODVertex* maybeKink = track->GetProdVertex();
    if (maybeKink->GetType() == AliAODVertex::kKink) retval = kFALSE;
  }

  if (fOptPileup == kRejectTracksFromPileupVertex) {
    // to be implemented
    // we need either to have here the AOD Event,
    // or to have the pileup vertex object
  }

  if (fApplySPDDeadPbPb2011) {
    Bool_t deadSPDLay1PbPb2011[20][4] = {
      {kTRUE, kTRUE, kTRUE, kTRUE},
      {kTRUE, kTRUE, kTRUE, kTRUE},
      {kTRUE, kTRUE, kTRUE, kTRUE},
      {kTRUE, kTRUE, kTRUE, kTRUE},
      {kTRUE, kTRUE, kTRUE, kTRUE},
      {kFALSE, kFALSE, kTRUE, kTRUE},
      {kTRUE, kTRUE, kFALSE, kFALSE},
      {kTRUE, kTRUE, kTRUE, kTRUE},
      {kFALSE, kFALSE, kTRUE, kTRUE},
      {kTRUE, kTRUE, kTRUE, kTRUE},
      {kTRUE, kTRUE, kFALSE, kFALSE},
      {kTRUE, kTRUE, kTRUE, kTRUE},
      {kFALSE, kFALSE, kFALSE, kFALSE},
      {kFALSE, kFALSE, kTRUE, kTRUE},
      {kFALSE, kFALSE, kFALSE, kFALSE},
      {kFALSE, kFALSE, kFALSE, kFALSE},
      {kTRUE, kTRUE, kTRUE, kTRUE},
      {kTRUE, kTRUE, kFALSE, kFALSE},
      {kFALSE, kFALSE, kFALSE, kFALSE},
      {kFALSE, kFALSE, kFALSE, kFALSE}
    };
    Bool_t deadSPDLay2PbPb2011[40][4] = {
      {kTRUE, kTRUE, kTRUE, kTRUE},
      {kTRUE, kTRUE, kTRUE, kTRUE},
      {kTRUE, kTRUE, kTRUE, kTRUE},
      {kTRUE, kTRUE, kTRUE, kTRUE},
      {kTRUE, kTRUE, kTRUE, kTRUE},
      {kTRUE, kTRUE, kTRUE, kTRUE},
      {kTRUE, kTRUE, kTRUE, kTRUE},
      {kTRUE, kTRUE, kTRUE, kTRUE},
      {kTRUE, kTRUE, kTRUE, kTRUE},
      {kTRUE, kTRUE, kTRUE, kTRUE},
      {kTRUE, kTRUE, kTRUE, kTRUE},
      {kTRUE, kTRUE, kTRUE, kTRUE},
      {kFALSE, kFALSE, kFALSE, kFALSE},
      {kFALSE, kFALSE, kTRUE, kTRUE},
      {kTRUE, kTRUE, kTRUE, kTRUE},
      {kTRUE, kTRUE, kTRUE, kTRUE},
      {kTRUE, kTRUE, kFALSE, kFALSE},
      {kTRUE, kTRUE, kTRUE, kTRUE},
      {kTRUE, kTRUE, kTRUE, kTRUE},
      {kTRUE, kTRUE, kTRUE, kTRUE},
      {kFALSE, kFALSE, kFALSE, kFALSE},
      {kFALSE, kFALSE, kFALSE, kFALSE},
      {kTRUE, kTRUE, kTRUE, kTRUE},
      {kTRUE, kTRUE, kTRUE, kTRUE},
      {kFALSE, kFALSE, kFALSE, kFALSE},
      {kFALSE, kFALSE, kFALSE, kFALSE},
      {kTRUE, kTRUE, kTRUE, kTRUE},
      {kTRUE, kTRUE, kTRUE, kTRUE},
      {kFALSE, kFALSE, kFALSE, kFALSE},
      {kFALSE, kFALSE, kFALSE, kFALSE},
      {kFALSE, kFALSE, kFALSE, kFALSE},
      {kFALSE, kFALSE, kFALSE, kFALSE},
      {kTRUE, kTRUE, kTRUE, kTRUE},
      {kTRUE, kTRUE, kTRUE, kTRUE},
      {kTRUE, kTRUE, kFALSE, kFALSE},
      {kTRUE, kTRUE, kTRUE, kTRUE},
      {kFALSE, kFALSE, kFALSE, kFALSE},
      {kFALSE, kFALSE, kFALSE, kFALSE},
      {kFALSE, kFALSE, kFALSE, kFALSE},
      {kFALSE, kFALSE, kFALSE, kFALSE}
    };
    Double_t xyz1[3], xyz2[3];
    esdTrack.GetXYZAt(3.9, 0., xyz1);
    esdTrack.GetXYZAt(7.6, 0., xyz2);
    Double_t phi1 = TMath::ATan2(xyz1[1], xyz1[0]);
    if (phi1 < 0) phi1 += 2 * TMath::Pi();
    Int_t lad1 = (Int_t)(phi1 / (2.*TMath::Pi() / 20.));
    Double_t phi2 = TMath::ATan2(xyz2[1], xyz2[0]);
    if (phi2 < 0) phi2 += 2 * TMath::Pi();
    Int_t lad2 = (Int_t)(phi2 / (2.*TMath::Pi() / 40.));
    Int_t mod1 = (Int_t)((xyz1[2] + 14) / 7.);
    Int_t mod2 = (Int_t)((xyz2[2] + 14) / 7.);
    Bool_t lay1ok = kFALSE;
    if (mod1 >= 0 && mod1 < 4 && lad1 < 20) {
      lay1ok = deadSPDLay1PbPb2011[lad1][mod1];
    }
    Bool_t lay2ok = kFALSE;
    if (mod2 >= 0 && mod2 < 4 && lad2 < 40) {
      lay2ok = deadSPDLay2PbPb2011[lad2][mod2];
    }
    if (!lay1ok && !lay2ok) retval = kFALSE;
  }

  if (fApplySPDMisalignedPP2012) {
    // Cut tracks crossing the SPD at 5.6<phi<2pi
    Double_t xyz1[3], xyz2[3];
    esdTrack.GetXYZAt(3.9, 0., xyz1);
    esdTrack.GetXYZAt(7.6, 0., xyz2);
    Double_t phi1 = TMath::ATan2(xyz1[1], xyz1[0]);
    if (phi1 < 0) phi1 += 2 * TMath::Pi();
    Double_t phi2 = TMath::ATan2(xyz2[1], xyz2[0]);
    if (phi2 < 0) phi2 += 2 * TMath::Pi();
    Bool_t lay1ok = kTRUE;
    if (phi1 > 5.6 && phi1 < 2.*TMath::Pi()) lay1ok = kFALSE;
    Bool_t lay2ok = kTRUE;
    if (phi2 > 5.6 && phi2 < 2.*TMath::Pi()) lay2ok = kFALSE;
    if (!lay1ok || !lay2ok) return kFALSE;
  }

  return retval;
}

//---------------------------------------------------------------------------
void AliRDHFJetsCuts::SetUseCentrality(Int_t flag)
{
  //
  // set centrality estimator
  //
  fUseCentrality = flag;
  if (fUseCentrality < kCentOff || fUseCentrality >= kCentInvalid) AliWarning("Centrality estimator not valid");

  return;
}



//---------------------------------------------------------------------------
void AliRDHFJetsCuts::PrintAll() const
{
  //
  // print all cuts values
  //

  printf("Minimum vtx type %d\n", fMinVtxType);
  printf("Minimum vtx contr %d\n", fMinVtxContr);
  printf("Max vtx red chi2 %f\n", fMaxVtxRedChi2);
  printf("Min SPD mult %d\n", fMinSPDMultiplicity);
  //  printf("Use PID %d  OldPid=%d\n",(Int_t)fUsePID,fPidHF ? fPidHF->GetOldPid() : -1);
  printf("Remove daughters from vtx %d\n", (Int_t)fRemoveDaughtersFromPrimary);
  printf("Recompute primary vertex %d\n", (Int_t)fRecomputePrimVertex);
  printf("Physics selection: %s\n", fUsePhysicsSelection ? "Yes" : "No");
  printf("Pileup rejection: %s\n", (fOptPileup > 0) ? "Yes" : "No");
  if (fOptPileup == 1) printf(" -- Reject pileup event");
  if (fOptPileup == 2) printf(" -- Reject tracks from pileup vtx");
  if (fUseCentrality > 0) {
    TString estimator = "";
    if (fUseCentrality == 1) estimator = "V0";
    if (fUseCentrality == 2) estimator = "Tracks";
    if (fUseCentrality == 3) estimator = "Tracklets";
    if (fUseCentrality == 4) estimator = "SPD clusters outer";
    printf("Centrality class considered: %.1f-%.1f, estimated with %s", fMinCentrality, fMaxCentrality, estimator.Data());
  }

  return;
}

//--------------------------------------------------------------------------
void AliRDHFJetsCuts::PrintTrigger() const
{
  // print the trigger selection

  printf("Selected trigger classes: %s %s\n", fTriggerClass[0].Data(), fTriggerClass[1].Data());

  cout << " Trigger selection pattern: ";
  if ( fTriggerMask & AliVEvent::kAny ) cout << " kAny ";
  if ( fTriggerMask & AliVEvent::kAnyINT ) cout << " kAnyINT ";
  if ( fTriggerMask & AliVEvent::kINT7 ) cout << " kINT7 ";
  if ( fTriggerMask & AliVEvent::kMB ) cout << " kMB ";
  if ( fTriggerMask & AliVEvent::kCINT5 ) cout << " kCINT5 ";
  if ( fTriggerMask & AliVEvent::kCentral ) cout << " kCentral ";
  if ( fTriggerMask & AliVEvent::kSemiCentral ) cout << " kSemiCentral ";
  if ( fTriggerMask & AliVEvent::kEMCEGA ) cout << " kEMCEGA ";
  if ( fTriggerMask & AliVEvent::kHighMult ) cout << " kHighMult ";
  if ( fTriggerMask & AliVEvent::kFastOnly ) cout << " kFastOnly ";
  cout << endl << endl;

}

//-------------------------------------------------------------------
Float_t AliRDHFJetsCuts::GetCentrality(AliAODEvent* aodEvent, AliRDHFJetsCuts::ECentrality estimator)
{
  //
  // Get centrality percentile
  //

  TClonesArray* mcArray = (TClonesArray*)((AliAODEvent*)aodEvent)->GetList()->FindObject(AliAODMCParticle::StdBranchName());
  if (mcArray) {fUseAOD049 = kFALSE;}

  AliAODHeader* header = static_cast <AliAODHeader*> (aodEvent->GetHeader());
  AliCentrality* centrality = header->GetCentralityP();
  Float_t cent = -999.;
  Bool_t isSelRun = kFALSE;
  Int_t selRun[5] = {138364, 138826, 138828, 138836, 138871};
  if (!centrality) return cent;
  else {
    if (estimator == kCentV0M) {
      cent = (Float_t)(centrality->GetCentralityPercentile("V0M"));
      if (cent < 0) {
        Int_t quality = centrality->GetQuality();
        if (quality <= 1) { // fQuality==1 means rejected by zVertex cut that we apply a part and we want to keep separate (Giacomo)
          cent = (Float_t)centrality->GetCentralityPercentileUnchecked("V0M");
        } else {
          Int_t runnum = aodEvent->GetRunNumber();
          for (Int_t ir = 0; ir < 5; ir++) {
            if (runnum == selRun[ir]) {
              isSelRun = kTRUE;
              break;
            }
          }
          if ((quality == 8 || quality == 9) && isSelRun)cent = (Float_t)centrality->GetCentralityPercentileUnchecked("V0M");
        }
      }

      //temporary fix for AOD049 outliers
      if (fUseAOD049 && cent >= 0) {
        Float_t v0 = 0;
        AliAODVZERO* aodV0 = aodEvent->GetVZEROData();
        v0 += aodV0->GetMTotV0A();
        v0 += aodV0->GetMTotV0C();
        if (cent == 0 && v0 < 19500)return -1; //filtering issue
        Float_t tkl = (Float_t)(aodEvent->GetTracklets()->GetNumberOfTracklets());
        Float_t val = 1.30552 +  0.147931 * v0;
        Float_t tklSigma[101] = {176.644, 156.401, 153.789, 153.015, 142.476, 137.951, 136.127, 129.852, 127.436, 124.86, 120.788, 115.611, 113.172, 110.496, 109.127, 104.421, 102.479, 99.9766, 97.5152, 94.0654, 92.4602, 89.3364, 87.1342, 83.3497, 82.6216, 81.1084, 78.0793, 76.1234, 72.9434, 72.1334, 68.0056, 68.2755, 66.0376, 62.9666, 62.4274, 59.65, 58.3776, 56.6361, 54.5184, 53.4224, 51.932, 50.8922, 48.2848, 47.912, 46.5717, 43.4114, 43.2083, 41.3065, 40.1863, 38.5255, 37.2851, 37.5396, 34.4949, 33.8366, 31.8043, 31.7412, 30.8392, 30.0274, 28.8793, 27.6398, 26.6488, 25.0183, 25.1489, 24.4185, 22.9107, 21.2002, 21.6977, 20.1242, 20.4963, 19.0235, 19.298, 17.4103, 16.868, 15.2939, 15.2939, 16.0295, 14.186, 14.186, 15.2173, 12.9504, 12.9504, 12.9504, 15.264, 12.3674, 12.3674, 12.3674, 12.3674, 12.3674, 18.3811, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544};
        if ( TMath::Abs(tkl - val) > 6.*tklSigma[(Int_t)cent] )return -1; //outlier
      }
    } else {
      if (estimator == kCentTRK) {
        cent = (Float_t)(centrality->GetCentralityPercentile("TRK"));
        if (cent < 0) {
          Int_t quality = centrality->GetQuality();
          if (quality <= 1) {
            cent = (Float_t)centrality->GetCentralityPercentileUnchecked("TRK");
          } else {
            Int_t runnum = aodEvent->GetRunNumber();
            for (Int_t ir = 0; ir < 5; ir++) {
              if (runnum == selRun[ir]) {
                isSelRun = kTRUE;
                break;
              }
            }
            if ((quality == 8 || quality == 9) && isSelRun)cent = (Float_t)centrality->GetCentralityPercentileUnchecked("TRK");
          }
        }
      } else {
        if (estimator == kCentTKL) {
          cent = (Float_t)(centrality->GetCentralityPercentile("TKL"));
          if (cent < 0) {
            Int_t quality = centrality->GetQuality();
            if (quality <= 1) {
              cent = (Float_t)centrality->GetCentralityPercentileUnchecked("TKL");
            } else {
              Int_t runnum = aodEvent->GetRunNumber();
              for (Int_t ir = 0; ir < 5; ir++) {
                if (runnum == selRun[ir]) {
                  isSelRun = kTRUE;
                  break;
                }
              }
              if ((quality == 8 || quality == 9) && isSelRun)cent = (Float_t)centrality->GetCentralityPercentileUnchecked("TKL");
            }
          }
        } else {
          if (estimator == kCentCL1) {
            cent = (Float_t)(centrality->GetCentralityPercentile("CL1"));
            if (cent < 0) {
              Int_t quality = centrality->GetQuality();
              if (quality <= 1) {
                cent = (Float_t)centrality->GetCentralityPercentileUnchecked("CL1");
              } else {
                Int_t runnum = aodEvent->GetRunNumber();
                for (Int_t ir = 0; ir < 5; ir++) {
                  if (runnum == selRun[ir]) {
                    isSelRun = kTRUE;
                    break;
                  }
                }
                if ((quality == 8 || quality == 9) && isSelRun)cent = (Float_t)centrality->GetCentralityPercentileUnchecked("CL1");
              }
            }
          } else {
            AliWarning("Centrality estimator not valid");

          }
        }
      }
    }
  }
  return cent;
}
//-------------------------------------------------------------------
Bool_t AliRDHFJetsCuts::CompareCuts(const AliRDHFJetsCuts* obj) const
{
  //
  // Compare two cuts objects
  //

  Bool_t areEqual = kTRUE;

  if (fMinVtxType != obj->fMinVtxType) { printf("Minimum vtx type %d  %d\n", fMinVtxType, obj->fMinVtxType); areEqual = kFALSE;}

  if (fMinVtxContr != obj->fMinVtxContr) { printf("Minimum vtx contr %d  %d\n", fMinVtxContr, obj->fMinVtxContr); areEqual = kFALSE;}

  if (TMath::Abs(fMaxVtxRedChi2 - obj->fMaxVtxRedChi2) > 1.e-10) {   printf("Max vtx red chi2 %f  %f\n", fMaxVtxRedChi2, obj->fMaxVtxRedChi2); areEqual = kFALSE;}

  if (fMinSPDMultiplicity != obj->fMinSPDMultiplicity) {  printf("Min SPD mult %d\n  %d", fMinSPDMultiplicity, obj->fMinSPDMultiplicity); areEqual = kFALSE;}



  if (fRemoveDaughtersFromPrimary != obj->fRemoveDaughtersFromPrimary) {printf("Remove daughters from vtx %d  %d\n", (Int_t)fRemoveDaughtersFromPrimary, (Int_t)obj->fRemoveDaughtersFromPrimary); areEqual = kFALSE;}
  if (fTrackCuts) {
    if (fTrackCuts->GetMinNClusterTPC() != obj->fTrackCuts->GetMinNClusterTPC()) {printf("MinNClsTPC %d  %d\n", fTrackCuts->GetMinNClusterTPC(), obj->fTrackCuts->GetMinNClusterTPC()); areEqual = kFALSE;}

    if (fTrackCuts->GetMinNClustersITS() != obj->fTrackCuts->GetMinNClustersITS()) {printf("MinNClsITS %d  %d\n", fTrackCuts->GetMinNClustersITS(), obj->fTrackCuts->GetMinNClustersITS()); areEqual = kFALSE;}

    if (TMath::Abs(fTrackCuts->GetMaxChi2PerClusterTPC() - obj->fTrackCuts->GetMaxChi2PerClusterTPC()) > 1.e-10) {printf("MaxChi2ClsTPC %f  %f\n", fTrackCuts->GetMaxChi2PerClusterTPC(), obj->fTrackCuts->GetMaxChi2PerClusterTPC()); areEqual = kFALSE;}

    if (fTrackCuts->GetClusterRequirementITS(AliESDtrackCuts::kSPD) != obj->fTrackCuts->GetClusterRequirementITS(AliESDtrackCuts::kSPD)) {printf("ClusterReq SPD %d  %d\n", fTrackCuts->GetClusterRequirementITS(AliESDtrackCuts::kSPD), obj->fTrackCuts->GetClusterRequirementITS(AliESDtrackCuts::kSPD)); areEqual = kFALSE;}
  }


  return areEqual;
}

//--------------------------------------------------------------------------
Bool_t AliRDHFJetsCuts::RecalcOwnPrimaryVtx(AliAODRecoDecayHF* d,
                                            AliAODEvent* aod) const
{
  //
  // Recalculate primary vertex without daughters
  //

  if (!aod) {
    AliError("Can not remove daughters from vertex without AOD event");
    return 0;
  }

  AliAODVertex* recvtx = d->RemoveDaughtersFromPrimaryVtx(aod);
  if (!recvtx) {
    AliDebug(2, "Removal of daughter tracks failed");
    return kFALSE;
  }


  //set recalculed primary vertex
  d->SetOwnPrimaryVtx(recvtx);
  delete recvtx;

  return kTRUE;
}
//--------------------------------------------------------------------------
Bool_t AliRDHFJetsCuts::SetMCPrimaryVtx(AliAODRecoDecayHF* d, AliAODEvent* aod) const
{
  //
  // Recalculate primary vertex without daughters
  //

  if (!aod) {
    AliError("Can not get MC vertex without AOD event");
    return kFALSE;
  }

  // load MC header
  AliAODMCHeader* mcHeader =
    (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
  if (!mcHeader) {
    AliError("Can not get MC vertex without AODMCHeader event");
    return kFALSE;
  }
  Double_t pos[3];
  Double_t covmatrix[6] = {0., 0., 0., 0., 0., 0.};
  mcHeader->GetVertex(pos);
  AliAODVertex* recvtx = new AliAODVertex(pos, covmatrix);

  if (!recvtx) {
    AliDebug(2, "Removal of daughter tracks failed");
    return kFALSE;
  }

  //set recalculed primary vertex
  d->SetOwnPrimaryVtx(recvtx);

  d->RecalculateImpPars(recvtx, aod);

  delete recvtx;

  return kTRUE;
}
//--------------------------------------------------------------------------
void AliRDHFJetsCuts::CleanOwnPrimaryVtx(AliAODRecoDecayHF* d,
                                         AliAODEvent* aod,
                                         AliAODVertex* origownvtx) const
{
  //
  // Clean-up own primary vertex if needed
  //

  if (fRemoveDaughtersFromPrimary || fUseMCVertex) {
    d->UnsetOwnPrimaryVtx();
    if (origownvtx) {
      d->SetOwnPrimaryVtx(origownvtx);
      delete origownvtx; origownvtx = NULL;
    }
    d->RecalculateImpPars(d->GetPrimaryVtx(), aod);
  } else {
    if (origownvtx) {
      delete origownvtx; origownvtx = NULL;
    }
  }
  return;
}
//--------------------------------------------------------------------------
Bool_t AliRDHFJetsCuts::IsSignalMC(AliAODRecoDecay* d, AliAODEvent* aod, Int_t pdg) const
{
  //
  // Checks if this candidate is matched to MC signal
  //

  if (!aod) return kFALSE;

  // get the MC array
  TClonesArray* mcArray = (TClonesArray*)((AliAODEvent*)aod)->GetList()->FindObject(AliAODMCParticle::StdBranchName());

  if (!mcArray) return kFALSE;

  // try to match
  Int_t label = d->MatchToMC(pdg, mcArray);

  if (label >= 0) {
    //printf("MATCH!\n");
    return kTRUE;
  }

  return kFALSE;
}


//--------------------------------------------------------------------------
Bool_t AliRDHFJetsCuts::RecomputePrimaryVertex(AliAODEvent* event) const
{
  // recompute event primary vertex from AOD tracks

  AliVertexerTracks* vertexer = new AliVertexerTracks(event->GetMagneticField());
  vertexer->SetITSMode();
  vertexer->SetMinClusters(3);

  AliAODVertex* pvtx = event->GetPrimaryVertex();
  if (strstr(pvtx->GetTitle(), "VertexerTracksWithConstraint")) {
    Float_t diamondcovxy[3];
    event->GetDiamondCovXY(diamondcovxy);
    Double_t pos[3] = {event->GetDiamondX(), event->GetDiamondY(), 0.};
    Double_t cov[6] = {diamondcovxy[0], diamondcovxy[1], diamondcovxy[2], 0., 0., 10.*10.};
    AliESDVertex* diamond = new AliESDVertex(pos, cov, 1., 1);
    vertexer->SetVtxStart(diamond);
    delete diamond; diamond = NULL;
  }

  AliESDVertex* vertexESD = (AliESDVertex*)vertexer->FindPrimaryVertex(event);
  if (!vertexESD) return kFALSE;
  if (vertexESD->GetNContributors() <= 0) {
    //AliDebug(2,"vertexing failed");
    delete vertexESD; vertexESD = NULL;
    return kFALSE;
  }
  delete vertexer; vertexer = NULL;

  // convert to AliAODVertex
  Double_t pos[3], cov[6], chi2perNDF;
  vertexESD->GetXYZ(pos); // position
  vertexESD->GetCovMatrix(cov); //covariance matrix
  chi2perNDF = vertexESD->GetChi2toNDF();
  delete vertexESD; vertexESD = NULL;

  pvtx->SetPosition(pos[0], pos[1], pos[2]);
  pvtx->SetChi2perNDF(chi2perNDF);
  pvtx->SetCovMatrix(cov);

  return kTRUE;
}
