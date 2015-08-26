/**************************************************************************
 * Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
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
#include <vector>
#include <map>

#include <TArrayD.h>
#include <TMath.h>
#include <THashList.h>
#include <TString.h>

#include "AliAnalysisUtils.h"
#include "AliAODTrack.h"
#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliInputEventHandler.h"
#include "AliVVertex.h"

#include "AliEmcalTriggerPatchInfo.h"
#include "AliEMCalHistoContainer.h"
#include "AliAnalysisTaskChargedParticlesRef.h"

/// \cond CLASSIMP
ClassImp(EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRef)
/// \endcond

namespace EMCalTriggerPtAnalysis {

/**
 * Dummy constructor
 */
AliAnalysisTaskChargedParticlesRef::AliAnalysisTaskChargedParticlesRef() :
    AliAnalysisTaskSE(),
    fTrackCuts(NULL),
    fAnalysisUtil(NULL),
    fHistos(NULL),
    fTriggerStringFromPatches(kFALSE),
    fYshift(0.465),
    fEtaSign(1)
{
}

/**
 * Main constructor
 * @param name Name of the task
 */
AliAnalysisTaskChargedParticlesRef::AliAnalysisTaskChargedParticlesRef(const char *name) :
    AliAnalysisTaskSE(name),
    fTrackCuts(NULL),
    fAnalysisUtil(NULL),
    fHistos(NULL),
    fTriggerStringFromPatches(kFALSE),
    fYshift(0.465),
    fEtaSign(1)
{
  DefineOutput(1, TList::Class());
}

/**
 * Destuctor
 */
AliAnalysisTaskChargedParticlesRef::~AliAnalysisTaskChargedParticlesRef() {
  if(fTrackCuts) delete fTrackCuts;
  if(fAnalysisUtil) delete fAnalysisUtil;
  if(fHistos) delete fHistos;
}

/**
 * Create the output histograms
 */
void AliAnalysisTaskChargedParticlesRef::UserCreateOutputObjects() {
  fAnalysisUtil = new AliAnalysisUtils;

  fTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(true, 1);
  fTrackCuts->SetName("Standard Track cuts");
  fTrackCuts->SetMinNCrossedRowsTPC(120);
  fTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");

  TArrayD oldbinning, newbinning;
  CreateOldPtBinning(oldbinning);
  CreateNewPtBinning(newbinning);

  fHistos = new AliEMCalHistoContainer("Ref");
  // Exclusive means classes without higher trigger classes:
  // EG2excl means EG2 && !EG1
  // EJ2excl means EJ2 && !EJ1
  // MBExcl means MinBias && !EMCAL trigger
  TString triggers[8] = {"MB", "EJ1", "EJ2", "EG1", "EG2", "EG2excl", "EJ2excl", "MBexcl"};
  Double_t ptcuts[5] = {1., 2., 5., 10., 20.};
  for(TString *trg = triggers; trg < triggers+8; trg++){
    fHistos->CreateTH1(Form("hEventCount%s", trg->Data()), Form("Event Counter for trigger class %s", trg->Data()), 1, 0.5, 1.5);
    fHistos->CreateTH1(Form("hVertexBefore%s", trg->Data()), Form("Vertex distribution before z-cut for trigger class %s", trg->Data()), 500, -50, 50);
    fHistos->CreateTH1(Form("hVertexAfter%s", trg->Data()), Form("Vertex distribution after z-cut for trigger class %s", trg->Data()), 100, -10, 10);
    fHistos->CreateTH1(Form("hPtEtaAllOldBinning%s", trg->Data()), Form("Charged particle pt distribution all eta old binning trigger %s", trg->Data()), oldbinning);
    fHistos->CreateTH1(Form("hPtEtaCentOldBinning%s", trg->Data()), Form("Charged particle pt distribution central eta old binning trigger %s", trg->Data()), oldbinning);
    fHistos->CreateTH1(Form("hPtEtaAllNewBinning%s", trg->Data()), Form("Charged particle pt distribution all eta new binning trigger %s", trg->Data()), newbinning);
    fHistos->CreateTH1(Form("hPtEtaCentNewBinning%s", trg->Data()), Form("Charged particle pt distribution central eta new binning trigger %s", trg->Data()), newbinning);
    fHistos->CreateTH1(Form("hPtEMCALEtaAllOldBinning%s", trg->Data()), Form("Charged particle in EMCAL pt distribution all eta old binning trigger %s", trg->Data()), oldbinning);
    fHistos->CreateTH1(Form("hPtEMCALEtaCentOldBinning%s", trg->Data()), Form("Charged particle in EMCAL pt distribution central eta old binning trigger %s", trg->Data()), oldbinning);
    fHistos->CreateTH1(Form("hPtEMCALEtaAllNewBinning%s", trg->Data()), Form("Charged particle in EMCAL pt distribution all eta new binning trigger %s", trg->Data()), newbinning);
    fHistos->CreateTH1(Form("hPtEMCALEtaCentNewBinning%s", trg->Data()), Form("Charged particle in EMCAL pt distribution central eta new binning trigger %s", trg->Data()), newbinning);
    for(int ipt = 0; ipt < 5; ipt++){
      fHistos->CreateTH1(
          Form("hEtaLabDistAllPt1%s", trg->Data()),
          Form("Eta (lab) distribution without etacut for tracks with Pt above 1 GeV/c trigger %s", trg->Data()),
          100,
          -1.,
          1.
          );
      fHistos->CreateTH1(
          Form("hEtaLabDistCutPt1%s", trg->Data()),
          Form("Eta (lab) distribution with etacut for tracks with Pt above 1 GeV/c trigger %s", trg->Data()),
          100,
          -1.,
          1.
          );
      fHistos->CreateTH1(
          Form("hEtaCentDistAllPt1%s", trg->Data()),
          Form("Eta (cent) distribution without etacut for tracks with Pt above 1 GeV/c trigger %s",
              trg->Data()),
              160,
              -1.3,
              1.3
              );
      fHistos->CreateTH1(
          Form("hEtaCentDistCutPt1%s", trg->Data()),
          Form("Eta (cent) distribution with etacut for tracks with Pt above 1 GeV/c trigger %s", trg->Data()),
          160,
          -1.3,
          1.3
          );
      fHistos->CreateTH1(
          Form("hEtaLabDistAllEMCALPt%d%s", static_cast<Int_t>(ptcuts[ipt]), trg->Data()),
          Form("Eta (lab) distribution without etacut for tracks in EMCAL with Pt above %.1f GeV/c trigger %s", ptcuts[ipt], trg->Data()),
          100,
          -1.,
          1.
          );
      fHistos->CreateTH1(
          Form("hEtaLabDistCutEMCALPt%d%s", static_cast<Int_t>(ptcuts[ipt]), trg->Data()),
          Form("Eta (lab) distribution with etacut for tracks in EMCAL with Pt above %.1f GeV/c trigger %s", ptcuts[ipt], trg->Data()),
          100,
          -1.,
          1.
          );
      fHistos->CreateTH1(
          Form("hEtaCentDistAllEMCALPt%d%s", static_cast<Int_t>(ptcuts[ipt]), trg->Data()),
          Form("Eta (cent) distribution without etacut for tracks in EMCAL with Pt above %.1f GeV/c trigger %s", ptcuts[ipt], trg->Data()),
          160,
          -1.3,
          1.3
          );
      fHistos->CreateTH1(
          Form("hEtaCentDistCutEMCALPt%d%s", static_cast<Int_t>(ptcuts[ipt]), trg->Data()),
          Form("Eta (cent) distribution with etacut for tracks in EMCAL with Pt above %.1f GeV/c trigger %s", ptcuts[ipt], trg->Data()),
          160,
          -1.3,
          1.3
          );
      fHistos->CreateTH1(
          Form("hPhiDistAllPt%d%s", static_cast<Int_t>(ptcuts[ipt]), trg->Data()),
          Form("#phi distribution of particles with Pt above %.1f GeV/c trigger %s", ptcuts[ipt], trg->Data()),
          300,
          0.,
          2*TMath::Pi()
          );
    }
  }
  PostData(1, fHistos->GetListOfHistograms());
}

/**
 * Simple unit test framework
 * - Select event using AliAnalysisUtil
 * - Assing trigger type (Request INT7, EJ*, EG*)
 * - Loop over tracks, select particles
 * - Fill distributions
 * @param option Not used
 */
void AliAnalysisTaskChargedParticlesRef::UserExec(Option_t*) {
  // Select event
  TString triggerstring = "";
  if(fTriggerStringFromPatches){
    triggerstring = GetFiredTriggerClassesFromPatches(dynamic_cast<TClonesArray *>(fInputEvent->FindListObject("EmcalTriggers")));
  } else {
    triggerstring = fInputEvent->GetFiredTriggerClasses();
  }
  Bool_t isMinBias = fInputHandler->IsEventSelected() & AliVEvent::kINT7,
      isEJ1 = triggerstring.Contains("EJ1"),
      isEJ2 = triggerstring.Contains("EJ2"),
      isEG1 = triggerstring.Contains("EG1"),
      isEG2 = triggerstring.Contains("EG2");
  if(!(isMinBias || isEG1 || isEG2 || isEJ1 || isEJ2)) return;
  const AliVVertex *vtx = fInputEvent->GetPrimaryVertex();
  //if(!fInputEvent->IsPileupFromSPD(3, 0.8, 3., 2., 5.)) return;         // reject pileup event
  if(vtx->GetNContributors() < 1) return;
  // Fill reference distribution for the primary vertex before any z-cut
  if(isMinBias) fHistos->FillTH1("hVertexBeforeMB", vtx->GetZ());
  if(isEJ1) fHistos->FillTH1("hVertexBeforeEJ1", vtx->GetZ());
  if(isEJ2) fHistos->FillTH1("hVertexBeforeEJ2", vtx->GetZ());
  if(isEG1) fHistos->FillTH1("hVertexBeforeEG1", vtx->GetZ());
  if(isEG2) fHistos->FillTH1("hVertexBeforeEG2", vtx->GetZ());
  if(!fAnalysisUtil->IsVertexSelected2013pA(fInputEvent)) return;       // Apply new vertex cut
  if(fAnalysisUtil->IsPileUpEvent(fInputEvent)) return;       // Apply new vertex cut
  // Apply vertex z cut
  if(vtx->GetZ() < -10. || vtx->GetZ() > 10.) return;

  // Fill Event counter and reference vertex distributions for the different trigger classes
  if(isMinBias){
    fHistos->FillTH1("hEventCountMB", 1);
    fHistos->FillTH1("hVertexAfterMB", vtx->GetZ());
    // Check for exclusive classes
    if(!(isEG1 || isEG2 || isEJ1 || isEJ2)){
      fHistos->FillTH1("hEventCountMBexcl", 1);
      fHistos->FillTH1("hVertexAfterMBexcl", vtx->GetZ());
    }
  }
  if(isEJ1){
    fHistos->FillTH1("hEventCountEJ1", 1);
    fHistos->FillTH1("hVertexAfterEJ1", vtx->GetZ());
  }
  if(isEJ2){
    fHistos->FillTH1("hEventCountEJ2", 1);
    fHistos->FillTH1("hVertexAfterEJ2", vtx->GetZ());
    // Check for exclusive classes
    if(!isEJ1){
      fHistos->FillTH1("hEventCountEJ2excl", 1);
      fHistos->FillTH1("hVertexAfterEJ2excl", vtx->GetZ());
    }
  }
  if(isEG1){
    fHistos->FillTH1("hEventCountEG1", 1);
    fHistos->FillTH1("hVertexAfterEG1", vtx->GetZ());
  }
  if(isEG2){
    fHistos->FillTH1("hEventCountEG2", 1);
    fHistos->FillTH1("hVertexAfterEG2", vtx->GetZ());
    // Check for exclusive classes
    if(!isEG1){
      fHistos->FillTH1("hEventCountEG2excl", 1);
      fHistos->FillTH1("hVertexAfterEG2excl", vtx->GetZ());
    }
  }

  // Loop over tracks, fill select particles
  // Histograms
  // - Full eta_{lab} (-0.8, 0.8), new binning
  // - Full eta_{lab} (-0.8, 0.8), old binning
  // - Eta distribution for tracks above 1, 2, 5, 10 GeV/c without eta cut
  // - Central eta_{cms} (-0.3, -0.2), new binning,
  // - Central eta_{cms} (-0.8, -0.2), old binning,
  // - Eta distribution for tracks above 1, 2, 5, 10 GeV/c
  // - Eta distribution for tracks above 1, 2, 5, 10 GeV/c with eta cut
  AliVTrack *checktrack(NULL);
  int ptmin[5] = {1,2,5,10,20}; // for eta distributions
  Bool_t isEMCAL(kFALSE);
  for(int itrk = 0; itrk < fInputEvent->GetNumberOfTracks(); ++itrk){
    checktrack = dynamic_cast<AliVTrack *>(fInputEvent->GetTrack(itrk));
    if(!checktrack) continue;
    if(TMath::Abs(checktrack->Eta()) > 0.8) continue;
    if(TMath::Abs(checktrack->Pt()) < 0.1) continue;
    isEMCAL = checktrack->Eta() > 1.5 && checktrack->Eta() < 3.1;

    // Calculate eta in cms frame according
    // EPJC74 (2014) 3054:
    // eta_cms = - eta_lab - |yshift|
    Double_t etacent = -1. * checktrack->Eta() - TMath::Abs(fYshift);
    etacent *= fEtaSign;

    // Distinguish track selection for ESD and AOD tracks
    AliESDtrack *esdtrack(NULL);
    AliAODTrack *aodtrack(NULL);
    if((esdtrack = dynamic_cast<AliESDtrack *>(checktrack))){
      if(!TrackSelectionESD(esdtrack)) continue;
    } else if((aodtrack = dynamic_cast<AliAODTrack *>(checktrack))){
      if(!TrackSelectionAOD(aodtrack)) continue;
    } else {
      continue;
    }

    // fill histograms allEta
    if(isMinBias){
      FillTrackHistos("MB", checktrack->Pt(), checktrack->Eta() * fEtaSign, checktrack->Phi(), etacent, TMath::Abs(etacent) < 0.3, isEMCAL);
      // check for exclusive classes
      if(!(isEG1 || isEG2 || isEJ1 || isEJ2)){
        FillTrackHistos("MBexcl", checktrack->Pt(), checktrack->Eta() * fEtaSign, checktrack->Phi(), etacent, TMath::Abs(etacent) < 0.3, isEMCAL);
      }
    }
    if(isEJ1){
      FillTrackHistos("EJ1", checktrack->Pt(), checktrack->Eta() * fEtaSign, checktrack->Phi(), etacent, TMath::Abs(etacent) < 0.3, isEMCAL);
    }
    if(isEJ2){
      FillTrackHistos("EJ2", checktrack->Pt(), checktrack->Eta() * fEtaSign, checktrack->Phi(), etacent, TMath::Abs(etacent) < 0.3, isEMCAL);
      // check for exclusive classes
      if(!isEJ1){
        FillTrackHistos("EJ2excl", checktrack->Pt(), checktrack->Eta() * fEtaSign, checktrack->Phi(), etacent, TMath::Abs(etacent) < 0.3, isEMCAL);
      }
    }
    if(isEG1){
      FillTrackHistos("EG1", checktrack->Pt(), checktrack->Eta() * fEtaSign, checktrack->Phi(), etacent, TMath::Abs(etacent) < 0.3, isEMCAL);
    }
    if(isEG2){
      FillTrackHistos("EG2", checktrack->Pt(), checktrack->Eta() * fEtaSign, checktrack->Phi(), etacent, TMath::Abs(etacent) < 0.3, isEMCAL);
      // check for exclusive classes
      if(!isEG1){
        FillTrackHistos("EG2excl", checktrack->Pt(), checktrack->Eta() * fEtaSign, checktrack->Phi(), etacent, TMath::Abs(etacent) < 0.3, isEMCAL);
      }
    }
  }
  PostData(1, fHistos->GetListOfHistograms());
}

/**
 * Fill track histograms
 * @param eventclass Trigger class fired
 * @param pt track \f$ p_{t} \f$
 * @param etalab Track \f$ \eta \f$ in lab frame
 * @param etacent Track \f$ \eta \f$ in cms frame
 * @param phi Track \f$ \eta \f$ in lab frame
 * @param etacut Track accepted by \f$ \eta \f$ cut
 * @param inEmcal Track in EMCAL \f$ \phi \f$ acceptance
 */
void AliAnalysisTaskChargedParticlesRef::FillTrackHistos(
    const char *eventclass,
    Double_t pt,
    Double_t etalab,
    Double_t etacent,
    Double_t phi,
    Bool_t etacut,
    Bool_t inEmcal
    )
{
  fHistos->FillTH1(Form("hPtEtaAllNewBinning%s", eventclass), TMath::Abs(pt));
  fHistos->FillTH1(Form("hPtEtaAllOldBinning%s", eventclass), TMath::Abs(pt));
  if(inEmcal){
    fHistos->FillTH1(Form("hPtEMCALEtaAllNewBinning%s", eventclass), TMath::Abs(pt));
    fHistos->FillTH1(Form("hPtEMCALEtaAllOldBinning%s", eventclass), TMath::Abs(pt));
  }

  int ptmin[5] = {1,2,5,10,20}; // for eta distributions
  for(int icut = 0; icut < 5; icut++){
    if(TMath::Abs(pt) > static_cast<double>(ptmin[icut])){
      fHistos->FillTH1(Form("hPhiDistAllPt%d%s", ptmin[icut], eventclass), phi);
      fHistos->FillTH1(Form("hEtaLabDistAllPt%d%s", ptmin[icut], eventclass), etalab);
      fHistos->FillTH1(Form("hEtaCentDistAllPt%d%s", ptmin[icut], eventclass), etacent);
      if(inEmcal){
        fHistos->FillTH1(Form("hEtaLabDistAllEMCALPt%d%s", ptmin[icut], eventclass), etalab);
        fHistos->FillTH1(Form("hEtaCentDistAllEMCALPt%d%s", ptmin[icut], eventclass), etacent);
      }
    }
  }

  if(etacut){
    fHistos->FillTH1(Form("hPtEtaCentNewBinning%s", eventclass), TMath::Abs(pt));
    fHistos->FillTH1(Form("hPtEtaCentOldBinning%s", eventclass), TMath::Abs(pt));
    if(inEmcal){
      fHistos->FillTH1(Form("hPtEMCALEtaCentNewBinning%s", eventclass), TMath::Abs(pt));
      fHistos->FillTH1(Form("hPtEMCALEtaCentOldBinning%s", eventclass), TMath::Abs(pt));
    }
    for(int icut = 0; icut < 5; icut++){
      if(TMath::Abs(pt) > static_cast<double>(ptmin[icut])){
        fHistos->FillTH1(Form("hEtaLabDistCutPt%d%s", ptmin[icut], eventclass), etalab);
        fHistos->FillTH1(Form("hEtaCentDistCutPt%d%s", ptmin[icut], eventclass), etacent);
        if(inEmcal){
          fHistos->FillTH1(Form("hEtaLabDistCutEMCALPt%d%s", ptmin[icut], eventclass), etalab);
          fHistos->FillTH1(Form("hEtaCentDistCutEMCALPt%d%s", ptmin[icut], eventclass), etacent);
        }
      }
    }
  }
}


/**
 * Create old pt binning
 * @param binning
 */
void AliAnalysisTaskChargedParticlesRef::CreateOldPtBinning(TArrayD &binning) const{
 std::vector<double> mybinning;
 std::map<double,double> definitions;
 definitions.insert(std::pair<double,double>(2.5, 0.1));
 definitions.insert(std::pair<double,double>(7., 0.25));
 definitions.insert(std::pair<double,double>(15., 0.5));
 definitions.insert(std::pair<double,double>(25., 1.));
 definitions.insert(std::pair<double,double>(40., 2.5));
 definitions.insert(std::pair<double,double>(50., 5.));
 definitions.insert(std::pair<double,double>(100., 10.));
 double currentval = 0;
 for(std::map<double,double>::iterator id = definitions.begin(); id != definitions.end(); ++id){
   double limit = id->first, binwidth = id->second;
   while(currentval < limit){
     currentval += binwidth;
     mybinning.push_back(currentval);
    }
  }
  binning.Set(mybinning.size());
  int ib = 0;
  for(std::vector<double>::iterator it = mybinning.begin(); it != mybinning.end(); ++it)
    binning[ib++] = *it;
}

/**
 * Create new Pt binning
 * @param binning
 */
void AliAnalysisTaskChargedParticlesRef::CreateNewPtBinning(TArrayD& binning) const {
  std::vector<double> mybinning;
  std::map<double,double> definitions;
  definitions.insert(std::pair<double, double>(1, 0.05));
  definitions.insert(std::pair<double, double>(2, 0.1));
  definitions.insert(std::pair<double, double>(4, 0.2));
  definitions.insert(std::pair<double, double>(7, 0.5));
  definitions.insert(std::pair<double, double>(16, 1));
  definitions.insert(std::pair<double, double>(36, 2));
  definitions.insert(std::pair<double, double>(40, 4));
  definitions.insert(std::pair<double, double>(50, 5));
  definitions.insert(std::pair<double, double>(100, 10));
  definitions.insert(std::pair<double, double>(200, 20));
  double currentval = 0.;
  mybinning.push_back(currentval);
  for(std::map<double,double>::iterator id = definitions.begin(); id != definitions.end(); ++id){
    double limit = id->first, binwidth = id->second;
    while(currentval < limit){
      currentval += binwidth;
      mybinning.push_back(currentval);
    }
  }
  binning.Set(mybinning.size());
  int ib = 0;
  for(std::vector<double>::iterator it = mybinning.begin(); it != mybinning.end(); ++it)
    binning[ib++] = *it;
}

/**
 * Run track selection for ESD tracks
 * @param track The track to check
 * @return True if the track is selected, false otherwise
 */
Bool_t AliAnalysisTaskChargedParticlesRef::TrackSelectionESD(AliESDtrack* track) {
  return fTrackCuts->AcceptTrack(track);
}

/**
 * Run track selection for AOD tracks
 * @param track The track to check
 * @return True if the track is selected, false otherwise
 */
Bool_t AliAnalysisTaskChargedParticlesRef::TrackSelectionAOD(AliAODTrack* track) {
  if(!track->TestFilterBit(AliAODTrack::kTrkGlobal)) return false;
  if(track->GetTPCNCrossedRows() < 120) return false;
  return true;
}

/**
 * Apply trigger selection using offline patches and trigger thresholds based on offline ADC Amplitude
 * @param triggerpatches Trigger patches found by the trigger maker
 * @return String with EMCAL trigger decision
 */
TString AliAnalysisTaskChargedParticlesRef::GetFiredTriggerClassesFromPatches(const TClonesArray* triggerpatches) const {
  TString triggerstring = "";
  Int_t nEJ1 = 0, nEJ2 = 0, nEG1 = 0, nEG2 = 0;
  double  minADC_EJ1 = 260.,
          minADC_EJ2 = 127.,
          minADC_EG1 = 140.,
          minADC_EG2 = 89.;
  for(TIter patchIter = TIter(triggerpatches).Begin(); patchIter != TIter::End(); ++patchIter){
    AliEmcalTriggerPatchInfo *patch = dynamic_cast<AliEmcalTriggerPatchInfo *>(*patchIter);
    if(!patch->IsOfflineSimple()) continue;
    if(patch->IsJetHighSimple() && patch->GetADCOfflineAmp() > minADC_EJ1) nEJ1++;
    if(patch->IsJetLowSimple() && patch->GetADCOfflineAmp() > minADC_EJ2) nEJ2++;
    if(patch->IsGammaHighSimple() && patch->GetADCOfflineAmp() > minADC_EG1) nEG1++;
    if(patch->IsGammaLowSimple() && patch->GetADCOfflineAmp() > minADC_EG2) nEG2++;
  }
  if(nEJ1) triggerstring += "EJ1";
  if(nEJ2){
    if(triggerstring.Length()) triggerstring += ",";
    triggerstring += "EJ2";
  }
  if(nEG1){
    if(triggerstring.Length()) triggerstring += ",";
    triggerstring += "EG1";
  }
  if(nEG2){
    if(triggerstring.Length()) triggerstring += ",";
    triggerstring += "EG2";
  }
  return triggerstring;
}

} /* namespace EMCalTriggerPtAnalysis */
