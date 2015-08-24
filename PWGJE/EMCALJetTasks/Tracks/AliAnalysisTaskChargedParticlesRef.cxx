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
  for(TString *trg = triggers; trg < triggers+8; trg++){
    fHistos->CreateTH1(Form("hEventCount%s", trg->Data()), Form("Event Counter for trigger class %s", trg->Data()), 1, 0.5, 1.5);
    fHistos->CreateTH1(Form("hVertexBefore%s", trg->Data()), Form("Vertex distribution before z-cut for trigger class %s", trg->Data()), 500, -50, 50);
    fHistos->CreateTH1(Form("hVertexAfter%s", trg->Data()), Form("Vertex distribution after z-cut for trigger class %s", trg->Data()), 100, -10, 10);
    fHistos->CreateTH1(Form("hPtEtaAllOldBinning%s", trg->Data()), Form("Charged particle pt distribution all eta old binning trigger %s", trg->Data()), oldbinning);
    fHistos->CreateTH1(Form("hPtEtaCentOldBinning%s", trg->Data()), Form("Charged particle pt distribution central eta old binning trigger %s", trg->Data()), oldbinning);
    fHistos->CreateTH1(Form("hPtEtaAllNewBinning%s", trg->Data()), Form("Charged particle pt distribution all eta new binning trigger %s", trg->Data()), newbinning);
    fHistos->CreateTH1(Form("hPtEtaCentNewBinning%s", trg->Data()), Form("Charged particle pt distribution central eta new binning trigger %s", trg->Data()), newbinning);
    fHistos->CreateTH1(Form("hEtaLabDistAllPt1%s", trg->Data()), Form("Eta (lab) distribution without etacut for tracks with Pt above 1 GeV/c trigger %s", trg->Data()), 100, -1., 1.);
    fHistos->CreateTH1(Form("hEtaLabDistAllPt2%s", trg->Data()), Form("Eta (lab) distribution without etacut for tracks with Pt above 2 GeV/c trigger %s", trg->Data()), 100, -1., 1.);
    fHistos->CreateTH1(Form("hEtaLabDistAllPt5%s", trg->Data()), Form("Eta (lab) distribution without etacut for tracks with Pt above 5 GeV/c trigger %s", trg->Data()), 100, -1., 1.);
    fHistos->CreateTH1(Form("hEtaLabDistAllPt10%s", trg->Data()), Form("Eta (lab) distribution without etacut for tracks with Pt above 10 GeV/c trigger %s", trg->Data()), 100, -1., 1.);
    fHistos->CreateTH1(Form("hEtaLabDistAllPt20%s", trg->Data()), Form("Eta (lab) distribution without etacut for tracks with Pt above 20 GeV/c trigger %s", trg->Data()), 100, -1., 1.);
    fHistos->CreateTH1(Form("hEtaLabDistCutPt1%s", trg->Data()), Form("Eta (lab) distribution with etacut for tracks with Pt above 1 GeV/c trigger %s", trg->Data()), 100, -1., 1.);
    fHistos->CreateTH1(Form("hEtaLabDistCutPt2%s", trg->Data()), Form("Eta (lab) distribution with etacut for tracks with Pt above 2 GeV/c trigger %s", trg->Data()), 100, -1., 1.);
    fHistos->CreateTH1(Form("hEtaLabDistCutPt5%s", trg->Data()), Form("Eta (lab) distribution with etacut for tracks with Pt above 5 GeV/c trigger %s", trg->Data()), 100, -1., 1.);
    fHistos->CreateTH1(Form("hEtaLabDistCutPt10%s", trg->Data()), Form("Eta (lab) distribution with etacut for tracks with Pt above 10 GeV/c trigger %s", trg->Data()), 100, -1., 1.);
    fHistos->CreateTH1(Form("hEtaLabDistCutPt20%s", trg->Data()), Form("Eta (lab) distribution with etacut for tracks with Pt above 20 GeV/c trigger %s", trg->Data()), 100, -1., 1.);
    fHistos->CreateTH1(Form("hEtaCentDistAllPt1%s", trg->Data()), Form("Eta (cent) distribution without etacut for tracks with Pt above 1 GeV/c trigger %s", trg->Data()), 160, -1.3, 1.3);
    fHistos->CreateTH1(Form("hEtaCentDistAllPt2%s", trg->Data()), Form("Eta (cent) distribution without etacut for tracks with Pt above 2 GeV/c trigger %s", trg->Data()), 160, -1.3, 1.3);
    fHistos->CreateTH1(Form("hEtaCentDistAllPt5%s", trg->Data()), Form("Eta (cent) distribution without etacut for tracks with Pt above 5 GeV/c trigger %s", trg->Data()), 160, -1.3, 1.3);
    fHistos->CreateTH1(Form("hEtaCentDistAllPt10%s", trg->Data()), Form("Eta (cent) distribution without etacut for tracks with Pt above 10 GeV/c trigger %s", trg->Data()), 160, -1.3, 1.3);
    fHistos->CreateTH1(Form("hEtaCentDistAllPt20%s", trg->Data()), Form("Eta (cent) distribution without etacut for tracks with Pt above 20 GeV/c trigger %s", trg->Data()), 160, -1.3, 1.3);
    fHistos->CreateTH1(Form("hEtaCentDistCutPt1%s", trg->Data()), Form("Eta (cent) distribution with etacut for tracks with Pt above 1 GeV/c trigger %s", trg->Data()), 160, -1.3, 1.3);
    fHistos->CreateTH1(Form("hEtaCentDistCutPt2%s", trg->Data()), Form("Eta (cent) distribution with etacut for tracks with Pt above 2 GeV/c trigger %s", trg->Data()), 160, -1.3, 1.3);
    fHistos->CreateTH1(Form("hEtaCentDistCutPt5%s", trg->Data()), Form("Eta (cent) distribution with etacut for tracks with Pt above 5 GeV/c trigger %s", trg->Data()), 160, -1.3, 1.3);
    fHistos->CreateTH1(Form("hEtaCentDistCutPt10%s", trg->Data()), Form("Eta (cent) distribution with etacut for tracks with Pt above 10 GeV/c trigger %s", trg->Data()), 160, -1.3, 1.3);
    fHistos->CreateTH1(Form("hEtaCentDistCutPt20%s", trg->Data()), Form("Eta (cent) distribution with etacut for tracks with Pt above 20 GeV/c trigger %s", trg->Data()), 160, -1.3, 1.3);
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
  TString triggerstring = fInputEvent->GetFiredTriggerClasses();
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
  for(int itrk = 0; itrk < fInputEvent->GetNumberOfTracks(); ++itrk){
    checktrack = dynamic_cast<AliVTrack *>(fInputEvent->GetTrack(itrk));
    if(!checktrack) continue;
    if(TMath::Abs(checktrack->Eta()) > 0.8) continue;
    if(TMath::Abs(checktrack->Pt()) < 0.1) continue;

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
      fHistos->FillTH1("hPtEtaAllNewBinningMB", TMath::Abs(checktrack->Pt()));
      fHistos->FillTH1("hPtEtaAllOldBinningMB", TMath::Abs(checktrack->Pt()));
      for(int icut = 0; icut < 5; icut++){
        if(TMath::Abs(checktrack->Pt()) > static_cast<double>(ptmin[icut])){
          fHistos->FillTH1(Form("hEtaLabDistAllPt%dMB", ptmin[icut]), checktrack->Eta());
          fHistos->FillTH1(Form("hEtaCentDistAllPt%dMB", ptmin[icut]), etacent);
        }
      }
      // check for exclusive classes
      if(!(isEG1 || isEG2 || isEJ1 || isEJ2)){
        fHistos->FillTH1("hPtEtaAllNewBinningMBexcl", TMath::Abs(checktrack->Pt()));
        fHistos->FillTH1("hPtEtaAllOldBinningMBexcl", TMath::Abs(checktrack->Pt()));
        for(int icut = 0; icut < 5; icut++){
          if(TMath::Abs(checktrack->Pt()) > static_cast<double>(ptmin[icut])){
            fHistos->FillTH1(Form("hEtaLabDistAllPt%dMBexcl", ptmin[icut]), checktrack->Eta());
            fHistos->FillTH1(Form("hEtaCentDistAllPt%dMBexcl", ptmin[icut]), etacent);
          }
        }
      }
    }
    if(isEJ1){
      fHistos->FillTH1("hPtEtaAllNewBinningEJ1", TMath::Abs(checktrack->Pt()));
      fHistos->FillTH1("hPtEtaAllOldBinningEJ1", TMath::Abs(checktrack->Pt()));
      for(int icut = 0; icut < 5; icut++){
        if(TMath::Abs(checktrack->Pt()) > static_cast<double>(ptmin[icut])){
          fHistos->FillTH1(Form("hEtaLabDistAllPt%dEJ1", ptmin[icut]), checktrack->Eta());
          fHistos->FillTH1(Form("hEtaCentDistAllPt%dEJ1", ptmin[icut]), etacent);
        }
      }
    }
    if(isEJ2){
      fHistos->FillTH1("hPtEtaAllNewBinningEJ2", TMath::Abs(checktrack->Pt()));
      fHistos->FillTH1("hPtEtaAllOldBinningEJ2", TMath::Abs(checktrack->Pt()));
      for(int icut = 0; icut < 5; icut++){
        if(TMath::Abs(checktrack->Pt()) > static_cast<double>(ptmin[icut])){
          fHistos->FillTH1(Form("hEtaLabDistAllPt%dEJ2", ptmin[icut]), checktrack->Eta());
          fHistos->FillTH1(Form("hEtaCentDistAllPt%dEJ2", ptmin[icut]), etacent);
        }
      }
      // check for exclusive classes
      if(!isEJ1){
        fHistos->FillTH1("hPtEtaAllNewBinningEJ2excl", TMath::Abs(checktrack->Pt()));
        fHistos->FillTH1("hPtEtaAllOldBinningEJ2excl", TMath::Abs(checktrack->Pt()));
        for(int icut = 0; icut < 5; icut++){
          if(TMath::Abs(checktrack->Pt()) > static_cast<double>(ptmin[icut])){
            fHistos->FillTH1(Form("hEtaLabDistAllPt%dEJ2excl", ptmin[icut]), checktrack->Eta());
            fHistos->FillTH1(Form("hEtaCentDistAllPt%dEJ2excl", ptmin[icut]), etacent);
          }
        }
      }
    }
    if(isEG1){
      fHistos->FillTH1("hPtEtaAllNewBinningEG1", TMath::Abs(checktrack->Pt()));
      fHistos->FillTH1("hPtEtaAllOldBinningEG1", TMath::Abs(checktrack->Pt()));
      for(int icut = 0; icut < 5; icut++){
        if(TMath::Abs(checktrack->Pt()) > static_cast<double>(ptmin[icut])){
          fHistos->FillTH1(Form("hEtaLabDistAllPt%dEG1", ptmin[icut]), checktrack->Eta());
          fHistos->FillTH1(Form("hEtaCentDistAllPt%dEG1", ptmin[icut]), etacent);
        }
      }
    }
    if(isEG2){
      fHistos->FillTH1("hPtEtaAllNewBinningEG2", TMath::Abs(checktrack->Pt()));
      fHistos->FillTH1("hPtEtaAllOldBinningEG2", TMath::Abs(checktrack->Pt()));
      for(int icut = 0; icut < 5; icut++){
        if(TMath::Abs(checktrack->Pt()) > static_cast<double>(ptmin[icut])){
          fHistos->FillTH1(Form("hEtaLabDistAllPt%dEG2", ptmin[icut]), checktrack->Eta());
          fHistos->FillTH1(Form("hEtaCentDistAllPt%dEG2", ptmin[icut]), etacent);
        }
      }
      // check for exclusive classes
      if(!isEG1){
        fHistos->FillTH1("hPtEtaAllNewBinningEG2excl", TMath::Abs(checktrack->Pt()));
        fHistos->FillTH1("hPtEtaAllOldBinningEG2excl", TMath::Abs(checktrack->Pt()));
        for(int icut = 0; icut < 5; icut++){
          if(TMath::Abs(checktrack->Pt()) > static_cast<double>(ptmin[icut])){
            fHistos->FillTH1(Form("hEtaLabDistAllPt%dEG2excl", ptmin[icut]), checktrack->Eta());
            fHistos->FillTH1(Form("hEtaCentDistAllPt%dEG2excl", ptmin[icut]), etacent);
          }
        }
      }
    }

    if(TMath::Abs(etacent) < 0.3){
      // Fill Histograms in central eta
      if(isMinBias){
        fHistos->FillTH1("hPtEtaCentNewBinningMB", TMath::Abs(checktrack->Pt()));
        fHistos->FillTH1("hPtEtaCentOldBinningMB", TMath::Abs(checktrack->Pt()));
        for(int icut = 0; icut < 5; icut++){
          if(TMath::Abs(checktrack->Pt()) > static_cast<double>(ptmin[icut])){
            fHistos->FillTH1(Form("hEtaLabDistCutPt%dMB", ptmin[icut]), checktrack->Eta());
            fHistos->FillTH1(Form("hEtaCentDistCutPt%dMB", ptmin[icut]), etacent);
          }
        }
        // check for exclusive classes
        if(!(isEG1 || isEG2 || isEJ1 || isEJ2)){
          fHistos->FillTH1("hPtEtaCentNewBinningMBexcl", TMath::Abs(checktrack->Pt()));
          fHistos->FillTH1("hPtEtaCentOldBinningMBexcl", TMath::Abs(checktrack->Pt()));
          for(int icut = 0; icut < 5; icut++){
            if(TMath::Abs(checktrack->Pt()) > static_cast<double>(ptmin[icut])){
              fHistos->FillTH1(Form("hEtaLabDistCutPt%dMBexcl", ptmin[icut]), checktrack->Eta());
              fHistos->FillTH1(Form("hEtaCentDistCutPt%dMBexcl", ptmin[icut]), etacent);
            }
          }
        }
      }
      if(isEJ1){
        fHistos->FillTH1("hPtEtaCentNewBinningEJ1", TMath::Abs(checktrack->Pt()));
        fHistos->FillTH1("hPtEtaCentOldBinningEJ1", TMath::Abs(checktrack->Pt()));
        for(int icut = 0; icut < 5; icut++){
          if(TMath::Abs(checktrack->Pt()) > static_cast<double>(ptmin[icut])){
            fHistos->FillTH1(Form("hEtaLabDistCutPt%dEJ1", ptmin[icut]), checktrack->Eta());
            fHistos->FillTH1(Form("hEtaCentDistCutPt%dEJ1", ptmin[icut]), etacent);
          }
        }
      }
      if(isEJ2){
        fHistos->FillTH1("hPtEtaCentNewBinningEJ2", TMath::Abs(checktrack->Pt()));
        fHistos->FillTH1("hPtEtaCentOldBinningEJ2", TMath::Abs(checktrack->Pt()));
        for(int icut = 0; icut < 5; icut++){
          if(TMath::Abs(checktrack->Pt()) > static_cast<double>(ptmin[icut])){
            fHistos->FillTH1(Form("hEtaLabDistCutPt%dEJ2", ptmin[icut]), checktrack->Eta());
            fHistos->FillTH1(Form("hEtaCentDistCutPt%dEJ2", ptmin[icut]), etacent);
          }
        }
        // check for exclusive classes
        if(!isEJ1){
          fHistos->FillTH1("hPtEtaCentNewBinningEJ2excl", TMath::Abs(checktrack->Pt()));
          fHistos->FillTH1("hPtEtaCentOldBinningEJ2excl", TMath::Abs(checktrack->Pt()));
          for(int icut = 0; icut < 5; icut++){
            if(TMath::Abs(checktrack->Pt()) > static_cast<double>(ptmin[icut])){
              fHistos->FillTH1(Form("hEtaLabDistCutPt%dEJ2excl", ptmin[icut]), checktrack->Eta());
              fHistos->FillTH1(Form("hEtaCentDistCutPt%dEJ2excl", ptmin[icut]), etacent);
            }
          }
        }
      }
      if(isEG1){
        fHistos->FillTH1("hPtEtaCentNewBinningEG1", TMath::Abs(checktrack->Pt()));
        fHistos->FillTH1("hPtEtaCentOldBinningEG1", TMath::Abs(checktrack->Pt()));
        for(int icut = 0; icut < 5; icut++){
          if(TMath::Abs(checktrack->Pt()) > static_cast<double>(ptmin[icut])){
            fHistos->FillTH1(Form("hEtaLabDistCutPt%dEG1", ptmin[icut]), checktrack->Eta());
            fHistos->FillTH1(Form("hEtaCentDistCutPt%dEG1", ptmin[icut]), checktrack->Eta());
          }
        }
      }
      if(isEG2){
        fHistos->FillTH1("hPtEtaCentNewBinningEG2", TMath::Abs(checktrack->Pt()));
        fHistos->FillTH1("hPtEtaCentOldBinningEG2", TMath::Abs(checktrack->Pt()));
        for(int icut = 0; icut < 5; icut++){
          if(TMath::Abs(checktrack->Pt()) > static_cast<double>(ptmin[icut])){
            fHistos->FillTH1(Form("hEtaLabDistCutPt%dEG2", ptmin[icut]), checktrack->Eta());
            fHistos->FillTH1(Form("hEtaCentDistCutPt%dEG2", ptmin[icut]), checktrack->Eta());
          }
        }
        // check for exclusive classes
        if(!isEG1){
          fHistos->FillTH1("hPtEtaCentNewBinningEG2excl", TMath::Abs(checktrack->Pt()));
          fHistos->FillTH1("hPtEtaCentOldBinningEG2excl", TMath::Abs(checktrack->Pt()));
          for(int icut = 0; icut < 5; icut++){
            if(TMath::Abs(checktrack->Pt()) > static_cast<double>(ptmin[icut])){
              fHistos->FillTH1(Form("hEtaLabDistCutPt%dEG2excl", ptmin[icut]), checktrack->Eta());
              fHistos->FillTH1(Form("hEtaCentDistCutPt%dEG2excl", ptmin[icut]), checktrack->Eta());
            }
          }
        }
      }
    }
  }
  PostData(1, fHistos->GetListOfHistograms());
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

} /* namespace EMCalTriggerPtAnalysis */
