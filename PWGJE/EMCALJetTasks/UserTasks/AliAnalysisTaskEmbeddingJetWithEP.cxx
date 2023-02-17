/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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
#include <fstream>
#include <sstream>
#include <iostream>

#include <TGrid.h>
#include <TChain.h>

#include <TClonesArray.h>
#include <TFile.h>
#include <THashList.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TList.h>
#include <TCanvas.h>
#include <TString.h>

#include <TRandom3.h>
#include <TMath.h>

#include "AliOADBContainer.h"
#include "AliAODVertex.h"
#include "AliAODHandler.h"

#include <AliAnalysisTask.h>
#include <AliAnalysisManager.h>
#include <AliVEventHandler.h>
#include <AliVCluster.h>
#include <AliVParticle.h>
#include <AliLog.h>

#include "AliMultSelection.h"

#include "AliTLorentzVector.h"
#include "AliEmcalJet.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"

#include "AliJEQnVectorHandler.h"

#include "AliLocalRhoParameter.h"
#include "AliRhoParameter.h"

#include "AliAnalysisTaskFlowVectorCorrectionsPWGJE.h"
#include "AliAnalysisTaskEmcalEmbeddingHelper.h"

#include "AliAnalysisTaskEmbeddingJetWithEP.h"

#include "AliDataFile.h"

class AliAnalysisTaskEmbeddingJetWithEP;

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskEmbeddingJetWithEP);
/// \endcond

/**
 * Default constructor. Needed by ROOT I/O
 */
AliAnalysisTaskEmbeddingJetWithEP::AliAnalysisTaskEmbeddingJetWithEP() :
  AliAnalysisTaskEmcalJet(),
  fAOD(nullptr),
  fOutputList(nullptr),
  fEventCuts(),
  fYAMLConfig(),
  fHistManager(),
  fUseRunList(),
  fEmbeddingQA(),
  fDoJetMatchingGeom(kFALSE),
  fDoJetMatchingMCFraction(kFALSE),
  fDoDifferentialRM(kFALSE),
  fMCJetContainer(nullptr),
  fDetJetContainer(nullptr),
  fDetJetContainerPPIntermediate(nullptr),
  fRequireMatchedJetAcc(kFALSE),
  fJetMatchingR(0.),
  fMinSharedPFraction(0.),
  fMCJetMinMatchingPt(0.),
  fDetJetMinMatchingPt(0.),
  fPlotJetMatchCandThresh(0.),
  fMinPt(0.),
  fMaxPt(0.),
  bUseJetCont2Acc(kFALSE),
  fCalibRefFileName(""),
  fCalibRefFile(nullptr),
  fCalibV0Ref(nullptr),
  fCalibRefObjList(nullptr),
  fPileupCut(kFALSE),
  fTPCQnMeasure(kFALSE),
  fPileupCutQA(kFALSE),
  fCalibQA(kFALSE),
  fGainCalibQA(kFALSE),
  fReCentCalibQA(kFALSE),
  fEPQA(kFALSE),
  fTrackQA(kFALSE),
  fBkgQA(kFALSE),
  fCalibType(0),
  fNormMethod(0),
  fV0Combin(kFALSE),
  fQ2VecHandler(nullptr),
  fQ3VecHandler(nullptr),
  fHCorrV0ChWeghts(NULL),
  fHCorrQ2xV0C(NULL),
  fHCorrQ2yV0C(NULL),
  fHCorrQ2xV0A(NULL),
  fHCorrQ2yV0A(NULL),
  fHCorrQ3xV0C(NULL),
  fHCorrQ3yV0C(NULL),
  fHCorrQ3xV0A(NULL),
  fHCorrQ3yV0A(NULL),
  fV0CutPU(NULL),
  fSPDCutPU(NULL),
  fMultCutPU(NULL),
  fCenCutLowPU(NULL),
  fCenCutHighPU(NULL),
  fFitModulationType(kNoFit), fFitModulation(nullptr), hBkgTracks(nullptr),
  fV2ResoV0(0.), fV3ResoV0(0.),
  fQaEventNum(-1)
{
    
    for(Int_t i(0); i < 2; i++){
      q2VecV0M[i] = 0.;
      q2VecV0C[i] = 0.;
      q2VecV0A[i] = 0.;
      q3VecV0M[i] = 0.;
      q3VecV0C[i] = 0.;
      q3VecV0A[i] = 0.;

      q2VecTpcM[i] = 0.;
      q2VecTpcP[i] = 0.;
      q2VecTpcN[i] = 0.;
      q3VecTpcM[i] = 0.;
      q3VecTpcP[i] = 0.;
      q3VecTpcN[i] = 0.;
    }

    for(Int_t i(0); i < 3; i++){
      q2V0[i] = 0.;
      q3V0[i] = 0.;
      q2Tpc[i] = 0.;
      q3Tpc[i] = 0.;

      psi2V0[i] = 0.;
      psi3V0[i] = 0.;
      psi2Tpc[i] = 0.;
      psi3Tpc[i] = 0.;
    }
  
}

/**
 * Standard constructor. Should be used by the user.
 * @param[in] name Name of the task
 */
AliAnalysisTaskEmbeddingJetWithEP::AliAnalysisTaskEmbeddingJetWithEP(const char *name):
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fAOD(nullptr),
  fOutputList(nullptr),
  fEventCuts(),
  fYAMLConfig(),
  fHistManager(name),
  fUseRunList(),
  fEmbeddingQA(),
  fDoJetMatchingGeom(kFALSE),
  fDoJetMatchingMCFraction(kFALSE),
  fDoDifferentialRM(kFALSE),
  fMCJetContainer(nullptr),
  fDetJetContainer(nullptr),
  fDetJetContainerPPIntermediate(nullptr),
  fRequireMatchedJetAcc(kFALSE),
  fJetMatchingR(0.),
  fMinSharedPFraction(0.),
  fMCJetMinMatchingPt(0.),
  fDetJetMinMatchingPt(0.),
  fPlotJetMatchCandThresh(0.),
  fMinPt(0.),
  fMaxPt(0.),
  bUseJetCont2Acc(kFALSE),
  fCalibRefFileName(""),
  fCalibRefFile(nullptr),
  fCalibRefObjList(nullptr),
  fCalibV0Ref(nullptr),
  fPileupCut(kFALSE),
  fTPCQnMeasure(kFALSE),
  fPileupCutQA(kFALSE),
  fCalibQA(kFALSE),
  fGainCalibQA(kFALSE),
  fReCentCalibQA(kFALSE),
  fEPQA(kFALSE),
  fTrackQA(kFALSE),
  fBkgQA(kFALSE),
  fCalibType(0),
  fNormMethod(0),
  fV0Combin(kFALSE),
  fQ2VecHandler(nullptr),
  fQ3VecHandler(nullptr),
  fHCorrV0ChWeghts(NULL),
  fHCorrQ2xV0C(NULL),
  fHCorrQ2yV0C(NULL),
  fHCorrQ2xV0A(NULL),
  fHCorrQ2yV0A(NULL),
  fHCorrQ3xV0C(NULL),
  fHCorrQ3yV0C(NULL),
  fHCorrQ3xV0A(NULL),
  fHCorrQ3yV0A(NULL),
  fV0CutPU(NULL),
  fSPDCutPU(NULL),
  fMultCutPU(NULL),
  fCenCutLowPU(NULL),
  fCenCutHighPU(NULL),
  fFitModulationType(kNoFit), fFitModulation(nullptr), hBkgTracks(nullptr),
  fV2ResoV0(0.), fV3ResoV0(0.),
  fQaEventNum(-1)
{
  
    for(Int_t i(0); i < 2; i++){
      q2VecV0M[i] = 0.;
      q2VecV0C[i] = 0.;
      q2VecV0A[i] = 0.;
      q3VecV0M[i] = 0.;
      q3VecV0C[i] = 0.;
      q3VecV0A[i] = 0.;

      q2VecTpcM[i] = 0.;
      q2VecTpcP[i] = 0.;
      q2VecTpcN[i] = 0.;
      q3VecTpcM[i] = 0.;
      q3VecTpcP[i] = 0.;
      q3VecTpcN[i] = 0.;
    }

    for(Int_t i(0); i < 3; i++){
      q2V0[i] = 0.;
      q3V0[i] = 0.;
      q2Tpc[i] = 0.;
      q3Tpc[i] = 0.;

      psi2V0[i] = 0.;
      psi3V0[i] = 0.;
      psi2Tpc[i] = 0.;
      psi3Tpc[i] = 0.;
    }
  
  
  
  if(fLocalRhoName=="") fLocalRhoName = Form("LocalRhoFrom_%s", GetName());
  SetMakeGeneralHistograms(kTRUE);
}

/**
 * Destructor
 */
AliAnalysisTaskEmbeddingJetWithEP::~AliAnalysisTaskEmbeddingJetWithEP()
{
  
  if(fCalibRefFile) {fCalibRefFile->Close(); fCalibRefFile = 0x0;}
  if(fCalibV0Ref)   {delete fCalibV0Ref;   fCalibV0Ref = 0x0;}
  if(fQ2VecHandler) {delete fQ2VecHandler; fQ2VecHandler = 0x0;}
  if(fQ3VecHandler) {delete fQ3VecHandler; fQ3VecHandler = 0x0;}
  if(fFitModulation) {delete fFitModulation; fFitModulation = 0x0;}
  if(hBkgTracks) {delete hBkgTracks; hBkgTracks = 0x0;}
}

void AliAnalysisTaskEmbeddingJetWithEP::SetRunList(bool removeDummyTask)
{
  
  fYAMLConfig.AddConfiguration(fRunListFileName, "runlist");
  fYAMLConfig.Initialize();
  fYAMLConfig.GetProperty("runlist", fUseRunList);
}


/**
 * Performing run-independent initialization.
 * Here the histograms should be instantiated.
 */
void AliAnalysisTaskEmbeddingJetWithEP::UserCreateOutputObjects()
{
  
  fOutputList = new TList();
  fOutputList->SetOwner(true);

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  fEventCuts.AddQAplotsToList(fOutput);
  fEventCuts.OverrideAutomaticTriggerSelection(\
    AliVEvent::kAnyINT | AliVEvent::kCentral | AliVEvent::kSemiCentral);
  
  for (Int_t i=0; i<3; i++) {
    auto jetCont = GetJetContainer(i);
    TString jetContName = jetCont->GetName();
    if (jetContName.Contains("mcparticles"))   fMCJetContainer = jetCont;
    else if (jetContName.Contains("Combined")) fDetJetContainer = jetCont;
    else  fDetJetContainerPPIntermediate = jetCont;
    
  }

  if (!fMCJetContainer) Printf("No MC jet container found!");
  Printf("mcJetContainer: %s", fMCJetContainer->GetName());

  if (!fDetJetContainer) Printf("No det-level jet container found!");
  Printf("det-level JetContainer: %s", fDetJetContainer->GetName());
  
  if (!fDetJetContainerPPIntermediate) {
    Printf("No intermediate pp det-level jet container found, despite MC-fraction matching enabled!");
  }
  Printf("Intermediate pp det-level JetContainer: %s", fDetJetContainerPPIntermediate->GetName());
  

  // == s == Set Out put Hist grams  ###########################################
  AllocateJetHistograms();
  AllocateMatchedJetHistograms();
  // == e == Set Out put Hist grams  ###########################################
  

  // == s == Set embedding Helper   ############################################
  const AliAnalysisTaskEmcalEmbeddingHelper * embeddingHelper = AliAnalysisTaskEmcalEmbeddingHelper::GetInstance();
  if (embeddingHelper) {
    bool res = fEmbeddingQA.Initialize();
    if(res) fEmbeddingQA.AddQAPlotsToList(fOutput);
  }
  // == e == Set embedding Helper   ############################################

  // == s == Add Objects into output file  #####################################
  TIter next(fHistManager.GetListOfHistograms());
  TObject* obj = 0;
  while ((obj = next())) {
    fOutput->Add(obj);
  }
  // == e == Add Objects into output file  #####################################
  

  // == s == Calib root file include  ==============================--------===-
  TString tempCalibFileName = AliDataFile::GetFileName(fCalibRefFileName.Data());
  TString tempCalibLocalFileName;
  // std::cout << fCalibRefFileName << std::endl;
  tempCalibLocalFileName = fCalibRefFileName;
  // Check access to CVMFS (will only be displayed locally)
  if(fCalibRefFileName.BeginsWith("alien://") && !gGrid){
    AliInfo("Trying to connect to AliEn ...");
    TGrid::Connect("alien://");
  }
  if(!tempCalibFileName.IsNull()) fCalibRefFile = TFile::Open(tempCalibFileName.Data());
  if(tempCalibFileName.IsNull())  fCalibRefFile = TFile::Open(tempCalibLocalFileName.Data());
  if(!fCalibRefFile) {
    AliWarning("V0-TPC Gain calibration file cannot be opened\n");
    return;
  }
  // == e == Calib root file include  ==============================--------===-


  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}

/*
 * This function allocates the histograms for basic jet QA.
 * A set of histograms (pT, eta, phi, area, number of jets, corrected pT) is allocated
 * per each jet container and per each centrality bin.
 */
void AliAnalysisTaskEmbeddingJetWithEP::AllocateJetHistograms()
{
  TString histName;
  TString histtitle;
  TString groupName;
  AliJetContainer* jetCont = 0;
  TIter next(&fJetCollArray);
  while ((jetCont = static_cast<AliJetContainer*>(next()))) {
    groupName = jetCont->GetName();
    // Protect against creating the histograms twice
    if (fHistManager.FindObject(groupName)) {
      AliWarning(TString::Format("%s: Found groupName %s in hist manager. The jet containers will be filled into the same histograms.", GetName(), groupName.Data()));
      continue;
    }
    fHistManager.CreateHistoGroup(groupName);

    TString rhoGroupName = TString::Format("%s/Rho", groupName.Data());
    if (fHistManager.FindObject(rhoGroupName)) {
      AliWarning(TString::Format("%s: Found groupName %s in hist manager. The jet containers will be filled into the same histograms.", GetName(), rhoGroupName.Data()));
      continue;
    }
    fHistManager.CreateHistoGroup(rhoGroupName);


    TString InclusiveGroupName = TString::Format("%s/Inclusive", groupName.Data());
    if (fHistManager.FindObject(InclusiveGroupName)) {
      AliWarning(TString::Format("%s: Found groupName %s in hist manager. The jet containers will be filled into the same histograms.", GetName(), InclusiveGroupName.Data()));
      continue;
    }
    fHistManager.CreateHistoGroup(InclusiveGroupName);

    TString IPlaneGroupName = TString::Format("%s/InPlane", groupName.Data());
    if (fHistManager.FindObject(IPlaneGroupName)) {
      AliWarning(TString::Format("%s: Found groupName %s in hist manager. The jet containers will be filled into the same histograms.", GetName(), IPlaneGroupName.Data()));
      continue;
    }
    fHistManager.CreateHistoGroup(IPlaneGroupName);

    TString OPlaneGroupName = TString::Format("%s/OutOfPlane", groupName.Data());
    if (fHistManager.FindObject(OPlaneGroupName)) {
      AliWarning(TString::Format("%s: Found groupName %s in hist manager. The jet containers will be filled into the same histograms.", GetName(), OPlaneGroupName.Data()));
      continue;
    }
    fHistManager.CreateHistoGroup(OPlaneGroupName);

    // A vs. pT    
    histName = TString::Format("%s/CentVsPtVsPtdence", rhoGroupName.Data());
    histtitle = histName + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});#it{A}_{jet}/#pi#it{R}{2}";

    for (Int_t cent = 0; cent < fNcentBins; cent++) {

      histName = TString::Format("%s/hNJets_%d", groupName.Data(), cent);
      histtitle = TString::Format("%s;number of jets;events", histName.Data());
      fHistManager.CreateTH1(histName, histtitle, 500, 0, 500);

      // histograms for jet angle relative to the event plane
      histName = TString::Format("%s/hJetPhiMinusPsi2_%d", groupName.Data(), cent);
      histtitle = "Jet phi minus psi2";
      fHistManager.CreateTH1(histName, histtitle, fNbins / 2, 0, TMath::TwoPi());
      histName = TString::Format("%s/hJetPhiMinusPsi3_%d", groupName.Data(), cent);
      histtitle = "Jet phi minus psi3";
      fHistManager.CreateTH1(histName, histtitle, fNbins / 2, 0, TMath::TwoPi());


      histName = TString::Format("%s/hDeltaPt_%d", groupName.Data(), cent);
      histtitle = "delta pt without EP distribution";
      fHistManager.CreateTH1(histName, histtitle, 300, -50, 250);
      histName = TString::Format("%s/hDeltaPt_Local_%d", groupName.Data(), cent);
      histtitle = "delta pt according EP distribution";
      fHistManager.CreateTH1(histName, histtitle, 300, -50, 250);
      histName = TString::Format("%s/hPhiVsDeltaPt_Local_%d", groupName.Data(), cent);
      histtitle = "phi vs delta pt according EP distribution";
      fHistManager.CreateTH2(histName, histtitle, fNbins/2, 0, TMath::TwoPi(), 300, -50, 250);
      
      // Rho histograms
      // histName = TString::Format("%s/hRhoVsCent", rhoGroupName.Data());
      // histtitle = histName + ";Centrality (%);#rho (GeV/#it{c});counts";
      // fHistManager.CreateTH2(histName, histtitle.Data(), 50, 0, 100, 100, 0, 500);
      histName = TString::Format("%s/hJetRho_%d", rhoGroupName.Data(), cent);
      histtitle = "Rho";
      fHistManager.CreateTH1(histName, histtitle, fNbins, 0.0, 300.0);
      histName = TString::Format("%s/hJetRhoLocal_%d", rhoGroupName.Data(), cent);
      histtitle = "Rho Local";
      fHistManager.CreateTH1(histName, histtitle, fNbins, 0.0, 300.0);

      histName = TString::Format("%s/hJetLocalRhoVsAverageRho_%d", rhoGroupName.Data(), cent);
      histtitle = "Local rho versus average rho";
      fHistManager.CreateTH2(histName, histtitle, fNbins, 0.0, 300.0, fNbins, 0.0, 300.0);
      histName = TString::Format("%s/hJetCorrPtLocalVsJetCorrPt_%d", rhoGroupName.Data(), cent);
      histtitle = "Local rho adjusted jet pT versus average rho adjusted jet pT";
      fHistManager.CreateTH2(histName, histtitle, fNbins, 0.0, 100.0, fNbins, 0.0, 100.0);

      // histo local rho vs delta phi
      histName = TString::Format("%s/hJetRhoVsDeltaPhi_%d", rhoGroupName.Data(), cent);
      histtitle = "Local rho versus angle relative to event plane";
      fHistManager.CreateTH2(histName, histtitle, fNbins, 0.0, TMath::TwoPi(), fNbins, 0.0, 300.0);

      //inlucive
      histName = TString::Format("%s/hJetPt_%d", InclusiveGroupName.Data(), cent);
      histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c});counts", histName.Data());
      fHistManager.CreateTH1(histName, histtitle, fNbins, fMinBinPt, fMaxBinPt);
      histName = TString::Format("%s/hJetCorrPt_%d", InclusiveGroupName.Data(), cent);
      histtitle = TString::Format("%s;#it{p}_{T,jet}^{corr} (GeV/#it{c});counts", histName.Data());
      fHistManager.CreateTH1(histName, histtitle, fNbins, -fMaxBinPt / 2, fMaxBinPt / 2);
      histName = TString::Format("%s/hJetCorrPtLocal_%d", InclusiveGroupName.Data(), cent);
      histtitle = TString::Format("%s;#it{p}_{T,jet}^{corr} local (GeV/#it{c});counts", histName.Data());
      fHistManager.CreateTH1(histName, histtitle, fNbins, -fMaxBinPt / 2, fMaxBinPt / 2);
      
      histName = TString::Format("%s/hJetArea_%d", InclusiveGroupName.Data(), cent);
      histtitle = TString::Format("%s;#it{A}_{jet};counts", histName.Data());
      fHistManager.CreateTH1(histName, histtitle, fNbins / 2, 0, 3);

      histName = TString::Format("%s/hJetPhi_%d", InclusiveGroupName.Data(), cent);
      histtitle = TString::Format("%s;#it{#phi}_{jet};counts", histName.Data());
      fHistManager.CreateTH1(histName, histtitle, fNbins / 2, 0, TMath::TwoPi());

      histName = TString::Format("%s/hJetEta_%d", InclusiveGroupName.Data(), cent);
      histtitle = TString::Format("%s;#it{#eta}_{jet};counts", histName.Data());
      fHistManager.CreateTH1(histName, histtitle, fNbins / 6, -1, 1);
        

      //v2 in plane
      histName = TString::Format("%s/hJetPt_%d", IPlaneGroupName.Data(), cent);
      histtitle = "Jet yeild of in-plane (v2)";
      fHistManager.CreateTH1(histName, histtitle, 300, -50, 250);
      histName = TString::Format("%s/hCorrJetPt_%d", IPlaneGroupName.Data(), cent);
      histtitle = "corr Jet yeild of in-plane (v2)";
      fHistManager.CreateTH1(histName, histtitle, 300, -50, 250);

      histName = TString::Format("%s/hJetArea_%d", IPlaneGroupName.Data(), cent);
      histtitle = TString::Format("%s;#it{A}_{jet};counts", histName.Data());
      fHistManager.CreateTH1(histName, histtitle, fNbins / 2, 0, 3);
      histName = TString::Format("%s/hJetPhi_%d", IPlaneGroupName.Data(), cent);
      histtitle = TString::Format("%s;#it{#phi}_{jet};counts", histName.Data());
      fHistManager.CreateTH1(histName, histtitle, fNbins / 2, 0, TMath::TwoPi());
      histName = TString::Format("%s/hJetEta_%d", IPlaneGroupName.Data(), cent);
      histtitle = TString::Format("%s;#it{#eta}_{jet};counts", histName.Data());
      fHistManager.CreateTH1(histName, histtitle, fNbins / 6, -1, 1);

      //v2 out of plane
      histName = TString::Format("%s/hJetPt_%d", OPlaneGroupName.Data(), cent);
      histtitle = "Jet yeild of out-plane (v2)";
      fHistManager.CreateTH1(histName, histtitle, 300, -50, 250);
      histName = TString::Format("%s/hCorrhJetPt_%d", OPlaneGroupName.Data(), cent);
      histtitle = "corr Jet yeild of in-plane (v2)";
      fHistManager.CreateTH1(histName, histtitle, 300, -50, 250);

      histName = TString::Format("%s/hJetArea_%d", OPlaneGroupName.Data(), cent);
      histtitle = TString::Format("%s;#it{A}_{jet};counts", histName.Data());
      fHistManager.CreateTH1(histName, histtitle, fNbins / 2, 0, 3);
      histName = TString::Format("%s/hJetPhi_%d", OPlaneGroupName.Data(), cent);
      histtitle = TString::Format("%s;#it{#phi}_{jet};counts", histName.Data());
      fHistManager.CreateTH1(histName, histtitle, fNbins / 2, 0, TMath::TwoPi());
      histName = TString::Format("%s/hJetEta_%d", OPlaneGroupName.Data(), cent);
      histtitle = TString::Format("%s;#it{#eta}_{jet};counts", histName.Data());
      fHistManager.CreateTH1(histName, histtitle, fNbins / 6, -1, 1);
      
    }
  }
}


/*
 * This function allocates histograms for matched truth-det jets in the case of embedding.
 */
void AliAnalysisTaskEmbeddingJetWithEP::AllocateMatchedJetHistograms()
{
  TString histname;
  Double_t jetR=0;
  auto jetCont1 = GetJetContainer(0);
  jetR = jetCont1->GetJetRadius();

  fHistManager.CreateHistoGroup("MatchedJetHistograms");

  TString title;
  TString groupName = "MatchedJetHistograms";

  Int_t nPtBins1 = TMath::CeilNint(fMaxPt-fMinPt);
  Int_t nPtBinsTruth2 = TMath::CeilNint(fMaxPt/2);
    
  // Response matrix, (centrality, pT-truth, pT-det)
  Int_t nbinsx = 9; Int_t minx = 0; Int_t maxx = 90;
  Int_t nbinsy = fMaxPt; Int_t miny = 0; Int_t maxy = fMaxPt;
  Int_t nbinsz = nPtBins1; Int_t minz = fMinPt; Double_t maxz = fMaxPt;

  if (fDoJetMatchingGeom) {
    // Matching distance, (pT-det, pT-truth, deltaR)
    nbinsx = nPtBins1; minx = fMinPt; maxx = fMaxPt;
    nbinsy = nPtBinsTruth2; miny = 0; maxy = fMaxPt;
    nbinsz = 15; minz = 0; maxz = 1.5*jetR;
    histname = TString::Format("%s/hMatchingDistance",groupName.Data());
    title = histname + ";#it{p}_{T}^{det} (GeV/#it{c});#it{p}_{T}^{truth} (GeV/#it{c});R";
    fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, fMaxPt, miny, maxy, nbinsz, minz, maxz);

    // Jet matching QA (copied from AliAnalysisTaskEmcalJetHCorrelations.cxx)
    histname = TString::Format("%s/fHistJetMatchingQA",groupName.Data());
    title = histname;
    std::vector<std::string> binLabels = {"noMatch", "matchedJet", "uniqueMatch","jetDistance", "passedAllCuts"};
    auto histMatchedJetCuts = fHistManager.CreateTH1(histname.Data(), title.Data(),binLabels.size(), 0, binLabels.size());
    // Set label names
    for (unsigned int i = 1; i <= binLabels.size(); i++) {
      histMatchedJetCuts->GetXaxis()->SetBinLabel(i, binLabels.at(i-1).c_str());
    }
    histMatchedJetCuts->GetYaxis()->SetTitle("Number of jets");
  }
  if (fDoJetMatchingMCFraction) {
    // (det pT, shared MC fraction, deltaR) of closest jets
    nbinsx = nPtBins1; minx = fMinPt; maxx = fMaxPt;
    nbinsy = 20; miny = 0; maxy = 1.;
    nbinsz = 15; minz = 0; maxz = 1.5*jetR;
    histname = TString::Format("%s/hMatchingDistanceVsMCFraction",groupName.Data());
    title = histname + ";#it{p}_{T}^{det} (GeV/#it{c});MC fraction;#DeltaR";
    fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);

    // Jet matching QA (copied from AliAnalysisTaskEmcalJetHCorrelations.cxx)
    histname = TString::Format("%s/fHistJetMatchingQA",groupName.Data());
    title = histname;
    std::vector<std::string> binLabels = {"noMatch", "matchedJet", "sharedMomentumFraction", "partLevelMatchedJet", "jetDistancePPdet", "jetDistancePPtruth", "passedAllCuts"};
    auto histMatchedJetCuts = fHistManager.CreateTH1(histname.Data(), title.Data(), binLabels.size(), 0, binLabels.size());
    // Set label names
    for (unsigned int i = 1; i <= binLabels.size(); i++) {
      histMatchedJetCuts->GetXaxis()->SetBinLabel(i, binLabels.at(i-1).c_str());
    }
    histMatchedJetCuts->GetYaxis()->SetTitle("Number of jets");
  }

  std::vector<TString> fHistJetEPName = {"Inclusive", "InPlane", "OutOfPlane"}; ///<
  for(Int_t jetEPBin = 0; jetEPBin < 3; jetEPBin++){
    TString subGroupName = fHistJetEPName.at(jetEPBin);
    groupName = TString::Format("MatchedJetHistograms/%s",subGroupName.Data());
    fHistManager.CreateHistoGroup(groupName);
    
    //This is a 5-dim RM with information on the angularity and matching distance
    if (fDoDifferentialRM) {
      //setup the THnSparse
      Int_t nCentBins=20;

      TString titleThn[6]= {"#it{p}_{T}^{truth} (GeV/#it{c})", "#it{p}_{T,corr}^{det} (GeV/#it{c})", "#Delta#it{R}", "shared mom fraction" ,"angularity", "Centrality (%)"};
      Int_t nbinsThn[6]  = {(Int_t)fMaxPt, (Int_t)fMaxPt, 15, 100, 100, nCentBins};
      Double_t minThn[6] = {0., 0., 0., 0.,0., 0.};
      Double_t maxThn[6] = {fMaxPt, fMaxPt, 1.5*jetR, 1, jetR, 100};
      histname = TString::Format("%s/hResponseMatrixDiff",groupName.Data());
      // (1) pt part LvL, (2) pt det LvL, (3) Matching distance (4) shared momentum fraction (5) angularity (6) centrality
      THnSparse* thn = fHistManager.CreateTHnSparse(histname.Data(), histname.Data(), 6, nbinsThn, minThn, maxThn);
      for (Int_t i = 0; i < 6; i++) {
        thn->GetAxis(i)->SetTitle(titleThn[i]);
        //thn->SetBinEdges(i, binEdges[i]);
      }
    }
    //This is a 3D RM for PbPb and a 2D RM for pp
    else {
      histname = TString::Format("%s/hResponseMatrix",groupName.Data());
      title = histname + ";Centrality (%);#it{p}_{T}^{truth} (GeV/#it{c});#it{p}_{T,corr}^{det} (GeV/#it{c})";
      fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
    }

    
    // JES shift, (centrality, pT-truth, (pT-det - pT-truth) / pT-truth)
    nbinsx = 20; minx = 0; maxx = 100;
    nbinsy = nPtBinsTruth2; miny = 0; maxy = fMaxPt;
    nbinsz = 250; minz = -5.; maxz = 5.;
    histname = TString::Format("%s/hJESshift",groupName.Data());
    title = histname + ";Centrality (%);#it{p}_{T}^{truth} (GeV/#it{c});#frac{#it{p}_{T,corr}^{det} - #it{p}_{T}^{truth}}{#it{p}_{T}^{truth}}";
    fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
    

    // NEF of det-level matched jets, (centrality, pT-truth, NEF)
    nbinsx = 20; minx = 0; maxx = 100;
    nbinsy = fMaxPt; miny = 0; maxy = fMaxPt;
    nbinsz = 50; minz = 0; maxz = 1.;
    histname = TString::Format("%s/hNEFVsPt",groupName.Data());
    title = histname + ";Centrality (%);#it{p}_{T}^{truth} (GeV/#it{c});Calo energy fraction";
    fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
    
    // z-leading (charged) of det-level matched jets, (centrality, pT-truth, z-leading)
    nbinsx = 20; minx = 0; maxx = 100;
    nbinsy = nPtBinsTruth2; miny = 0; maxy = fMaxPt;
    nbinsz = 50; minz = 0; maxz = 1.;
    histname = TString::Format("%s/hZLeadingVsPt",groupName.Data());
    title = histname + ";Centrality (%);#it{p}_{T}^{truth} (GeV/#it{c});#it{z}_{leading}";
    fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
  }
}


/**
 * This function is executed automatically for the first event.
 * Some extra initialization can be performed here.
 */
void AliAnalysisTaskEmbeddingJetWithEP::ExecOnce()
{
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  
  if(!fLocalRho) {
    fLocalRho = new AliLocalRhoParameter(fLocalRhoName.Data(), 0); 
      if(!(InputEvent()->FindListObject(fLocalRho->GetName()))) {
        InputEvent()->AddObject(fLocalRho);
      } else {
        AliFatal(Form("%s: Container with name %s already present. Aborting", \
        GetName(), fLocalRho->GetName()));
      }
  }
  
  AliAnalysisTaskEmcalJet::ExecOnce();
  if(!GetJetContainer()) AliFatal(Form("%s: Couldn't find jet container. Aborting !", GetName()));


}


/**
 * Run analysis code here, if needed.
 * It will be executed before FillHistograms().
 * If this function return kFALSE, FillHistograms() will *not*
 * be executed for the current event
 * @return Always kTRUE
 */
Bool_t AliAnalysisTaskEmbeddingJetWithEP::Run()
{
  // std::cout << "ChecKuma Run Number === " << CheckRunNum << "==================" << std::endl;
  if(!fEventCuts.AcceptEvent(InputEvent())) return kFALSE;
  
  if(fPileupCut){
    SetupPileUpRemovalFunctions();
    Bool_t kPileupCutEvent = CheckEventIsPileUp2018();
    if(kPileupCutEvent) return kFALSE;
  }

  DoEventPlane();
  SetModulationRhoFit();
  // std::cout << "Fomula = " << fFitModulation->GetExpFormula() << std::endl;
  MeasureBkg();
  DoJetLoop();

  FillMatchedJetHistograms();

  // Only fill the embedding qa plots if:
  //  - We are using the embedding helper
  //  - The class has been initialized
  //  - Both jet collections are available
  if (fEmbeddingQA.IsInitialized()) {
    fEmbeddingQA.RecordEmbeddedEventProperties();
  }
  return kTRUE;
}

void AliAnalysisTaskEmbeddingJetWithEP::DoEventPlane(){

  if (!fAOD && AODEvent() && IsStandardAOD()) {
      // In case there is an AOD handler writing a standard AOD, use the AOD
      // event in memory rather than the input (ESD) event.
      fAOD = dynamic_cast<AliAODEvent*> (AODEvent());
  }
  if (!fAOD) {
      AliWarning("AliAnalysisTaskJetQnVectors::Exec(): bad AOD");
      return;
  }
  AliAODHandler* aodHandler = static_cast<AliAODHandler*>((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
  if(!aodHandler) {
    AliWarning("AliAnalysisTaskJetQnVectors::Exec(): No AliInputEventHandler!");
    return;
  }
  
  //== s == qn Calibration  111111111111111111111111111111111111111111111111111
  // std::cout << "bef calib (qx,qy) = " << q2VecV0M[0] << "," << q2VecV0M[1] << ")" << std::endl;
  QnGainCalibration();
  if(0 ){
    std::cout << "recent calibA (qx,qy) = " \
      << q2VecV0A[0] << "," << q2VecV0A[1] << ")" << std::endl;
    std::cout << "recent calibM (qx,qy) = " \
      << q2VecV0M[0] << "," << q2VecV0M[1] << ")" << std::endl;
  }
  // std::cout << "gain calib (qx,qy) = " << q2VecV0M[0] << "," << q2VecV0M[1] << ")" << std::endl;
  QnRecenteringCalibration();

  //== s == combin V0C and V0A  ################################################
  Double_t q2ChiV0C = 0.;
  Double_t q2ChiV0A = 0.;
  Double_t q3ChiV0C = 0.;
  Double_t q3ChiV0A = 0.;
  if(fV0Combin){
    Double_t psiReso = 0;
    q2ChiV0A = CalculateEventPlaneChi(psiReso);
    q2ChiV0C = CalculateEventPlaneChi(psiReso);
    q3ChiV0A = CalculateEventPlaneChi(psiReso);
    q3ChiV0C = CalculateEventPlaneChi(psiReso);
  }else {
    q2ChiV0A = 1.0;
    q3ChiV0A = 1.0;
  }
  q2VecV0M[0] = q2ChiV0C*q2ChiV0C*q2VecV0C[0] + q2ChiV0A*q2ChiV0A*q2VecV0A[0];
  q2VecV0M[1] = q2ChiV0C*q2ChiV0C*q2VecV0C[1] + q2ChiV0A*q2ChiV0A*q2VecV0A[1];
  q3VecV0M[0] = q3ChiV0C*q2ChiV0C*q3VecV0C[0] + q3ChiV0A*q2ChiV0A*q3VecV0A[0];
  q3VecV0M[1] = q3ChiV0C*q2ChiV0C*q3VecV0C[1] + q3ChiV0A*q2ChiV0A*q3VecV0A[1];
  if(0){
    std::cout << "recent calibA (qx,qy) = " \
      << q2VecV0A[0] << "," << q2VecV0A[1] << ")" << std::endl;
    std::cout << "recent calibM (qx,qy) = " \
      << q2VecV0M[0] << "," << q2VecV0M[1] << ")" << std::endl;
  }
  //== s == combin V0C and V0A  ################################################

  psi2V0[0] = CalcEPAngle(q2VecV0M[0], q2VecV0M[1]);
  psi2V0[1] = CalcEPAngle(q2VecV0C[0], q2VecV0C[1]);
  psi2V0[2] = CalcEPAngle(q2VecV0A[0], q2VecV0A[1]);

  psi3V0[0] = CalcEPAngle(q3VecV0M[0], q3VecV0M[1]);
  psi3V0[1] = CalcEPAngle(q3VecV0C[0], q3VecV0C[1]);
  psi3V0[2] = CalcEPAngle(q3VecV0A[0], q3VecV0A[1]);
  //== e == qn Calibration  111111111111111111111111111111111111111111111111111
}


void AliAnalysisTaskEmbeddingJetWithEP::SetModulationRhoFit() 
{
  // set modulation fit
  TString histName;

  if(fFitModulation) delete fFitModulation;
  fFitModulation = 0x0;
  const char * fitFunction = "[0]*(1.+2.*([1]*TMath::Cos(2.*(x-[2]))+[3]*TMath::Cos(3.*(x-[4]))))";
  switch (fFitModulationType)  {
    case kNoFit : { fFitModulation = new TF1("fix_kNoFit", "[0]", 0, TMath::TwoPi()); } break;
    case kV2 : {
      fitFunction = "[0]*([1]+[2]*[3]*TMath::Cos([2]*(x-[4])))";
      fFitModulation = new TF1("fit_kV2", fitFunction, 0, TMath::TwoPi());
      fFitModulation->SetParameter(0, 0.);  // normalization
      fFitModulation->SetParameter(1, 0.2); // v2
    } break;
    case kCombined: {
      fitFunction = "[0]*(1.+2.*([1]*TMath::Cos(2.*(x-[2]))+[3]*TMath::Cos(3.*(x-[4]))))";
      fFitModulation = new TF1("fit_kCombined", fitFunction, 0, TMath::TwoPi());
      fFitModulation->SetParameter(0, 0.);       // normalization
      fFitModulation->SetParameter(1, 0.2);      // v2
      fFitModulation->SetParameter(3, 0.2);      // v3
    } break;
    default : { // for the combined fit, the 'direct fourier series' or the user supplied vn values we use v2 and v3
      fitFunction = "[0]*(1.+2.*([1]*TMath::Cos(2.*(x-[2]))+[3]*TMath::Cos(3.*(x-[4]))))";
      fFitModulation = new TF1("fit_kCombined", fitFunction, 0, TMath::TwoPi());
      fFitModulation->SetParameter(0, 0.);       // normalization
      fFitModulation->SetParameter(1, 0.2);      // v2
      fFitModulation->SetParameter(3, 0.2);      // v3
    } break;
  }

  

  if(hBkgTracks) delete hBkgTracks;
  histName = "hBkgTracks";
  // hBkgTracks = new TH1F(histName, histName, 100, 0.0, TMath::TwoPi());
  hBkgTracks = new TH1F(histName, histName, 25, 0.0, TMath::TwoPi());
}


void AliAnalysisTaskEmbeddingJetWithEP::MeasureBkg(){
  TString histName;
  
  UInt_t sumAcceptedTracks = 0;
  AliParticleContainer* partCont = 0;
  TIter next(&fParticleCollArray);
  while ((partCont = static_cast<AliParticleContainer*>(next()))) {
    // groupname = partCont->GetName();
    for(auto part : partCont->accepted()) {
      if (!part) continue;
      if (partCont->GetLoadedClass()->InheritsFrom("AliVTrack")) {
        const AliVTrack* track = static_cast<const AliVTrack*>(part);

        Float_t trackDeltaPhi = track->Phi() - psi2V0[0];
        if (trackDeltaPhi < 0.0) trackDeltaPhi += TMath::TwoPi();
        hBkgTracks->Fill(trackDeltaPhi);
      }
    }
  }
  
  fLocalRho->SetVal(fRho->GetVal());
  fFitModulation->SetParameter(0, fLocalRho->GetVal());
  fFitModulation->FixParameter(2, psi2V0[0]);
  fFitModulation->FixParameter(4, psi3V0[0]);
  
  hBkgTracks->Fit(fFitModulation, "N0Q"); 
  
  if(0){
    TCanvas *cBkgRhoFit = new TCanvas("cBkgRhoFit", "cBkgRhoFit", 2000, 1500);
    
    TH1F* hBkgTracks_Event = (TH1F*) hBkgTracks->Clone("hnew");
    histName = hBkgTracks->GetName();
    // histName = hBkgTracks->GetName() + std::to_string(CheckRunNum);
    hBkgTracks_Event->SetName(histName);
    hBkgTracks_Event->Draw("E");
    hBkgTracks->Fit(fFitModulation, "N0Q");
    fFitModulation->SetLineColor(632);
    fFitModulation->Draw("same");

    TF1* rhoFitV2Com = new TF1("rhoFitV2Com", "[0]*(1.+2.*([1]*TMath::Cos(2.*(x-[2]))))", 0.0, TMath::TwoPi());
    rhoFitV2Com->SetParameter(0, fFitModulation->GetParameter(0));
    rhoFitV2Com->SetParameter(1, fFitModulation->GetParameter(1));//v2
    rhoFitV2Com->SetParameter(2, fFitModulation->GetParameter(2));//psi2
    rhoFitV2Com->SetLineColor(808);
    rhoFitV2Com->Draw("same");

    TF1* rhoFitV3Com = new TF1("rhoFitV3Com", "[0]*(1.+2.*([1]*TMath::Cos(3.*(x-[2]))))", 0.0, TMath::TwoPi());
    rhoFitV3Com->SetParameter(0, fFitModulation->GetParameter(0));
    rhoFitV3Com->SetParameter(1, fFitModulation->GetParameter(1));//v3
    rhoFitV3Com->SetParameter(2, fFitModulation->GetParameter(2));//psi3
    rhoFitV3Com->SetLineColor(824);
    rhoFitV3Com->Draw("same");

    histName = "checkOutput/cBkgRhoFit_Cent" + std::to_string(fCentBin) +".root";
    cBkgRhoFit->SaveAs(histName);
    
    // histName = "checkOutput/hBkgTracks_Event" + std::to_string(CheckRunNum) +".root";
    hBkgTracks_Event->SaveAs(histName);
    delete cBkgRhoFit;
    delete hBkgTracks_Event;
  }
  
  // fV2ResoV0 = CalcEPReso(2, psi2V0[0], psi2Tpc[1], psi2Tpc[2]);
  // fV3ResoV0 = CalcEPReso(3, psi3V0[0], psi3Tpc[1], psi3Tpc[2]);
  // std::cout << "v2Reso = " << fV2ResoV0 << ", v3Reso = " << fV3ResoV0 << std::endl;

  fLocalRho->SetLocalRho(fFitModulation);
  // fLocalRho->SetVal(fRho->GetVal());

}


Double_t AliAnalysisTaskEmbeddingJetWithEP::CalcEPReso(Int_t n, \
  Double_t &psiA, Double_t &psiB, Double_t &psiC){
  
  Double_t vnReso = -999.;
  vnReso = TMath::Sqrt((TMath::Abs(TMath::Cos(n*(psiA - psiB))) \
                          * TMath::Abs(TMath::Cos(n*(psiA - psiC)))) \
                        / TMath::Abs(TMath::Cos(n*(psiB - psiC))));

  return vnReso;
}



/**
 * This function performs a loop over the reconstructed jets
 * in the current event and fills the relevant histograms.
 */
void AliAnalysisTaskEmbeddingJetWithEP::DoJetLoop()
{
  TString histName;
  TString groupName;
  AliJetContainer* jetCont = 0;
  TIter next(&fJetCollArray);
  
  while ((jetCont = static_cast<AliJetContainer*>(next()))) {
    
    groupName = jetCont->GetName();
    UInt_t count = 0;
    
    Double_t jetR = jetCont->GetJetRadius();

    Double_t rhoVal = 0;
    Double_t leadingJetEta = -999.;
    Double_t leadingJetPhi = -999.;
    Double_t leadingJetPt  = -999.;

    for(auto jet : jetCont->accepted()) {
      if (!jet) continue;
      count++;

      // Filling histos for angle relative to event plane
      Double_t phiMinusPsi2 = jet->Phi() - psi2V0[0];
      if (phiMinusPsi2 < 0.0) phiMinusPsi2 += TMath::TwoPi();
      histName = TString::Format("%s/hJetPhiMinusPsi2_%d", groupName.Data(), fCentBin);
      fHistManager.FillTH1(histName, phiMinusPsi2);
      Double_t phiMinusPsi3 = jet->Phi() - psi3V0[0];
      if (phiMinusPsi3 < 0.0) phiMinusPsi3 += TMath::TwoPi();
      histName = TString::Format("%s/hJetPhiMinusPsi3_%d", groupName.Data(), fCentBin);
      fHistManager.FillTH1(histName, phiMinusPsi3);

      //if (jetCont->GetRhoParameter()) {
      Double_t localRhoVal = 0.0;
      Double_t localRhoValScaled = 0.0;
      Double_t jetPtCorr = 0.0;
      Double_t jetPtCorrLocal = 0.0;
      Double_t deltaPhiJetEP = -999.0;

      deltaPhiJetEP = jet->Phi() - psi2V0[0];
      if (jet->Phi() - psi2V0[0] < 0.0) deltaPhiJetEP += TMath::TwoPi();

      localRhoVal = fLocalRho->GetLocalVal(jet->Phi(), GetJetContainer()->GetJetRadius(), fLocalRho->GetVal());
      
      // localRhoValScaled = localRhoVal * jetCont->GetRhoVal() / p0;
      jetPtCorr = jet->Pt() - jetCont->GetRhoVal() * jet->Area();
      localRhoValScaled = fLocalRho->GetLocalVal(\
        jet->Phi(), GetJetContainer()->GetJetRadius(), fLocalRho->GetVal());
      jetPtCorrLocal = jet->Pt() - localRhoValScaled * jet->Area();
      
      if(0){
        std::cout << "(jetPt, jetPhi, jetA, globalRho, localRho, localRhoValScal) = (" \
          << jet->Pt() << ", " << jet->Phi() << ", " << jet->Area() << ", " \
          << jetCont->GetRhoVal() * jet->Area() << ", " << localRhoValScaled * jet->Area() \
          << ", " << localRhoValScaled 
          << ")" << std::endl;
      }
      
      TString rhoGroupName = TString::Format("%s/Rho", groupName.Data());
      // histName = TString::Format("%s/hRhoVsCent", rhoGroupName.Data());
      // fHistManager.FillTH2(histName.Data(), fCent, jetCont->GetRhoVal());
      histName = TString::Format("%s/hJetRho_%d", rhoGroupName.Data(), fCentBin);
      fHistManager.FillTH1(histName, jetCont->GetRhoVal());
      histName = TString::Format("%s/hJetRhoLocal_%d", rhoGroupName.Data(), fCentBin);
      fHistManager.FillTH1(histName, localRhoVal); // trying out local rho val

      histName = TString::Format("%s/hJetLocalRhoVsAverageRho_%d",rhoGroupName.Data(), fCentBin);
      fHistManager.FillTH2(histName, jetCont->GetRhoVal(), localRhoValScaled);
      histName = TString::Format("%s/hJetCorrPtLocalVsJetCorrPt_%d",rhoGroupName.Data(), fCentBin);
      fHistManager.FillTH2(histName, jetPtCorr, jetPtCorrLocal);
      histName = TString::Format("%s/hJetRhoVsDeltaPhi_%d",rhoGroupName.Data(),fCentBin);
      fHistManager.FillTH2(histName, deltaPhiJetEP, localRhoValScaled);
      
      
      //inclusive Jet
      TString InclusiveGroupName = TString::Format("%s/Inclusive", groupName.Data());
      histName = TString::Format("%s/hJetArea_%d", InclusiveGroupName.Data(), fCentBin);
      fHistManager.FillTH1(histName, jet->Area());
      histName = TString::Format("%s/hJetPhi_%d", InclusiveGroupName.Data(), fCentBin);
      fHistManager.FillTH1(histName, jet->Phi());
      histName = TString::Format("%s/hJetEta_%d", InclusiveGroupName.Data(), fCentBin);
      fHistManager.FillTH1(histName, jet->Eta());
      histName = TString::Format("%s/hJetPt_%d", InclusiveGroupName.Data(), fCentBin);

      fHistManager.FillTH1(histName, jet->Pt());
      histName = TString::Format("%s/hJetCorrPt_%d", InclusiveGroupName.Data(), fCentBin);
      fHistManager.FillTH1(histName, jetPtCorr);
      histName = TString::Format("%s/hJetCorrPtLocal_%d", InclusiveGroupName.Data(), fCentBin);
      fHistManager.FillTH1(histName, jetPtCorrLocal);
      
      //V2 In plane Jet
      if ((phiMinusPsi2 < TMath::Pi()/4) || (phiMinusPsi2 >= 7*TMath::Pi()/4)\
      || (phiMinusPsi2 >= 3*TMath::Pi()/4 && phiMinusPsi2 < 5*TMath::Pi()/4)) {
        TString IPlaneGroupName = TString::Format("%s/InPlane", groupName.Data());

        histName = TString::Format("%s/hJetPt_%d", IPlaneGroupName.Data(), fCentBin);
        fHistManager.FillTH1(histName, jet->Pt());
        histName = TString::Format("%s/hCorrJetPt_%d", IPlaneGroupName.Data(), fCentBin);
        fHistManager.FillTH1(histName, jetPtCorrLocal);
      }
      else {
        TString OPlaneGroupName = TString::Format("%s/InPlane", groupName.Data());
        histName = TString::Format("%s/hJetPt_%d", OPlaneGroupName.Data(), fCentBin);
        fHistManager.FillTH1(histName, jet->Pt());
        histName = TString::Format("%s/hCorrJetPt_%d", OPlaneGroupName.Data(), fCentBin);
        fHistManager.FillTH1(histName, jetPtCorrLocal);
      }
      
      if(leadingJetPt < jetPtCorrLocal){
        leadingJetEta = jet->Eta();
        leadingJetPhi = jet->Phi();
        leadingJetPt  = jetPtCorrLocal;
      }

    }

    Double_t rcPt = 0., rcEta = 0., rcPhi = 0.;
    CalcRandomCone(rcPt, rcEta, rcPhi, leadingJetEta, leadingJetPhi, jetR);    
    Double_t rcLocalRhoValScaled = fLocalRho->GetLocalVal(rcPhi, jetR, fLocalRho->GetVal());
    Double_t deltaLoacalPt = rcPt - rcLocalRhoValScaled*jetR*jetR*TMath::Pi();
    Double_t deltaGlobalPt = rcPt - fLocalRho->GetVal()*jetR*jetR*TMath::Pi();
    histName = TString::Format("%s/hDeltaPt_%d", groupName.Data(), fCentBin);
    fHistManager.FillTH1(histName, deltaGlobalPt);
    histName = TString::Format("%s/hDeltaPt_Local_%d", groupName.Data(), fCentBin);
    fHistManager.FillTH1(histName, deltaLoacalPt);
    histName = TString::Format("%s/hPhiVsDeltaPt_Local_%d", groupName.Data(), fCentBin);
    fHistManager.FillTH2(histName, rcPhi, deltaLoacalPt);

    histName = TString::Format("%s/hNJets_%d", groupName.Data(), fCentBin);
    fHistManager.FillTH1(histName, count);
  }
}


/**
 * This function performs matching of det-level jets to truth-level jets, and fills relevant histograms.
 * There are two matching approaches implemented, each of which expect a certain set of jet containers to be attached:
 * (1) fDoJetMatchingGeometrical: Do geometrical matching, with two jet containers fMCJetContainer and fDetJetContainer.
 *     This is appropriate for pp and p-Pb.
 * (2) fDoJetMatchingMCFraction: Do MC fraction based jet matching, with pp-truth, pp-det, and combined-det jet containers.
 *     This is appropriate for Pb-Pb.
 */
void AliAnalysisTaskEmbeddingJetWithEP::FillMatchedJetHistograms()
{
  // Loop over all jets and fill the ClosestJet(), i.e. the matching candidate.
  // Note: Allow truth jets to be outside of EMCALfid or fail 5 GeV requirement, since these can still contribute accepted det-jets
  //       (but for the jet reconstruction efficiency denominator, the criteria should be enforced).
  if (fDoJetMatchingGeom) {
    if (fRequireMatchedJetAcc) {
      fMCJetContainer->SetJetAcceptanceType(AliEmcalJet::kTPC);
      ComputeJetMatches(fDetJetContainer, fMCJetContainer, kTRUE);
    }
    else {
      ComputeJetMatches(fDetJetContainer, fMCJetContainer, kFALSE);
    }
  }
  else if (fDoJetMatchingMCFraction) {
    // First match PbPb-det to pp-det
    ComputeJetMatches(fDetJetContainer, fDetJetContainerPPIntermediate, kTRUE);
    
    // Then match pp-det to pp-truth
    if (fRequireMatchedJetAcc) { // Require pp-truth be accepted (i.e. leading track req), but still allow geometrical acceptance migration
      fMCJetContainer->SetJetAcceptanceType(AliEmcalJet::kTPC);
      ComputeJetMatches(fDetJetContainerPPIntermediate, fMCJetContainer, kTRUE);
    }
    else{ // Don't require pp-truth jet to be accepted
      ComputeJetMatches(fDetJetContainerPPIntermediate, fMCJetContainer, kFALSE);
    }

  }
  
  // Loop through accepted det-level jets, and retrieve matching candidate.
  // It match passes criteria (i.e. matching distance, uniqueness, MC fraction), fill matching histos.
  Double_t rhoVal = 0;
  if (fDetJetContainer->GetRhoParameter()) rhoVal = fDetJetContainer->GetRhoVal();
  
  for (auto jet : fDetJetContainer->accepted()) {

    // == s == Estimate Jet angle from EP   ====================================
    std::vector<TString> fHistJetEPName = {"Inclusive", "InPlane", "OutOfPlane"}; ///<
    TString groupName = "";

    Double_t phiMinusPsi2 = jet->Phi() - psi2V0[0];
    if (phiMinusPsi2 < 0.0) phiMinusPsi2 += TMath::TwoPi();
    Double_t phiMinusPsi3 = jet->Phi() - psi3V0[0];
    if (phiMinusPsi3 < 0.0) phiMinusPsi3 += TMath::TwoPi();

    //V2 In plane Jet
    if ((phiMinusPsi2 < TMath::Pi()/4) || (phiMinusPsi2 >= 7*TMath::Pi()/4)\
    || (phiMinusPsi2 >= 3*TMath::Pi()/4 && phiMinusPsi2 < 5*TMath::Pi()/4)) {
      groupName = "MatchedJetHistograms/InPlane";
    }else groupName = "MatchedJetHistograms/OutOfPlane";
    // == e == Estimate Jet angle from EP   ====================================
    

    Float_t jetPtCorr = jet->Pt() - fDetJetContainer->GetRhoVal() * jet->Area();
    Float_t localRhoValScaled = fLocalRho->GetLocalVal(\
    jet->Phi(), GetJetContainer()->GetJetRadius(), fLocalRho->GetVal());
    Float_t jetPtCorrLocal = jet->Pt() - localRhoValScaled * jet->Area();
    // Float_t detPt = GetJetPt(jet, rhoVal);
    Float_t detPt = jetPtCorrLocal;
    
    // Get the matched part-level jet
    const AliEmcalJet* matchedPartLevelJet = GetMatchedPartLevelJet(jet, detPt);
    if (!matchedPartLevelJet) continue;
  
    Float_t truthPt = matchedPartLevelJet->Pt();
    
    // Fill response matrix (centrality, pT-truth, pT-det)
    TString histname;
    //This is a 5-dim RM with information on the angularity and matching distance
    if (fDoDifferentialRM) {
      // (1) pt part LvL, (2) pt det LvL, (3) Matching distance (4) angularity (5) centrality
      Double_t angularity     = GetAngularity(matchedPartLevelJet);
      Double_t matchDistance  = matchedPartLevelJet->ClosestJetDistance();
      Double_t sharedFraction = fDetJetContainer->GetFractionSharedPt(jet, nullptr);
      Double_t x[6] = {truthPt, detPt, matchDistance, sharedFraction, angularity, fCent};

      histname = "MatchedJetHistograms/Inclusive/hResponseMatrixDiff";
      fHistManager.FillTHnSparse(histname, x);
      histname = TString::Format("%s/hResponseMatrixDiff",groupName.Data());
      fHistManager.FillTHnSparse(histname, x);
    }
    //This is a 3D RM for PbPb and a 2D RM for pp
    else {
      histname = "MatchedJetHistograms/Inclusive/hResponseMatrix";
      fHistManager.FillTH3(histname, fCent, truthPt, detPt);
      histname = TString::Format("%s/hResponseMatrix",groupName.Data());
      fHistManager.FillTH3(histname, fCent, truthPt, detPt);
    }

    // Fill JES shift (centrality, pT-truth, (pT-det - pT-truth) / pT-truth)
    histname = "MatchedJetHistograms/Inclusive/hJESshift";
    fHistManager.FillTH3(histname, fCent, truthPt, (detPt-truthPt)/truthPt );
    histname = TString::Format("%s/hJESshift",groupName.Data());
    fHistManager.FillTH3(histname, fCent, truthPt, (detPt-truthPt)/truthPt);
    
    // Fill NEF of det-level matched jets (centrality, pT-truth, NEF)
    histname = "MatchedJetHistograms/Inclusive/hNEFVsPt";
    fHistManager.FillTH3(histname, fCent, truthPt, jet->NEF());
    histname = TString::Format("%s/hNEFVsPt",groupName.Data());
    fHistManager.FillTH3(histname, fCent, truthPt, jet->NEF());
    
    // Fill z-leading (charged) of det-level matched jets (centrality, pT-truth, z-leading)
    TLorentzVector leadPart;
    fDetJetContainer->GetLeadingHadronMomentum(leadPart, jet);
    Double_t z = GetParallelFraction(leadPart.Vect(), jet);
    if(z == 1 || (z>1 && z-1 < 1e-3)) z = 0.999; // so that it will contribute to the bin <1
    histname = "MatchedJetHistograms/Inclusive/hZLeadingVsPt";
    fHistManager.FillTH3(histname, fCent, truthPt, z);
    histname = TString::Format("%s/hZLeadingVsPt",groupName.Data());
    fHistManager.FillTH3(histname, fCent, truthPt, z);
    
  } //jet loop
}


Bool_t  AliAnalysisTaskEmbeddingJetWithEP::QnGainCalibration(){
  TString histName;
  TString groupName;
  groupName="EventPlane";
  
  fCalibQA = kTRUE;
  if(fCalibRefFile){
    TList *tempCalibRefObj = (TList *)fCalibRefFile->Get("fWgtsV0ZDC");
    fCalibRefObjList = tempCalibRefObj;

    //V0 Channel Gains:
    fHCorrV0ChWeghts = (TH2F *)fCalibRefObjList->FindObject(Form("hWgtV0ChannelsvsVzRun%d",fRunNumber));
    if(fHCorrV0ChWeghts){
      printf("\n ===========> Info:: V0 Channel Weights Found for Run %d \n ",fRunNumber);
    }
  }
  AliAODVZERO* fAodV0 = dynamic_cast<AliAODVZERO*>(fAOD->GetVZEROData());

  AliVEvent *fVevent = dynamic_cast<AliVEvent*>(InputEvent());
  const AliVVertex *pointVtx = fVevent->GetPrimaryVertex();
  Double_t fVtxZ = -999;
  fVtxZ  = pointVtx->GetZ();
  
  Int_t ibinV0 = 0;
  Double_t fSumMV0A = 0.;
  Double_t fSumMV0C = 0.;
  Double_t fSumMV0M = 0.;
  Double_t fV0chGain = 0.;
  Double_t fMultV0 = 0.;
  
  for(int iV0 = 0; iV0 < 64; iV0++) { //0-31 is V0C, 32-63 VOA
    fMultV0 = fAodV0->GetMultiplicity(iV0);
    
    /// V0 Channel Gain Correction:
    if(fHCorrV0ChWeghts){
      ibinV0    = fHCorrV0ChWeghts->FindBin(fVtxZ,iV0);
      fV0chGain = fHCorrV0ChWeghts->GetBinContent(ibinV0);
    }
    
    fMultV0 = fMultV0*fV0chGain;   //Corrected Multiplicity
    // std::cout << "(fV0chGain, fMultV0) = (" << fV0chGain << ","  << fMultV0 << ")" << std::endl;

    Double_t fPhiV0  = TMath::PiOver4()*(0.5 + iV0 % 8);
    
    q2VecV0M[0] += TMath::Cos(2*fPhiV0) * fMultV0;
    q2VecV0M[1] += TMath::Sin(2*fPhiV0) * fMultV0;
    q3VecV0M[0] += TMath::Cos(3*fPhiV0) * fMultV0;
    q3VecV0M[1] += TMath::Sin(3*fPhiV0) * fMultV0;
    fSumMV0M += fMultV0;

    if(iV0 < 32){
      q2VecV0C[0] += TMath::Cos(2*fPhiV0) * fMultV0;
      q2VecV0C[1] += TMath::Sin(2*fPhiV0) * fMultV0;
      q3VecV0C[0] += TMath::Cos(3*fPhiV0) * fMultV0;
      q3VecV0C[1] += TMath::Sin(3*fPhiV0) * fMultV0;
      fSumMV0C += fMultV0;
    }
    else if(iV0 >= 32){
      q2VecV0A[0] += TMath::Cos(2*fPhiV0) * fMultV0;
      q2VecV0A[1] += TMath::Sin(2*fPhiV0) * fMultV0;
      q3VecV0A[0] += TMath::Cos(3*fPhiV0) * fMultV0;
      q3VecV0A[1] += TMath::Sin(3*fPhiV0) * fMultV0;
      fSumMV0A += fMultV0;
    }

    
  }///V0 Channel loop
  
  /// Now the q vectors:
  if(fSumMV0A<=1e-4 || fSumMV0C<=1e-4){
    q2VecV0M[0] = 0.;
    q2VecV0M[1] = 0.;
    q3VecV0M[0] = 0.;
    q3VecV0M[1] = 0.;
    q2VecV0C[0] = 0.;
    q2VecV0C[1] = 0.;
    q3VecV0C[0] = 0.;
    q3VecV0C[1] = 0.;
    q2VecV0A[0] = 0.;
    q2VecV0A[1] = 0.;
    q3VecV0A[0] = 0.;
    q3VecV0A[1] = 0.;
    
    return kFALSE;       
  }
  else{
    q2VecV0M[0] = q2VecV0M[0]/fSumMV0M;
    q2VecV0M[1] = q2VecV0M[1]/fSumMV0M;
    q3VecV0M[0] = q3VecV0M[0]/fSumMV0M;
    q3VecV0M[1] = q3VecV0M[1]/fSumMV0M;
    q2VecV0C[0] = q2VecV0C[0]/fSumMV0C;
    q2VecV0C[1] = q2VecV0C[1]/fSumMV0C;
    q3VecV0C[0] = q3VecV0C[0]/fSumMV0C;
    q3VecV0C[1] = q3VecV0C[1]/fSumMV0C;
    q2VecV0A[0] = q2VecV0A[0]/fSumMV0A;
    q2VecV0A[1] = q2VecV0A[1]/fSumMV0A;
    q3VecV0A[0] = q3VecV0A[0]/fSumMV0A;
    q3VecV0A[1] = q3VecV0A[1]/fSumMV0A;
    
    return kTRUE;  
  }



}

Bool_t  AliAnalysisTaskEmbeddingJetWithEP::QnRecenteringCalibration(){
  TList *tempCalibRefObj = (TList *)fCalibRefFile->Get("fWgtsV0ZDC");
  fCalibRefObjList = tempCalibRefObj;

  //Get V0A, V0C <Q> Vectors:
  fHCorrQ2xV0C = (TH1D *)fCalibRefObjList->FindObject(Form("fHisAvgQNxvsCentV0CRun%d",fRunNumber));
  fHCorrQ2yV0C = (TH1D *)fCalibRefObjList->FindObject(Form("fHisAvgQNyvsCentV0CRun%d",fRunNumber));    
  fHCorrQ2xV0A = (TH1D *)fCalibRefObjList->FindObject(Form("fHisAvgQNxvsCentV0ARun%d",fRunNumber));
  fHCorrQ2yV0A = (TH1D *)fCalibRefObjList->FindObject(Form("fHisAvgQNyvsCentV0ARun%d",fRunNumber));
	
  fHCorrQ3xV0C = (TH1D *)fCalibRefObjList->FindObject(Form("fHisAvgQ3xvsCentV0CRun%d",fRunNumber));
  fHCorrQ3yV0C = (TH1D *)fCalibRefObjList->FindObject(Form("fHisAvgQ3yvsCentV0CRun%d",fRunNumber));    
  fHCorrQ3xV0A = (TH1D *)fCalibRefObjList->FindObject(Form("fHisAvgQ3xvsCentV0ARun%d",fRunNumber));
  fHCorrQ3yV0A = (TH1D *)fCalibRefObjList->FindObject(Form("fHisAvgQ3yvsCentV0ARun%d",fRunNumber));    
  if(fHCorrQ2xV0C && fHCorrQ2yV0C && fHCorrQ2xV0A && fHCorrQ2yV0A){
    printf(" ===========> Info:: V0A,V0C <Q> Found for Run %d \n ",fRunNumber);
  }

  Int_t icentbin = 0;
  Double_t avgqx=0,avgqy=0; 
  //cout<<" => Before qnxV0C "<<qnxV0C<<"\tqnyV0C "<<qnyV0C<<"\tqnxV0A "<<qnxV0A<<"\tqnyV0A "<<qnyV0A<<endl;
  
  if(fHCorrQ2xV0C && fHCorrQ2yV0C){
    icentbin = fHCorrQ2xV0C->FindBin(fCent);
    avgqx = fHCorrQ2xV0C->GetBinContent(icentbin);
    avgqy = fHCorrQ2yV0C->GetBinContent(icentbin);
    q2VecV0C[0] -= avgqx;
    q2VecV0C[1] -= avgqy;	
    //cout<<" V0C PsiN: "<<gPsiN<<" Cent: "<<fCent<<"\t <qx> "<<avgqx<<"\t <qy> "<<avgqy<<endl;
  }
  if(fHCorrQ2xV0A && fHCorrQ2yV0A){
    icentbin = fHCorrQ2xV0A->FindBin(fCent);
    avgqx = fHCorrQ2xV0A->GetBinContent(icentbin);
    avgqy = fHCorrQ2yV0A->GetBinContent(icentbin);
    q2VecV0A[0] -= avgqx;
    q2VecV0A[1] -= avgqy;
    //cout<<" V0A PsiN: "<<gPsiN<<" Cent: "<<fCent<<"\t <qx> "<<avgqx<<"\t <qy> "<<avgqy<<endl;
  }
  //cout<<" => After qnxV0C "<<qnxV0C<<"\tqnyV0C "<<qnyV0C<<" qnxV0A"<<qnxV0A<<"\tqnyV0A"<<qnyV0A<<endl;
  if(fHCorrQ3xV0C && fHCorrQ3yV0C){
    icentbin = fHCorrQ3xV0C->FindBin(fCent);
    avgqx = fHCorrQ3xV0C->GetBinContent(icentbin);
    avgqy = fHCorrQ3yV0C->GetBinContent(icentbin);
    q3VecV0C[0] -= avgqx;
    q3VecV0C[1] -= avgqy;      
    //cout<<" V0C PsiN: "<<gPsiN<<" Cent: "<<fCent<<"\t <qx> "<<avgqx<<"\t <qy> "<<avgqy<<endl;
  }
  if(fHCorrQ3xV0A && fHCorrQ3yV0A){
    icentbin = fHCorrQ3xV0A->FindBin(fCent);
    avgqx = fHCorrQ3xV0A->GetBinContent(icentbin);
    avgqy = fHCorrQ3yV0A->GetBinContent(icentbin);
    q3VecV0A[0] -= avgqx;
    q3VecV0A[1] -= avgqy;           
    //cout<<" V0A PsiN: "<<gPsiN<<" Cent: "<<fCent<<"\t <qx> "<<avgqx<<"\t <qy> "<<avgqy<<endl;
  }
  //cout<<" => After qnxV0C "<<qnxV0C<<"\tqnyV0C "<<qnyV0C<<" qnxV0A "<<qnxV0A<<"\tqnyV0A "<<qnyV0A<<endl;
  
  return kTRUE;
}

//_____________________________________________________________________________
TH1F* AliAnalysisTaskEmbeddingJetWithEP::GetResoFromOutputFile(detectorType det, Int_t h, TArrayD* cen)
{
    if(!fOutputList) {
        printf(" > Please add fOutputList first < \n");
        return 0x0;
    }
    TH1F* r(0x0);
    (cen) ? r = new TH1F("R", "R", cen->GetSize()-1, cen->GetArray()) : r = new TH1F("R", "R", 10, 0, 10);
    if(!cen) r->GetXaxis()->SetTitle("number of centrality bin");
    r->GetYaxis()->SetTitle(Form("Resolution #Psi_{%i}", h));
    for(Int_t i(0); i < 10; i++) {
        TProfile* temp((TProfile*)fOutputList->FindObject(Form("fProfV%iResolution_%i", h, i)));
        if(!temp) break;
        Double_t a(temp->GetBinContent(3)); //cos(2[psi_V0A - psi_V0C])
        Double_t b(temp->GetBinContent(5)); //cos(2[psi_TPC - psi_V0C])
        Double_t c(temp->GetBinContent(7)); //cos(2[psi_TPC - psi_V0A])

        Double_t d(temp->GetBinContent(9));  //cos(2[psi_V0M - psi_TPCnega])
        Double_t e(temp->GetBinContent(10)); //cos(2[psi_V0M - psi_TPCposi])
        Double_t f(temp->GetBinContent(11)); //cos(2[psi_TPCnega - psi_TPCposi])

        Double_t _a(temp->GetBinError(3)), _b(temp->GetBinError(5)), _c(temp->GetBinError(7));
        Double_t _d(temp->GetBinError(9)), _e(temp->GetBinError(10)), _f(temp->GetBinError(11));
        Double_t error(0);
        if(a <= 0 || b <= 0 || c <= 0 || d <= 0 || e <= 0 || f <= 0) continue;
        switch (det) {
            case kVZEROA : {
                r->SetBinContent(1+i, TMath::Sqrt((a*b)/c));
                if(i==0) r->SetNameTitle("VZEROA resolution", "VZEROA resolution");
                error = TMath::Power((2.*a*TMath::Sqrt((a*b)/c))/3.,2.)*_a*_a+TMath::Power((2.*b*TMath::Sqrt((a*b)/c))/3.,2.)*_b*_b+TMath::Power(2.*c*TMath::Sqrt((a*b)/c),2.)*_c*_c;
                if(error > 0.) error = TMath::Sqrt(error);
                r->SetBinError(1+i, error);
            } break;
            case kVZEROC : {
                r->SetBinContent(1+i, TMath::Sqrt((a*c)/b));
                error = TMath::Power((2.*a*TMath::Sqrt((a*c)/b))/3.,2.)*_a*_a+TMath::Power((2.*b*TMath::Sqrt((a*c)/b)),2.)*_b*_b+TMath::Power(2.*c*TMath::Sqrt((a*c)/b)/3.,2.)*_c*_c;
                if(error > 0.) error = TMath::Sqrt(error);
                if(i==0) r->SetNameTitle("VZEROC resolution", "VZEROC resolution");
                r->SetBinError(1+i, error);
            } break;
            case kTPC : {
                r->SetBinContent(1+i, TMath::Sqrt((b*c)/a));
                if(i==0) r->SetNameTitle("TPC resolution", "TPC resolution");
                r->SetBinError(1+i, TMath::Sqrt(_a*_a+_b*_b+_c*_c));
            } break;
            case kVZEROComb : {
                r->SetBinContent(1+i, TMath::Sqrt((d*e)/f));
                if(i==0) r->SetNameTitle("VZEROComb resolution", "VZEROComb resolution");
                r->SetBinError(1+i, TMath::Sqrt(_d*_d+_e*_e+_f*_f));
            } break;
            default : break;
        }
    }
    return r;
}

//_____________________________________________________________________________
Double_t AliAnalysisTaskEmbeddingJetWithEP::CalculateEventPlaneChi(Double_t res)
{
    // return chi for given resolution to combine event plane estimates from two subevents
    // see Phys. Rev. C no. CS6346 (http://arxiv.org/abs/nucl-ex/9805001)
    Double_t chi(2.), delta(1.), con((TMath::Sqrt(TMath::Pi()))/(2.*TMath::Sqrt(2)));
    for (Int_t i(0); i < 15; i++) {
        chi = ((con*chi*TMath::Exp(-chi*chi/4.)*(TMath::BesselI0(chi*chi/4.)+TMath::BesselI1(chi*chi/4.))) < res) ? chi + delta : chi - delta;
        delta = delta / 2.;
    }
    return chi;
}


//_____________________________________________________________________________
void AliAnalysisTaskEmbeddingJetWithEP::CalcRandomCone(Double_t &pt, Double_t &eta, Double_t &phi,
  Double_t &leadingJetEta, Double_t &leadingJetPhi, Double_t &jetR) const
{
    // get a random cone
    pt = 0; eta = 0; phi = 0;
    Double_t dJet(999);// no jet: same as jet very far away

    // the random cone acceptance has to equal the jet acceptance
    // this also insures safety when runnnig on the semi-good tpc runs for 11h data,
    // where jet acceptance is adjusted to reduced acceptance - hence random cone acceptance as well
    Float_t minPhi(GetJetContainer()->GetJetPhiMin()), maxPhi(GetJetContainer()->GetJetPhiMax());
    if(maxPhi > TMath::TwoPi()) maxPhi = TMath::TwoPi();
    if(minPhi < 0 ) minPhi = 0.;
    
    
    // construct a random cone and see if it's far away enough from the leading jet
    Int_t attempts(1000);
    while(kTRUE) {
        attempts--;
        eta = gRandom->Uniform(GetJetContainer()->GetJetEtaMin(), GetJetContainer()->GetJetEtaMax());
        phi = gRandom->Uniform(minPhi, maxPhi);

        dJet = TMath::Sqrt((leadingJetEta-eta)*(leadingJetEta-eta)\
          +(leadingJetPhi-phi)*(leadingJetPhi-phi));
        if(dJet > 2*jetR) break;
        else if (attempts == 0) {
            printf(" > No random cone after 1000 tries, giving up ... !\n");
            return;
        }
    }
    // get the charged energy (if tracks are provided)
    AliParticleContainer* tracksCont = 0;
    TIter next(&fParticleCollArray);
    while ((tracksCont = static_cast<AliParticleContainer*>(next()))) {
        // std::cout << "particle Contaner in RC jet Estimate: " << tracksCont->GetName() << std::endl;
        tracksCont->ResetCurrentID();
        AliVParticle* track = tracksCont->GetNextAcceptParticle();

        for(auto track : tracksCont->accepted()) {
            Float_t etaTrack(track->Eta()), phiTrack(track->Phi());
            // get distance from cone
            if(TMath::Abs(phiTrack-phi) > TMath::Abs(phiTrack - phi + TMath::TwoPi()))\
              phiTrack+=TMath::TwoPi();
            if(TMath::Abs(phiTrack-phi) > TMath::Abs(phiTrack - phi - TMath::TwoPi()))\
              phiTrack-=TMath::TwoPi();
            
            Float_t rangeR = TMath::Sqrt(TMath::Abs((etaTrack-eta)*(etaTrack-eta)\
              +(phiTrack-phi)*(phiTrack-phi)));
            if(rangeR <= jetR) pt += track->Pt();
        }
    }
    
}

void AliAnalysisTaskEmbeddingJetWithEP::SetupPileUpRemovalFunctions(){
  
  ////==========> LHC18q/r PileUp Removal Functions: ---- Do not Remove them !!! -----
  Double_t parV0[8] = {43.8011, 0.822574, 8.49794e-02, 1.34217e+02, 7.09023e+00, 4.99720e-02, -4.99051e-04, 1.55864e-06};
  fV0CutPU  = new TF1("fV0CutPU", "[0]+[1]*x - 6.*[2]*([3] + [4]*sqrt(x) + [5]*x + [6]*x*sqrt(x) + [7]*x*x)", 0, 100000);
  fV0CutPU->SetParameters(parV0);
  
  fSPDCutPU = new TF1("fSPDCutPU", "400. + 4.*x", 0, 10000);

  Double_t parFB32[8] = {2093.36, -66.425, 0.728932, -0.0027611, 1.01801e+02, -5.23083e+00, -1.03792e+00, 5.70399e-03};
  fMultCutPU = new TF1("fMultCutPU", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 6.*([4]+[5]*sqrt(x)+[6]*x+[7]*x*x)", 0, 90);
  fMultCutPU->SetParameters(parFB32);
  
  Double_t parV0CL0[6] = {0.320462, 0.961793, 1.02278, 0.0330054, -0.000719631, 6.90312e-06};
  fCenCutLowPU  = new TF1("fCenCutLowPU", "[0]+[1]*x - 6.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)",  0, 100);
  fCenCutLowPU->SetParameters(parV0CL0);
  
  fCenCutHighPU = new TF1("fCenCutHighPU", "[0]+[1]*x + 5.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
  fCenCutHighPU->SetParameters(parV0CL0);
  //--------------------------------------------------------------------------------------

}


Bool_t AliAnalysisTaskEmbeddingJetWithEP::CheckEventIsPileUp2018(){
  
  /// Todo Rihan: I can check for PileUp and get TPC event Plane in Same Function
  /// Utilizing same track loop. This method would save time..
  if (!fAOD && AODEvent() && IsStandardAOD()) {
    fAOD = dynamic_cast<AliAODEvent*> (AODEvent());
  }
  if (!fAOD) {
    AliWarning("AliAnalysisTaskJetQnVectors::Exec(): bad AOD");
    return kFALSE;
  }
  
  Bool_t BisPileup = kFALSE;

  Double_t centrV0M=-99.0;
  Double_t centrCL1=-99.0;
  Double_t centrCL0=-99.0;

  AliMultSelection* fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");
  if(!fMultSelection) {
    printf("\n\n **WARNING** ::UserExec() AliMultSelection object not found.\n\n");
    exit(111);
  }
  
  centrV0M = fMultSelection->GetMultiplicityPercentile("V0M");
  centrCL1 = fMultSelection->GetMultiplicityPercentile("CL1");
  centrCL0 = fMultSelection->GetMultiplicityPercentile("CL0");

  Int_t nITSClsLy0 = fAOD->GetNumberOfITSClusters(0);
  Int_t nITSClsLy1 = fAOD->GetNumberOfITSClusters(1);
  Int_t nITSCls = nITSClsLy0 + nITSClsLy1;

  AliAODTracklets* aodTrkl = (AliAODTracklets*)fAOD->GetTracklets();
  Int_t nITSTrkls = aodTrkl->GetNumberOfTracklets();
  const Int_t nTracks = fAOD->GetNumberOfTracks();
  Int_t multTrk = 0;
  
  for (Int_t it = 0; it < nTracks; it++) {
    AliAODTrack* aodTrk = (AliAODTrack*)fAOD->GetTrack(it);

    if (!aodTrk){
      delete aodTrk;
      continue;
    }

    if(aodTrk->TestFilterBit(32)){
      if((TMath::Abs(aodTrk->Eta()) < 0.8)&&(aodTrk->GetTPCNcls() >= 70)&&(aodTrk->Pt() >= 0.2))
      multTrk++;
    }
  }
  
  AliAODVZERO* aodV0 = fAOD->GetVZEROData();
  Float_t multV0a = aodV0->GetMTotV0A();
  Float_t multV0c = aodV0->GetMTotV0C();
  Float_t multV0Tot = multV0a + multV0c;
  UShort_t multV0aOn = aodV0->GetTriggerChargeA();
  UShort_t multV0cOn = aodV0->GetTriggerChargeC();
  UShort_t multV0On = multV0aOn + multV0cOn;
  
  Int_t tpcClsTot = fAOD->GetNumberOfTPCClusters();
  Float_t nclsDif = Float_t(tpcClsTot) \
    - (60932.9 + 69.2897*multV0Tot - 0.000217837*multV0Tot*multV0Tot);
  
  
  if(centrCL0 < fCenCutLowPU->Eval(centrV0M)) BisPileup=kTRUE;
  if(centrCL0 > fCenCutHighPU->Eval(centrV0M)) BisPileup=kTRUE;
  if(Float_t(nITSCls) > fSPDCutPU->Eval(nITSTrkls)) BisPileup=kTRUE;
  if(multV0On < fV0CutPU->Eval(multV0Tot)) BisPileup=kTRUE;
  if(Float_t(multTrk) < fMultCutPU->Eval(centrV0M)) BisPileup=kTRUE;
  if(((AliAODHeader*)fAOD->GetHeader())->GetRefMultiplicityComb08() < 0) BisPileup=kTRUE;
  if(fAOD->IsIncompleteDAQ()) BisPileup=kTRUE;
  //if (nclsDif > 200000)//can be increased to 200000
  // BisPileup=kTRUE;

  Int_t multEsd = ((AliAODHeader*)fAOD->GetHeader())->GetNumberOfESDTracks();

  if(fPileupCutQA){
    TString histName;
    TString groupName;
    groupName="PileupCutQA";

    histName = TString::Format("%s/fHistCentCL0VsV0MBefore", groupName.Data());
    fHistManager.FillTH2(histName, centrV0M,centrCL0);
    histName = TString::Format("%s/fHistTPCVsESDTrkBefore", groupName.Data());
    fHistManager.FillTH2(histName, multTrk,multEsd);
    histName = TString::Format("%s/fHistTPConlyVsCL1Before", groupName.Data());
    fHistManager.FillTH2(histName, centrCL1,multTrk);
    histName = TString::Format("%s/fHistTPConlyVsV0MBefore", groupName.Data());
    fHistManager.FillTH2(histName, centrV0M,multTrk);

    if(!BisPileup){
      histName = TString::Format("%s/fHistCentCL0VsV0MAfter", groupName.Data());
      fHistManager.FillTH2(histName, centrV0M,centrCL0);
      histName = TString::Format("%s/fHistTPCVsESDTrkAfter", groupName.Data());
      fHistManager.FillTH2(histName, multTrk,multEsd);
      histName = TString::Format("%s/fHistTPConlyVsCL1After", groupName.Data());
      fHistManager.FillTH2(histName, centrCL1,multTrk);
      histName = TString::Format("%s/fHistTPConlyVsV0MAfter", groupName.Data());
      fHistManager.FillTH2(histName, centrV0M,multTrk);
    }
  }
  
  return BisPileup;
}


/*
 * Loop over jets of two specified jet collections, and fill the ClosestJet(), i.e. the matching candidate.
 * The first collection always uses the container acceptance criteria.
 * The second collection can be configured to use the container acceptance criteria or not.
 */
void AliAnalysisTaskEmbeddingJetWithEP::ComputeJetMatches(AliJetContainer* jetCont1, AliJetContainer* jetCont2, Bool_t bUseJetCont2Acc) {

  for (auto jet1 : jetCont1->all()) jet1->ResetMatching();
  for (auto jet2 : jetCont2->all()) jet2->ResetMatching();

  for (auto jet1 : jetCont1->accepted()) {//detector level
    if (jet1->Pt() < fDetJetMinMatchingPt) continue;
    
    if (bUseJetCont2Acc) {
      for (auto jet2 : jetCont2->accepted()) SetJetClosestCandidate(jet1, jet2);
    }else {
      for (auto jet2 : jetCont2->all()) {//truth level
        if (jet2->Pt() < fMCJetMinMatchingPt) continue;
        SetJetClosestCandidate(jet1, jet2);
      }
    }
  }

}


/*
 * Given two jets, set them as closest if they are closer than the current closest jets.
 */
void AliAnalysisTaskEmbeddingJetWithEP::SetJetClosestCandidate(AliEmcalJet* jet1, AliEmcalJet* jet2) {
  
  Double_t deltaR = jet1->DeltaR(jet2);
  if (deltaR > 0.) {
    if (deltaR < jet1->ClosestJetDistance()) jet1->SetClosestJet(jet2, deltaR);
    if (deltaR < jet2->ClosestJetDistance()) jet2->SetClosestJet(jet1, deltaR);
  }
}

/*
 * Return a pointer to the matched truth-level jet, if it passes the matching criteria.
 * For fDoJetMatchingGeometrical, this means (1) within R = fJetMatchingR, (2) unique match.
 * For fDoJetMatchingMCFraction, this means also shared MC fraction requirement. That is, if the jet (combined jet)
 * is matched to a pp det-level jet, which is matched to a pp truth-level jet -- it must satisfy
 * (1) The shared momentum fraction being larger than some minimum value fMinSharedMomentumFraction, and
 * (2) Their matched distance being below the max matching distance fMaxMatchedJetDistance
 *
 * @param[in] detJet det-level jet to be checked for a successful truth-level match.
 * @return Pointer to truth-level matched jet, if it exists. False otherwise.
*/
const AliEmcalJet* AliAnalysisTaskEmbeddingJetWithEP::GetMatchedPartLevelJet(const AliEmcalJet* detJet, Double_t detJetPt) {

  // Track in histogram how many matches pass each distinct matching criteria
  TString histNameQA = "MatchedJetHistograms/fHistJetMatchingQA";
  
  bool returnValue = false;
  const AliEmcalJet* partLevelJet = nullptr;
  
  // First, check if combined jet has a pp det-level match assigned
  if (detJet->ClosestJet()) {
    fHistManager.FillTH1(histNameQA, "matchedJet");
    returnValue = true;
    
    // Check shared momentum fraction.
    double sharedFraction = fDetJetContainer->GetFractionSharedPt(detJet, nullptr);
    if (sharedFraction < fMinSharedPFraction) returnValue = false;
    else fHistManager.FillTH1(histNameQA, "sharedPFraction");
  
    // Check that the combined jet has a particle-level match
    AliEmcalJet * detLevelJetPP = detJet->ClosestJet();
    partLevelJet = detLevelJetPP->ClosestJet();
    if (!partLevelJet) returnValue = false;
    else fHistManager.FillTH1(histNameQA, "partLevelMatchedJet");
    
    // Check the matching distance between the combined and pp det-level jets
    double matchedJetDistance = detJet->ClosestJetDistance();
    if (matchedJetDistance > fJetMatchingR) returnValue = false;
    else fHistManager.FillTH1(histNameQA, "jetDistancePPdet");
    
    
    // Check the matching distance between the combined and pp truth-level jets
    if (partLevelJet) {
      Double_t deltaR = detJet->DeltaR(partLevelJet);
      if (deltaR > fJetMatchingR) returnValue = false;
      else fHistManager.FillTH1(histNameQA, "jetDistancePPtruth");
    }
    
    // Record all cuts passed
    if (returnValue == true) fHistManager.FillTH1(histNameQA, "passedAllCuts");
    
    // Fill (det pT, shared MC fraction, deltaR) of closest jets
    TString histname = "MatchedJetHistograms/hMatchingDistanceVsMCFraction";
    fHistManager.FillTH3(histname, detJetPt, sharedFraction, matchedJetDistance);
    
  }
  else {
    fHistManager.FillTH1(histNameQA, "noMatch");
    returnValue = false;
  }
  
  if (returnValue) return partLevelJet;

  return 0;
}

/*
 * Compute the Angularity of a jet - based on particle tracks (for particle level info)
 */
Double_t AliAnalysisTaskEmbeddingJetWithEP::GetAngularity(const AliEmcalJet* jet)
{
  //
  Double_t angularity=-1;

  if (jet->GetNumberOfTracks()== 0) return 0;

  Double_t den =0.;
  Double_t num = 0.;
  AliVParticle *vp1 = 0x0;
  //loop over all tracks in the jet
  for(UInt_t i = 0; i < jet->GetNumberOfTracks(); i++) {
    vp1 = static_cast<AliVParticle*>(jet->Track(i));

    if (!vp1){
      Printf("AliVParticle associated to constituent not found");
      continue;
    }

    Double_t dphi = GetRelativePhi(vp1->Phi(),jet->Phi());
    Double_t dr2  = (vp1->Eta()-jet->Eta())*(vp1->Eta()-jet->Eta()) + dphi*dphi;
    Double_t dr   = TMath::Sqrt(dr2);
    num=num+vp1->Pt()*dr;
    den=den+vp1->Pt();
  }
  if (den>0) angularity=num/den;

  return angularity;
}

/*
 * Compute dPhi distance
 */
Double_t AliAnalysisTaskEmbeddingJetWithEP::GetRelativePhi(Double_t mphi, Double_t vphi){
  if (vphi < -1*TMath::Pi()) vphi += (2*TMath::Pi());
  else if (vphi > TMath::Pi()) vphi -= (2*TMath::Pi());
  if (mphi < -1*TMath::Pi()) mphi += (2*TMath::Pi());
  else if (mphi > TMath::Pi()) mphi -= (2*TMath::Pi());
  double dphi = mphi-vphi;
  if (dphi < -1*TMath::Pi()) dphi += (2*TMath::Pi());
  else if (dphi > TMath::Pi()) dphi -= (2*TMath::Pi());
  return dphi;//dphi in [-Pi, Pi]
}


/**
 * This function is called once at the end of the analysis.
 */
void AliAnalysisTaskEmbeddingJetWithEP::Terminate(Option_t *) 
{
}


/**
 * This function adds the task to the analysis manager. Often, this function is called
 * by an AddTask C macro. However, by compiling the code, it ensures that we do not
 * have to deal with difficulties caused by CINT.
 */
AliAnalysisTaskEmbeddingJetWithEP * AliAnalysisTaskEmbeddingJetWithEP::AddTaskEmbeddingJetWithEP(
  const char *ntracks, const char *nclusters, const char* ncells, const char *suffix)
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmbeddingJetWithEP", "No analysis manager to connect to.");
    return 0;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler)
  {
    ::Error("AddTaskEmbeddingJetWithEP", "This task requires an input event handler");
    return 0;
  }
  
  enum EDataType_t {kUnknown, kESD, kAOD};


  EDataType_t dataType = kAOD;
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------
  TString trackName(ntracks);
  if (trackName == "usedefault") trackName = "tracks";

  TString name("AliAnalysisTaskEmbeddingJetWithEP");
  if (strcmp(suffix,"") != 0) {
    name += "_";
    name += suffix;
  }

  AliAnalysisTaskEmbeddingJetWithEP* rawJetTask = new AliAnalysisTaskEmbeddingJetWithEP(name);
  // rawJetTask->LoadSpliForqnPerce(qnSplineFileName); //new
  rawJetTask->SetVzRange(-10,10);



  if (trackName == "mcparticles") rawJetTask->AddMCParticleContainer(trackName);
  else if (trackName == "tracks") rawJetTask->AddTrackContainer(trackName);
  else if (!trackName.IsNull()) rawJetTask->AddParticleContainer(trackName);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  mgr->AddTask(rawJetTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  TString contname(name);
  contname += "_histos";
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(),
      TList::Class(),AliAnalysisManager::kOutputContainer,
      Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput  (rawJetTask, 0,  cinput1 );
  mgr->ConnectOutput (rawJetTask, 1, coutput1 );

  return rawJetTask;
}













