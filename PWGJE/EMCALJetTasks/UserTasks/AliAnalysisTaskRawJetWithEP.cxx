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

#include <TGrid.h>
#include <TChain.h>

#include <TClonesArray.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TCanvas.h>
#include <TString.h>

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

#include "AliTLorentzVector.h"
#include "AliEmcalJet.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"

#include "AliJEQnVectorHandler.h"

#include "AliLocalRhoParameter.h"
#include "AliRhoParameter.h"

#include "AliAnalysisTaskFlowVectorCorrectionsPWGJE.h"
#include "AliAnalysisTaskRawJetWithEP1.h"

#include "AliDataFile.h"

class AliAnalysisTaskRawJetWithEP1;

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskRawJetWithEP1);
/// \endcond

/**
 * Default constructor. Needed by ROOT I/O
 */
AliAnalysisTaskRawJetWithEP1::AliAnalysisTaskRawJetWithEP1() :
  AliAnalysisTaskEmcalJet(),
  fAOD(nullptr),
  fOutputList(nullptr),
  fEventCuts(),
  fOADBFile(nullptr),
  fCalibRefFileName(""),
  fCalibRefFile(nullptr),
  fCalibRefObjList(nullptr),
  fCalibV0Ref(nullptr),
  fQ2VecHandler(nullptr),
  fQ3VecHandler(nullptr),
  fV0Q2VectTask(0x0),
  fHistManager(),
  fBefCalibQ2_Cent0(0x0),
  fAftCalibQ2_Cent0(0x0),
  fHCorrV0ChWeghts(NULL),
  fHCorrQ2xV0C(NULL),
  fHCorrQ2yV0C(NULL),
  fHCorrQ2xV0A(NULL),
  fHCorrQ2yV0A(NULL),
  fHCorrQ3xV0C(NULL),
  fHCorrQ3yV0C(NULL),
  fHCorrQ3xV0A(NULL),
  fHCorrQ3yV0A(NULL), 
  fFitModulationType(kNoFit), fFitModulation(0), hBkgTracks(0),
  fV2ResoV0(0.), fV3ResoV0(0.), CheckRunNum(0),
  fV0Q2Vector(0.), fV0Ep2Angle(0.)
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
AliAnalysisTaskRawJetWithEP1::AliAnalysisTaskRawJetWithEP1(const char *name):
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fAOD(nullptr),
  fOutputList(nullptr),
  fEventCuts(),
  fOADBFile(nullptr),
  fCalibRefFileName(""),
  fCalibRefFile(nullptr),
  fCalibRefObjList(nullptr),
  fCalibV0Ref(nullptr),
  fQ2VecHandler(nullptr),
  fQ3VecHandler(nullptr),
  fV0Q2VectTask(0x0),
  fHistManager(name),
  fBefCalibQ2_Cent0(0x0),
  fAftCalibQ2_Cent0(0x0),
  fHCorrV0ChWeghts(NULL),
  fHCorrQ2xV0C(NULL),
  fHCorrQ2yV0C(NULL),
  fHCorrQ2xV0A(NULL),
  fHCorrQ2yV0A(NULL),
  fHCorrQ3xV0C(NULL),
  fHCorrQ3yV0C(NULL),
  fHCorrQ3xV0A(NULL),
  fHCorrQ3yV0A(NULL), 
  fFitModulationType(kNoFit),  fFitModulation(0), hBkgTracks(0),
  fV2ResoV0(0.), fV3ResoV0(0.), CheckRunNum(0),
  fV0Q2Vector(0.), fV0Ep2Angle(0.)
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
AliAnalysisTaskRawJetWithEP1::~AliAnalysisTaskRawJetWithEP1()
{
  
  if(fOADBFile) {fOADBFile->Close(); fOADBFile = 0x0;}
  if(fCalibRefFile) {fCalibRefFile->Close(); fCalibRefFile = 0x0;}
  if(fCalibV0Ref)   {delete fCalibV0Ref;   fCalibV0Ref = 0x0;}
  if(fCalibRefObjList) {delete fCalibRefObjList; fCalibRefObjList = 0x0;}
  if(fQ2VecHandler) {delete fQ2VecHandler; fQ2VecHandler = 0x0;}
  if(fQ3VecHandler) {delete fQ3VecHandler; fQ3VecHandler = 0x0;}
  if(fFitModulation) {delete fFitModulation; fFitModulation = 0x0;}
  if(hBkgTracks) {delete hBkgTracks; hBkgTracks = 0x0;}
}

/**
 * Performing run-independent initialization.
 * Here the histograms should be instantiated.
 */
void AliAnalysisTaskRawJetWithEP1::UserCreateOutputObjects()
{
  
  fOutputList = new TList();
  fOutputList->SetOwner(true);

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  fEventCuts.AddQAplotsToList(fOutput);
  fEventCuts.OverrideAutomaticTriggerSelection(AliVEvent::kAnyINT | AliVEvent::kCentral | AliVEvent::kSemiCentral);
  
  // == s == Set Out put Hist grams  ###########################################
  AllocateEventPlaneHistograms();
  AllocateTrackHistograms();
  AllocateJetHistograms();
  // == e == Set Out put Hist grams  ###########################################
  
  // == s == Add Objects into output file  #####################################
  TIter next(fHistManager.GetListOfHistograms());
  TObject* obj = 0;
  while ((obj = next())) {
    fOutput->Add(obj);
  }
  // == e == Add Objects into output file  #####################################
  

  // == s == Calib root file include  ==============================--------===-
  TString pathToFileCMVFNS = AliDataFile::GetFileName(fOADBFileName.Data());
  TString pathToFileLocal;
  std::cout << fOADBFileName << std::endl;
  pathToFileLocal = fOADBFileName;
  // Check access to CVMFS (will only be displayed locally)
  if(fOADBFileName.BeginsWith("alien://") && !gGrid){
    AliInfo("Trying to connect to AliEn ...");
    TGrid::Connect("alien://");
  }
  if(!pathToFileCMVFNS.IsNull()) fOADBFile = TFile::Open(pathToFileCMVFNS.Data());
  if(pathToFileCMVFNS.IsNull())  fOADBFile = TFile::Open(pathToFileLocal.Data());
  if(!fOADBFile) {
    AliWarning("OADB V0-TPC calibration file cannot be opened\n");
    return;
  }

  TString tempCalibFileName = AliDataFile::GetFileName(fCalibRefFileName.Data());
  TString tempCalibLocalFileName;
  std::cout << fCalibRefFileName << std::endl;
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

void AliAnalysisTaskRawJetWithEP1::AllocateEventPlaneHistograms()
{
  TString histName;
  TString histtitle;
  TString groupName;
  groupName="EventPlane";
  fHistManager.CreateHistoGroup(groupName);

  histName  = TString::Format("%s/hCentrality", groupName.Data());
  histtitle = TString::Format("%s;Centrality;counts", histName.Data());
  fHistManager.CreateTH1(histName, histtitle, 100, 0, 100);

  // == s == gain calibration QA ============================================
  histName  = TString::Format("%s/hV0CellChGains", groupName.Data());
  histtitle = TString::Format("%s;cell ch number;CorrGain", histName.Data());
  fHistManager.CreateTH1(histName, histtitle, 64, 0, 64);
  histName  = TString::Format("%s/hV0CellChRatio", groupName.Data());
  histtitle = TString::Format("%s;cell ch number;Ch ratio", histName.Data());
  fHistManager.CreateTH1(histName, histtitle, 64, 0, 64);

  histName  = TString::Format("%s/hV0CellChVsMultBefEq", groupName.Data());
  histtitle = TString::Format("%s;cell ch number;multiplisity before equalization", histName.Data());
  fHistManager.CreateTH1(histName, histtitle, 64, 0, 64);

  histName  = TString::Format("%s/hV0CellChVsMultAftEq", groupName.Data());
  histtitle = TString::Format("%s;cell ch number;multiplisity before equalization", histName.Data());
  fHistManager.CreateTH1(histName, histtitle, 64, 0, 64);
  // == e == gain calibration QA ============================================

  // == s == recentrering calibration QA =======================================
  histName  = TString::Format("%s/hCentVsQxV0ABefCalib", groupName.Data());
  histtitle = TString::Format("%s;centrality;Q_{x} V0A", histName.Data());
  fHistManager.CreateTH1(histName, histtitle, 10, 0, 10);
  histName  = TString::Format("%s/hCentVsQyV0ABefCalib", groupName.Data());
  histtitle = TString::Format("%s;centrality;Q_{y} V0A", histName.Data());
  fHistManager.CreateTH1(histName, histtitle, 10, 0, 10);

  histName  = TString::Format("%s/hCentVsQxV0AAftCalib", groupName.Data());
  histtitle = TString::Format("%s;centrality;Q_{x} V0A", histName.Data());
  fHistManager.CreateTH1(histName, histtitle, 10, 0, 10);
  histName  = TString::Format("%s/hCentVsQyV0AAftCalib", groupName.Data());
  histtitle = TString::Format("%s;centrality;Q_{y} V0A", histName.Data());
  fHistManager.CreateTH1(histName, histtitle, 10, 0, 10);

  histName  = TString::Format("%s/hCentVsQxV0CBefCalib", groupName.Data());
  histtitle = TString::Format("%s;centrality;Q_{x} V0C", histName.Data());
  fHistManager.CreateTH1(histName, histtitle, 10, 0, 10);
  histName  = TString::Format("%s/hCentVsQyV0CBefCalib", groupName.Data());
  histtitle = TString::Format("%s;centrality;Q_{y} V0C", histName.Data());
  fHistManager.CreateTH1(histName, histtitle, 10, 0, 10);

  histName  = TString::Format("%s/hCentVsQxV0CAftCalib", groupName.Data());
  histtitle = TString::Format("%s;centrality;Q_{x} V0C", histName.Data());
  fHistManager.CreateTH1(histName, histtitle, 10, 0, 10);
  histName  = TString::Format("%s/hCentVsQyV0CAftCalib", groupName.Data());
  histtitle = TString::Format("%s;centrality;Q_{y} V0C", histName.Data());
  fHistManager.CreateTH1(histName, histtitle, 10, 0, 10);
  // == e == recentrering calibration QA =======================================


  for (Int_t cent = 0; cent < fNcentBins; cent++) {
      // == s == Event plane angle histograms Setting
      histName = TString::Format("%s/hPsi2V0M_%d", groupName.Data(), cent);
      histtitle = "Psi2 from VZEROCA";
      //fHistManager.CreateTH1(histName, histtitle, fNbins / 2, 0, TMath::TwoPi());
      fHistManager.CreateTH1(histName, histtitle, fNbins / 2, -TMath::TwoPi(), 2*TMath::TwoPi());
      histName = TString::Format("%s/hPsi3V0M_%d", groupName.Data(), cent);
      histtitle = "Psi3 from V0M";
      //fHistManager.CreateTH1(histName, histtitle, fNbins / 2, 0, TMath::TwoPi());
      fHistManager.CreateTH1(histName, histtitle, fNbins / 2, -TMath::TwoPi(), 2*TMath::TwoPi());

      histName = TString::Format("%s/hPsi2V0C_%d", groupName.Data(), cent);
      histtitle = "Psi3 from V0A";
      //fHistManager.CreateTH1(histName, histtitle, fNbins / 2, 0, TMath::TwoPi());
      fHistManager.CreateTH1(histName, histtitle, fNbins / 2, -TMath::TwoPi(), 2*TMath::TwoPi());
      histName = TString::Format("%s/hPsi3V0C_%d", groupName.Data(), cent);
      histtitle = "Psi2 from V0C";
      //fHistManager.CreateTH1(histName, histtitle, fNbins / 2, 0, TMath::TwoPi());
      fHistManager.CreateTH1(histName, histtitle, fNbins / 2, -TMath::TwoPi(), 2*TMath::TwoPi());

      histName = TString::Format("%s/hPsi2V0A_%d", groupName.Data(), cent);
      histtitle = "Psi3 from V0A";
      //fHistManager.CreateTH1(histName, histtitle, fNbins / 2, 0, TMath::TwoPi());
      fHistManager.CreateTH1(histName, histtitle, fNbins / 2, -TMath::TwoPi(), 2*TMath::TwoPi());

      histName = TString::Format("%s/hPsi3V0A_%d", groupName.Data(), cent);
      histtitle = "Psi3 from V0A";
      //fHistManager.CreateTH1(histName, histtitle, fNbins / 2, 0, TMath::TwoPi());
      fHistManager.CreateTH1(histName, histtitle, fNbins / 2, -TMath::TwoPi(), 2*TMath::TwoPi());

      histName = TString::Format("%s/hPsi2TPCN_%d", groupName.Data(), cent);
      histtitle = "Psi2 from TPC eta negative";
      //fHistManager.CreateTH1(histName, histtitle, fNbins / 2, 0, TMath::TwoPi());
      fHistManager.CreateTH1(histName, histtitle, fNbins / 2, -TMath::TwoPi(), 2*TMath::TwoPi());
      histName = TString::Format("%s/hPsi3TPCN_%d", groupName.Data(), cent);
      histtitle = "Psi3 from TPC eta negative";
      //fHistManager.CreateTH1(histName, histtitle, fNbins / 2, 0, TMath::TwoPi());
      fHistManager.CreateTH1(histName, histtitle, fNbins / 2, -TMath::TwoPi(), 2*TMath::TwoPi());

      histName = TString::Format("%s/hPsi2TPCP_%d", groupName.Data(), cent);
      histtitle = "Psi2 from TPC eta positive";
      //fHistManager.CreateTH1(histName, histtitle, fNbins / 2, 0, TMath::TwoPi());
      fHistManager.CreateTH1(histName, histtitle, fNbins / 2, -TMath::TwoPi(), 2*TMath::TwoPi());
      histName = TString::Format("%s/hPsi3TPCP_%d", groupName.Data(), cent);
      histtitle = "Psi3 from TPC eta positive";
      //fHistManager.CreateTH1(histName, histtitle, fNbins / 2, 0, TMath::TwoPi());
      fHistManager.CreateTH1(histName, histtitle, fNbins / 2, -TMath::TwoPi(), 2*TMath::TwoPi());
  }

  // cdf and pdf of chisquare distribution
  histName = TString::Format("%s/hPvalueCDF", groupName.Data());
  histtitle = TString::Format("%s; CDF #chi^{2}", histName.Data());
  fHistManager.CreateTH1(histName, histtitle, 50, 0, 1);
  
  histName = TString::Format("%s/hPvalueCDFCent", groupName.Data());
  histtitle = TString::Format("%s; centrality; p-value", histName.Data());
  fHistManager.CreateTH2(histName, histtitle, 40, 0, 100, 40, 0, 1);
  
  histName = TString::Format("%s/hChi2Cent", groupName.Data());
  histtitle = TString::Format("%s; CDF #chi^{2}; #tilde{#chi^{2}}", histName.Data());
  fHistManager.CreateTH2(histName, histtitle, 100, 0, 100, 100, 0, 5);

  histName = TString::Format("%s/hPChi2", groupName.Data());
  histtitle = TString::Format("%s; p-value; #tilde{#chi^{2}}", histName.Data());
  fHistManager.CreateTH2(histName, histtitle, 1000, 0, 1, 100, 0, 5);

  histName = TString::Format("%s/hKolmogorovTest", groupName.Data());
  histtitle = TString::Format("%s;CDF #chi^{2}", histName.Data());
  fHistManager.CreateTH1(histName, histtitle, 50, 0, 1);

  histName = TString::Format("%s/hKolmogorovTestCent", groupName.Data());
  histtitle = TString::Format("%s; centrality; Kolmogorov p", histName.Data());
  fHistManager.CreateTH2(histName, histtitle, 100, 0, 100, 45, 0, 1);

  histName = TString::Format("%s/hPvalueCDFROOT", groupName.Data());
  histtitle = TString::Format("%s; CDF #chi^{2} ROOT", histName.Data());
  fHistManager.CreateTH1(histName, histtitle, 50, 0, 1);

  histName = TString::Format("%s/hPvalueCDFROOTCent", groupName.Data());
  histtitle = TString::Format("%s; centrality; p-value ROOT", histName.Data());
  fHistManager.CreateTH2(histName, histtitle, 40, 0, 100, 45, 0, 1);

  histName = TString::Format("%s/hChi2ROOTCent", groupName.Data());
  histtitle = TString::Format("%s; p-value; #tilde{#chi^{2}}", histName.Data());
  fHistManager.CreateTH2(histName, histtitle, 100, 0, 100, 45, 0, 1);

  histName = TString::Format("%s/hPChi2Root", groupName.Data());
  histtitle = TString::Format("%s;CDF #chi^{2}; #tilde{#chi^{2}} ROOT", histName.Data());
  fHistManager.CreateTH2(histName, histtitle, 1000, 0, 1, 100, 0, 5);

  histName = TString::Format("%s/hPKolmogorov", groupName.Data());
  histtitle = TString::Format("%s; p-value; kolmogorov p", histName.Data());
  fHistManager.CreateTH2(histName, histtitle, 40, 0, 1, 40, 0, 1);
  
  
  histName = TString::Format("%s/hRhoVsMult", groupName.Data());
  histtitle = TString::Format("%s; multiplicity; #rho [GeV/c]", histName.Data());
  fHistManager.CreateTH2(histName, histtitle, 100, 0, 4000, 100, 0, 250);
  
  histName = TString::Format("%s/hRhoVsCent", groupName.Data());
  histtitle = TString::Format("%s; centrality; #rho [GeV/c]", histName.Data());
  fHistManager.CreateTH2(histName, histtitle, 100, 0, 100, 100, 0, 250);
  
  histName = TString::Format("%s/hRhoAVsMult", groupName.Data());
  histtitle = TString::Format("%s; multiplicity; #rho * A (jet) [GeV/c]", histName.Data());
  fHistManager.CreateTH2(histName, histtitle, 100, 0, 4000, 100, 0, 50);
  
  histName = TString::Format("%s/hRhoAVsCent", groupName.Data());
  histtitle = TString::Format("%s; centrality; #rho * A (jet) [GeV/c]", histName.Data());
  fHistManager.CreateTH2(histName, histtitle, 100, 0, 100, 100, 0, 50);

}

/*
 * This function allocates the histograms for basic tracking QA.
 * A set of histograms (pT, eta, phi, difference between kinematic properties
 * at the vertex and at the EMCal surface, number of tracks) is allocated
 * per each particle container and per each centrality bin.
 */
void AliAnalysisTaskRawJetWithEP1::AllocateTrackHistograms()
{
  TString histName;
  TString histtitle;
  TString groupName;
  AliParticleContainer* partCont = 0;
  TIter next(&fParticleCollArray);
  while ((partCont = static_cast<AliParticleContainer*>(next()))) {
    groupName = partCont->GetName();
    // Protect against creating the histograms twice
    if (fHistManager.FindObject(groupName)) {
      AliWarning(TString::Format("%s: Found groupName %s in hist manager. The track containers will be filled into the same histograms.", GetName(), groupName.Data()));
      continue;
    }
    fHistManager.CreateHistoGroup(groupName);

    // adding histo for counting events
    histName = TString::Format("Hist nEvents");
    histtitle = TString::Format("Number of Events");
    fHistManager.CreateTH1(histName, histtitle, 1, 0.0, 1.0);

    for (Int_t cent = 0; cent < fNcentBins; cent++) {
      histName = TString::Format("%s/hTrackPt_%d", groupName.Data(), cent);
      histtitle = TString::Format("%s;#it{p}_{T,track} (GeV/#it{c});counts", histName.Data());
      fHistManager.CreateTH1(histName, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt / 2);

      histName = TString::Format("%s/hTrackPhi_%d", groupName.Data(), cent);
      histtitle = TString::Format("%s;#it{#phi}_{track};counts", histName.Data());
      fHistManager.CreateTH1(histName, histtitle, fNbins / 2, 0, TMath::TwoPi());

      histName = TString::Format("%s/hTrackEta_%d", groupName.Data(), cent);
      histtitle = TString::Format("%s;#it{#eta}_{track};counts", histName.Data());
      fHistManager.CreateTH1(histName, histtitle, fNbins / 6, -1, 1);

      histName = TString::Format("%s/hNTracks_%d", groupName.Data(), cent);
      histtitle = TString::Format("%s;number of tracks;events", histName.Data());
      if (fForceBeamType != kpp) fHistManager.CreateTH1(histName, histtitle, 500, 0, 5000);
      else fHistManager.CreateTH1(histName, histtitle, 200, 0, 200);
      



      // == e == Event plane angle histograms Setting
    }
  }

  histName = "fHistSumNTracks";
  histtitle = TString::Format("%s;Sum of n tracks;events", histName.Data());
  if (fForceBeamType != kpp) fHistManager.CreateTH1(histName, histtitle, 500, 0, 5000);
  else fHistManager.CreateTH1(histName, histtitle, 200, 0, 200);
}

/*
 * This function allocates the histograms for basic jet QA.
 * A set of histograms (pT, eta, phi, area, number of jets, corrected pT) is allocated
 * per each jet container and per each centrality bin.
 */
void AliAnalysisTaskRawJetWithEP1::AllocateJetHistograms()
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

    // A vs. pT
    if (fForceBeamType == kAA) {
      histName = TString::Format("%s/hRhoVsCent", groupName.Data());
      histtitle = histName + ";Centrality (%);#rho (GeV/#it{c});counts";
    	fHistManager.CreateTH2(histName, histtitle.Data(), 50, 0, 100, 100, 0, 500);

      histName = TString::Format("%s/CentVsPtVsPtdence", groupName.Data());
      histtitle = histName + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});#it{A}_{jet}/#pi#it{R}^{2}";
      // fHistManager.CreateTH3(histName, histtitle.Data(), 10, 0, 100, 300, 0, 250, 75, 0, 3);
    }

    for (Int_t cent = 0; cent < fNcentBins; cent++) {
      histName = TString::Format("%s/hJetPt_%d", groupName.Data(), cent);
      histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c});counts", histName.Data());
      fHistManager.CreateTH1(histName, histtitle, fNbins, fMinBinPt, fMaxBinPt);

      histName = TString::Format("%s/hJetArea_%d", groupName.Data(), cent);
      histtitle = TString::Format("%s;#it{A}_{jet};counts", histName.Data());
      fHistManager.CreateTH1(histName, histtitle, fNbins / 2, 0, 3);

      histName = TString::Format("%s/hJetPhi_%d", groupName.Data(), cent);
      histtitle = TString::Format("%s;#it{#phi}_{jet};counts", histName.Data());
      fHistManager.CreateTH1(histName, histtitle, fNbins / 2, 0, TMath::TwoPi());

      histName = TString::Format("%s/hJetEta_%d", groupName.Data(), cent);
      histtitle = TString::Format("%s;#it{#eta}_{jet};counts", histName.Data());
      fHistManager.CreateTH1(histName, histtitle, fNbins / 6, -1, 1);

      histName = TString::Format("%s/hNJets_%d", groupName.Data(), cent);
      histtitle = TString::Format("%s;number of jets;events", histName.Data());
      if (fForceBeamType != kpp) {
        fHistManager.CreateTH1(histName, histtitle, 500, 0, 500);
      }
      else {
        fHistManager.CreateTH1(histName, histtitle, 100, 0, 100);
      }

      // histograms for jet angle relative to the event plane
      histName = TString::Format("%s/hJetPhiMinusPsi2_%d", groupName.Data(), cent);
      histtitle = "Jet phi minus psi2";
      fHistManager.CreateTH1(histName, histtitle, fNbins / 2, 0, TMath::TwoPi());
      histName = TString::Format("%s/hJetPhiMinusPsi3_%d", groupName.Data(), cent);
      histtitle = "Jet phi minus psi3";
      fHistManager.CreateTH1(histName, histtitle, fNbins / 2, 0, TMath::TwoPi());

      // histograms for in-plane vs out-of-plane jets
      histName = TString::Format("%s/hJetsInPlaneOutOfPlaneV2_%d", groupName.Data(), cent);
      histtitle = "In-plane vs out-of-plane jets (v2)";
      fHistManager.CreateTH2(histName, histtitle, 2, 0, 2, 10, 0, 100);
      histName = TString::Format("%s/hJetsInPlaneOutOfPlaneV3_%d", groupName.Data(), cent);
      histtitle = "In-plane vs out-of-plane jets (v2)";
      fHistManager.CreateTH2(histName, histtitle, 2, 0, 2, 10, 0, 100);


      if (!jetCont->GetRhoName().IsNull()) {
        
        // Rho histograms
        histName = TString::Format("%s/hJetRho_%d", groupName.Data(), cent);
        
        histtitle = "Rho";
        fHistManager.CreateTH1(histName, histtitle, fNbins, 0.0, 300.0);
        histName = TString::Format("%s/hJetRhoLocal_%d", groupName.Data(), cent);
        histtitle = "Rho Local";
        fHistManager.CreateTH1(histName, histtitle, fNbins, 0.0, 300.0);
        histName = TString::Format("%s/hJetCorrPt_%d", groupName.Data(), cent);
        histtitle = TString::Format("%s;#it{p}_{T,jet}^{corr} (GeV/#it{c});counts", histName.Data());
        fHistManager.CreateTH1(histName, histtitle, fNbins, -fMaxBinPt / 2, fMaxBinPt / 2);
        histName = TString::Format("%s/hJetCorrPtLocal_%d", groupName.Data(), cent);
        histtitle = TString::Format("%s;#it{p}_{T,jet}^{corr} local (GeV/#it{c});counts", histName.Data());
        fHistManager.CreateTH1(histName, histtitle, fNbins, -fMaxBinPt / 2, fMaxBinPt / 2);
        
        // = s = Create histograms of Jet Yeild ====================================================
        //v2 inlucive
        histName = TString::Format("%s/hJetsYieldV2_%d", groupName.Data(), cent);
        histtitle = "Jet yeild (v2)";
        fHistManager.CreateTH1(histName, histtitle, 300, -50, 250);
        histName = TString::Format("%s/hCorrJetsYieldV2_%d", groupName.Data(), cent);
        histtitle = "corr Jet yeild (v2)";
        fHistManager.CreateTH1(histName, histtitle, 300, -50, 250);
        histName = TString::Format("%s/hLocalCorrJetsYieldV2_%d", groupName.Data(), cent);
        histtitle = "local corr Jet yeild (v2)";
        fHistManager.CreateTH1(histName, histtitle, 300, -50, 250);

        //v2 in plane
        histName = TString::Format("%s/hJetsYieldInPlaneV2_%d", groupName.Data(), cent);
        histtitle = "Jet yeild of in-plane (v2)";
        fHistManager.CreateTH1(histName, histtitle, 300, -50, 250);
        histName = TString::Format("%s/hCorrJetsYieldInPlaneV2_%d", groupName.Data(), cent);
        histtitle = "corr Jet yeild of in-plane (v2)";
        fHistManager.CreateTH1(histName, histtitle, 300, -50, 250);
        histName = TString::Format("%s/hLocalCorrJetsYieldInPlaneV2_%d", groupName.Data(), cent);
        histtitle = "local corr Jet yeild of in-plane (v2)";
        fHistManager.CreateTH1(histName, histtitle, 300, -50, 250);

        //v2 out of plane
        histName = TString::Format("%s/hJetsYieldOutPlaneV2_%d", groupName.Data(), cent);
        histtitle = "Jet yeild of out-plane (v2)";
        fHistManager.CreateTH1(histName, histtitle, 300, -50, 250);
        histName = TString::Format("%s/hCorrJetsYieldOutPlaneV2_%d", groupName.Data(), cent);
        histtitle = "corr Jet yeild of in-plane (v2)";
        fHistManager.CreateTH1(histName, histtitle, 300, -50, 250);
        histName = TString::Format("%s/hLocalCorrJetsYieldOutPlaneV2_%d", groupName.Data(), cent);
        histtitle = "local corr Jet yeild of out-plane (v2)";
        fHistManager.CreateTH1(histName, histtitle, 300, -50, 250);


        //v3 inlucive
        histName = TString::Format("%s/hJetsYieldV3_%d", groupName.Data(), cent);
        histtitle = "Jet yeild (v3)";
        fHistManager.CreateTH1(histName, histtitle, 300, -50, 250);
        histName = TString::Format("%s/hCorrJetsYieldV3_%d", groupName.Data(), cent);
        histtitle = "corr Jet yeild (v3)";
        fHistManager.CreateTH1(histName, histtitle, 300, -50, 250);
        histName = TString::Format("%s/hLocalCorrJetsYieldV3_%d", groupName.Data(), cent);
        histtitle = "local corr Jet yeild (v3)";
        fHistManager.CreateTH1(histName, histtitle, 300, -50, 250);

        //v3 in plane
        histName = TString::Format("%s/hJetsYieldInPlaneV3_%d", groupName.Data(), cent);
        histtitle = "Jet yeild of in-plane (v3)";
        fHistManager.CreateTH1(histName, histtitle, 300, -50, 250);
        histName = TString::Format("%s/hCorrJetsYieldInPlaneV3_%d", groupName.Data(), cent);
        histtitle = "corr Jet yeild of in-plane (v3)";
        fHistManager.CreateTH1(histName, histtitle, 300, -50, 250);
        histName = TString::Format("%s/hLocalCorrJetsYieldInPlaneV3_%d", groupName.Data(), cent);
        histtitle = "local corr Jet yeild of in-plane (v3)";
        fHistManager.CreateTH1(histName, histtitle, 300, -50, 250);

        //v3 out of plane
        histName = TString::Format("%s/hJetsYieldOutPlaneV3_%d", groupName.Data(), cent);
        histtitle = "Jet yeild of out-plane (v3)";
        fHistManager.CreateTH1(histName, histtitle, 300, -50, 250);
        histName = TString::Format("%s/hCorrJetsYieldOutPlaneV3_%d", groupName.Data(), cent);
        histtitle = "corr Jet yeild of in-plane (v3)";
        fHistManager.CreateTH1(histName, histtitle, 300, -50, 250);
        histName = TString::Format("%s/hLocalCorrJetsYieldOutPlaneV3_%d", groupName.Data(), cent);
        histtitle = "local corr Jet yeild of out-plane (v3)";
        fHistManager.CreateTH1(histName, histtitle, 300, -50, 250);
        // = e = Create histograms of Jet Yeild ====================================================
        


        histName = TString::Format("%s/hJetLocalRhoVsAverageRho_%d", groupName.Data(), cent);
        histtitle = "Local rho versus average rho";
        fHistManager.CreateTH2(histName, histtitle, fNbins, 0.0, 300.0, fNbins, 0.0, 300.0);
        histName = TString::Format("%s/hJetCorrPtLocalVsJetCorrPt_%d", groupName.Data(), cent);
        histtitle = "Local rho adjusted jet pT versus average rho adjusted jet pT";
        fHistManager.CreateTH2(histName, histtitle, fNbins, 0.0, 100.0, fNbins, 0.0, 100.0);

        // histo local rho vs delta phi
        histName = TString::Format("%s/hJetRhoVsDeltaPhi_%d", groupName.Data(), cent);
        histtitle = "Local rho versus angle relative to event plane";
        fHistManager.CreateTH2(histName, histtitle, fNbins, 0.0, TMath::TwoPi(), fNbins, 0.0, 300.0);
        
      }
    }
  }
}


/**
 * This function is executed automatically for the first event.
 * Some extra initialization can be performed here.
 */
void AliAnalysisTaskRawJetWithEP1::ExecOnce()
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
Bool_t AliAnalysisTaskRawJetWithEP1::Run()
{
  std::cout << "ChecKuma Run Number === " << CheckRunNum << "==================" << std::endl;
  CheckRunNum++;
  if(!fEventCuts.AcceptEvent(InputEvent())) return kFALSE;
  
  // if(fPreRunNum != fRunNumber){
  //   fPreRunNum = fRunNumber;
  //   fRunOrder += 1;
  //   fRunNumList.push_back(fRunNumber);
  // }
  std::cout << "fCentBin = " << fCentBin << ", fCent = " << fCent << "  CheKumaaaaa" << std::endl;

  DoEventPlane();
  SetModulationRhoFit();
  // std::cout << "Fomula = " << fFitModulation->GetExpFormula() << std::endl;
  // MeasureTpcEPQA();
  MeasureBkg();
  DoJetLoop();
  
  return kTRUE;
}

void AliAnalysisTaskRawJetWithEP1::DoEventPlane(){

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

  QnGainCalibration();
  QnRecenteringCalibration();

  QnJEHandlarEPGet();
  // if(CheckRunNum == fQaEventNum) VzeroGainCalibQA();
  
  TString histName;
  TString groupName;
  groupName="EventPlane";

  histName = TString::Format("%s/hCentrality", groupName.Data());
  fHistManager.FillTH1(histName, fCent);

  histName = TString::Format("%s/hPsi2V0M_%d", groupName.Data(), fCentBin);
  fHistManager.FillTH1(histName, psi2V0[0]);
  histName = TString::Format("%s/hPsi3V0M_%d", groupName.Data(), fCentBin);
  fHistManager.FillTH1(histName, psi3V0[0]);
  
  histName = TString::Format("%s/hPsi2V0C_%d", groupName.Data(), fCentBin);
  fHistManager.FillTH1(histName, psi2V0[1]);
  histName = TString::Format("%s/hPsi3V0C_%d", groupName.Data(), fCentBin);
  fHistManager.FillTH1(histName, psi3V0[1]);

  histName = TString::Format("%s/hPsi2V0A_%d", groupName.Data(), fCentBin);
  fHistManager.FillTH1(histName, psi2V0[2]);
  histName = TString::Format("%s/hPsi3V0A_%d", groupName.Data(), fCentBin);
  fHistManager.FillTH1(histName, psi3V0[2]);

  histName = TString::Format("%s/hPsi2TPCN_%d", groupName.Data(), fCentBin);
  fHistManager.FillTH1(histName, psi2Tpc[2]);
  histName = TString::Format("%s/hPsi3TPCN_%d", groupName.Data(), fCentBin);
  fHistManager.FillTH1(histName, psi3Tpc[2]);
  
  histName = TString::Format("%s/hPsi2TPCP_%d", groupName.Data(), fCentBin);
  fHistManager.FillTH1(histName, psi2Tpc[1]);
  histName = TString::Format("%s/hPsi2TPCP_%d", groupName.Data(), fCentBin);
  fHistManager.FillTH1(histName, psi3Tpc[1]);
  
}


void AliAnalysisTaskRawJetWithEP1::SetModulationRhoFit() 
{
  // set modulation fit
  TString histName;

  if(fFitModulation) delete fFitModulation;
  fFitModulation = 0x0;
  const char * fitFunction = "[0]*([1]+[2]*([3]*TMath::Cos([2]*(x-[4]))+[7]*TMath::Cos([5]*(x-[6]))))";
  switch (fFitModulationType)  {
    case kNoFit : { fFitModulation = new TF1("fix_kNoFit", "[0]", 0, TMath::TwoPi()); } break;
    case kV2 : {
      fitFunction = "[0]*([1]+[2]*[3]*TMath::Cos([2]*(x-[4])))";
      fFitModulation = new TF1("fit_kV2", fitFunction, 0, TMath::TwoPi());
      fFitModulation->SetParameter(0, 0.);        // normalization
      fFitModulation->SetParameter(3, 0.2);       // v2
      fFitModulation->FixParameter(1, 1.);        // constant
      fFitModulation->FixParameter(2, 2.);        // constant
    } break;
    case kCombined: {
      fitFunction = "[0]*([1]+[2]*([3]*TMath::Cos([2]*(x-[4]))+[7]*TMath::Cos([5]*(x-[6]))))";
      fFitModulation = new TF1("fit_kCombined", fitFunction, 0, TMath::TwoPi());
      fFitModulation->SetParameter(0, 0.);       // normalization
      fFitModulation->SetParameter(3, 0.2);      // v2
      fFitModulation->FixParameter(1, 1.);       // constant
      fFitModulation->FixParameter(2, 2.);       // constant
      fFitModulation->FixParameter(5, 3.);       // constant
      fFitModulation->SetParameter(7, 0.2);      // v3
    } break;
    default : { // for the combined fit, the 'direct fourier series' or the user supplied vn values we use v2 and v3
      fitFunction = "[0]*([1]+[2]*([3]*TMath::Cos([2]*(x-[4]))+[7]*TMath::Cos([5]*(x-[6]))))";
      fFitModulation = new TF1("fit_kCombined", fitFunction, 0, TMath::TwoPi());
      fFitModulation->SetParameter(0, 0.);       // normalization
      fFitModulation->SetParameter(3, 0.2);      // v2
      fFitModulation->FixParameter(1, 1.);       // constant
      fFitModulation->FixParameter(2, 2.);       // constant
      fFitModulation->FixParameter(5, 3.);       // constant
      fFitModulation->SetParameter(7, 0.2);      // v3
    } break;
  }

  

  if(hBkgTracks) delete hBkgTracks;
  histName = "hBkgTracks";
  // hBkgTracks = new TH1F(histName, histName, 100, 0.0, TMath::TwoPi());
  hBkgTracks = new TH1F(histName, histName, 25, 0.0, TMath::TwoPi());
}


void AliAnalysisTaskRawJetWithEP1::MeasureBkg(){
  TString histName;
  
  UInt_t sumAcceptedTracks = 0;
  AliParticleContainer* partCont = 0;
  TIter next(&fParticleCollArray);
  TH1F *_tempkBkgTracks = new TH1F("_tempSwap", "_tempSwap", 30, 0., TMath::TwoPi());
  while ((partCont = static_cast<AliParticleContainer*>(next()))) {
    // groupname = partCont->GetName();
    for(auto part : partCont->accepted()) {
      if (!part) continue;
      if (partCont->GetLoadedClass()->InheritsFrom("AliVTrack")) {
        const AliVTrack* track = static_cast<const AliVTrack*>(part);

        Float_t trackDeltaPhi = track->Phi() - psi2V0[0];
        if (trackDeltaPhi < 0.0) trackDeltaPhi += TMath::TwoPi();
        hBkgTracks->Fill(trackDeltaPhi);
        _tempkBkgTracks->Fill(trackDeltaPhi);
      }
    }
  }
  
  fLocalRho->SetVal(fRho->GetVal());
  
  fFitModulation->SetParameter(0, fLocalRho->GetVal()); //Fix ChecKuma????
  fFitModulation->FixParameter(4, psi2V0[0]);
  fFitModulation->FixParameter(6, psi3V0[0]);
  
  hBkgTracks->Fit(fFitModulation, "N0Q"); //??????????????????????????????????????
  

  if(1){ // ChecKuma fit parameters
    std::cout << "psi2V0 = " << psi2V0[0] << ", psi3V0 = " << psi3V0[0] << std::endl;
    
    if(CheckRunNum == 24) hBkgTracks->SaveAs("checkOutput/checkhBkgTracks.root");
    std::cout << "  p0 = " << fFitModulation->GetParameter(0)\
              << ", p1 = " << fFitModulation->GetParameter(1)\
              << ", p2 = " << fFitModulation->GetParameter(2)\
              << ", p3 = " << fFitModulation->GetParameter(3)\
              << ", p4 = " << fFitModulation->GetParameter(4)\
              << ", p5 = " << fFitModulation->GetParameter(5)\
              << ", p6 = " << fFitModulation->GetParameter(6)\
              << ", p7 = " << fFitModulation->GetParameter(7)\
              << std::endl;
    
  }


  // if(CheckRunNum == 5){
  if(0){
    TCanvas *cBkgRhoFit = new TCanvas("cBkgRhoFit", "cBkgRhoFit", 2000, 1500);
    
    TH1F* hBkgTracks_Event = (TH1F*) hBkgTracks->Clone("hnew");
    histName = hBkgTracks->GetName() + std::to_string(CheckRunNum);
    hBkgTracks_Event->SetName(histName);
    hBkgTracks_Event->Draw("E");
    hBkgTracks->Fit(fFitModulation, "Q");
    fFitModulation->SetLineColor(632);
    fFitModulation->Draw("same");

    TF1* rhoFitV2Com = new TF1("rhoFitV2Com", "[0]*([1]+[2]*([3]*TMath::Cos([2]*(x-[4]))))", 0.0, TMath::TwoPi());
    rhoFitV2Com->SetParameter(0, fFitModulation->GetParameter(0));
    rhoFitV2Com->SetParameter(1, fFitModulation->GetParameter(1));
    rhoFitV2Com->SetParameter(2, fFitModulation->GetParameter(2));
    rhoFitV2Com->SetParameter(3, fFitModulation->GetParameter(3));
    rhoFitV2Com->SetParameter(4, fFitModulation->GetParameter(4));
    rhoFitV2Com->SetLineColor(808);
    rhoFitV2Com->Draw("same");

    TF1* rhoFitV3Com = new TF1("rhoFitV3Com", "[0]*([1]+[2]*([5]*TMath::Cos([3]*(x-[4]))))", 0.0, TMath::TwoPi());
    rhoFitV3Com->SetParameter(0, fFitModulation->GetParameter(0));
    rhoFitV3Com->SetParameter(1, fFitModulation->GetParameter(1));
    rhoFitV3Com->SetParameter(3, fFitModulation->GetParameter(5));
    rhoFitV3Com->SetParameter(4, fFitModulation->GetParameter(6));
    rhoFitV3Com->SetParameter(5, fFitModulation->GetParameter(7));
    rhoFitV3Com->SetLineColor(824);
    rhoFitV3Com->Draw("same");

    histName = "checkOutput/cBkgRhoFit_Cent" + std::to_string(fCentBin) +".root";
    cBkgRhoFit->SaveAs(histName);
    
    // histName = "checkOutput/hBkgTracks_Event" + std::to_string(CheckRunNum) +".root";
    hBkgTracks_Event->SaveAs(histName);
    delete cBkgRhoFit;
    delete hBkgTracks_Event;
    !chcekDrawBkg++;
  }
  
  fV2ResoV0 = CalcEPReso(2, psi2V0[0], psi2Tpc[1], psi2Tpc[2]);
  fV3ResoV0 = CalcEPReso(3, psi3V0[0], psi3Tpc[1], psi3Tpc[2]);
  std::cout << "v2Reso = " << fV2ResoV0 << ", v3Reso = " << fV3ResoV0 << std::endl;

  
  Double_t v2ObjV0 = -999., v3ObjV0= -999.;
  v2ObjV0 = fFitModulation->GetParameter(3);
  v3ObjV0 = fFitModulation->GetParameter(7);
  
  std::cout << "v2obj = " << v2ObjV0 << ", v3Obj = " << v3ObjV0 << std::endl;
  fLocalRho->SetLocalRho(fFitModulation);
  

  // fLocalRho->SetVal(fRho->GetVal());
  // std::cout << "fRho = " << fRho->GetVal() << std::endl;

  // the quality of the fit is evaluated from 1 - the cdf of the chi square distribution
  // three methods are available, all with their drawbacks. all are stored, one is selected to do the cut
  // Int_t numOfFreePara = 2; //v2, v3
  // Int_t NDF(fFitModulation->GetXaxis()->GetNbins()-numOfFreePara);
  // std::cout << "nbis: " << fFitModulation->GetXaxis()->GetNbins() << std::endl;
  // if(NDF == 0 || (float)NDF <= 0.) return;
  
  // Double_t CDF(1.-ChiSquareCDF(NDF, ChiSquare(fFitModulation, fFitModulation)));
  // Double_t CDFROOT(1.-ChiSquareCDF(NDF, fFitModulation->GetChisquare()));
  
}


Double_t AliAnalysisTaskRawJetWithEP1::CalcEPReso(Int_t n, \
  Double_t &psiA, Double_t &psiB, Double_t &psiC){
  
  Double_t vnReso = -999.;
  vnReso = TMath::Sqrt((TMath::Abs(TMath::Cos(n*(psiA - psiB))) \
                          * TMath::Abs(TMath::Cos(n*(psiA - psiC)))) \
                        / TMath::Abs(TMath::Cos(n*(psiB - psiC))));

  return vnReso;
}


void AliAnalysisTaskRawJetWithEP1::MeasureTpcEPQA(){
  TString histName;
  TString groupName;

  Double_t EtaAcc = 0.9;
  UInt_t sumAcceptedTracks = 0;
  AliParticleContainer* partCont = 0;
  
  TIter next(&fParticleCollArray);
  
  while ((partCont = static_cast<AliParticleContainer*>(next()))) {
    
    groupName = partCont->GetName();
    std::cout << "partCont: " << partCont << std::endl; //checkuma
    // counting events
    histName = TString::Format("Hist nEvents");
    fHistManager.FillTH1(histName, 0.5);
    
    UInt_t count = 0;
    for(auto part : partCont->accepted()) {
      if (!part) continue;
      count++;

      histName = TString::Format("%s/hTrackPt_%d", groupName.Data(), fCentBin);
      fHistManager.FillTH1(histName, part->Pt());

      histName = TString::Format("%s/hTrackPhi_%d", groupName.Data(), fCentBin);
      fHistManager.FillTH1(histName, part->Phi());

      histName = TString::Format("%s/hTrackEta_%d", groupName.Data(), fCentBin);
      fHistManager.FillTH1(histName, part->Eta());

      if (partCont->GetLoadedClass()->InheritsFrom("AliVTrack")) {
        const AliVTrack* track = static_cast<const AliVTrack*>(part);

        // Filling histos for angle relative to event plane
        Double_t phiMinusPsi2 = track->Phi() - psi2V0[0];
        Double_t phiMinusPsi3 = track->Phi() - psi3V0[0];
        if (phiMinusPsi2 < 0.0) phiMinusPsi2 += TMath::TwoPi();
        if (phiMinusPsi3 < 0.0) phiMinusPsi3 += TMath::TwoPi();
        // histName = TString::Format("%s/hTrackPhiMinusPsi2_%d", groupName.Data(), fCentBin);
        // fHistManager.FillTH1(histName, phiMinusPsi2);
        // histName = TString::Format("%s/hTrackPhiMinusPsi3_%d", groupName.Data(), fCentBin);
        // fHistManager.FillTH1(histName, phiMinusPsi3);
      }
    }
    sumAcceptedTracks += count;

    histName = TString::Format("%s/hNTracks_%d", groupName.Data(), fCentBin);
    fHistManager.FillTH1(histName, count);
    
    
  }

  histName = "fHistSumNTracks";
  fHistManager.FillTH1(histName, sumAcceptedTracks);
  
}


//________________________________________________________________________
void AliAnalysisTaskRawJetWithEP1::LoadSpliForqnPerce()
{
    // load splines for qn percentiles
    std::cout << "splinesFilePath: " << fSplinesFileName << std::endl;

    TString listnameTPC[3] = {"SplineListq2TPC", "SplineListq2TPCPosEta", "SplineListq2TPCNegEta"};
    TString listnameV0[3] = {"SplineListq2V0", "SplineListq2V0A", "SplineListq2V0C"};

    TString pathToFileCMVFNS = AliDataFile::GetFileName(fSplinesFileName.Data());
    std::cout << "pathToFileCMVFNS: " << pathToFileCMVFNS << std::endl;
    // Check access to CVMFS (will only be displayed locally)
    if (pathToFileCMVFNS.IsNull())
        AliFatal("Cannot access data files from CVMFS: please export ALICE_DATA=root://eospublic.cern.ch//eos/experiment/alice/analysis-data and run again");
    TFile* splinesfile = TFile::Open(pathToFileCMVFNS.Data());
    if(!splinesfile)
        AliFatal("File with splines for qn percentiles not found!");


    // for(int iDet=0; iDet<3; iDet++) {
    //     fSplineListqnPercTPC[iDet] = (TList*)splinesfile->Get(listnameTPC[iDet].Data());
    //     if(!fSplineListqnPercTPC[iDet])
    //         AliFatal("TList with splines for qnTPC percentiles not found in the spline file!");
    //     fSplineListqnPercV0[iDet] = (TList*)splinesfile->Get(listnameV0[iDet].Data());
    //     if(!fSplineListqnPercV0[iDet])
    //         AliFatal("TList with splines for qnV0 percentiles not found in the spline file!");

    //     fSplineListqnPercTPC[iDet]->SetOwner();
    //     fSplineListqnPercV0[iDet]->SetOwner();
    // }
    splinesfile->Close();
}


/**
 * This function performs a loop over the reconstructed jets
 * in the current event and fills the relevant histograms.
 */
void AliAnalysisTaskRawJetWithEP1::DoJetLoop()
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
    if (jetCont->GetRhoParameter()) { //kuma ??
      rhoVal = jetCont->GetRhoVal();
      histName = TString::Format("%s/hRhoVsCent", groupName.Data());
      fHistManager.FillTH2(histName.Data(), fCent, rhoVal);
    }

    
    for(auto jet : jetCont->accepted()) {
      if (!jet) continue;
      count++;
      
      histName = TString::Format("%s/hJetPt_%d", groupName.Data(), fCentBin);
      fHistManager.FillTH1(histName, jet->Pt());
      histName = TString::Format("%s/hJetArea_%d", groupName.Data(), fCentBin);
      fHistManager.FillTH1(histName, jet->Area());
      histName = TString::Format("%s/hJetPhi_%d", groupName.Data(), fCentBin);
      fHistManager.FillTH1(histName, jet->Phi());
      histName = TString::Format("%s/hJetEta_%d", groupName.Data(), fCentBin);
      fHistManager.FillTH1(histName, jet->Eta());

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


      // localRhoVal = fLocalRho->GetLocalVal(jet->Phi(), GetJetContainer()->GetJetRadius(), fLocalRho->GetVal());
      
      // localRhoValScaled = localRhoVal * jetCont->GetRhoVal() / p0;
      jetPtCorr = jet->Pt() - jetCont->GetRhoVal() * jet->Area();
      jetPtCorrLocal = jet->Pt() - localRhoValScaled * jet->Area();

      /*
      std::cout << "=s01= kumaaaaaaaaaakumaaaaaaaaaakumaaaaaaaaaakumaaaaaaaaaakumaaaaaaaaaa" << std::endl;
      std::cout << "(jetPhi, corrPhi2V0, p0, p1, localRho, jetA, localRhoValScal, phiMinusPsi2, phiMinusPsi3) = (" \
      << jet->Phi() << ", " << psi2V0[0] << ", " << p0 << ", " << p1 << ", "\
      << localRhoVal << ", " << jet->Area() << ", " << localRhoValScaled << ", " \
      << phiMinusPsi2 << ", " << phiMinusPsi3 \
      << ")" << std::endl;
      std::cout << "=e01= kumaaaaaaaaaakumaaaaaaaaaakumaaaaaaaaaakumaaaaaaaaaakumaaaaaaaaaa" << std::endl;
      */

      histName = TString::Format("%s/hJetRho_%d", groupName.Data(), fCentBin);
      fHistManager.FillTH1(histName, jetCont->GetRhoVal());
      histName = TString::Format("%s/hJetRhoLocal_%d", groupName.Data(), fCentBin);
      fHistManager.FillTH1(histName, localRhoVal); // trying out local rho val

      histName = TString::Format("%s/hJetCorrPt_%d", groupName.Data(), fCentBin);
      fHistManager.FillTH1(histName, jetPtCorr);
      histName = TString::Format("%s/hJetCorrPtLocal_%d", groupName.Data(), fCentBin);
      fHistManager.FillTH1(histName, jetPtCorrLocal);
      //fHistManager.FillTH1(histName, jet->Pt() - localRhoValScaled * 0.4);

      histName = TString::Format("%s/hJetLocalRhoVsAverageRho_%d", groupName.Data(), fCentBin);
      fHistManager.FillTH2(histName, jetCont->GetRhoVal(), localRhoValScaled);
      histName = TString::Format("%s/hJetCorrPtLocalVsJetCorrPt_%d", groupName.Data(), fCentBin);
      fHistManager.FillTH2(histName, jetPtCorr, jetPtCorrLocal);
      histName = TString::Format("%s/hJetRhoVsDeltaPhi_%d", groupName.Data(), fCentBin);
      fHistManager.FillTH2(histName, deltaPhiJetEP, localRhoValScaled);
      
      
      //V2 inclusive Jet
      histName = TString::Format("%s/hJetsYieldV2_%d", groupName.Data(), fCentBin);
      fHistManager.FillTH1(histName, jet->Pt());
      histName = TString::Format("%s/hCorrJetsYieldV2_%d", groupName.Data(), fCentBin);
      fHistManager.FillTH1(histName, jetPtCorr);
      histName = TString::Format("%s/hLocalCorrJetsYieldV2_%d", groupName.Data(), fCentBin);
      fHistManager.FillTH1(histName, jetPtCorrLocal);
      
      //V2 In plane Jet
      if ((phiMinusPsi2 < TMath::Pi()/4) || (phiMinusPsi2 >= 7*TMath::Pi()/4)\
      || (phiMinusPsi2 >= 3*TMath::Pi()/4 && phiMinusPsi2 < 5*TMath::Pi()/4)) {
        histName = TString::Format("%s/hJetsYieldInPlaneV2_%d", groupName.Data(), fCentBin);
        fHistManager.FillTH1(histName, jet->Pt());
        histName = TString::Format("%s/hCorrJetsYieldInPlaneV2_%d", groupName.Data(), fCentBin);
        fHistManager.FillTH1(histName, jetPtCorr);
        histName = TString::Format("%s/hLocalCorrJetsYieldInPlaneV2_%d", groupName.Data(), fCentBin);
        fHistManager.FillTH1(histName, jetPtCorrLocal);
      }
      else {
        histName = TString::Format("%s/hJetsYieldOutPlaneV2_%d", groupName.Data(), fCentBin);
        fHistManager.FillTH1(histName, jet->Pt());
        histName = TString::Format("%s/hCorrJetsYieldOutPlaneV2_%d", groupName.Data(), fCentBin);
        fHistManager.FillTH1(histName, jetPtCorr);
        histName = TString::Format("%s/hLocalCorrJetsYieldOutPlaneV2_%d", groupName.Data(), fCentBin);
        fHistManager.FillTH1(histName, jetPtCorrLocal);
      }
      
      /*
      //V3 In-plane
      histName = TString::Format("%s/hJetsYieldV3_%d", groupName.Data(), fCentBin);
      fHistManager.FillTH1(histName, jet->Pt());
      histName = TString::Format("%s/hCorrJetsYieldV3_%d", groupName.Data(), fCentBin);
      fHistManager.FillTH1(histName, jetPtCorr);
      histName = TString::Format("%s/hLocalCorrJetsYieldV3_%d", groupName.Data(), fCentBin);
      if ((phiMinusPsi3 < TMath::Pi()/6) \
      || (phiMinusPsi3 >= TMath::Pi()/2 && phiMinusPsi3 < 5*TMath::Pi()/6)\
      || (phiMinusPsi3 >= 7*TMath::Pi()/6 && phiMinusPsi3 < 3*TMath::Pi()/2)\
      || (phiMinusPsi3 >= 11*TMath::Pi()/6)) {
        histName = TString::Format("%s/hJetsYieldInPlaneV3_%d", groupName.Data(), fCentBin);
        fHistManager.FillTH1(histName, jet->Pt());
        histName = TString::Format("%s/hCorrJetsYieldInPlaneV3_%d", groupName.Data(), fCentBin);
        fHistManager.FillTH1(histName, jetPtCorr);
        histName = TString::Format("%s/hLocalCorrJetsYieldInPlaneV3_%d", groupName.Data(), fCentBin);
        fHistManager.FillTH1(histName, jetPtCorrLocal);
      }
      else {
        histName = TString::Format("%s/hJetsYieldOutPlaneV3_%d", groupName.Data(), fCentBin);
        fHistManager.FillTH1(histName, jet->Pt());
        histName = TString::Format("%s/hCorrJetsYieldOutPlaneV3_%d", groupName.Data(), fCentBin);
        fHistManager.FillTH1(histName, jetPtCorr);
        histName = TString::Format("%s/hLocalCorrJetsYieldOutPlaneV3_%d", groupName.Data(), fCentBin);
        fHistManager.FillTH1(histName, jetPtCorrLocal);
      }
      */
    }

    histName = TString::Format("%s/hNJets_%d", groupName.Data(), fCentBin);
    fHistManager.FillTH1(histName, count);
  }
}

Bool_t  AliAnalysisTaskRawJetWithEP1::QnJEHandlarEPGet(){
  fQ2VecHandler = new AliJEQnVectorHandler(0,1,2,fOADBFileName);
  fQ3VecHandler = new AliJEQnVectorHandler(0,1,3,fOADBFileName);
  
  fQ2VecHandler->ResetAODEvent();
  fQ2VecHandler->SetAODEvent(fAOD);
  fQ2VecHandler->ComputeCalibratedQnVectorTPC();
  fQ2VecHandler->ComputeCalibratedQnVectorV0();

  fQ3VecHandler->ResetAODEvent();
  fQ3VecHandler->SetAODEvent(fAOD);
  fQ3VecHandler->ComputeCalibratedQnVectorTPC();
  fQ3VecHandler->ComputeCalibratedQnVectorV0();

  //fill histos with EP angle
  fQ2VecHandler->GetEventPlaneAngleTPC(psi2Tpc[0],psi2Tpc[2],psi2Tpc[1]);
  fQ2VecHandler->GetEventPlaneAngleV0(psi2V0[0],psi2V0[2],psi2V0[1]);
  std::cout << "(Psi2FullTPC,Psi2PosTPC,Psi2NegTPC) = ("\
    << psi2Tpc[0] << "," << psi2Tpc[2] << "," << psi2Tpc[1] << ")" << std::endl;
  std::cout << "(Psi2FullV0,Psi2V0A,Psi2V0C) = ("\
    << psi2V0[0] << "," << psi2V0[2] << "," << psi2V0[1] << ")" << std::endl;

	//fill histos for q2 spline calibration
  fQ2VecHandler->GetQnVecTPC(q2VecTpcM, q2VecTpcP, q2VecTpcN);
  fQ2VecHandler->GetQnVecV0(q2VecV0M, q2VecV0A, q2VecV0C);
  fQ2VecHandler->GetqnTPC(q2Tpc[0],q2Tpc[2],q2Tpc[1]);
  fQ2VecHandler->GetqnV0(q2V0[0],q2Tpc[2],q2Tpc[1]);
  std::cout << "(q2FullTPC,q2PosTPC,q2NegTPC) = ("\
    << q2Tpc[0] << "," << q2Tpc[2] << "," << q2Tpc[1] << ")" << std::endl;
  std::cout << "(q2FullV0,q2V0A,q2V0C) = ("\
    << q2V0[0] << "," << q2V0[2] << "," << q2V0[1] << ")" << std::endl;


  //fill histos with EP angle
  fQ3VecHandler->GetEventPlaneAngleTPC(psi3Tpc[0],psi3Tpc[2],psi3Tpc[1]);
  fQ3VecHandler->GetEventPlaneAngleV0(psi3V0[0],psi3V0[2],psi3V0[1]);
  std::cout << "(Psi3FullTPC,Psi3PosTPC,Psi3NegTPC) = ("\
    << psi3Tpc[0] << "," << psi3Tpc[2] << "," << psi3Tpc[1] << ")" << std::endl;
  std::cout << "(Psi3FullV0,Psi3V0A,Psi3V0C) = ("\
    << psi3V0[0] << "," << psi3V0[2] << "," << psi3V0[1] << ")" << std::endl;

	//fill histos for q3 spline calibration
  fQ3VecHandler->GetQnVecTPC(q3VecTpcM, q3VecTpcP, q3VecTpcN);
  fQ3VecHandler->GetQnVecV0(q3VecV0M, q3VecV0A, q3VecV0C);
  fQ3VecHandler->GetqnTPC(q3Tpc[0],q3Tpc[2],q3Tpc[1]);
  fQ3VecHandler->GetqnV0(q3V0[0],q3V0[2],q3V0[1]);
  std::cout << "(q3FullTPC,q3PosTPC,q3NegTPC) = ("\
    << q3Tpc[0] << "," << q3Tpc[2] << "," << q3Tpc[1] << ")" << std::endl;
  std::cout << "(q3FullV0,q3V0A,q3V0C) = ("\
    << q3V0[0] << "," << q3V0[2] << "," << q3V0[1] << ")" << std::endl;
}

Bool_t  AliAnalysisTaskRawJetWithEP1::QnGainCalibration(){
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
  
  AliOADBContainer* cont = (AliOADBContainer*) fOADBFile->Get("hMultV0BefCorPfpx");
  TH1D* fHistMultV0 = ((TH1D*) cont->GetObject(fRunNumber));

  AliVEvent *fVevent = dynamic_cast<AliVEvent*>(InputEvent());
  const AliVVertex *pointVtx = fVevent->GetPrimaryVertex();
  Double_t fVtxZ = -999;
  fVtxZ  = pointVtx->GetZ();
  
  Int_t ibinV0 = 0;
  Double_t fSumMV0A = 0;
  Double_t fSumMV0C = 0;
  Double_t fV0chGain = 0.;
  Double_t fMultV0 = 0.;
  for(int iV0 = 0; iV0 < 64; iV0++) { //0-31 is V0C, 32-63 VOA
    fMultV0 = fAodV0->GetMultiplicity(iV0);
    
    if(CheckRunNum == fQaEventNum){
      histName = TString::Format("%s/hV0CellChVsMultBefEq", groupName.Data());
      fHistManager.FillTH1(histName, iV0, fMultV0);
    }

    /// V0 Channel Gain Correction:
    if(fHCorrV0ChWeghts){
      ibinV0    = fHCorrV0ChWeghts->FindBin(fVtxZ,iV0);
      fV0chGain = fHCorrV0ChWeghts->GetBinContent(ibinV0);
    }
    if(CheckRunNum == fQaEventNum){
      histName  = TString::Format("%s/hV0CellChGains", groupName.Data());
      fHistManager.FillTH1(histName, iV0, fV0chGain);
      std::cout << "Bef fMultV0 = " << fMultV0 << std::endl;
    }
    
    fMultV0 = fMultV0*fV0chGain;   //Corrected Multiplicity
    // std::cout << "fV0chGain = " << fV0chGain << std::endl;

    if(0){ //jet channel baias correction
      Double_t tagV0Mult = 0.;
      Double_t refV0Mult = 0.;
      Double_t refGainRatio = 0.;
      tagV0Mult = fHistMultV0->GetBinContent(iV0+1);
      // if (iV0 < 8) refV0Mult = fAodV0->GetMultiplicity(1);
      if (iV0 < 8) refV0Mult = fHistMultV0->GetBinContent(1);
      else if (iV0 >= 8  && iV0 < 16) refV0Mult = fHistMultV0->GetBinContent(9);
      else if (iV0 >= 16 && iV0 < 24) refV0Mult = fHistMultV0->GetBinContent(17);
      else if (iV0 >= 24 && iV0 < 32) refV0Mult = fHistMultV0->GetBinContent(25);
      else if (iV0 >= 32 && iV0 < 40) refV0Mult = fHistMultV0->GetBinContent(33);
      else if (iV0 >= 40 && iV0 < 48) refV0Mult = fHistMultV0->GetBinContent(41);
      else if (iV0 >= 48 && iV0 < 56) refV0Mult = fHistMultV0->GetBinContent(49);
      else if (iV0 >= 56 && iV0 < 64) refV0Mult = fHistMultV0->GetBinContent(57);

      refGainRatio = tagV0Mult/refV0Mult;
      fMultV0 = fMultV0/refGainRatio;   //Corrected Multiplicity

      if(CheckRunNum == fQaEventNum){
        // histName  = TString::Format("%s/hV0CellChRatio", groupName.Data());
        // fHistManager.FillTH1(histName, iV0, refGainRatio);
        histName = TString::Format("%s/hV0CellChVsMultAftEq", groupName.Data());
        fHistManager.FillTH1(histName, iV0, fMultV0);
        std::cout << "Aft fMultV0 = " << fMultV0 << std::endl;
        std::cout << "refMultV0 = " << refV0Mult << std::endl;
        std::cout << "ratio = " << refGainRatio << std::endl;
      }
    }

    Double_t fPhiV0  = TMath::PiOver4()*(0.5 + iV0 % 8);
    
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

Bool_t  AliAnalysisTaskRawJetWithEP1::QnRecenteringCalibration(){
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


/**
 * This function is called once at the end of the analysis.
 */
void AliAnalysisTaskRawJetWithEP1::Terminate(Option_t *) 
{
}


/**
 * This function adds the task to the analysis manager. Often, this function is called
 * by an AddTask C macro. However, by compiling the code, it ensures that we do not
 * have to deal with difficulties caused by CINT.
 */
AliAnalysisTaskRawJetWithEP1 * AliAnalysisTaskRawJetWithEP1::AddTaskRawJetWithEP1(
  const char *ntracks, const char *nclusters, const char* ncells, const char *suffix)
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskRawJetWithEP1", "No analysis manager to connect to.");
    return 0;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler)
  {
    ::Error("AddTaskRawJetWithEP1", "This task requires an input event handler");
    return 0;
  }
  
  enum EDataType_t {kUnknown, kESD, kAOD};
  EDataType_t dataType = kAOD;
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------
  TString trackName(ntracks);
  if (trackName == "usedefault") trackName = "tracks";

  TString name("AliAnalysisTaskRawJetWithEP1");
  if (strcmp(suffix,"") != 0) {
    name += "_";
    name += suffix;
  }

  AliAnalysisTaskRawJetWithEP1* rawJetTask = new AliAnalysisTaskRawJetWithEP1(name);
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













