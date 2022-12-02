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
#include "AliAnalysisTaskRawJetWithEP.h"

#include "AliDataFile.h"

class AliAnalysisTaskRawJetWithEP;

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskRawJetWithEP);
/// \endcond

/**
 * Default constructor. Needed by ROOT I/O
 */
AliAnalysisTaskRawJetWithEP::AliAnalysisTaskRawJetWithEP() :
  AliAnalysisTaskEmcalJet(),
  fAOD(nullptr),
  fUseAliEventCuts(kTRUE),
  fEventCuts(0),
  fEventCutList(0),
  fUseManualEventCuts(kFALSE),
  fYAMLConfig(),
  fUseRunList(),
  fOADBFileName(""),
  fOADBFile(nullptr),
  fRunListFileName(""),
  fSplinesFileName(""),
  fCalibRefFileName(""),
  fCalibRefFile(nullptr),
  fCalibRefObjList(nullptr),
  fCalibV0Ref(nullptr),
  fExLJetFromFit(kTRUE),
  fLeadingJet(0),
  fLeadingJetAfterSub(0),
  fFitModulation(0),
  fPileupCut(kFALSE),
  fTPCQnMeasure(kFALSE),
  fPileupCutQA(kFALSE),
  fCalibQA(kFALSE),
  fGainCalibQA(kFALSE),
  fReCentCalibQA(kFALSE),
  fEPQA(kFALSE),
  fTrackQA(kFALSE),
  fBkgQA(kFALSE),
  fSepEP(kFALSE),
  fCalibType(0),
  fNormMethod(0),
  fV0Combin(kFALSE),
  fQnVCalibType(kOrig),
  fQ2VecHandler(nullptr),
  fQ3VecHandler(nullptr),
  fHistManager(),
  fHCorrV0ChWeghts(NULL),
  fHCorrQ2xV0C(NULL),
  fHCorrQ2yV0C(NULL),
  fHCorrQ2xV0A(NULL),
  fHCorrQ2yV0A(NULL),
  fHCorrQ3xV0C(NULL),
  fHCorrQ3yV0C(NULL),
  fHCorrQ3xV0A(NULL),
  fHCorrQ3yV0A(NULL),
  fFitModulationType(kNoFit),
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
AliAnalysisTaskRawJetWithEP::AliAnalysisTaskRawJetWithEP(const char *name):
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fAOD(nullptr),
  fUseAliEventCuts(kTRUE),
  fEventCuts(0),
  fEventCutList(0),
  fUseManualEventCuts(kFALSE),
  fYAMLConfig(),
  fUseRunList(),
  fOADBFileName(""),
  fOADBFile(nullptr),
  fRunListFileName(""),
  fSplinesFileName(""),
  fCalibRefFileName(""),
  fCalibRefFile(nullptr),
  fCalibRefObjList(nullptr),
  fCalibV0Ref(nullptr),
  fExLJetFromFit(kTRUE),
  fLeadingJet(0),
  fLeadingJetAfterSub(0),
  fFitModulation(0),
  fPileupCut(kFALSE),
  fTPCQnMeasure(kFALSE),
  fPileupCutQA(kFALSE),
  fCalibQA(kFALSE),
  fGainCalibQA(kFALSE),
  fReCentCalibQA(kFALSE),
  fEPQA(kFALSE),
  fTrackQA(kFALSE),
  fBkgQA(kFALSE),
  fSepEP(kFALSE),
  fCalibType(0),
  fNormMethod(0),
  fV0Combin(kFALSE),
  fQnVCalibType(kOrig),
  fQ2VecHandler(nullptr),
  fQ3VecHandler(nullptr),
  fHistManager(name),
  fHCorrV0ChWeghts(NULL),
  fHCorrQ2xV0C(NULL),
  fHCorrQ2yV0C(NULL),
  fHCorrQ2xV0A(NULL),
  fHCorrQ2yV0A(NULL),
  fHCorrQ3xV0C(NULL),
  fHCorrQ3yV0C(NULL),
  fHCorrQ3xV0A(NULL),
  fHCorrQ3yV0A(NULL),
  fFitModulationType(kNoFit),
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
AliAnalysisTaskRawJetWithEP::~AliAnalysisTaskRawJetWithEP()
{
  
  if(fOADBFile) {fOADBFile->Close(); fOADBFile = 0x0;}
  if(fCalibRefFile) {fCalibRefFile->Close(); fCalibRefFile = 0x0;}
  if(fCalibV0Ref)   {delete fCalibV0Ref;   fCalibV0Ref = 0x0;}
  if(fCalibRefObjList) {delete fCalibRefObjList; fCalibRefObjList = 0x0;}
  if(fQ2VecHandler) {delete fQ2VecHandler; fQ2VecHandler = 0x0;}
  if(fQ3VecHandler) {delete fQ3VecHandler; fQ3VecHandler = 0x0;}
  if(fFitModulation) {delete fFitModulation; fFitModulation = 0x0;}
}

void AliAnalysisTaskRawJetWithEP::SetRunList(bool removeDummyTask)
{
  
  fYAMLConfig.AddConfiguration(fRunListFileName, "runlist");
  fYAMLConfig.Initialize();
  fYAMLConfig.GetProperty("runlist", fUseRunList);
}


/**
 * Performing run-independent initialization.
 * Here the histograms should be instantiated.
 */
void AliAnalysisTaskRawJetWithEP::UserCreateOutputObjects()
{
  
  // fOutputList = new TList();
  // fOutputList->SetOwner(true);

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();
  
  // Intialize AliEventCuts
  if (fUseAliEventCuts) {
    fEventCutList = new TList();
    fEventCutList ->SetOwner();
    fEventCutList ->SetName("EventCutOutput");
    
    fEventCuts.OverrideAutomaticTriggerSelection(AliVEvent::kAnyINT | AliVEvent::kCentral | AliVEvent::kSemiCentral);
    fEventCuts.SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE,1);
    if(fUseManualEventCuts==1)
    {
      fEventCuts.SetManualMode();
      fEventCuts.fMC = false;
      fEventCuts.SetupPbPb2018();
      fEventCuts.fUseVariablesCorrelationCuts = true;
    }
    fEventCuts.AddQAplotsToList(fEventCutList);
    fOutput->Add(fEventCutList);
  }
  
  // fEventCuts.AddQAplotsToList(fOutput);

  // == s == Set Out put Hist grams  ###########################################
  if(fEPQA)          AllocateEventPlaneHistograms();
  if(fBkgQA)         AllocateBkgHistograms();
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


void AliAnalysisTaskRawJetWithEP::AllocateEventPlaneHistograms()
{
  TString histName;
  TString histtitle;
  TString groupName;
  groupName="EventPlane";
  fHistManager.CreateHistoGroup(groupName);

  histName  = TString::Format("%s/hCentrality", groupName.Data());
  histtitle = TString::Format("%s;Centrality;counts", histName.Data());
  histtitle = TString::Format("%s;cell ch number;CorrGain", histName.Data());
  fHistManager.CreateTH1(histName, histtitle, 100, 0, 100);

  for (Int_t cent = 0; cent < fNcentBins; cent++) {
      // == s == Event plane angle histograms Setting
      histName = TString::Format("%s/hPsi2V0M_%d", groupName.Data(), cent);
      histtitle = "Psi2 from VZEROCA";
      //fHistManager.CreateTH1(histName, histtitle, 50, 0, TMath::TwoPi());
      fHistManager.CreateTH1(histName, histtitle, 50, -TMath::TwoPi(), 2*TMath::TwoPi());
      histName = TString::Format("%s/hPsi3V0M_%d", groupName.Data(), cent);
      histtitle = "Psi3 from V0M";
      //fHistManager.CreateTH1(histName, histtitle, 50, 0, TMath::TwoPi());
      fHistManager.CreateTH1(histName, histtitle, 50, -TMath::TwoPi(), 2*TMath::TwoPi());

      histName = TString::Format("%s/hPsi2V0C_%d", groupName.Data(), cent);
      histtitle = "Psi3 from V0A";
      //fHistManager.CreateTH1(histName, histtitle, 50, 0, TMath::TwoPi());
      fHistManager.CreateTH1(histName, histtitle, 50, -TMath::TwoPi(), 2*TMath::TwoPi());
      histName = TString::Format("%s/hPsi3V0C_%d", groupName.Data(), cent);
      histtitle = "Psi2 from V0C";
      //fHistManager.CreateTH1(histName, histtitle, 50, 0, TMath::TwoPi());
      fHistManager.CreateTH1(histName, histtitle, 50, -TMath::TwoPi(), 2*TMath::TwoPi());

      histName = TString::Format("%s/hPsi2V0A_%d", groupName.Data(), cent);
      histtitle = "Psi3 from V0A";
      //fHistManager.CreateTH1(histName, histtitle, 50, 0, TMath::TwoPi());
      fHistManager.CreateTH1(histName, histtitle, 50, -TMath::TwoPi(), 2*TMath::TwoPi());

      histName = TString::Format("%s/hPsi3V0A_%d", groupName.Data(), cent);
      histtitle = "Psi3 from V0A";
      //fHistManager.CreateTH1(histName, histtitle, 50, 0, TMath::TwoPi());
      fHistManager.CreateTH1(histName, histtitle, 50, -TMath::TwoPi(), 2*TMath::TwoPi());

      histName = TString::Format("%s/hPsi2TPCN_%d", groupName.Data(), cent);
      histtitle = "Psi2 from TPC eta negative";
      //fHistManager.CreateTH1(histName, histtitle, 50, 0, TMath::TwoPi());
      fHistManager.CreateTH1(histName, histtitle, 50, -TMath::TwoPi(), 2*TMath::TwoPi());
      histName = TString::Format("%s/hPsi3TPCN_%d", groupName.Data(), cent);
      histtitle = "Psi3 from TPC eta negative";
      //fHistManager.CreateTH1(histName, histtitle, 50, 0, TMath::TwoPi());
      fHistManager.CreateTH1(histName, histtitle, 50, -TMath::TwoPi(), 2*TMath::TwoPi());

      histName = TString::Format("%s/hPsi2TPCP_%d", groupName.Data(), cent);
      histtitle = "Psi2 from TPC eta positive";
      //fHistManager.CreateTH1(histName, histtitle, 50, 0, TMath::TwoPi());
      fHistManager.CreateTH1(histName, histtitle, 50, -TMath::TwoPi(), 2*TMath::TwoPi());
      histName = TString::Format("%s/hPsi3TPCP_%d", groupName.Data(), cent);
      histtitle = "Psi3 from TPC eta positive";
      //fHistManager.CreateTH1(histName, histtitle, 50, 0, TMath::TwoPi());
      fHistManager.CreateTH1(histName, histtitle, 50, -TMath::TwoPi(), 2*TMath::TwoPi());

      histName = TString::Format("%s/hPsi2V0MVsTPCP_%d", groupName.Data(), cent);
      histtitle = "Psi2 from V0M vs TPC postitive";
      fHistManager.CreateTH2(histName, histtitle, 50, -TMath::TwoPi(), 2*TMath::TwoPi(),\
        50, -TMath::TwoPi(), 2*TMath::TwoPi());
      histName = TString::Format("%s/hPsi2V0MVsTPCN_%d", groupName.Data(), cent);
      histtitle = "Psi2 from V0M vs TPC negative";
      fHistManager.CreateTH2(histName, histtitle, 50, -TMath::TwoPi(), 2*TMath::TwoPi(),\
        50, -TMath::TwoPi(), 2*TMath::TwoPi());
      histName = TString::Format("%s/hPsi2V0AVsV0C_%d", groupName.Data(), cent);
      histtitle = "Psi2 from V0A vs V0C";
      fHistManager.CreateTH2(histName, histtitle, 50, -TMath::TwoPi(), 2*TMath::TwoPi(),\
        50, -TMath::TwoPi(), 2*TMath::TwoPi());
      histName = TString::Format("%s/hPsi2TPCPVsTPCN_%d", groupName.Data(), cent);
      histtitle = "Psi2 from TPC posi vs nega";
      fHistManager.CreateTH2(histName, histtitle, 50, -TMath::TwoPi(), 2*TMath::TwoPi(),\
        50, -TMath::TwoPi(), 2*TMath::TwoPi());


      // profiles for all correlator permutations which are necessary to calculate each second and third order event plane resolution
      TProfile *tempHist;
      THashList *parent;
      histName = TString::Format("hProfV2Resolution_%d", cent);
      histtitle = histName;
      parent = static_cast<THashList *>(fHistManager.FindObject(groupName.Data()));
      tempHist = new TProfile(histName, histtitle, 11, -0.5, 10.5);
      tempHist->GetXaxis()->SetBinLabel(3, "<cos(2(#Psi_{VZEROA} - #Psi_{VZEROC}))>");
      tempHist->GetXaxis()->SetBinLabel(4, "<cos(2(#Psi_{VZEROC} - #Psi_{VZEROA}))>");
      tempHist->GetXaxis()->SetBinLabel(5, "<cos(2(#Psi_{VZEROA} - #Psi_{TPC}))>");
      tempHist->GetXaxis()->SetBinLabel(6, "<cos(2(#Psi_{TPC} - #Psi_{VZEROA}))>");
      tempHist->GetXaxis()->SetBinLabel(7, "<cos(2(#Psi_{VZEROC} - #Psi_{TPC}))>");
      tempHist->GetXaxis()->SetBinLabel(8, "<cos(2(#Psi_{TPC} - #Psi_{VZEROC}))>");
      tempHist->GetXaxis()->SetBinLabel(9, "<cos(2(#Psi_{VZERO} - #Psi_{TPC_A}))>");
      tempHist->GetXaxis()->SetBinLabel(10, "<cos(2(#Psi_{VZERO} - #Psi_{TPC_B}))>");
      tempHist->GetXaxis()->SetBinLabel(11, "<cos(2(#Psi_{TPC_A} - #Psi_{TPC_B}))>");
      parent->Add(tempHist);

      histName = TString::Format("hProfV3Resolution_%d", cent);
      histtitle = histName;
      parent = static_cast<THashList *>(fHistManager.FindObject(groupName.Data()));
      tempHist = new TProfile(histName, histtitle, 11, -0.5, 10.5);
      tempHist->GetXaxis()->SetBinLabel(3, "<cos(3(#Psi_{VZEROA} - #Psi_{VZEROC}))>");
      tempHist->GetXaxis()->SetBinLabel(4, "<cos(3(#Psi_{VZEROC} - #Psi_{VZEROA}))>");
      tempHist->GetXaxis()->SetBinLabel(5, "<cos(3(#Psi_{VZEROA} - #Psi_{TPC}))>");
      tempHist->GetXaxis()->SetBinLabel(6, "<cos(3(#Psi_{TPC} - #Psi_{VZEROA}))>");
      tempHist->GetXaxis()->SetBinLabel(7, "<cos(3(#Psi_{VZEROC} - #Psi_{TPC}))>");
      tempHist->GetXaxis()->SetBinLabel(8, "<cos(3(#Psi_{TPC} - #Psi_{VZEROC}))>");
      tempHist->GetXaxis()->SetBinLabel(9, "<cos(3(#Psi_{VZERO} - #Psi_{TPC_A}))>");
      tempHist->GetXaxis()->SetBinLabel(10, "<cos(3(#Psi_{VZERO} - #Psi_{TPC_B}))>");
      tempHist->GetXaxis()->SetBinLabel(11, "<cos(3(#Psi_{TPC_A} - #Psi_{TPC_B}))>");
      parent->Add(tempHist);
  }

  histName  = TString::Format("%s/v2", groupName.Data());
  histtitle = TString::Format("%s;centrality;v2", histName.Data());
  fHistManager.CreateTProfile(histName, histtitle, 10, 0, 10);
  histName  = TString::Format("%s/v3", groupName.Data());
  histtitle = TString::Format("%s;centrality;v3", histName.Data());
  fHistManager.CreateTProfile(histName, histtitle, 10, 0, 10);

}

void AliAnalysisTaskRawJetWithEP::AllocateBkgHistograms()
{
  TString histName;
  TString histtitle;
  TString groupName;
  groupName="BackgroundFit";
  fHistManager.CreateHistoGroup(groupName);

  // == s == cdf and pdf of chisquare distribution #############################
  //  = v2 and v3 combind fit local rho =
  histName = TString::Format("%s/hPvalueCDF_lRhoCombinFit", groupName.Data());
  histtitle = TString::Format("%s; CDF #chi^{2}", histName.Data());
  fHistManager.CreateTH1(histName, histtitle, 50, 0, 1);
  
  histName = TString::Format("%s/hPvalueCDFCent_lRhoCombinFit", groupName.Data());
  histtitle = TString::Format("%s; centrality; p-value", histName.Data());
  fHistManager.CreateTH2(histName, histtitle, 40, 0, 100, 40, 0, 1);
  
  histName = TString::Format("%s/hChi2Cent_lRhoCombinFit", groupName.Data());
  histtitle = TString::Format("%s; CDF #chi^{2}; #tilde{#chi^{2}}", histName.Data());
  fHistManager.CreateTH2(histName, histtitle, 100, 0, 100, 100, 0, 5);

  histName = TString::Format("%s/hPChi2_lRhoCombinFit", groupName.Data());
  histtitle = TString::Format("%s; p-value; #tilde{#chi^{2}}", histName.Data());
  fHistManager.CreateTH2(histName, histtitle, 100, 0, 1, 100, 0, 5);

  histName = TString::Format("%s/hPvalueCDFROOT_lRhoCombinFit", groupName.Data());
  histtitle = TString::Format("%s; CDF #chi^{2} ROOT", histName.Data());
  fHistManager.CreateTH1(histName, histtitle, 50, 0, 1);

  histName = TString::Format("%s/hPvalueCDFROOTCent_lRhoCombinFit", groupName.Data());
  histtitle = TString::Format("%s; centrality; p-value ROOT", histName.Data());
  fHistManager.CreateTH2(histName, histtitle, 40, 0, 100, 45, 0, 1);

  histName = TString::Format("%s/hChi2ROOTCent_lRhoCombinFit", groupName.Data());
  histtitle = TString::Format("%s; p-value; #tilde{#chi^{2}}", histName.Data());
  fHistManager.CreateTH2(histName, histtitle, 100, 0, 100, 45, 0, 1);

  histName = TString::Format("%s/hPChi2ROOT_lRhoCombinFit", groupName.Data());
  histtitle = TString::Format("%s;CDF #chi^{2}; #tilde{#chi^{2}} ROOT", histName.Data());
  fHistManager.CreateTH2(histName, histtitle, 100, 0, 1, 100, 0, 5);

  //  = v2  fit local rho =
  histName = TString::Format("%s/hPvalueCDF_lRhoV2Fit", groupName.Data());
  histtitle = TString::Format("%s; CDF #chi^{2}", histName.Data());
  fHistManager.CreateTH1(histName, histtitle, 50, 0, 1);
  
  histName = TString::Format("%s/hPvalueCDFCent_lRhoV2Fit", groupName.Data());
  histtitle = TString::Format("%s; centrality; p-value", histName.Data());
  fHistManager.CreateTH2(histName, histtitle, 40, 0, 100, 40, 0, 1);
  
  histName = TString::Format("%s/hChi2Cent_lRhoV2Fit", groupName.Data());
  histtitle = TString::Format("%s; CDF #chi^{2}; #tilde{#chi^{2}}", histName.Data());
  fHistManager.CreateTH2(histName, histtitle, 100, 0, 100, 100, 0, 5);

  histName = TString::Format("%s/hPChi2_lRhoV2Fit", groupName.Data());
  histtitle = TString::Format("%s; p-value; #tilde{#chi^{2}}", histName.Data());
  fHistManager.CreateTH2(histName, histtitle, 100, 0, 1, 100, 0, 5);

  histName = TString::Format("%s/hPvalueCDFROOT_lRhoV2Fit", groupName.Data());
  histtitle = TString::Format("%s; CDF #chi^{2} ROOT", histName.Data());
  fHistManager.CreateTH1(histName, histtitle, 50, 0, 1);

  histName = TString::Format("%s/hPvalueCDFROOTCent_lRhoV2Fit", groupName.Data());
  histtitle = TString::Format("%s; centrality; p-value ROOT", histName.Data());
  fHistManager.CreateTH2(histName, histtitle, 40, 0, 100, 45, 0, 1);

  histName = TString::Format("%s/hChi2ROOTCent_lRhoV2Fit", groupName.Data());
  histtitle = TString::Format("%s; p-value; #tilde{#chi^{2}}", histName.Data());
  fHistManager.CreateTH2(histName, histtitle, 100, 0, 100, 45, 0, 1);

  histName = TString::Format("%s/hPChi2ROOT_lRhoV2Fit", groupName.Data());
  histtitle = TString::Format("%s;CDF #chi^{2}; #tilde{#chi^{2}} ROOT", histName.Data());
  fHistManager.CreateTH2(histName, histtitle, 100, 0, 1, 100, 0, 5);

    //  = fit global rho (rho0) =
  histName = TString::Format("%s/hPvalueCDF_gRhoFit", groupName.Data());
  histtitle = TString::Format("%s; CDF #chi^{2}", histName.Data());
  fHistManager.CreateTH1(histName, histtitle, 50, 0, 1);
  
  histName = TString::Format("%s/hPvalueCDFCent_gRhoFit", groupName.Data());
  histtitle = TString::Format("%s; centrality; p-value", histName.Data());
  fHistManager.CreateTH2(histName, histtitle, 40, 0, 100, 40, 0, 1);
  
  histName = TString::Format("%s/hChi2Cent_gRhoFit", groupName.Data());
  histtitle = TString::Format("%s; CDF #chi^{2}; #tilde{#chi^{2}}", histName.Data());
  fHistManager.CreateTH2(histName, histtitle, 100, 0, 100, 100, 0, 5);

  histName = TString::Format("%s/hPChi2_gRhoFit", groupName.Data());
  histtitle = TString::Format("%s; p-value; #tilde{#chi^{2}}", histName.Data());
  fHistManager.CreateTH2(histName, histtitle, 100, 0, 1, 100, 0, 5);

  histName = TString::Format("%s/hPvalueCDFROOT_gRhoFit", groupName.Data());
  histtitle = TString::Format("%s; CDF #chi^{2} ROOT", histName.Data());
  fHistManager.CreateTH1(histName, histtitle, 50, 0, 1);

  histName = TString::Format("%s/hPvalueCDFROOTCent_gRhoFit", groupName.Data());
  histtitle = TString::Format("%s; centrality; p-value ROOT", histName.Data());
  fHistManager.CreateTH2(histName, histtitle, 40, 0, 100, 45, 0, 1);

  histName = TString::Format("%s/hChi2ROOTCent_gRhoFit", groupName.Data());
  histtitle = TString::Format("%s; p-value; #tilde{#chi^{2}}", histName.Data());
  fHistManager.CreateTH2(histName, histtitle, 100, 0, 100, 45, 0, 1);

  histName = TString::Format("%s/hPChi2ROOT_gRhoFit", groupName.Data());
  histtitle = TString::Format("%s;CDF #chi^{2}; #tilde{#chi^{2}} ROOT", histName.Data());
  fHistManager.CreateTH2(histName, histtitle, 100, 0, 1, 100, 0, 5);
  // == e == cdf and pdf of chisquare distribution #############################


  
  histName = TString::Format("%s/hRhoVsMult", groupName.Data());
  histtitle = TString::Format("%s; multiplicity; #rho [GeV/c]", histName.Data());
  fHistManager.CreateTH2(histName, histtitle, 100, 0, 4000, 100, 0, 250);
  
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
void AliAnalysisTaskRawJetWithEP::AllocateTrackHistograms()
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
      fHistManager.CreateTH1(histName, histtitle, 50, fMinBinPt, fMaxBinPt / 2);

      histName = TString::Format("%s/hTrackPhi_%d", groupName.Data(), cent);
      histtitle = TString::Format("%s;#it{#phi}_{track};counts", histName.Data());
      fHistManager.CreateTH1(histName, histtitle, 50, 0, TMath::TwoPi());

      histName = TString::Format("%s/hTrackEta_%d", groupName.Data(), cent);
      histtitle = TString::Format("%s;#it{#eta}_{track};counts", histName.Data());
      fHistManager.CreateTH1(histName, histtitle, fNbins / 6, -1, 1);

      histName = TString::Format("%s/hNTracks_%d", groupName.Data(), cent);
      histtitle = TString::Format("%s;number of tracks;events", histName.Data());
      fHistManager.CreateTH1(histName, histtitle, 500, 0, 5000);
      
      
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
void AliAnalysisTaskRawJetWithEP::AllocateJetHistograms()
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

    TString GenGroupName = TString::Format("%s/General", groupName.Data());
    if (fHistManager.FindObject(GenGroupName)) {
      AliWarning(TString::Format("%s: Found groupName %s in hist manager. The jet containers will be filled into the same histograms.", GetName(), GenGroupName.Data()));
      continue;
    }
    fHistManager.CreateHistoGroup(GenGroupName);

    TString RhoGroupName = TString::Format("%s/Rho", groupName.Data());
    if (fHistManager.FindObject(RhoGroupName)) {
      AliWarning(TString::Format("%s: Found groupName %s in hist manager. The jet containers will be filled into the same histograms.", GetName(), RhoGroupName.Data()));
      continue;
    }
    fHistManager.CreateHistoGroup(RhoGroupName);

    TString InclusiveGroupName = TString::Format("%s/Inclusive", groupName.Data());
    if (fHistManager.FindObject(InclusiveGroupName)) {
      AliWarning(TString::Format("%s: Found groupName %s in hist manager. The jet containers will be filled into the same histograms.", GetName(), InclusiveGroupName.Data()));
      continue;
    }
    fHistManager.CreateHistoGroup(InclusiveGroupName);

    TString IPlaneGroupName = TString::Format("%s/InPlane", groupName.Data());
    TString OPlaneGroupName = TString::Format("%s/OutOfPlane", groupName.Data());
    if(fSepEP){
      if (fHistManager.FindObject(IPlaneGroupName)) {
        AliWarning(TString::Format("%s: Found groupName %s in hist manager. The jet containers will be filled into the same histograms.", GetName(), IPlaneGroupName.Data()));
        continue;
      }
      fHistManager.CreateHistoGroup(IPlaneGroupName);
      
      if (fHistManager.FindObject(OPlaneGroupName)) {
        AliWarning(TString::Format("%s: Found groupName %s in hist manager. The jet containers will be filled into the same histograms.", GetName(), OPlaneGroupName.Data()));
        continue;
      }
      fHistManager.CreateHistoGroup(OPlaneGroupName);
    }


    // A vs. pT
    histName = TString::Format("%s/hRhoVsCent", RhoGroupName.Data());
    histtitle = histName + ";Centrality (%);#rho (GeV/#it{c});counts";
    fHistManager.CreateTH2(histName, histtitle.Data(), 50, 0, 100, 100, 0, 500);

    
    for (Int_t cent = 0; cent < fNcentBins; cent++) {
      histName = TString::Format("%s/hJetArea_%d", GenGroupName.Data(), cent);
      histtitle = TString::Format("%s;#it{A}_{jet};counts", histName.Data());
      fHistManager.CreateTH1(histName, histtitle, 50, 0, 3);

      histName = TString::Format("%s/hJetPhi_%d", GenGroupName.Data(), cent);
      histtitle = TString::Format("%s;#it{#phi}_{jet};counts", histName.Data());
      fHistManager.CreateTH1(histName, histtitle, 50, 0, TMath::TwoPi());

      histName = TString::Format("%s/hJetEta_%d", GenGroupName.Data(), cent);
      histtitle = TString::Format("%s;#it{#eta}_{jet};counts", histName.Data());
      fHistManager.CreateTH1(histName, histtitle, 50, -1, 1);

      histName = TString::Format("%s/hNJets_%d", GenGroupName.Data(), cent);
      histtitle = TString::Format("%s;number of jets;events", histName.Data());
      fHistManager.CreateTH1(histName, histtitle, 500, 0, 500);

      // histograms for jet angle relative to the event plane
      histName = TString::Format("%s/hJetPhiMinusPsi2_%d", GenGroupName.Data(), cent);
      histtitle = "Jet phi minus psi2";
      fHistManager.CreateTH1(histName, histtitle, 50, 0, TMath::TwoPi());
      histName = TString::Format("%s/hJetPhiMinusPsi3_%d", GenGroupName.Data(), cent);
      histtitle = "Jet phi minus psi3";
      fHistManager.CreateTH1(histName, histtitle, 50, 0, TMath::TwoPi());

      histName = TString::Format("%s/hDeltaPt_%d", GenGroupName.Data(), cent);
      histtitle = "delta pt without EP distribution";
      fHistManager.CreateTH1(histName, histtitle, 100, -50, 50);
      histName = TString::Format("%s/hDeltaPt_Local_%d", GenGroupName.Data(), cent);
      histtitle = "delta pt according EP distribution";
      fHistManager.CreateTH1(histName, histtitle, 100, -50, 50);
      histName = TString::Format("%s/hPhiVsDeltaPt_Global_%d", GenGroupName.Data(), cent);
      histtitle = "phi vs delta pt according EP distribution";
      fHistManager.CreateTH2(histName, histtitle, 50, 0, TMath::TwoPi(), 100, -50, 50);
      histName = TString::Format("%s/hPhiVsDeltaPt_Local_%d", GenGroupName.Data(), cent);
      histtitle = "phi vs delta pt according EP distribution";
      fHistManager.CreateTH2(histName, histtitle, 50, 0, TMath::TwoPi(), 100, -50, 50);


      if (!jetCont->GetRhoName().IsNull()) {
        
        // == s ==  Rho histograms   11111111111111111111111111111111111111111111111111111111111111
        histtitle = "GlobalRho";
        histName = TString::Format("%s/hJetRho_%d", RhoGroupName.Data(), cent);
        fHistManager.CreateTH1(histName, histtitle, 60, 0.0, 300.0);
        histtitle = "LocalRho";
        histName = TString::Format("%s/hJetRhoLocal_%d", RhoGroupName.Data(), cent);
        fHistManager.CreateTH1(histName, histtitle, 60, 0.0, 300.0);
        
        histName = TString::Format("%s/hJetLRhoVsAveRho_%d", RhoGroupName.Data(), cent);
        histtitle = "Local rho versus average rho";
        fHistManager.CreateTH2(histName, histtitle, 60, 0.0, 300.0, 60, 0.0, 300.0);

        // histo local rho vs delta phi
        histName = TString::Format("%s/hJetGRhoVsDeltaPhi_%d", RhoGroupName.Data(), cent);
        histtitle = "Global rho versus angle relative to event plane";
        fHistManager.CreateTH2(histName, histtitle, 60, 0.0, TMath::TwoPi(), 60, 0.0, 300.0);
        histName = TString::Format("%s/hJetLRhoVsDeltaPhi_%d", RhoGroupName.Data(), cent);
        histtitle = "Local rho versus angle relative to event plane";
        fHistManager.CreateTH2(histName, histtitle, 60, 0.0, TMath::TwoPi(), 60, 0.0, 300.0);
        // == e ==  Rho histograms   11111111111111111111111111111111111111111111111111111111111111


        // = s = Create histograms of Jet Yeild ====================================================
        //inlucive
        histName = TString::Format("%s/hJetPt_%d", InclusiveGroupName.Data(), cent);
        histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c});counts", histName.Data());
        fHistManager.CreateTH1(histName, histtitle, 60, -50, 250);
        histName = TString::Format("%s/hJetCorrPt_%d", InclusiveGroupName.Data(), cent);
        histtitle = TString::Format("%s;#it{p}_{T,jet}^{corr} (GeV/#it{c});counts", histName.Data());
        fHistManager.CreateTH1(histName, histtitle, 60, -50, 250);
        
        if(fSepEP){
          //v2 in plane
          histName = TString::Format("%s/hJetPt_%d", IPlaneGroupName.Data(), cent);
          histtitle = "Jet yeild of in-plane (v2)";
          fHistManager.CreateTH1(histName, histtitle, 60, -50, 250);
          histName = TString::Format("%s/hJetCorrPtLocal_%d", IPlaneGroupName.Data(), cent);
          histtitle = "corr Jet yeild of in-plane (v2)";
          fHistManager.CreateTH1(histName, histtitle, 60, -50, 250);

          //v2 out of plane
          histName = TString::Format("%s/hJetPt_%d", OPlaneGroupName.Data(), cent);
          histtitle = "Jet yeild of out-plane (v2)";
          fHistManager.CreateTH1(histName, histtitle, 60, -50, 250);
          histName = TString::Format("%s/hJetCorrPtLocal_%d", OPlaneGroupName.Data(), cent);
          histtitle = "corr Jet yeild of in-plane (v2)";
          fHistManager.CreateTH1(histName, histtitle, 60, -50, 250);
        }
        // = e = Create histograms of Jet Yeild ===================================================
        
      }
    }
  }
}


/**
 * This function is executed automatically for the first event.
 * Some extra initialization can be performed here.
 */
void AliAnalysisTaskRawJetWithEP::ExecOnce()
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
 * This function (overloading the base class) uses AliEventCuts to perform event selection.
 */
Bool_t AliAnalysisTaskRawJetWithEP::IsEventSelected()
{
  if (fUseAliEventCuts) {
    if (!fEventCuts.AcceptEvent(InputEvent()))
    {
      PostData(1, fOutput);
      return kFALSE;
    }
  }
  else {
    Bool_t answer = AliAnalysisTaskEmcal::IsEventSelected();
    return answer;
  }
  return kTRUE;
}

/**
 * Run analysis code here, if needed.
 * It will be executed before FillHistograms().
 * If this function return kFALSE, FillHistograms() will *not*
 * be executed for the current event
 * @return Always kTRUE
 */
Bool_t AliAnalysisTaskRawJetWithEP::Run()
{
  
  if(!fEventCuts.AcceptEvent(InputEvent())) return kFALSE;
  
  if(fPileupCut){
    SetupPileUpRemovalFunctions();
    Bool_t kPileupCutEvent = CheckEventIsPileUp2018();
    if(kPileupCutEvent) return kFALSE;
  }
  
  if(!fLocalRho) {
    AliWarning(Form("%s: No LocalRho object found, attempting to get it from Event based on name!",GetName()));
    fLocalRho = GetLocalRhoFromEvent(fLocalRhoName);
  }

  DoEventPlane();
  if(fTPCQnMeasure) MeasureTpcEPQA();
  MeasureBkg();
  DoJetLoop();

  // if(fLocalRho) delete fLocalRho;
  return kTRUE;
}

void AliAnalysisTaskRawJetWithEP::DoEventPlane(){
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
  switch (fQnVCalibType) {
    case kOrig : {
      // std::cout << "bef calib (qx,qy) = " << q2VecV0M[0] << "," << q2VecV0M[1] << ")" << std::endl;
      QnGainCalibration();
      if(0){
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
    } break;
    case kJeHand : {
      QnJEHandlarEPGet(); 
    } break;
    default : break;
  }
  //== e == qn Calibration  111111111111111111111111111111111111111111111111111
  
  
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

  histName = TString::Format("%s/hPsi2V0MVsTPCP_%d", groupName.Data(), fCentBin);
  fHistManager.FillTH2(histName, psi2V0[0], psi2Tpc[1]);
  histName = TString::Format("%s/hPsi2V0MVsTPCN_%d", groupName.Data(), fCentBin);
  fHistManager.FillTH2(histName, psi2V0[0], psi2Tpc[2]);
  histName = TString::Format("%s/hPsi2V0AVsV0C_%d", groupName.Data(), fCentBin);
  fHistManager.FillTH2(histName, psi2V0[2], psi2V0[1]);
  histName = TString::Format("%s/hPsi2TPCPVsTPCN_%d", groupName.Data(), fCentBin);
  fHistManager.FillTH2(histName, psi2V0[1], psi2Tpc[2]);


  histName = TString::Format("%s/hProfV2Resolution_%d", groupName.Data(), fCentBin);
  fHistManager.FillProfile(histName, 2., TMath::Cos(2.*(psi2V0[2] - psi2V0[1])));
  fHistManager.FillProfile(histName, 3., TMath::Cos(2.*(psi2V0[1] - psi2V0[2])));
  fHistManager.FillProfile(histName, 4., TMath::Cos(2.*(psi2V0[2] - psi2Tpc[0])));
  fHistManager.FillProfile(histName, 5., TMath::Cos(2.*(psi2Tpc[2] - psi2V0[2])));
  fHistManager.FillProfile(histName, 6., TMath::Cos(2.*(psi2V0[1] - psi2Tpc[0])));
  fHistManager.FillProfile(histName, 7., TMath::Cos(2.*(psi2Tpc[0] - psi2V0[1])));
  fHistManager.FillProfile(histName, 8., TMath::Cos(2.*(psi2V0[0] - psi2Tpc[1])));
  fHistManager.FillProfile(histName, 9., TMath::Cos(2.*(psi2V0[0] - psi2Tpc[2])));
  fHistManager.FillProfile(histName, 10., TMath::Cos(2.*(psi2Tpc[1] - psi2Tpc[2])));
    
  histName = TString::Format("%s/hProfV3Resolution_%d", groupName.Data(), fCentBin);
  fHistManager.FillProfile(histName, 2., TMath::Cos(3.*(psi2V0[2] - psi2V0[1])));
  fHistManager.FillProfile(histName, 3., TMath::Cos(3.*(psi2V0[1] - psi2V0[2])));
  fHistManager.FillProfile(histName, 4., TMath::Cos(3.*(psi2V0[2] - psi2Tpc[0])));
  fHistManager.FillProfile(histName, 5., TMath::Cos(3.*(psi2Tpc[2] - psi2V0[2])));
  fHistManager.FillProfile(histName, 6., TMath::Cos(3.*(psi2V0[1] - psi2Tpc[0])));
  fHistManager.FillProfile(histName, 7., TMath::Cos(3.*(psi2Tpc[0] - psi2V0[1])));
  fHistManager.FillProfile(histName, 8., TMath::Cos(3.*(psi2V0[0] - psi2Tpc[1])));
  fHistManager.FillProfile(histName, 9., TMath::Cos(3.*(psi2V0[0] - psi2Tpc[2])));
  fHistManager.FillProfile(histName, 10., TMath::Cos(3.*(psi2Tpc[1] - psi2Tpc[2])));
}


Bool_t AliAnalysisTaskRawJetWithEP::MeasureBkg(){
  TString groupName;
  TString histName;

  // AliAnalysisTaskJetV2 ==================================================================
  Int_t iTracks(fTracks->GetEntriesFast());
  Double_t excludeInEta = -999;
  Double_t excludeInPhi = -999;
  Double_t excludeInPt  = -999;
  // if(fLocalRho->GetVal() <= 0 ) return kFALSE;   // no use fitting an empty event ...
  if(fExLJetFromFit) {
      if(fLeadingJet) {
          excludeInEta = fLeadingJet->Eta();
          excludeInPhi = fLeadingJet->Phi();
          excludeInPt = fLeadingJet->Pt();
      }
  }

  // check the acceptance of the track selection that will be used
  // if one uses e.g. semi-good tpc tracks, accepance in phi is reduced to 0 < phi < 4
  // the defaults (-10 < phi < 10) which accept all, are then overwritten
  Double_t lowBound(0.), upBound(TMath::TwoPi());     // bounds for fit
  if(GetParticleContainer()->GetParticlePhiMin() > lowBound){
    lowBound = GetParticleContainer()->GetParticlePhiMin();
  }
  if(GetParticleContainer()->GetParticlePhiMax() < upBound){
    upBound = GetParticleContainer()->GetParticlePhiMax();
  }
  
  TH1F _tempSwap;     // on stack for quick access
  TH1F _tempSwapN;    // on stack for quick access, bookkeeping histogram

  // non poissonian error when using pt weights
  Double_t sumPt(0.), sumPt2(0.), trackN(0.);
  Double_t tempJetR = 0.;
  tempJetR = GetJetContainer()->GetJetRadius();

  AliParticleContainer* partCont = 0;
  TIter next(&fParticleCollArray);
  while ((partCont = static_cast<AliParticleContainer*>(next()))) {
      for(auto part : partCont->accepted()) {
        if (!part) continue;
        if (partCont->GetLoadedClass()->InheritsFrom("AliVTrack")) {
          const AliVTrack* track = static_cast<const AliVTrack*>(part);
          
          if(( (TMath::Abs(track->Eta() - excludeInEta) < tempJetR) \
          || (TMath::Abs(track->Eta()) - tempJetR - GetJetContainer()->GetJetEtaMax() ) > 0 )) continue;
          
          _tempSwap.Fill(track->Phi(), track->Pt());
          
          sumPt += track->Pt();
          sumPt2 += track->Pt()*track->Pt();
          trackN += 1;
          _tempSwapN.Fill(track->Phi());
          
        }
      }
  }
  
  // in the case of pt weights overwrite the poissonian error estimate
  // which is assigned by root by a more sophisticated appraoch
  // the assumption here is that the bin error will be dominated 
  // by the uncertainty in the mean pt in a bin and in the uncertainty
  // of the number of tracks in a bin, the first of which will be estimated 
  // from the sample standard deviation of all tracks in the 
  // event, for the latter use a poissonian estimate. 
  // the two contrubitions are assumed to be uncorrelated
  // not one track passes the cuts > 2 avoids possible division by 0 later on
  if(trackN < 2) return kFALSE; 
  for(Int_t l = 0; l < _tempSwap.GetNbinsX(); l++) {
      if(_tempSwapN.GetBinContent(l+1) == 0) {
          _tempSwap.SetBinContent(l+1,0);
          _tempSwap.SetBinError(l+1,0);
      }
      else {
          Double_t vartimesnsq = sumPt2*trackN - sumPt*sumPt;
          Double_t variance = vartimesnsq/(trackN*(trackN-1.));
          Double_t SDOMSq = variance / _tempSwapN.GetBinContent(l+1);
          Double_t SDOMSqOverMeanSq \
            = SDOMSq * _tempSwapN.GetBinContent(l+1) * _tempSwapN.GetBinContent(l+1) \
              / (_tempSwapN.GetBinContent(l+1) * _tempSwapN.GetBinContent(l+1));
          Double_t poissonfrac = 1./_tempSwapN.GetBinContent(l+1);
          Double_t vartotalfrac = SDOMSqOverMeanSq + poissonfrac;
          Double_t vartotal \
            = vartotalfrac * _tempSwap.GetBinContent(l+1) * _tempSwap.GetBinContent(l+1);
          if(vartotal > 0.0001) _tempSwap.SetBinError(l+1,TMath::Sqrt(vartotal));
          else {
              _tempSwap.SetBinContent(l+1,0);
              _tempSwap.SetBinError(l+1,0);
          }
      }
  }

  /// === s === determine background fit function   ###############################################
  // TF1* fFitModulation = 0x0;
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
  /// === e === determine background fit function   ###############################################

  fLocalRho->SetVal(fRho->GetVal());
  fFitModulation->SetParameter(0, fLocalRho->GetVal());
  fFitModulation->FixParameter(2, psi2V0[0]);
  fFitModulation->FixParameter(4, psi3V0[0]);

  _tempSwap.Fit(fFitModulation, "N0Q", "", lowBound, upBound);

  Double_t tempV2 = -999.9;
  Double_t tempV3 = -999.9;
  tempV2 = fFitModulation->GetParameter(1);
  tempV3 = fFitModulation->GetParameter(3);
  
  groupName="EventPlane";
  histName = TString::Format("%s/v2", groupName.Data());
  fHistManager.FillProfile(histName, fCentBin, tempV2);
  histName = TString::Format("%s/v3", groupName.Data());
  fHistManager.FillProfile(histName, fCentBin, tempV3);
  
  // AliAnalysisTaskJetV2 ===========================================================================
  
  if(0){
    TCanvas *cBkgRhoFit = new TCanvas("cBkgRhoFit", "cBkgRhoFit", 2000, 1500);
    
    TH1F* hBkgTracks_Event = (TH1F*) _tempSwap.Clone("hnew");
    hBkgTracks_Event->SetName("hBkgTracks");
    hBkgTracks_Event->Draw("E");
    _tempSwap.Fit(fFitModulation, "N0Q");
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
    rhoFitV3Com->SetParameter(1, fFitModulation->GetParameter(3));//v3
    rhoFitV3Com->SetParameter(2, fFitModulation->GetParameter(4));//psi3
    rhoFitV3Com->SetLineColor(824);
    rhoFitV3Com->Draw("same");

    histName = "checkOutput/cBkgRhoFit_Cent" + std::to_string(fCentBin) +".root";
    cBkgRhoFit->SaveAs(histName);

    histName = "checkOutput/fFitModulation_Cent" + std::to_string(fCentBin) +".root";
    fFitModulation->SaveAs(histName);
    histName = "checkOutput/v2Fit_Cent" + std::to_string(fCentBin) +".root";
    rhoFitV2Com->SaveAs(histName);
    histName = "checkOutput/v3Fit_Cent" + std::to_string(fCentBin) +".root";
    rhoFitV3Com->SaveAs(histName);
    
    // histName = "checkOutput/hBkgTracks_Event" + std::to_string(CheckRunNum) +".root";
    hBkgTracks_Event->SaveAs(histName);
    delete cBkgRhoFit;
    delete hBkgTracks_Event;
  }
  
  // fV2ResoV0 = CalcEPReso(2, psi2V0[0], psi2Tpc[1], psi2Tpc[2]);
  // fV3ResoV0 = CalcEPReso(3, psi3V0[0], psi3Tpc[1], psi3Tpc[2]);
  // std::cout << "v2Reso = " << fV2ResoV0 << ", v3Reso = " << fV3ResoV0 << std::endl;

  fLocalRho->SetLocalRho(fFitModulation);
  BkgFitEvaluation(&_tempSwap, fFitModulation);

  // delete fFitModulation;
  
  return kTRUE;
}


Double_t AliAnalysisTaskRawJetWithEP::CalcEPReso(Int_t n, \
  Double_t &psiA, Double_t &psiB, Double_t &psiC){
  
  Double_t vnReso = -999.;
  vnReso = TMath::Sqrt((TMath::Abs(TMath::Cos(n*(psiA - psiB))) \
                          * TMath::Abs(TMath::Cos(n*(psiA - psiC)))) \
                        / TMath::Abs(TMath::Cos(n*(psiB - psiC))));

  return vnReso;
}


void AliAnalysisTaskRawJetWithEP::MeasureTpcEPQA(){
  TString histName;
  TString groupName;

  Double_t EtaAcc = 0.9;
  UInt_t sumAcceptedTracks = 0;
  AliParticleContainer* partCont = 0;
  
  TIter next(&fParticleCollArray);
  
  while ((partCont = static_cast<AliParticleContainer*>(next()))) {
    
    groupName = partCont->GetName();
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
void AliAnalysisTaskRawJetWithEP::LoadSpliForqnPerce()
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
void AliAnalysisTaskRawJetWithEP::DoJetLoop()
{
  TString histName;
  TString groupName;
  AliJetContainer* jetCont = 0;
  TIter next(&fJetCollArray);
  
  while ((jetCont = static_cast<AliJetContainer*>(next()))) {
    groupName = jetCont->GetName();
    TString GenGroupName = TString::Format("%s/General", groupName.Data());
    TString RhoGroupName = TString::Format("%s/Rho", groupName.Data());
    TString InclusiveGroupName = TString::Format("%s/Inclusive", groupName.Data());
    TString IPlaneGroupName = TString::Format("%s/InPlane", groupName.Data());
    TString OPlaneGroupName = TString::Format("%s/OutOfPlane", groupName.Data());
    UInt_t count = 0;
    
    Double_t jetR = jetCont->GetJetRadius();
    Double_t rhoVal = 0;
    if (jetCont->GetRhoParameter()) { //kuma ??
      rhoVal = jetCont->GetRhoVal();
      histName = TString::Format("%s/hRhoVsCent", RhoGroupName.Data());
      fHistManager.FillTH2(histName.Data(), fCent, rhoVal);
    }

    Double_t leadingJetEta = -999.;
    Double_t leadingJetPhi = -999.;
    Double_t leadingJetPt  = -999.;
    
    for(auto jet : jetCont->accepted()) {
      if (!jet) continue;
      count++;
      
      histName = TString::Format("%s/hNJets_%d", GenGroupName.Data(), fCentBin);
      fHistManager.FillTH1(histName, count);
      histName = TString::Format("%s/hJetArea_%d", GenGroupName.Data(), fCentBin);
      fHistManager.FillTH1(histName, jet->Area());
      histName = TString::Format("%s/hJetPhi_%d", GenGroupName.Data(), fCentBin);
      fHistManager.FillTH1(histName, jet->Phi());
      histName = TString::Format("%s/hJetEta_%d", GenGroupName.Data(), fCentBin);
      fHistManager.FillTH1(histName, jet->Eta());

      // Filling histos for angle relative to event plane
      Double_t phiMinusPsi2 = jet->Phi() - psi2V0[0];
      if (phiMinusPsi2 < 0.0) phiMinusPsi2 += TMath::TwoPi();
      histName = TString::Format("%s/hJetPhiMinusPsi2_%d", GenGroupName.Data(), fCentBin);
      fHistManager.FillTH1(histName, phiMinusPsi2);
      Double_t phiMinusPsi3 = jet->Phi() - psi3V0[0];
      if (phiMinusPsi3 < 0.0) phiMinusPsi3 += TMath::TwoPi();
      histName = TString::Format("%s/hJetPhiMinusPsi3_%d", GenGroupName.Data(), fCentBin);
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
      
      histName = TString::Format("%s/hJetRho_%d", RhoGroupName.Data(), fCentBin);
      fHistManager.FillTH1(histName, jetCont->GetRhoVal());
      histName = TString::Format("%s/hJetRhoLocal_%d", RhoGroupName.Data(), fCentBin);
      fHistManager.FillTH1(histName, localRhoVal); // trying out local rho val

      histName = TString::Format("%s/hJetLRhoVsAveRho_%d", RhoGroupName.Data(), fCentBin);
      fHistManager.FillTH2(histName, jetCont->GetRhoVal(), localRhoValScaled);
      histName = TString::Format("%s/hJetGRhoVsDeltaPhi_%d", RhoGroupName.Data(), fCentBin);
      fHistManager.FillTH2(histName, deltaPhiJetEP, jetCont->GetRhoVal());
      histName = TString::Format("%s/hJetLRhoVsDeltaPhi_%d", RhoGroupName.Data(), fCentBin);
      fHistManager.FillTH2(histName, deltaPhiJetEP, localRhoValScaled);
      
      Double_t rcPt = 0., rcEta = 0., rcPhi = 0.;
      CalcRandomCone(rcPt, rcEta, rcPhi, leadingJetEta, leadingJetPhi, jetR);    
      Double_t rcLocalRhoValScaled = fLocalRho->GetLocalVal(rcPhi, jetR, fLocalRho->GetVal());
      Double_t deltaLoacalPt = rcPt - rcLocalRhoValScaled*jetR*jetR*TMath::Pi();
      Double_t deltaGlobalPt = rcPt - fLocalRho->GetVal()*jetR*jetR*TMath::Pi();
      histName = TString::Format("%s/hDeltaPt_%d", GenGroupName.Data(), fCentBin);
      fHistManager.FillTH1(histName, deltaGlobalPt);
      histName = TString::Format("%s/hDeltaPt_Local_%d", GenGroupName.Data(), fCentBin);
      fHistManager.FillTH1(histName, deltaLoacalPt);
      histName = TString::Format("%s/hPhiVsDeltaPt_Global_%d", GenGroupName.Data(), fCentBin);
      fHistManager.FillTH2(histName, rcPhi, deltaGlobalPt);
      histName = TString::Format("%s/hPhiVsDeltaPt_Local_%d", GenGroupName.Data(), fCentBin);
      fHistManager.FillTH2(histName, rcPhi, deltaLoacalPt);
      
      //inclusive Jet
      histName = TString::Format("%s/hJetPt_%d", InclusiveGroupName.Data(), fCentBin);
      fHistManager.FillTH1(histName, jet->Pt());
      histName = TString::Format("%s/hJetCorrPt_%d", InclusiveGroupName.Data(), fCentBin);
      fHistManager.FillTH1(histName, jetPtCorr);
      
      //V2 In plane Jet
      if(fSepEP){
        if ((phiMinusPsi2 < TMath::Pi()/4) || (phiMinusPsi2 >= 7*TMath::Pi()/4)\
        || (phiMinusPsi2 >= 3*TMath::Pi()/4 && phiMinusPsi2 < 5*TMath::Pi()/4)) {
          histName = TString::Format("%s/hJetPt_%d", IPlaneGroupName.Data(), fCentBin);
          fHistManager.FillTH1(histName, jet->Pt());
          histName = TString::Format("%s/hJetCorrPtLocal_%d", IPlaneGroupName.Data(), fCentBin);
          fHistManager.FillTH1(histName, jetPtCorrLocal);
        }
        else {
          histName = TString::Format("%s/hJetPt_%d", OPlaneGroupName.Data(), fCentBin);
          fHistManager.FillTH1(histName, jet->Pt());
          histName = TString::Format("%s/hJetCorrPtLocal_%d", OPlaneGroupName.Data(), fCentBin);
          fHistManager.FillTH1(histName, jetPtCorrLocal);
        }
      }
      
    }
  }
}

void AliAnalysisTaskRawJetWithEP::QnJEHandlarEPGet()
{
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
  if(0){
    std::cout << "(Psi2FullTPC,Psi2PosTPC,Psi2NegTPC) = ("\
      << psi2Tpc[0] << "," << psi2Tpc[2] << "," << psi2Tpc[1] << ")" << std::endl;
    std::cout << "(Psi2FullV0,Psi2V0A,Psi2V0C) = ("\
      << psi2V0[0] << "," << psi2V0[2] << "," << psi2V0[1] << ")" << std::endl;
  }

	//fill histos for q2 spline calibration
  fQ2VecHandler->GetQnVecTPC(q2VecTpcM, q2VecTpcP, q2VecTpcN);
  fQ2VecHandler->GetQnVecV0(q2VecV0M, q2VecV0A, q2VecV0C);
  fQ2VecHandler->GetqnTPC(q2Tpc[0],q2Tpc[2],q2Tpc[1]);
  fQ2VecHandler->GetqnV0(q2V0[0],q2V0[2],q2V0[1]);
  if(0){
    std::cout << "recent calibA (qx,qy) = " \
      << q2VecV0A[0] << "," << q2VecV0A[1] << ")" << std::endl;
    std::cout << "recent calibM (qx,qy) = " \
      << q2VecV0M[0] << "," << q2VecV0M[1] << ")" << std::endl;
    std::cout << "(q2FullTPC,q2PosTPC,q2NegTPC) = ("\
      << q2Tpc[0] << "," << q2Tpc[2] << "," << q2Tpc[1] << ")" << std::endl;
    std::cout << "(q2FullV0,q2V0A,q2V0C) = ("\
      << q2V0[0] << "," << q2V0[2] << "," << q2V0[1] << ")" << std::endl;    
  }

  //fill histos with EP angle
  fQ3VecHandler->GetEventPlaneAngleTPC(psi3Tpc[0],psi3Tpc[2],psi3Tpc[1]);
  fQ3VecHandler->GetEventPlaneAngleV0(psi3V0[0],psi3V0[2],psi3V0[1]);
  if(0){
    std::cout << "(Psi3FullTPC,Psi3PosTPC,Psi3NegTPC) = ("\
      << psi3Tpc[0] << "," << psi3Tpc[2] << "," << psi3Tpc[1] << ")" << std::endl;
    std::cout << "(Psi3FullV0,Psi3V0A,Psi3V0C) = ("\
      << psi3V0[0] << "," << psi3V0[2] << "," << psi3V0[1] << ")" << std::endl;
  }

	//fill histos for q3 spline calibration
  fQ3VecHandler->GetQnVecTPC(q3VecTpcM, q3VecTpcP, q3VecTpcN);
  fQ3VecHandler->GetQnVecV0(q3VecV0M, q3VecV0A, q3VecV0C);
  fQ3VecHandler->GetqnTPC(q3Tpc[0],q3Tpc[2],q3Tpc[1]);
  fQ3VecHandler->GetqnV0(q3V0[0],q3V0[2],q3V0[1]);
  if(0){
    std::cout << "(q3FullTPC,q3PosTPC,q3NegTPC) = ("\
      << q3Tpc[0] << "," << q3Tpc[2] << "," << q3Tpc[1] << ")" << std::endl;
    std::cout << "(q3FullV0,q3V0A,q3V0C) = ("\
      << q3V0[0] << "," << q3V0[2] << "," << q3V0[1] << ")" << std::endl;
  }
}

Bool_t  AliAnalysisTaskRawJetWithEP::QnGainCalibration(){
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
  
  // AliOADBContainer* cont = (AliOADBContainer*) fOADBFile->Get("hMultV0BefCorPfpx");
  // TH1D* fHistMultV0 = ((TH1D*) cont->GetObject(fRunNumber));

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

Bool_t  AliAnalysisTaskRawJetWithEP::QnRecenteringCalibration(){
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

// //_____________________________________________________________________________
// TH1F* AliAnalysisTaskRawJetWithEP::GetResoFromOutputFile(detectorType det, Int_t h, TArrayD* cen)
// {
//     if(!fOutputList) {
//         printf(" > Please add fOutputList first < \n");
//         return 0x0;
//     }
//     TH1F* r(0x0);
//     (cen) ? r = new TH1F("R", "R", cen->GetSize()-1, cen->GetArray()) : r = new TH1F("R", "R", 10, 0, 10);
//     if(!cen) r->GetXaxis()->SetTitle("number of centrality bin");
//     r->GetYaxis()->SetTitle(Form("Resolution #Psi_{%i}", h));
//     for(Int_t i(0); i < 10; i++) {
//         TProfile* temp((TProfile*)fOutputList->FindObject(Form("fProfV%iResolution_%i", h, i)));
//         if(!temp) break;
//         Double_t a(temp->GetBinContent(3)); //cos(2[psi_V0A - psi_V0C])
//         Double_t b(temp->GetBinContent(5)); //cos(2[psi_TPC - psi_V0C])
//         Double_t c(temp->GetBinContent(7)); //cos(2[psi_TPC - psi_V0A])

//         Double_t d(temp->GetBinContent(9));  //cos(2[psi_V0M - psi_TPCnega])
//         Double_t e(temp->GetBinContent(10)); //cos(2[psi_V0M - psi_TPCposi])
//         Double_t f(temp->GetBinContent(11)); //cos(2[psi_TPCnega - psi_TPCposi])

//         Double_t _a(temp->GetBinError(3)), _b(temp->GetBinError(5)), _c(temp->GetBinError(7));
//         Double_t _d(temp->GetBinError(9)), _e(temp->GetBinError(10)), _f(temp->GetBinError(11));
//         Double_t error(0);
//         if(a <= 0 || b <= 0 || c <= 0 || d <= 0 || e <= 0 || f <= 0) continue;
//         switch (det) {
//             case kVZEROA : {
//                 r->SetBinContent(1+i, TMath::Sqrt((a*b)/c));
//                 if(i==0) r->SetNameTitle("VZEROA resolution", "VZEROA resolution");
//                 error = TMath::Power((2.*a*TMath::Sqrt((a*b)/c))/3.,2.)*_a*_a+TMath::Power((2.*b*TMath::Sqrt((a*b)/c))/3.,2.)*_b*_b+TMath::Power(2.*c*TMath::Sqrt((a*b)/c),2.)*_c*_c;
//                 if(error > 0.) error = TMath::Sqrt(error);
//                 r->SetBinError(1+i, error);
//             } break;
//             case kVZEROC : {
//                 r->SetBinContent(1+i, TMath::Sqrt((a*c)/b));
//                 error = TMath::Power((2.*a*TMath::Sqrt((a*c)/b))/3.,2.)*_a*_a+TMath::Power((2.*b*TMath::Sqrt((a*c)/b)),2.)*_b*_b+TMath::Power(2.*c*TMath::Sqrt((a*c)/b)/3.,2.)*_c*_c;
//                 if(error > 0.) error = TMath::Sqrt(error);
//                 if(i==0) r->SetNameTitle("VZEROC resolution", "VZEROC resolution");
//                 r->SetBinError(1+i, error);
//             } break;
//             case kTPC : {
//                 r->SetBinContent(1+i, TMath::Sqrt((b*c)/a));
//                 if(i==0) r->SetNameTitle("TPC resolution", "TPC resolution");
//                 r->SetBinError(1+i, TMath::Sqrt(_a*_a+_b*_b+_c*_c));
//             } break;
//             case kVZEROComb : {
//                 r->SetBinContent(1+i, TMath::Sqrt((d*e)/f));
//                 if(i==0) r->SetNameTitle("VZEROComb resolution", "VZEROComb resolution");
//                 r->SetBinError(1+i, TMath::Sqrt(_d*_d+_e*_e+_f*_f));
//             } break;
//             default : break;
//         }
//     }
//     return r;
// }

//_____________________________________________________________________________
Double_t AliAnalysisTaskRawJetWithEP::CalculateEventPlaneChi(Double_t res)
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
void AliAnalysisTaskRawJetWithEP::BkgFitEvaluation(TH1F* hBkgTracks, TF1* fFitModulation)
{
  // the quality of the fit is evaluated from 1 - the cdf of the chi square distribution
  // three methods are available, all with their drawbacks. 
  // all are stored, one is selected to do the cut
  Int_t numOfFreePara = 2; //v2, v3
  Int_t NDF = 1;
  NDF = (Int_t)fFitModulation->GetXaxis()->GetNbins() - numOfFreePara;
  if(NDF == 0 || (float)NDF <= 0.) return;
  
  Double_t ChiSqr = 999.;
  Double_t CDF = 1.;
  Double_t CDFROOT = 1.;
  ChiSqr = ChiSquare(*hBkgTracks, fFitModulation);
  // CDF = 1. - ChiSquareCDF(NDF, ChiSqr);
  // CDFROOT = 1.-ChiSquareCDF(NDF, fFitModulation->GetChisquare());
  
  CDF = 1. - ChiSquarePDF(NDF, ChiSqr);  
  CDFROOT = 1.-ChiSquarePDF(NDF, fFitModulation->GetChisquare());
  
  // std::cout << "(ChiSqr, ROOTChi, CDF, CDFROOT) = (" << ChiSqr << ", " << fFitModulation->GetChisquare() << ", " << CDF << ", " << CDFROOT << ")" << std::endl;

  
  TString histName;
  TString groupName;
  groupName="BackgroundFit";
  
  // == v2 and v3 combind fit ==
  histName = TString::Format("%s/hPvalueCDF_lRhoCombinFit", groupName.Data());
  fHistManager.FillTH1(histName, CDF);
  histName = TString::Format("%s/hPvalueCDFCent_lRhoCombinFit", groupName.Data());
  fHistManager.FillTH2(histName, fCent, CDF);
  histName = TString::Format("%s/hChi2Cent_lRhoCombinFit", groupName.Data());
  fHistManager.FillTH2(histName, fCent, ChiSqr/((float)NDF));
  histName = TString::Format("%s/hPChi2_lRhoCombinFit", groupName.Data());
  fHistManager.FillTH2(histName, CDF, ChiSqr/((float)NDF));

  histName = TString::Format("%s/hPvalueCDFROOT_lRhoCombinFit", groupName.Data());
  fHistManager.FillTH1(histName, CDFROOT);
  histName = TString::Format("%s/hPvalueCDFROOTCent_lRhoCombinFit", groupName.Data());
  fHistManager.FillTH2(histName, fCent, CDFROOT);
  histName = TString::Format("%s/hChi2ROOTCent_lRhoCombinFit", groupName.Data());
  fHistManager.FillTH2(histName, fCent, ChiSqr/((float)NDF));
  histName = TString::Format("%s/hPChi2ROOT_lRhoCombinFit", groupName.Data());
  fHistManager.FillTH2(histName, CDFROOT, ChiSqr/((float)NDF));

  // std::cout <<"Comb Fit (ChiSqr, CDF) = ("<< ChiSqr/((float)NDF) << ", "<< CDF <<")"<< std::endl; 
  
  // == v2 fit ==
  TF1* tempV2Fit = new TF1("tempRhoFitV2", "[0]*(1.+2.*([1]*TMath::Cos(2.*(x-[2]))))", \
    0.0, TMath::TwoPi());
  tempV2Fit->SetParameter(0, fRho->GetVal());
  tempV2Fit->FixParameter(2, psi2V0[0]);
  hBkgTracks->Fit(tempV2Fit, "N0Q");
  numOfFreePara = 1; //v2
  NDF = tempV2Fit->GetXaxis()->GetNbins() - numOfFreePara;
  if(NDF == 0 || (float)NDF <= 0.) return;
  ChiSqr = 999.;
  ChiSqr = ChiSquare(*hBkgTracks, tempV2Fit);
  // CDF = 1. - ChiSquareCDF(NDF, ChiSqr);
  CDF = 1. - ChiSquarePDF(NDF, ChiSqr);
  // CDFROOT = 1. - ChiSquareCDF(NDF, tempV2Fit->GetChisquare());
  CDFROOT = 1. - ChiSquarePDF(NDF, tempV2Fit->GetChisquare());

  histName = TString::Format("%s/hPvalueCDF_lRhoV2Fit", groupName.Data());
  fHistManager.FillTH1(histName, CDF);
  histName = TString::Format("%s/hPvalueCDFCent_lRhoV2Fit", groupName.Data());
  fHistManager.FillTH2(histName, fCent, CDF);
  histName = TString::Format("%s/hChi2Cent_lRhoV2Fit", groupName.Data());
  fHistManager.FillTH2(histName, fCent, ChiSqr/((float)NDF));
  histName = TString::Format("%s/hPChi2_lRhoV2Fit", groupName.Data());
  fHistManager.FillTH2(histName, CDF, ChiSqr/((float)NDF));

  histName = TString::Format("%s/hPvalueCDFROOT_lRhoV2Fit", groupName.Data());
  fHistManager.FillTH1(histName, CDFROOT);
  histName = TString::Format("%s/hPvalueCDFROOTCent_lRhoV2Fit", groupName.Data());
  fHistManager.FillTH2(histName, fCent, CDFROOT);
  histName = TString::Format("%s/hChi2ROOTCent_lRhoV2Fit", groupName.Data());
  fHistManager.FillTH2(histName, fCent, ChiSqr/((float)NDF));
  histName = TString::Format("%s/hPChi2ROOT_lRhoV2Fit", groupName.Data());
  fHistManager.FillTH2(histName, CDFROOT, ChiSqr/((float)NDF));
  // std::cout << "v2Fit (ChiSqr, CDF) = ("<< ChiSqr/((float)NDF) << ", "<< CDF <<")"<< std::endl;
  delete tempV2Fit;

  
  // == global rho fit ==
  TF1* tempGlobalFit = new TF1("tempGlobalRhoFit", "[0]", 0.0, TMath::TwoPi());
  tempGlobalFit->FixParameter(0, fRho->GetVal());
  tempGlobalFit->SetParameter(0, fRho->GetVal());
  tempGlobalFit->FixParameter(2, psi2V0[0]);
  hBkgTracks->Fit(tempGlobalFit, "N0Q");
  numOfFreePara = 0; //v2
  NDF = tempGlobalFit->GetXaxis()->GetNbins() - numOfFreePara;
  if(NDF == 0 || (float)NDF <= 0.) return;
  ChiSqr = 999.;
  ChiSqr = ChiSquare(*hBkgTracks, tempGlobalFit);
  // CDF = 1. - ChiSquareCDF(NDF, ChiSqr);
  // CDFROOT = 1.-ChiSquareCDF(NDF, tempGlobalFit->GetChisquare());
  CDF = 1. - ChiSquarePDF(NDF, ChiSqr);
  CDFROOT = 1.-ChiSquarePDF(NDF, tempGlobalFit->GetChisquare());

  histName = TString::Format("%s/hPvalueCDF_gRhoFit", groupName.Data());
  fHistManager.FillTH1(histName, CDF);
  histName = TString::Format("%s/hPvalueCDFCent_gRhoFit", groupName.Data());
  fHistManager.FillTH2(histName, fCent, CDF);
  histName = TString::Format("%s/hChi2Cent_gRhoFit", groupName.Data());
  fHistManager.FillTH2(histName, fCent, ChiSqr/((float)NDF));
  histName = TString::Format("%s/hPChi2_gRhoFit", groupName.Data());
  fHistManager.FillTH2(histName, CDF, ChiSqr/((float)NDF));

  histName = TString::Format("%s/hPvalueCDFROOT_gRhoFit", groupName.Data());
  fHistManager.FillTH1(histName, CDFROOT);
  histName = TString::Format("%s/hPvalueCDFROOTCent_gRhoFit", groupName.Data());
  fHistManager.FillTH2(histName, fCent, CDFROOT);
  histName = TString::Format("%s/hChi2ROOTCent_gRhoFit", groupName.Data());
  fHistManager.FillTH2(histName, fCent, ChiSqr/((float)NDF));
  histName = TString::Format("%s/hPChi2ROOT_gRhoFit", groupName.Data());
  fHistManager.FillTH2(histName, CDFROOT, ChiSqr/((float)NDF));
  
  // std::cout << "NoFit (ChiSqr, CDF) = (" << ChiSqr/((float)NDF) << ", "<< CDF <<")"<< std::endl;
  delete tempGlobalFit;
  
}

//_____________________________________________________________________________
void AliAnalysisTaskRawJetWithEP::CalcRandomCone(Double_t &pt, Double_t &eta, Double_t &phi,
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

AliEmcalJet* AliAnalysisTaskRawJetWithEP::GetLeadingJet(AliLocalRhoParameter* localRho) {
    // return pointer to the highest pt jet (before background subtraction) within acceptance
    // only rudimentary cuts are applied on this level, hence the implementation outside of
    // the framework
    Int_t iJets(fJets->GetEntriesFast());
    Double_t pt(0);
    AliEmcalJet* leadingJet(0x0);
    if(!localRho) {
        for(Int_t i(0); i < iJets; i++) {
            AliEmcalJet* jet = static_cast<AliEmcalJet*>(fJets->At(i));
            // if(!PassesSimpleCuts(jet)) continue;
            if(jet->Pt() > pt) {
              leadingJet = jet;
              pt = leadingJet->Pt();
            }
        }
        return leadingJet;
    } else {
        // return leading jet after background subtraction
        Double_t rho(0);
        for(Int_t i(0); i < iJets; i++) {
            AliEmcalJet* jet = static_cast<AliEmcalJet*>(fJets->At(i));
            // if(!PassesSimpleCuts(jet)) continue;
            rho = localRho->GetLocalVal(jet->Phi(), GetJetContainer()->GetJetRadius(), localRho->GetVal());
            // if(fUse2DIntegration) rho = localRho->GetLocalValInEtaPhi(jet->Phi(), GetJetContainer()->GetJetRadius(), localRho->GetVal());
            if((jet->Pt()-jet->Area()*rho) > pt) {
              leadingJet = jet;
              pt = (leadingJet->Pt()-jet->Area()*rho);
            }
        }
        return leadingJet;
    }
    return 0x0;
}


// void AliAnalysisTaskRawJetWithEP::SetupPileUpRemovalFunctions(){
  
//   ////==========> LHC18q/r PileUp Removal Functions: ---- Do not Remove them !!! -----
//   Double_t parV0[8] = {43.8011, 0.822574, 8.49794e-02, 1.34217e+02, 7.09023e+00, 4.99720e-02, -4.99051e-04, 1.55864e-06};
//   fV0CutPU  = new TF1("fV0CutPU", "[0]+[1]*x - 6.*[2]*([3] + [4]*sqrt(x) + [5]*x + [6]*x*sqrt(x) + [7]*x*x)", 0, 100000);
//   fV0CutPU->SetParameters(parV0);
  
//   fSPDCutPU = new TF1("fSPDCutPU", "400. + 4.*x", 0, 10000);

//   Double_t parFB32[8] = {2093.36, -66.425, 0.728932, -0.0027611, 1.01801e+02, -5.23083e+00, -1.03792e+00, 5.70399e-03};
//   fMultCutPU = new TF1("fMultCutPU", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 6.*([4]+[5]*sqrt(x)+[6]*x+[7]*x*x)", 0, 90);
//   fMultCutPU->SetParameters(parFB32);
  
//   Double_t parV0CL0[6] = {0.320462, 0.961793, 1.02278, 0.0330054, -0.000719631, 6.90312e-06};
//   fCenCutLowPU  = new TF1("fCenCutLowPU", "[0]+[1]*x - 6.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)",  0, 100);
//   fCenCutLowPU->SetParameters(parV0CL0);
  
//   fCenCutHighPU = new TF1("fCenCutHighPU", "[0]+[1]*x + 5.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
//   fCenCutHighPU->SetParameters(parV0CL0);
//   //--------------------------------------------------------------------------------------

// }


Bool_t AliAnalysisTaskRawJetWithEP::CheckEventIsPileUp2018(){
  
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
  
  
  // if(centrCL0 < fCenCutLowPU->Eval(centrV0M)) BisPileup=kTRUE;
  // if(centrCL0 > fCenCutHighPU->Eval(centrV0M)) BisPileup=kTRUE;
  // if(Float_t(nITSCls) > fSPDCutPU->Eval(nITSTrkls)) BisPileup=kTRUE;
  // if(multV0On < fV0CutPU->Eval(multV0Tot)) BisPileup=kTRUE;
  // if(Float_t(multTrk) < fMultCutPU->Eval(centrV0M)) BisPileup=kTRUE;
  // if(((AliAODHeader*)fAOD->GetHeader())->GetRefMultiplicityComb08() < 0) BisPileup=kTRUE;
  // if(fAOD->IsIncompleteDAQ()) BisPileup=kTRUE;
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

/**
 * This function is called once at the end of the analysis.
 */
void AliAnalysisTaskRawJetWithEP::Terminate(Option_t *) 
{
}


/**
 * This function adds the task to the analysis manager. Often, this function is called
 * by an AddTask C macro. However, by compiling the code, it ensures that we do not
 * have to deal with difficulties caused by CINT.
 */
AliAnalysisTaskRawJetWithEP * AliAnalysisTaskRawJetWithEP::AddTaskRawJetWithEP(
  const char *ntracks, const char *nclusters, const char* ncells, const char *suffix)
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskRawJetWithEP", "No analysis manager to connect to.");
    return 0;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler)
  {
    ::Error("AddTaskRawJetWithEP", "This task requires an input event handler");
    return 0;
  }
  
  enum EDataType_t {kUnknown, kESD, kAOD};


  EDataType_t dataType = kAOD;
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------
  TString trackName(ntracks);
  if (trackName == "usedefault") trackName = "tracks";

  TString name("AliAnalysisTaskRawJetWithEP");
  if (strcmp(suffix,"") != 0) {
    name += "_";
    name += suffix;
  }

  AliAnalysisTaskRawJetWithEP* rawJetTask = new AliAnalysisTaskRawJetWithEP(name);
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













