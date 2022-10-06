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
#include "AliAnalysisTaskEPCalibForJet.h"

#include "AliDataFile.h"

class AliAnalysisTaskEPCalibForJet;

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskEPCalibForJet);
/// \endcond

//Default constructor. Needed by ROOT I/O
AliAnalysisTaskEPCalibForJet::AliAnalysisTaskEPCalibForJet() :
  AliAnalysisTaskEmcalJet(),
  fHistManager(),
  fEventCuts(),
  fYAMLConfig(),
  fUseRunList(),
  fAOD(nullptr),
  fOutputList(nullptr),
  fCalibRefFileName(""),
  fCalibRefFile(nullptr),
  fCalibRefObjList(nullptr),
  fCalibV0Ref(nullptr),
  fQ2VecHandler(nullptr),
  fQ3VecHandler(nullptr),
  fPileupCut(kFALSE),
  fTPCQnMeasure(kFALSE),
  fPileupCutQA(kFALSE),
  fCalibQA(kFALSE),
  fGainCalibQA(kFALSE),
  fReCentCalibQA(kFALSE),
  fEPQA(kFALSE),
  fTrackQA(kFALSE),
  fBkgQA(kFALSE),
  fV0Combin(kFALSE),
  fCalibType(0),
  fQnVCalibType(kOrig),
  fNormMethod(0),
  fV2ResoV0(0.), fV3ResoV0(0.),
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
  CheckRunNum(0), fQaEventNum(-1)
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

// Standard constructor. Should be used by the user. @param[in] name Name of the task
AliAnalysisTaskEPCalibForJet::AliAnalysisTaskEPCalibForJet(const char *name):
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fHistManager(name),
  fEventCuts(),
  fYAMLConfig(),
  fUseRunList(),
  fAOD(nullptr),
  fOutputList(nullptr),
  fCalibRefFileName(""),
  fCalibRefFile(nullptr),
  fCalibRefObjList(nullptr),
  fCalibV0Ref(nullptr),
  fQ2VecHandler(nullptr),
  fQ3VecHandler(nullptr),
  fPileupCut(kFALSE),
  fTPCQnMeasure(kFALSE),
  fPileupCutQA(kFALSE),
  fCalibQA(kFALSE),
  fGainCalibQA(kFALSE),
  fReCentCalibQA(kFALSE),
  fEPQA(kFALSE),
  fTrackQA(kFALSE),
  fBkgQA(kFALSE),
  fV0Combin(kFALSE),
  fCalibType(0),
  fQnVCalibType(kOrig),
  fNormMethod(0),
  fV2ResoV0(0.), fV3ResoV0(0.),
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
  CheckRunNum(0), fQaEventNum(-1)
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

// Destructor
AliAnalysisTaskEPCalibForJet::~AliAnalysisTaskEPCalibForJet()
{
  if(fCalibRefFile) {fCalibRefFile->Close(); fCalibRefFile = 0x0;}
  if(fCalibV0Ref)   {delete fCalibV0Ref;   fCalibV0Ref = 0x0;}
  if(fCalibRefObjList) {delete fCalibRefObjList; fCalibRefObjList = 0x0;}
  if(fQ2VecHandler) {delete fQ2VecHandler; fQ2VecHandler = 0x0;}
  if(fQ3VecHandler) {delete fQ3VecHandler; fQ3VecHandler = 0x0;}
  if(fFitModulation) {delete fFitModulation; fFitModulation = 0x0;}
  if(hBkgTracks) {delete hBkgTracks; hBkgTracks = 0x0;}
}


void AliAnalysisTaskEPCalibForJet::SetRunList(bool removeDummyTask)
{
  
  fYAMLConfig.AddConfiguration(fRunListFileName, "runlist");
  fYAMLConfig.Initialize();
  fYAMLConfig.GetProperty("runlist", fUseRunList);
}


/**
 * Performing run-independent initialization.
 * Here the histograms should be instantiated.
 */
void AliAnalysisTaskEPCalibForJet::UserCreateOutputObjects()
{
  
  fOutputList = new TList();
  fOutputList->SetOwner(true);

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();
  
  fEventCuts.AddQAplotsToList(fOutput);
  fEventCuts.OverrideAutomaticTriggerSelection(AliVEvent::kAnyINT | AliVEvent::kCentral | AliVEvent::kSemiCentral);
  
  // == s == Set Out put Hist grams  ###########################################
  if(fPileupCutQA) AllocatePileupCutHistograms();
  if(fGainCalibQA) AllocateGainCalibHistograms();
  if(fReCentCalibQA) AllocateReCentCalibHistograms();
  if(fEPQA) AllocateEventPlaneHistograms();
  if(fTrackQA) AllocateTrackHistograms();
  if(fBkgQA) AllocateBkgHistograms();
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

void AliAnalysisTaskEPCalibForJet::AllocatePileupCutHistograms(){
  TString histName;
  TString histtitle;
  TString groupName;
  groupName="PileupCutQA";

  Int_t gMaxGlobalmult  = 4000;
  Int_t gMaxTPCcorrmult = 5000;
  Int_t gMaxESDtracks   = 20000;

  fHistManager.CreateHistoGroup(groupName);

  histName  = TString::Format("%s/fHistCentCL0VsV0MBefore", groupName.Data());
  histtitle = TString::Format("%s;Cent(V0M); Cent(CL0)", histName.Data());
  fHistManager.CreateTH2(histName, histtitle, 100,0,100,100,0,100);
  histName  = TString::Format("%s/fHistTPCVsESDTrkBefore", groupName.Data());
  histtitle = TString::Format("%s;TPC1; ESD trk", histName.Data());
  fHistManager.CreateTH2(histName, histtitle, 100,0,5000,200,0,20000);
  histName  = TString::Format("%s/fHistTPConlyVsCL1Before", groupName.Data());
  histtitle = TString::Format("%s;Cent(CL1); TPC(FB128)", histName.Data());
  fHistManager.CreateTH2(histName, histtitle, 100,0,100,250,0,gMaxTPCcorrmult);
  histName  = TString::Format("%s/fHistTPConlyVsV0MBefore", groupName.Data());
  histtitle = TString::Format("%s;Cent(V0M); TPC(FB128)", histName.Data());
  fHistManager.CreateTH2(histName, histtitle, 100,0,100,250,0,gMaxTPCcorrmult);

  histName  = TString::Format("%s/fHistCentCL0VsV0MAfter", groupName.Data());
  histtitle = TString::Format("%s;Cent(V0M); Cent(CL0)", histName.Data());
  fHistManager.CreateTH2(histName, histtitle, 100,0,100,100,0,100);
  histName  = TString::Format("%s/fHistTPCVsESDTrkAfter", groupName.Data());
  histtitle = TString::Format("%s;TPC1; ESD trk", histName.Data());
  fHistManager.CreateTH2(histName, histtitle, 100,0,5000,200,0,20000);
  histName  = TString::Format("%s/fHistTPConlyVsCL1After", groupName.Data());
  histtitle = TString::Format("%s;Cent(CL1); TPC(FB128)", histName.Data());
  fHistManager.CreateTH2(histName, histtitle, 100,0,100,250,0,gMaxTPCcorrmult);
  histName  = TString::Format("%s/fHistTPConlyVsV0MAfter", groupName.Data());
  histtitle = TString::Format("%s;Cent(V0M); TPC(FB128)", histName.Data());
  fHistManager.CreateTH2(histName, histtitle, 100,0,100,250,0,gMaxTPCcorrmult);

}

void AliAnalysisTaskEPCalibForJet::AllocateGainCalibHistograms(){
  TString histName;
  TString histtitle;
  TString groupName;
  groupName="GainCalib";

  fHistManager.CreateHistoGroup(groupName);
  THashList *parent = static_cast<THashList *>(fHistManager.FindObject(groupName.Data()));
  
  for(Int_t eventNumBin = 0; eventNumBin < fUseRunList.size(); eventNumBin++){
    Int_t runEventNum = stoi(fUseRunList.at(eventNumBin));
    histName = TString::Format("hAvgV0ChannelsvsVz_%d", runEventNum);
    histtitle = histName;
    
    TProfile2D *tempHist = new TProfile2D(histName, histtitle, 10,-10,10,64,0,64);
    parent->Add(tempHist);
    
    histName = TString::Format("GainCalib/hTPCPosiTrkVzPhiEta_%d", runEventNum);
    fHistManager.CreateTH3(histName, histtitle,10,-10,10,50,0,6.283185,16,-0.8,0.8);
    histName = TString::Format("GainCalib/hTPCNegaTrkVzPhiEta_%d", runEventNum);
    fHistManager.CreateTH3(histName, histtitle,10,-10,10,50,0,6.283185,16,-0.8,0.8);
    
  }

}

void AliAnalysisTaskEPCalibForJet::AllocateReCentCalibHistograms(){
  TString histName;
  TString histtitle;
  TString groupName;
  groupName="ReCentCalib";

  fHistManager.CreateHistoGroup(groupName);
  THashList *parent = static_cast<THashList *>(fHistManager.FindObject(groupName.Data()));
  
  
  for(Int_t eventNumBin = 0; eventNumBin < fUseRunList.size(); eventNumBin++){
    Int_t runEventNum = stoi(fUseRunList.at(eventNumBin));

    // == s ==  V0 recentering Measurement  ################################################
    //  Q2
    histName = TString::Format("%s/hAvgQ2XvsCentV0CBef_%d", groupName.Data(), runEventNum);
    histtitle = histtitle = TString::Format("%s;Centrality;<Q2_{x}>^{V0C}", histName.Data());
    fHistManager.CreateTProfile(histName, histtitle,90,0,90);
    histName = TString::Format("%s/hAvgQ2YvsCentV0CBef_%d", groupName.Data(), runEventNum);
    histtitle = histtitle = TString::Format("%s;Centrality;<Q2_{y}>^{V0C}", histName.Data());
    fHistManager.CreateTProfile(histName, histtitle,90,0,90);
    histName = TString::Format("%s/hAvgQ2XvsCentV0ABef_%d", groupName.Data(), runEventNum);
    histtitle = histtitle = TString::Format("%s;Centrality;<Q2_{x}>^{V0A}", histName.Data());
    fHistManager.CreateTProfile(histName, histtitle,90,0,90);
    histName = TString::Format("%s/hAvgQ2YvsCentV0ABef_%d", groupName.Data(), runEventNum);
    histtitle = histtitle = TString::Format("%s;Centrality;<Q2_{y}>^{V0A}", histName.Data());
    fHistManager.CreateTProfile(histName, histtitle,90,0,90);
    
    //  Q3
    histName = TString::Format("%s/hAvgQ3XvsCentV0CBef_%d", groupName.Data(), runEventNum);
    histtitle = histtitle = TString::Format("%s;Centrality;<Q_3{x}>^{V0C}", histName.Data());
    fHistManager.CreateTProfile(histName, histtitle,90,0,90);
    histName = TString::Format("%s/hAvgQ3YvsCentV0CBef_%d", groupName.Data(), runEventNum);
    histtitle = histtitle = TString::Format("%s;Centrality;<Q3_{y}>^{V0C}", histName.Data());
    fHistManager.CreateTProfile(histName, histtitle,90,0,90);
    histName = TString::Format("%s/hAvgQ3XvsCentV0ABef_%d", groupName.Data(), runEventNum);
    histtitle = histtitle = TString::Format("%s;Centrality;<Q3_{x}>^{V0A}", histName.Data());
    fHistManager.CreateTProfile(histName, histtitle,90,0,90);
    histName = TString::Format("%s/hAvgQ3YvsCentV0ABef_%d", groupName.Data(), runEventNum);
    histtitle = histtitle = TString::Format("%s;Centrality;<Q3_{y}>^{V0A}", histName.Data());
    fHistManager.CreateTProfile(histName, histtitle,90,0,90);

    // Calib After
    //  Q2
    histName = TString::Format("%s/hAvgQ2XvsCentV0CAft_%d", groupName.Data(), runEventNum);
    histtitle = histtitle = TString::Format("%s;Centrality;<Q2_{x}>^{V0C}", histName.Data());
    fHistManager.CreateTProfile(histName, histtitle,90,0,90);
    histName = TString::Format("%s/hAvgQ2YvsCentV0CAft_%d", groupName.Data(), runEventNum);
    histtitle = histtitle = TString::Format("%s;Centrality;<Q2_{y}>^{V0C}", histName.Data());
    fHistManager.CreateTProfile(histName, histtitle,90,0,90);
    histName = TString::Format("%s/hAvgQ2XvsCentV0AAft_%d", groupName.Data(), runEventNum);
    histtitle = histtitle = TString::Format("%s;Centrality;<Q2_{x}>^{V0A}", histName.Data());
    fHistManager.CreateTProfile(histName, histtitle,90,0,90);
    histName = TString::Format("%s/hAvgQ2YvsCentV0AAft_%d", groupName.Data(), runEventNum);
    histtitle = histtitle = TString::Format("%s;Centrality;<Q2_{y}>^{V0A}", histName.Data());
    fHistManager.CreateTProfile(histName, histtitle,90,0,90);
    
    //  Q3
    histName = TString::Format("%s/hAvgQ3XvsCentV0CAft_%d", groupName.Data(), runEventNum);
    histtitle = histtitle = TString::Format("%s;Centrality;<Q_3{x}>^{V0C}", histName.Data());
    fHistManager.CreateTProfile(histName, histtitle,90,0,90);
    histName = TString::Format("%s/hAvgQ3YvsCentV0CAft_%d", groupName.Data(), runEventNum);
    histtitle = histtitle = TString::Format("%s;Centrality;<Q3_{y}>^{V0C}", histName.Data());
    fHistManager.CreateTProfile(histName, histtitle,90,0,90);
    histName = TString::Format("%s/hAvgQ3XvsCentV0AAft_%d", groupName.Data(), runEventNum);
    histtitle = histtitle = TString::Format("%s;Centrality;<Q3_{x}>^{V0A}", histName.Data());
    fHistManager.CreateTProfile(histName, histtitle,90,0,90);
    histName = TString::Format("%s/hAvgQ3YvsCentV0AAft_%d", groupName.Data(), runEventNum);
    histtitle = histtitle = TString::Format("%s;Centrality;<Q3_{y}>^{V0A}", histName.Data());
    fHistManager.CreateTProfile(histName, histtitle,90,0,90);
    // == e ==  V0 recentering Measurement  ################################################


    // == s ==  TPC recentering Measurement  ###############################################
    //  Q2
    histName = TString::Format("%s/hAvgQ2XvsCentTPCPBef_%d", groupName.Data(), runEventNum);
    histtitle = histtitle = TString::Format("%s;Centrality;<Q2_{x}>^{TPC_Posi}", histName.Data());
    fHistManager.CreateTProfile(histName, histtitle,90,0,90);
    histName = TString::Format("%s/hAvgQ2YvsCentTPCPBef_%d", groupName.Data(), runEventNum);
    histtitle = histtitle = TString::Format("%s;Centrality;<Q2_{y}>^{TPC_Posi}", histName.Data());
    fHistManager.CreateTProfile(histName, histtitle,90,0,90);
    histName = TString::Format("%s/hAvgQ2XvsCentTPCNBef_%d", groupName.Data(), runEventNum);
    histtitle = histtitle = TString::Format("%s;Centrality;<Q2_{x}>^{TPC_Nega}", histName.Data());
    fHistManager.CreateTProfile(histName, histtitle,90,0,90);
    histName = TString::Format("%s/hAvgQ2YvsCentTPCNBef_%d", groupName.Data(), runEventNum);
    histtitle = histtitle = TString::Format("%s;Centrality;<Q2_{y}>^{TPC_Nega}", histName.Data());
    fHistManager.CreateTProfile(histName, histtitle,90,0,90);

    //  Q3
    histName = TString::Format("%s/hAvgQ3XvsCentTPCPBef_%d", groupName.Data(), runEventNum);
    histtitle = histtitle = TString::Format("%s;Centrality;<Q_3{x}>^{TPC_Posi}", histName.Data());
    fHistManager.CreateTProfile(histName, histtitle,90,0,90);
    histName = TString::Format("%s/hAvgQ3YvsCentTPCPBef_%d", groupName.Data(), runEventNum);
    histtitle = histtitle = TString::Format("%s;Centrality;<Q3_{y}>^{TPC_Posi}", histName.Data());
    fHistManager.CreateTProfile(histName, histtitle,90,0,90);
    histName = TString::Format("%s/hAvgQ3XvsCentTPCNBef_%d", groupName.Data(), runEventNum);
    histtitle = histtitle = TString::Format("%s;Centrality;<Q3_{x}>^{TPC_Nega}", histName.Data());
    fHistManager.CreateTProfile(histName, histtitle,90,0,90);
    histName = TString::Format("%s/hAvgQ3YvsCentTPCNBef_%d", groupName.Data(), runEventNum);
    histtitle = histtitle = TString::Format("%s;Centrality;<Q3_{y}>^{TPC_Nega}", histName.Data());
    fHistManager.CreateTProfile(histName, histtitle,90,0,90);


    //  Q2
    histName = TString::Format("%s/hAvgQ2XvsCentTPCPAft_%d", groupName.Data(), runEventNum);
    histtitle = histtitle = TString::Format("%s;Centrality;<Q2_{x}>^{TPC_Posi}", histName.Data());
    fHistManager.CreateTProfile(histName, histtitle,90,0,90);
    histName = TString::Format("%s/hAvgQ2YvsCentTPCPAft_%d", groupName.Data(), runEventNum);
    histtitle = histtitle = TString::Format("%s;Centrality;<Q2_{y}>^{TPC_Posi}", histName.Data());
    fHistManager.CreateTProfile(histName, histtitle,90,0,90);
    histName = TString::Format("%s/hAvgQ2XvsCentTPCNAft_%d", groupName.Data(), runEventNum);
    histtitle = histtitle = TString::Format("%s;Centrality;<Q2_{x}>^{TPC_Nega}", histName.Data());
    fHistManager.CreateTProfile(histName, histtitle,90,0,90);
    histName = TString::Format("%s/hAvgQ2YvsCentTPCNAft_%d", groupName.Data(), runEventNum);
    histtitle = histtitle = TString::Format("%s;Centrality;<Q2_{y}>^{TPC_Nega}", histName.Data());
    fHistManager.CreateTProfile(histName, histtitle,90,0,90);

    //  Q3
    histName = TString::Format("%s/hAvgQ3XvsCentTPCPAft_%d", groupName.Data(), runEventNum);
    histtitle = histtitle = TString::Format("%s;Centrality;<Q_3{x}>^{TPC_Posi}", histName.Data());
    fHistManager.CreateTProfile(histName, histtitle,90,0,90);
    histName = TString::Format("%s/hAvgQ3YvsCentTPCPAft_%d", groupName.Data(), runEventNum);
    histtitle = histtitle = TString::Format("%s;Centrality;<Q3_{y}>^{TPC_Posi}", histName.Data());
    fHistManager.CreateTProfile(histName, histtitle,90,0,90);
    histName = TString::Format("%s/hAvgQ3XvsCentTPCNAft_%d", groupName.Data(), runEventNum);
    histtitle = histtitle = TString::Format("%s;Centrality;<Q3_{x}>^{TPC_Nega}", histName.Data());
    fHistManager.CreateTProfile(histName, histtitle,90,0,90);
    histName = TString::Format("%s/hAvgQ3YvsCentTPCNAft_%d", groupName.Data(), runEventNum);
    histtitle = histtitle = TString::Format("%s;Centrality;<Q3_{y}>^{TPC_Nega}", histName.Data());
    fHistManager.CreateTProfile(histName, histtitle,90,0,90);

    // == e ==  TPC recentering Measurement  ###############################################
  }

}


void AliAnalysisTaskEPCalibForJet::AllocateEventPlaneHistograms()
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
  fHistManager.CreateTH1(histName, histtitle, 64, 0, 64);
  histName  = TString::Format("%s/hV0CellChRatio", groupName.Data());
  histtitle = TString::Format("%s;cell ch number;Ch ratio", histName.Data());
  fHistManager.CreateTH1(histName, histtitle, 64, 0, 64);

  histName  =TString::Format("%s/hV0CellChVsMultBefEq", groupName.Data());
  histtitle =TString::Format("%s;cell ch number;multiplisity before equalization",histName.Data());
  fHistManager.CreateTH1(histName, histtitle, 64, 0, 64);

  histName  =TString::Format("%s/hV0CellChVsMultAftEq", groupName.Data());
  histtitle =TString::Format("%s;cell ch number;multiplisity before equalization",histName.Data());
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

}

void AliAnalysisTaskEPCalibForJet::AllocateBkgHistograms()
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
  fHistManager.CreateTH2(histName, histtitle, 1000, 0, 1, 100, 0, 5);

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
  fHistManager.CreateTH2(histName, histtitle, 1000, 0, 1, 100, 0, 5);

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
  fHistManager.CreateTH2(histName, histtitle, 1000, 0, 1, 100, 0, 5);

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
  fHistManager.CreateTH2(histName, histtitle, 1000, 0, 1, 100, 0, 5);

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
  fHistManager.CreateTH2(histName, histtitle, 1000, 0, 1, 100, 0, 5);

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
  fHistManager.CreateTH2(histName, histtitle, 1000, 0, 1, 100, 0, 5);
  // == e == cdf and pdf of chisquare distribution #############################


  
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
void AliAnalysisTaskEPCalibForJet::AllocateTrackHistograms()
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


/**
 * This function is executed automatically for the first event.
 * Some extra initialization can be performed here.
 */
void AliAnalysisTaskEPCalibForJet::ExecOnce()
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
Bool_t AliAnalysisTaskEPCalibForJet::Run()
{
  // std::cout << "ChecKuma Run Number === " << CheckRunNum << "==================" << std::endl;
  
  CheckRunNum++;
  if(!fEventCuts.AcceptEvent(InputEvent())) return kFALSE;
  
  if(fPileupCut){
    SetupPileUpRemovalFunctions();
    Bool_t kPileupCutEvent = CheckEventIsPileUp2018();
    if(kPileupCutEvent) return kFALSE;
  }
  if(fGainCalibQA) DoMeasureChGainDiff();
  if(fEPQA) DoEventPlane();
  if(fBkgQA){
    SetModulationRhoFit();
    // std::cout << "Fomula = " << fFitModulation->GetExpFormula() << std::endl;
    //  MeasureTpcEPQA();
    MeasureBkg();
  }
  
  return kTRUE;
}

void AliAnalysisTaskEPCalibForJet::DoMeasureChGainDiff(){
  TString histName;
  TString groupName;
  groupName="GainCalib";
  
  AliAODVZERO* fAodV0 = dynamic_cast<AliAODVZERO*>(fAOD->GetVZEROData());
  AliVEvent *fVevent = dynamic_cast<AliVEvent*>(InputEvent());
  const AliVVertex *pointVtx = fVevent->GetPrimaryVertex();
  Double_t fVtxZ = -999;
  fVtxZ  = pointVtx->GetZ();
  
  Int_t ibinV0 = 0;
  Double_t fMultV0 = 0.;
  
  
  for(int iV0 = 0; iV0 < 64; iV0++) { //0-31 is V0C, 32-63 VOA
    fMultV0 = fAodV0->GetMultiplicity(iV0);
    
    histName = TString::Format("%s/hAvgV0ChannelsvsVz_%d", groupName.Data(), fRunNumber);
    TProfile2D* tempProfile2DHist = dynamic_cast<TProfile2D *>(fHistManager.FindObject(histName));
    if(!tempProfile2DHist){
      Fatal("THistManager::FillTProfile", "Histogram %s not found in parent group %s", histName.Data(), groupName.Data());
    }
    tempProfile2DHist->Fill(fVtxZ, iV0, fMultV0);//?????
    // fHistManager.FillProfile(histName, fVtxZ, iV0, fMultV0);
  }
  
}


void AliAnalysisTaskEPCalibForJet::DoEventPlane(){

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
  QnV0GainCalibration();
  if(0){
    std::cout << "recent calibA (qx,qy) = " \
      << q2VecV0A[0] << "," << q2VecV0A[1] << ")" << std::endl;
    std::cout << "recent calibM (qx,qy) = " \
      << q2VecV0M[0] << "," << q2VecV0M[1] << ")" << std::endl;
  }
  
  if(fTPCQnMeasure) MeasureQnTPC();
  
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


void AliAnalysisTaskEPCalibForJet::SetModulationRhoFit() 
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


void AliAnalysisTaskEPCalibForJet::MeasureBkg(){
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

  fLocalRho->SetLocalRho(fFitModulation);
  // fLocalRho->SetVal(fRho->GetVal());

  BkgFitEvaluation();
}


Double_t AliAnalysisTaskEPCalibForJet::CalcEPReso(Int_t n, \
  Double_t &psiA, Double_t &psiB, Double_t &psiC){
  
  Double_t vnReso = -999.;
  vnReso = TMath::Sqrt((TMath::Abs(TMath::Cos(n*(psiA - psiB))) \
                          * TMath::Abs(TMath::Cos(n*(psiA - psiC)))) \
                        / TMath::Abs(TMath::Cos(n*(psiB - psiC))));

  return vnReso;
}


void AliAnalysisTaskEPCalibForJet::MeasureTpcEPQA(){
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

void  AliAnalysisTaskEPCalibForJet::MeasureQnTPC(){
  TString histName;
  TString groupName;

  Double_t centrV0M=-99.0;
  Double_t centrCL1=-99.0;
  Double_t centrCL0=-99.0;

  
  AliMultSelection* fMultSelection \
    = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");

  if(!fMultSelection) {
    printf("\n\n **WARNING** ::UserExec() AliMultSelection object not found.\n\n");
    exit(111);
  }
  
  
  AliVEvent *fVevent = dynamic_cast<AliVEvent*>(InputEvent());
  const AliVVertex *pointVtx = fVevent->GetPrimaryVertex();
  Double_t fVtxZ = -999.;
  fVtxZ  = pointVtx->GetZ();

  centrV0M = fMultSelection->GetMultiplicityPercentile("V0M");
  centrCL1 = fMultSelection->GetMultiplicityPercentile("CL1");
  centrCL0 = fMultSelection->GetMultiplicityPercentile("CL0");

  AliAODTracklets* aodTrkl = (AliAODTracklets*)fAOD->GetTracklets();
  Int_t nITSTrkls = aodTrkl->GetNumberOfTracklets();

  const Int_t nTracks = fAOD->GetNumberOfTracks();
  Int_t multTrk = 0;
  
  Int_t    fTPCclustMinforEP= 70;      // Fixed for EP calculation
  Int_t    trkChrg=0,trkNClus=0;     //
  Double_t fMinPtCutforEP   = 0.2;    // Fixed for EP calculation
  Double_t fMaxPtCutforEP   = 2.0;    // Fixed for EP calculation
  Double_t fEtaGapPosforEP  = 0.1;    // could be made variable in AddTask Macro.
  Double_t fEtaGapNegforEP  =-0.1;    // could be made variable in AddTask Macro.

  Double_t trkPt=0, trkPhi=0, trkEta=0, trkChi2=0, trkdEdx=0, trkWgt=1.0;  
  Double_t SumQ2xTPCPos = 0., SumQ2yTPCPos = 0., SumQ2xTPCNeg = 0., SumQ2yTPCNeg = 0;
  Double_t SumQ3xTPCPos = 0., SumQ3yTPCPos = 0., SumQ3xTPCNeg = 0., SumQ3yTPCNeg = 0;
  
  Double_t fWgtMultTPCPos=0., fWgtMultTPCNeg=0;
  
  for (Int_t it = 0; it < nTracks; it++) {
    AliAODTrack* aodTrk = (AliAODTrack*)fAOD->GetTrack(it);
    if (!aodTrk) continue;
    
    if (aodTrk->TestFilterBit(32)){
      if((TMath::Abs(aodTrk->Eta()) < 0.8)&&(aodTrk->GetTPCNcls() >= 70)&&(aodTrk->Pt() >= 0.2))
        multTrk++;
    }
    
    /// Now get TPC Q-Vectors:
    // if(!aodTrk->TestFilterBit(fFilterBit))    continue;  //// Only use FB tracks. 
    trkPt    = aodTrk->Pt();
    trkPhi   = aodTrk->Phi();
    trkEta   = aodTrk->Eta();
    trkChrg  = aodTrk->Charge();
    trkChi2  = aodTrk->Chi2perNDF();
    trkNClus = aodTrk->GetTPCNcls();
    // trkdEdx  = aodTrk->GetDetPid()->GetTPCsignal(); //?????????????????????????
    
    
    //Apply track cuts for EP  here:
    if(trkNClus< fTPCclustMinforEP)        continue;
    if((trkEta < -0.8) || (trkEta > 0.8)) continue;    
    if(trkPt < fMinPtCutforEP) continue;
    if(trkPt > fMaxPtCutforEP) continue;
    if(trkChi2 < 0.1)          continue;
    if(trkChi2 > 4.0)          continue;
    // if(trkdEdx < 10)           continue;
    if(!TMath::Abs(trkChrg))   continue;
    
    groupName="GainCalib";
    if(trkChrg > 0){ ///+ve Ch done
      histName = TString::Format("%s/hTPCPosiTrkVzPhiEta_%d", groupName.Data(), fRunNumber);
      fHistManager.FillTH3(histName, fVtxZ,trkPhi,trkEta);
    }else{  //-Ve charge
      histName = TString::Format("%s/hTPCNegaTrkVzPhiEta_%d", groupName.Data(), fRunNumber);
      fHistManager.FillTH3(histName, fVtxZ,trkPhi,trkEta);
    }
    

    Int_t trkID = aodTrk->GetID();
    // trkWgt = GetNUAWeightForTrack(fVertexZEvent,trkPhi,trkEta,trkChrg); //???????
    
    ///Used Pt as weight for Better resolution:
    if(trkEta >= fEtaGapPosforEP){
      SumQ2xTPCPos += trkWgt*TMath::Cos(2*trkPhi);
      SumQ2yTPCPos += trkWgt*TMath::Sin(2*trkPhi);
      SumQ3xTPCPos += trkWgt*TMath::Cos(3*trkPhi);
      SumQ3yTPCPos += trkWgt*TMath::Sin(3*trkPhi);
      
      fWgtMultTPCPos += trkWgt; 
      // vecPosEPTrkID.push_back(trkID);
      // vecPosEPTrkPhi.push_back(trkPhi);
      // vecPosEPTrkNUAWgt.push_back(trkWgt);
    }
    else if(trkEta <= fEtaGapNegforEP){
      SumQ2xTPCNeg += trkWgt*TMath::Cos(2*trkPhi);
      SumQ2yTPCNeg += trkWgt*TMath::Sin(2*trkPhi);
      SumQ3xTPCNeg += trkWgt*TMath::Cos(3*trkPhi);
      SumQ3yTPCNeg += trkWgt*TMath::Sin(3*trkPhi);
      
      fWgtMultTPCNeg += trkWgt;
      // vecNegEPTrkID.push_back(trkID);
      // vecNegEPTrkPhi.push_back(trkPhi);
      // vecNegEPTrkNUAWgt.push_back(trkWgt);
    }
    
  }//AOD track loop
  
    /// Set The q vector values for Event Plane: 
  if(fWgtMultTPCPos<0.1 || fWgtMultTPCNeg<0.1){        /// this means there is not enough tracks in this event!!
    q2VecTpcP[0] = 0;
    q2VecTpcP[1] = 0;
    q2VecTpcN[0] = 0;
    q2VecTpcN[1] = 0;

    q3VecTpcP[0] = 0;
    q3VecTpcP[1] = 0;
    q3VecTpcN[0] = 0;
    q3VecTpcN[1] = 0;
  }
  else{
    q2VecTpcP[0] = SumQ2xTPCPos;
    q2VecTpcP[1] = SumQ2yTPCPos;
    q2VecTpcN[0] = SumQ2xTPCNeg;
    q2VecTpcN[1] = SumQ2yTPCNeg;

    q3VecTpcP[0] = SumQ3xTPCPos;
    q3VecTpcP[1] = SumQ3yTPCPos;
    q3VecTpcN[0] = SumQ3xTPCNeg;
    q3VecTpcN[1] = SumQ3yTPCNeg;
  }

  // std::cout << "SumQ2xTPCPos:SumQ2xTPCNeg = " << SumQ2xTPCPos <<" : "<< SumQ2xTPCNeg << std::endl;

  if(fReCentCalibQA){
    groupName="ReCentCalib";

    histName = TString::Format("%s/hAvgQ2XvsCentTPCPBef_%d", groupName.Data(), fRunNumber);
    fHistManager.FillProfile(histName, fCent, q2VecTpcP[0]);
    histName = TString::Format("%s/hAvgQ2YvsCentTPCPBef_%d", groupName.Data(), fRunNumber);
    fHistManager.FillProfile(histName, fCent, q2VecTpcP[1]);
    histName = TString::Format("%s/hAvgQ2XvsCentTPCNBef_%d", groupName.Data(), fRunNumber);
    fHistManager.FillProfile(histName, fCent, q2VecTpcN[0]);
    histName = TString::Format("%s/hAvgQ2YvsCentTPCNBef_%d", groupName.Data(), fRunNumber);
    fHistManager.FillProfile(histName, fCent, q2VecTpcN[1]);

    histName = TString::Format("%s/hAvgQ3XvsCentTPCPBef_%d", groupName.Data(), fRunNumber);
    fHistManager.FillProfile(histName, fCent, q3VecTpcP[0]);
    histName = TString::Format("%s/hAvgQ3YvsCentTPCPBef_%d", groupName.Data(), fRunNumber);
    fHistManager.FillProfile(histName, fCent, q3VecTpcP[1]);
    histName = TString::Format("%s/hAvgQ3XvsCentTPCNBef_%d", groupName.Data(), fRunNumber);
    fHistManager.FillProfile(histName, fCent, q3VecTpcN[0]);
    histName = TString::Format("%s/hAvgQ3YvsCentTPCNBef_%d", groupName.Data(), fRunNumber);
    fHistManager.FillProfile(histName, fCent, q3VecTpcN[1]);

  }
  
}


Bool_t  AliAnalysisTaskEPCalibForJet::QnV0GainCalibration(){
  TString histName;
  TString groupName;
  groupName="EventPlane";
  
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
  Double_t fMultV0 = 0.;
  
  for(int iV0 = 0; iV0 < 64; iV0++) { //0-31 is V0C, 32-63 VOA
    fMultV0 = fAodV0->GetMultiplicity(iV0);

    Double_t fV0chGain = 1.;
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


Bool_t  AliAnalysisTaskEPCalibForJet::QnRecenteringCalibration(){
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
  // if(fHCorrQ2xV0C && fHCorrQ2yV0C && fHCorrQ2xV0A && fHCorrQ2yV0A){
  //   printf(" ===========> Info:: V0A,V0C <Q> Found for Run %d \n ",fRunNumber);
  // }


  if(fReCentCalibQA){
    TString histName;
    TString groupName;
    groupName="ReCentCalib";

    histName = TString::Format("%s/hAvgQ2XvsCentV0CBef_%d", groupName.Data(), fRunNumber);
    fHistManager.FillProfile(histName, fCent, q2VecV0C[0]);
    histName = TString::Format("%s/hAvgQ2YvsCentV0CBef_%d", groupName.Data(), fRunNumber);
    fHistManager.FillProfile(histName, fCent, q2VecV0C[1]);
    histName = TString::Format("%s/hAvgQ2XvsCentV0ABef_%d", groupName.Data(), fRunNumber);
    fHistManager.FillProfile(histName, fCent, q2VecV0A[0]);
    histName = TString::Format("%s/hAvgQ2YvsCentV0ABef_%d", groupName.Data(), fRunNumber);
    fHistManager.FillProfile(histName, fCent, q2VecV0A[1]);

    histName = TString::Format("%s/hAvgQ3XvsCentV0CBef_%d", groupName.Data(), fRunNumber);
    fHistManager.FillProfile(histName, fCent, q3VecV0C[0]);
    histName = TString::Format("%s/hAvgQ3YvsCentV0CBef_%d", groupName.Data(), fRunNumber);
    fHistManager.FillProfile(histName, fCent, q3VecV0C[1]);
    histName = TString::Format("%s/hAvgQ3XvsCentV0ABef_%d", groupName.Data(), fRunNumber);
    fHistManager.FillProfile(histName, fCent, q3VecV0A[0]);
    histName = TString::Format("%s/hAvgQ3YvsCentV0ABef_%d", groupName.Data(), fRunNumber);
    fHistManager.FillProfile(histName, fCent, q3VecV0A[1]);

  }


  Int_t icentbin = 0;
  Double_t avgqx=0,avgqy=0; 

  // std::cout<<"before V0C q2x:q2y = "<< q2VecV0C[0] <<" : "<< q2VecV0C[1]<<std::endl;
  // std::cout<<"before V0A q2x:q2y = "<< q2VecV0A[0] <<" : "<< q2VecV0A[1]<<std::endl;
  
  if(fHCorrQ2xV0C && fHCorrQ2yV0C){
    icentbin = fHCorrQ2xV0C->FindBin(fCent);
    avgqx = fHCorrQ2xV0C->GetBinContent(icentbin);
    avgqy = fHCorrQ2yV0C->GetBinContent(icentbin);
    q2VecV0C[0] -= avgqx;
    q2VecV0C[1] -= avgqy;
    
    // std::cout << "before q2x:aveQ2x:q2y:aveQ2y = " << q2VecV0C[0] <<" : "<< avgqx\
    <<" : "<< q2VecV0C[1] <<" : "<< avgqy << std::endl;
  }
  if(fHCorrQ2xV0A && fHCorrQ2yV0A){
    icentbin = fHCorrQ2xV0A->FindBin(fCent);
    avgqx = fHCorrQ2xV0A->GetBinContent(icentbin);
    avgqy = fHCorrQ2yV0A->GetBinContent(icentbin);
    q2VecV0A[0] -= avgqx;
    q2VecV0A[1] -= avgqy;
  }

  if(fHCorrQ3xV0C && fHCorrQ3yV0C){
    icentbin = fHCorrQ3xV0C->FindBin(fCent);
    avgqx = fHCorrQ3xV0C->GetBinContent(icentbin);
    avgqy = fHCorrQ3yV0C->GetBinContent(icentbin);
    q3VecV0C[0] -= avgqx;
    q3VecV0C[1] -= avgqy;      

  }
  if(fHCorrQ3xV0A && fHCorrQ3yV0A){
    icentbin = fHCorrQ3xV0A->FindBin(fCent);
    avgqx = fHCorrQ3xV0A->GetBinContent(icentbin);
    avgqy = fHCorrQ3yV0A->GetBinContent(icentbin);
    q3VecV0A[0] -= avgqx;
    q3VecV0A[1] -= avgqy;           

  }
  
  if(fReCentCalibQA){
    TString histName;
    TString groupName;
    groupName="ReCentCalib";

    histName = TString::Format("%s/hAvgQ2XvsCentV0CAft_%d", groupName.Data(), fRunNumber);
    fHistManager.FillProfile(histName, fCent, q2VecV0C[0]);
    histName = TString::Format("%s/hAvgQ2YvsCentV0CAft_%d", groupName.Data(), fRunNumber);
    fHistManager.FillProfile(histName, fCent, q2VecV0C[1]);
    histName = TString::Format("%s/hAvgQ2XvsCentV0AAft_%d", groupName.Data(), fRunNumber);
    fHistManager.FillProfile(histName, fCent, q2VecV0A[0]);
    histName = TString::Format("%s/hAvgQ2YvsCentV0AAft_%d", groupName.Data(), fRunNumber);
    fHistManager.FillProfile(histName, fCent, q2VecV0A[1]);

    histName = TString::Format("%s/hAvgQ3XvsCentV0CAft_%d", groupName.Data(), fRunNumber);
    fHistManager.FillProfile(histName, fCent, q3VecV0C[0]);
    histName = TString::Format("%s/hAvgQ3YvsCentV0CAft_%d", groupName.Data(), fRunNumber);
    fHistManager.FillProfile(histName, fCent, q3VecV0C[1]);
    histName = TString::Format("%s/hAvgQ3XvsCentV0AAft_%d", groupName.Data(), fRunNumber);
    fHistManager.FillProfile(histName, fCent, q3VecV0A[0]);
    histName = TString::Format("%s/hAvgQ3YvsCentV0AAft_%d", groupName.Data(), fRunNumber);
    fHistManager.FillProfile(histName, fCent, q3VecV0A[1]);
  }

  return kTRUE;
}

//_____________________________________________________________________________
TH1F* AliAnalysisTaskEPCalibForJet::GetResoFromOutputFile(detectorType det, Int_t h, TArrayD* cen)
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
Double_t AliAnalysisTaskEPCalibForJet::CalculateEventPlaneChi(Double_t res)
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

void AliAnalysisTaskEPCalibForJet::SetupPileUpRemovalFunctions(){
  
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

//_____________________________________________________________________________
void AliAnalysisTaskEPCalibForJet::BkgFitEvaluation()
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
  // std::cout << "CDF = " << ChiSquarePDF(NDF, ChiSqr) << std::endl;

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



Bool_t AliAnalysisTaskEPCalibForJet::CheckEventIsPileUp2018(){
  
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


/**
 * This function is called once at the end of the analysis.
 */
void AliAnalysisTaskEPCalibForJet::Terminate(Option_t *) 
{
}


/**
 * This function adds the task to the analysis manager. Often, this function is called
 * by an AddTask C macro. However, by compiling the code, it ensures that we do not
 * have to deal with difficulties caused by CINT.
 */
AliAnalysisTaskEPCalibForJet * AliAnalysisTaskEPCalibForJet::AddTaskEPCalibForJet(
  const char *ntracks, const char *nclusters, const char* ncells, const char *suffix)
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEPCalibForJet", "No analysis manager to connect to.");
    return 0;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler)
  {
    ::Error("AddTaskEPCalibForJet", "This task requires an input event handler");
    return 0;
  }
  
  enum EDataType_t {kUnknown, kESD, kAOD};


  EDataType_t dataType = kAOD;
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------
  TString trackName(ntracks);
  if (trackName == "usedefault") trackName = "tracks";

  TString name("AliAnalysisTaskEPCalibForJet");
  if (strcmp(suffix,"") != 0) {
    name += "_";
    name += suffix;
  }

  AliAnalysisTaskEPCalibForJet* rawJetTask = new AliAnalysisTaskEPCalibForJet(name);
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










