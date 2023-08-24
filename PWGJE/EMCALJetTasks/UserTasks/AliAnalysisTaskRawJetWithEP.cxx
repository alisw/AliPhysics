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
 * 
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include <TRandom3.h>
#include <TMath.h>
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

#include "AliAODHandler.h"
#include "AliAnalysisManager.h"
#include "AliMultSelection.h"
#include "AliAODVertex.h"
#include "AliDataFile.h"
#include "AliLog.h"
#include "TGrid.h"


#include <AliAnalysisTask.h>
#include <AliVEventHandler.h>
#include <THistManager.h>

#include <AliVCluster.h>
#include <AliVParticle.h>
#include <AliLog.h>

#include "AliTLorentzVector.h"
#include "AliEmcalJet.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"

#include "AliLocalRhoParameter.h"
#include "AliRhoParameter.h"

#include "AliAnalysisTaskRawJetWithEP.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskRawJetWithEP);
/// \endcond

/**
 * Default constructor. Needed by ROOT I/O
 */
AliAnalysisTaskRawJetWithEP::AliAnalysisTaskRawJetWithEP() : 
    AliAnalysisTaskEmcalJet(),
    fOutputList(nullptr),
    fHistNEvents(nullptr),
    fCalibType(AliAnalysisTaskRawJetWithEP::kQnCalib),
    fNormMethod(AliAnalysisTaskRawJetWithEP::kQoverM),
    fOADBFileName(""),
    fMesLev(0),
    fAOD(nullptr),
    fPrevEventRun(-1),
    fTriggerClass(""),
    fTriggerMask(AliVEvent::kAny),
    fRejectTPCPileup(false),
    fUseAliEventCuts(kTRUE),
    fUseManualEventCuts(kFALSE),
    fRunListFileName(""),
    fCalibRefFileName(""),
    fCalibRefObjList(nullptr),
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
    fOwnEventCut(kFALSE),
    fDoEP(kFALSE),
    fDoTrack(kFALSE),
    fDoBkg(kFALSE),
    fDoJet(kFALSE),
    fEPQA(kFALSE),
    fTrackQA(kFALSE),
    fBkgQA(kFALSE),
    fJetQA(kFALSE),
    fSepEP(kFALSE),
    fV0Combin(kFALSE),
    fV0KindForBkg(0),
    fQnVCalibType("kOrig"),
    fHistManager(),
    fHCorrV0ChWeghts(nullptr),
    fHCorrQ2xV0C(nullptr),
    fHCorrQ2yV0C(nullptr),
    fHCorrQ2xV0A(nullptr),
    fHCorrQ2yV0A(nullptr),
    fHCorrQ3xV0C(nullptr),
    fHCorrQ3yV0C(nullptr),
    fHCorrQ3xV0A(nullptr),
    fHCorrQ3yV0A(nullptr),
    fQaEventNum(-1),
    fV2ResoV0(0.),
    fV3ResoV0(0.),
    fFitModulationType(kNoFit),
    fEventCuts(""),
    fEtaMinTPC(0.),
    fEtaMaxTPC(0.8),
    fPtMinTPC(0.2),
    fPtMaxTPC(5.),
    fUsedTrackPosIDs(),
    fUsedTrackNegIDs(),
    fFractionOfTracksForQnTPC(1.1),
    fV0(nullptr),
    fRun(0),
    fZvtx(-9999.),
    fCentrality(-9999.),
    fIsOADBFileOpen(false),
    fMultV0BefCorPfpx(0),
    fOADBzArray_contQx2am(0),
    fOADBzArray_contQy2am(0),
    fOADBzArray_contQx2as(0),
    fOADBzArray_contQy2as(0),
    fOADBzArray_contQx2cm(0),
    fOADBzArray_contQy2cm(0),
    fOADBzArray_contQx2cs(0),
    fOADBzArray_contQy2cs(0),
    fOADBzArray_contQx3am(0),
    fOADBzArray_contQy3am(0),
    fOADBzArray_contQx3as(0),
    fOADBzArray_contQy3as(0),
    fOADBzArray_contQx3cm(0),
    fOADBzArray_contQy3cm(0),
    fOADBzArray_contQx3cs(0),
    fOADBzArray_contQy3cs(0),
    fOADBcentArray_contTPCposEta(0),
    fOADBcentArray_contTPCnegEta(0),
    fCalibObjRun(-9999),
    fHistMultV0(nullptr),
    fV0CalibZvtxDiff(true),
    fEnablePhiDistrHistos(false)
{

    for(Int_t i(0); i < 8; i++) V0MultForAngle[i] = 0.;

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

    fUsedTrackPosIDs = TBits(1000);
    fUsedTrackNegIDs = TBits(1000);

    for(int iZvtx = 0; iZvtx < 14; iZvtx++) {
        fQx2mV0A[iZvtx] = nullptr;
        fQy2mV0A[iZvtx] = nullptr;
        fQx2sV0A[iZvtx] = nullptr;
        fQy2sV0A[iZvtx] = nullptr;
        fQx2mV0C[iZvtx] = nullptr;
        fQy2mV0C[iZvtx] = nullptr;
        fQx2sV0C[iZvtx] = nullptr;
        fQy2sV0C[iZvtx] = nullptr;

        fQx3mV0A[iZvtx] = nullptr;
        fQy3mV0A[iZvtx] = nullptr;
        fQx3sV0A[iZvtx] = nullptr;
        fQy3sV0A[iZvtx] = nullptr;
        fQx3mV0C[iZvtx] = nullptr;
        fQy3mV0C[iZvtx] = nullptr;
        fQx3sV0C[iZvtx] = nullptr;
        fQy3sV0C[iZvtx] = nullptr;  
    }

    for(int iCent = 0; iCent < 11; iCent++) {
        fWeightsTPCPosEta[iCent] = nullptr;
        fWeightsTPCNegEta[iCent] = nullptr;
    }

    fPhiVsCentrTPC[0]=nullptr;
    fPhiVsCentrTPC[1]=nullptr;

    fOADBzArray_contQx2am = new TObjArray();
    fOADBzArray_contQy2am = new TObjArray();
    fOADBzArray_contQx2as = new TObjArray();
    fOADBzArray_contQy2as = new TObjArray();
    fOADBzArray_contQx2cm = new TObjArray();
    fOADBzArray_contQy2cm = new TObjArray();
    fOADBzArray_contQx2cs = new TObjArray();
    fOADBzArray_contQy2cs = new TObjArray();
    
    fOADBzArray_contQx3am = new TObjArray();
    fOADBzArray_contQy3am = new TObjArray();
    fOADBzArray_contQx3as = new TObjArray();
    fOADBzArray_contQy3as = new TObjArray();
    fOADBzArray_contQx3cm = new TObjArray();
    fOADBzArray_contQy3cm = new TObjArray();
    fOADBzArray_contQx3cs = new TObjArray();
    fOADBzArray_contQy3cs = new TObjArray();

    fOADBcentArray_contTPCposEta = new TObjArray();
    fOADBcentArray_contTPCnegEta = new TObjArray();
}

/**
 * Standard constructor. Should be used by the user.
 *
 * @param[in] name Name of the task
 */
AliAnalysisTaskRawJetWithEP::AliAnalysisTaskRawJetWithEP(const char *name) : 
    AliAnalysisTaskEmcalJet(name, kTRUE),
    fOutputList(nullptr),
    fHistNEvents(nullptr),
    fCalibType(AliAnalysisTaskRawJetWithEP::kQnCalib),
    fNormMethod(AliAnalysisTaskRawJetWithEP::kQoverSqrtM),
    fOADBFileName(""),
    fMesLev(0),
    fAOD(nullptr),
    fPrevEventRun(-1),
    fTriggerClass(""),
    fTriggerMask(AliVEvent::kAny),
    fUseAliEventCuts(kTRUE),
    fUseManualEventCuts(kFALSE),
    fRunListFileName(""),
    fCalibRefFileName(""),
    fCalibRefObjList(nullptr),
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
    fOwnEventCut(kFALSE),
    fDoEP(kFALSE),
    fDoTrack(kFALSE),
    fDoBkg(kFALSE),
    fDoJet(kFALSE),
    fEPQA(kFALSE),
    fTrackQA(kFALSE),
    fBkgQA(kFALSE),
    fSepEP(kFALSE),
    fV0Combin(kFALSE),
    fV0KindForBkg(0),
    fQnVCalibType("kOrig"),
    fHistManager(name),
    fHCorrV0ChWeghts(nullptr),
    fHCorrQ2xV0C(nullptr),
    fHCorrQ2yV0C(nullptr),
    fHCorrQ2xV0A(nullptr),
    fHCorrQ2yV0A(nullptr),
    fHCorrQ3xV0C(nullptr),
    fHCorrQ3yV0C(nullptr),
    fHCorrQ3xV0A(nullptr),
    fHCorrQ3yV0A(nullptr),
    fQaEventNum(-1),
    fV2ResoV0(0.),
    fV3ResoV0(0.),
    fFitModulationType(kNoFit),
    fEventCuts(""),
    fEtaMinTPC(0.),
    fEtaMaxTPC(0.8),
    fPtMinTPC(0.2),
    fPtMaxTPC(5.),
    fUsedTrackPosIDs(),
    fUsedTrackNegIDs(),
    fFractionOfTracksForQnTPC(1.1),
    fV0(nullptr),
    fRun(0),
    fZvtx(-9999.),
    fCentrality(-9999.),
    fIsOADBFileOpen(false),
    fMultV0BefCorPfpx(0),
    fOADBzArray_contQx2am(0),
    fOADBzArray_contQy2am(0),
    fOADBzArray_contQx2as(0),
    fOADBzArray_contQy2as(0),
    fOADBzArray_contQx2cm(0),
    fOADBzArray_contQy2cm(0),
    fOADBzArray_contQx2cs(0),
    fOADBzArray_contQy2cs(0),
    fOADBzArray_contQx3am(0),
    fOADBzArray_contQy3am(0),
    fOADBzArray_contQx3as(0),
    fOADBzArray_contQy3as(0),
    fOADBzArray_contQx3cm(0),
    fOADBzArray_contQy3cm(0),
    fOADBzArray_contQx3cs(0),
    fOADBzArray_contQy3cs(0),
    fOADBcentArray_contTPCposEta(0),
    fOADBcentArray_contTPCnegEta(0),
    fCalibObjRun(-9999),
    fHistMultV0(nullptr),
    fV0CalibZvtxDiff(true),
    fEnablePhiDistrHistos(false)
{

    for(Int_t i(0); i < 8; i++) V0MultForAngle[i] = 0.;

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
        V0Mult2[i] = 0.;
        V0Mult3[i] = 0.;
        TpcMult2[i] = 0.;
        TpcMult3[i] = 0.;
        
        q2V0[i] = 0.;
        q3V0[i] = 0.;
        q2Tpc[i] = 0.;
        q3Tpc[i] = 0.;
    
        psi2V0[i] = 0.;
        psi3V0[i] = 0.;
        psi2Tpc[i] = 0.;
        psi3Tpc[i] = 0.;

        q2NormV0[i] = -9999.;
        q3NormV0[i] = -9999.;
        q2NormTpc[i] = -9999.;
        q3NormTpc[i] = -9999.;
    }

    if(fLocalRhoName=="") fLocalRhoName = Form("LocalRhoFrom_%s", GetName());
    SetMakeGeneralHistograms(kTRUE);

    fUsedTrackPosIDs = TBits(1000);
    fUsedTrackNegIDs = TBits(1000);

    for(int iZvtx = 0; iZvtx < 14; iZvtx++) {
        fQx2mV0A[iZvtx] = nullptr;
        fQy2mV0A[iZvtx] = nullptr;
        fQx2sV0A[iZvtx] = nullptr;
        fQy2sV0A[iZvtx] = nullptr;
        fQx2mV0C[iZvtx] = nullptr;
        fQy2mV0C[iZvtx] = nullptr;
        fQx2sV0C[iZvtx] = nullptr;
        fQy2sV0C[iZvtx] = nullptr; 

        fQx3mV0A[iZvtx] = nullptr;
        fQy3mV0A[iZvtx] = nullptr;
        fQx3sV0A[iZvtx] = nullptr;
        fQy3sV0A[iZvtx] = nullptr;
        fQx3mV0C[iZvtx] = nullptr;
        fQy3mV0C[iZvtx] = nullptr;
        fQx3sV0C[iZvtx] = nullptr;
        fQy3sV0C[iZvtx] = nullptr;  
    }

    for(int iCent = 0; iCent < 11; iCent++) {
        fWeightsTPCPosEta[iCent] = nullptr;
        fWeightsTPCNegEta[iCent] = nullptr;
    }

    fPhiVsCentrTPC[0]=nullptr;
    fPhiVsCentrTPC[1]=nullptr;

    fOADBzArray_contQx2am = new TObjArray();
    fOADBzArray_contQy2am = new TObjArray();
    fOADBzArray_contQx2as = new TObjArray();
    fOADBzArray_contQy2as = new TObjArray();
    fOADBzArray_contQx2cm = new TObjArray();
    fOADBzArray_contQy2cm = new TObjArray();
    fOADBzArray_contQx2cs = new TObjArray();
    fOADBzArray_contQy2cs = new TObjArray();

    fOADBzArray_contQx3am = new TObjArray();
    fOADBzArray_contQy3am = new TObjArray();
    fOADBzArray_contQx3as = new TObjArray();
    fOADBzArray_contQy3as = new TObjArray();
    fOADBzArray_contQx3cm = new TObjArray();
    fOADBzArray_contQy3cm = new TObjArray();
    fOADBzArray_contQx3cs = new TObjArray();
    fOADBzArray_contQy3cs = new TObjArray();

    fOADBcentArray_contTPCposEta = new TObjArray();
    fOADBcentArray_contTPCnegEta = new TObjArray();

    // Load OADB information from file during constructor
    // bool LoadedCalibrations = LoadOADBCalibrations();
    // if (!LoadedCalibrations) {
    //   AliError("Calibrations failed to load");
    // } else {
    //   AliInfo("Calibrations loaded correctly!\n");
    // }
}

/**
 * Destructor
 */
AliAnalysisTaskRawJetWithEP::~AliAnalysisTaskRawJetWithEP()
{
    // standard destructor
    if(fOutputList) {
		if(!fOutputList->IsOwner()) {
			delete fHistNEvents;
		}
		delete fOutputList;
    }

    if(fCalibRefObjList) {delete fCalibRefObjList; fCalibRefObjList = 0x0;}
    if(fFitModulation) {delete fFitModulation; fFitModulation = 0x0;}
}


//________________________________________________________________________
void AliAnalysisTaskRawJetWithEP::UserCreateOutputObjects()
{
    AliAnalysisTaskEmcalJet::UserCreateOutputObjects();
    
    // fOutputList = new TList();
    // fOutputList->SetOwner(true);
    
    /*
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
        fOutputList->Add(fEventCutList);
    } 
    */
    
    if(fOwnEventCut){
        fHistNEvents = new TH1F("fHistNEvents","Number of processed events;;Number of events",5,0.5,5.5);
        fHistNEvents->Sumw2();
        fHistNEvents->SetMinimum(0);
        fHistNEvents->GetXaxis()->SetBinLabel(1,"Read from AOD");
        fHistNEvents->GetXaxis()->SetBinLabel(2,"Pass Phys. Sel. + Trig");
        fHistNEvents->GetXaxis()->SetBinLabel(3,"Centrality > 90%");
        fHistNEvents->GetXaxis()->SetBinLabel(4,"No vertex");
        fHistNEvents->GetXaxis()->SetBinLabel(5,"Fail Pileup");
        // fOutputList->Add(fHistNEvents);
        fOutput->Add(fHistNEvents);
    }

    // == s == Set Out put Hist grams  ###########################################
    if(fEPQA)   AllocateEventPlaneHistograms();
    if(fBkgQA)  AllocateBkgHistograms();
    if(fJetQA)  AllocateJetHistograms();
    // == e == Set Out put Hist grams  ###########################################
    
    // == s == Add Objects into output file  ##################################### // previously make error for Run()
    TIter next(fHistManager.GetListOfHistograms());
    TObject* obj = 0;
    while ((obj = next())) {
        // fOutputList->Add(obj);
        fOutput->Add(obj);
    }
    // == e == Add Objects into output file  #####################################
    
    // post data
    if(fMesLev==3)std::cout << "UserCreateOutputObjects: post data ========="<< std::endl;
    
    // PostData(1, fOutputList);
    PostData(1, fOutput);
    
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

    for (Int_t cent = 0; cent < 11; cent++) {
        // == s == Event plane angle histograms Setting
        histName = TString::Format("%s/hPsi2V0AVsV0C_%d", groupName.Data(), cent);
        histtitle = "Psi2 from V0A vs V0C";
        fHistManager.CreateTH2(histName, histtitle, 50, 0.0, TMath::Pi(), 50, 0.0,TMath::Pi());
        histName = TString::Format("%s/hPsi2V0MVsV0C_%d", groupName.Data(), cent);
        histtitle = "Psi2 from V0M vs V0C";
        fHistManager.CreateTH2(histName, histtitle, 50, 0.0, TMath::Pi(), 50, 0.0,TMath::Pi());
        histName = TString::Format("%s/hPsi2V0MVsV0A_%d", groupName.Data(), cent);
        histtitle = "Psi2 from V0M vs V0A";
        fHistManager.CreateTH2(histName, histtitle, 50, 0.0,TMath::Pi(), 50, 0.0,TMath::Pi());

        histName = TString::Format("%s/hPsi2V0MVsTPCP_%d", groupName.Data(), cent);
        histtitle = "Psi2 from V0M vs TPC postitive";
        fHistManager.CreateTH2(histName, histtitle, 50, 0.0,TMath::Pi(), 50, 0.0,TMath::Pi());
        histName = TString::Format("%s/hPsi2V0MVsTPCN_%d", groupName.Data(), cent);
        histtitle = "Psi2 from V0M vs TPC negative";
        fHistManager.CreateTH2(histName, histtitle, 50, 0.0,TMath::Pi(), 50, 0.0,TMath::Pi());
        histName = TString::Format("%s/hPsi2TPCPVsTPCN_%d", groupName.Data(), cent);
        histtitle = "Psi2 from TPC posi vs nega";
        fHistManager.CreateTH2(histName, histtitle, 50, 0.0,TMath::Pi(), 50, 0.0,TMath::Pi());


        // profiles for all correlator permutations which are necessary to calculate each second and third order event plane resolution
        TProfile *tempHist;
        THashList *parent;
        histName = TString::Format("hProfV2Resolution_%d", cent);
        histtitle = histName;
        parent = static_cast<THashList *>(fHistManager.FindObject(groupName.Data()));
        tempHist = new TProfile(histName, histtitle, 12, -0.5, 11.5);
        tempHist->GetXaxis()->SetBinLabel(2,  "<cos(2(#Psi_{V0C} - #Psi_{V0A}))>");
        tempHist->GetXaxis()->SetBinLabel(3,  "<cos(2(#Psi_{V0C} - #Psi_{TPC}))>");
        tempHist->GetXaxis()->SetBinLabel(4,  "<cos(2(#Psi_{V0C} - #Psi_{TPC_N}))>");
        tempHist->GetXaxis()->SetBinLabel(5,  "<cos(2(#Psi_{V0C} - #Psi_{TPC_P}))>");
        tempHist->GetXaxis()->SetBinLabel(6,  "<cos(2(#Psi_{V0A} - #Psi_{TPC}))>");
        tempHist->GetXaxis()->SetBinLabel(7,  "<cos(2(#Psi_{V0A} - #Psi_{TPC_N}))>");
        tempHist->GetXaxis()->SetBinLabel(8,  "<cos(2(#Psi_{V0A} - #Psi_{TPC_P}))>");
        tempHist->GetXaxis()->SetBinLabel(9,  "<cos(2(#Psi_{V0} - #Psi_{TPC_N}))>");
        tempHist->GetXaxis()->SetBinLabel(10, "<cos(2(#Psi_{V0} - #Psi_{TPC_P}))>");
        tempHist->GetXaxis()->SetBinLabel(11, "<cos(2(#Psi_{TPC_N} - #Psi_{TPC_P}))>");
        parent->Add(tempHist);

        histName = TString::Format("hProfV3Resolution_%d", cent);
        histtitle = histName;
        parent = static_cast<THashList *>(fHistManager.FindObject(groupName.Data()));
        tempHist = new TProfile(histName, histtitle, 12, -0.5, 11.5);
        tempHist->GetXaxis()->SetBinLabel(2,  "<cos(3(#Psi_{V0C} - #Psi_{V0A}))>");
        tempHist->GetXaxis()->SetBinLabel(3,  "<cos(3(#Psi_{V0C} - #Psi_{TPC}))>");
        tempHist->GetXaxis()->SetBinLabel(4,  "<cos(3(#Psi_{V0C} - #Psi_{TPC_N}))>");
        tempHist->GetXaxis()->SetBinLabel(5,  "<cos(3(#Psi_{V0C} - #Psi_{TPC_P}))>");
        tempHist->GetXaxis()->SetBinLabel(6,  "<cos(3(#Psi_{V0A} - #Psi_{TPC}))>");
        tempHist->GetXaxis()->SetBinLabel(7,  "<cos(3(#Psi_{V0A} - #Psi_{TPC_N}))>");
        tempHist->GetXaxis()->SetBinLabel(8,  "<cos(3(#Psi_{V0A} - #Psi_{TPC_P}))>");
        tempHist->GetXaxis()->SetBinLabel(9,  "<cos(3(#Psi_{V0} - #Psi_{TPC_N}))>");
        tempHist->GetXaxis()->SetBinLabel(10, "<cos(3(#Psi_{V0} - #Psi_{TPC_P}))>");
        tempHist->GetXaxis()->SetBinLabel(11, "<cos(3(#Psi_{TPC_N} - #Psi_{TPC_P}))>");
        parent->Add(tempHist);
    }

    histName  = TString::Format("%s/v2", groupName.Data());
    histtitle = TString::Format("%s;centrality;v2", histName.Data());
    fHistManager.CreateTProfile(histName, histtitle, 10, 0, 100);
    histName  = TString::Format("%s/v3", groupName.Data());
    histtitle = TString::Format("%s;centrality;v3", histName.Data());
    fHistManager.CreateTProfile(histName, histtitle, 10, 0, 100);

    histName  = TString::Format("%s/CentQ2x_V0M", groupName.Data());
    histtitle = TString::Format("%s;centrality;Q2x_V0M", histName.Data());
    fHistManager.CreateTH2(histName, histtitle, 100, 0, 100, 100, 0, 1000);
    histName  = TString::Format("%s/CentQ2y_V0M", groupName.Data());
    histtitle = TString::Format("%s;centrality;Q2y_V0M", histName.Data());
    fHistManager.CreateTH2(histName, histtitle, 100, 0, 100, 100, 0, 1000);
    histName  = TString::Format("%s/CentQ2x_V0C", groupName.Data());
    histtitle = TString::Format("%s;centrality;Q2x_V0C", histName.Data());
    fHistManager.CreateTH2(histName, histtitle, 100, 0, 100, 100, 0, 1000);
    histName  = TString::Format("%s/CentQ2y_V0C", groupName.Data());
    histtitle = TString::Format("%s;centrality;Q2y_V0C", histName.Data());
    fHistManager.CreateTH2(histName, histtitle, 100, 0, 100, 100, 0, 1000);
    histName  = TString::Format("%s/CentQ2x_V0A", groupName.Data());
    histtitle = TString::Format("%s;centrality;Q2x_V0A", histName.Data());
    fHistManager.CreateTH2(histName, histtitle, 100, 0, 100, 100, 0, 1000);
    histName  = TString::Format("%s/CentQ2y_V0A", groupName.Data());
    histtitle = TString::Format("%s;centrality;Q2y_V0A", histName.Data());
    fHistManager.CreateTH2(histName, histtitle, 100, 0, 100, 100, 0, 1000);
    
    histName  = TString::Format("%s/CentQ3x_V0M", groupName.Data());
    histtitle = TString::Format("%s;centrality;Q3x_V0M", histName.Data());
    fHistManager.CreateTH2(histName, histtitle, 100, 0, 100, 100, 0, 1000);
    histName  = TString::Format("%s/CentQ3y_V0M", groupName.Data());
    histtitle = TString::Format("%s;centrality;Q3y_V0M", histName.Data());
    fHistManager.CreateTH2(histName, histtitle, 100, 0, 100, 100, 0, 1000);
    histName  = TString::Format("%s/CentQ3x_V0C", groupName.Data());
    histtitle = TString::Format("%s;centrality;Q3x_V0C", histName.Data());
    fHistManager.CreateTH2(histName, histtitle, 100, 0, 100, 100, 0, 1000);
    histName  = TString::Format("%s/CentQ3y_V0C", groupName.Data());
    histtitle = TString::Format("%s;centrality;Q3y_V0C", histName.Data());
    fHistManager.CreateTH2(histName, histtitle, 100, 0, 100, 100, 0, 1000);
    histName  = TString::Format("%s/CentQ3x_V0A", groupName.Data());
    histtitle = TString::Format("%s;centrality;Q3x_V0A", histName.Data());
    fHistManager.CreateTH2(histName, histtitle, 100, 0, 100, 100, 0, 1000);
    histName  = TString::Format("%s/CentQ3y_V0A", groupName.Data());
    histtitle = TString::Format("%s;centrality;Q3y_V0A", histName.Data());
    fHistManager.CreateTH2(histName, histtitle, 100, 0, 100, 100, 0, 1000);

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
    fHistManager.CreateTH2(histName, histtitle, 100, 0, 100, 40, 0, 1);
    
    histName = TString::Format("%s/hChi2Cent_lRhoCombinFit", groupName.Data());
    histtitle = TString::Format("%s; centrality; #tilde{#chi^{2}}", histName.Data());
    fHistManager.CreateTH2(histName, histtitle, 100, 0, 100, 100, 0, 5);

    histName = TString::Format("%s/hPChi2_lRhoCombinFit", groupName.Data());
    histtitle = TString::Format("%s; p-value; #tilde{#chi^{2}}", histName.Data());
    fHistManager.CreateTH2(histName, histtitle, 100, 0, 1, 100, 0, 5);

    histName = TString::Format("%s/hPvalueCDFROOT_lRhoCombinFit", groupName.Data());
    histtitle = TString::Format("%s; CDF #chi^{2} ROOT", histName.Data());
    fHistManager.CreateTH1(histName, histtitle, 50, 0, 1);

    histName = TString::Format("%s/hPvalueCDFROOTCent_lRhoCombinFit", groupName.Data());
    histtitle = TString::Format("%s; centrality; p-value ROOT", histName.Data());
    fHistManager.CreateTH2(histName, histtitle, 50, 0, 100, 45, 0, 1);

    histName = TString::Format("%s/hChi2ROOTCent_lRhoCombinFit", groupName.Data());
    histtitle = TString::Format("%s; p-value; #tilde{#chi^{2}}", histName.Data());
    fHistManager.CreateTH2(histName, histtitle, 100, 0, 1, 100, 0, 5);

    histName = TString::Format("%s/hPChi2ROOT_lRhoCombinFit", groupName.Data());
    histtitle = TString::Format("%s;CDF #chi^{2}; #tilde{#chi^{2}} ROOT", histName.Data());
    fHistManager.CreateTH2(histName, histtitle, 100, 0, 1, 100, 0, 5);

    //  = v2  fit local rho =
    histName = TString::Format("%s/hPvalueCDF_lRhoV2Fit", groupName.Data());
    histtitle = TString::Format("%s; CDF #chi^{2}", histName.Data());
    fHistManager.CreateTH1(histName, histtitle, 50, 0, 1);
    
    histName = TString::Format("%s/hPvalueCDFCent_lRhoV2Fit", groupName.Data());
    histtitle = TString::Format("%s; centrality; p-value", histName.Data());
    fHistManager.CreateTH2(histName, histtitle, 50, 0, 100, 100, 0, 1);
    
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
    fHistManager.CreateTH2(histName, histtitle, 50, 0, 100, 100, 0, 1);

    histName = TString::Format("%s/hChi2ROOTCent_lRhoV2Fit", groupName.Data());
    histtitle = TString::Format("%s; p-value; #tilde{#chi^{2}}", histName.Data());
    fHistManager.CreateTH2(histName, histtitle, 100, 0, 1, 100, 0, 5);

    histName = TString::Format("%s/hPChi2ROOT_lRhoV2Fit", groupName.Data());
    histtitle = TString::Format("%s;CDF #chi^{2}; #tilde{#chi^{2}} ROOT", histName.Data());
    fHistManager.CreateTH2(histName, histtitle, 100, 0, 1, 100, 0, 5);

        //  = fit global rho (rho0) =
    histName = TString::Format("%s/hPvalueCDF_gRhoFit", groupName.Data());
    histtitle = TString::Format("%s; CDF #chi^{2}", histName.Data());
    fHistManager.CreateTH1(histName, histtitle, 50, 0, 1);
    
    histName = TString::Format("%s/hPvalueCDFCent_gRhoFit", groupName.Data());
    histtitle = TString::Format("%s; centrality; p-value", histName.Data());
    fHistManager.CreateTH2(histName, histtitle, 50, 0, 100, 100, 0, 1);
    
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
    fHistManager.CreateTH2(histName, histtitle, 50, 0, 100, 100, 0, 1);

    histName = TString::Format("%s/hChi2ROOTCent_gRhoFit", groupName.Data());
    histtitle = TString::Format("%s; p-value; #tilde{#chi^{2}}", histName.Data());
    fHistManager.CreateTH2(histName, histtitle, 100, 0, 1, 100, 0, 5);

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
        histName = TString::Format("nEvents");
        histtitle = TString::Format("Number of Events");
        fHistManager.CreateTH1(histName, histtitle, 1, 0.0, 1.0);

        for (Int_t cent = 0; cent < 11; cent++) {
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

        
        for (Int_t cent = 0; cent < 11; cent++) {
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
            fHistManager.CreateTH1(histName, histtitle, 20, 0, 20);

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
            histName = TString::Format("%s/hPhiVsDeltaPt_Global_%d",GenGroupName.Data(),cent);
            histtitle = "phi vs delta pt according EP distribution";
            fHistManager.CreateTH2(histName, histtitle, 50, 0, TMath::TwoPi(), 100, -50, 50);
            histName = TString::Format("%s/hPhiVsDeltaPt_Local_%d", GenGroupName.Data(), cent);
            histtitle = "phi vs delta pt according EP distribution";
            fHistManager.CreateTH2(histName, histtitle, 50, 0, TMath::TwoPi(), 100, -50, 50);


            if (!jetCont->GetRhoName().IsNull()) {
                    
                    // == s ==  Rho histograms   11111111111111111111111111111111111111111111111111111111111111
                    histtitle = "GlobalRho";
                    histName = TString::Format("%s/hJetRho_%d", RhoGroupName.Data(), cent);
                    fHistManager.CreateTH1(histName, histtitle, 500, 0.0, 250.0);
                    histtitle = "LocalRho";
                    histName = TString::Format("%s/hJetRhoLocal_%d", RhoGroupName.Data(), cent);
                    fHistManager.CreateTH1(histName, histtitle, 500, 0.0, 250.0);
                    
                    histName = TString::Format("%s/hJetLRhoVsAveRho_%d", RhoGroupName.Data(), cent);
                    histtitle = "Local rho versus average rho";
                    fHistManager.CreateTH2(histName, histtitle, 500, 0.0, 250.0, 500, 0.0, 250.0);

                    // histo local rho vs delta phi
                    histName = TString::Format("%s/hJetGRhoVsDeltaPhi_%d", RhoGroupName.Data(), cent);
                    histtitle = "Global rho versus angle relative to event plane";
                    fHistManager.CreateTH2(histName, histtitle, 60, 0.0, TMath::TwoPi(), 500, 0.0, 250.0);
                    histName = TString::Format("%s/hJetLRhoVsDeltaPhi_%d", RhoGroupName.Data(), cent);
                    histtitle = "Local rho versus angle relative to event plane";
                    fHistManager.CreateTH2(histName, histtitle, 60, 0.0, TMath::TwoPi(), 5000, 0.0, 250.0);
                    
                    if(cent==0){
                        histName = TString::Format("%s/hCentJetPtCorrAPerCone", RhoGroupName.Data());
                        histtitle = histName + ";Centrality (%);#it{p}_{T, ch jet}^{raw}; #it{A}_{jet}/(#pi*#it{R^{2}})";
                        fHistManager.CreateTH3(histName.Data(), histtitle.Data(), 100, 0, 100, 250, -50, 200, 100, 0., 1.5);

                        histName = TString::Format("%s/hCentJetPtRawGRhoA", RhoGroupName.Data());
                        histtitle = histName + ";Centrality (%);#it{p}_{T, ch jet}^{raw}; median #rho * A [GeV/#it{c}]";
                        fHistManager.CreateTH3(histName.Data(), histtitle.Data(), 100, 0, 100, 250, 0, 250, 100, 0., 50.);
                        histName = TString::Format("%s/hCentJetPtRawLRhoA", RhoGroupName.Data());
                        histtitle = histName + ";Centrality (%);#it{p}_{T, ch jet}^{raw}; local #rho * A [GeV/#it{c}]";
                        fHistManager.CreateTH3(histName.Data(), histtitle.Data(), 100, 0, 100, 250, 0, 250, 100, 0., 50.);
                        

                        histName = "hCentGRhoLRho";
                        histtitle = histName + ";Centrality (%);#rho [GeV/#it{c} / (m^{2})]; local #rho [GeV/#it{c} / (m^{2})]";
                        fHistManager.CreateTH3(histName.Data(), histtitle.Data(), 100, 0, 100, 250, 0, 250, 250, 0., 250.);
                    }
                    // == e ==  Rho histograms   11111111111111111111111111111111111111111111111111111111111111


                    // = s = Create histograms of Jet Yeild ====================================================
                    //inlucive
                    histName = TString::Format("%s/hJetPt_%d", InclusiveGroupName.Data(), cent);
                    histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c});counts", histName.Data());
                    fHistManager.CreateTH1(histName, histtitle, 60, -50, 250);
                    histName = TString::Format("%s/hJetCorrPt_%d", InclusiveGroupName.Data(), cent);
                    histtitle = TString::Format("%s;#it{p}_{T,jet}^{corr} (GeV/#it{c});counts", histName.Data());
                    fHistManager.CreateTH1(histName, histtitle, 60, -50, 250);

                    histName = TString::Format("%s/hAzAngleJetPt_%d", InclusiveGroupName.Data(), cent);
                    histtitle = TString::Format("%s;#phi_{jet} - #Psi_{2, EP} [rad];#it{p}_{T,jet}^{corr} [GeV/#it{c}]", histName.Data());
                    fHistManager.CreateTH2(histName, histtitle, 100, 0.0, TMath::TwoPi()/4, 60, -50, 250);
                    histName = TString::Format("%s/hAzAngleJetCorrPt_%d", InclusiveGroupName.Data(), cent);
                    histtitle = TString::Format("%s;#phi_{jet} - #Psi_{2, EP} [rad];#it{p}_{T,jet}^{corr} [GeV/#it{c}]", histName.Data());
                    fHistManager.CreateTH2(histName, histtitle, 100, 0.0, TMath::TwoPi()/4, 60, -50, 250);
                    
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
        histName = TString::Format("%s/hMultJetPt", InclusiveGroupName.Data());
        histtitle = TString::Format("%s;Multiplicity;#it{p}_{T,jet}^{corr} [GeV/#it{c}]", histName.Data());
        fHistManager.CreateTH2(histName, histtitle, 2500, 0., 5000., 60, -50, 250);
        histName = TString::Format("%s/hMultJetCorrPt", InclusiveGroupName.Data());
        histtitle = TString::Format("%s;Multiplicity;#it{p}_{T,jet}^{corr} [GeV/#it{c}]", histName.Data());
        fHistManager.CreateTH2(histName, histtitle, 2500, 0., 5000., 60, -50, 250);

        histName = TString::Format("%s/hJetPtVsEP2AngleVsCent", InclusiveGroupName.Data());
        histtitle = TString::Format("%s;#it{p}_{T,jet}^{raw} [GeV/#it{c}]; #phi^{jet}-#Psi_{EP, 2} [rad];Centrality [%]", histName.Data());
        fHistManager.CreateTH3(histName, histtitle, 60, -50, 250, 100, 0., TMath::TwoPi(), 20, 0., 100.);
        histName = TString::Format("%s/hJetCorrPtVsEP2AngleVsCent", InclusiveGroupName.Data());
        histtitle = TString::Format("%s;#it{p}_{T,jet}^{corr} [GeV/#it{c}];#phi^{jet}-#Psi_{EP, 2} [rad];Centrality [%]", histName.Data());
        fHistManager.CreateTH3(histName, histtitle, 60, -50, 250, 100, 0., TMath::TwoPi(), 20, 0., 100.);
        // fHistManager.CreateTH3(histName, histtitle, 60, -50, 250, 100, 0., TMath::Pi(), 20, 0., 100.);
        histName = TString::Format("%s/hLocalRhoAVsEP2AngleVsCent", InclusiveGroupName.Data());
        histtitle = TString::Format("%s;#it{#rho}_{ch} (#it{#phi}) [GeV/#it{c}/ (m^{2})];#phi^{jet}-#Psi_{EP, 2} [rad];Centrality [%]", histName.Data());
        fHistManager.CreateTH3(histName, histtitle, 125, -20, 230, 100, 0., TMath::TwoPi(), 20, 0., 100.);

        histName = TString::Format("%s/hJetPtVsEP3AngleVsCent", InclusiveGroupName.Data());
        histtitle = TString::Format("%s;#it{p}_{T,jet}^{raw} [GeV/#it{c}];#phi^{jet}-#Psi_{EP, 3} [rad];Centrality [%]", histName.Data());
        fHistManager.CreateTH3(histName, histtitle, 60, -50, 250, 100, 0., TMath::TwoPi(), 20, 0., 100.);
        histName = TString::Format("%s/hJetCorrPtVsEP3AngleVsCent", InclusiveGroupName.Data());
        histtitle = TString::Format("%s;#it{p}_{T,jet}^{corr} [GeV/#it{c}];#phi^{jet}-#Psi_{EP, 3} [rad];Centrality [%]", histName.Data());
        fHistManager.CreateTH3(histName, histtitle, 60, -50, 250, 100, 0., TMath::TwoPi(), 20, 0., 100.);
        // fHistManager.CreateTH3(histName, histtitle, 60, -50, 250, 100, 0., TMath::TwoPi()*2/3, 20, 0., 100.);
        histName = TString::Format("%s/hLocalRhoAVsEP3AngleVsCent", InclusiveGroupName.Data());
        histtitle = TString::Format("%s;#it{#rho}_{ch} (#it{#phi}) [GeV/#it{c}/ (m^{2})];#phi^{jet}-#Psi_{EP, 3} [rad];Centrality [%]", histName.Data());
        fHistManager.CreateTH3(histName, histtitle, 125, -20, 230, 100, 0., TMath::TwoPi(), 20, 0., 100.);
    }
}


/**
 * This function is executed automatically for the first event.
 * Some extra initialization can be performed here.
 */
void AliAnalysisTaskRawJetWithEP::ExecOnce()
{
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    
    if(fMesLev==3)std::cout << "Start setup reference calibration values!!" << std::endl;
    TString EPCalibRefFileName = "";
    if(fQnVCalibType == "kOrig"){
        TList *lCalibRefHists = nullptr;
        
        if (fAOD->GetRunNumber() < 295584) EPCalibRefFileName = "alien:///alice/cern.ch/user/t/tkumaoka/CalibV0GainCorrectionLHC15o_Oct2021.root"; //LHC15o
        else if (fAOD->GetRunNumber() <= 296623) EPCalibRefFileName = "alien:///alice/cern.ch/user/t/tkumaoka/CalibV0GainCorrectionLHC18q_Oct2021.root";
        else if (fAOD->GetRunNumber() >  296623) EPCalibRefFileName = "alien:///alice/cern.ch/user/t/tkumaoka/CalibV0GainCorrectionLHC18r_Oct2021.root";
        TString tempCalibFileName = AliDataFile::GetFileName(EPCalibRefFileName.Data());
        TString tempCalibLocalFileName = EPCalibRefFileName;
        
        // Check access to CVMFS (will only be displayed locally)
        if(EPCalibRefFileName.BeginsWith("alien://") && !gGrid){
            TGrid::Connect("alien://");
        }
        TFile* EPCalibRefFile = nullptr;
        if(!tempCalibFileName.IsNull()) EPCalibRefFile = TFile::Open(tempCalibFileName.Data());
        if(tempCalibFileName.IsNull())  EPCalibRefFile = TFile::Open(tempCalibLocalFileName.Data());
        lCalibRefHists = (TList *)EPCalibRefFile->Get("fWgtsV0ZDC");
        
        SetCalibOrigRefObjList(lCalibRefHists);
    }
    else if(fQnVCalibType == "kJeHand"){
        if (fAOD->GetRunNumber() < 295584) EPCalibRefFileName = "alien:///alice/cern.ch/user/t/tkumaoka/calibV0TPCRun2Vtx10P115oPass2.root"; //LHC15o
        else if (fAOD->GetRunNumber() <= 296623) EPCalibRefFileName = "alien:///alice/cern.ch/user/t/tkumaoka/calibV0TPCRun2Vtx10P118qPass3.root"; //LHC18q
        else if (fAOD->GetRunNumber() > 296623) EPCalibRefFileName = "alien:///alice/cern.ch/user/t/tkumaoka/calibV0TPCRun2Vtx10P118rPass3.root"; //LHC18r

        TString pathToFileLocal = EPCalibRefFileName;

        TString tempCalibFileName = AliDataFile::GetFileName(EPCalibRefFileName.Data());
        TString tempCalibLocalFileName = EPCalibRefFileName;
        
        // Check access to CVMFS (will only be displayed locally)
        if(EPCalibRefFileName.BeginsWith("alien://") && !gGrid){
            // AliInfo("Trying to connect to AliEn ...");
            TGrid::Connect("alien://");
        }
        TFile* EPCalibRefFile = nullptr;
        if(!tempCalibFileName.IsNull()) EPCalibRefFile = TFile::Open(tempCalibFileName.Data());
        if(tempCalibFileName.IsNull())  EPCalibRefFile = TFile::Open(tempCalibLocalFileName.Data());

        AliOADBContainer *lRefMultV0BefCorPfpx = new AliOADBContainer();
        TObjArray *lRefQx2am = new TObjArray();
        TObjArray *lRefQy2am = new TObjArray();
        TObjArray *lRefQx2as = new TObjArray();
        TObjArray *lRefQy2as = new TObjArray();
        TObjArray *lRefQx3am = new TObjArray();
        TObjArray *lRefQy3am = new TObjArray();
        TObjArray *lRefQx3as = new TObjArray();
        TObjArray *lRefQy3as = new TObjArray();
        TObjArray *lRefQx2cm = new TObjArray();
        TObjArray *lRefQy2cm = new TObjArray();
        TObjArray *lRefQx2cs = new TObjArray();
        TObjArray *lRefQy2cs = new TObjArray();
        TObjArray *lRefQx3cm = new TObjArray();
        TObjArray *lRefQy3cm = new TObjArray(); 
        TObjArray *lRefQx3cs = new TObjArray(); 
        TObjArray *lRefQy3cs = new TObjArray();
        TObjArray *lRefTPCposEta = new TObjArray();
        TObjArray *lRefTPCnegEta = new TObjArray();


        lRefMultV0BefCorPfpx = (AliOADBContainer *) EPCalibRefFile->Get("hMultV0BefCorPfpx");

        bool LoadedCaliRef = ExtractRecentPara(EPCalibRefFile, lRefQx2am, lRefQy2am, lRefQx2as, lRefQy2as, lRefQx3am, lRefQy3am, lRefQx3as, lRefQy3as, lRefQx2cm, lRefQy2cm, lRefQx2cs, lRefQy2cs, lRefQx3cm,lRefQy3cm, lRefQx3cs, lRefQy3cs, lRefTPCposEta, lRefTPCnegEta);
        if (!LoadedCaliRef) {
            std::cout << "Calibrations failed to load!\n" << std::endl;
        } else {
            std::cout << "Calibrations loaded correctly!\n" << std::endl;
        }
        
        SetLRefMultV0BefCorPfpx(lRefMultV0BefCorPfpx);
        SetLRefQx2am(lRefQx2am);
        SetLRefQy2am(lRefQy2am);
        SetLRefQx2as(lRefQx2as);
        SetLRefQy2as(lRefQy2as);
        SetLRefQx3am(lRefQx3am);
        SetLRefQy3am(lRefQy3am);
        SetLRefQx3as(lRefQx3as);
        SetLRefQy3as(lRefQy3as);
        SetLRefQx2cm(lRefQx2cm);
        SetLRefQy2cm(lRefQy2cm);
        SetLRefQx2cs(lRefQx2cs);
        SetLRefQy2cs(lRefQy2cs);
        SetLRefQx3cm(lRefQx3cm);
        SetLRefQy3cm(lRefQy3cm);
        SetLRefQx3cs(lRefQx3cs);
        SetLRefQy3cs(lRefQy3cs);
        SetLRefTPCposEta(lRefTPCposEta);
        SetLRefTPCnegEta(lRefTPCnegEta);
        
        if(lRefMultV0BefCorPfpx) delete lRefMultV0BefCorPfpx;
        if(lRefQx2am) delete lRefQx2am;
        if(lRefQy2am) delete lRefQy2am;
        if(lRefQx2as) delete lRefQx2as;
        if(lRefQy2as) delete lRefQy2as;
        if(lRefQx3am) delete lRefQx3am;
        if(lRefQy3am) delete lRefQy3am;
        if(lRefQx3as) delete lRefQx3as;
        if(lRefQy3as) delete lRefQy3as;
        if(lRefQx2cm) delete lRefQx2cm;
        if(lRefQy2cm) delete lRefQy2cm;
        if(lRefQx2cs) delete lRefQx2cs;
        if(lRefQy2cs) delete lRefQy2cs;
        if(lRefQx3cm) delete lRefQx3cm;
        if(lRefQy3cm) delete lRefQy3cm;
        if(lRefQx3cs) delete lRefQx3cs;
        if(lRefQy3cs) delete lRefQy3cs;
        if(lRefTPCposEta) delete lRefTPCposEta;
        if(lRefTPCnegEta) delete lRefTPCnegEta;
    }
    if(fMesLev==3)std::cout << "Finish setup reference calibration values!!" << std::endl;

    if(fMesLev==3)std::cout << "ExecOnce: Start load fLocalRho ================" << std::endl;
    if(!fLocalRho) {
        fLocalRho = new AliLocalRhoParameter(fLocalRhoName.Data(), 0); 
        if(!(InputEvent()->FindListObject(fLocalRho->GetName()))) {
            InputEvent()->AddObject(fLocalRho);
        } else {
            AliFatal(Form("%s: Container with name %s already present. Aborting", \
            GetName(), fLocalRho->GetName()));
        }
    }
    
    if(fMesLev==3)std::cout << "ExecOnce: Start AliAnalysisTaskEmcalJet::ExecOnce() ="<< std::endl;
    
    AliAnalysisTaskEmcalJet::ExecOnce();
    if(!GetJetContainer()) AliFatal(Form("%s: Couldn't find jet container. Aborting !", GetName()));
    
}

/**
 * This function (overloading the base class) uses AliEventCuts to perform event selection.
 */
/*
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
*/


//________________________________________________________________________
Bool_t AliAnalysisTaskRawJetWithEP::Run()
{
    // main event loop
    if(fMesLev==3)std::cout << "Start Run() function #########################" << std::endl;
    if(fOwnEventCut){
        // main event loop
        fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
        if (!fAOD && AODEvent() && IsStandardAOD()) {
            // In case there is an AOD handler writing a standard AOD, use the AOD
            // event in memory rather than the input (ESD) event.
            fAOD = dynamic_cast<AliAODEvent*> (AODEvent());
        }
        if (!fAOD) {
            AliWarning("AliAnalysisTaskRawJetWithEP::Exec(): bad AOD");
            return kFALSE;
        }

        fHistNEvents->Fill(1);

        if(TMath::Abs(fAOD->GetMagneticField())<0.001) return kFALSE;

        AliAODHandler* aodHandler = static_cast<AliAODHandler*>((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
        if(!aodHandler) {
            AliWarning("AliAnalysisTaskRawJetWithEP::Exec(): No AliInputEventHandler!");
            return kFALSE;
        }

        unsigned int maskPhysSel = aodHandler->IsEventSelected();
        TString firedTriggerClasses = fAOD->GetFiredTriggerClasses();
        if((fAOD->GetRunNumber()<136851 || fAOD->GetRunNumber()>139517)) {
            if(!(firedTriggerClasses.Contains(fTriggerClass.Data()))) return kFALSE;
        }
        if((maskPhysSel & fTriggerMask)==0.) {
            return kFALSE;
        }

        AliMultSelection *multSelection = dynamic_cast<AliMultSelection*>(fAOD->FindListObject("MultSelection"));
        if(!multSelection){
            AliWarning("AliMultSelection could not be found in the aod event list of objects");
            return kFALSE;
        }
        fHistNEvents->Fill(2);

        float cent = multSelection->GetMultiplicityPercentile("V0M");
        if(cent>90) {
            fHistNEvents->Fill(3);
            return kFALSE;
        }

        const AliAODVertex *vertex = fAOD->GetPrimaryVertex();
        if(!vertex || TMath::Abs(vertex->GetZ())>10. || vertex->GetNContributors()<=0) {
            fHistNEvents->Fill(4);
            return kFALSE;
        }
        if (fRejectTPCPileup)   {
            fEventCuts.SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE,1);
            Bool_t acceptEventCuts = fEventCuts.AcceptEvent(InputEvent());
            if(!acceptEventCuts)   { 
                fHistNEvents->Fill(5);
                return kFALSE;}
        }
        int runnumber = fAOD->GetRunNumber();
        fPrevEventRun = runnumber;
    }

    if(fDoEP) DoEventPlane();
    if(fDoTrack) MeasureTpcEPQA();
    if(fDoJet) DoJetLoop();

    // Post output data
    if(fMesLev==3)std::cout << "Run(): Post output data ================" << std::endl;
    
    // PostData(1, fOutputList);
    PostData(1, fOutput);

    return kTRUE;
}


Bool_t AliAnalysisTaskRawJetWithEP::DoEventPlane(){
    if (!fAOD && AODEvent() && IsStandardAOD()) {
        // In case there is an AOD handler writing a standard AOD, use the AOD
        // event in memory rather than the input (ESD) event.
        fAOD = dynamic_cast<AliAODEvent*> (AODEvent());
    }
    if (!fAOD) {
        AliWarning("AliAnalysisTaskRawJetWithEP::Exec(): bad AOD");
        return kFALSE;
    }
    AliAODHandler* aodHandler = static_cast<AliAODHandler*>((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
    if(!aodHandler) {
        AliWarning("AliAnalysisTaskJetQnVectors::Exec(): No AliInputEventHandler!");
        return kFALSE;
    }

    //== s == qn Calibration  111111111111111111111111111111111111111111111111111
    if(fQnVCalibType == "kOrig"){
        if(!fCalibRefObjList){
            AliWarning("AliAnalysisTaskRawJetWithEP::: No fCalibRefObjList!!!!");
            return kFALSE;
        }
        QnV0GainCalibration();
        QnRecenteringCalibration();

        //== s == combin V0C and V0A  ################################################
        Double_t q2ChiV0C = 0.;
        Double_t q2ChiV0A = 0.;
        Double_t q3ChiV0C = 0.;
        Double_t q3ChiV0A = 0.;
        if(fV0Combin){
            if(fCent < 15.){
            q2ChiV0C = 0.78059*TMath::Exp(-0.5*((fCent-21.4249)/19.1914)*((fCent-21.4249)/19.1914));
            }
            else if((fCent >= 15.)&&(fCent < 55.)){
                q2ChiV0C = 0.76684*TMath::Exp(-0.5*((fCent-24.8316)/35.6943)*((fCent-24.8316)/35.6943));
            }
            else if(fCent > 55.){
                q2ChiV0C = 0.60019*TMath::Exp(-0.5*((fCent-44.9328)/19.7624)*((fCent-44.9328)/19.7624));
            }

            if(fCent < 15.){
                q2ChiV0A = 0.6914*TMath::Exp(-0.5*((fCent-21.457)/17.647)*((fCent-21.457)/17.647));
            }
            else if((fCent >= 15.)&&(fCent < 55.)){
                q2ChiV0A = 0.6772*TMath::Exp(-0.5*((fCent-24.758)/32.234)*((fCent-24.758)/32.234));
            }
            else if(fCent > 55.){
                q2ChiV0A = 0.5153*TMath::Exp(-0.5*((fCent-43.046)/20.165)*((fCent-43.046)/20.165));
            }
            
            if(fCent < 45.){
                q3ChiV0C = 0.3233*TMath::Exp(-0.5*((fCent-9.8289)/33.896)*((fCent-9.8289)/33.896));
            }
            else if(fCent > 45.){
                q3ChiV0C = 0.21312*TMath::Exp(-0.5*((fCent-36.463)/17.963)*((fCent-36.463)/17.963));
            }

            if(fCent < 55.){
                q3ChiV0A = 0.2325*TMath::Exp(-0.5*((fCent-11.805)/30.671)*((fCent-11.805)/30.671));
            }
            else if(fCent > 55.){
                q3ChiV0A = 0.1942*TMath::Exp(-0.5*((fCent-24.688)/24.008)*((fCent-24.688)/24.008));
            }
        }else {
            q2ChiV0C = 1.0;
            q3ChiV0C = 1.0;
            q2ChiV0A = 1.0;
            q3ChiV0A = 1.0;
        }
        q2VecV0M[0] = q2ChiV0C*q2ChiV0C*q2VecV0C[0] + q2ChiV0A*q2ChiV0A*q2VecV0A[0];
        q2VecV0M[1] = q2ChiV0C*q2ChiV0C*q2VecV0C[1] + q2ChiV0A*q2ChiV0A*q2VecV0A[1];
        q3VecV0M[0] = q3ChiV0C*q2ChiV0C*q3VecV0C[0] + q3ChiV0A*q2ChiV0A*q3VecV0A[0];
        q3VecV0M[1] = q3ChiV0C*q2ChiV0C*q3VecV0C[1] + q3ChiV0A*q2ChiV0A*q3VecV0A[1];
        //== s == combin V0C and V0A  ################################################

        psi2V0[0] = CalcEPAngle(q2VecV0M[0], q2VecV0M[1]);
        psi2V0[1] = CalcEPAngle(q2VecV0C[0], q2VecV0C[1]);
        psi2V0[2] = CalcEPAngle(q2VecV0A[0], q2VecV0A[1]);

        psi3V0[0] = CalcEPAngle(q3VecV0M[0], q3VecV0M[1]);
        psi3V0[1] = CalcEPAngle(q3VecV0C[0], q3VecV0C[1]);
        psi3V0[2] = CalcEPAngle(q3VecV0A[0], q3VecV0A[1]);
    }
    else if(fQnVCalibType == "kJeHand"){
        if(!fMultV0BefCorPfpx){
            AliWarning("AliAnalysisTaskRawJetWithEP::: No fMultV0BefCorPfpx!!!!!!");
            return kFALSE;
        }
        QnJEHandlarEPGet();
    }
    if(fQnVCalibType == "kNon"){
        QnCalcWOCalib();

        q2VecV0M[0] = q2VecV0C[0] + q2VecV0A[0];
        q2VecV0M[1] = q2VecV0C[1] + q2VecV0A[1];
        q3VecV0M[0] = q3VecV0C[0] + q3VecV0A[0];
        q3VecV0M[1] = q3VecV0C[1] + q3VecV0A[1];

        psi2V0[0] = CalcEPAngle(q2VecV0M[0], q2VecV0M[1]);
        psi2V0[1] = CalcEPAngle(q2VecV0C[0], q2VecV0C[1]);
        psi2V0[2] = CalcEPAngle(q2VecV0A[0], q2VecV0A[1]);

        psi3V0[0] = CalcEPAngle(q3VecV0M[0], q3VecV0M[1]);
        psi3V0[1] = CalcEPAngle(q3VecV0C[0], q3VecV0C[1]);
        psi3V0[2] = CalcEPAngle(q3VecV0A[0], q3VecV0A[1]);
    }
    //== e == qn Calibration  111111111111111111111111111111111111111111111111111
    
    if(fEPQA){
        TString histName;
        TString groupName;
        groupName="EventPlane";
    
        histName = TString::Format("%s/hCentrality", groupName.Data());
        fHistManager.FillTH1(histName, fCent);

        Int_t iCentBin = GetCentBin(); // fCentBin

        histName = TString::Format("%s/hPsi2V0AVsV0C_%d", groupName.Data(), iCentBin);
        fHistManager.FillTH2(histName, psi2V0[2], psi2V0[1]);
        histName = TString::Format("%s/hPsi2V0MVsV0C_%d", groupName.Data(), iCentBin);
        fHistManager.FillTH2(histName, psi2V0[0], psi2V0[1]);
        histName = TString::Format("%s/hPsi2V0MVsV0A_%d", groupName.Data(), iCentBin);
        fHistManager.FillTH2(histName, psi2V0[0], psi2V0[2]);
        histName = TString::Format("%s/hPsi2V0MVsTPCP_%d", groupName.Data(),iCentBin);
        fHistManager.FillTH2(histName, psi2V0[0], psi2Tpc[1]);
        histName = TString::Format("%s/hPsi2V0MVsTPCN_%d", groupName.Data(), iCentBin);
        fHistManager.FillTH2(histName, psi2V0[0], psi2Tpc[2]);
        histName = TString::Format("%s/hPsi2TPCPVsTPCN_%d", groupName.Data(), iCentBin);
        fHistManager.FillTH2(histName, psi2Tpc[1], psi2Tpc[2]);


        histName = TString::Format("%s/hProfV2Resolution_%d", groupName.Data(), iCentBin);
        fHistManager.FillProfile(histName, 1., TMath::Cos(2.*(psi2V0[1] - psi2V0[2])));
        fHistManager.FillProfile(histName, 2., TMath::Cos(2.*(psi2V0[1] - psi2Tpc[0])));
        fHistManager.FillProfile(histName, 3., TMath::Cos(2.*(psi2V0[1] - psi2Tpc[1])));
        fHistManager.FillProfile(histName, 4., TMath::Cos(2.*(psi2V0[1] - psi2Tpc[2])));
        fHistManager.FillProfile(histName, 5., TMath::Cos(2.*(psi2V0[2] - psi2Tpc[0])));
        fHistManager.FillProfile(histName, 6., TMath::Cos(2.*(psi2V0[2] - psi2Tpc[1])));
        fHistManager.FillProfile(histName, 7., TMath::Cos(2.*(psi2V0[2] - psi2Tpc[2])));
        fHistManager.FillProfile(histName, 8., TMath::Cos(2.*(psi2V0[0] - psi2Tpc[1])));
        fHistManager.FillProfile(histName, 9., TMath::Cos(2.*(psi2V0[0] - psi2Tpc[2])));
        fHistManager.FillProfile(histName, 10., TMath::Cos(2.*(psi2Tpc[1] - psi2Tpc[2])));
            
        histName = TString::Format("%s/hProfV3Resolution_%d", groupName.Data(), iCentBin);
        fHistManager.FillProfile(histName, 1., TMath::Cos(2.*(psi3V0[1] - psi3V0[2])));
        fHistManager.FillProfile(histName, 2., TMath::Cos(2.*(psi3V0[1] - psi3Tpc[0])));
        fHistManager.FillProfile(histName, 3., TMath::Cos(2.*(psi3V0[1] - psi3Tpc[1])));
        fHistManager.FillProfile(histName, 4., TMath::Cos(2.*(psi3V0[1] - psi3Tpc[2])));
        fHistManager.FillProfile(histName, 5., TMath::Cos(2.*(psi3V0[2] - psi3Tpc[0])));
        fHistManager.FillProfile(histName, 6., TMath::Cos(2.*(psi3V0[2] - psi3Tpc[1])));
        fHistManager.FillProfile(histName, 7., TMath::Cos(2.*(psi3V0[2] - psi3Tpc[2])));
        fHistManager.FillProfile(histName, 8., TMath::Cos(2.*(psi3V0[0] - psi3Tpc[1])));
        fHistManager.FillProfile(histName, 9., TMath::Cos(2.*(psi3V0[0] - psi3Tpc[2])));
        fHistManager.FillProfile(histName, 10., TMath::Cos(2.*(psi3Tpc[1] - psi3Tpc[2])));
    }

    return kTRUE;
}


Bool_t AliAnalysisTaskRawJetWithEP::MeasureBkg(Double_t baseJetRho){
    if(fMesLev==3)std::cout << "=== Start MeasureBkg()  =====================" << std::endl;
    
    if(fLocalRho) {
        if(fMesLev==3)std::cout << "Run(): SetVal of fLocalRho" << std::endl;
        // if(!fRho)std::cout << "Don't Find fRho" << std::endl;
        if(fMesLev==3)std::cout << "Run(): fRho->GetVal() = " << baseJetRho<< std::endl;
        // fLocalRho->SetVal(fRho->GetVal());
        fLocalRho->SetVal(baseJetRho);
    }
    else std::cout << "Run(): Cannot find fLocalRho" << std::endl;

    TString groupName;
    TString histName;

    // AliAnalysisTaskJetV2 ==================================================================
    // Int_t iTracks(fTracks->GetEntriesFast()); //????
    Double_t excludeInEta = -999;
    Double_t excludeInPhi = -999;
    Double_t excludeInPt  = -999;
    if(fLocalRho->GetVal() <= 0 ) return kFALSE;   // no use fitting an empty event ...
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
    
    Int_t fNAcceptedTracks = 0;
    AliParticleContainer * partContForBKG = 0;
    AliTLorentzVector trackVec;
    AliVParticle* track;
    Double_t etaTPC = 0.9;
    Double_t trackEta = 0.;
    Double_t trackPhi = 0.;
    Double_t trackPt = 0.;
    Double_t deltaR = 0.;
    Double_t trackID = -1.;
    TIter nextPartCont(&fParticleCollArray);
    while ((partContForBKG = static_cast<AliParticleContainer*>(nextPartCont()))) {
        TString partContName = partContForBKG->GetName();
        
        for (auto trackIterator : partContForBKG->accepted_momentum() ) {
            trackVec.Clear();
            trackVec = trackIterator.first;
            trackEta = trackVec.Eta();
            trackPhi = trackVec.Phi_0_2pi();
            trackPt  = trackVec.Pt();
            //To get track ID and particle pointer
            track = trackIterator.second;
            TClonesArray* fTracksContArray = partContForBKG->GetArray();
            trackID = fTracksContArray->IndexOf(track);
            if (TMath::Abs(trackEta) < etaTPC) fNAcceptedTracks++;
        }
        
    }

    TH1F _tempSwap;     // on stack for quick access
    TH1F _tempSwapN;    // on stack for quick access, bookkeeping histogram

    if(fNAcceptedTracks < 49) fNAcceptedTracks = 49;       // avoid aliasing effects
    _tempSwap = TH1F("_tempSwap", "_tempSwap", TMath::CeilNint(TMath::Sqrt(fNAcceptedTracks)), lowBound, upBound);
    _tempSwapN = TH1F("_tempSwapN", "_tempSwapN", TMath::CeilNint(TMath::Sqrt(fNAcceptedTracks)), lowBound, upBound);
    _tempSwap.Sumw2();

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
            fitFunction = "[0]*(1.+2.*[1]*TMath::Cos(2.*(x-[2])))";
            fFitModulation = new TF1("fit_kV2", fitFunction, 0, TMath::TwoPi());
            fFitModulation->SetParameter(0, 0.);  // normalization
            fFitModulation->SetParameter(1, 0.2); // v2
            fFitModulation->FixParameter(2, psi2V0[fV0KindForBkg]);
        } break;
        case kCombined: {
            fitFunction = "[0]*(1.+2.*([1]*TMath::Cos(2.*(x-[2]))+[3]*TMath::Cos(3.*(x-[4]))))";
            fFitModulation = new TF1("fit_kCombined", fitFunction, 0, TMath::TwoPi());
            fFitModulation->SetParameter(0, 0.);       // normalization
            fFitModulation->SetParameter(1, 0.2);      // v2
            fFitModulation->SetParameter(3, 0.2);      // v3

            fFitModulation->SetParameter(0, baseJetRho);
            fFitModulation->FixParameter(2, psi2V0[fV0KindForBkg]);
            fFitModulation->FixParameter(4, psi3V0[fV0KindForBkg]);
        } break;
        default : { // for the combined fit, the 'direct fourier series' or the user supplied vn values we use v2 and v3
            fitFunction = "[0]*(1.+2.*([1]*TMath::Cos(2.*(x-[2]))+[3]*TMath::Cos(3.*(x-[4]))))";
            fFitModulation = new TF1("fit_kCombined", fitFunction, 0, TMath::TwoPi());
            fFitModulation->SetParameter(0, 0.);       // normalization
            fFitModulation->SetParameter(1, 0.2);      // v2
            fFitModulation->SetParameter(3, 0.2);      // v3

            fFitModulation->SetParameter(0, baseJetRho);
            fFitModulation->FixParameter(2, psi2V0[fV0KindForBkg]);
            fFitModulation->FixParameter(4, psi3V0[fV0KindForBkg]);
        } break;
    }
    /// === e === determine background fit function   ###############################################

    
    // fFitModulation->SetParameter(0, fLocalRho->GetVal());
    _tempSwap.Fit(fFitModulation, "N0Q", "", lowBound, upBound);
    Double_t tempV2 = -999.9;
    Double_t tempV3 = -999.9;
    tempV2 = fFitModulation->GetParameter(1);
    tempV3 = fFitModulation->GetParameter(3);
    
    if(fEPQA){
        groupName="EventPlane";
        histName = TString::Format("%s/v2", groupName.Data());
        fHistManager.FillProfile(histName, fCent, tempV2);
        histName = TString::Format("%s/v3", groupName.Data());
        fHistManager.FillProfile(histName, fCent, tempV3);
    }

    // AliAnalysisTaskJetV2 ===========================================================================
    
    if(0){
        TCanvas *cBkgRhoFit = new TCanvas("cBkgRhoFit", "cBkgRhoFit", 1000, 700);
        
        TH1F* hBkgTracks_Event = (TH1F*) _tempSwap.Clone("hnew");
        hBkgTracks_Event->SetName("hBkgTracks");
        hBkgTracks_Event->Draw("E");
        _tempSwap.Fit(fFitModulation, "N0Q");
        fFitModulation->SetLineColor(632);
        fFitModulation->SetLineWidth(5); 
        fFitModulation->Draw("same");

        TF1* rho0Fit = new TF1("rho0Fit", "[0]", 0.0, TMath::TwoPi());
        rho0Fit->SetParameter(0, fFitModulation->GetParameter(0));
        rho0Fit->SetLineColor(12);
        rho0Fit->SetLineStyle(9);
        rho0Fit->SetLineWidth(5); 
        rho0Fit->Draw("same");

        TF1* rhoFitV2Com = new TF1("rhoFitV2Com", "[0]*(1.+2.*([1]*TMath::Cos(2.*(x-[2]))))", 0.0, TMath::TwoPi());
        rhoFitV2Com->SetParameter(0, fFitModulation->GetParameter(0));
        rhoFitV2Com->SetParameter(1, fFitModulation->GetParameter(1));//v2
        rhoFitV2Com->SetParameter(2, fFitModulation->GetParameter(2));//psi2
        rhoFitV2Com->SetLineColor(867);
        rhoFitV2Com->SetLineStyle(9);
        rhoFitV2Com->SetLineWidth(5); 
        rhoFitV2Com->Draw("same");

        TF1* rhoFitV3Com = new TF1("rhoFitV3Com", "[0]*(1.+2.*([1]*TMath::Cos(3.*(x-[2]))))", 0.0, TMath::TwoPi());
        rhoFitV3Com->SetParameter(0, fFitModulation->GetParameter(0));
        rhoFitV3Com->SetParameter(1, fFitModulation->GetParameter(3));//v3
        rhoFitV3Com->SetParameter(2, fFitModulation->GetParameter(4));//psi3
        rhoFitV3Com->SetLineColor(418);
        rhoFitV3Com->SetLineStyle(9);
        rhoFitV3Com->SetLineWidth(5); 
        rhoFitV3Com->Draw("same");

        histName = "checkOutput/cBkgRhoFit_Cent" + std::to_string(fCent) +".root";
        cBkgRhoFit->SaveAs(histName);


        delete cBkgRhoFit;
        delete hBkgTracks_Event;
    }
    
    fLocalRho->SetLocalRho(fFitModulation);
    BkgFitEvaluation(baseJetRho, &_tempSwap, fFitModulation);
    
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
    Int_t iCentBin = GetCentBin(); // fCentBin

    while ((partCont = static_cast<AliParticleContainer*>(next()))) {
        
        groupName = partCont->GetName();
        // counting events
        // histName = TString::Format("HistnEvents");
        // fHistManager.FillTH1(histName, 0.5);
        
        UInt_t count = 0;
        for(auto part : partCont->accepted()) {
            if (!part) continue;
            count++;
            
            histName = TString::Format("%s/hTrackPt_%d", groupName.Data(), iCentBin);
            fHistManager.FillTH1(histName, part->Pt());

            histName = TString::Format("%s/hTrackPhi_%d", groupName.Data(), iCentBin);
            fHistManager.FillTH1(histName, part->Phi());

            histName = TString::Format("%s/hTrackEta_%d", groupName.Data(), iCentBin);
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

        histName = TString::Format("%s/hNTracks_%d", groupName.Data(), iCentBin);
        fHistManager.FillTH1(histName, count);
        
        
    }

    histName = "fHistSumNTracks";
    fHistManager.FillTH1(histName, sumAcceptedTracks);
    
}


/**
 * This function performs a loop over the reconstructed jets
 * in the current event and fills the relevant histograms.
 */
void AliAnalysisTaskRawJetWithEP::DoJetLoop()
{
    if(fMesLev==3)std::cout << "=== Start DoJetLoop()  =====================" << std::endl;
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
        if (jetCont->GetRhoParameter()) {
            rhoVal = jetCont->GetRhoVal();
            histName = TString::Format("%s/hRhoVsCent", RhoGroupName.Data());
            fHistManager.FillTH2(histName.Data(), fCent, rhoVal);

            MeasureBkg(rhoVal);
        }
        
        Double_t leadingJetEta = -999.;
        Double_t leadingJetPhi = -999.;
        Double_t leadingJetPt  = -999.;
        
        for(auto jet : jetCont->accepted()) {
            if (!jet) continue;
            count++;

            Int_t iCentBin = GetCentBin(); //fCentBin

            histName = TString::Format("%s/hNJets_%d", GenGroupName.Data(), iCentBin);
            fHistManager.FillTH1(histName, count);
            histName = TString::Format("%s/hJetArea_%d", GenGroupName.Data(), iCentBin);
            fHistManager.FillTH1(histName, jet->Area());
            histName = TString::Format("%s/hJetPhi_%d", GenGroupName.Data(), iCentBin);
            fHistManager.FillTH1(histName, jet->Phi());
            histName = TString::Format("%s/hJetEta_%d", GenGroupName.Data(), iCentBin);
            fHistManager.FillTH1(histName, jet->Eta());

            // Filling histos for angle relative to event plane
            Double_t phiMinusPsi2 = jet->Phi() - psi2V0[0];
            
            if (phiMinusPsi2 < 0.0) phiMinusPsi2 += TMath::TwoPi();
            histName = TString::Format("%s/hJetPhiMinusPsi2_%d", GenGroupName.Data(), iCentBin);
            fHistManager.FillTH1(histName, phiMinusPsi2);
            Double_t phiMinusPsi3 = jet->Phi() - psi3V0[0];
            if (phiMinusPsi3 < 0.0) phiMinusPsi3 += TMath::TwoPi();
            histName = TString::Format("%s/hJetPhiMinusPsi3_%d", GenGroupName.Data(), iCentBin);
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
            localRhoVal = fLocalRho->GetLocalValInEtaPhi(jet->Phi(), GetJetContainer()->GetJetRadius(), fLocalRho->GetVal());
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
            
            histName = TString::Format("%s/hJetRho_%d", RhoGroupName.Data(), iCentBin);
            fHistManager.FillTH1(histName, jetCont->GetRhoVal());
            histName = TString::Format("%s/hJetRhoLocal_%d", RhoGroupName.Data(), iCentBin);
            fHistManager.FillTH1(histName, localRhoVal); // trying out local rho val

            histName = TString::Format("%s/hJetLRhoVsAveRho_%d", RhoGroupName.Data(), iCentBin);
            fHistManager.FillTH2(histName, jetCont->GetRhoVal(), localRhoValScaled);
            histName = TString::Format("%s/hJetGRhoVsDeltaPhi_%d", RhoGroupName.Data(), iCentBin);
            fHistManager.FillTH2(histName, deltaPhiJetEP, jetCont->GetRhoVal());
            histName = TString::Format("%s/hJetLRhoVsDeltaPhi_%d", RhoGroupName.Data(), iCentBin);
            fHistManager.FillTH2(histName, deltaPhiJetEP, localRhoValScaled);

            histName = TString::Format("%s/hCentJetPtCorrAPerCone", RhoGroupName.Data());
            fHistManager.FillTH3(histName, fCent, jetPtCorrLocal, jet->Area()/(TMath::Pi()*jetR*jetR));
            histName = TString::Format("%s/hCentJetPtRawGRhoA", RhoGroupName.Data());
            fHistManager.FillTH3(histName, fCent, jet->Pt(), jetCont->GetRhoVal() * jet->Area());
            histName = TString::Format("%s/hCentJetPtRawLRhoA", RhoGroupName.Data());
            fHistManager.FillTH3(histName, fCent, jet->Pt(), localRhoValScaled * jet->Area());
            
            
            Double_t rcPt = 0., rcEta = 0., rcPhi = 0.;
            CalcRandomCone(rcPt, rcEta, rcPhi, leadingJetEta, leadingJetPhi, jetR);    
            Double_t rcLocalRhoValScaled = fLocalRho->GetLocalVal(rcPhi, jetR, fLocalRho->GetVal());
            Double_t deltaLoacalPt = rcPt - rcLocalRhoValScaled*jetR*jetR*TMath::Pi();
            Double_t deltaGlobalPt = rcPt - fLocalRho->GetVal()*jetR*jetR*TMath::Pi();
            // std::cout << "delta L pT, rcPt, lRho, gRho, A = " << deltaLoacalPt << ", "\
            //     << rcPt << ", " << rcLocalRhoValScaled << ", " << fLocalRho->GetVal() << ", "\
            //     << jetR*jetR*TMath::Pi() << std::endl;

            histName = TString::Format("%s/hDeltaPt_%d", GenGroupName.Data(), iCentBin);
            fHistManager.FillTH1(histName, deltaGlobalPt);
            histName = TString::Format("%s/hDeltaPt_Local_%d", GenGroupName.Data(), iCentBin);
            fHistManager.FillTH1(histName, deltaLoacalPt);
            Double_t tempDPhi = -1.;
            tempDPhi = rcPhi - psi2V0[0];
            if(tempDPhi>TMath::Pi()) tempDPhi -= TMath::Pi();
            if(tempDPhi<0) tempDPhi += TMath::Pi();
            histName = TString::Format("%s/hPhiVsDeltaPt_Global_%d", GenGroupName.Data(), iCentBin);
            fHistManager.FillTH2(histName, tempDPhi, deltaGlobalPt);
            histName = TString::Format("%s/hPhiVsDeltaPt_Local_%d", GenGroupName.Data(), iCentBin);
            fHistManager.FillTH2(histName, tempDPhi, deltaLoacalPt);
            
            //inclusive Jet
            histName = TString::Format("%s/hJetPt_%d", InclusiveGroupName.Data(), iCentBin);
            fHistManager.FillTH1(histName, jet->Pt());
            histName = TString::Format("%s/hJetCorrPt_%d", InclusiveGroupName.Data(), iCentBin);
            fHistManager.FillTH1(histName, jetPtCorr);
            
            Double_t jetEP2Angle = jet->Phi() - psi2V0[0];
            if(jetEP2Angle < 0.) jetEP2Angle += TMath::TwoPi();
            // if((jetEP2Angle>TMath::Pi()/2)&&(jetEP2Angle<=TMath::Pi()/2)){
            //     jetEP2Angle = TMath::Pi() - jetEP2Angle;
            // } else if((jetEP2Angle>TMath::Pi())&&(jetEP2Angle<=TMath::Pi()*3/2)){
            //     jetEP2Angle = jetEP2Angle - TMath::Pi();
            // } else if((jetEP2Angle>TMath::Pi()*3/2)&&(jetEP2Angle<=TMath::TwoPi())){
            //     jetEP2Angle = TMath::TwoPi() - jetEP2Angle;
            // }
            
            
            histName = TString::Format("%s/hAzAngleJetPt_%d", InclusiveGroupName.Data(), iCentBin);
            fHistManager.FillTH2(histName, jetEP2Angle, jet->Pt());
            histName = TString::Format("%s/hAzAngleJetCorrPt_%d", InclusiveGroupName.Data(), iCentBin);
            fHistManager.FillTH2(histName, jetEP2Angle, jetPtCorrLocal);

            Int_t v0Angle = (Int_t) jet->Phi() / (TMath::Pi()/4);
            Double_t tempV0Mult = V0MultForAngle[v0Angle];
            histName = TString::Format("%s/hMultJetCorrPt", InclusiveGroupName.Data());
            fHistManager.FillTH2(histName, tempV0Mult, jet->Pt());
            histName = TString::Format("%s/hMultJetCorrPt", InclusiveGroupName.Data());
            fHistManager.FillTH2(histName, tempV0Mult, jetPtCorrLocal);
            
            histName = TString::Format("%s/hJetPtVsEP2AngleVsCent", InclusiveGroupName.Data()); 
            fHistManager.FillTH3(histName, jet->Pt(), jetEP2Angle, fCent);
            histName = TString::Format("%s/hJetCorrPtVsEP2AngleVsCent", InclusiveGroupName.Data());
            fHistManager.FillTH3(histName, jetPtCorrLocal, jetEP2Angle, fCent);
            histName = TString::Format("%s/hLocalRhoAVsEP2AngleVsCent", InclusiveGroupName.Data());
            fHistManager.FillTH3(histName, localRhoValScaled * jet->Area(), jetEP2Angle, fCent);
            

            Double_t jetEP3Angle = jet->Phi() - psi3V0[0];
            if(jetEP3Angle < 0.) jetEP3Angle += TMath::TwoPi();
            // // move to 0-(2/3)pi
            // if((jetEP3Angle>TMath::TwoPi()/3)&&(jetEP3Angle<=TMath::TwoPi()*2/3)){
            //     jetEP3Angle = jetEP3Angle - TMath::TwoPi()/3;
            // } else if((jetEP3Angle>TMath::TwoPi()*2/3)&&(jetEP3Angle<TMath::TwoPi())){
            //     jetEP3Angle = jetEP3Angle - TMath::TwoPi()*2/3;
            // } 
            // // move to 0-(1/3)pi
            // if(jetEP3Angle > TMath::Pi()/3) jetEP3Angle = TMath::TwoPi()/3 - jetEP3Angle;

            histName = TString::Format("%s/hJetPtVsEP3AngleVsCent", InclusiveGroupName.Data()); 
            fHistManager.FillTH3(histName, jet->Pt(), jetEP3Angle, fCent);
            histName = TString::Format("%s/hJetCorrPtVsEP3AngleVsCent", InclusiveGroupName.Data());
            fHistManager.FillTH3(histName, jetPtCorrLocal, jetEP3Angle, fCent);
            histName = TString::Format("%s/hLocalRhoAVsEP3AngleVsCent", InclusiveGroupName.Data());
            fHistManager.FillTH3(histName, localRhoValScaled * jet->Area(), jetEP3Angle, fCent);
            
            //V2 In plane Jet
            if(fSepEP){
                if ((phiMinusPsi2 < TMath::Pi()/4) || (phiMinusPsi2 >= 7*TMath::Pi()/4)\
                    || (phiMinusPsi2 >= 3*TMath::Pi()/4 && phiMinusPsi2 < 5*TMath::Pi()/4)) {
                    histName = TString::Format("%s/hJetPt_%d", IPlaneGroupName.Data(), iCentBin);
                    fHistManager.FillTH1(histName, jet->Pt());
                    histName = TString::Format("%s/hJetCorrPtLocal_%d", IPlaneGroupName.Data(), iCentBin);
                    fHistManager.FillTH1(histName, jetPtCorrLocal);
                }
                else {
                    histName = TString::Format("%s/hJetPt_%d", OPlaneGroupName.Data(), iCentBin);
                    fHistManager.FillTH1(histName, jet->Pt());
                    histName = TString::Format("%s/hJetCorrPtLocal_%d", OPlaneGroupName.Data(), iCentBin);
                    fHistManager.FillTH1(histName, jetPtCorrLocal);
                }
            }
            
        }
    }
}

Bool_t AliAnalysisTaskRawJetWithEP::QnJEHandlarEPGet()
{
    ResetAODEvent();
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    SetAODEvent(fAOD); 
    if(!fIsOADBFileOpen || fCalibObjRun!=fRun) {
        fIsOADBFileOpen = OpenInfoCalbration();
        if(!fIsOADBFileOpen)
            return kFALSE;
        fCalibObjRun = fRun;
    }
    
    //== Q2 Vector ######################################## 
    Double_t harmonic = 2.;
    ComputeQvecV0(q2VecV0M, q2VecV0C, q2VecV0A, q2V0, V0Mult2, harmonic);
    ComputeQvecTpc(q2VecTpcM, q2VecTpcN, q2VecTpcP, q2Tpc, TpcMult2, harmonic);
    
    // == s == ComputeEventPlaneAngle ======================
    // Inisialize
    for(Int_t i = 0; i<3; i++){
        psi2V0[i] = -1;
        psi2Tpc[i] = -1;
    }
    
    psi2V0[0] = ComputeEventPlaneAngle(q2VecV0M, harmonic);
    psi2V0[1] = ComputeEventPlaneAngle(q2VecV0C, harmonic);
    psi2V0[2] = ComputeEventPlaneAngle(q2VecV0A, harmonic);
    
    psi2Tpc[0] = ComputeEventPlaneAngle(q2VecTpcM, harmonic);
    psi2Tpc[1] = ComputeEventPlaneAngle(q2VecTpcN, harmonic);
    psi2Tpc[2] = ComputeEventPlaneAngle(q2VecTpcP, harmonic);
    
    //== Q3 Vector ######################################## 
    harmonic = 3.;
    ComputeQvecV0(q3VecV0M, q3VecV0C, q3VecV0A, q3V0, V0Mult3, harmonic);
    ComputeQvecTpc(q3VecTpcM, q3VecTpcN, q3VecTpcP, q3Tpc, TpcMult3, harmonic);
    
    // Inisialize
    for(Int_t i = 0; i<3; i++){
        psi3V0[i] = -1;
        psi3Tpc[i] = -1;
    }
    

    psi3V0[0] = ComputeEventPlaneAngle(q3VecV0M, harmonic);
    psi3V0[1] = ComputeEventPlaneAngle(q3VecV0C, harmonic);
    psi3V0[2] = ComputeEventPlaneAngle(q3VecV0A, harmonic);

    psi3Tpc[0] = ComputeEventPlaneAngle(q3VecTpcM, harmonic);
    psi3Tpc[1] = ComputeEventPlaneAngle(q3VecTpcN, harmonic);
    psi3Tpc[2] = ComputeEventPlaneAngle(q3VecTpcP, harmonic);

    return kTRUE;
}


Bool_t AliAnalysisTaskRawJetWithEP::QnV0GainCalibration(){
    //V0 Channel Gains:
    fHCorrV0ChWeghts = (TH2F *)fCalibRefObjList->FindObject(Form("hWgtV0ChannelsvsVzRun%d",fRunNumber));

    AliAODVZERO* fAodV0 = dynamic_cast<AliAODVZERO*>(fAOD->GetVZEROData());
    AliVEvent *fVevent = dynamic_cast<AliVEvent*>(InputEvent());
    const AliVVertex *pointVtx = fVevent->GetPrimaryVertex();
    Double_t fVtxZ = -999;
    fVtxZ  = pointVtx->GetZ();
    
    Double_t fMultV0 = 0.;
    Int_t ibinV0 = 0;
    Double_t fSumMV0A = 0.;
    Double_t fSumMV0C = 0.;
    Double_t fSumMV0M = 0.;
    for(int i = 0; i < 2; i++){
        q2VecV0M[i] = 0.;
        q2VecV0C[i] = 0.;
        q2VecV0A[i] = 0.;
        q3VecV0M[i] = 0.;
        q3VecV0C[i] = 0.;
        q3VecV0A[i] = 0.;
    }
    for(int i = 0; i < 3; i++){
        V0Mult2[i] = 0.;
    }


    for(int iV0 = 0; iV0 < 64; iV0++) { //0-31 is V0C, 32-63 VOA
        fMultV0 = fAodV0->GetMultiplicity(iV0);
        Double_t fV0chGain = 1.;
        Double_t fPhiV0  = TMath::PiOver4()*(0.5 + iV0 % 8);

        /// V0 Channel Gain Correction:
        if(fHCorrV0ChWeghts){
            ibinV0    = fHCorrV0ChWeghts->FindBin(fVtxZ,iV0);
            fV0chGain = fHCorrV0ChWeghts->GetBinContent(ibinV0);
        }
        fMultV0 = fMultV0*fV0chGain;
        
        q2VecV0M[0] += TMath::Cos(2*fPhiV0) * fMultV0;
        q2VecV0M[1] += TMath::Sin(2*fPhiV0) * fMultV0;
        q3VecV0M[0] += TMath::Cos(3*fPhiV0) * fMultV0;
        q3VecV0M[1] += TMath::Sin(3*fPhiV0) * fMultV0;
        V0Mult2[0] += fMultV0;

        if(iV0 < 32){
            q2VecV0C[0] += TMath::Cos(2*fPhiV0) * fMultV0;
            q2VecV0C[1] += TMath::Sin(2*fPhiV0) * fMultV0;
            q3VecV0C[0] += TMath::Cos(3*fPhiV0) * fMultV0;
            q3VecV0C[1] += TMath::Sin(3*fPhiV0) * fMultV0;
            V0Mult2[1] += fMultV0;
        }
        else if(iV0 >= 32){
            q2VecV0A[0] += TMath::Cos(2*fPhiV0) * fMultV0;
            q2VecV0A[1] += TMath::Sin(2*fPhiV0) * fMultV0;
            q3VecV0A[0] += TMath::Cos(3*fPhiV0) * fMultV0;
            q3VecV0A[1] += TMath::Sin(3*fPhiV0) * fMultV0;
            V0Mult2[2] += fMultV0;
        }
    
    }///V0 Channel loop

    /// Now the q vectors:
    if(V0Mult2[1]<=1e-4 || V0Mult2[2]<=1e-4){
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
    return kTRUE;  

}


Bool_t AliAnalysisTaskRawJetWithEP::QnRecenteringCalibration(){
    //Get V0A, V0C <Q> Vectors:
    fHCorrQ2xV0C = (TH1D *)fCalibRefObjList->FindObject(Form("fHisAvgQ2xvsCentV0CRun%d",fRunNumber));
    fHCorrQ2yV0C = (TH1D *)fCalibRefObjList->FindObject(Form("fHisAvgQ2yvsCentV0CRun%d",fRunNumber));    
    fHCorrQ2xV0A = (TH1D *)fCalibRefObjList->FindObject(Form("fHisAvgQ2xvsCentV0ARun%d",fRunNumber));
    fHCorrQ2yV0A = (TH1D *)fCalibRefObjList->FindObject(Form("fHisAvgQ2yvsCentV0ARun%d",fRunNumber));
    
    fHCorrQ3xV0C = (TH1D *)fCalibRefObjList->FindObject(Form("fHisAvgQ3xvsCentV0CRun%d",fRunNumber));
    fHCorrQ3yV0C = (TH1D *)fCalibRefObjList->FindObject(Form("fHisAvgQ3yvsCentV0CRun%d",fRunNumber));    
    fHCorrQ3xV0A = (TH1D *)fCalibRefObjList->FindObject(Form("fHisAvgQ3xvsCentV0ARun%d",fRunNumber));
    fHCorrQ3yV0A = (TH1D *)fCalibRefObjList->FindObject(Form("fHisAvgQ3yvsCentV0ARun%d",fRunNumber));    

    Int_t icentbin = 0;
    Double_t avgqx=0,avgqy=0; 

    if(fHCorrQ2xV0C && fHCorrQ2yV0C){
        icentbin = fHCorrQ2xV0C->FindBin(fCent);
        avgqx = fHCorrQ2xV0C->GetBinContent(icentbin);
        avgqy = fHCorrQ2yV0C->GetBinContent(icentbin);
        q2VecV0C[0] -= avgqx;
        q2VecV0C[1] -= avgqy;
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
    
    if(V0Mult2[0] != 0.){
        q2VecV0M[0] /= V0Mult2[0];
        q2VecV0M[1] /= V0Mult2[0];
        q3VecV0M[0] /= V0Mult2[0];
        q3VecV0M[1] /= V0Mult2[0];
    }

    if(V0Mult2[1] != 0.){
        q2VecV0C[0] /= V0Mult2[1];
        q2VecV0C[1] /= V0Mult2[1];
        q3VecV0C[0] /= V0Mult2[1];
        q3VecV0C[1] /= V0Mult2[1];
    }

    if(V0Mult2[2] != 0.){
        q2VecV0A[0] /= V0Mult2[2];
        q2VecV0A[1] /= V0Mult2[2];
        q3VecV0A[0] /= V0Mult2[2];
        q3VecV0A[1] /= V0Mult2[2];
    }

    q2V0[0] = TMath::Sqrt(q2VecV0M[0]*q2VecV0M[0]+q2VecV0M[1]*q2VecV0M[1])*TMath::Sqrt(V0Mult2[0]);
    q2V0[1] = TMath::Sqrt(q2VecV0C[0]*q2VecV0C[0]+q2VecV0C[1]*q2VecV0C[1])*TMath::Sqrt(V0Mult2[1]);
    q2V0[2] = TMath::Sqrt(q2VecV0A[0]*q2VecV0A[0]+q2VecV0A[1]*q2VecV0A[1])*TMath::Sqrt(V0Mult2[2]);

    q3V0[0] = TMath::Sqrt(q3VecV0M[0]*q3VecV0M[0]+q3VecV0M[1]*q3VecV0M[1])*TMath::Sqrt(V0Mult2[0]);
    q3V0[1] = TMath::Sqrt(q3VecV0C[0]*q3VecV0C[0]+q3VecV0C[1]*q3VecV0C[1])*TMath::Sqrt(V0Mult2[1]);
    q3V0[2] = TMath::Sqrt(q3VecV0A[0]*q3VecV0A[0]+q3VecV0A[1]*q3VecV0A[1])*TMath::Sqrt(V0Mult2[2]);

    return kTRUE;
}


Bool_t AliAnalysisTaskRawJetWithEP::QnCalcWOCalib(){
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
void AliAnalysisTaskRawJetWithEP::BkgFitEvaluation(Double_t baseJetRho, TH1F* hBkgTracks, TF1* fFitModulation)
{
    // the quality of the fit is evaluated from 1 - the cdf of the chi square distribution
    // three methods are available, all with their drawbacks. 
    // all are stored, one is selected to do the cut
    Int_t numOfFreePara = 2; //v2, v3
    switch (fFitModulationType)  {
        case kNoFit : { numOfFreePara = 0; } break;
        case kV2 : { numOfFreePara = 1;} break;
        case kCombined: { numOfFreePara = 2; } break;
        default : { numOfFreePara = 2;} break;
    }
    
    Int_t NDF = 1;
    NDF = (Int_t)fFitModulation->GetXaxis()->GetNbins() - numOfFreePara;
    if(NDF == 0 || (float)NDF <= 0.) return;
    
    Double_t ChiSqr = 999.;
    Double_t CDF = 1.;
    Double_t CDFROOT = 1.;
    ChiSqr = ChiSquare(*hBkgTracks, fFitModulation);
    CDF = 1. - ChiSquareCDF(NDF, ChiSqr);
    CDFROOT = 1.-ChiSquareCDF(NDF, fFitModulation->GetChisquare());
    
    // CDF = 1. - ChiSquarePDF(NDF, ChiSqr);
    // CDFROOT = 1.-ChiSquarePDF(NDF, fFitModulation->GetChisquare());
    
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
    // tempV2Fit->SetParameter(0, fRho->GetVal());
    tempV2Fit->SetParameter(0, baseJetRho);
    tempV2Fit->FixParameter(2, psi2V0[0]);
    hBkgTracks->Fit(tempV2Fit, "N0Q");
    numOfFreePara = 2; //v2
    NDF = tempV2Fit->GetXaxis()->GetNbins() - numOfFreePara;
    if(NDF == 0 || (float)NDF <= 0.) return;
    ChiSqr = 999.;
    ChiSqr = ChiSquare(*hBkgTracks, tempV2Fit);
    // CDF = 1. - ChiSquareCDF(NDF, ChiSqr);
    CDF = 1. - ChiSquarePDF(NDF, ChiSqr);
    CDFROOT = 1. - ChiSquareCDF(NDF, tempV2Fit->GetChisquare());
    // CDFROOT = 1. - ChiSquarePDF(NDF, tempV2Fit->GetChisquare());

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
    // tempGlobalFit->FixParameter(0, fRho->GetVal());
    tempGlobalFit->FixParameter(0, baseJetRho);
    hBkgTracks->Fit(tempGlobalFit, "N0Q");
    numOfFreePara = 1; //v2
    NDF = tempGlobalFit->GetXaxis()->GetNbins() - numOfFreePara;
    if(NDF == 0 || (float)NDF <= 0.) return;
    ChiSqr = 999.;
    ChiSqr = ChiSquare(*hBkgTracks, tempGlobalFit);
    CDF = 1. - ChiSquareCDF(NDF, ChiSqr);
    CDFROOT = 1.-ChiSquareCDF(NDF, tempGlobalFit->GetChisquare());
    // CDF = 1. - ChiSquarePDF(NDF, ChiSqr);
    // CDFROOT = 1.-ChiSquarePDF(NDF, tempGlobalFit->GetChisquare());

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
//
//   ////==========> LHC18q/r PileUp Removal Functions: ---- Do not Remove them !!! -----
//   Double_t parV0[8] = {43.8011, 0.822574, 8.49794e-02, 1.34217e+02, 7.09023e+00, 4.99720e-02, -4.99051e-04, 1.55864e-06};
//   fV0CutPU  = new TF1("fV0CutPU", "[0]+[1]*x - 6.*[2]*([3] + [4]*sqrt(x) + [5]*x + [6]*x*sqrt(x) + [7]*x*x)", 0, 100000);
//   fV0CutPU->SetParameters(parV0);
//   fSPDCutPU = new TF1("fSPDCutPU", "400. + 4.*x", 0, 10000);
//
//   Double_t parFB32[8] = {2093.36, -66.425, 0.728932, -0.0027611, 1.01801e+02, -5.23083e+00, -1.03792e+00, 5.70399e-03};
//   fMultCutPU = new TF1("fMultCutPU", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 6.*([4]+[5]*sqrt(x)+[6]*x+[7]*x*x)", 0, 90);
//   fMultCutPU->SetParameters(parFB32);
//
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


// Creates the QnVector Handlers. Needs to be run in add task so that calibration
/**
 * This function is called once at the end of the analysis.
 */
void AliAnalysisTaskRawJetWithEP::Terminate(Option_t *) 
{
}


/// ================================================================================================
//________________________________________________________________
bool AliAnalysisTaskRawJetWithEP::SetAODEvent(AliAODEvent* event)
{
    if(!event) {
        AliWarning("Event not found!");
        return false;
    }
    fRun = event->GetRunNumber();
    const AliAODVertex* trkVtx = dynamic_cast<AliAODVertex*>(event->GetPrimaryVertex());
    if (!trkVtx || trkVtx->GetNContributors()<=0) {
        AliWarning("Primary vertex not found!");
        return false;
    }
    else fZvtx = trkVtx->GetZ();

    AliMultSelection* MultSelection = dynamic_cast<AliMultSelection*>(event->FindListObject("MultSelection"));
    if(!MultSelection) {
        AliWarning("AliMultSelection object not found!");
        return false;
    }    
    else fCentrality = MultSelection->GetMultiplicityPercentile("V0M");

    fV0 = dynamic_cast<AliAODVZERO*>(event->GetVZEROData());
    if(!fV0) {
        AliWarning("V0 info not found!");
        return false;
    }

    //only if everything ok get event
    fAOD = event;

    return true;
}

//________________________________________________________________
void AliAnalysisTaskRawJetWithEP::ResetAODEvent() 
{
    fRun            = -9999;
    fCentrality     = -9999.;
    fZvtx           = -9999.;
    fV0             = nullptr;

    fUsedTrackPosIDs.ResetAllBits();
    fUsedTrackNegIDs.ResetAllBits();
}


//________________________________________________________________
void AliAnalysisTaskRawJetWithEP::GetQnVec(double QnVecFull[2], double QnVecA[2], double QnVecC[2], Double_t QnNorm[3], Double_t Multi[3])
{
    switch(fNormMethod) {
        case 0: //kQoverQlength
            for(int iComp=0; iComp<2; iComp++) {
                QnNorm[0] > 0 ? QnVecFull[iComp] = QnVecFull[iComp] / QnNorm[0] : QnVecFull[iComp] = QnVecFull[iComp];
                QnNorm[2] > 0 ? QnVecA[iComp]  = QnVecA[iComp] / QnNorm[2] : QnVecA[iComp] = QnVecA[iComp];
                QnNorm[1] > 0 ? QnVecC[iComp]  = QnVecC[iComp] / QnNorm[1] : QnVecC[iComp] = QnVecC[iComp];
            }
        break;
        case 1: //kQoverM
            for(int iComp=0; iComp<2; iComp++) {
                Multi[0] > 0 ? QnVecFull[iComp] = QnVecFull[iComp] / Multi[0] : QnVecFull[iComp] = QnVecFull[iComp];
                Multi[2] > 0 ? QnVecA[iComp]  = QnVecA[iComp] / Multi[2] : QnVecA[iComp] = QnVecA[iComp];
                Multi[1] > 0 ? QnVecC[iComp]  = QnVecC[iComp] / Multi[1] : QnVecC[iComp] = QnVecC[iComp];
            }
        break;
        case 2: //kQoverSqrtM
            for(int iComp=0; iComp<2; iComp++) {
                Multi[0] > 0 ? QnVecFull[iComp] = QnVecFull[iComp] / TMath::Sqrt(Multi[0]) : QnVecFull[iComp] = QnVecFull[iComp];
                Multi[2] > 0 ? QnVecA[iComp]  = QnVecA[iComp] / TMath::Sqrt(Multi[2]) : QnVecA[iComp] = QnVecA[iComp];
                Multi[1] > 0 ? QnVecC[iComp]  = QnVecC[iComp] / TMath::Sqrt(Multi[1]) : QnVecC[iComp] = QnVecC[iComp];
            }
        break;
        case 3: //kNone
            for(int iComp=0; iComp<2; iComp++) {
                QnVecFull[iComp] = QnVecFull[iComp];
                QnVecA[iComp]  = QnVecA[iComp];
                QnVecC[iComp]  = QnVecC[iComp];
            }
        break;
    }
}


//________________________________________________________________
void AliAnalysisTaskRawJetWithEP::Getqn(Double_t Qn[3], Double_t QnNorm[3], Double_t Multi[3])
{
    Multi[0] > 0 ? Qn[0] = QnNorm[0] / TMath::Sqrt(Multi[0]) : Qn[0] = 0;
    Multi[1] > 0 ? Qn[1] = QnNorm[1] / TMath::Sqrt(Multi[1]) : Qn[1] = 0;
    Multi[2] > 0 ? Qn[2] = QnNorm[2] / TMath::Sqrt(Multi[2]) : Qn[2] = 0;
}


//________________________________________________________________
void AliAnalysisTaskRawJetWithEP::EnablePhiDistrHistos() 
{
    fEnablePhiDistrHistos=true;
    TString histonames[2] = {"TPCPosEta","TPCNegEta"};
    for(int iHisto=0; iHisto<2; iHisto++) {
        if(fPhiVsCentrTPC[iHisto]) {
            delete fPhiVsCentrTPC[iHisto];
            fPhiVsCentrTPC[iHisto]=nullptr;
        }
        fPhiVsCentrTPC[iHisto] = new TH2F(Form("fPhiVsCentrTPC%s",histonames[iHisto].Data()),";centrality (%);#varphi;Entries",10,0.,100.,180,0.,2*TMath::Pi());
    }
}

//__________________________________________________________
void AliAnalysisTaskRawJetWithEP::ComputeQvecTpc(Double_t QnVecTpcM[2],Double_t QnVecTpcN[2],Double_t QnVecTpcP[2], Double_t QnNorm[3], Double_t Multi[3], unsigned int harmonic) 
{
    Int_t centbin = GetCentBin();
    
    //initialise Q vectors
    for(int iComp=0; iComp<2; iComp++) {
        QnVecTpcM[iComp] = 0.;
        QnVecTpcP[iComp]  = 0.;
        QnVecTpcN[iComp]  = 0.;
    }
    for(int i=0; i<3; i++) {
        QnNorm[3] = 1.;
        Multi[3] = 0.;
    }
    
    fUsedTrackPosIDs.ResetAllBits();
    fUsedTrackNegIDs.ResetAllBits();
    
    //reset phi distributions
    if(fEnablePhiDistrHistos) {
        fPhiVsCentrTPC[0]->Reset();
        fPhiVsCentrTPC[1]->Reset();
    }

    int nTracks=fAOD->GetNumberOfTracks();
    for(int iTrack=0; iTrack<nTracks; iTrack++) {
        AliAODTrack* track=dynamic_cast<AliAODTrack*>(fAOD->GetTrack(iTrack));
        double pt=track->Pt();
        double eta=track->Eta();
        if(!track || !IsTrackSelected(track)) continue;
        double phi=track->Phi();
        double pseudoRand = pt*1000.-(long)(pt*1000);
        if(pseudoRand>fFractionOfTracksForQnTPC) continue;
        if(track->GetID()>=0) fUsedTrackPosIDs.SetBitNumber(track->GetID());
        else fUsedTrackNegIDs.SetBitNumber(TMath::Abs(track->GetID()));
        double qx = TMath::Cos(harmonic*phi);
        double qy = TMath::Sin(harmonic*phi);
        double weight = 1.;
        if(eta>0) {
            if(fWeightsTPCPosEta[centbin]) {
                int phibin = fWeightsTPCPosEta[centbin]->GetXaxis()->FindBin(phi);
                weight = 1./fWeightsTPCPosEta[centbin]->GetBinContent(phibin);
            }
            QnVecTpcM[0] += weight*qx;
            QnVecTpcM[1] += weight*qy;
            Multi[0]     += weight;
            QnVecTpcP[0]  += weight*qx;
            QnVecTpcP[1]  += weight*qy;
            Multi[2]      += weight;
            if(fEnablePhiDistrHistos && fPhiVsCentrTPC[0])
                fPhiVsCentrTPC[0]->Fill(fCentrality,phi);
        }
        else {
            if(fWeightsTPCNegEta[centbin]) {
                int phibin = fWeightsTPCNegEta[centbin]->GetXaxis()->FindBin(phi);
                weight = 1./fWeightsTPCNegEta[centbin]->GetBinContent(phibin);
            }
            QnVecTpcM[0] += weight*qx;
            QnVecTpcM[1] += weight*qy;
            Multi[0]     += weight;
            QnVecTpcN[0]  += weight*qx;
            QnVecTpcN[1]  += weight*qy;
            Multi[1]     += weight;
            if(fEnablePhiDistrHistos && fPhiVsCentrTPC[1])
                fPhiVsCentrTPC[1]->Fill(fCentrality,phi);
        }
    }
    
    QnNorm[0] = TMath::Sqrt(QnVecTpcM[0]*QnVecTpcM[0]+QnVecTpcM[1]*QnVecTpcM[1]);
    QnNorm[2]  = TMath::Sqrt(QnVecTpcP[0]*QnVecTpcP[0]+QnVecTpcP[1]*QnVecTpcP[1]);
    QnNorm[1]  = TMath::Sqrt(QnVecTpcN[0]*QnVecTpcN[0]+QnVecTpcN[1]*QnVecTpcN[1]);
}

//__________________________________________________________
void AliAnalysisTaskRawJetWithEP::ComputeQvecV0(Double_t QnVecV0M[2],Double_t QnVecV0C[2],Double_t QnVecV0A[2], Double_t QnNorm[3], Double_t Multi[3], unsigned int harmonic)
{
    TString histName;
    TString groupName;
    groupName="EventPlane";
    if(fMesLev==3)std::cout << "ComputeQvecV0(): Start Initialize Qn values  ===============" << std::endl;
    
    if(harmonic == 2){
        for(int i=0; i<8; i++) V0MultForAngle[i] = 0;
    }
    
    
    for(int iComp=0; iComp<2; iComp++) {
        QnVecV0M[iComp] = 0.;
        QnVecV0A[iComp] = 0.;
        QnVecV0C[iComp] = 0.;
    }

    for(int i=0; i<3; i++) {
        QnNorm[i] = 1.;
        Multi[i] = 0.;
    }
    

    short zvtxbin = GetVertexZbin();
    for (int iCh = 0; iCh < 64; iCh++) {
        
        double phiCh = TMath::PiOver4()*(0.5 + iCh % 8);
        double multv0 = fV0->GetMultiplicity(iCh);
        
        
        if(harmonic == 2){
            int iV0ChAngle  = iCh % 8;
            V0MultForAngle[iV0ChAngle] += multv0;
        }
        
        if (iCh < 32) { // V0C side
            double multCorC = -10;
            if (iCh < 8)
                multCorC = multv0/fHistMultV0->GetBinContent(iCh+1)*fHistMultV0->GetBinContent(1);
            else if (iCh >= 8 && iCh < 16)
                multCorC = multv0/fHistMultV0->GetBinContent(iCh+1)*fHistMultV0->GetBinContent(9);
            else if (iCh >= 16 && iCh < 24)
                multCorC = multv0/fHistMultV0->GetBinContent(iCh+1)*fHistMultV0->GetBinContent(17);
            else if (iCh >= 24 && iCh < 32)
                multCorC = multv0/fHistMultV0->GetBinContent(iCh+1)*fHistMultV0->GetBinContent(25);
            if (multCorC < 0) {
                AliWarning("Problem with multiplicity in V0C");
                continue;
            }
            
            QnVecV0C[0] += TMath::Cos(harmonic*phiCh) * multCorC;
            QnVecV0C[1] += TMath::Sin(harmonic*phiCh) * multCorC;
            QnVecV0M[0] += TMath::Cos(harmonic*phiCh) * multCorC;
            QnVecV0M[1] += TMath::Sin(harmonic*phiCh) * multCorC;

            Multi[1] += multCorC;  
            Multi[0] += multCorC;  
        } 
        else { // V0A side
            double multCorA = -10;
            
            if (iCh >= 32 && iCh < 40)
                multCorA = multv0/fHistMultV0->GetBinContent(iCh+1)*fHistMultV0->GetBinContent(33);
            else if (iCh >= 40 && iCh < 48)
                multCorA = multv0/fHistMultV0->GetBinContent(iCh+1)*fHistMultV0->GetBinContent(41);
            else if (iCh >= 48 && iCh < 56)
                multCorA = multv0/fHistMultV0->GetBinContent(iCh+1)*fHistMultV0->GetBinContent(49);
            else if (iCh >= 56 && iCh < 64)
                multCorA = multv0/fHistMultV0->GetBinContent(iCh+1)*fHistMultV0->GetBinContent(57);
            
            if (multCorA < 0) {
                AliWarning("Problem with multiplicity in V0A");
                continue;
            }
            
            QnVecV0A[0] += TMath::Cos(harmonic*phiCh) * multCorA;
            QnVecV0A[1] += TMath::Sin(harmonic*phiCh) * multCorA;
            QnVecV0M[0] += TMath::Cos(harmonic*phiCh) * multCorA;
            QnVecV0M[1] += TMath::Sin(harmonic*phiCh) * multCorA;

            Multi[2] += multCorA;  
            Multi[0] += multCorA;              
        }
    }

    // int iCentBin = static_cast<int>(fCentrality)+1;
    Int_t iCentBin = GetCentBin();
    
    //only recentering and not width equalisation to preserve multiplicity dependence (needed for qn)
    if(harmonic == 2){
        
        QnVecV0A[0] = (QnVecV0A[0] - fQx2mV0A[zvtxbin]->GetBinContent(iCentBin));///fQx2sV0A[zvtxbin]->GetBinContent(iCentBin);
        QnVecV0A[1] = (QnVecV0A[1] - fQy2mV0A[zvtxbin]->GetBinContent(iCentBin));///fQy2sV0A[zvtxbin]->GetBinContent(iCentBin);
        QnVecV0C[0] = (QnVecV0C[0] - fQx2mV0C[zvtxbin]->GetBinContent(iCentBin));///fQx2sV0C[zvtxbin]->GetBinContent(iCentBin);
        QnVecV0C[1] = (QnVecV0C[1] - fQy2mV0C[zvtxbin]->GetBinContent(iCentBin));///fQy2sV0C[zvtxbin]->GetBinContent(iCentBin);
        
        // Double_t avgqxC = fQx2mV0C[zvtxbin]->GetBinContent(iCentBin);
        // Double_t avgqyC = fQy2mV0C[zvtxbin]->GetBinContent(iCentBin);
        // Double_t avgqxA = fQx2mV0A[zvtxbin]->GetBinContent(iCentBin);
        // Double_t avgqyA = fQy2mV0A[zvtxbin]->GetBinContent(iCentBin);
        
    }
    else if(harmonic == 3){
        
        QnVecV0A[0] = (QnVecV0A[0] - fQx3mV0A[zvtxbin]->GetBinContent(iCentBin));///fQx2sV0A[zvtxbin]->GetBinContent(iCentBin);
        QnVecV0A[1] = (QnVecV0A[1] - fQy3mV0A[zvtxbin]->GetBinContent(iCentBin));///fQy2sV0A[zvtxbin]->GetBinContent(iCentBin);
        QnVecV0C[0] = (QnVecV0C[0] - fQx3mV0C[zvtxbin]->GetBinContent(iCentBin));///fQx2sV0C[zvtxbin]->GetBinContent(iCentBin);
        QnVecV0C[1] = (QnVecV0C[1] - fQy3mV0C[zvtxbin]->GetBinContent(iCentBin));///fQy2sV0C[zvtxbin]->GetBinContent(iCentBin);
    }

    
    QnNorm[0] = TMath::Sqrt(QnVecV0M[0]*QnVecV0M[0]+QnVecV0M[1]*QnVecV0M[1]);
    QnNorm[1] = TMath::Sqrt(QnVecV0C[0]*QnVecV0C[0]+QnVecV0C[1]*QnVecV0C[1]);
    QnNorm[2] = TMath::Sqrt(QnVecV0A[0]*QnVecV0A[0]+QnVecV0A[1]*QnVecV0A[1]);
    
}

//__________________________________________________________
short AliAnalysisTaskRawJetWithEP::GetVertexZbin() const
{
    if(!fV0CalibZvtxDiff)
        return 0; //if it is not Zvtx differential, always first bin

    short zvtxbin = -10;
    
    if (fZvtx >= -10. && fZvtx < -8.)     zvtxbin = 0;
    else if (fZvtx >= -8. && fZvtx < -6.) zvtxbin = 1;
    else if (fZvtx >= -6. && fZvtx < -4.) zvtxbin = 2;
    else if (fZvtx >= -4. && fZvtx < -3.) zvtxbin = 3;
    else if (fZvtx >= -3. && fZvtx < -2.) zvtxbin = 4;
    else if (fZvtx >= -2. && fZvtx < -1.) zvtxbin = 5;
    else if (fZvtx >= -1. && fZvtx < 0)   zvtxbin = 6;
    else if (fZvtx >= 0 && fZvtx < 1.)    zvtxbin = 7;
    else if (fZvtx >= 1. && fZvtx < 2.)   zvtxbin = 8;
    else if (fZvtx >= 2. && fZvtx < 3.)   zvtxbin = 9;
    else if (fZvtx >= 3. && fZvtx < 4.)   zvtxbin = 10;
    else if (fZvtx >= 4. && fZvtx < 6.)   zvtxbin = 11;
    else if (fZvtx >= 6. && fZvtx < 8.)   zvtxbin = 12;
    else if (fZvtx >= 8. && fZvtx <= 10.) zvtxbin = 13;
    
    return zvtxbin;
}

//__________________________________________________________
Int_t AliAnalysisTaskRawJetWithEP::GetCentBin()
{
    Int_t centbin = 10;
    
    if (fCent >= 0. && fCent < 5.)        centbin = 0;
    else if (fCent >= 5. && fCent < 10.)  centbin = 1;
    else if (fCent >= 10. && fCent < 20.) centbin = 2;
    else if (fCent >= 20. && fCent < 30.) centbin = 3;
    else if (fCent >= 30. && fCent < 40.) centbin = 4;
    else if (fCent >= 40. && fCent < 50.) centbin = 5;
    else if (fCent >= 50. && fCent < 60.) centbin = 6;
    else if (fCent >= 60. && fCent < 70.) centbin = 7;
    else if (fCent >= 70. && fCent < 80.) centbin = 8;
    else if (fCent >= 80. && fCent < 90.) centbin = 9;

    return centbin;
}

//__________________________________________________________
bool AliAnalysisTaskRawJetWithEP::IsTrackSelected(AliAODTrack* track) {

    if(fRun>=244918 && fRun<=246994) {//PbPb2015
        if(fCalibType==kQnCalib) {
            if(!track->TestFilterBit(BIT(8)) && !track->TestFilterBit(BIT(9)))
                return false;
            if(track->Pt()<fPtMinTPC || track->Pt()>fPtMaxTPC || TMath::Abs(track->Eta())>fEtaMaxTPC || TMath::Abs(track->Eta())<fEtaMinTPC)
                return false;
        }
        else if(fCalibType==kQnFrameworkCalib) {
            if(!track->TestFilterBit(BIT(8)) && !track->TestFilterBit(BIT(9)))
                return false;
            if(track->Pt()<0.2 || track->Pt()>5 || TMath::Abs(track->Eta())>0.8 || TMath::Abs(track->Eta())<0.)
                return false;
        }
    }
    else if(fRun>=295581 && fRun<=297317) { //PbPb2018
        if(fCalibType==kQnCalib) {
            if(!track->TestFilterBit(BIT(8)) && !track->TestFilterBit(BIT(9)))
                return false;
            if(track->Pt()<fPtMinTPC || track->Pt()>fPtMaxTPC || TMath::Abs(track->Eta())>fEtaMaxTPC || TMath::Abs(track->Eta())<fEtaMinTPC)
                return false;
        }
        else if(fCalibType==kQnFrameworkCalib) {
            if(!track->TestFilterBit(BIT(8)) && !track->TestFilterBit(BIT(9)))
                return false;
            if(track->Pt()<0.2 || track->Pt()>5 || TMath::Abs(track->Eta())>0.8 || TMath::Abs(track->Eta())<0.)
                return false;
        }
    }
    else { //default
        if(fCalibType==kQnCalib) {
            if(!track->TestFilterBit(BIT(8)) && !track->TestFilterBit(BIT(9)))
                return false;
            if(track->Pt()<fPtMinTPC || track->Pt()>fPtMaxTPC || TMath::Abs(track->Eta())>fEtaMaxTPC || TMath::Abs(track->Eta())<fEtaMinTPC)
                return false;
        }
        else if(fCalibType==kQnFrameworkCalib) {
            if(!track->TestFilterBit(BIT(8)) && !track->TestFilterBit(BIT(9)))
                return false;
            if(track->Pt()<0.2 || track->Pt()>5 || TMath::Abs(track->Eta())>0.8 || TMath::Abs(track->Eta())<0.)
                return false;
        }
    }

    return true;
}


/// ==========================================================================================
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

    std::cout << "Success add AliAnalysisTaskRawJetWithEP!!"<< std::endl;

    return rawJetTask;
}


bool AliAnalysisTaskRawJetWithEP::ExtractRecentPara(TFile *EPCalibRefFile, TObjArray *lRefQx2am, TObjArray *lRefQy2am, TObjArray *lRefQx2as, TObjArray *lRefQy2as, TObjArray *lRefQx3am, TObjArray *lRefQy3am, TObjArray *lRefQx3as, TObjArray *lRefQy3as, TObjArray *lRefQx2cm, TObjArray *lRefQy2cm, TObjArray *lRefQx2cs, TObjArray *lRefQy2cs, TObjArray *lRefQx3cm,TObjArray *lRefQy3cm, TObjArray *lRefQx3cs, TObjArray *lRefQy3cs, TObjArray *lRefTPCposEta, TObjArray *lRefTPCnegEta) {
    for(int iZvtx = 0; iZvtx < 14; iZvtx++) {
        // V0 A-side #################################################################
        // Mean Qx correction
        // Includes check if Zvtx is differential
        AliOADBContainer* contQx2am = (AliOADBContainer*) EPCalibRefFile->Get(Form("fqxa2m_%d", iZvtx));
        if(!contQx2am) { //check if it is not Zvtx differential
            contQx2am = (AliOADBContainer*) EPCalibRefFile->Get("fqxa2m");
            if(contQx2am) fV0CalibZvtxDiff = false;
        }
        if(!contQx2am) {
            AliWarning("OADB object fqxa2m is not available in the file\n");
            return false;
        }
        
        lRefQx2am->Add(contQx2am);
        
        AliOADBContainer* contQy2am = nullptr; // Mean Qy correction
        AliOADBContainer* contQx2as = nullptr; // Sigma Qx correction
        AliOADBContainer* contQy2as = nullptr; // Sigma Qy correction
        if(fV0CalibZvtxDiff){
            contQy2am = (AliOADBContainer*) EPCalibRefFile->Get(Form("fqya2m_%d", iZvtx));
            contQx2as = (AliOADBContainer*) EPCalibRefFile->Get(Form("fqxa2s_%d", iZvtx));
            contQy2as = (AliOADBContainer*) EPCalibRefFile->Get(Form("fqya2s_%d", iZvtx));
        }
        else {
            contQy2am = (AliOADBContainer*) EPCalibRefFile->Get("fqya2m");
            contQx2as = (AliOADBContainer*) EPCalibRefFile->Get("fqxa2s");
            contQy2as = (AliOADBContainer*) EPCalibRefFile->Get("fqya2s");
        }
        if(!contQy2am) {
            AliWarning("OADB object fqya2m is not available in the file\n");
            return false;
        }
        lRefQy2am->Add(contQy2am);
        if(!contQx2as) {
            AliWarning("OADB object fqxa2s is not available in the file\n");
            return false;
        }
        lRefQx2as->Add(contQx2as);
        if(!contQy2as) {
            AliWarning("OADB object fqya2s is not available in the file\n");
            return false;
        }
        lRefQy2as->Add(contQy2as);


        AliOADBContainer* contQx3am = (AliOADBContainer*) EPCalibRefFile->Get(Form("fqxa3m_%d", iZvtx));
        if(!contQx3am) { //check if it is not Zvtx differential
            contQx3am = (AliOADBContainer*) EPCalibRefFile->Get("fqxa3m");
            if(contQx3am) fV0CalibZvtxDiff = false;
        }
        if(!contQx3am) {
            AliWarning("OADB object fqxa3m is not available in the file\n");
            return false;
        }
        lRefQx3am->Add(contQx3am);

        
        AliOADBContainer* contQy3am = nullptr; // Mean Qy correction
        AliOADBContainer* contQx3as = nullptr; // Sigma Qx correction
        AliOADBContainer* contQy3as = nullptr; // Sigma Qy correction
        if(fV0CalibZvtxDiff){
            contQy3am = (AliOADBContainer*) EPCalibRefFile->Get(Form("fqya3m_%d", iZvtx));
            contQx3as = (AliOADBContainer*) EPCalibRefFile->Get(Form("fqxa3s_%d", iZvtx));
            contQy3as = (AliOADBContainer*) EPCalibRefFile->Get(Form("fqya3s_%d", iZvtx));
        }
        else {
            contQy3am = (AliOADBContainer*) EPCalibRefFile->Get("fqya3m");
            contQx3as = (AliOADBContainer*) EPCalibRefFile->Get("fqxa3s");
            contQy3as = (AliOADBContainer*) EPCalibRefFile->Get("fqya3s");
        }
        if(!contQy3am) {
            AliWarning("OADB object fqya3m is not available in the file\n");
            return false;
        }
        lRefQy3am->Add(contQy3am);

        if(!contQx3as) {
            AliWarning("OADB object fqxa3s is not available in the file\n");
            return false;
        }
        lRefQx3as->Add(contQx3as);

        if(!contQy3as) {
            AliWarning("OADB object fqya3s is not available in the file\n");
            return false;
        }
        lRefQy3as->Add(contQy3as);


        // V0 C-side  #################################################################
        AliOADBContainer* contQx2cm = nullptr; // Mean Qx correction
        AliOADBContainer* contQy2cm = nullptr; // Mean Qy correction
        AliOADBContainer* contQx2cs = nullptr; // Sigma Qx correction
        AliOADBContainer* contQy2cs = nullptr; // Sigma Qy correction

        AliOADBContainer* contQx3cm = nullptr; // Mean Qx correction
        AliOADBContainer* contQy3cm = nullptr; // Mean Qy correction
        AliOADBContainer* contQx3cs = nullptr; // Sigma Qx correction
        AliOADBContainer* contQy3cs = nullptr; // Sigma Qy correction
        if(fV0CalibZvtxDiff){
            contQx2cm = (AliOADBContainer*) EPCalibRefFile->Get(Form("fqxc2m_%d", iZvtx));
            contQy2cm = (AliOADBContainer*) EPCalibRefFile->Get(Form("fqyc2m_%d", iZvtx));
            contQx2cs = (AliOADBContainer*) EPCalibRefFile->Get(Form("fqxc2s_%d", iZvtx));
            contQy2cs = (AliOADBContainer*) EPCalibRefFile->Get(Form("fqyc2s_%d", iZvtx));
            
            contQx3cm = (AliOADBContainer*) EPCalibRefFile->Get(Form("fqxc3m_%d", iZvtx));
            contQy3cm = (AliOADBContainer*) EPCalibRefFile->Get(Form("fqyc3m_%d", iZvtx));
            contQx3cs = (AliOADBContainer*) EPCalibRefFile->Get(Form("fqxc3s_%d", iZvtx));
            contQy3cs = (AliOADBContainer*) EPCalibRefFile->Get(Form("fqyc3s_%d", iZvtx));
        }
        else{
            contQx2cm = (AliOADBContainer*) EPCalibRefFile->Get("fqxc2m");
            contQy2cm = (AliOADBContainer*) EPCalibRefFile->Get("fqyc2m");
            contQx2cs = (AliOADBContainer*) EPCalibRefFile->Get("fqxc2s");
            contQy2cs = (AliOADBContainer*) EPCalibRefFile->Get("fqyc2s");

            contQx3cm = (AliOADBContainer*) EPCalibRefFile->Get("fqxc3m");
            contQy3cm = (AliOADBContainer*) EPCalibRefFile->Get("fqyc3m");
            contQx3cs = (AliOADBContainer*) EPCalibRefFile->Get("fqxc3s");
            contQy3cs = (AliOADBContainer*) EPCalibRefFile->Get("fqyc3s");
        }

        if(!contQx2cm) {
            AliWarning("OADB object fqxc2m is not available in the file\n");
            return false;
        }
        lRefQx2cm->Add(contQx2cm);

        if(!contQy2cm) {
            AliWarning("OADB object fqyc2m is not available in the file\n");
            return false;
        }
        lRefQy2cm->Add(contQy2cm);

        if(!contQx2cs) {
            AliWarning("OADB object fqxc2s is not available in the file\n");
            return false;
        }
        lRefQx2cs->Add(contQx2cs);

        if(!contQy2cs) {
            AliWarning("OADB object fqyc2s is not available in the file\n");
            return false;
        }
        lRefQy2cs->Add(contQy2cs);

        
        if(!contQx3cm) {
            AliWarning("OADB object fqxc3m is not available in the file\n");
            return false;
        }
        lRefQx3cm->Add(contQx3cm);

        if(!contQy3cm) {
            AliWarning("OADB object fqyc3m is not available in the file\n");
            return false;
        }
        lRefQy3cm->Add(contQy3cm);

        if(!contQx3cs) {
            AliWarning("OADB object fqxc2s is not available in the file\n");
            return false;
        }
        lRefQx3cs->Add(contQx3cs);

        if(!contQy3cs) {
            AliWarning("OADB object fqyc2s is not available in the file\n");
            return false;
        }
        lRefQy3cs->Add(contQy3cs);
    }
    
    
    //load TPC calibrations (not mandatory)
    for(int iCent = 0; iCent < 9; iCent++) {
        AliOADBContainer* contTPCposEta = (AliOADBContainer*) EPCalibRefFile->Get(Form("fphidistr_poseta_%d_%d", iCent*10, (iCent+1)*10));
        if(!contTPCposEta) {
            AliWarning("OADB object fphidistr_poseta (TPC Calibration) is not available in the file\n");
        }
        lRefTPCposEta->Add(contTPCposEta);

        AliOADBContainer* contTPCnegEta = (AliOADBContainer*) EPCalibRefFile->Get(Form("fphidistr_negeta_%d_%d", iCent*10, (iCent+1)*10));
        if(!contTPCnegEta) {
            AliWarning("OADB object fphidistr_negeta (TPC Calibration) is not available in the file\n");
            fWeightsTPCNegEta[iCent] = nullptr;
        }
        lRefTPCnegEta->Add(contTPCnegEta);
    }
    
    return true;
}

bool AliAnalysisTaskRawJetWithEP::OpenInfoCalbration() 
{
    if(!fMultV0BefCorPfpx) {
        AliWarning("OADB object hMultV0BefCorPfpx is not available\n");
        return false;
    }
    
    if(!(fMultV0BefCorPfpx->GetObject(fRun))) {
        AliWarning(Form("OADB object hMultV0BefCorPfpx is not available for run %i\n", fRun));
        return false;
    }

    fHistMultV0 = ((TH1D*) fMultV0BefCorPfpx->GetObject(fRun));
    
    for(int iZvtx = 0; iZvtx < 14; iZvtx++) {
        AliOADBContainer* contQx2am = 0;
        // If we do not have z-vertex differential objects, then only the first index is 
        // in the OADBContainer array
        if(fV0CalibZvtxDiff) contQx2am = (AliOADBContainer* ) fOADBzArray_contQx2am->At(iZvtx);
        else contQx2am = (AliOADBContainer* ) fOADBzArray_contQx2am->At(0);
        if(!contQx2am) {
            AliWarning("OADB object fqxa2m is not available\n");
            return false;
        }
        if(!(contQx2am->GetObject(fRun))) {
            AliWarning(Form("OADB object fqxa2m is not available for run %i\n", fRun));
            return false;
        }
        fQx2mV0A[iZvtx] = ((TH1D*) contQx2am->GetObject(fRun));
        
        AliOADBContainer* contQy2am = 0;
        if(fV0CalibZvtxDiff) contQy2am = (AliOADBContainer* ) fOADBzArray_contQy2am->At(iZvtx);
        else contQy2am = (AliOADBContainer* ) fOADBzArray_contQy2am->At(0);
        if(!contQy2am) {
            AliWarning("OADB object fqya2m is not available\n");
            return false;
        }
        if(!(contQy2am->GetObject(fRun))) {
            AliWarning(Form("OADB object fqya2m is not available for run %i\n", fRun));
            return false;
        }
        fQy2mV0A[iZvtx] = ((TH1D*) contQy2am->GetObject(fRun));
        
        AliOADBContainer* contQx2as = 0;
        if(fV0CalibZvtxDiff) contQx2as = (AliOADBContainer* ) fOADBzArray_contQx2as->At(iZvtx);
        else contQx2as = (AliOADBContainer* ) fOADBzArray_contQx2as->At(0);
        if(!contQx2as) {
            AliWarning("OADB object fqxa2s is not available\n");
            return false;
        }
        if(!(contQx2as->GetObject(fRun))) {
            AliWarning(Form("OADB object fqxa2s is not available for run %i\n", fRun));
            return false;
        }
        fQx2sV0A[iZvtx] = ((TH1D*) contQx2as->GetObject(fRun));
        
        AliOADBContainer* contQy2as = 0;
        if(fV0CalibZvtxDiff) contQy2as = (AliOADBContainer* ) fOADBzArray_contQy2as->At(iZvtx);
        else contQy2as = (AliOADBContainer* ) fOADBzArray_contQy2as->At(0);
        if(!contQy2as) {
            AliWarning("OADB object fqya2s is not available\n");
            return false;
        }
        if(!(contQy2as->GetObject(fRun))) {
            AliWarning(Form("OADB object fqya2s is not available for run %i\n", fRun));
            return false;
        }
        fQy2sV0A[iZvtx] = ((TH1D*) contQy2as->GetObject(fRun));
        
        AliOADBContainer* contQx2cm = 0;
        if(fV0CalibZvtxDiff) contQx2cm = (AliOADBContainer* ) fOADBzArray_contQx2cm->At(iZvtx);
        else contQx2cm = (AliOADBContainer* ) fOADBzArray_contQx2cm->At(0);
        if(!contQx2cm) {
            AliWarning("OADB object fqxc2m is not available\n");
            return false;
        }
        if(!(contQx2cm->GetObject(fRun))) {
            AliWarning(Form("OADB object fqxc2m is not available for run %i\n", fRun));
            return false;
        }
        fQx2mV0C[iZvtx] = ((TH1D*) contQx2cm->GetObject(fRun));
        
        AliOADBContainer* contQy2cm = 0;
        if(fV0CalibZvtxDiff) contQy2cm = (AliOADBContainer* ) fOADBzArray_contQy2cm->At(iZvtx);
        else contQy2cm = (AliOADBContainer* ) fOADBzArray_contQy2cm->At(0);
        if(!contQy2cm) {
            AliWarning("OADB object fqyc2m is not available\n");
            return false;
        }
        if(!(contQy2cm->GetObject(fRun))) {
            AliWarning(Form("OADB object fqyc2m is not available for run %i\n", fRun));
            return false;
        }
        fQy2mV0C[iZvtx] = ((TH1D*) contQy2cm->GetObject(fRun));

        AliOADBContainer* contQx2cs = 0;
        if(fV0CalibZvtxDiff) contQx2cs = (AliOADBContainer* ) fOADBzArray_contQx2cs->At(iZvtx);
        else contQx2cs = (AliOADBContainer* ) fOADBzArray_contQx2cs->At(0);
        if(!contQx2cs) {
            AliWarning("OADB object fqxc%ds is not available\n");
            return false;
        }
        if(!(contQx2cs->GetObject(fRun))) {
            AliWarning(Form("OADB object fqxc2s is not available for run %i\n",  fRun));
            return false;
        }
        fQx2sV0C[iZvtx] = ((TH1D*) contQx2cs->GetObject(fRun));
        
        AliOADBContainer* contQy2cs = nullptr;
        if(fV0CalibZvtxDiff) contQy2cs = (AliOADBContainer* ) fOADBzArray_contQy2cs->At(iZvtx);
        else contQy2cs = (AliOADBContainer* ) fOADBzArray_contQy2cs->At(0);
        if(!contQy2cs) {
            AliWarning("OADB object fqyc%ds is not available\n");
            return false;
        }
        if(!(contQy2cs->GetObject(fRun))) {
            AliWarning(Form("OADB object fqyc2s is not available for run %i\n",  fRun));
            return false;
        }
        fQy2sV0C[iZvtx] = ((TH1D*) contQy2cs->GetObject(fRun));

        // == Q3 Vector
        AliOADBContainer* contQx3am = 0;
        // If we do not have z-vertex differential objects, then only the first index is 
        // in the OADBContainer array
        if(fV0CalibZvtxDiff) contQx3am = (AliOADBContainer* ) fOADBzArray_contQx3am->At(iZvtx);
        else contQx3am = (AliOADBContainer* ) fOADBzArray_contQx3am->At(0);
        if(!contQx3am) {
            AliWarning("OADB object fqxa3m is not available\n");
            return false;
        }
        if(!(contQx3am->GetObject(fRun))) {
            AliWarning(Form("OADB object fqxa3m is not available for run %i\n", fRun));
            return false;
        }
        fQx3mV0A[iZvtx] = ((TH1D*) contQx3am->GetObject(fRun));
        
        AliOADBContainer* contQy3am = 0;
        if(fV0CalibZvtxDiff) contQy3am = (AliOADBContainer* ) fOADBzArray_contQy3am->At(iZvtx);
        else contQy3am = (AliOADBContainer* ) fOADBzArray_contQy3am->At(0);
        if(!contQy3am) {
            AliWarning("OADB object fqya3m is not available\n");
            return false;
        }
        if(!(contQy3am->GetObject(fRun))) {
            AliWarning(Form("OADB object fqya3m is not available for run %i\n", fRun));
            return false;
        }
        fQy3mV0A[iZvtx] = ((TH1D*) contQy3am->GetObject(fRun));
        
        AliOADBContainer* contQx3as = 0;
        if(fV0CalibZvtxDiff) contQx3as = (AliOADBContainer* ) fOADBzArray_contQx3as->At(iZvtx);
        else contQx3as = (AliOADBContainer* ) fOADBzArray_contQx3as->At(0);
        if(!contQx3as) {
            AliWarning("OADB object fqxa2s is not available\n");
            return false;
        }
        if(!(contQx3as->GetObject(fRun))) {
            AliWarning(Form("OADB object fqxa2s is not available for run %i\n", fRun));
            return false;
        }
        fQx3sV0A[iZvtx] = ((TH1D*) contQx3as->GetObject(fRun));
        
        AliOADBContainer* contQy3as = 0;
        if(fV0CalibZvtxDiff) contQy3as = (AliOADBContainer* ) fOADBzArray_contQy3as->At(iZvtx);
        else contQy3as = (AliOADBContainer* ) fOADBzArray_contQy3as->At(0);
        if(!contQy3as) {
            AliWarning("OADB object fqya2s is not available\n");
            return false;
        }
        if(!(contQy3as->GetObject(fRun))) {
            AliWarning(Form("OADB object fqya2s is not available for run %i\n", fRun));
            return false;
        }
        fQy3sV0A[iZvtx] = ((TH1D*) contQy3as->GetObject(fRun));
        
        AliOADBContainer* contQx3cm = 0;
        if(fV0CalibZvtxDiff) contQx3cm = (AliOADBContainer* ) fOADBzArray_contQx3cm->At(iZvtx);
        else contQx3cm = (AliOADBContainer* ) fOADBzArray_contQx3cm->At(0);
        if(!contQx3cm) {
            AliWarning("OADB object fqxc2m is not available\n");
            return false;
        }
        if(!(contQx3cm->GetObject(fRun))) {
            AliWarning(Form("OADB object fqxc2m is not available for run %i\n", fRun));
            return false;
        }
        fQx3mV0C[iZvtx] = ((TH1D*) contQx3cm->GetObject(fRun));
        
        AliOADBContainer* contQy3cm = 0;
        if(fV0CalibZvtxDiff) contQy3cm = (AliOADBContainer* ) fOADBzArray_contQy3cm->At(iZvtx);
        else contQy3cm = (AliOADBContainer* ) fOADBzArray_contQy3cm->At(0);
        if(!contQy3cm) {
            AliWarning("OADB object fqyc2m is not available\n");
            return false;
        }
        if(!(contQy3cm->GetObject(fRun))) {
            AliWarning(Form("OADB object fqyc2m is not available for run %i\n", fRun));
            return false;
        }
        fQy3mV0C[iZvtx] = ((TH1D*) contQy3cm->GetObject(fRun));

        AliOADBContainer* contQx3cs = 0;
        if(fV0CalibZvtxDiff) contQx3cs = (AliOADBContainer* ) fOADBzArray_contQx3cs->At(iZvtx);
        else contQx3cs = (AliOADBContainer* ) fOADBzArray_contQx3cs->At(0);
        if(!contQx3cs) {
            AliWarning("OADB object fqxc%ds is not available\n");
            return false;
        }
        if(!(contQx3cs->GetObject(fRun))) {
            AliWarning(Form("OADB object fqxc2s is not available for run %i\n",  fRun));
            return false;
        }
        fQx3sV0C[iZvtx] = ((TH1D*) contQx3cs->GetObject(fRun));
        
        AliOADBContainer* contQy3cs = nullptr;
        if(fV0CalibZvtxDiff) contQy3cs = (AliOADBContainer* ) fOADBzArray_contQy3cs->At(iZvtx);
        else contQy3cs = (AliOADBContainer* ) fOADBzArray_contQy3cs->At(0);
        if(!contQy3cs) {
            AliWarning("OADB object fqyc%ds is not available\n");
            return false;
        }
        if(!(contQy3cs->GetObject(fRun))) {
            AliWarning(Form("OADB object fqyc2s is not available for run %i\n",  fRun));
            return false;
        }
        fQy3sV0C[iZvtx] = ((TH1D*) contQy3cs->GetObject(fRun));

        if(!fV0CalibZvtxDiff) //assign only first element of array if it is not Zvtx differential
            break;
    }
    
    //load TPC calibrations (not mandatory)
    for(int iCent = 0; iCent < 9; iCent++) {
        AliOADBContainer* contTPCposEta = 0;
        contTPCposEta = (AliOADBContainer* ) fOADBcentArray_contTPCposEta->At(iCent);
        if(!contTPCposEta) {
            AliWarning("OADB object fphidistr_poseta (TPC Calibration) is not available\n");
            fWeightsTPCPosEta[iCent] = nullptr;
        }
        else {
            if(!(contTPCposEta->GetObject(fRun))) {
                AliWarning(Form("OADB object fphidistr_poseta (TPC Calibration) is not available for run %i\n", fRun));
                fWeightsTPCPosEta[iCent] = nullptr;
            }
            else {
                fWeightsTPCPosEta[iCent] = ((TH1D*) contTPCposEta->GetObject(fRun));   
            }
        }

        AliOADBContainer* contTPCnegEta = 0;
        contTPCnegEta = (AliOADBContainer* ) fOADBcentArray_contTPCnegEta->At(iCent);
        if(!contTPCnegEta) {
            AliWarning("OADB object fphidistr_negeta (TPC Calibration) is not available in the file\n");
            fWeightsTPCNegEta[iCent] = nullptr;
            return true;
        }
        else {        
            if(!(contTPCnegEta->GetObject(fRun))) {
                AliWarning(Form("OADB object fphidistr_negeta (TPC Calibration) is not available for run %i\n", fRun));
                fWeightsTPCNegEta[iCent] = nullptr;
            }
            else {
                fWeightsTPCNegEta[iCent] = ((TH1D*) contTPCnegEta->GetObject(fRun));   
            }
        }
    }

    return true;
}




