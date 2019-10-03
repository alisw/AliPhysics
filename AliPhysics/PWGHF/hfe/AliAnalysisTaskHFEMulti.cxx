/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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
//
// The analysis task:
// Fork of the AliAnalysisTaskHFE to do studies in multiplicity
// 
// Author:
// Jan Wagner - j.wagner@cern.ch
//
//
#include <TAxis.h>
#include <TBits.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH3D.h>
#include <TIterator.h>
#include <TList.h>
#include <TLegend.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TParticle.h>
#include <TProfile.h>
#include <TString.h>
#include <TF1.h>
#include <TTree.h>
#include <TRandom3.h>

#include "AliESDtrackCuts.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisUtils.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliAODTrack.h"
#include "AliAODVertex.h"
#include "AliCentrality.h"
#include "AliCFContainer.h"
#include "AliCFManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMCParticle.h"
#include "AliMultiplicity.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliOADBContainer.h"
#include "AliStack.h"
#include "AliTriggerAnalysis.h"
#include "AliVVertex.h"

#include "AliHFEcollection.h"
#include "AliHFEcontainer.h"
#include "AliHFEcuts.h"
#include "AliHFEelecbackground.h"
#include "AliHFENonPhotonicElectron.h"
#include "AliHFEmcQA.h"
#include "AliHFEpairs.h"
#include "AliHFEpid.h"
#include "AliHFEpidQAmanager.h"
#include "AliHFEsecVtxs.h"
#include "AliHFEsecVtx.h"
#include "AliHFEsignalCuts.h"
#include "AliHFEtaggedTrackAnalysis.h"
#include "AliHFEtools.h"
#include "AliHFEparamBag.h"
#include "AliHFEV0taginfo.h"
#include "AliHFEvarManager.h"
#include "AliAnalysisTaskHFEMulti.h"
#include "AliAODMCHeader.h"
#include "TClonesArray.h"

ClassImp(AliAnalysisTaskHFEMulti)

//____________________________________________________________
AliAnalysisTaskHFEMulti::AliAnalysisTaskHFEMulti():
AliAnalysisTaskSE("PID efficiency Analysis")
    , fAODMCHeader(NULL)
    , fAODArrayMCInfo(NULL)
    , fQAlevel(0)
    , fPlugins(0)
    , fCollisionSystem(3)
    , fFillSignalOnly(kTRUE)
    , fRejectMCFakeTracks(kFALSE)
    , fFillNoCuts(kFALSE)
    , fBackGroundFactorApply(kFALSE)
    , fRemovePileUp(kFALSE)
    , fIdentifiedAsPileUp(kFALSE)
    , fIdentifiedAsOutInz(kFALSE)
    , fPassTheEventCut(kFALSE)
    , fRejectKinkMother(kFALSE)
    , fisppMultiBin(kFALSE)
    , fPbPbUserCentralityBinning(kFALSE)
    , fRemoveFirstEvent(kFALSE)
    , fisNonHFEsystematics(kFALSE)
    , fSpecialTrigger(NULL)
    , fCentralityF(-1)
    , fCentralityPercent(-1)
    , fCentralityEstimator("V0M")
    , fContributors(0.5)
    , fSPDtracklets(0)
    , fSPDtrackletsCorr(0.0)
    , fSPDtrkF(-1)
    , fWeightBackGround(0.)
    , fVz(0.0)
    , fContainer(NULL)
    , fVarManager(NULL)
    , fSignalCuts(NULL)
    , fCFM(NULL)
    , fPID(NULL)
    , fPIDqa(NULL)
    , fPIDpreselect(NULL)
    , fCuts(NULL)
    , fAnalysisUtils(NULL)
    , fCutspreselect(NULL)
    , fMCQA(NULL)
    , fExtraCuts(NULL)
    , fBackgroundSubtraction(NULL)
    , fMCNtrWeight(NULL)
    , fRefMulti(-1.)
    , fMultiEstimatorSystem(kNtrk10)
    , fPIDResponse(NULL)
    , fQA(NULL)
    , fOutput(NULL)
    , fParams(NULL)
    , fHistMCQA(NULL)
    , fHistSECVTX(NULL)
    , fHistELECBACKGROUND(NULL)
    , fQACollection(NULL)
{
    //
    // Dummy constructor
    //
    memset(fElecBackgroundFactor, 0, sizeof(Double_t) * kElecBgSpecies * kBgPtBins * kCentBins * kBgLevels);
    memset(fkBackGroundFactorArray, 0, sizeof(TF1 *) * 12);
    memset(fBinLimit, 0, sizeof(Double_t) * (kBgPtBins+1));
    memset(&fisppMultiBin, kFALSE, sizeof(fisppMultiBin));
    memset(fCentralityLimits, 0, sizeof(Float_t) * 12);
    memset(fMultEstimatorAvg,0,sizeof(TProfile *) *4);

    SetppAnalysis();
}

//____________________________________________________________
AliAnalysisTaskHFEMulti::AliAnalysisTaskHFEMulti(const char * name):
    AliAnalysisTaskSE(name)
    , fAODMCHeader(NULL)
    , fAODArrayMCInfo(NULL)
    , fQAlevel(0)
    , fPlugins(0)
    , fCollisionSystem(3)
    , fFillSignalOnly(kTRUE)
    , fRejectMCFakeTracks(kFALSE)
    , fFillNoCuts(kFALSE)
    , fBackGroundFactorApply(kFALSE)
    , fRemovePileUp(kFALSE)
    , fIdentifiedAsPileUp(kFALSE)
    , fIdentifiedAsOutInz(kFALSE)
    , fPassTheEventCut(kFALSE)  
    , fRejectKinkMother(kFALSE)
    , fisppMultiBin(kFALSE)
    , fPbPbUserCentralityBinning(kFALSE)
    , fRemoveFirstEvent(kFALSE)
    , fisNonHFEsystematics(kFALSE)
    , fSpecialTrigger(NULL)
    , fCentralityF(-1)
    , fCentralityPercent(-1)
    , fCentralityEstimator("V0M")
    , fContributors(0.5)
    , fSPDtracklets(0)
    , fSPDtrackletsCorr(0.0)
    , fSPDtrkF(-1)
    , fWeightBackGround(0.)
    , fVz(0.0)
    , fContainer(NULL)
    , fVarManager(NULL)
    , fSignalCuts(NULL)
    , fCFM(NULL)
    , fPID(NULL)
    , fPIDqa(NULL)
    , fPIDpreselect(NULL)
    , fCuts(NULL)
    , fAnalysisUtils(NULL)
    , fCutspreselect(NULL)
    , fMCQA(NULL)
    , fExtraCuts(NULL)
    , fBackgroundSubtraction(NULL)
    , fMCNtrWeight(NULL)
    , fRefMulti(-1.)
    , fMultiEstimatorSystem(kNtrk10)
    , fPIDResponse(NULL)
    , fQA(NULL)
    , fOutput(NULL)
    , fParams(NULL)
    , fHistMCQA(NULL)
    , fHistSECVTX(NULL)
    , fHistELECBACKGROUND(NULL)
    , fQACollection(0x0)
{
    //
    // Default constructor
    // 
    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class());
    DefineOutput(3, TObject::Class());

    fPID = new AliHFEpid("hfePid");
    fPIDqa = new AliHFEpidQAmanager;
    fVarManager = new AliHFEvarManager("hfeVarManager");
    fAnalysisUtils = new AliAnalysisUtils;

    memset(fElecBackgroundFactor, 0, sizeof(Double_t) * kElecBgSpecies * kBgPtBins * kCentBins * kBgLevels);
    memset(fkBackGroundFactorArray, 0, sizeof(TF1 *) * 12);
    memset(fBinLimit, 0, sizeof(Double_t) * (kBgPtBins+1));
    memset(&fisppMultiBin, kFALSE, sizeof(fisppMultiBin));
    memset(fCentralityLimits, 0, sizeof(Float_t) * 12);
    memset(fMultEstimatorAvg,0,sizeof(TProfile *) *4);

    SetppAnalysis();
}

//____________________________________________________________
AliAnalysisTaskHFEMulti::AliAnalysisTaskHFEMulti(const AliAnalysisTaskHFEMulti &ref):
    AliAnalysisTaskSE(ref)
    , fAODMCHeader(NULL)
    , fAODArrayMCInfo(NULL)
    , fQAlevel(0)
    , fPlugins(0)
    , fCollisionSystem(ref.fCollisionSystem)
    , fFillSignalOnly(ref.fFillSignalOnly)
    , fRejectMCFakeTracks(ref.fRejectMCFakeTracks)
    , fFillNoCuts(ref.fFillNoCuts)
    , fBackGroundFactorApply(ref.fBackGroundFactorApply)
    , fRemovePileUp(ref.fRemovePileUp)
    , fIdentifiedAsPileUp(ref.fIdentifiedAsPileUp)
    , fIdentifiedAsOutInz(ref.fIdentifiedAsOutInz)
    , fPassTheEventCut(ref.fPassTheEventCut)
    , fRejectKinkMother(ref.fRejectKinkMother)
    , fisppMultiBin(ref.fisppMultiBin)
    , fPbPbUserCentralityBinning(ref.fPbPbUserCentralityBinning)
    , fRemoveFirstEvent(ref.fRemoveFirstEvent)
    , fisNonHFEsystematics(ref.fisNonHFEsystematics)
    , fSpecialTrigger(ref.fSpecialTrigger)
    , fCentralityF(ref.fCentralityF)
    , fCentralityPercent(ref.fCentralityPercent)
    , fCentralityEstimator(ref.fCentralityEstimator)
    , fContributors(ref.fContributors)
    , fSPDtracklets(ref.fSPDtracklets)
    , fSPDtrackletsCorr(ref.fSPDtrackletsCorr)
    , fSPDtrkF(ref.fSPDtrkF)
    , fWeightBackGround(ref.fWeightBackGround)
    , fVz(ref.fVz)
    , fContainer(NULL)
    , fVarManager(NULL)
    , fSignalCuts(NULL)
    , fCFM(NULL)
    , fPID(NULL)
    , fPIDqa(NULL)
    , fPIDpreselect(NULL)
    , fCuts(NULL)
    , fAnalysisUtils(NULL)
    , fCutspreselect(NULL)
    , fMCQA(NULL)
    , fExtraCuts(NULL)
    , fBackgroundSubtraction(NULL)
    , fMCNtrWeight(ref.fMCNtrWeight)
    , fRefMulti(-1.)
    , fMultiEstimatorSystem(kNtrk10)
    , fPIDResponse(ref.fPIDResponse)
    , fQA(NULL)
    , fOutput(NULL)
    , fParams(NULL)
    , fHistMCQA(NULL)
    , fHistSECVTX(NULL)
    , fHistELECBACKGROUND(NULL)
    , fQACollection(NULL)
{
    //
    // Copy Constructor
    //
    ref.Copy(*this);
}

//____________________________________________________________
AliAnalysisTaskHFEMulti &AliAnalysisTaskHFEMulti::operator=(const AliAnalysisTaskHFEMulti &ref){
    //
    // Assignment operator
    //
    if(this == &ref) 
        ref.Copy(*this);
    return *this;
}

//____________________________________________________________
void AliAnalysisTaskHFEMulti::Copy(TObject &o) const {
    // 
    // Copy into object o
    //
    AliAnalysisTaskHFEMulti &target = dynamic_cast<AliAnalysisTaskHFEMulti &>(o);
    target.fAODMCHeader = fAODMCHeader;
    target.fAODArrayMCInfo = fAODArrayMCInfo;
    target.fQAlevel = fQAlevel;
    target.fPlugins = fPlugins;
    target.fCollisionSystem = fCollisionSystem;
    target.fFillSignalOnly = fFillSignalOnly;
    target.fRejectMCFakeTracks = fRejectMCFakeTracks;
    target.fFillNoCuts = fFillNoCuts;
    target.fBackGroundFactorApply = fBackGroundFactorApply;
    target.fRemovePileUp = fRemovePileUp;
    target.fIdentifiedAsPileUp = fIdentifiedAsPileUp;
    target.fIdentifiedAsOutInz = fIdentifiedAsOutInz;
    target.fPassTheEventCut = fPassTheEventCut;
    target.fRejectKinkMother = fRejectKinkMother;
    target.fisppMultiBin =   fisppMultiBin;
    target.fPbPbUserCentralityBinning = fPbPbUserCentralityBinning;
    target.fRemoveFirstEvent = fRemoveFirstEvent;
    target.fisNonHFEsystematics = fisNonHFEsystematics;
    target.fSpecialTrigger = fSpecialTrigger;
    target.fCentralityF = fCentralityF;
    target.fCentralityPercent = fCentralityPercent;
    target.fCentralityEstimator = fCentralityEstimator;
    target.fContributors = fContributors;
    target.fSPDtracklets = fSPDtracklets;
    target.fSPDtrackletsCorr = fSPDtrackletsCorr;
    target.fSPDtrkF = fSPDtrkF;
    target.fWeightBackGround = fWeightBackGround;
    target.fVz = fVz;
    target.fContainer = fContainer;
    target.fVarManager = fVarManager;
    target.fSignalCuts = fSignalCuts;
    target.fCFM = fCFM;
    target.fPID = fPID;
    target.fPIDqa = fPIDqa;
    target.fPIDpreselect = fPIDpreselect;
    target.fCuts = fCuts;
    target.fAnalysisUtils = fAnalysisUtils;
    target.fCutspreselect = fCutspreselect;
    target.fMCQA = fMCQA;
    target.fExtraCuts = fExtraCuts;
    target.fBackgroundSubtraction = fBackgroundSubtraction;
    target.fMCNtrWeight = fMCNtrWeight;
    target.fRefMulti = fRefMulti;
    target.fMultiEstimatorSystem = fMultiEstimatorSystem;
    target.fPIDResponse = fPIDResponse;
    target.fQA = fQA;
    target.fOutput = fOutput;
    target.fParams = fParams;
    target.fHistMCQA = fHistMCQA;
    target.fHistSECVTX = fHistSECVTX;
    target.fHistELECBACKGROUND = fHistELECBACKGROUND;
    target.fQACollection = fQACollection;
}

//____________________________________________________________
AliAnalysisTaskHFEMulti::~AliAnalysisTaskHFEMulti(){
    //
    // Destructor
    //
    if(fPID) delete fPID;
    if(fPIDpreselect) delete fPIDpreselect;
    if(fVarManager) delete fVarManager;
    if(fCFM) delete fCFM;
    if(fSignalCuts) delete fSignalCuts;
    if(fMCQA) delete fMCQA;
    if(fBackgroundSubtraction) delete fBackgroundSubtraction;
    if(fSpecialTrigger) delete fSpecialTrigger;
    if(fAnalysisUtils) delete fAnalysisUtils;
    // Delete output objects only if we are not running in PROOF mode because otherwise this produces a crash during merging
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if(mgr && mgr->GetAnalysisType() != AliAnalysisManager::kProofAnalysis){
        if(fPIDqa) delete fPIDqa;
        if(fOutput) delete fOutput;
        if(fParams) delete fParams;
        if(fQA) delete fQA;
    }
}

//____________________________________________________________
void AliAnalysisTaskHFEMulti::UserCreateOutputObjects(){
    //
    // Creating output container and output objects
    // Here we also Initialize the correction framework container and 
    // the objects for
    // - PID
    // - MC QA
    // - SecVtx
    // QA histograms are created if requested
    // Called once per worker
    //
    AliDebug(3, "Creating Output Objects");

    // Make lists for Output
    if(!fQA) fQA = new TList;
    fQA->SetOwner();
    if(!fOutput) fOutput = new TList;
    fOutput->SetOwner();
    if(!fParams) {
        AliError("Parameter class not set, using default");
        fParams = new AliHFEparamBag("default");
    }

    // Automatic determination of the analysis mode
    AliVEventHandler *inputHandler = dynamic_cast<AliVEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    if(!TString(inputHandler->IsA()->GetName()).CompareTo("AliAODInputHandler")){
        SetAODAnalysis();
    } else {
        SetESDAnalysis();
        if(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())
            SetHasMCData();
    }
    printf("Analysis Mode: %s Analysis\n", IsAODanalysis() ? "AOD" : "ESD");
    printf("MC Data available %s\n", HasMCData() ? "Yes" : "No");

    // Initialize PIDResponse handler
    fPIDResponse = dynamic_cast<AliInputEventHandler *>(inputHandler)->GetPIDResponse();

    // First Part: Make QA histograms
    fQACollection = new AliHFEcollection("TaskQA", "QA histos from the Electron Task");
    //fQACollection->CreateTH1F("nElectronTracksEvent", "Number of Electron Candidates", 100, 0, 100);
    //fQACollection->CreateTH1F("nElectron", "Number of electrons", 100, 0, 100);
    fQACollection->CreateTH2F("NtrVsZ","Ntracklets vs Z vertex",300,-15,15,375,-0.5,374.5);
    fQACollection->CreateTH2F("NtrCorrVsZ","Ntracklets after correction vs Z vertex",300,-15,15,375,-0.5,374.5);
    fQACollection->CreateTH2F("NchCorrVsNtr","Ntracklets after correction vs Ncharged;N_{tr}^{cor};N_{ch}",375,-0.5,374.5,375,-0.5,374.5);
    fQACollection->CreateTH2F("NtrCorrVsNch","Ntracklets after correction vs Ncharged;N_{ch};N_{tr}^{cor}",375,-0.5,374.5,375,-0.5,374.5);
    fQACollection->CreateTH2F("VZEROVsZ","VZERO vs Z vertex",300,-15,15,400,0,800);
    fQACollection->CreateTH2F("VZEROCorrVsZ","VZERO after correction vs Z vertex",300,-15,15,400,0,800);
    fQA->Add(fQACollection);

    for(int i=0;i<4;i++){
        if(fMultEstimatorAvg[i]){
            TProfile* hprof=new TProfile(*fMultEstimatorAvg[i]);
            hprof->SetName(Form("Profile%d",i));
            fQA->Add(hprof);
        }
    }

    // Initialize PID
    fPID->SetHasMCData(HasMCData());
    if(!fPID->GetNumberOfPIDdetectors()) fPID->AddDetector("TPC", 0);
    if(IsQAOn(kPIDqa)){
        AliInfo("PID QA switched on");
        fPIDqa->Initialize(fPID);
        fQA->Add(fPIDqa->MakeList("HFEpidQA"));
    }
    fPID->SortDetectors();

    // Background subtraction-------------------------------------------------------------------
    if (GetPlugin(kNonPhotonicElectron)){
        if(!fBackgroundSubtraction) fBackgroundSubtraction = new AliHFENonPhotonicElectron();
        if(IsAODanalysis()) fBackgroundSubtraction->SetAOD(kTRUE);
        fBackgroundSubtraction->Init();
        fOutput->Add(fBackgroundSubtraction->GetListOutput());
    }
    //------------------------------------------------------------------------------------------


    // Initialize correction Framework and Cuts
    const Int_t kNcutSteps = AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kNcutStepsRecTrack + AliHFEcuts::kNcutStepsDETrack + AliHFEcuts::kNcutStepsSecvtxTrack;
    fCFM = new AliCFManager;
    fCFM->SetNStepParticle(kNcutSteps);
    MakeParticleContainer();
    MakeEventContainer();
    // Temporary fix: Initialize particle cuts with NULL
    for(Int_t istep = 0; istep < kNcutSteps; istep++)
        fCFM->SetParticleCutsList(istep, NULL);
    if(!fCuts){
        AliWarning("Cuts not available. Default cuts will be used");
        fCuts = new AliHFEcuts;
        fCuts->CreateStandardCuts();
    }
    if(IsAODanalysis()) fCuts->SetAOD();
    // Make clone for V0 tagging step
    fCuts->SetPIDResponse(fPIDResponse);
    fCuts->Initialize(fCFM);
    if(fCuts->IsQAOn()) fQA->Add(fCuts->GetQAhistograms());
    fSignalCuts = new AliHFEsignalCuts("HFEsignalCuts", "HFE MC Signal definition");
    fVarManager->SetSignalCuts(fSignalCuts);

    // add output objects to the List
    fOutput->AddAt(fContainer, 0);
    fOutput->AddAt(fCFM->GetEventContainer(), 1);

    // mcQA----------------------------------
    if (HasMCData() && IsQAOn(kMCqa)) {
        AliInfo("MC QA on");
        if(!fMCQA) fMCQA = new AliHFEmcQA;
        if(!fHistMCQA) fHistMCQA = new TList();
        fHistMCQA->SetOwner();
        if(IsPbPb()) fMCQA->SetPbPb();
        if(fisppMultiBin) fMCQA->SetPPMultiBin();
        if(TestBit(kTreeStream)){
            fMCQA->EnableDebugStreamer();
        }
        if(TestBit(kWeightHist)){
            fMCQA->EnableGetWeightHist();
        }
        fMCQA->CreatDefaultHistograms(fHistMCQA);
        fMCQA->SetBackgroundWeightFactor(fElecBackgroundFactor[0][0][0],fBinLimit);
        fQA->Add(fHistMCQA);
    } 


    //fQA->Print();

    PrintStatus();
    // Done!!!
    PostData(1, fOutput);
    PostData(2, fQA);
    PostData(3, fParams);
}

//____________________________________________________________
void AliAnalysisTaskHFEMulti::UserExec(Option_t *){
    //
    // Run the analysis
    // 

    //printf("test00\n");

    AliDebug(3, "Starting Single Event Analysis");
    if(!fInputEvent){
        AliError("Reconstructed Event not available");
        return;
    }
    if(HasMCData() && IsESDanalysis()){
        AliDebug(4, Form("MC Event: %p", fMCEvent));
        if(!fMCEvent){
            AliError("No MC Event, but MC Data required");
            return;
        }
    }
    if(!fCuts){
        AliError("HFE cuts not available");
        return;
    }
    if(!fPID->IsInitialized()){
        // Initialize PID with the given run number
        fPID->InitializePID(fInputEvent->GetRunNumber());
    }

    if(fRemoveFirstEvent){
        if(fAnalysisUtils->IsFirstEventInChunk(fInputEvent)) return;
    }

    //AliESDEvent *ev = dynamic_cast<AliESDEvent *>(fInputEvent);

    if(IsESDanalysis() && HasMCData()){
        // Protect against missing MC trees
        AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
        if(!mcH){ 
            AliError("No MC Event Handler available");
            return;
        }
        if(!mcH->InitOk()) return;
        if(!mcH->TreeK()) return;
        if(!mcH->TreeTR()) return;

        // Background subtraction-------------------------------------------------------------------
        if(GetPlugin(kNonPhotonicElectron)) fBackgroundSubtraction->SetMCEvent(fMCEvent);
        //------------------------------------------------------------------------------------------
    }

    if(IsAODanalysis() && HasMCData()){
        // take MC info
        AliAODEvent *aodE = dynamic_cast<AliAODEvent *>(fInputEvent);
        if(!aodE){ 
            AliError("No AOD Event");
            return;
        }
        fAODMCHeader = dynamic_cast<AliAODMCHeader *>(fInputEvent->FindListObject(AliAODMCHeader::StdBranchName()));
        if(!fAODMCHeader){ 
            AliError("No AliAODMCHeader");
            return;
        }
        fAODArrayMCInfo = dynamic_cast<TClonesArray *>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
        if(!fAODArrayMCInfo){ 
            AliError("No AOD MC particles");
            return;
        }
        fSignalCuts->SetMCAODInfo(fAODArrayMCInfo);
        // Background subtraction-------------------------------------------------------------------
        if (GetPlugin(kNonPhotonicElectron)) fBackgroundSubtraction->SetAODArrayMCInfo(fAODArrayMCInfo);
        //------------------------------------------------------------------------------------------
    }

    //printf("test2\n");

    fCentralityF = -1;
    if(!ReadCentrality()) fCentralityF = -1;

    // See if pile up and z in the range
    RejectionPileUpVertexRangeEventCut();

    // Protect agains missing 
    if(HasMCData()){
        //printf("Has MC data\n");
        fSignalCuts->SetMCEvent(fMCEvent);
        ProcessMC();  // Run the MC loop + MC QA in case MC Data are available
    }

    AliPIDResponse *pidResponse = fPIDResponse;
    if(!pidResponse){
        AliDebug(1, "Using default PID Response");
        pidResponse = AliHFEtools::GetDefaultPID(HasMCData(), fInputEvent->IsA() == AliESDEvent::Class());
    }
    fPID->SetPIDResponse(pidResponse);
    if(fPIDpreselect) fPIDpreselect->SetPIDResponse(pidResponse);

    // Background subtraction-------------------------------------------------------------------
    if(GetPlugin(kNonPhotonicElectron)) fBackgroundSubtraction->InitRun(fInputEvent,pidResponse);
    //------------------------------------------------------------------------------------------

    // Event loop
    if(IsAODanalysis()){
        ProcessAOD();
    } else {
        //ESD analysis TODO
    }
    // Done!!!
    PostData(1, fOutput);
    PostData(2, fQA);
}

//____________________________________________________________
void AliAnalysisTaskHFEMulti::Terminate(Option_t *){
    //
    // Terminate not implemented at the moment
    //
}

//_______________________________________________________________
Bool_t AliAnalysisTaskHFEMulti::IsEventInBinZero() {
    //
    //
    //

    //printf("test in IsEventInBinZero\n");
    if(!fInputEvent){
        AliError("Reconstructed Event not available");
        return kFALSE;
    }

    // check vertex
    const AliVVertex *vertex = fInputEvent->GetPrimaryVertex();
    if(!vertex) return kTRUE;
    //if(vertex) return kTRUE;

    // check tracks
    if(fInputEvent->GetNumberOfTracks()<=0) return kTRUE;
    //if(fInputEvent->GetNumberOfTracks()>0) return kTRUE;


    return kFALSE;

}
//____________________________________________________________
void AliAnalysisTaskHFEMulti::ProcessMC(){
    //
    // Runs the MC Loop (filling the container for the MC Cut Steps with the observables pt, eta and phi)
    // In case MC QA is on also MC QA loop is done
    //
    AliDebug(3, "Processing MC Information");
    Double_t eventContainer [5] = {0., 0., 0., 0., 0.};
    if(IsESDanalysis()) eventContainer[0] = fMCEvent->GetPrimaryVertex()->GetZ();
    else eventContainer[0] = fAODMCHeader->GetVtxZ();
    eventContainer[2] = fCentralityF;
    eventContainer[3] = fContributors;
    eventContainer[4] = fSPDtrkF;
    fVz = eventContainer[0];
    //printf("z position is %f\n",eventContainer[0]);
    //if(fCFM->CheckEventCuts(AliHFEcuts::kEventStepGenerated, fMCEvent)) 
    fCFM->GetEventContainer()->Fill(eventContainer,AliHFEcuts::kEventStepGenerated);
    Int_t nElectrons = 0;
    if(IsESDanalysis()){
        if(!((fIdentifiedAsPileUp) || (TMath::Abs(fVz) > fCuts->GetVertexRange()) || (fCentralityF < 0))){ //kStepMCGeneratedZOutNoPileUpCentralityFine
            if (HasMCData() && IsQAOn(kMCqa)) {
                AliDebug(2, "Running MC QA");

                if(fMCEvent->Stack()){
                    fMCQA->SetMCEvent(fMCEvent);
                    fMCQA->SetGenEventHeader(fMCEvent->GenEventHeader());
                    fMCQA->SetCentrality(fCentralityF);
                    fMCQA->SetPercentrality(static_cast<Int_t>(fCentralityPercent));
                    if(IsPbPb()) { fMCQA->SetPbPb();}
                    else
                    {
                        if(fisppMultiBin) fMCQA->SetPPMultiBin();
                        else fMCQA->SetPP();
                    }
                    fMCQA->Init();

                    fMCQA->GetMesonKine();

                    // loop over all tracks for decayed electrons
                    for (Int_t igen = 0; igen < fMCEvent->GetNumberOfTracks(); igen++){
                        TParticle* mcpart = fMCEvent->Stack()->Particle(igen);
                        if(!mcpart) continue;
                        fMCQA->GetQuarkKine(mcpart, igen, AliHFEmcQA::kCharm);
                        fMCQA->GetQuarkKine(mcpart, igen, AliHFEmcQA::kBeauty);
                        fMCQA->GetHadronKine(mcpart, AliHFEmcQA::kCharm);
                        fMCQA->GetHadronKine(mcpart, AliHFEmcQA::kBeauty);
                        fMCQA->GetDecayedKine(mcpart, AliHFEmcQA::kCharm,  AliHFEmcQA::kElectronPDG); // no accept cut
                        fMCQA->GetDecayedKine(mcpart, AliHFEmcQA::kBeauty, AliHFEmcQA::kElectronPDG); // no accept cut
                        fMCQA->GetDecayedKine(mcpart, AliHFEmcQA::kOthers, AliHFEmcQA::kElectronPDG); // no accept cut
                    }
                    //fMCQA->EndOfEventAna(AliHFEmcQA::kCharm);
                    //fMCQA->EndOfEventAna(AliHFEmcQA::kBeauty);
                }

            } // end of MC QA loop
        }
        // -----------------------------------------------------------------
        fCFM->SetMCEventInfo(fMCEvent);
        // fCFM->CheckEventCuts(AliCFManager::kEvtRecCuts, fESD);
    } else {
        fMCQA->SetMCArray(fAODArrayMCInfo);

        if(!((fIdentifiedAsPileUp) || (TMath::Abs(fVz) > fCuts->GetVertexRange()) || (fCentralityF < 0))){ //kStepMCGeneratedZOutNoPileUpCentralityFine
            if (HasMCData() && IsQAOn(kMCqa)) {
                AliDebug(2, "Running MC QA");

                fMCQA->SetCentrality(fCentralityF);
                fMCQA->SetPercentrality(static_cast<Int_t>(fCentralityPercent));

                if(IsPbPb()) { fMCQA->SetPbPb();}
                else
                {
                    if(fisppMultiBin) fMCQA->SetPPMultiBin();
                    else fMCQA->SetPP();
                }
                fMCQA->Init();

                //fMCQA->GetMesonKine();

                // loop over all tracks for decayed electrons
                AliAODMCParticle * mcpart;
                for (Int_t igen = 0; igen < fAODArrayMCInfo->GetEntriesFast(); igen++){
                    mcpart = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(igen));
                    if(!mcpart) continue;
                    fMCQA->GetDecayedKine(mcpart, AliHFEmcQA::kCharm,  AliHFEmcQA::kElectronPDG); // no accept cut
                    fMCQA->GetDecayedKine(mcpart, AliHFEmcQA::kBeauty, AliHFEmcQA::kElectronPDG); // no accept cut
                    fMCQA->GetDecayedKine(mcpart, AliHFEmcQA::kOthers, AliHFEmcQA::kElectronPDG); // no accept cut
                }

            } // end of MC QA loop
        }

        fCFM->SetMCEventInfo(fInputEvent);
    }
    // Run MC loop
    AliVParticle *mctrack = NULL;
    Int_t numberofmctracks = 0;
    if(IsESDanalysis()){
        numberofmctracks = fMCEvent->GetNumberOfTracks();
    }
    else {
        numberofmctracks = fAODArrayMCInfo->GetEntriesFast();
    }
    AliDebug(3, Form("Number of Tracks: %d",numberofmctracks));
    //printf("Number of MC track %d\n",numberofmctracks);
    for(Int_t imc = 0; imc <numberofmctracks; imc++){
        if(IsESDanalysis()) {
            if(!(mctrack = fMCEvent->GetTrack(imc))) continue;
        }
        else {
            if(!(mctrack = (AliVParticle *) fAODArrayMCInfo->At(imc))) continue;
        }
        //printf("Test in ProcessMC\n");
        AliDebug(4, "Next MC Track");
        if(ProcessMCtrack(mctrack)) nElectrons++;
    }

    // fCFM->CheckEventCuts(AliCFManager::kEvtRecCuts, fESD);
    //fQACollection->Fill("nElectron", nElectrons);
}

//____________________________________________________________
void AliAnalysisTaskHFEMulti::ProcessAOD(){
    //
    // Run Analysis in AOD Mode
    // Function is still in development
    //
    //printf("Process AOD\n");
    AliDebug(3, "Processing AOD Event");
    Double_t eventContainer[5];
    eventContainer[0] = 0.0;
    if(HasMCData()) eventContainer[0] = fVz;
    else {
        if(fInputEvent->GetPrimaryVertex()) eventContainer[0] = fInputEvent->GetPrimaryVertex()->GetZ();
    }
    eventContainer[1] = 1.; // No Information available in AOD analysis, assume all events have V0AND
    eventContainer[2] = fCentralityF; 
    eventContainer[3] = fContributors; 
    eventContainer[4] = fSPDtrkF; 

    //printf("value event container %f, %f, %f, %f\n",eventContainer[0],eventContainer[1],eventContainer[2],eventContainer[3]);

    AliAODEvent *fAOD = dynamic_cast<AliAODEvent *>(fInputEvent);
    if(!fAOD){
        AliError("AOD Event required for AOD Analysis");
        return;
    }
    if(!fAOD->GetPrimaryVertex()||TMath::Abs(fAOD->GetMagneticField())<0.001) return;

    //printf("Will fill\n");
    //
    fCFM->GetEventContainer()->Fill(eventContainer, AliHFEcuts::kEventStepRecNoCut);
    //printf("Fill\n");
    //
    if(fIdentifiedAsPileUp) return; 
    fCFM->GetEventContainer()->Fill(eventContainer, AliHFEcuts::kEventStepRecNoPileUp);

    //
    if(fIdentifiedAsOutInz) return;
    fCFM->GetEventContainer()->Fill(eventContainer, AliHFEcuts::kEventStepZRange);  

    //
    if(!fPassTheEventCut) return;
    fCFM->GetEventContainer()->Fill(eventContainer, AliHFEcuts::kEventStepReconstructed);
    //printf("pass\n");

    fContainer->NewEvent();

    Int_t nch = GetNcharged();

    fCFM->SetRecEventInfo(fAOD);
    

    if(!fExtraCuts){
        fExtraCuts = new AliHFEextraCuts("hfeExtraCuts","HFE Extra Cuts");
    }
    fExtraCuts->SetRecEventInfo(fAOD);

    // Get Number of contributors to the primary vertex for multiplicity-dependent correction
    Int_t ncontribVtx = 0;
    AliAODVertex *priVtx = fAOD->GetPrimaryVertex();
    if(priVtx){
        ncontribVtx = priVtx->GetNContributors();
    }

    // Look for kink mother
    Int_t numberofvertices = fAOD->GetNumberOfVertices();
    Double_t listofmotherkink[numberofvertices];
    Int_t numberofmotherkink = 0;
    for(Int_t ivertex=0; ivertex < numberofvertices; ivertex++) {
        AliAODVertex *aodvertex = fAOD->GetVertex(ivertex);
        if(!aodvertex) continue;
        if(aodvertex->GetType()==AliAODVertex::kKink) {
            AliAODTrack *mother = (AliAODTrack *) aodvertex->GetParent();
            if(!mother) continue;
            Int_t idmother = mother->GetID();
            listofmotherkink[numberofmotherkink] = idmother;
            //printf("ID %d\n",idmother);
            numberofmotherkink++;
        }
    }
    //printf("Number of kink mother in the events %d\n",numberofmotherkink);


    if(ncontribVtx > 0){
        fQACollection->Fill("NtrVsZ", fAOD->GetPrimaryVertex()->GetZ(), static_cast<Double_t>(fSPDtracklets));
        fQACollection->Fill("NtrCorrVsZ", fAOD->GetPrimaryVertex()->GetZ(), fSPDtrackletsCorr);
        fQACollection->Fill("NchCorrVsNtr",fSPDtrackletsCorr, static_cast<Double_t>(nch));
        fQACollection->Fill("NtrCorrVsNch",static_cast<Double_t>(nch),fSPDtrackletsCorr);
        fQACollection->Fill("VZEROVsZ", fAOD->GetPrimaryVertex()->GetZ(), static_cast<Double_t>(fSPDtracklets));
        fQACollection->Fill("VZEROCorrVsZ", fAOD->GetPrimaryVertex()->GetZ(), fSPDtrackletsCorr);
    }
    
    // Background subtraction-------------------------------------------------------------------
    if(fRefMulti>0){
        if (GetPlugin(kNonPhotonicElectron)) fBackgroundSubtraction->FillPoolAssociatedTracks(fInputEvent, fSPDtrkF);
    }else{
        if (GetPlugin(kNonPhotonicElectron)) fBackgroundSubtraction->FillPoolAssociatedTracks(fInputEvent, fCentralityF);
    }
    //------------------------------------------------------------------------------------------

    // Loop over tracks
    AliAODTrack *track = NULL;
    AliAODMCParticle *mctrack = NULL;
    Int_t nElectronCandidates = 0;
    Bool_t signal;

    //printf("Number of track %d\n",(Int_t) fAOD->GetNumberOfTracks());
    Bool_t kinkmother(kFALSE), kinkdaughter(kFALSE); Double_t kinkstatus(0);
//BEGIN TRACK LOOP
    for(Int_t itrack = 0; itrack < fAOD->GetNumberOfTracks(); itrack++){
        kinkmother=kFALSE;
        kinkdaughter=kFALSE;
        kinkstatus = 0.;
        track = (AliAODTrack *) fAOD->GetTrack(itrack); mctrack = NULL;
        if(!track) continue;

        for(int ivx = 0; ivx < numberofmotherkink; ivx++){
            if(track->GetID() == listofmotherkink[ivx]){
                kinkmother = kTRUE;
                break;
            }
        }
        AliAODVertex *pvx = track->GetProdVertex();
        if(pvx && (pvx->GetType() == AliAODVertex::kKink)) kinkdaughter = kTRUE;
        kinkstatus = 0.;
        if(kinkmother) kinkstatus = 1.;
        else if(kinkdaughter) kinkstatus = 2.;

        signal = kTRUE;
        if(HasMCData()){
            Int_t label = TMath::Abs(track->GetLabel());
            if(label && label < fAODArrayMCInfo->GetEntriesFast())
                mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(label));
            if(fFillSignalOnly && !fCFM->CheckParticleCuts(AliHFEcuts::kStepMCGenerated, mctrack)) signal = kFALSE;
            if(fRejectMCFakeTracks && IsMCFakeTrack(track)) signal = kFALSE;
        }

        if(fRefMulti>0){
            fVarManager->NewTrack(track, mctrack, fSPDtrkF, -1, signal);
        }else{
            fVarManager->NewTrack(track, mctrack, fCentralityF, -1, signal);
        }

        if(fFillNoCuts) {
            if(signal || !fFillSignalOnly){
                fVarManager->FillContainer(fContainer, "recTrackContReco", AliHFEcuts::kStepRecNoCut, kFALSE);
                fVarManager->FillContainer(fContainer, "recTrackContMC", AliHFEcuts::kStepRecNoCut, kTRUE);
            }
        }


        // RecKine: ITSTPC cuts  
        if(!ProcessCutStep(AliHFEcuts::kStepRecKineITSTPC, track)) continue;

        // Reject kink mother
        if(fRejectKinkMother) {
            Bool_t kinkmotherpass = kTRUE;
            for(Int_t ikinkmother = 0; ikinkmother < numberofmotherkink; ikinkmother++) {
                if(track->GetID() == listofmotherkink[ikinkmother]) {
                    kinkmotherpass = kFALSE;
                    continue;
                }
            }
            if(!kinkmotherpass) continue;
        }       

        // RecPrim
        if(!ProcessCutStep(AliHFEcuts::kStepRecPrim, track)) continue;

        // HFEcuts: ITS layers cuts
        if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsITS, track)) continue;

        // HFE cuts: TOF PID and mismatch flag
        if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsTOF, track)) continue;

        // HFE cuts: TPC PID cleanup
        if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsTPC, track)) continue;

        // HFEcuts: Nb of tracklets TRD0
        if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsTRD, track)) continue;

        // Fill correlation maps before PID
        if(signal && fContainer->GetCorrelationMatrix("correlationstepbeforePID")) {
            //printf("Fill correlation maps before PID\n");
            fVarManager->FillCorrelationMatrix(fContainer->GetCorrelationMatrix("correlationstepbeforePID"));
        }

        //printf("Will process to PID\n");

        // track accepted, do PID
        AliHFEpidObject hfetrack;
        hfetrack.SetAnalysisType(AliHFEpidObject::kAODanalysis);
        hfetrack.SetRecTrack(track);
        if(HasMCData()) hfetrack.SetMCTrack(mctrack);
        if(fRefMulti>0){
            hfetrack.SetCentrality(fSPDtrkF);
        }else{
            hfetrack.SetCentrality(fCentralityF);
        }
        hfetrack.SetMulitplicity(ncontribVtx); // for correction
        if(IsPbPb()) hfetrack.SetPbPb();
        else{
            if(IspPb()) hfetrack.SetpPb();
            else hfetrack.SetPP();
        }
        fPID->SetVarManager(fVarManager);
        if(!fPID->IsSelected(&hfetrack, fContainer, "recTrackCont", fPIDqa)) continue;   
        // we will do PID here as soon as possible

        // Background subtraction----------------------------------------------------------------------------------------------
        if (GetPlugin(kNonPhotonicElectron)) {
            Int_t indexmother = -1;
            Int_t mcsource = -1;  
            Int_t mcQAsource = -1;
            Double_t weightNonPhotonicFactor = 1.; 
            //printf("weight %f \n",weightNonPhotonicFactor);
            if(HasMCData() && mctrack){  
                mcsource = fBackgroundSubtraction->FindMother(TMath::Abs(track->GetLabel()),indexmother);
                if(fBackgroundSubtraction->GetLevelBack()>=0) {
                    if(fMCQA) {
                        fMCQA->SetCentrality(fCentralityF);
                        fMCQA->SetPercentrality(static_cast<Int_t>(fCentralityPercent));
                        mcQAsource = fMCQA->GetElecSource(mctrack, kTRUE);
                        fMCQA->SetContainerStep(2);
                        weightNonPhotonicFactor = TMath::Abs(fMCQA->GetWeightFactor(mctrack, fBackgroundSubtraction->GetLevelBack())); // positive:conversion e, negative: nonHFE 
                    }
                }
            }
            if(fRefMulti>0)
                fBackgroundSubtraction->LookAtNonHFE(itrack, track, fInputEvent, weightNonPhotonicFactor, fSPDtrkF, -1,mcsource, indexmother,mcQAsource);
            else
                fBackgroundSubtraction->LookAtNonHFE(itrack, track, fInputEvent, weightNonPhotonicFactor, fCentralityF, -1,mcsource, indexmother,mcQAsource);
        }
        //---------------------------------------------------------------------------------------------------------------------

        // Apply weight for background contamination
        //Double_t weightBackGround = 1.0;
        if(signal) {
            // Apply weight for background contamination
            if(fBackGroundFactorApply) {
                if(IsPbPb() && fCentralityF >= 0) fWeightBackGround =  fkBackGroundFactorArray[fCentralityF >= 0 ? fCentralityF : 0]->Eval(TMath::Abs(track->P()));
                else    fWeightBackGround =  fkBackGroundFactorArray[0]->Eval(TMath::Abs(track->P())); // pp case

                if(fWeightBackGround < 0.0) fWeightBackGround = 0.0;
                else if(fWeightBackGround > 1.0) fWeightBackGround = 1.0;
                // weightBackGround as special weight
                fVarManager->FillContainer(fContainer, "hadronicBackground", 1, kFALSE, fWeightBackGround);
            }
            fVarManager->FillCorrelationMatrix(fContainer->GetCorrelationMatrix("correlationstepafterPID"));
        }

        nElectronCandidates++;    

    }
//END TRACK LOOP

    // Background subtraction-------------------------------------------------------------------
    if(fRefMulti>0){
        if (GetPlugin(kNonPhotonicElectron)) fBackgroundSubtraction->CountPoolAssociated(fInputEvent, fSPDtrkF);
    }else{
        if (GetPlugin(kNonPhotonicElectron)) fBackgroundSubtraction->CountPoolAssociated(fInputEvent, fCentralityF);
    }
    //------------------------------------------------------------------------------------------

    //fQACollection->Fill("nElectronTracksEvent", nElectronCandidates);
}

//____________________________________________________________
Bool_t AliAnalysisTaskHFEMulti::ProcessMCtrack(AliVParticle *track){
    //
    // Filter the Monte Carlo Track
    // Additionally Fill a THnSparse for Signal To Background Studies
    // Works for AOD and MC analysis Type
    //
    if(fRefMulti>0){
        fVarManager->NewTrack(track, NULL, fSPDtrkF, -1, kTRUE);
    }else{
        fVarManager->NewTrack(track, NULL, fCentralityF, -1, kTRUE);
    }
    //printf("Is primary %d\n",((Int_t)track->IsPrimary()));


    Double_t vertex[3] = {0.,0.,0.}; // Production vertex cut to mask gammas which are NOT supposed to have hits in the first ITS layer(s)
    if(IsESDanalysis()){
        AliMCParticle *mctrack = dynamic_cast<AliMCParticle *>(track);
        if(mctrack){
            vertex[0] = mctrack->Particle()->Vx();
            vertex[1] = mctrack->Particle()->Vy();
        }
    } else {
        AliAODMCParticle *aodmctrack = dynamic_cast<AliAODMCParticle *>(track);
        if(aodmctrack) aodmctrack->XvYvZv(vertex);
    }

    //printf("MC Generated\n");
    if(!fCFM->CheckParticleCuts(AliHFEcuts::kStepMCGenerated, track)) return kFALSE;
    //printf("MC Generated pass\n");
    fVarManager->FillContainer(fContainer, "MCTrackCont", AliHFEcuts::kStepMCGenerated, kFALSE);

    // Step GeneratedZOutNoPileUp
    if((fIdentifiedAsPileUp) || (TMath::Abs(fVz) > fCuts->GetVertexRange()) || (fCentralityF < 0)) return kFALSE;
    fVarManager->FillContainer(fContainer, "MCTrackCont", AliHFEcuts::kStepMCGeneratedZOutNoPileUpCentralityFine, kFALSE);
    //printf("In ProcessMCtrack %f\n",fCentralityF);

    // Step Generated Event Cut
    if(!fPassTheEventCut) return kFALSE;
    fVarManager->FillContainer(fContainer, "MCTrackCont", AliHFEcuts::kStepMCGeneratedEventCut, kFALSE);

    if(IsESDanalysis()){
        if(!fCFM->CheckParticleCuts(AliHFEcuts::kStepMCInAcceptance, track)) return kFALSE;
        fVarManager->FillContainer(fContainer, "MCTrackCont", AliHFEcuts::kStepMCInAcceptance, kFALSE);
    }
    return kTRUE;
}

//____________________________________________________________
Bool_t AliAnalysisTaskHFEMulti::PreSelectTrack(AliESDtrack *track) const {
    //
    // Preselect tracks
    //


    Bool_t survived = kTRUE;

    if(fCutspreselect) {
        //printf("test preselect\n");
        if(!fCutspreselect->IsSelected(track)) survived=kFALSE;
    }
    //printf("survived %d\n",(Int_t)survived);

    if(survived && fPIDpreselect){
        // Apply PID
        AliHFEpidObject hfetrack;
        hfetrack.SetAnalysisType(AliHFEpidObject::kESDanalysis);
        hfetrack.SetRecTrack(track);
        if(!fPIDpreselect->IsSelected(&hfetrack)) {
            //printf("Did not pass AliHFEcuts::kPID\n");
            survived = kFALSE;
        }
        //else printf("Pass AliHFEcuts::kPID\n");
    }

    return survived; 

}
//____________________________________________________________
void AliAnalysisTaskHFEMulti::MakeEventContainer(){
    //
    // Create the event container for the correction framework and link it
    // 1st bin: Vertex z-position
    // 2nd bin: V0AND decision (normalization to sigma_inel)
    // 3rd bin: Centrality class (for pp defined as number of contributors in vertex.)
    // 4th bin: Number of contributors > 0
    //

    const Int_t kNvar = 5;  // number of variables on the grid: 
    Int_t nBins[kNvar] =     {100,  2,   11,  2,  11};
    Double_t binMin[kNvar] = {-10., 0.,  0., 0.,   0.};
    Double_t binMax[kNvar] = { 10., 2., 11., 2., 11.};

    AliCFContainer *evCont = new AliCFContainer("eventContainer", "Container for events", AliHFEcuts::kNcutStepsEvent, kNvar, nBins);

    Double_t *vertexBins = AliHFEtools::MakeLinearBinning(nBins[0], binMin[0], binMax[0]);
    Double_t *v0andBins = AliHFEtools::MakeLinearBinning(nBins[1], binMin[1], binMax[1]);
    Double_t *centralityBins = AliHFEtools::MakeLinearBinning(nBins[2], binMin[2], binMax[2]);
    Double_t *contributorsBins = AliHFEtools::MakeLinearBinning(nBins[3], binMin[3], binMax[3]);
    Double_t *trackletsBins = AliHFEtools::MakeLinearBinning(nBins[4], binMin[4], binMax[4]);
    evCont->SetBinLimits(0, vertexBins);
    evCont->SetBinLimits(1, v0andBins);
    evCont->SetBinLimits(2, centralityBins);
    evCont->SetBinLimits(3, contributorsBins);
    evCont->SetBinLimits(4, trackletsBins);
    delete[] vertexBins; delete[] v0andBins; delete[] centralityBins; delete[] contributorsBins; delete[] trackletsBins;

    fCFM->SetEventContainer(evCont);
}

//____________________________________________________________
void AliAnalysisTaskHFEMulti::MakeParticleContainer(){
    //
    // Create the particle container for the correction framework manager and 
    // link it
    //
    if(!fContainer) fContainer = new AliHFEcontainer("trackContainer");
    fVarManager->DefineVariables(fContainer);

    // Create Correction Framework containers
    fContainer->CreateContainer("MCTrackCont", "Track Container filled with MC information", AliHFEcuts::kNcutStepsMCTrack);
    fContainer->CreateContainer("recTrackContReco", "Track Container filled with MC information", AliHFEcuts::kNcutStepsRecTrack + fPID->GetNumberOfPIDdetectors());
    fContainer->CreateContainer("recTrackContMC", "Track Container filled with MC information", AliHFEcuts::kNcutStepsRecTrack + fPID->GetNumberOfPIDdetectors());

    fContainer->CreateContainer("hadronicBackground", "Container for Hadronic Background", 3);


    fContainer->CreateCorrelationMatrix("correlationstepafterPID","THnSparse with correlations");
    if(!fVarManager->IsVariableDefined("centrality")) {
        fContainer->CreateCorrelationMatrix("correlationstepbeforePID","THnSparse with correlations");
        fContainer->CreateCorrelationMatrix("correlationstepafterTOF","THnSparse with correlations");
    }

    // Define the step names
    for(UInt_t istep = 0; istep < AliHFEcuts::kNcutStepsMCTrack; istep++){
        fContainer->SetStepTitle("MCTrackCont", AliHFEcuts::MCCutName(istep), istep);
    }
    for(UInt_t istep = 0; istep < AliHFEcuts::kNcutStepsRecTrack; istep++){
        fContainer->SetStepTitle("recTrackContReco", AliHFEcuts::RecoCutName(istep), istep);
        fContainer->SetStepTitle("recTrackContMC", AliHFEcuts::RecoCutName(istep), istep);
    }
    for(UInt_t ipid = 0; ipid < fPID->GetNumberOfPIDdetectors(); ipid++){
        fContainer->SetStepTitle("recTrackContReco", fPID->SortedDetectorName(ipid), AliHFEcuts::kNcutStepsRecTrack + ipid);
        fContainer->SetStepTitle("recTrackContMC", fPID->SortedDetectorName(ipid), AliHFEcuts::kNcutStepsRecTrack + ipid);
    }
}

//____________________________________________________________
void AliAnalysisTaskHFEMulti::PrintStatus() const {
    //
    // Print Analysis status
    //
    printf("\n\tAnalysis Settings\n\t========================================\n\n");
    printf("\tPrimary Vertex resolution: %s\n", GetPlugin(kPriVtx) ? "YES" : "NO");
    printf("\n");
    printf("\tParticle Identification Detectors:\n");
    fPID->PrintStatus();
    printf("\n");
    printf("\tQA: \n");
    printf("\t\tPID: %s\n", IsQAOn(kPIDqa) ? "YES" :  "NO");
    printf("\t\tCUTS: %s\n", (fCuts != NULL && fCuts->IsQAOn()) ? "YES" : "NO");
    printf("\t\tMC: %s\n", IsQAOn(kMCqa) ? "YES" : "NO");
    printf("\n");
}

//____________________________________________________________
Bool_t AliAnalysisTaskHFEMulti::FillProductionVertex(const AliVParticle * const track) const{
    //
    // Find the production vertex of the associated MC track
    //
    if(!fMCEvent) return kFALSE;
    const AliVParticle *mctrack = NULL;
    TString objectType = track->IsA()->GetName();
    if(objectType.CompareTo("AliESDtrack") == 0 || objectType.CompareTo("AliAODTrack") == 0){
        // Reconstructed track
        mctrack = fMCEvent->GetTrack(TMath::Abs(track->GetLabel()));
    } else {
        // MCParticle
        mctrack = track;
    }

    if(!mctrack) return kFALSE;

    Double_t xv = 0.0;
    Double_t yv = 0.0;

    if(TString(mctrack->IsA()->GetName()).CompareTo("AliMCParticle") == 0){
        // case MCParticle
        const AliMCParticle *mcpart = dynamic_cast<const AliMCParticle *>(mctrack);
        if(mcpart){
            xv =  mcpart->Xv();
            yv =  mcpart->Yv();
        }
    } else {
        // case AODMCParticle
        const AliAODMCParticle *mcpart = dynamic_cast<const AliAODMCParticle *>(mctrack);
        if(mcpart){
            xv =  mcpart->Xv();
            yv =  mcpart->Yv();
        }
    }

    //printf("xv %f, yv %f\n",xv,yv);
    fQACollection->Fill("radius", TMath::Abs(xv),TMath::Abs(yv));

    return kTRUE;

}
//__________________________________________
void AliAnalysisTaskHFEMulti::SwitchOnPlugin(Int_t plug){
    //
    // Switch on Plugin
    // Available:
    //  - Primary vertex studies
    //  - Secondary vertex Studies
    //  - Post Processing
    //
    switch(plug){
        case kPriVtx: SETBIT(fPlugins, plug); break;
        case kNonPhotonicElectron: SETBIT(fPlugins, plug); break; 
        default: AliError("Unknown Plugin");
    };
}
//__________________________________________
Bool_t AliAnalysisTaskHFEMulti::ProcessCutStep(Int_t cutStep, AliVParticle *track){
    //
    // Check single track cuts for a given cut step
    // Fill the particle container
    //
    const Int_t kMCOffset = AliHFEcuts::kNcutStepsMCTrack;
    if(!fCFM->CheckParticleCuts(cutStep + kMCOffset, track)) return kFALSE;
    if(fVarManager->IsSignalTrack()) {
        fVarManager->FillContainer(fContainer, "recTrackContReco", cutStep, kFALSE);
        fVarManager->FillContainer(fContainer, "recTrackContMC", cutStep, kTRUE);
    }
    return kTRUE;
}
//___________________________________________________
Bool_t AliAnalysisTaskHFEMulti::ReadCentrality() {
    //
    // Recover the centrality of the event from ESD or AOD
    //

    Float_t fCentralityLimitstemp[12];
    Float_t fCentralityLimitsdefault[12]= {0.,5.,10., 20., 30., 40., 50., 60.,70.,80., 90., 100.00001};
    if(!fPbPbUserCentralityBinning) memcpy(fCentralityLimitstemp,fCentralityLimitsdefault,sizeof(fCentralityLimitsdefault));
    else memcpy(fCentralityLimitstemp,fCentralityLimits,sizeof(fCentralityLimitsdefault));


    Int_t bin = -1;
    if(IsPbPb()||IspPb()) {
        // Centrality
        AliCentrality *centrality = fInputEvent->GetCentrality();
        fCentralityPercent = centrality->GetCentralityPercentile(fCentralityEstimator.Data());
        //printf("centrality %f\n",fCentralityPercent);

        for(Int_t ibin = 0; ibin < 11; ibin++){
            if(fCentralityPercent >= fCentralityLimitstemp[ibin] && fCentralityPercent < fCentralityLimitstemp[ibin+1]){
                bin = ibin;
                //printf("test bin %f, low %f, high %f, %d\n",fCentralityPercent,fCentralityLimitstemp[ibin],fCentralityLimitstemp[ibin+1],ibin);
                break;
            }
        }

        if(bin == -1) bin = 11; // Overflow
    } else {
        // PP: Tracklet multiplicity, use common definition
        Int_t itsMultiplicity = GetITSMultiplicity(fInputEvent);
        Int_t multiplicityLimits[8] = {0, 1, 9, 17, 25, 36, 60, 500};
        for(Int_t ibin = 0; ibin < 7; ibin++){  
            if(itsMultiplicity >= multiplicityLimits[ibin] && itsMultiplicity < multiplicityLimits[ibin + 1]){
                bin = ibin;
                break;
            }
        }
        if(bin == -1) bin = 7;  // Overflow
    }
    fCentralityF = bin;
    AliDebug(2, Form("Centrality class %d\n", fCentralityF));


    // contributors, to be outsourced
    const AliVVertex *vtx;
    if(IsAODanalysis()){
        AliAODEvent *fAOD = dynamic_cast<AliAODEvent *>(fInputEvent);
        if(!fAOD){
            AliError("AOD Event required for AOD Analysis");
            return kFALSE;
        }
        vtx = fAOD->GetPrimaryVertex();
    } else {
        AliESDEvent *fESD = dynamic_cast<AliESDEvent *>(fInputEvent);
        if(!fESD){
            AliError("ESD Event required for ESD Analysis");
            return kFALSE;
        }
        vtx = fESD->GetPrimaryVertex() ;
    }
    if(!vtx){ 
        fContributors = 0.5;
        return kFALSE;
    }
    else {
        Int_t contributorstemp = vtx->GetNContributors();
        if( contributorstemp <=  0) {
            fContributors =  0.5;
            //printf("Number of contributors %d and vz %f\n",contributorstemp,vtx->GetZ());
        }
        else fContributors = 1.5;
        //printf("Number of contributors %d\n",contributorstemp);
    }




    fSPDtracklets=GetITSMultiplicity(fInputEvent);
    if(IspPb() && fRefMulti>0){
        Int_t period = (fInputEvent->GetRunNumber() < 195484)?0:1; //TODO get correct period with nicer function
        const Int_t local_mbins = 8;
        Int_t multiplicityLimits[local_mbins] = {0, 1, 22, 29, 35, 44, 70, 200};

        if(fMultiEstimatorSystem == kVZERO){
            period += 2;
            Int_t tmpMulti[local_mbins] = {0, 91, 133, 173, 227, 789, 800,801};
            memcpy(multiplicityLimits,tmpMulti,sizeof(multiplicityLimits));
            fSPDtracklets = GetVZEROMultiplicity(fInputEvent);
        }
        //fSPDtrackletsCorr = GetCorrectedNtracklets(fMultEstimatorAvg[period], fSPDtracklets, vtx->GetZ(), fRefMulti); 
        fSPDtrackletsCorr =static_cast<Int_t>(GetCorrectedNtracklets(fMultEstimatorAvg[period], fSPDtracklets, vtx->GetZ(), fRefMulti));
        //Weight MC Ntr for difference to data
        if(HasMCData() && fMCNtrWeight){
            //Double_t loc_weight = fMCNtrWeight->GetBinContent(fMCNtrWeight->FindBin(fSPDtrackletsCorr));
            //fSPDtrackletsCorr *= loc_weight;
        }

        if(fSPDtrackletsCorr < multiplicityLimits[0]){
            //store in underflow bin
            fSPDtrkF=-1;
            return kTRUE;
        }

        for(Int_t ibin = 0; ibin < local_mbins-1; ibin++){  
            if(fSPDtrackletsCorr >= multiplicityLimits[ibin] && fSPDtrackletsCorr < multiplicityLimits[ibin + 1]){
                bin = ibin;
                break;
            }else bin = 11; //overflow
        }
        fSPDtrkF=bin;
    }

    return kTRUE;
}

//___________________________________________________
Int_t AliAnalysisTaskHFEMulti::GetITSMultiplicity(AliVEvent *ev){
    //
    // Definition of the Multiplicity according to the JPSI group (F. Kramer)
    //
    Int_t nTracklets = 0;
    Int_t nAcc = 0;
    Double_t etaRange = 1.0;

    if (ev->IsA() == AliAODEvent::Class()) {
        AliAODTracklets *tracklets = ((AliAODEvent*)ev)->GetTracklets();
        nTracklets = tracklets->GetNumberOfTracklets();
        for (Int_t nn = 0; nn < nTracklets; nn++) {
            Double_t theta = tracklets->GetTheta(nn);
            Double_t eta = -TMath::Log(TMath::Tan(theta/2.0));
            if (TMath::Abs(eta) < etaRange) nAcc++;
        }
    } else if (ev->IsA() == AliESDEvent::Class()) {
        nTracklets = ((AliESDEvent*)ev)->GetMultiplicity()->GetNumberOfTracklets();
        for (Int_t nn = 0; nn < nTracklets; nn++) {
            Double_t eta = ((AliESDEvent*)ev)->GetMultiplicity()->GetEta(nn);
            if (TMath::Abs(eta) < etaRange) nAcc++;
        }
    } else return -1;

    return nAcc;
}


//___________________________________________________
Int_t AliAnalysisTaskHFEMulti::GetVZEROMultiplicity(AliVEvent *ev){
    //
    // Definition of the Multiplicity according to the D2H group (Zaida)
    //
    if (ev->IsA() == AliAODEvent::Class()) {
        AliAODVZERO *vzeroAOD = dynamic_cast<AliAODVZERO *>( dynamic_cast<AliAODEvent *>(ev)->GetVZEROData());
        return static_cast<Int_t>(vzeroAOD->GetMTotV0A());
    } else if (ev->IsA() == AliESDEvent::Class()) {
        return -1;
    } else return -1;

}

//___________________________________________________
void AliAnalysisTaskHFEMulti::RejectionPileUpVertexRangeEventCut() {
    //
    // Recover the centrality of the event from ESD or AOD
    //
    if(IsAODanalysis()){

        AliAODEvent *fAOD = dynamic_cast<AliAODEvent *>(fInputEvent);
        if(!fAOD){
            AliError("AOD Event required for AOD Analysis");
            return;
        }
        // PileUp
        fIdentifiedAsPileUp = kFALSE;
        if(fRemovePileUp && fAOD->IsPileupFromSPD()) fIdentifiedAsPileUp = kTRUE; 
        // Z vertex
        fIdentifiedAsOutInz = kFALSE;
        //printf("Z vertex %f and out %f\n",fAOD->GetPrimaryVertex()->GetZ(),fCuts->GetVertexRange());
        if(TMath::Abs(fAOD->GetPrimaryVertex()->GetZ()) > fCuts->GetVertexRange()) fIdentifiedAsOutInz = kTRUE;
        // Event Cut
        fPassTheEventCut = kTRUE;
        if(!fCFM->CheckEventCuts(AliHFEcuts::kEventStepReconstructed, fAOD)) fPassTheEventCut = kFALSE; 


    } else {

        AliDebug(3, "Processing ESD Centrality");
        AliESDEvent *fESD = dynamic_cast<AliESDEvent *>(fInputEvent);
        if(!fESD){
            AliError("ESD Event required for ESD Analysis");
            return;
        }
        // PileUp
        fIdentifiedAsPileUp = kFALSE;
        if(fRemovePileUp && fESD->IsPileupFromSPD()) fIdentifiedAsPileUp = kTRUE; 



        // Z vertex
        fIdentifiedAsOutInz = kFALSE;
        Bool_t findvertex = kTRUE;
        const AliESDVertex* vtxESD = fESD->GetPrimaryVertex();
        if((!vtxESD) || (vtxESD->GetNContributors() <= 0)) findvertex = kFALSE;
        if(findvertex) {
            if(TMath::Abs(vtxESD->GetZ()) > fCuts->GetVertexRange()) fIdentifiedAsOutInz = kTRUE;
        }

        //Event Cut
        fPassTheEventCut = kTRUE;
        if(!fCFM->CheckEventCuts(AliHFEcuts::kEventStepReconstructed, fESD)) fPassTheEventCut = kFALSE;   

    }

}


//___________________________________________________
Bool_t AliAnalysisTaskHFEMulti::IsMCFakeTrack(const AliVTrack *const trk) const {
    //
    // Check whether track is MC Fake track using the sign of the track label
    //
    return trk->GetLabel() < 0;
}
//____________________________________________________________________________
Double_t AliAnalysisTaskHFEMulti::GetCorrectedNtracklets(TProfile* estimatorAvg, Double_t uncorrectedNacc, Double_t vtxZ, Double_t refMult) {
//Example from D2H:
  
   //Correct the number of accepted tracklets based on a TProfile Hist
  
  

  if(TMath::Abs(vtxZ)>10.0){
    //    printf("ERROR: Z vertex out of range for correction of multiplicity\n");
    return uncorrectedNacc;
  }

  if(!estimatorAvg){
    printf("ERROR: Missing TProfile for correction of multiplicity\n");
    return uncorrectedNacc;
  }

  Double_t localAvg = estimatorAvg->GetBinContent(estimatorAvg->FindBin(vtxZ));

  Double_t deltaM = uncorrectedNacc*(refMult/localAvg - 1);

  Double_t correctedNacc = uncorrectedNacc + (deltaM>0 ? 1 : -1) * gRandom->Poisson(TMath::Abs(deltaM));

  return correctedNacc;
}

//____________________________________________________________________________
Int_t AliAnalysisTaskHFEMulti::GetNcharged(){
    //counts all tracks in eta<1 with charge!=0

    Int_t Nch = 0;

    if(!HasMCData()) return Nch; // if no MC info return 0

    // loop over all tracks 
    for (Int_t igen = 0; igen < fAODArrayMCInfo->GetEntriesFast(); igen++){
        AliAODMCParticle *mctrack=(AliAODMCParticle*)fAODArrayMCInfo->UncheckedAt(igen);
        Int_t charge = mctrack->Charge();
        Double_t eta = mctrack->Eta();
        Bool_t isPhysPrim = mctrack->IsPhysicalPrimary();
        if(charge!=0){
            if(eta > -1.0 && eta < 1.0){
                if(isPhysPrim){
                    Nch++;
                }
            }
        }
    }
    return Nch;
}
