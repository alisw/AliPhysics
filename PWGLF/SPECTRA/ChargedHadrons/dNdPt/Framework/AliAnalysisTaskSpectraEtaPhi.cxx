#include <iostream>
#include "TChain.h"
#include "TH1F.h"
#include "TList.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TGeoGlobalMagField.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliAODEvent.h"
#include "AliHeader.h"
#include "AliMCEvent.h"
#include "AliGenEventHeader.h"
#include "AliESDtrackCuts.h"
#include "AlidNdPtTools.h"
#include "AliAnalysisTaskMKBase.h"
#include "AliAnalysisTaskSpectraEtaPhi.h"

class AliAnalysisTaskSpectraEtaPhi;

using namespace std;

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSpectraEtaPhi)
/// \endcond
//_____________________________________________________________________________

AliAnalysisTaskSpectraEtaPhi::AliAnalysisTaskSpectraEtaPhi()
    : AliAnalysisTaskMKBase()
    , fHistEffContNCluster(0)
    , fHistEffContZ(0)
    , fHistEffContEta(0)
    , fHistEffContPhi(0)
    , fHistTrackNCluster(0)
    , fHistTrackZ(0)
    , fHistTrackEta(0)
    , fHistTrackPhi(0)
    , fHistEvent(0)
{
    // default contructor
}

//_____________________________________________________________________________

AliAnalysisTaskSpectraEtaPhi::AliAnalysisTaskSpectraEtaPhi(const char* name)
    : AliAnalysisTaskMKBase(name)
    , fHistEffContNCluster(0)
    , fHistEffContZ(0)
    , fHistEffContEta(0)
    , fHistEffContPhi(0)
    , fHistTrackNCluster(0)
    , fHistTrackZ(0)
    , fHistTrackEta(0)
    , fHistTrackPhi(0)
    , fHistEvent(0)
{
    // constructor
}

//_____________________________________________________________________________

AliAnalysisTaskSpectraEtaPhi::~AliAnalysisTaskSpectraEtaPhi()
{
    // destructor
}

//_____________________________________________________________________________

void AliAnalysisTaskSpectraEtaPhi::AddOutput()
{
    // MC
    AddAxis("cent");
    AddAxis("MCpT","ptfewpatrick");
    AddAxis("NClusterPID", "N_{Cluster, PID}", 50, -0.5, 200.5);
    AddAxis("ESDTrackCut", 9, -1.5,7.5);//default, nogeo, nogeonogold, geoCut116, geoCut117, geoCut118 geoCut119, tpconlyminimal
    AddAxis("MCinfo",4,-0.5,3.5);  // 0=prim, 1=decay 2=material, 3=genprim
    fHistEffContNCluster = (THnF*)CreateHist("fHistEffContNCluster");
    fOutputList->Add(fHistEffContNCluster);
    
    AddAxis("cent");
    AddAxis("MCpT","ptfewpatrick");
    AddAxis("z","Z",60,-30.5,+29.5);
    AddAxis("ESDTrackCut", 9, -1.5,7.5);//default, nogeo, nogeonogold, geoCut116, geoCut117, geoCut118 geoCut119, tpconlyminimal
    AddAxis("MCinfo",4,-0.5,3.5);  // 0=prim, 1=decay 2=material, 3=genprim
    fHistEffContZ = (THnF*)CreateHist("fHistEffContZ");
    fOutputList->Add(fHistEffContZ);
    
    AddAxis("cent");
    AddAxis("MCpT","ptfewpatrick");
    AddAxis("eta","#eta",10,-1.,+1.);
    AddAxis("ESDTrackCut", 9, -1.5,7.5);//default, nogeo, nogeonogold, geoCut116, geoCut117, geoCut118 geoCut119, tpconlyminimal
    AddAxis("MCinfo",4,-0.5,3.5);  // 0=prim, 1=decay 2=material, 3=genprim
    AddAxis("MCQ",3,-1.5,1.5);
    fHistEffContEta = (THnF*)CreateHist("fHistEffContEta");
    fOutputList->Add(fHistEffContEta);
    
    AddAxis("cent");
    AddAxis("MCpT","ptfewpatrick");
    AddAxis("phi","#phi",72,0.,2*TMath::Pi());
    AddAxis("ESDTrackCut", 9, -1.5,7.5);//default, nogeo, nogeonogold, geoCut116, geoCut117, geoCut118 geoCut119, tpconlyminimal
    AddAxis("MCinfo",4,-0.5,3.5);  // 0=prim, 1=decay 2=material, 3=genprim
    fHistEffContPhi = (THnF*)CreateHist("fHistEffContPhi");
    fOutputList->Add(fHistEffContPhi);
    
    // data
    AddAxis("cent");
    AddAxis("pt", "ptfewpatrick");
    AddAxis("NClusterPID", "N_{Cluster, PID}", 50, -0.5, 200.5);
    AddAxis("ESDTrackCut", 9, -1.5,7.5);//default, nogeo, nogeonogold, geoCut116, geoCut117, geoCut118 geoCut119, tpconlyminimal
    fHistTrackNCluster = (THnF*)CreateHist("fHistTrackNCluster");
    fOutputList->Add(fHistTrackNCluster);
    
    AddAxis("cent");
    AddAxis("pt", "ptfewpatrick");
    AddAxis("z","Z",60,-30.5,+29.5);
    AddAxis("ESDTrackCut", 9, -1.5,7.5);//default, nogeo, nogeonogold, geoCut116, geoCut117, geoCut118 geoCut119, tpconlyminimal
    fHistTrackZ = (THnF*)CreateHist("fHistTrackZ");
    fOutputList->Add(fHistTrackZ);
    
    AddAxis("cent");
    AddAxis("pt", "ptfewpatrick");
    AddAxis("eta","#eta",10,-1.,+1.);
    AddAxis("ESDTrackCut", 9, -1.5,7.5);//default, nogeo, nogeonogold, geoCut116, geoCut117, geoCut118 geoCut119, tpconlyminimal
    AddAxis("Q",3,-1.5,1.5);
    fHistTrackEta = (THnF*)CreateHist("fHistTrackEta");
    fOutputList->Add(fHistTrackEta);
    
    AddAxis("cent");
    AddAxis("pt", "ptfewpatrick");
    AddAxis("phi","#phi",72,0.,2*TMath::Pi());
    AddAxis("ESDTrackCut", 9, -1.5,7.5);//default, nogeo, nogeonogold, geoCut116, geoCut117, geoCut118 geoCut119, tpconlyminimal
    fHistTrackPhi = (THnF*)CreateHist("fHistTrackPhi");
    fOutputList->Add(fHistTrackPhi);
    
    // Event hists
    AddAxis("cent");
    AddAxis("nAcc","mult6kfine");
    AddAxis("zV",8,-20,20);
    fHistEvent = (THnF*)CreateHist("fHistEvent");
    fOutputList->Add(fHistEvent);
    
}

//_____________________________________________________________________________

Bool_t AliAnalysisTaskSpectraEtaPhi::IsEventSelected()
{
    return fIsAcceptedAliEventCuts;
}

//_____________________________________________________________________________

void AliAnalysisTaskSpectraEtaPhi::AnaEvent()
{
   
   LoopOverAllTracks();
   if (fIsMC) LoopOverAllParticles();
   
   FillHist(fHistEvent, fMultPercentileV0M, static_cast<double>(fNTracksAcc), static_cast<double>(fZv));
   
}

//_____________________________________________________________________________

void AliAnalysisTaskSpectraEtaPhi::AnaTrack(Int_t flag)
{
//    if (!fAcceptTrackM) return;
    for(int i=0; i<8;++i){
        if(fAcceptTrack[i])
        {
            FillHist(fHistTrackNCluster, fMultPercentileV0M, fPt, fTPCSignalN, static_cast<Double_t>(i));
            FillHist(fHistTrackZ, fMultPercentileV0M, fPt, fZInner, static_cast<Double_t>(i));
            FillHist(fHistTrackEta, fMultPercentileV0M, fPt, fEta, static_cast<Double_t>(i), static_cast<Double_t>(fChargeSign));
            FillHist(fHistTrackPhi, fMultPercentileV0M, fPt, fPhi, static_cast<Double_t>(i));
        }
    }
}

//_____________________________________________________________________________

void AliAnalysisTaskSpectraEtaPhi::AnaTrackMC(Int_t flag)
{
//    if (!fAcceptTrackM) return;
     
    if (fMCParticleType==AlidNdPtTools::kOther) { Log("RecTrack.PDG.",fMCPDGCode); }
    if (TMath::Abs(fMCQ > 1)) { Log("RecTrack.Q>1.PDG.",fMCPDGCode); }
    
    for(int i=0; i<8;++i){
        if(fAcceptTrack[i])
        {
            FillHist(fHistEffContNCluster, fMultPercentileV0M, fPt, fTPCSignalN, static_cast<Double_t>(i), fMCProdcutionType);
            FillHist(fHistEffContZ, fMultPercentileV0M, fPt, fZInner, static_cast<Double_t>(i), fMCProdcutionType);
            FillHist(fHistEffContEta, fMultPercentileV0M, fPt, fEta, static_cast<Double_t>(i), fMCProdcutionType, static_cast<Double_t>(fChargeSign));
            FillHist(fHistEffContPhi, fMultPercentileV0M, fPt, fPhi, static_cast<Double_t>(i), fMCProdcutionType);
        }
    }

}

//_____________________________________________________________________________

void AliAnalysisTaskSpectraEtaPhi::AnaParticleMC(Int_t flag)
{
    if (!fMCisPrim) return;
    if (!fMCIsCharged) return;
    if (TMath::Abs(fMCEta) > 0.8) return;
    
    if (fMCParticleType==AlidNdPtTools::kOther) { Log("GenPrim.PDG.",fMCPDGCode); }
    if (TMath::Abs(fMCQ > 1)) { Log("GenPrim.Q>1.PDG.",fMCPDGCode); }
    
    FillHist(fHistEffContNCluster, fMultPercentileV0M, fPt, fTPCSignalN, -1., 3.);
    FillHist(fHistEffContZ, fMultPercentileV0M, fPt, fZInner, -1., 3.);
    FillHist(fHistEffContEta, fMultPercentileV0M, fPt, fEta, -1., 3., static_cast<double>(fChargeSign));
    FillHist(fHistEffContPhi, fMultPercentileV0M, fPt, fPhi, -1., 3.);

}

//_____________________________________________________________________________

AliAnalysisTaskSpectraEtaPhi* AliAnalysisTaskSpectraEtaPhi::AddTaskSpectra(const char* name, const char* outfile)
{
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskSpectraEtaPhi", "No analysis manager to connect to.");
        return 0;
    }

    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddTaskSpectraEtaPhi", "This task requires an input event handler");
        return NULL;
    }
    
    // Setup output file
    //===========================================================================
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":";
    fileName += name;  // create a subfolder in the file
    if (outfile) { // if a finename is given, use that one
        fileName = TString(outfile);
    }
    

    // create the task
    //===========================================================================
    AliAnalysisTaskSpectraEtaPhi *task = new AliAnalysisTaskSpectraEtaPhi(name);
    if (!task) { return 0; }
    
    // configure the task
    //===========================================================================
    task->SelectCollisionCandidates(AliVEvent::kAnyINT);
    task->SetESDtrackCutsM(AlidNdPtTools::CreateESDtrackCuts("defaultEta08"));
    //default, nogeo, nogeonogold, geoCut116, geoCut117, geoCut118 geoCut119, tpconlyminimal
    task->SetESDtrackCuts(0,AlidNdPtTools::CreateESDtrackCuts("defaultEta08"));
    task->SetESDtrackCuts(1,AlidNdPtTools::CreateESDtrackCuts("tpcitsnogeoEta08"));
    task->SetESDtrackCuts(2,AlidNdPtTools::CreateESDtrackCuts("tpcitsnogeonogoldEta08"));
    task->SetESDtrackCuts(3,AlidNdPtTools::CreateESDtrackCuts("defaultEta08", 116));
    task->SetESDtrackCuts(4,AlidNdPtTools::CreateESDtrackCuts("defaultEta08", 117));
    task->SetESDtrackCuts(5,AlidNdPtTools::CreateESDtrackCuts("defaultEta08", 118));
    task->SetESDtrackCuts(6,AlidNdPtTools::CreateESDtrackCuts("defaultEta08", 119));
    task->SetESDtrackCuts(7,AlidNdPtTools::CreateESDtrackCuts("tpconlyminimalEta08"));
    task->SetNeedEventMult(kTRUE);
    task->SetNeedTrackIP(kTRUE);
    
    // attach the task to the manager and configure in and ouput
    //===========================================================================
    mgr->AddTask(task);
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task,1,mgr->CreateContainer(name, TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    
  return task;
}

