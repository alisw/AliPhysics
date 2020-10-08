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
#include "AlidNdPtEventCuts.h"
#include "AlidNdPt.h"
#include "AliAnalysisTaskMKBase.h"
#include "AliAnalysisTaskSpectraINEL0.h"
class AliAnalysisTaskSpectraINEL0;

using namespace std;

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSpectraINEL0)
/// \endcond
//_____________________________________________________________________________

AliAnalysisTaskSpectraINEL0::AliAnalysisTaskSpectraINEL0()
    : AliAnalysisTaskMKBase()
    , fHistEffCont(0)
    , fHistTrack(0)
    , fHistEvent(0)
    , fHistVtxInfo(0)
    , fHistTrackINEL0(0)
{
    // default contructor
}

//_____________________________________________________________________________

AliAnalysisTaskSpectraINEL0::AliAnalysisTaskSpectraINEL0(const char* name)
    : AliAnalysisTaskMKBase(name)
    , fHistEffCont(0)
    , fHistTrack(0)
    , fHistEvent(0)
    , fHistVtxInfo(0)
    , fHistTrackINEL0(0)
{
    // constructor
}

//_____________________________________________________________________________

AliAnalysisTaskSpectraINEL0::~AliAnalysisTaskSpectraINEL0()
{
    // destructor
}

//_____________________________________________________________________________

void AliAnalysisTaskSpectraINEL0::AddOutput()
{
    AddAxis("cent");
    AddAxis("nAcc","mult6kcoarse");
    AddAxis("MCpT","pt");
    AddAxis("MCQ",3,-1.5,1.5);
    AddAxis("MCpid",10,-0.5,9.5);  // 0=e, 1=mu, 2=pi, 3=K, 4=p, 6=sigmaP, 7=sigmaM, 8=xi, 9=omega, 5=other
    AddAxis("MCinfo",4,-0.5,3.5);  // 0=prim, 1=decay 2=material, 3=genprim
    fHistEffCont = CreateHist("fHistEffCont");
    fOutputList->Add(fHistEffCont);

    AddAxis("cent");
    AddAxis("nAcc","mult6kcoarse");
    AddAxis("pt");
    AddAxis("Q",3,-1.5,1.5);
    fHistTrack = CreateHist("fHistTrack");
    fOutputList->Add(fHistTrack);

    AddAxis("cent");
    AddAxis("nAcc","mult6kfine");
    AddAxis("zV",8,-20,20);
    fHistEvent = CreateHist("fHistEvent");
    fOutputList->Add(fHistEvent);

    AddAxis("VtxInfo", 2, -0.5, 1.5); // 0 == with Vertex; 1== without Vertex
    AddAxis("nAcc","mult6kfine");
    AddAxis("zV",8,-30,30);
    fHistVtxInfo = CreateHist("fHistVtxInfo");
    fOutputList->Add(fHistVtxInfo);

    AddAxis("cent");
    AddAxis("nAcc","mult6kcoarse");
    AddAxis("MCpT","pt");
    AddAxis("EventSelection", 2, -0.5, 1.5); // 0 == No Selection; 1 == with selection
    fHistTrackINEL0 = CreateHist("fHistTrackINEL0");
    fOutputList->Add(fHistTrackINEL0);

}

//_____________________________________________________________________________

Bool_t AliAnalysisTaskSpectraINEL0::IsEventSelected()
{
  return fIsAcceptedAliEventCuts;
  //return oldEventCuts();
}

Bool_t AliAnalysisTaskSpectraINEL0::oldEventCuts()
{
  return kFALSE;
}

//_____________________________________________________________________________

void AliAnalysisTaskSpectraINEL0::FillDefaultHistograms(Int_t step) {
    // call function from mother class
    AliAnalysisTaskMKBase::FillDefaultHistograms(step);

    if(step==0){
        if(fVtxStatus) FillHist(fHistVtxInfo, 0, fNTracksAcc, fZv);
        else FillHist(fHistVtxInfo, 1, fNTracksAcc, fZv);
    }
    if(step==1){
        LoopOverAllTracks(step);
        FillHist(fHistEvent, fMultPercentileV0M, fNTracksAcc, fZv);
    }


    if (fIsMC) LoopOverAllParticles(step);
}

//_____________________________________________________________________________

void AliAnalysisTaskSpectraINEL0::AnaTrack(Int_t flag)
{
    if(0==flag) return; // data only after event selection
    if (!fAcceptTrackM) return;

    FillHist(fHistTrack, fMultPercentileV0M, fNTracksAcc, fPt, fChargeSign);
}

//_____________________________________________________________________________

void AliAnalysisTaskSpectraINEL0::AnaTrackMC(Int_t flag)
{
    if (!fAcceptTrackM) return;

    if (fMCParticleType==AlidNdPtTools::kOther) { Log("RecTrack.PDG.",fMCPDGCode); }
    if (TMath::Abs(fMCQ > 1)) { Log("RecTrack.Q>1.PDG.",fMCPDGCode); }

    if (flag == 0) {

    } else if (flag == 1) {

        FillHist(fHistEffCont, fMultPercentileV0M, fNTracksAcc, fMCPt, fMCChargeSign, fMCParticleType, fMCProdcutionType);
    } else {
        Err("AliAnalysisTaskMKBase::FillDefaultHistograms:InvalidStep");
    }
}

//_____________________________________________________________________________

void AliAnalysisTaskSpectraINEL0::AnaParticleMC(Int_t flag)
{
    if (!fMCisPrim) return;
    if (!fMCIsCharged) return;
    if (TMath::Abs(fMCEta) > 0.8) return;

    if (fMCParticleType==AlidNdPtTools::kOther) { Log("GenPrim.PDG.",fMCPDGCode); }
    if (TMath::Abs(fMCQ > 1)) { Log("GenPrim.Q>1.PDG.",fMCPDGCode); }

    if (flag == 0) {
        if(fZv < 10. && fZv > -10.) FillHist(fHistTrackINEL0, fMultPercentileV0M, fNTracksAcc, fMCPt, 0);
    } else if (flag == 1) {
        FillHist(fHistTrackINEL0, fMultPercentileV0M, fNTracksAcc, fMCPt, 1);

        FillHist(fHistEffCont, fMultPercentileV0M, fNTracksAcc, fMCPt, fMCChargeSign, fMCParticleType, 3);
    } else {
        Err("AliAnalysisTaskMKBase::FillDefaultHistograms:InvalidStep");
    }
}

//_____________________________________________________________________________

AliAnalysisTaskSpectraINEL0* AliAnalysisTaskSpectraINEL0::AddTaskSpectraINEL0(const char* name, const char* outfile)
{
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskSpectraINEL0", "No analysis manager to connect to.");
        return 0;
    }

    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddTaskSpectraINEL0", "This task requires an input event handler");
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
    AliAnalysisTaskSpectraINEL0 *task = new AliAnalysisTaskSpectraINEL0(name);
    if (!task) { return 0; }

    // configure the task
    //===========================================================================
    task->SelectCollisionCandidates(AliVEvent::kAnyINT);
    task->SetESDtrackCutsM(AlidNdPtTools::CreateESDtrackCuts("defaultEta08"));
//     task->SetESDtrackCuts(0,AlidNdPtTools::CreateESDtrackCuts("defaultEta08"));

    // attach the task to the manager and configure in and ouput
    //===========================================================================
    mgr->AddTask(task);
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task,1,mgr->CreateContainer(name, TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));

  return task;
}
