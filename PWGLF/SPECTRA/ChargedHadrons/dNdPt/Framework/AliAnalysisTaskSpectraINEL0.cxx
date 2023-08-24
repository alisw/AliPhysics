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
#include "AliAnalysisTaskSpectraINEL0.h"
class AliAnalysisTaskSpectraINEL0;

// namespace

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSpectraINEL0)
/// \endcond
//_____________________________________________________________________________

AliAnalysisTaskSpectraINEL0::AliAnalysisTaskSpectraINEL0()
    : AliAnalysisTaskMKBase()
{
    // default contructor
    AliAnalysisTaskMKBase::fUseBaseOutput = kTRUE;
    AliAnalysisTaskMKBase::fNeedEventVertex = kTRUE;
    AliAnalysisTaskMKBase::fNeedEventMult = kTRUE;
}

//_____________________________________________________________________________

AliAnalysisTaskSpectraINEL0::AliAnalysisTaskSpectraINEL0(const char* name)
    : AliAnalysisTaskMKBase(name)
{
    // constructor
    AliAnalysisTaskMKBase::fUseBaseOutput = kTRUE;
    AliAnalysisTaskMKBase::fNeedEventVertex = kTRUE;
    AliAnalysisTaskMKBase::fNeedEventMult = kTRUE;
}

//_____________________________________________________________________________

AliAnalysisTaskSpectraINEL0::~AliAnalysisTaskSpectraINEL0()
{
    // destructor
}

//_____________________________________________________________________________

void AliAnalysisTaskSpectraINEL0::AddOutput()
{
    std::vector<double> ptbins = {0.15,  0.2,   0.25,  0.3,   0.35,  0.4,   0.45, 0.5,  0.55,
                                0.6,  0.65,  0.7,   0.75,  0.8,   0.85,  0.9,   0.95,  1.0,   1.1,  1.2,  1.3,
                                1.4,  1.5,   1.6,   1.7,   1.8,   1.9,   2.0,   2.2,   2.4,   2.6,  2.8,  3.0,
                                3.2,  3.4,   3.6,   3.8,   4.0,   4.5,   5.0,   5.5,   6.0,   6.5,  7.0,  8.0,
                                9.0,  10.0,  11.0,  12.0,  13.0,  14.0,  15.0,  16.0,  18.0,  20.0, 22.0, 24.0,
                                26.0, 28.0,  30.0,  32.0,  34.0,  36.0,  40.0,  45.0,  50.0,  60.0, 70.0, 80.0,
                                90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 180.0, 200.0 };
    std::vector<double> inverseptbins;
    for(int i = ptbins.size() - 1; i > 0 ; i--){
      inverseptbins.push_back(1./ptbins[i]);
    }

    std::vector<double> centBins = {0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100. };

    const Int_t nmultbins = 241;
    std::vector<double> multBins;
    multBins.reserve(nmultbins + 1);
    multBins.push_back(-0.5);
    {
        int i = 0;
        for (; i <= 100; i++) {
            multBins.push_back(multBins.back() + 1.);
        }
        for (; i <= 100 + 90; i++) {
            multBins.push_back(multBins.back() + 10.);
        }
        for (; i <= 10 + 90 + 50; i++) {
            multBins.push_back(multBins.back() + 100.);
        }
    }

    Axis centAxis = {"cent", "centrality", centBins};
    Axis multAxis = {"mult", "#it{N}_{ch}", multBins};
    Axis zVrtAxis = {"zV", "vertex_{Z}", {-20, 20}, 8};
    Axis ptAxis = {"pt", "#it{p}_{T} (GeV/c)", ptbins};
    Axis QAxis = {"Q", "Q", {-1.5,1.5}, 3};
    Axis relResoAxis = {"rel. Reso", "#sigma(#it{p}_{T}) / #it{p}_{T}", {0.0,0.3}, 300};
    Axis Sigma1ptAxis = {"Sigma1ptAxis", "#sigma(1/#it{p}_{T})", {0.0,0.3}, 300};
    Axis inverseptAxis = {"1/pt", "1/#it{p}_{T} (GeV/c)", {0.005, 6.677}, 60};
    //Axis MCpidAxis = {"MCpid", "MCpid", {-0.5, 9.5}, 10}; // 0=e, 1=mu, 2=pi, 3=K, 4=p, 6=sigmaP, 7=sigmaM, 8=xi, 9=omega, 5=other
    Axis MCinfoAxis = {"MCinfo", "MCinfo", {-0.5,3.5}, 4}; // 0=prim, 1=decay 2=material, 3=genprim
    //Axis MCQAxis = {"MCQ", "MCQ", {-1.5,1.5}, 3};
    Axis MCptAxis = {"MCpt", "#it{p}_{T} (GeV/c)", ptbins};
    Axis EvtSelectionAxis = {"EvtSele", "EvtSele", {-0.5,1.5}, 2}; // 0 == No Selection; 1 == with selection
    Axis VtxInfoAxis = {"VtxInfo", "VtxInfo", {-0.5,1.5}, 2}; // 0 == with Vertex; 1== without Vertex


    fHistEffCont.AddAxis(centAxis);
    fHistEffCont.AddAxis(multAxis);
    fHistEffCont.AddAxis(MCptAxis);
    fHistEffCont.AddAxis(MCinfoAxis); // 0=prim, 1=decay 2=material, 3=genprim
    fOutputList->Add(fHistEffCont.GenerateHist("fHistEffCont"));

    fHistTrack.AddAxis(centAxis);
    fHistTrack.AddAxis(multAxis);
    fHistTrack.AddAxis(ptAxis);
    fHistTrack.AddAxis(QAxis);
    fOutputList->Add(fHistTrack.GenerateHist("fHistTrack"));

    fHistTrackINEL0.AddAxis(centAxis);
    fHistTrackINEL0.AddAxis(multAxis);
    fHistTrackINEL0.AddAxis(MCptAxis);
    fHistTrackINEL0.AddAxis(EvtSelectionAxis);
    fOutputList->Add(fHistTrackINEL0.GenerateHist("fHistTrackINEL0"));

    fHistRelResoFromCov.AddAxis(ptAxis);
    fHistRelResoFromCov.AddAxis(relResoAxis);
    fOutputList->Add(fHistRelResoFromCov.GenerateHist("fHistRelResoFromCov"));

    fHistSigma1pt.AddAxis(inverseptAxis);
    fHistSigma1pt.AddAxis(Sigma1ptAxis);
    fOutputList->Add(fHistSigma1pt.GenerateHist("fHistSigma1pt"));

    fHistEvent.AddAxis(centAxis);
    fHistEvent.AddAxis(multAxis);
    fHistEvent.AddAxis(zVrtAxis);
    fOutputList->Add(fHistEvent.GenerateHist("fHistEvent"));

    fHistVtxInfo.AddAxis(VtxInfoAxis);// 0 == with Vertex; 1== without Vertex
    fHistVtxInfo.AddAxis(multAxis);
    fHistVtxInfo.AddAxis(zVrtAxis);
    fOutputList->Add(fHistVtxInfo.GenerateHist("fHistVtxInfo"));


}

//_____________________________________________________________________________

Bool_t AliAnalysisTaskSpectraINEL0::IsEventSelected()
{
  return fIsAcceptedAliEventCuts;
}

//_____________________________________________________________________________

void AliAnalysisTaskSpectraINEL0::FillDefaultHistograms(Int_t step) {
    // call function from mother class
    AliAnalysisTaskMKBase::FillDefaultHistograms(step);

    if(step==0){
        if(fVtxStatus) fHistVtxInfo.Fill(0, fNTracksAcc, fZv);
        else fHistVtxInfo.Fill(1, fNTracksAcc, fZv);
    }
    if(step==1){
        LoopOverAllTracks(step);
        fHistEvent.Fill(fMultPercentileV0M, fNTracksAcc, fZv);
    }


    if (fIsMC) LoopOverAllParticles(step);
}

//_____________________________________________________________________________

void AliAnalysisTaskSpectraINEL0::AnaTrack(Int_t flag)
{
    if(0==flag) return; // data only after event selection
    if (!fAcceptTrackM) return;

    fHistTrack.Fill(fMultPercentileV0M, fNTracksAcc, fPt, fChargeSign);
    fHistRelResoFromCov.Fill(fPt, fSigma1Pt*fPt);
    fHistSigma1pt.Fill(f1Pt, fSigma1Pt);
}

//_____________________________________________________________________________

void AliAnalysisTaskSpectraINEL0::AnaTrackMC(Int_t flag)
{
    if(fMCPileUpTrack) return; // reject pileup tracks
    if (!fAcceptTrackM) return;

    if (fMCParticleType==AlidNdPtTools::kOther) { Log("RecTrack.PDG.",fMCPDGCode); }
    if (TMath::Abs(fMCQ > 1)) { Log("RecTrack.Q>1.PDG.",fMCPDGCode); }

    if (flag == 0) {

    } else if (flag == 1) {

        fHistEffCont.Fill(fMultPercentileV0M, fNTracksAcc, fMCPt, fMCProdcutionType);
    } else {
        Err("AliAnalysisTaskMKBase::FillDefaultHistograms:InvalidStep");
    }
}

//_____________________________________________________________________________

void AliAnalysisTaskSpectraINEL0::AnaParticleMC(Int_t flag)
{
    if(fMCPileUpTrack) return; // reject pileup tracks

    if (!fMCisPrim) return;
    if (!fMCIsCharged) return;
    if (TMath::Abs(fMCEta) > 0.8) return;

    if (fMCParticleType==AlidNdPtTools::kOther) { Log("GenPrim.PDG.",fMCPDGCode); }
    if (TMath::Abs(fMCQ > 1)) { Log("GenPrim.Q>1.PDG.",fMCPDGCode); }

    if (flag == 0) {
        if(fZv < 10. && fZv > -10.) fHistTrackINEL0.Fill(fMultPercentileV0M, fNTracksAcc, fMCPt, 0);
    } else if (flag == 1) {
        fHistTrackINEL0.Fill(fMultPercentileV0M, fNTracksAcc, fMCPt, 1);

        fHistEffCont.Fill(fMultPercentileV0M, fNTracksAcc, fMCPt, 3);
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
