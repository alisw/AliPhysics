#include "AliAnalysisTaskDCArStudy.h"

#include <iostream>
#include "TChain.h"
#include "TH1F.h"
#include "TList.h"
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
#include "AliMCSpectraWeights.h"


class AliAnalysisTaskDCArStudy;

namespace {
using namespace Hist;
using namespace std;
} // namespace

ClassImp(AliAnalysisTaskDCArStudy)

//_____________________________________________________________________________

AliAnalysisTaskDCArStudy::AliAnalysisTaskDCArStudy()
    : AliAnalysisTaskMKBase()
    , fHistDCA{}
    , fHistDCAPCC{}
    , fMCSpectraWeights(nullptr)
{
    // default contructor

}

//_____________________________________________________________________________

AliAnalysisTaskDCArStudy::AliAnalysisTaskDCArStudy(const char* name)
    : AliAnalysisTaskMKBase(name)
    , fHistDCA{}
    , fHistDCAPCC{}
    , fMCSpectraWeights(nullptr)
{
    // constructor

}

//_____________________________________________________________________________

AliAnalysisTaskDCArStudy::~AliAnalysisTaskDCArStudy()
{
    // destructor
}

//_____________________________________________________________________________

void AliAnalysisTaskDCArStudy::AddOutput()
{
    //dcar:pt:mult:mcinfo
    auto const DCAbins = 500;
    auto const DCAbinWidth = 2./DCAbins;
    std::vector<double> centBins = {0.,  10., 20., 30., 40.,
        50., 60., 70., 80., 90.};
    std::vector<double> ptBins = {0.0, 0.1,0.12,0.14,0.16,0.18,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.2,2.4,2.6,2.8,3,3.2,3.4,3.6,3.8,4,4.5,5,5.5,6,6.5,7,8,10,13,20, 30, 50, 80, 100, 200};
    std::vector<double> multBins;
    {
        // variable mult binning total 0-6000
        // 1-width 0-10              11
        // 10-width 10-100            9
        // 100-width 100-1000         9
        // 200-width 1000-6000       25
        multBins.push_back(-0.5);
        int i = 0;
        for(; i <= 10; i++)
        {
            multBins.push_back(multBins.back() + 1);
        }
        for(; i <= 9; i++)
        {
            multBins.push_back(multBins.back() + 10);
        }
        for(; i <= 9; i++)
        {
            multBins.push_back(multBins.back() + 100);
        }
        for(; i <= 25; i++)
        {
            multBins.push_back(multBins.back() + 200);
        }
    }

    Axis DCAaxis{"DCA_{xy}", "DCAxy",  {-1-DCAbinWidth/2.,1+DCAbinWidth/2}, DCAbins+1};
    Axis multAxisNch{"Nch", "#it{N}_{ch}", multBins};
    Axis centAxis{"cent", "centrality", centBins};
    Axis ptAxis{"pt", "#it{p}_{T} (GeV/c)", ptBins};
    // 0=prim, 1=decay 2=material, 3=genprim
    Axis mcInfoAxis{"mcInfo", "mcInfo", {-0.5, 3.5}, 4};
    // 0=none, 1=weighted 2=weightedRandom, 3=weightSys,
    // 4=weightSysRandom
    //    Axis mcWeightAxis{"weight", "weight", {-0.5, 4.5}, 5};
    Axis weightValues{"weight", "weight", {0.,5.}, 50};
    Axis secLambdaKaon{"SecDecayType", "SecDecayType", {-1.5, 2.5}, 4};//0=lambda, 1=kaon, 2=other

    // ------
    // hists
    // -----
    double requiredMemory = 0.;

    fHistDCA.AddAxis(DCAaxis);
    fHistDCA.AddAxis(ptAxis);
    fHistDCA.AddAxis(multAxisNch);
    fHistDCA.AddAxis(centAxis);
    fHistDCA.AddAxis(mcInfoAxis);
    fOutputList->Add(fHistDCA.GenerateHist("fHistDCA"));
    requiredMemory += fHistDCA.GetSize();

    fHistDCAPCC.AddAxis(DCAaxis);
    fHistDCAPCC.AddAxis(ptAxis);
    fHistDCAPCC.AddAxis(multAxisNch);
    fHistDCAPCC.AddAxis(centAxis);
    fHistDCAPCC.AddAxis(mcInfoAxis);
    fOutputList->Add(fHistDCAPCC.GenerateHist("fHistDCAPCC"));
    requiredMemory += fHistDCAPCC.GetSize();


    fHistSecWeights.AddAxis(ptAxis);
    fHistSecWeights.AddAxis(weightValues);
    fHistSecWeights.AddAxis(secLambdaKaon);
    fOutputList->Add(fHistSecWeights.GenerateHist("fHistSecWeights"));
    requiredMemory += fHistSecWeights.GetSize();

    AliError(Form("Estimated memory usage of histograms: %.0f Bytes (%f MiB)",
                  requiredMemory, requiredMemory / 1048576));

}

//_____________________________________________________________________________

Bool_t AliAnalysisTaskDCArStudy::IsEventSelected()
{
    return fIsAcceptedAliEventCuts;
}

//_____________________________________________________________________________


void AliAnalysisTaskDCArStudy::AnaEvent()
{
    LoopOverAllTracks();
}

void AliAnalysisTaskDCArStudy::AnaEventMC() {

    AliMCSpectraWeightsHandler* mcWeightsHandler = static_cast<AliMCSpectraWeightsHandler*>(fEvent->FindListObject("fMCSpectraWeights"));
    fMCSpectraWeights = (mcWeightsHandler) ? mcWeightsHandler->fMCSpectraWeight : nullptr;
    
//    LoopOverAllParticles();
    LoopOverAllTracks();
}

//_____________________________________________________________________________

void AliAnalysisTaskDCArStudy::AnaTrackMC(Int_t flag)
{
    if (fAcceptTrack[0]) {
        double fMCweight = 1.0;
        if(fMCSpectraWeights && 0==fMCPrimSec && !fMCPileUpTrack && fMCParticle->Particle()){ // only for primary particles
            fMCweight = fMCSpectraWeights->GetMCSpectraWeight(fMCParticle->Particle(), 0);
        }
        if(fMCSpectraWeights && 1==fMCPrimSec && !fMCPileUpTrack && fMCParticle->Particle()){ // only for secondaries from decay
            fMCweight = fMCSpectraWeights->GetWeightForSecondaryParticle(fMCParticle->Particle());

            N.Fill(fPt, fMCweight, fMCSpectraWeights->IdentifySecondaryType(fMCParticle->Particle()));

        }
        fHistDCA.Fill(fDCAr, fPt, fNTracksAcc, fMultPercentileV0M, fMCPrimSec);
        fHistDCAPCC.FillWeight(fMCweight, fDCAr, fPt, fNTracksAcc, fMultPercentileV0M, fMCPrimSec);



    }
}

//_____________________________________________________________________________

void AliAnalysisTaskDCArStudy::AnaTrackDATA(Int_t flag)
{
    if (fAcceptTrack[0]) {
        fHistDCA.Fill(fDCAr, fPt, fNTracksAcc, fMultPercentileV0M, -1);
        fHistDCAPCC.Fill(fDCAr, fPt, fNTracksAcc, fMultPercentileV0M, -1);
    }
}


//_____________________________________________________________________________

AliAnalysisTaskDCArStudy* AliAnalysisTaskDCArStudy::AddTaskDCArStudy(const char* name, const char* outfile)
{
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskDCArStudy", "No analysis manager to connect to.");
        return 0;
    }

    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddTaskDCArStudy", "This task requires an input event handler");
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
    AliAnalysisTaskDCArStudy *task = new AliAnalysisTaskDCArStudy(name);
    if (!task) { return 0; }

    // configure the task
    //===========================================================================
    task->SelectCollisionCandidates(AliVEvent::kAnyINT);
    task->SetESDtrackCutsM(AlidNdPtTools::CreateESDtrackCuts("defaultEta08"));
    task->SetESDtrackCuts(0,AlidNdPtTools::CreateESDtrackCuts("TPCITSforDCArStudyEta08"));
    task->SetNeedEventMult(kTRUE);
    task->SetNeedEventVertex(kTRUE);
    task->SetNeedTrackTPC(kTRUE);

    // attach the task to the manager and configure in and ouput
    //===========================================================================
    mgr->AddTask(task);
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task,1,mgr->CreateContainer(name, TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));

  return task;
}
