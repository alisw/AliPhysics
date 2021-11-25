#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TList.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency.h"
#include "AliAODMCParticle.h"
#include "AliMultSelection.h"

using namespace std;

ClassImp(AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency)

    AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency::AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency() : AliAnalysisTaskSE(),
                                                     fAOD(0), fOutputList(0)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency::AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency(const char *name) : AliAnalysisTaskSE(name),
                                                                 fAOD(0), fOutputList(0)
{
    // constructor
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency::~AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency()
{
    // destructor
    if (fOutputList)
    {
        delete fOutputList;
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency::UserCreateOutputObjects()
{
    fOutputList = new TList();
    fOutputList->SetOwner(kTRUE);

    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    if (man)
    {
        AliInputEventHandler *inputHandler = dynamic_cast<AliInputEventHandler *>(man->GetInputEventHandler());
        if (inputHandler)
            fPIDResponse = inputHandler->GetPIDResponse();
        else
            AliFatal("Input handler needed");
    }
    PostData(1, fOutputList);
}
//_____________________________________________________________________________
void AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency::UserExec(Option_t *)
{
    fAOD = dynamic_cast<AliAODEvent *>(InputEvent());

    if (!fAOD)
    {
        PostData(1, fOutputList);
        return;
    }
    PostData(1, fOutputList);
}
//_____________________________________________________________________________
void AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________
