//
// Class AliRsnAnalysisTask2ndStep
//
// TODO
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//
#include <TTree.h>
#include <TChain.h>
#include <TH1F.h>
#include <TList.h>
#include <TObjArray.h>

#include "AliLog.h"

#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"

#include "AliRsnEvent.h"
#include "AliRsnPair.h"
#include "AliRsnPairMgr.h"

#include "AliRsnAnalysisTask2ndStep.h"

ClassImp(AliRsnAnalysisTask2ndStep)

//_____________________________________________________________________________
AliRsnAnalysisTask2ndStep::AliRsnAnalysisTask2ndStep(const char *name) :
  AliAnalysisTaskSE(name),
  fOutList(0x0),
  fPairMgrs(0x0),
  fEventBuffer(0x0),
  fRsnHandlerAOD(0x0),
  fAnalysisMgr(0x0)
{
//
// Default constructor
//

    DefineOutput(1, TList::Class());
}

//_____________________________________________________________________________
void AliRsnAnalysisTask2ndStep::CreateHandlers(AliAnalysisManager* am)
{
//
// Sets the handlers.
//

    fAnalysisMgr = am;
    if (!fAnalysisMgr) {
        AliWarning("Passed a NULL AnalysisManager!");
        return;
    }
    
    // this task is made only for building histograms from
    // a set of RSN events saved as non-standard branches in AOD tree
    fRsnHandlerAOD = new AliAODInputHandler();
    if (fRsnHandlerAOD) fAnalysisMgr->SetInputEventHandler(fRsnHandlerAOD);
}

//_____________________________________________________________________________
void AliRsnAnalysisTask2ndStep::UserCreateOutputObjects()
{
//
// Creates output objects, as a TList of all histograms
//
    AliInfo("--->");
    
    //OpenFile (1);
    if (fOutList) {
        fOutList->Clear();
        delete fOutList;
    }
    fOutList = new TList();
    fOutList->SetOwner();
    AliRsnPairMgr *pairMgr = 0x0;
    AliRsnPair *pair = 0x0;
    
    if (!fPairMgrs) {
        AliError("No pair managers defined!");
        return;
    }
    
    TObjArrayIter next(fPairMgrs);
    while ( (pairMgr = (AliRsnPairMgr*)next()) ) {
        TList *listTmp = new TList();
        listTmp->SetName(pairMgr->GetName());
        TObjArrayIter nextPair(pairMgr->GetPairs());
        while ( (pair = (AliRsnPair*)nextPair()) ) {
            pair->Init();
            //pair->Print();
            pair->GenerateHistograms(pairMgr->GetName(), listTmp);
            //listTmp->Add(pair->GenerateHistograms(pairMgr->GetName()));
            //fOutList->Add(listTmp);
        }
        fOutList->Add(listTmp);
    }
    fEventBuffer = new AliRsnEventBuffer(1000);
    
    AliInfo("<---");
}

//_____________________________________________________________________________
void AliRsnAnalysisTask2ndStep::UserExec(Option_t *)
{
//
// Executes the task
//
    
    if (fEntry++ % 1000 == 0) AliInfo(Form("Event %d",-1));
    
    // retrieve AOD
    AliAODEvent *aod = dynamic_cast<AliAODEvent*>(fInputEvent);
    
    // find RSN event
    AliRsnEvent *rsn = (AliRsnEvent*)(aod->GetList()->FindObject("rsnEvents"));
    
    // execute analysis
    ProcessEventAnalysis(rsn);
    PostEventProcess();
    
    PostData(1, fOutList);
}

//_____________________________________________________________________________
void AliRsnAnalysisTask2ndStep::Terminate(Option_t *)
{
//
// A simple check on the output list
//

    fOutList = dynamic_cast<TList*>(GetOutputData(1));
    if (!fOutList) { AliError("--- Output list not available ---"); return; }
    fOutList->Print();
}

//_____________________________________________________________________________
void AliRsnAnalysisTask2ndStep::ProcessEventAnalysis(AliRsnEvent *curEvent)
{
//
// Process of one event
//

    // Adds event to Event Buffer
    fEventBuffer->AddEvent(curEvent);
    
    // loop over all Pair managers
    AliRsnPair *pair = 0;
    AliRsnPairMgr *mgr = 0;
    TObjArrayIter nextMgr(fPairMgrs);
    while ( (mgr = (AliRsnPairMgr*)nextMgr()) ) {
        TObjArrayIter nextPair(mgr->GetPairs());
        while ( (pair = (AliRsnPair*)nextPair()) ) {
            pair->ProcessPair(fEventBuffer);
        }
    }
}

//_____________________________________________________________________________
void AliRsnAnalysisTask2ndStep::PostEventProcess()
{
//
// Some work done after event processing
//

}

//_____________________________________________________________________________
void AliRsnAnalysisTask2ndStep::AddPairMgr(AliRsnPairMgr *pairmgr)
{
    if (!fPairMgrs) fPairMgrs = new TObjArray;
    fPairMgrs->Add(pairmgr);
}
