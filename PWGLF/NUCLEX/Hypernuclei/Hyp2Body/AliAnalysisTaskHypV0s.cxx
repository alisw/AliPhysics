#include "AliAnalysisTaskHypV0s.h"

#include <Riostream.h>
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliVertexerHyperTriton2Body.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskHypV0s);

AliAnalysisTaskHypV0s::AliAnalysisTaskHypV0s(std::string name)
    : AliAnalysisTaskSE(name.data()),
      fV0Vertexer{},
      fInputHandler{nullptr},
      fPIDResponse{nullptr}

{
}

AliAnalysisTaskHypV0s::~AliAnalysisTaskHypV0s()
{
}

void AliAnalysisTaskHypV0s::UserCreateOutputObjects()
{
    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    fInputHandler = (AliInputEventHandler *)(man->GetInputEventHandler());
    fPIDResponse = fInputHandler->GetPIDResponse();
    fV0Vertexer. SetV0VertexerDCAFirstToPV(0.);
    fV0Vertexer. SetV0VertexerDCASecondToPV(0.);
    fV0Vertexer. SetV0VertexerDCAV0Daughters(2);
    fV0Vertexer.SetV0VertexerCosinePA(-2);
    fV0Vertexer.SetV0VertexerMinRadius(0.);
    fV0Vertexer.SetV0VertexerMaxRadius(200);
    fV0Vertexer.SetV0VertexerMaxChisquare(1000);
}
void AliAnalysisTaskHypV0s::UserExec(Option_t *)
{

    AliESDEvent *esdEvent = dynamic_cast<AliESDEvent *>(InputEvent());
    AliMCEvent *mcEvent = MCEvent();
    std::vector<AliESDv0> V0Vector;
    V0Vector = fV0Vertexer.Tracks2V0vertices(esdEvent, fPIDResponse, mcEvent);
    esdEvent->ResetV0s();
    int nV0s = V0Vector.size();

    for (int iV0 = 0; iV0 < nV0s; iV0++)
    { // This is the begining of the V0 loop (we analyse only offline
        // V0s)
        AliESDv0 *v0 = &V0Vector[iV0];
        esdEvent->AddV0(v0);
    }
    // loop on tracklets to match them with mother MClabel1m
}

void AliAnalysisTaskHypV0s::Terminate(Option_t *) {}

AliAnalysisTaskHypV0s *AliAnalysisTaskHypV0s::AddTask(TString suffix)
{
    // Get the current analysis manager
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr)
    {
        ::Error("AddTaskHypV0s", "No analysis manager found.");
        return nullptr;
    }
    // Check the analysis type using the event handlers connected to the analysis
    // manager.

    if (!mgr->GetInputEventHandler())
    {
        ::Error("AddTaskHypV0s", "This task requires an input event handler");
        return nullptr;
    }

    TString tskname = "AliAnalysisTaskHypV0s";
    tskname.Append(suffix.Data());
    AliAnalysisTaskHypV0s *task = new AliAnalysisTaskHypV0s(tskname.Data());

    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

    return task;
}
