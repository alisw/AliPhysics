#ifndef __CINT__//|
#include <AliAnalysisManager.h>//|
#include <AliMultiInputEventHandler.h>//|
#include <EventMixing/EventMixing/AliMixEventPool.h>//|
#include <EventMixing/EventMixing/AliMixEventCutObj.h>//|
#include <EventMixing/EventMixing/AliMixInputEventHandler.h>//|
#include <AliVEvent.h>//|
#endif//|

void AddMixingHandler(Double_t centMin = 70, Double_t centMax = 80, Double_t centStep = 2, AliMultiInputEventHandler* multiInputHandler, Bool_t useMC = kFALSE, Bool_t usePhysSel = kFALSE,TString opts = "")
{

   if (!multiInputHandler) return;

   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   const Int_t bufferSize = 1;
   const Int_t mixNum = 5;
   AliMixInputEventHandler *mixHandler = new AliMixInputEventHandler(bufferSize, mixNum);
   mixHandler->SetInputHandlerForMixing(dynamic_cast<AliMultiInputEventHandler*>(mgr->GetInputEventHandler()));
   AliMixEventPool *evPool = new AliMixEventPool();

   //AliMixEventCutObj *multi = new AliMixEventCutObj(AliMixEventCutObj::kMultiplicity, 2, 10002, 10000);
   AliMixEventCutObj *zvertex = new AliMixEventCutObj(AliMixEventCutObj::kZVertex, -10, 10, 5);

   AliMixEventCutObj *centrality = new AliMixEventCutObj(AliMixEventCutObj::kCentrality, centMin, centMax, centStep, "V0M");

   evPool->AddCut(centrality);
   //evPool->AddCut(multi);
   evPool->AddCut(zvertex);

   // adds event pool (comment it and u will have default mixing)
   mixHandler->SetEventPool(evPool);

   // only use events with physics selection
   if (usePhysSel) mixHandler->SelectCollisionCandidates(AliVEvent::kMB);

   multiInputHandler->AddInputEventHandler(mixHandler);

}
