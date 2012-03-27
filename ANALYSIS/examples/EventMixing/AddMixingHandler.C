#ifndef __CINT__
#include <AliVEvent.h>
#endif
void AddMixingHandler ( AliMultiInputEventHandler *multiInputHandler,TString format = "esd", Bool_t useMC = kFALSE, TString opts = "" ) {

   const Int_t bufferSize = 1;
   const Int_t mixNum = 5;
   if ( !multiInputHandler ) return;

   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   AliMixInputEventHandler *mixHandler = new AliMixInputEventHandler ( bufferSize, mixNum );
   mixHandler->SetInputHandlerForMixing ( dynamic_cast<AliMultiInputEventHandler *> ( mgr->GetInputEventHandler() ) );
   AliMixEventPool *evPool = new AliMixEventPool();

   AliMixEventCutObj *multi = new AliMixEventCutObj ( AliMixEventCutObj::kMultiplicity, 2, 102, 10 );
   AliMixEventCutObj *zvertex = new AliMixEventCutObj ( AliMixEventCutObj::kZVertex, -10, 10, 1 );

//    AliMixEventCutObj *centrality = new AliMixEventCutObj(AliMixEventCutObj::kCentrality, 0, 20, 10, "V0M");
//    AliMixEventCutObj *multi = new AliMixEventCutObj(AliMixEventCutObj::kMultiplicity, 2, 102, 10);
//    AliMixEventCutObj *zvertex = new AliMixEventCutObj(AliMixEventCutObj::kZVertex, -5, 5, 1);
//
//    AliMixEventCutObj *multi = new AliMixEventCutObj(AliMixEventCutObj::kMultiplicity, 2, 5002, 100);
//    AliMixEventCutObj *zvertex = new AliMixEventCutObj(AliMixEventCutObj::kZVertex, -10, 10, 1);
//
//     AliMixEventCutObj *centrality = new AliMixEventCutObj ( AliMixEventCutObj::kCentrality, 0, 100, 10, "V0M" );
//     AliMixEventCutObj *zvertex = new AliMixEventCutObj ( AliMixEventCutObj::kZVertex, -10, 10, 1 );

//    evPool->AddCut(centrality);
   evPool->AddCut(multi);
   evPool->AddCut ( zvertex );

   // adds event pool (comment it and u will have default mixing)
   mixHandler->SetEventPool ( evPool );

   mixHandler->SelectCollisionCandidates(AliVEvent::kAny);
//    mixHandler->SelectCollisionCandidates(AliVEvent::kAnyINT);


//    mixHandler->DoMixIfNotEnoughEvents(kFALSE);

   multiInputHandler->AddInputEventHandler ( mixHandler );

   // adds mixing info task
   gROOT->LoadMacro ( "AddAnalysisTaskMixInfo.C" );
   AddAnalysisTaskMixInfo (opts);

}
