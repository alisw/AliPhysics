#ifndef __CINT__
#include <AliAnalysisManager.h>
#include <AliMultiInputEventHandler.h>
#include <ANALYSIS/EventMixing/AliMixEventPool.h>
#include <ANALYSIS/EventMixing/AliMixEventCutObj.h>
#include <PWG2/RESONANCES/AliRsnAnalysisTask.h>
#include <PWG2/RESONANCES/AliRsnMiniAnalysisTask.h>
#endif

void AddMixingHandler ( AliMultiInputEventHandler *multiInputHandler,AliAnalysisTaskSE *task, TString format = "esd", Bool_t useMC = kFALSE,Bool_t isRsnMini=kFALSE,const Int_t mixNum = 10, TString opts = "" ) {

   Bool_t valid = kTRUE;
   Int_t isPP = AliAnalysisManager::GetGlobalInt("rsnIsPP",valid);

   if (isRsnMini) {
      AliRsnMiniAnalysisTask *taskRsn = (AliRsnMiniAnalysisTask *) task;

      // settings
      if (isPP) {
         taskRsn->UseMultiplicity("QUALITY");
         taskRsn->SetMaxDiffMult(1.0);
      } else {
         taskRsn->UseCentrality("V0M");
         taskRsn->SetMaxDiffMult(10.0);
      }
      // set mixing
      taskRsn->UseContinuousMix();
      //task->UseBinnedMix();
      taskRsn->SetNMix(mixNum);
      taskRsn->SetMaxDiffVz(1.0);
//      taskRsn->SetMaxDiffMult(10.0);
      taskRsn->SetMaxDiffAngle(1E20);

   } else {
      if ( !multiInputHandler ) return;

      AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
      const Int_t bufferSize = 1;
      AliMixInputEventHandler *mixHandler = new AliMixInputEventHandler ( bufferSize, mixNum );
      mixHandler->SetInputHandlerForMixing ( dynamic_cast<AliMultiInputEventHandler *> ( mgr->GetInputEventHandler() ) );
      AliMixEventPool *evPool = new AliMixEventPool();

      AliMixEventCutObj *multi = 0;

      if (isPP) {
         multi = new AliMixEventCutObj ( AliMixEventCutObj::kMultiplicity, 2, 102, 10 );
      } else {
         multi = new AliMixEventCutObj ( AliMixEventCutObj::kCentrality, 0, 100, 10, "V0M" );
      }
      AliMixEventCutObj *zvertex = new AliMixEventCutObj ( AliMixEventCutObj::kZVertex, -10, 10, 1 );

      evPool->AddCut(multi);
      evPool->AddCut ( zvertex );

      // adds event pool (comment it and u will have default mixing)
      mixHandler->SetEventPool ( evPool );

//    mixHandler->SelectCollisionCandidates(AliVEvent::kAny);

//    mixHandler->DoMixIfNotEnoughEvents(kFALSE);

      multiInputHandler->AddInputEventHandler ( mixHandler );

      // adds mixing info task
      RsnLoadMacroFromMixHandler( "AddAnalysisTaskMixInfo.C" );
      AddAnalysisTaskMixInfo (opts );

   }
}

Bool_t RsnLoadMacroFromMixHandler(TString macro,TString path="") {

   Bool_t valid;
   TString lego_path = AliAnalysisManager::GetGlobalStr("rsnLegoTrainPath",valid);
   if (!valid) lego_path = "$ALICE_ROOT/PWG2/RESONANCES/macros/lego_train";

   if (!gSystem->AccessPathName(macro.Data())) {
      gROOT->LoadMacro(macro.Data());
      Printf("Macro loaded from %s/%s ...",gSystem->pwd(),macro.Data());
      return kTRUE;
   }

   if (!gSystem->AccessPathName(gSystem->ExpandPathName(Form("%s/%s",lego_path.Data(),macro.Data())))) {
      gROOT->LoadMacro(gSystem->ExpandPathName(Form("%s/%s",lego_path.Data(),macro.Data())));
      Printf("Macro loaded from %s ...",gSystem->ExpandPathName(Form("%s/%s",lego_path.Data(),macro.Data())));
      return kTRUE;
   }

   Printf("Error loading %s",macro.Data());

   return kFALSE;
}
