#ifndef __CINT__
#include <ANALYSIS/AliAnalysisTaskSE.h>
#include <PWG2/RESONANCES/AliRsnCutSet.h>
#include <PWG2/RESONANCES/AliRsnInputHandler.h>
#include <PWG2/RESONANCES/AliRsnCutKaonForPhi2010.h>
#include <PWG2/RESONANCES/AliRsnMiniAnalysisTask.h>
#include <PWG2/RESONANCES/AliRsnAnalysisTask.h>
#endif
Int_t AddRsnDaughterCutsPhi2010(AliPID::EParticleType type1,AliPID::EParticleType type2,TString opt,Bool_t isRsnMini=kFALSE,AliRsnInputHandler *rsnIH=0,AliAnalysisTaskSE *task=0)
{

   if (!rsnIH) return 0;

   // === USER HAS TO SET CORRECT NUMBER OF CUTS SETS =====
   Int_t numberOfCuts = 1;

   //---------------------------------------------
   //  Define single cuts
   //---------------------------------------------

   Printf("AddRsnDaughterCutsPhi2010 Option : %s",opt.Data());

   AliRsnCutKaonForPhi2010 *cut = new AliRsnCutKaonForPhi2010("cutKaonPhi2010",3.0,3.0,0.8);
   if (opt.Contains("qualityonly")) cut->SetMode(AliRsnCutKaonForPhi2010::kQuality);

   if (opt.Contains("tpconly")) {
      cut->SetMode(AliRsnCutKaonForPhi2010::kOnlyTPC);
      if (opt.Contains("sigma1")) cut->SetCutTPC(1.0);
      if (opt.Contains("sigma2")) cut->SetCutTPC(2.0);
      if (opt.Contains("sigma3")) cut->SetCutTPC(3.0);
   }

   if (opt.Contains("tofonly")) {
      cut->SetMode(AliRsnCutKaonForPhi2010::kOnlyTOF);
      if (opt.Contains("sigma1")) cut->SetCutTOF(1.0);
      if (opt.Contains("sigma2")) cut->SetCutTOF(2.0);
      if (opt.Contains("sigma3")) cut->SetCutTOF(3.0);
   }

   //---------------------------------------------
   //  Combine cuts
   //---------------------------------------------
   TString cutname = "kaonPhi2010";
   if (!opt.IsNull()) cutname += Form("_%s",opt.Data());
   AliRsnCutSet *cuts = new AliRsnCutSet(cutname.Data(), AliRsnTarget::kDaughter);
   cuts->AddCut(cut);
   cuts->SetCutScheme(cut->GetName());

   if (opt.Contains("mon")) {
      AddMonitorOutput(cuts->GetMonitorOutput());
   }
   if (isRsnMini) {
      AliRsnMiniAnalysisTask *taskRsnMini = dynamic_cast<AliRsnMiniAnalysisTask *>(task);
      if (taskRsnMini) {
         taskRsnMini->AddTrackCuts(cuts);
      }
   } else {
      AliRsnDaughterSelector *sel = rsnIH->GetSelector();
      sel->Add(cuts, kTRUE);
   }
   return numberOfCuts;

}
