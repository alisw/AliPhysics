#ifndef __CINT__
#include <PWG2/RESONANCES/AliRsnInputHandler.h>
#include <PWG2/RESONANCES/AliRsnCutSet.h>
#include <PWG2/RESONANCES/AliRsnCutDaughterKStar2010PP.h>
#endif
Int_t AddRsnDaughterCutsKStar2010(AliPID::EParticleType type1,AliPID::EParticleType type2,TString opt,Bool_t isRsnMini=kFALSE,AliRsnInputHandler *rsnIH=0,AliAnalysisTaskSE *task=0)
{

   if (!rsnIH) return 0;

   // === USER HAS TO SET CORRECT NUMBER OF CUTS SETS =====
   Int_t numberOfCuts = 2;

   Printf("AddRsnDaughterCutsKStar2010 Option : %s",opt.Data());

   // integrated kaon cut
   AliRsnCutDaughterKStar2010PP *cutK = new AliRsnCutDaughterKStar2010PP("cutKaonForKStar", type1);
   // cut set
   AliRsnCutSet *cutSetK = new AliRsnCutSet("KaonForKStar", AliRsnTarget::kDaughter);
   cutSetK->AddCut(cutK);
   cutSetK->SetCutScheme(cutK->GetName());

   // integrated proton cut
   AliRsnCutDaughterKStar2010PP *cutP = new AliRsnCutDaughterKStar2010PP("cutPionForKStar", type2);
   // cut set
   AliRsnCutSet *cutSetP = new AliRsnCutSet("PionForKStar", AliRsnTarget::kDaughter);
   cutSetP->AddCut(cutP);
   cutSetP->SetCutScheme(cutP->GetName());

   if (opt.Contains("mon")) {
      AddMonitorOutput(cutSetK->GetMonitorOutput());
      AddMonitorOutput(cutSetP->GetMonitorOutput());
   }
   if (isRsnMini) {
      AliRsnMiniAnalysisTask *taskRsnMini = dynamic_cast<AliRsnMiniAnalysisTask *>(task);
      if (taskRsnMini) {
         taskRsnMini->AddTrackCuts(cutSetK);
         taskRsnMini->AddTrackCuts(cutSetP);
      }
   } else {
      AliRsnDaughterSelector *sel = rsnIH->GetSelector();
      sel->Add(cutSetK, kTRUE);
      sel->Add(cutSetP, kTRUE);
   }


   return numberOfCuts;
}
