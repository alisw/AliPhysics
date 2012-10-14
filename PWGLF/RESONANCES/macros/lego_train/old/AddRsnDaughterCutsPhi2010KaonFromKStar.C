#ifndef __CINT__
#include <PWG2/RESONANCES/AliRsnInputHandler.h>
#include <PWG2/RESONANCES/AliRsnCutKaonForPhi2010.h>

#endif
Int_t AddRsnDaughterCutsPhi2010KaonFromKStar(AliPID::EParticleType type1,AliPID::EParticleType type2,TString opt,Bool_t isRsnMini=kFALSE,AliRsnInputHandler *rsnIH=0,AliAnalysisTaskSE *task=0)
{

   if (!rsnIH) return 0;

   // === USER HAS TO SET CORRECT NUMBER OF CUTS SETS =====
   Int_t numberOfCuts = 1;

   Printf("Option : %s",opt.Data());

   //---------------------------------------------
   //  Define single cuts
   //---------------------------------------------

   AliRsnCutDaughterKStar2010PP *cutK = new AliRsnCutDaughterKStar2010PP("cutPhiKaonForKStar", AliPID::kKaon);
   // cut set
   AliRsnCutSet *cutSetK = new AliRsnCutSet("PhiKaonForKStar", AliRsnTarget::kDaughter);
   cutSetK->AddCut(cutK);
   cutSetK->SetCutScheme(cutK->GetName());
   if (opt.Contains("mon")) {
      AddMonitorOutput(cuts->GetMonitorOutput());
   }
   if (isRsnMini) {
      AliRsnMiniAnalysisTask *taskRsnMini = dynamic_cast<AliRsnMiniAnalysisTask *>(task);
      if (taskRsnMini) {
         taskRsnMini->AddTrackCuts(cutSetK);
      }
   } else {
      AliRsnDaughterSelector *sel = rsnIH->GetSelector();
      sel->Add(cutSetK, kTRUE);
   }


   return numberOfCuts;

}
