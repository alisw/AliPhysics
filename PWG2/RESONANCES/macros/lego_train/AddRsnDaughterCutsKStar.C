#ifndef __CINT__
#include <PWG2/RESONANCES/AliRsnInputHandler.h>
#include <PWG2/RESONANCES/AliRsnCutSet.h>
#include <PWG2/RESONANCES/AliRsnCutDaughterKStar2010PP.h>
#endif
Int_t AddRsnDaughterCutsKStar(AliPID::EParticleType type1,AliPID::EParticleType type2,TString opt,Bool_t isRsnMini=kFALSE,AliRsnInputHandler *rsnIH=0,AliAnalysisTaskSE *task=0)
{

   if (!rsnIH) return 0;

   // === USER HAS TO SET CORRECT NUMBER OF CUTS SETS =====
   Int_t numberOfCuts = 2;

   Printf("AddRsnDaughterCutsKStar Option : %s",opt.Data());

   AliRsnCutPion *cutPi = 0;
   AliRsnCutKaon *cutK = 0;
   if (opt.Contains("TPCTOFpidDefaultKstarPP2010")) {
      cutPi = new AliRsnCutPion("cutPionTPCTOFpidDefaultKstarPP2010", AliRsnCutPion::kTPCTOFpidDefaultKstarPP2010);
      cutK = new AliRsnCutKaon("cutKaonTPCTOFpidDefaultKstarPP2010", AliRsnCutKaon::kTPCTOFpidDefaultKstarPP2010);
   } else if (opt.Contains("FastTPCpid1point5sigma")) {
      cutPi = new AliRsnCutPion("cutPionForKStarFastTPCpid1point5sigma", AliRsnCutPion::kFastTPCpid1point5sigma);
      cutK = new AliRsnCutKaon("cutKaonForKStarFastTPCpid1point5sigma", AliRsnCutKaon::kFastTPCpid1point5sigma);
   } else if (opt.Contains("FastTPCpid2sigma")) {
      cutPi = new AliRsnCutPion("cutPionForKStarFastTPCpid2sigma", AliRsnCutPion::kFastTPCpid2sigma);
      cutK = new AliRsnCutKaon("cutKaonForKStarFastTPCpid2sigma", AliRsnCutKaon::kFastTPCpid2sigma);
   } else {
      cutPi = new AliRsnCutPion("cutPionDefault");
      cutK = new AliRsnCutKaon("cutKaonDefault")
   }
   AliRsnCutSet *cutSetPi = new AliRsnCutSet(Form("set%s",cutPi->GetName()), AliRsnTarget::kDaughter);
   cutSetPi->AddCut(cutPi);
   cutSetPi->SetCutScheme(cutPi->GetName());

   // cut set
   AliRsnCutSet *cutSetK = new AliRsnCutSet(Form("set%s",cutK->GetName()), AliRsnTarget::kDaughter);
   cutSetK->AddCut(cutK);
   cutSetK->SetCutScheme(cutK->GetName());


   if (opt.Contains("mon")) {
      Printf("Monitoring cut AddRsnDaughterCutsKStar Option : %s",opt.Data());
      AddMonitorOutput(cutSetPi->GetMonitorOutput());
      AddMonitorOutput(cutSetK->GetMonitorOutput());
   }
   if (isRsnMini) {
      AliRsnMiniAnalysisTask *taskRsnMini = dynamic_cast<AliRsnMiniAnalysisTask *>(task);
      if (taskRsnMini) {
         taskRsnMini->AddTrackCuts(cutSetPi);
         taskRsnMini->AddTrackCuts(cutSetK);
      }
   } else {
      AliRsnDaughterSelector *sel = rsnIH->GetSelector();
      sel->Add(cutSetPi, kTRUE);
      sel->Add(cutSetK, kTRUE);
   }


   return numberOfCuts;
}
