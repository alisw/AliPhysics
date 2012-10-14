#ifndef __CINT__
#include <PWG2/RESONANCES/AliRsnCutPID.h>
#include <PWG2/RESONANCES/AliRsnInputHandler.h>
#include <PWG2/RESONANCES/AliRsnCutSet.h>
#endif

Int_t AddRsnDaughterCutsBPID(AliPID::EParticleType type1,AliPID::EParticleType type2,TString opt,Bool_t isRsnMini=kFALSE,AliRsnInputHandler *rsnIH=0,AliAnalysisTaskSE *task=0)
{

   if (!rsnIH) return 0;

   // === USER HAS TO SET CORRECT NUMBER OF CUTS SETS =====
   Int_t numberOfCuts = 1;

   // gets selector
   AliRsnDaughterSelector *sel = rsnIH->GetSelector();

   //---------------------------------------------
   //  Define single cuts
   //---------------------------------------------

   AliRsnCutPID *cut1 = new AliRsnCutPID(Form("cut%sBPID%s",AliPID::ParticleName(type1),opt.Data()),type1,0.0,kFALSE);
   AliRsnCutSet *cuts1 = new AliRsnCutSet(Form("%sBPID%s",AliPID::ParticleName(type1),opt.Data()), AliRsnTarget::kDaughter);
   cuts1->AddCut(cut1);
   cuts1->SetCutScheme(cut1->GetName());
   sel->Add(cuts1, kTRUE);


   AliRsnCutSet *cuts2 = 0;
   if (type1 != type2) {
      AliRsnCutPID *cut2 = new AliRsnCutPID(Form("cut%sBPID%s",AliPID::ParticleName(type2),opt.Data()),type2,0.0,kFALSE);
      cuts2 = new AliRsnCutSet(Form("%sBPID%s",AliPID::ParticleName(type2),opt.Data()), AliRsnTarget::kDaughter);
      cuts2->AddCut(cut2);
      cuts2->SetCutScheme(cut2->GetName());
      sel->Add(cuts2, kTRUE);
      numberOfCuts++;
   }
   if (opt.Contains("mon")) {
      AddMonitorOutput(cuts1->GetMonitorOutput(),opt);
      if (type1 != type2) {
         AddMonitorOutput(cuts2->GetMonitorOutput(),opt);
      }
   }
   if (isRsnMini) {
      AliRsnMiniAnalysisTask *taskRsnMini = dynamic_cast<AliRsnMiniAnalysisTask *>(task);
      if (taskRsnMini) {
         taskRsnMini->AddTrackCuts(cuts1);
         if (type1 != type2) taskRsnMini->AddTrackCuts(cuts2);
      }
   } else {
      AliRsnDaughterSelector *sel = rsnIH->GetSelector();
      sel->Add(cuts1, kTRUE);
      if (type1 != type2)  sel->Add(cuts2, kTRUE);
   }


   return numberOfCuts;
}
