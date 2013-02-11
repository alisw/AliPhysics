#ifndef __CINT__
#include <AliRsnCutPID.h>
#include <AliRsnInputHandler.h>
#include <AliRsnCutSet.h>
#include <AliRsnCutValue.h>
#endif

Int_t AddRsnDaughterCutsPDG(AliPID::EParticleType type1,AliPID::EParticleType type2,TString opt,AliRsnInputHandler *rsnIH=0,AliAnalysisTaskSE *task=0)
{

   if (!rsnIH) return 0;
   Bool_t valid;
   Int_t isRsnMini = AliRsnTrainManager::GetGlobalInt("IsRsnMini",valid);
   TString rsnQualityCut = AliRsnTrainManager::GetGlobalStr("RsnQualityCut",valid); 

   // === USER HAS TO SET CORRECT NUMBER OF CUTS SETS =====
   Int_t numberOfCuts = 1;

   // gets selector
   AliRsnDaughterSelector *sel = rsnIH->GetSelector();

   //---------------------------------------------
   //  Define single cuts
   //---------------------------------------------

   Bool_t useQuality = kFALSE;
   if (opt.Contains("quality")) {
      useQuality = kTRUE;
   }


   Double_t etaRange=0.8;

   AliRsnCutValue *cutEta;
   Bool_t useEta = kFALSE;
   if (opt.Contains("eta")) {
      Printf("Using ETA range (%.2f,%.2f)",-etaRange,etaRange);
      useEta = kTRUE;
   }

   AliRsnCutSet *cuts1 = new AliRsnCutSet(Form("%sPDG%s",AliPID::ParticleName(type1),opt.Data()), AliRsnTarget::kDaughter);

   Double_t nSigmaTPC=3.0;
   Double_t nSigmaTOF=3.0;
   Double_t ptTPCMax=0.8;

   TString scheme;

   AliRsnCutPID *cut1 = new AliRsnCutPID(Form("cut%sPDG%s",AliPID::ParticleName(type1),opt.Data()),type1,0.0,kTRUE);
   cuts1->AddCut(cut1);
   if (!scheme.IsNull()) scheme += "&";
   scheme += cut1->GetName();
   if (useEta) {
      AliRsnCutValue *cutEta1 = new AliRsnCutValue(Form("cut%sETA%s",AliPID::ParticleName(type1),opt.Data()),-etaRange,etaRange);
      AliRsnValueDaughter *valEta1 = new AliRsnValueDaughter(Form("val%sETA%s",AliPID::ParticleName(type1)),AliRsnValueDaughter::kEta);
      cutEta1->SetValueObj(valEta1);
      cuts1->AddCut(cutEta1);
      if (!scheme.IsNull()) scheme += "&";
      scheme += cutEta1->GetName();
   }
   if (useQuality) {
      AliRsnCutTrackQuality *qualityCut1 = new AliRsnCutTrackQuality("cutQuatityPDG1");
      if (!rsnQualityCut.IsNull()) {
         AliESDtrackCuts *esdTK = RsnQualityCut(rsnQualityCut.Data());
         qualityCut1->SetESDtrackCuts(esdTK);
      } else {
         qualityCut1->SetDefaults2010();
      }
      cuts1->AddCut(qualityCut1);
      if (!scheme.IsNull()) scheme += "&";
      scheme += qualityCut1->GetName();
   }
   cuts1->SetCutScheme(scheme.Data());
   sel->Add(cuts1, kTRUE);

   scheme = "";
   AliRsnCutSet *cuts2 = 0;
   if (type1 != type2) {
      AliRsnCutPID *cut2 = new AliRsnCutPID(Form("cut%sPDG%s",AliPID::ParticleName(type2),opt.Data()),type2,0.0,kTRUE);

      cuts2 = new AliRsnCutSet(Form("%sPDG%s",AliPID::ParticleName(type2),opt.Data()), AliRsnTarget::kDaughter);
      cuts2->AddCut(cut2);
      if (!scheme.IsNull()) scheme += "&";
      scheme += cut2->GetName();
      if (useQuality) {
         AliRsnCutTrackQuality *qualityCut2 = new AliRsnCutTrackQuality("cutQuatityPDG2");
         qualityCut2->SetDefaults2010();
         cuts2->AddCut(qualityCut2);
         if (!scheme.IsNull()) scheme += "&";
         scheme += qualityCut2->GetName();
      }
      if (useEta) {
         AliRsnCutValue *cutEta2 = new AliRsnCutValue(Form("cut%sETA%s",AliPID::ParticleName(type2),opt.Data()),-etaRange,etaRange);
         AliRsnValueDaughter *valEta2 = new AliRsnValueDaughter(Form("val%sETA%s",AliPID::ParticleName(type2)),AliRsnValueDaughter::kEta);
         cutEta2->SetValueObj(valEta2);
         cuts2->AddCut(cutEta2);
         if (!scheme.IsNull()) scheme += "&";
         scheme += cutEta2->GetName();
      }

      cuts2->SetCutScheme(scheme.Data());
      sel->Add(cuts2, kTRUE);
      numberOfCuts++;
   }
   if (opt.Contains("mon")) {
      AddMonitorOutput(cuts1->GetMonitorOutput(),opt);
      if (type1 != type2) AddMonitorOutput(cuts2->GetMonitorOutput());
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


