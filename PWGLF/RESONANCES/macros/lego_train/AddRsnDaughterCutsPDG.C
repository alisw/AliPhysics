#ifndef __CINT__
#include <PWG2/RESONANCES/AliRsnCutPID.h>
#include <PWG2/RESONANCES/AliRsnInputHandler.h>
#include <PWG2/RESONANCES/AliRsnCutSet.h>
#endif

Int_t AddRsnDaughterCutsPDG(AliPID::EParticleType type1,AliPID::EParticleType type2,TString opt,Bool_t isRsnMini=kFALSE,AliRsnInputHandler *rsnIH=0,AliAnalysisTaskSE *task=0)
{

   if (!rsnIH) return 0;

   // === USER HAS TO SET CORRECT NUMBER OF CUTS SETS =====
   Int_t numberOfCuts = 1;

   // gets selector
   AliRsnDaughterSelector *sel = rsnIH->GetSelector();

   //---------------------------------------------
   //  Define single cuts
   //---------------------------------------------

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
   AliRsnCutKaonForPhi2010 *cutQuality1 = new AliRsnCutKaonForPhi2010("cutKaonPhi2010",nSigmaTPC,nSigmaTOF,ptTPCMax);
   cutQuality1->SetMode(AliRsnCutKaonForPhi2010::kQuality);
   cuts1->AddCut(cutQuality1);

   AliRsnCutPID *cut1 = new AliRsnCutPID(Form("cut%sPDG%s",AliPID::ParticleName(type1),opt.Data()),type1,0.0,kTRUE);
   cuts1->AddCut(cut1);
   if (useEta) {
      AliRsnCutValue *cutEta1 = new AliRsnCutValue(Form("cut%sETA%s",AliPID::ParticleName(type1),opt.Data()),-etaRange,etaRange);
      AliRsnValueDaughter *valEta1 = new AliRsnValueDaughter(Form("val%sETA%s",AliPID::ParticleName(type1)),AliRsnValueDaughter::kEta);
      cutEta1->SetValueObj(valEta1);
      cuts1->AddCut(cutEta1);

      cuts1->SetCutScheme(Form("%s&%s&%s",cutQuality1->GetName(),cut1->GetName(),cutEta1->GetName()));
   } else {
      cuts1->SetCutScheme(Form("%s&%s",cutQuality1->GetName(),cut1->GetName()));
   }
   sel->Add(cuts1, kTRUE);

   AliRsnCutSet *cuts2 = 0;
   if (type1 != type2) {
      AliRsnCutPID *cut2 = new AliRsnCutPID(Form("cut%sPDG%s",AliPID::ParticleName(type2),opt.Data()),type2,0.0,kTRUE);
      AliRsnCutKaonForPhi2010 *cutQuality2 = new AliRsnCutKaonForPhi2010("cutKaonPhi2010",nSigmaTPC,nSigmaTOF,ptTPCMax);
      cutQuality2->SetMode(AliRsnCutKaonForPhi2010::kQuality);
      cuts2->AddCut(cutQuality2);

      cuts2 = new AliRsnCutSet(Form("%sPDG%s",AliPID::ParticleName(type2),opt.Data()), AliRsnTarget::kDaughter);
      cuts2->AddCut(cut2);
      if (useEta) {
         AliRsnCutValue *cutEta2 = new AliRsnCutValue(Form("cut%sETA%s",AliPID::ParticleName(type2),opt.Data()),-etaRange,etaRange);
         AliRsnValueDaughter *valEta2 = new AliRsnValueDaughter(Form("val%sETA%s",AliPID::ParticleName(type2)),AliRsnValueDaughter::kEta);
         cutEta2->SetValueObj(valEta2);
         cuts2->AddCut(cutEta2);
         cuts2->SetCutScheme(Form("%s&%s&%s",cutQuality2->GetName(),cut2->GetName(),cutEta2->GetName()));
      } else {
         cuts2->SetCutScheme(Form("%s&%s",cutQuality2->GetName(),cut2->GetName()));
      }
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


