#ifndef __CINT__
#include <AliRsnCutPIDNSigma.h>
#include <PWGLF/RESONANCES/AliRsnCutTrackQuality.h>
#endif
Int_t AddRsnDaughterCutsRhoNsigma(AliPID::EParticleType type1,AliPID::EParticleType type2,TString opt,AliRsnInputHandler *rsnIH=0,AliAnalysisTaskSE *task=0)
{

   if (!rsnIH) return 0;

   Bool_t valid;
   Int_t isRsnMini = AliRsnTrainManager::GetGlobalInt("IsRsnMini",valid);

   // === USER HAS TO SET CORRECT NUMBER OF CUTS SETS =====
   Int_t numberOfCuts = 1;

   //---------------------------------------------
   //  Define single cuts
   //---------------------------------------------

   Printf("AddRsnDaughterCutsRho Option : %s",opt.Data());

   Double_t nSigmaTPC=3.0;
   Double_t nSigmaTOF=3.0;
   Double_t etaRange=0.8;
   Bool_t useTPC_Pi=kFALSE;
   Bool_t useTOF_Pi=kFALSE;

   if (opt.Contains("qualityonly")) {
      useTPC_Pi=kFALSE;
      useTOF_Pi=kFALSE;
   } else if (!opt.Contains("nsig")) {
      useTPC_Pi=kTRUE;
      useTOF_Pi=kTRUE;
   }

   if (opt.Contains("PiTPCnsig"))  useTPC_Pi=kTRUE;
   if (opt.Contains("PiTOFnsig"))  useTOF_Pi=kTRUE;

   if (opt.Contains("PiTPCnsig10")) nSigmaTPC = 1.0;
   if (opt.Contains("PiTPCnsig15")) nSigmaTPC = 1.5;
   if (opt.Contains("PiTPCnsig20")) nSigmaTPC = 2.0;
   if (opt.Contains("PiTPCnsig25")) nSigmaTPC = 2.5;
   if (opt.Contains("PiTPCnsig30")) nSigmaTPC = 3.0;

   if (opt.Contains("PiTOFnsig10")) nSigmaTOF = 1.0;
   if (opt.Contains("PiTOFnsig15")) nSigmaTOF = 1.5;
   if (opt.Contains("PiTOFnsig20")) nSigmaTOF = 2.0;
   if (opt.Contains("PiTOFnsig25")) nSigmaTOF = 2.5;
   if (opt.Contains("PiTOFnsig30")) nSigmaTOF = 3.0;

   Bool_t usePDG=kFALSE;
   if (opt.Contains("pdg")) {
      Printf("Using PDG");
      usePDG = kTRUE;
   }

   Bool_t useEta = kFALSE;
   if (opt.Contains("eta")) {
      Printf("Using ETA range (%.2f,%.2f)",-etaRange,etaRange);
      useEta = kTRUE;
   }

//---------------------------------------------
//  Combine cuts
//---------------------------------------------

   TString cutname = "Pi_Rho";
   if (!opt.IsNull()) cutname += Form("_%s",opt.Data());
   AliRsnCutSet *cuts = new AliRsnCutSet(cutname.Data(), AliRsnTarget::kDaughter);

   TString scheme="";

   AliRsnCutTrackQuality *qualityCut = new AliRsnCutTrackQuality("cutQuatityPi");
   //qualityCut->SetDefaults2010();

   qualityCut->SetDCAZmax(0.2);
   qualityCut->SetDCARmax(0.02);
   qualityCut->AddStatusFlag(AliESDtrack::kTPCin   , kTRUE);
   qualityCut->AddStatusFlag(AliESDtrack::kTPCrefit, kTRUE);
   qualityCut->AddStatusFlag(AliESDtrack::kITSrefit, kTRUE);

   qualityCut->SetPtRange(0.2, 1E+20);
   qualityCut->SetEtaRange(-0.8, 0.8);
   qualityCut->SetSPDminNClusters(0);
   qualityCut->SetITSminNClusters(0);
   qualityCut->SetITSmaxChi2(1E+20);
   qualityCut->SetTPCminNClusters(70);
   qualityCut->SetTPCmaxChi2(4.0);
   qualityCut->SetRejectKinkDaughters();

   cuts->AddCut(qualityCut);
   if (!scheme.IsNull()) scheme += "&";
   scheme += qualityCut->GetName();

   if (useTPC_Pi) {
      AliRsnCutPIDNSigma *cutPiTPC = new AliRsnCutPIDNSigma("cutPIDNSigmaTPCPi",AliPID::kPion,AliRsnCutPIDNSigma::kTPC);
      cutPiTPC->SinglePIDRange(nSigmaTPC);
      cuts->AddCut(cutPiTPC);
      if (!scheme.IsNull()) scheme += "&";
      scheme += cutPiTPC->GetName();
   }
   if (useTOF_Pi) {
      AliRsnCutPIDNSigma *cutPiTOF = new AliRsnCutPIDNSigma("cutPIDNSigmaTOFPi",AliPID::kPion,AliRsnCutPIDNSigma::kTOF);
      cutPiTOF->SinglePIDRange(nSigmaTOF);
      cuts->AddCut(cutPiTOF);
      if (!scheme.IsNull()) scheme += "&";
      scheme += cutPiTOF->GetName();
   }
   if (useEta) {
      Printf("Adding ETA ...");
      AliRsnValueDaughter *valEta = new AliRsnValueDaughter(Form("val%sETA%s",AliPID::ParticleName(type1),opt.Data()),AliRsnValueDaughter::kEta);
      AliRsnCutValue *cutEta = new AliRsnCutValue(Form("cut%sETA%s",AliPID::ParticleName(type1),opt.Data()),-etaRange,etaRange);
      cutEta->SetTargetType(AliRsnTarget::kDaughter);
      cutEta->SetValueObj(valEta);
      cuts->AddCut(cutEta);
      if (!scheme.IsNull()) scheme += "&";
      scheme += cutEta->GetName();
   }
   if (usePDG) {
      Printf("Adding PDG ...");
      AliRsnCutPID *cutPDG = new AliRsnCutPID(Form("cut%sPDG%s",AliPID::ParticleName(type1),opt.Data()),type1,0.0,kTRUE);
      cuts->AddCut(cutPDG);
      if (!scheme.IsNull()) scheme += "&";
      scheme += cutPDG->GetName();
   }

   Printf ("CUT Scheme is '%s'",scheme.Data());
   cuts->SetCutScheme(scheme.Data());

   if (opt.Contains("mon")) {
      AddMonitorOutput(cuts->GetMonitorOutput(),opt);
   }
   if (isRsnMini) {
      AliRsnMiniAnalysisTask *taskRsnMini = dynamic_cast<AliRsnMiniAnalysisTask *>(task);
      if (taskRsnMini) {
         taskRsnMini->AddTrackCuts(cuts);
      }
   } else {
      AliRsnDaughterSelector *sel = rsnIH->GetSelector();
//       sel->SetLabelCheck(kFALSE);
      sel->Add(cuts, kTRUE);
   }
   return numberOfCuts;

}

