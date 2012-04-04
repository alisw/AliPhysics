#ifndef __CINT__
#endif
Int_t AddRsnDaughterCutsPhiNsigma(AliPID::EParticleType type1,AliPID::EParticleType type2,TString opt,Bool_t isRsnMini=kFALSE,AliRsnInputHandler *rsnIH=0,AliAnalysisTaskSE *task=0)
{

   if (!rsnIH) return 0;

   Bool_t valid = kTRUE;
   Int_t isPP = AliAnalysisManager::GetGlobalInt("rsnIsPP",valid);

   Bool_t usePPCut = kFALSE;

   if (isPP && (opt.Contains("usePP"))) usePPCut = kTRUE;


   // === USER HAS TO SET CORRECT NUMBER OF CUTS SETS =====
   Int_t numberOfCuts = 1;

   //---------------------------------------------
   //  Define single cuts
   //---------------------------------------------

   Printf("AliRsnCutPIDNSigma Option : %s",opt.Data());

   Double_t nSigmaTPC=3.0;
   Double_t nSigmaTOF=3.0;
   Double_t etaRange=0.8;

   Bool_t useTPC_K=kFALSE;
   Bool_t useTOF_K=kFALSE;

   if (opt.Contains("qualityonly")) {
      useTPC_K=kFALSE;
      useTOF_K=kFALSE;
   } else if (!opt.Contains("nsig")) {
      useTPC_K=kTRUE;
      useTOF_K=kTRUE;
   }

   if (opt.Contains("KTPCnsig"))  useTPC_K=kTRUE;
   if (opt.Contains("KTOFnsig"))  useTOF_K=kTRUE;

   if (opt.Contains("KTPCnsig10")) nSigmaTPC = 1.0;
   if (opt.Contains("KTPCnsig15")) nSigmaTPC = 1.5;
   if (opt.Contains("KTPCnsig20")) nSigmaTPC = 2.0;
   if (opt.Contains("KTPCnsig25")) nSigmaTPC = 2.5;
   if (opt.Contains("KTPCnsig30")) nSigmaTPC = 3.0;

   if (opt.Contains("KTOFnsig10")) nSigmaTOF = 1.0;
   if (opt.Contains("KTOFnsig15")) nSigmaTOF = 1.5;
   if (opt.Contains("KTOFnsig20")) nSigmaTOF = 2.0;
   if (opt.Contains("KTOFnsig25")) nSigmaTOF = 2.5;
   if (opt.Contains("KTOFnsig30")) nSigmaTOF = 3.0;

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

   TString cutname = "K_Phi";
   if (!opt.IsNull()) cutname += Form("_%s",opt.Data());
   AliRsnCutSet *cuts = new AliRsnCutSet(cutname.Data(), AliRsnTarget::kDaughter);

   TString scheme="";

   AliRsnCutTrackQuality *qualityCut = new AliRsnCutTrackQuality("cutQuatityK");
   qualityCut->SetDefaults2010();
   cuts->AddCut(qualityCut);
   if (!scheme.IsNull()) scheme += "&";
   scheme += qualityCut->GetName();


   if (useTPC_K) {
      AliRsnCutPIDNSigma *cutKTPC = new AliRsnCutPIDNSigma("cutPIDNSigmaTPCK",AliPID::kKaon,AliRsnCutPIDNSigma::kTPC);
      cutKTPC->SinglePIDRange(nSigmaTPC);
      cuts->AddCut(cutKTPC);
      if (!scheme.IsNull()) scheme += "&";
      scheme += cutKTPC->GetName();
   }

   if (useTOF_K) {
      AliRsnCutPIDNSigma *cutKTOF = new AliRsnCutPIDNSigma("cutPIDNSigmaTOFK",AliPID::kKaon,AliRsnCutPIDNSigma::kTOF);
      cutKTOF->SinglePIDRange(nSigmaTOF);
      cuts->AddCut(cutKTOF);
      if (!scheme.IsNull()) scheme += "&";
      scheme += cutKTOF->GetName();
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
      AliRsnMiniAnalysisTask *taskRsnMini = (AliRsnMiniAnalysisTask*)task;
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

