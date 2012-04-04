#ifndef __CINT__
#endif
Int_t AddRsnDaughterCutsLambdaNsigma(AliPID::EParticleType type1,AliPID::EParticleType type2,TString opt,Bool_t isRsnMini=kFALSE,AliRsnInputHandler *rsnIH=0,AliAnalysisTaskSE *task=0)
{

   if (!rsnIH) return 0;

   Bool_t valid = kTRUE;
   Int_t isPP = AliAnalysisManager::GetGlobalInt("rsnIsPP",valid);

   Bool_t usePPCut = kFALSE;

   if (isPP && (opt.Contains("usePP"))) usePPCut = kTRUE;


   // === USER HAS TO SET CORRECT NUMBER OF CUTS SETS =====
   Int_t numberOfCuts = 2;

   //---------------------------------------------
   //  Define single cutsP
   //---------------------------------------------

   Printf("AddRsnDaughterCutsLambda Option : %s",opt.Data());

   Double_t nSigmaTPC_P=3.0;
   Double_t nSigmaTPC_K=3.0;
   Double_t nSigmaTOF_P=3.0;
   Double_t nSigmaTOF_K=3.0;
   Double_t etaRange=0.8;

   // default values
   Double_t nSigmaTPC_P=3.0;
   Double_t nSigmaTPC_K=3.0;
   Double_t nSigmaTOF_P=3.0;
   Double_t nSigmaTOF_K=3.0;
   Double_t etaRange=0.8;

   Bool_t useTPC_P=kFALSE;
   Bool_t useTOF_P=kFALSE;
   Bool_t useTPC_K=kFALSE;
   Bool_t useTOF_K=kFALSE;

   if (opt.Contains("qualityonly")) {
      useTPC_P=kFALSE;
      useTOF_P=kFALSE;
      useTPC_K=kFALSE;
      useTOF_K=kFALSE;
   } else if (!opt.Contains("nsig")) {
      useTPC_P=kTRUE;
      useTOF_P=kTRUE;
      useTPC_K=kTRUE;
      useTOF_K=kTRUE;
   }

   if (opt.Contains("PTPCnsig")) useTPC_P=kTRUE;
   if (opt.Contains("PTOFnsig")) useTOF_P=kTRUE;
   if (opt.Contains("KTPCnsig"))  useTPC_K=kTRUE;
   if (opt.Contains("KTOFnsig"))  useTOF_K=kTRUE;

   if (opt.Contains("PTPCnsig10")) nSigmaTPC_P = 1.0;
   if (opt.Contains("PTPCnsig15")) nSigmaTPC_P = 1.5;
   if (opt.Contains("PTPCnsig20")) nSigmaTPC_P = 2.0;
   if (opt.Contains("PTPCnsig25")) nSigmaTPC_P = 2.5;
   if (opt.Contains("PTPCnsig30")) nSigmaTPC_P = 3.0;

   if (opt.Contains("KTPCnsig10")) nSigmaTPC_K = 1.0;
   if (opt.Contains("KTPCnsig15")) nSigmaTPC_K = 1.5;
   if (opt.Contains("KTPCnsig20")) nSigmaTPC_K = 2.0;
   if (opt.Contains("KTPCnsig25")) nSigmaTPC_K = 2.5;
   if (opt.Contains("KTPCnsig30")) nSigmaTPC_K = 3.0;

   if (opt.Contains("PTOFnsig10")) nSigmaTOF_P = 1.0;
   if (opt.Contains("PTOFnsig15")) nSigmaTOF_P = 1.5;
   if (opt.Contains("PTOFnsig20")) nSigmaTOF_P = 2.0;
   if (opt.Contains("PTOFnsig25")) nSigmaTOF_P = 2.5;
   if (opt.Contains("PTOFnsig30")) nSigmaTOF_P = 3.0;

   if (opt.Contains("KTOFnsig10")) nSigmaTOF_K = 1.0;
   if (opt.Contains("KTOFnsig15")) nSigmaTOF_K = 1.5;
   if (opt.Contains("KTOFnsig20")) nSigmaTOF_K = 2.0;
   if (opt.Contains("KTOFnsig25")) nSigmaTOF_K = 2.5;
   if (opt.Contains("KTOFnsig30")) nSigmaTOF_K = 3.0;


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

   // PROTON SETTINGS ===========================================

   TString scheme="";
   TString cutname = "p_Lambda";
   if (!opt.IsNull()) cutname += Form("_%s",opt.Data());
   AliRsnCutSet *cutsP = new AliRsnCutSet(cutname.Data(), AliRsnTarget::kDaughter);

   AliRsnCutTrackQuality *qualityCutP = new AliRsnCutTrackQuality("cutQuatity{");
   qualityCutP->SetDefaults2010();
   cutsP->AddCut(qualityCutP);
   if (!scheme.IsNull()) scheme += "&";
   scheme += qualityCutP->GetName();

   if (useTPC_P) {
      AliRsnCutPIDNSigma *cutPTPC = new AliRsnCutPIDNSigma("cutNSigmaTPCP",AliPID::kProton,AliRsnCutPIDNSigma::kTPC);
      cutPTPC->SinglePIDRange(nSigmaTPC_P);
      cutsP->AddCut(cutPTPC);
      if (!scheme.IsNull()) scheme += "&";
      scheme += cutPTPC->GetName();
   }

   if (useTOF_P) {
      AliRsnCutPIDNSigma *cutPTOF = new AliRsnCutPIDNSigma("cutNSigmaTOFP",AliPID::kProton,AliRsnCutPIDNSigma::kTOF);
      cutPTOF->SinglePIDRange(nSigmaTOF_P);
      cutsP->AddCut(cutPTOF);
      if (!scheme.IsNull()) scheme += "&";
      scheme += cutPTOF->GetName();
   }

   if (useEta) {
      AliRsnValueDaughter *valEtaP = new AliRsnValueDaughter(Form("val%sETA%s",AliPID::ParticleName(type1),opt.Data()),AliRsnValueDaughter::kEta);
      AliRsnCutValue *cutEtaP = new AliRsnCutValue(Form("cut%sETA%s",AliPID::ParticleName(type1),opt.Data()),-etaRange,etaRange);
      cutEtaP->SetTargetType(AliRsnTarget::kDaughter);
      cutEtaP->SetValueObj(valEtaP);
      cutsP->AddCut(cutEtaP);
      if (!scheme.IsNull()) scheme += "&";
      scheme += cutEtaP->GetName();
   }
   if (usePDG) {
      AliRsnCutPID *cutPDGP = new AliRsnCutPID(Form("cut%sPDG%s",AliPID::ParticleName(type1),opt.Data()),type1,0.0,kTRUE);
      cutsP->AddCut(cutPDGP);
      if (!scheme.IsNull()) scheme += "&";
      scheme += cutPDGP->GetName();
   }

   Printf ("CUT Scheme for PROTON is '%s'",scheme.Data());
   cutsP->SetCutScheme(scheme.Data());

   // END PROTON =======================================

   // KAON SETTINGS =======================================
   scheme="";
   cutname = "K_Lambda";
   if (!opt.IsNull()) cutname += Form("_%s",opt.Data());
   AliRsnCutSet *cutsK = new AliRsnCutSet(cutname.Data(), AliRsnTarget::kDaughter);

   AliRsnCutTrackQuality *qualityCutK = new AliRsnCutTrackQuality("cutQuatityK");
   qualityCutK->SetDefaults2010();
   cutsK->AddCut(qualityCutK);
   if (!scheme.IsNull()) scheme += "&";
   scheme += qualityCutK->GetName();


   if (useTPC_K) {
      AliRsnCutPIDNSigma *cutKTPC = new AliRsnCutPIDNSigma("cutNSigmaTPCK",AliPID::kKaon,AliRsnCutPIDNSigma::kTPC);
      cutKTPC->SinglePIDRange(nSigmaTPC_K);
      cutsK->AddCut(cutKTPC);
      if (!scheme.IsNull()) scheme += "&";
      scheme += cutKTPC->GetName();
   }

   if (useTOF_K) {
      AliRsnCutPIDNSigma *cutKTOF = new AliRsnCutPIDNSigma("cutNSigmaTOFK",AliPID::kKaon,AliRsnCutPIDNSigma::kTOF);
      cutKTOF->SinglePIDRange(nSigmaTOF_K);
      cutsK->AddCut(cutKTOF);
      if (!scheme.IsNull()) scheme += "&";
      scheme += cutKTOF->GetName();
   }
   
   if (useEta) {
      AliRsnValueDaughter *valEtaK = new AliRsnValueDaughter(Form("val%sETA%s",AliPID::ParticleName(type2),opt.Data()),AliRsnValueDaughter::kEta);
      AliRsnCutValue *cutEtaK = new AliRsnCutValue(Form("cut%sETA%s",AliPID::ParticleName(type2),opt.Data()),-etaRange,etaRange);
      cutEtaK->SetTargetType(AliRsnTarget::kDaughter);
      cutEtaK->SetValueObj(valEtaK);
      cutsK->AddCut(cutEtaK);
      if (!scheme.IsNull()) scheme += "&";
      scheme += cutEtaK->GetName();
   }
   if (usePDG) {
      AliRsnCutPID *cutPDGK = new AliRsnCutPID(Form("cut%sPDG%s",AliPID::ParticleName(type2),opt.Data()),type2,0.0,kTRUE);
      cutsK->AddCut(cutPDGK);
      if (!scheme.IsNull()) scheme += "&";
      scheme += cutPDGK->GetName();
   }

   Printf ("CUT Scheme for KAON is '%s'",scheme.Data());
   cutsK->SetCutScheme(scheme.Data());

   // END KAON =======================================

   if (opt.Contains("mon")) {
      AddMonitorOutput(cutsP->GetMonitorOutput(),opt);
      AddMonitorOutput(cutsK->GetMonitorOutput(),opt);
   }
   if (isRsnMini) {
      AliRsnMiniAnalysisTask *taskRsnMini = dynamic_cast<AliRsnMiniAnalysisTask *>(task);
      if (taskRsnMini) {
         taskRsnMini->AddTrackCuts(cutsP);
         taskRsnMini->AddTrackCuts(cutsK);
      }
   } else {
      AliRsnDaughterSelector *sel = rsnIH->GetSelector();
//       sel->SetLabelCheck(kFALSE);
      sel->Add(cutsP, kTRUE);
      sel->Add(cutsK, kTRUE);
   }
   return numberOfCuts;

}
