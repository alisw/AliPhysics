#ifndef __CINT__
#endif
Int_t AddRsnDaughterCutsLambdaNsigma(AliPID::EParticleType type1,AliPID::EParticleType type2,TString opt,AliRsnInputHandler *rsnIH=0,AliAnalysisTaskSE *task=0)
{

   if (!rsnIH) return 0;

   Bool_t valid = kTRUE;
   //   Int_t collisionType = AliRsnTrainManager::GetGlobalInt("IsCollisionType",valid);
   Int_t useCommonQualityCut = AliRsnTrainManager::GetGlobalInt("RsnCommonQualityCut",valid);
   TString rsnQualityCut = AliRsnTrainManager::GetGlobalStr("RsnQualityCut",valid);
   Int_t isMC = AliRsnTrainManager::GetGlobalInt("IsMC",valid);
   Int_t isRsnMini = AliRsnTrainManager::GetGlobalInt("IsRsnMini",valid);
   Int_t isMixing = AliRsnTrainManager::GetGlobalInt("IsMixing",valid);

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
   Double_t PtMin_P=0.15;
   Double_t PtMax_P=1.e10;
   Double_t PtMin_K=0.15;
   Double_t PtMax_K=1.e10;
   Double_t PMax_P=1.1;
   Double_t PMax_K=0.6;

   Bool_t useTPC_P=kFALSE;
   Bool_t useTOF_P=kFALSE;
   Bool_t rejectUnmatchedTOF_P=kTRUE;
   Bool_t useTPC_K=kFALSE;
   Bool_t useTOF_K=kFALSE;
   Bool_t rejectUnmatchedTOF_K=kTRUE;

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
   if (opt.Contains("PTPCnsig10")) nSigmaTPC_P = 1.0;
   if (opt.Contains("PTPCnsig15")) nSigmaTPC_P = 1.5;
   if (opt.Contains("PTPCnsig20")) nSigmaTPC_P = 2.0;
   if (opt.Contains("PTPCnsig25")) nSigmaTPC_P = 2.5;
   if (opt.Contains("PTPCnsig30")) nSigmaTPC_P = 3.0;
   if (opt.Contains("PTPCnsig40")) nSigmaTPC_P = 4.0;
   if (opt.Contains("PTPCnsig50")) nSigmaTPC_P = 5.0;
   if (opt.Contains("PTPCnsig1000")) nSigmaTPC_P = 100.0;

   if (opt.Contains("KTPCnsig"))  useTPC_K=kTRUE;
   if (opt.Contains("KTPCnsig10")) nSigmaTPC_K = 1.0;
   if (opt.Contains("KTPCnsig15")) nSigmaTPC_K = 1.5;
   if (opt.Contains("KTPCnsig20")) nSigmaTPC_K = 2.0;
   if (opt.Contains("KTPCnsig25")) nSigmaTPC_K = 2.5;
   if (opt.Contains("KTPCnsig30")) nSigmaTPC_K = 3.0;
   if (opt.Contains("KTPCnsig40")) nSigmaTPC_K = 4.0;
   if (opt.Contains("KTPCnsig50")) nSigmaTPC_K = 5.0;
   if (opt.Contains("KTPCnsig1000")) nSigmaTPC_K = 100.0;

   if (opt.Contains("PTOFnsig")) useTOF_P=kTRUE;
   if (opt.Contains("PTOFacceptUnmatched")) rejectUnmatchedTOF_P=kFALSE;
   if (opt.Contains("PTOFnsig10")) nSigmaTOF_P = 1.0;
   if (opt.Contains("PTOFnsig15")) nSigmaTOF_P = 1.5;
   if (opt.Contains("PTOFnsig20")) nSigmaTOF_P = 2.0;
   if (opt.Contains("PTOFnsig25")) nSigmaTOF_P = 2.5;
   if (opt.Contains("PTOFnsig30")) nSigmaTOF_P = 3.0;
   if (opt.Contains("PTOFnsig40")) nSigmaTOF_P = 4.0;
   if (opt.Contains("PTOFnsig50")) nSigmaTOF_P = 5.0;
   if (opt.Contains("PTOFnsig1000")) nSigmaTOF_P = 100.0;

   if (opt.Contains("KTOFnsig"))  useTOF_K=kTRUE;
   if (opt.Contains("KTOFacceptUnmatched")) rejectUnmatchedTOF_K=kFALSE;
   if (opt.Contains("KTOFnsig10")) nSigmaTOF_K = 1.0;
   if (opt.Contains("KTOFnsig15")) nSigmaTOF_K = 1.5;
   if (opt.Contains("KTOFnsig20")) nSigmaTOF_K = 2.0;
   if (opt.Contains("KTOFnsig25")) nSigmaTOF_K = 2.5;
   if (opt.Contains("KTOFnsig30")) nSigmaTOF_K = 3.0;
   if (opt.Contains("KTOFnsig40")) nSigmaTOF_K = 4.0;
   if (opt.Contains("KTOFnsig50")) nSigmaTOF_K = 5.0;
   if (opt.Contains("KTOFnsig1000")) nSigmaTOF_K = 100.0;


   Bool_t usePDG=kFALSE;
   if (opt.Contains("pdg")) {
      Printf("Using PDG");
      usePDG = kTRUE;
   }

   Bool_t useEta = kFALSE;
   if (opt.Contains("eta")) {
     for(int j=1;j<=9;j++) if(opt.Contains(Form("eta0%i",j))) etaRange=0.1*j;

      Printf("Using ETA range (%.2f,%.2f)",-etaRange,etaRange);
      useEta = kTRUE;
   }

   Bool_t usePPt=kFALSE;
   if(opt.Contains("PPt")){
     Printf("Using Proton pT range (%.2f,%.2f)",PtMin_P,PtMax_P);
     usePPt=kTRUE;
   }

   Bool_t useKPt=kFALSE;
   if(opt.Contains("KPt")){
     Printf("Using Kaon pT range (%.2f,%.2f)",PtMin_K,PtMax_K);
     useKPt=kTRUE;
   }

   Bool_t usePMax_P=kFALSE;
   if(opt.Contains("PPMax")){
     for(int j=1;j<=9;j++) if(opt.Contains(Form("PPMax0%i",j))) PMax_P=0.1*j;
     for(int j=10;j<=30;j++) if(opt.Contains(Form("PPMax%i",j))) PMax_P=0.1*j;
     Printf("Using Proton momentum range (0,%.2f)",PMax_P);
     usePMax_P=kTRUE;
   }

   Bool_t usePMax_K=kFALSE;
   if(opt.Contains("KPMax")){
     for(int j=1;j<=9;j++) if(opt.Contains(Form("KPMax0%i",j))) PMax_K=0.1*j;
     for(int j=10;j<=30;j++) if(opt.Contains(Form("KPMax%i",j))) PMax_K=0.1*j;
     Printf("Using Kaon momentum range (0,%.2f)",PMax_K);
     usePMax_K=kTRUE;
   }

   // PROTON SETTINGS ===========================================

   TString scheme="";
   TString cutname = "p_Lambda";
   if (!opt.IsNull()) cutname += Form("_%s",opt.Data());
   AliRsnCutSet *cutsP = new AliRsnCutSet(cutname.Data(), AliRsnTarget::kDaughter);

   AliRsnCutTrackQuality *qualityCutP = new AliRsnCutTrackQuality("cutQualityP");
   if (useCommonQualityCut>=0) {
      qualityCutP->SetAODTestFilterBit(useCommonQualityCut);
   } else {
      qualityCutP->SetDefaults2010();
   }
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
      if(rejectUnmatchedTOF_P){
	if (!scheme.IsNull()) scheme += "&";
	scheme += cutPTOF->GetName();
      }else{
	AliRsnCutTOFMatch *cutPTOFMatch = new AliRsnCutTOFMatch("cutPTOFMatch");
	cutsP->AddCut(cutPTOFMatch);
	if (!scheme.IsNull()) scheme += "&";
	scheme += Form("(%s|(!%s))",cutPTOF->GetName(),cutPTOFMatch->GetName());
      }
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

   if(usePPt){
     AliRsnValueDaughter *valPtP = new AliRsnValueDaughter(Form("val%sPt%s",AliPID::ParticleName(type1),opt.Data()),AliRsnValueDaughter::kPt);
     AliRsnCutValue *cutPtP = new AliRsnCutValue(Form("cut%sPt%s",AliPID::ParticleName(type1),opt.Data()),PtMin_P,PtMax_P);
     cutPtP->SetTargetType(AliRsnTarget::kDaughter);
     cutPtP->SetValueObj(valPtP);
     cutsP->AddCut(cutPtP);
     if (!scheme.IsNull()) scheme += "&";
     scheme += cutPtP->GetName();
   }

   if(usePMax_P){
     AliRsnValueDaughter *valPP = new AliRsnValueDaughter(Form("val%sP%s",AliPID::ParticleName(type1),opt.Data()),AliRsnValueDaughter::kP);
     AliRsnCutValue *cutPP = new AliRsnCutValue(Form("cut%sP%s",AliPID::ParticleName(type1),opt.Data()),0.,PMax_P);
     cutPP->SetTargetType(AliRsnTarget::kDaughter);
     cutPP->SetValueObj(valPP);
     cutsP->AddCut(cutPP);
     if (!scheme.IsNull()) scheme += "&";
     scheme += cutPP->GetName();
   }

   Printf ("CUT Scheme for PROTON is '%s'",scheme.Data());
   cutsP->SetCutScheme(scheme.Data());

   // END PROTON =======================================

   // KAON SETTINGS =======================================
   scheme="";
   cutname = "K_Lambda";
   if (!opt.IsNull()) cutname += Form("_%s",opt.Data());
   AliRsnCutSet *cutsK = new AliRsnCutSet(cutname.Data(), AliRsnTarget::kDaughter);

   AliRsnCutTrackQuality *qualityCutK = new AliRsnCutTrackQuality("cutQualityK");
   if (useCommonQualityCut>=0) {
      qualityCutK->SetAODTestFilterBit(useCommonQualityCut);
   } else {
      qualityCutK->SetDefaults2010();
   }
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
      if(rejectUnmatchedTOF_K){
	if (!scheme.IsNull()) scheme += "&";
	scheme += cutKTOF->GetName();
      }else{
	AliRsnCutTOFMatch *cutKTOFMatch = new AliRsnCutTOFMatch("cutKTOFMatch");
	cutsK->AddCut(cutKTOFMatch);
	if (!scheme.IsNull()) scheme += "&";
	scheme += Form("(%s|(!%s))",cutKTOF->GetName(),cutKTOFMatch->GetName());
      }
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

   if(useKPt){
     AliRsnValueDaughter *valPtK = new AliRsnValueDaughter(Form("val%sPt%s",AliPID::ParticleName(type2),opt.Data()),AliRsnValueDaughter::kPt);
     AliRsnCutValue *cutPtK = new AliRsnCutValue(Form("cut%sPt%s",AliPID::ParticleName(type2),opt.Data()),PtMin_K,PtMax_K);
     cutPtK->SetTargetType(AliRsnTarget::kDaughter);
     cutPtK->SetValueObj(valPtK);
     cutsK->AddCut(cutPtK);
     if (!scheme.IsNull()) scheme += "&";
     scheme += cutPtK->GetName();
   }

   if(usePMax_K){
     AliRsnValueDaughter *valPK = new AliRsnValueDaughter(Form("val%sP%s",AliPID::ParticleName(type2),opt.Data()),AliRsnValueDaughter::kP);
     AliRsnCutValue *cutPK = new AliRsnCutValue(Form("cut%sP%s",AliPID::ParticleName(type2),opt.Data()),0.,PMax_K);
     cutPK->SetTargetType(AliRsnTarget::kDaughter);
     cutPK->SetValueObj(valPK);
     cutsK->AddCut(cutPK);
     if (!scheme.IsNull()) scheme += "&";
     scheme += cutPK->GetName();
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
