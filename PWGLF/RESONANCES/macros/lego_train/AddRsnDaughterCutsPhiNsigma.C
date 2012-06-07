#ifndef __CINT__
#include <AliRsnCutTrackQuality.h>
#endif
Int_t AddRsnDaughterCutsPhiNsigma(AliPID::EParticleType type1,AliPID::EParticleType type2,TString opt,Bool_t isRsnMini=kFALSE,AliRsnInputHandler *rsnIH=0,AliAnalysisTaskSE *task=0)
{

   if (!rsnIH) return 0;

   Bool_t valid = kTRUE;
   Int_t isPP = AliAnalysisManager::GetGlobalInt("rsnIsPP",valid);
   Int_t useCommonQualityCut = AliAnalysisManager::GetGlobalInt("rsnCommonQualityCut",valid);


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
   Double_t trackPtMin=0.;
   Double_t trackPtMax=1.e10;

   Bool_t useTPC_K=kFALSE;
   Bool_t useTOF_K=kFALSE;
   Bool_t useTrackPtCut=kFALSE;

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

   if (opt.Contains("trackPt")) {
      useTrackPtCut = kTRUE;
      if (opt.Contains("trackPtMin02")) trackPtMin = 0.2;
      if (opt.Contains("trackPtMin05")) trackPtMin = 0.5;

      if (opt.Contains("trackPtMax18")) trackPtMax = 1.8;
      if (opt.Contains("trackPtMax20")) trackPtMax = 2.0;
      if (opt.Contains("trackPtMax25")) trackPtMax = 2.5;
   }

   Bool_t usePDG=kFALSE;
   if (opt.Contains("pdg")) {
      Printf("Using PDG");
      usePDG = kTRUE;
   }

   Bool_t useEta = kFALSE;
   if (opt.Contains("eta")) {
     if(opt.Contains("eta08")) etaRange=0.8;
     if(opt.Contains("eta07")) etaRange=0.7;
     if(opt.Contains("eta06")) etaRange=0.6;
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
   AliRsnCutTrackQuality *qualityCut = new AliRsnCutTrackQuality("cutQualityK");
   if (useCommonQualityCut>=0) {
      qualityCut->SetAODTestFilterBit(useCommonQualityCut);
   } else {
      qualityCut->SetDefaults2010();
   }
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

   if (useTrackPtCut) {
      Printf("Adding Pt min=%.3f max=%.3f ...",trackPtMin,trackPtMax);
      AliRsnValueDaughter *valTrackPt = new AliRsnValueDaughter(Form("val%sTrackPt%s",AliPID::ParticleName(type1),opt.Data()),AliRsnValueDaughter::kPt);

      AliRsnCutValue *cutTrackPt = new AliRsnCutValue(Form("cut%sTrackPt%s",AliPID::ParticleName(type1),opt.Data()),trackPtMin,trackPtMax);
      cutTrackPt->SetTargetType(AliRsnTarget::kDaughter);
      cutTrackPt->SetValueObj(valTrackPt);
      cuts->AddCut(cutTrackPt);
      if (!scheme.IsNull()) scheme += "&";
      scheme += cutTrackPt->GetName();
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
      AliRsnMiniAnalysisTask *taskRsnMini = (AliRsnMiniAnalysisTask *)task;
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

