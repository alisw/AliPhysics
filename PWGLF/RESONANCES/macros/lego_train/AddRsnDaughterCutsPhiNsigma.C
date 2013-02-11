#ifndef __CINT__
#include <AliRsnCutTrackQuality.h>
#endif
Int_t AddRsnDaughterCutsPhiNsigma(AliPID::EParticleType type1,AliPID::EParticleType type2,TString opt,AliRsnInputHandler *rsnIH=0,AliAnalysisTaskSE *task=0)
{

   if (!rsnIH) return 0;

   Bool_t valid = kTRUE;
//   Int_t collisionType = AliRsnTrainManager::GetGlobalInt("IsCollisionType",valid);
   Int_t useCommonQualityCut = AliRsnTrainManager::GetGlobalInt("RsnCommonQualityCut",valid);
   TString rsnQualityCut = AliRsnTrainManager::GetGlobalStr("RsnQualityCut",valid);
   Int_t isMC = AliRsnTrainManager::GetGlobalInt("IsMC",valid);
   Int_t isRsnMini = AliRsnTrainManager::GetGlobalInt("IsRsnMini",valid);
   Int_t isMixing = AliRsnTrainManager::GetGlobalInt("IsMixing",valid);

   // experts only (don't touch)
   Int_t isRsnDev = AliAnalysisManager::GetGlobalInt("rsnUseRSNParDev",valid);

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
   Int_t NclTPC=70;
   Char_t DCAxyFormula[100]="0.0182+0.035/pt^1.01";

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

   if (opt.Contains("KTPCnsig05")) nSigmaTPC = 0.5;
   if (opt.Contains("KTPCnsig08")) nSigmaTPC = 0.8;
   if (opt.Contains("KTPCnsig10")) nSigmaTPC = 1.0;
   if (opt.Contains("KTPCnsig15")) nSigmaTPC = 1.5;
   if (opt.Contains("KTPCnsig20")) nSigmaTPC = 2.0;
   if (opt.Contains("KTPCnsig25")) nSigmaTPC = 2.5;
   if (opt.Contains("KTPCnsig30")) nSigmaTPC = 3.0;
   if (opt.Contains("KTPCnsig40")) nSigmaTPC = 4.0;
   if (opt.Contains("KTPCnsig50")) nSigmaTPC = 5.0;
   if (opt.Contains("KTPCnsig1000")) nSigmaTPC = 100.0;

   if (opt.Contains("KTOFnsig10")) nSigmaTOF = 1.0;
   if (opt.Contains("KTOFnsig15")) nSigmaTOF = 1.5;
   if (opt.Contains("KTOFnsig20")) nSigmaTOF = 2.0;
   if (opt.Contains("KTOFnsig25")) nSigmaTOF = 2.5;
   if (opt.Contains("KTOFnsig30")) nSigmaTOF = 3.0;
   if (opt.Contains("KTOFnsig1000")) nSigmaTOF = 100.0;

   if (opt.Contains("trackPt")) {
      useTrackPtCut = kTRUE;
      if (opt.Contains("trackPtMin015")) trackPtMin = 0.15;
      if (opt.Contains("trackPtMin02")) trackPtMin = 0.2;
      if (opt.Contains("trackPtMin05")) trackPtMin = 0.5;
      if (opt.Contains("trackPtMin06")) trackPtMin = 0.6;

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
      if(opt.Contains("eta09")) etaRange=0.9;
      if(opt.Contains("eta08")) etaRange=0.8;
      if(opt.Contains("eta07")) etaRange=0.7;
      if(opt.Contains("eta06")) etaRange=0.6;
      Printf("Using ETA range (%.2f,%.2f)",-etaRange,etaRange);
      useEta = kTRUE;
   }

   Bool_t useNclTPC = kFALSE;
   if (opt.Contains("NclTPC")) {
      if (opt.Contains("NclTPC70")) NclTPC=70;
      if (opt.Contains("NclTPC75")) NclTPC=75;
      if (opt.Contains("NclTPC80")) NclTPC=80;
      if (opt.Contains("NclTPC85")) NclTPC=85;
      if (opt.Contains("NclTPC90")) NclTPC=90;
      useNclTPC = kTRUE;
   }

   Bool_t useDCAxy = kFALSE;
   if (opt.Contains("DCAxy")) {
      if (opt.Contains("DCAxyFormula7s")) sprintf(DCAxyFormula,"0.0182+0.035/pt^1.01");
      if (opt.Contains("DCAxyFormula6s")) sprintf(DCAxyFormula,"0.0156+0.03/pt^1.01");
      if (opt.Contains("DCAxyFormula5s")) sprintf(DCAxyFormula,"0.013+0.025/pt^1.01");
      useDCAxy = kTRUE;
   }

//---------------------------------------------
//  Combine cuts
//---------------------------------------------

   TString cutname = "K_Phi";
   if (!opt.IsNull()) cutname += Form("_%s",opt.Data());
   AliRsnCutSet *cuts = new AliRsnCutSet(cutname.Data(), AliRsnTarget::kDaughter);

   TString scheme="";
   AliRsnCutTrackQuality *qualityCut = new AliRsnCutTrackQuality("cutQualityK");
   if (!rsnQualityCut.IsNull()) {
      AliESDtrackCuts *esdTK = RsnQualityCut(rsnQualityCut.Data());
      if(useDCAxy) esdTK->SetMaxDCAToVertexXYPtDep(DCAxyFormula);
      qualityCut->SetESDtrackCuts(esdTK);
   } else {
      if (useCommonQualityCut>=0) {
         qualityCut->SetAODTestFilterBit(useCommonQualityCut);
         if(useDCAxy) {qualityCut->SetCheckOnlyFilterBit(kFALSE); qualityCut->SetDCARPtFormula(DCAxyFormula);}
      } else {
         qualityCut->SetDefaults2010();
         if(useDCAxy) qualityCut->SetDCARPtFormula(DCAxyFormula);
      }
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

   if (useNclTPC) {
      Printf("Adding NclTPC >= %i",NclTPC);
      AliRsnValueDaughter *valNclTPC = new AliRsnValueDaughter(Form("val%sNclTPC%s",AliPID::ParticleName(type1),opt.Data()),AliRsnValueDaughter::kNTPCclusters);
      AliRsnCutValue *cutNclTPC = new AliRsnCutValue(Form("cut%sNclTPC%s",AliPID::ParticleName(type1),opt.Data()),NclTPC-0.1,1000.);
      cutNclTPC->SetTargetType(AliRsnTarget::kDaughter);
      cutNclTPC->SetValueObj(valNclTPC);
      cuts->AddCut(cutNclTPC);
      if (!scheme.IsNull()) scheme += "&";
      scheme += cutNclTPC->GetName();
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
      if (isRsnDev>=0 && opt.Contains("pairPID")) {
         AliRsnActPostDaughterSelection *pairPID = new AliRsnActPostDaughterSelection();
         pairPID->SetID(0);

         const char *fn="rsnRange.txt";
         if (!gSystem->AccessPathName(fn)) {
            TString minStr = gSystem->GetFromPipe(TString::Format("head -n 1 %s").Data());
            TString maxStr = gSystem->GetFromPipe(TString::Format("tail -n 1 %s").Data());
            pairPID->SetMass(minStr.Atof(),maxStr.Atof());
         } else {
            //       pairPID->SetMass(1.01,1.03);
            pairPID->SetMass(1.015,1.025);
            pairPID->SetMass(1.019,1.021);
            pairPID->SetMass(1.0195,1.0205);
            pairPID->SetMass(1.1000,1.1005);
            //       pairPID->SetMass(1.1005,1.1010);
         }
         sel->AddAction(pairPID);
      }

   }
   return numberOfCuts;

}

