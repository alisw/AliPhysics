#ifndef __CINT__
#include <AliRsnCutPIDNSigma.h>
#include <PWGLF/RESONANCES/AliRsnCutTrackQuality.h>
#endif
Int_t AddRsnDaughterCutsRho(AliPID::EParticleType type1,AliPID::EParticleType type2,TString opt,Bool_t isRsnMini=kFALSE,AliRsnInputHandler *rsnIH=0,AliAnalysisTaskSE *task=0)
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

   if (opt.Contains("TPCsigma1")) nSigmaTPC = 1.0;
   if (opt.Contains("TPCsigma2")) nSigmaTPC = 2.0;
   if (opt.Contains("TPCsigma3")) nSigmaTPC = 3.0;

   if (opt.Contains("TOFsigma1")) nSigmaTOF = 1.0;
   if (opt.Contains("TOFsigma2")) nSigmaTOF = 2.0;
   if (opt.Contains("TOFsigma3")) nSigmaTOF = 3.0;

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

   Bool_t usetof = kFALSE;
   if(opt.Contains("tof")) {
      Printf("Using tof PID range (%.2f,%.2f)",0.0,1E+20);
      usetof = kTRUE;
   }

   Bool_t usetpc = kFALSE;
   if(opt.Contains("tpc")) {
      Printf("Using tpc PID range (%.2f,%.2f)",0.0,1E+20);
      usetpc = kTRUE;
   }

//---------------------------------------------
//  Combine cuts
//---------------------------------------------

   TString cutname = "pionRho";
   if (!opt.IsNull()) cutname += Form("_%s",opt.Data());
   AliRsnCutSet *cuts = new AliRsnCutSet(cutname.Data(), AliRsnTarget::kDaughter);

   TString scheme="";

   AliRsnCutTrackQuality *qualityCut = new AliRsnCutTrackQuality("cutQuatity");
   //qualityCut->SetDefaults2010();
   //qualityCut->SetDCAZmax(0.3);
   //qualityCut->SetDCARmax(0.05);

   qualityCut->SetDCAZmax(0.2);
   qualityCut->SetDCARmax(0.02);
   qualityCut->AddStatusFlag(AliESDtrack::kTPCin   , kTRUE);
   qualityCut->AddStatusFlag(AliESDtrack::kTPCrefit, kTRUE);
   qualityCut->AddStatusFlag(AliESDtrack::kITSrefit, kTRUE);


   qualityCut->SetPtRange(0.15, 1E+20);
   qualityCut->SetEtaRange(-0.8, 0.8);
   qualityCut->SetSPDminNClusters(0);
   qualityCut->SetITSminNClusters(0);
   qualityCut->SetITSmaxChi2(1E+20);
   qualityCut->SetTPCminNClusters(80);
   qualityCut->SetTPCmaxChi2(4.0);
   qualityCut->SetRejectKinkDaughters();
   cuts->AddCut(qualityCut);
   if (!scheme.IsNull()) scheme += "&";
   scheme += qualityCut->GetName();

   if (usetpc) {
      AliRsnCutPIDNSigma *cutPiTPC = new AliRsnCutPIDNSigma("cutPIDNSigmaTPC",AliPID::kPion,AliRsnCutPIDNSigma::kTPC);
      cutPiTPC->SinglePIDRange(nSigmaTPC);
      //cutPiTPC->AddPIDRange(nSigmaTPC,0.0,0.7);
      cuts->AddCut(cutPiTPC);
      if (!scheme.IsNull()) scheme += "&";
      scheme += cutPiTPC->GetName();
   }

   if(usetof) {
      AliRsnCutPIDNSigma *cutPiTOF = new AliRsnCutPIDNSigma("cutPIDNSigmaTOF",AliPID::kPion,AliRsnCutPIDNSigma::kTOF);
      cutPiTOF->SinglePIDRange(nSigmaTOF);
      //cutPiTOF->AddPIDRange(nSigmaTOF,0.7,1e20);
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

