#ifndef __CINT__
#include <PWGLF/RESONANCES/AliRsnCutSet.h>
#include <PWGLF/RESONANCES/AliRsnInputHandler.h>
#include <PWGLF/RESONANCES/AliRsnCutKaonForPhi2010.h>
#include <PWGLF/RESONANCES/AliRsnMiniAnalysisTask.h>
#include <PWGLF/RESONANCES/AliRsnAnalysisTask.h>
#include <PWGLF/RESONANCES/AliRsnValueDaughter.h>
#include <RESONANCES/AliRsnCutPID.h>
#include <RESONANCES/AliRsnCutKaonForPhi2010PP.h>
#include <RESONANCES/AliRsnCutValue.h>
#endif
Int_t AddRsnDaughterCutsPhi2010(AliPID::EParticleType type1,AliPID::EParticleType type2,TString opt,Bool_t isRsnMini=kFALSE,AliRsnInputHandler *rsnIH=0,AliAnalysisTaskSE *task=0)
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

   Printf("AddRsnDaughterCutsPhi2010 Option : %s",opt.Data());

   Bool_t useTrackPtCut = kFALSE;
   Double_t trackPtMin = 0.0;
   Double_t trackPtMax = 100000.0;
   Double_t nSigmaTPC=3.0;
   Double_t nSigmaTOF=3.0;
   Double_t ptTPCMax=0.6;
   Double_t etaRange=0.5;
//     etaRange=0.1;

   if (opt.Contains("trackPt")) {
      useTrackPtCut = kTRUE;
      if (opt.Contains("trackPtMax18")) trackPtMax = 1.8;
      if (opt.Contains("trackPtMax20")) trackPtMax = 2.0;
      if (opt.Contains("trackPtMax25")) trackPtMax = 2.5;
   }

   if (opt.Contains("TPCsigma1")) nSigmaTPC = 1.0;
   if (opt.Contains("TPCsigma2")) nSigmaTPC = 2.0;
   if (opt.Contains("TPCsigma3")) nSigmaTPC = 3.0;
   if (opt.Contains("TPCsigma1.5")) nSigmaTPC = 1.5;
   if (opt.Contains("TPCsigma2.5")) nSigmaTPC = 2.5;

   if (opt.Contains("TOFsigma1")) nSigmaTOF = 1.0;
   if (opt.Contains("TOFsigma2")) nSigmaTOF = 2.0;
   if (opt.Contains("TOFsigma3")) nSigmaTOF = 3.0;



   if (opt.Contains("tpcptMax")) {
      if (opt.Contains("tpcptMax04")) ptTPCMax=0.4;
      if (opt.Contains("tpcptMax05")) ptTPCMax=0.5;
      if (opt.Contains("tpcptMax06")) ptTPCMax=0.6;
      if (opt.Contains("tpcptMax07")) ptTPCMax=0.7;
      if (opt.Contains("tpcptMax08")) ptTPCMax=0.8;
      if (opt.Contains("tpcptMax09")) ptTPCMax=0.9;
      if (opt.Contains("tpcptMax10")) ptTPCMax=1.0;
   }

   AliRsnCutKaonForPhi2010PP *cutPP=0;
   AliRsnCutKaonForPhi2010 *cutPbPb=0;
   AliRsnCut *cut;
   if (usePPCut) {
      Printf("Using AliRsnCutKaonForPhi2010PP ...");
      AliRsnCutKaonForPhi2010PP *cutPP = new AliRsnCutKaonForPhi2010PP("cutKaonPhi2010PP");
      cutPP->SetTPCNSigmaLow(nSigmaTPC);
      cutPP->SetTPCNSigmaHigh(5.0);
      cutPP->SetTPCLimit(ptTPCMax);
      cutPP->SetTOFNSigma(nSigmaTOF);
      cut = cutPP;
   }
   else {
      Printf("Using AliRsnCutKaonForPhi2010 ...");
      AliRsnCutKaonForPhi2010 *cutPbPb = new AliRsnCutKaonForPhi2010("cutKaonPhi2010",nSigmaTPC,nSigmaTOF,ptTPCMax);
      if (opt.Contains("qualityonly")) cutPbPb->SetMode(AliRsnCutKaonForPhi2010::kQuality);
      if (opt.Contains("tofonly")) cutPbPb->SetMode(AliRsnCutKaonForPhi2010::kOnlyTOF);
      if (opt.Contains("tpconly")) cutPbPb->SetMode(AliRsnCutKaonForPhi2010::kOnlyTPC);

      cut = cutPbPb;
   }

   Bool_t usePDG=kFALSE;
   if (opt.Contains("pdg")) {
      Printf("Using PDG");
      usePDG = kTRUE;
      if (opt.Contains("pdgPI")) type1 = AliPID::kPion;
//       type1 = AliPID::kPion;

   }


   Bool_t useEta = kFALSE;
   if (opt.Contains("eta")) {
      Printf("Using ETA range (%.2f,%.2f)",-etaRange,etaRange);
      useEta = kTRUE;
   }

   //---------------------------------------------
   //  Combine cuts
   //---------------------------------------------
   TString cutname = "kaonPhi2010";
   if (!opt.IsNull()) cutname += Form("_%s",opt.Data());
   AliRsnCutSet *cuts = new AliRsnCutSet(cutname.Data(), AliRsnTarget::kDaughter);
   cuts->AddCut(cut);

   TString scheme="";
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

   if (!scheme.IsNull()) scheme += "&";
   scheme += cut->GetName();

   Printf ("------ scheme '%s'",scheme.Data());
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

