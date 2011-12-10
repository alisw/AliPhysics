#ifndef __CINT__
#include <PWG2/RESONANCES/AliRsnCutPrimaryVertex.h>
#include <PWG2/RESONANCES/AliRsnValuePair.h>
#include <PWG2/RESONANCES/AliRsnListOutput.h>
#include <PWG2/RESONANCES/AliRsnPairDef.h>
#include <PWG2/RESONANCES/AliRsnLoopPair.h>
#include <PWG2/RESONANCES/AliRsnAnalysisTask.h>
#include <PWG2/RESONANCES/AliRsnLoopDaughter.h>
#include <PWG2/RESONANCES/AliRsnValueDaughter.h>
#include <PWG2/RESONANCES/AliRsnMiniAnalysisTask.h>
#include <PWG2/RESONANCES/AliRsnCutMiniPair.h>
#include <PWG2/RESONANCES/AliRsnInputHandler.h>
#include <PWG2/RESONANCES/AliRsnMiniMonitor.h>
#include <ANALYSIS/AliAnalysisManager.h>
#include <PWG2/RESONANCES/AliRsnValueEvent.h>
#endif

Bool_t RsnConfig(AliAnalysisTaskSE *task,Bool_t isMC,Bool_t isMixing,AliRsnInputHandler *rsnIH=0,TList *listRsn=0) {

   if (!task) return kFALSE;

   if (!listRsn) {
      return kFALSE;
   }

   Int_t cutIndex = 0;
   Int_t numOfCuts = 0;

   // set commol eventCuts

   gROOT->LoadMacro("AddRsnCommonEventCuts.C");
   AliRsnCutSet *commonEventCuts = AddRsnCommonEventCuts(task);

   gROOT->LoadMacro("AddRsnCommonPairCuts.C");
   AliRsnCutSet *commonPairCuts = AddRsnCommonPairCuts();


   // TODO this is tmp hack
   if (!rsnIH) rsnIH = new AliRsnInputHandler();

   TIter next(listRsn);
   TNamed *rsnObj=0;
   TString rsnName,rsnNameOpt,rsnNameOptFull,rsnCutName,rsnCutOpt,rsnCutNameOptFull;
   while ((rsnObj = (TNamed *)next())) {
      GetOptionFromString(rsnObj->GetName(),rsnName,rsnNameOpt);
      rsnNameOptFull=rsnName; if (!rsnNameOpt.IsNull()) rsnNameOptFull += Form("_%s",rsnNameOpt.Data());

      GetOptionFromString(rsnObj->GetTitle(),rsnCutName,rsnCutOpt);
      rsnCutNameOptFull=rsnCutName; if (!rsnCutOpt.IsNull())rsnCutNameOptFull += Form("_%s",rsnCutOpt.Data());

      gROOT->LoadMacro(Form("AddRsnDaughterCuts%s.C",rsnCutName.Data()));
      gROOT->LoadMacro(Form("AddRsnPairs%s.C",rsnName.Data()));

      rsnNameOptFull.ToLower();
      rsnName.ToLower();
      // add cuts
      if (!rsnName.CompareTo("phi")) {
         numOfCuts = gROOT->ProcessLine(Form("AddRsnDaughterCuts%s(AliPID::kKaon,AliPID::kKaon,\"%s\",%d,(AliRsnInputHandler*)%p,(AliAnalysisTaskSE*)%p)",rsnCutName.Data(), rsnCutOpt.Data(),gRsnUseMiniPackage,rsnIH, task));
         if (numOfCuts) {
            if (rsnNameOpt.Contains("mon")) AddParticleMonitor(task,isMC,cutIndex,commonEventCuts,commonPairCuts,Form("%s_%s_K",rsnNameOptFull.Data(),rsnCutNameOptFull.Data()));
            AddRsnPairsPhi(task,isMC,isMixing,AliPID::kKaon,cutIndex,AliPID::kKaon,cutIndex,commonEventCuts,commonPairCuts,Form("%s_%s",rsnNameOptFull.Data(),rsnCutNameOptFull.Data()));
            cutIndex+=numOfCuts;
         }
      } else if (!rsnName.CompareTo("kstar")) {
         numOfCuts = gROOT->ProcessLine(Form("AddRsnDaughterCuts%s(AliPID::kKaon,AliPID::kPion,\"%s\",%d,(AliRsnInputHandler*)%p,(AliAnalysisTaskSE*)%p)",rsnCutName.Data(), rsnCutOpt.Data(),gRsnUseMiniPackage,rsnIH,task));
         if (numOfCuts) {
            if (rsnNameOpt.Contains("mon")) AddParticleMonitor(task,isMC,cutIndex,commonEventCuts,commonPairCuts,Form("%s_%s_K",rsnNameOptFull.Data(),rsnCutNameOptFull.Data()));
            if (rsnNameOpt.Contains("mon")) AddParticleMonitor(task,isMC,cutIndex+1,commonEventCuts,commonPairCuts,Form("%s_%s_pi",rsnNameOptFull.Data(),rsnCutNameOptFull.Data()));
            AddRsnPairsKStar(task,isMC,isMixing,AliPID::kKaon,cutIndex,AliPID::kPion,cutIndex+1,commonEventCuts,commonPairCuts,Form("%s_%s",rsnNameOptFull.Data(),rsnCutNameOptFull.Data()));
            cutIndex+=numOfCuts;
         }
      } else {
         Printf("Error : Particle %s is not supported !!!!",rsnName.Data());
         return kFALSE;
      }

   }

   return kTRUE;
}

void GetOptionFromString(TString str,TString &outStr1,TString &outStr2,TString d=":") {
//   TStringO
   TObjArray *tokens = str.Tokenize(d.Data());
   TObjString *objStr =  (TObjString *)tokens->At(0);
   if (!objStr) {
      outStr1 = "err";
      outStr2 = "";
      return;
   }
   outStr1 = objStr->GetString();

   objStr =  (TObjString *) tokens->At(1);
   if (!objStr) {
      outStr2 = "";
      return;
   }
   outStr2 = objStr->GetString();

}

Bool_t AddPair(AliAnalysisTaskSE *task, Bool_t isMC,Bool_t isMixing, AliPID::EParticleType pType1,Int_t listID1, AliPID::EParticleType pType2,Int_t listID2, Int_t pdgMother,Double_t massMother, AliRsnCutSet *commonEventCuts=0,AliRsnCutSet *commonPairCuts=0, TString name = "") {

   Bool_t typeSame = (pType1 == pType2);

//    Printf("------------- id1=%d id2=%d",listID1,listID2);

   TList *listLoops = new TList;

   // pair definition
   AliRsnPairDef  *pairDefPM         = new AliRsnPairDef((AliRsnDaughter::ESpecies)pType1, '+', (AliRsnDaughter::ESpecies)pType2, '-', pdgMother, massMother);
   AliRsnPairDef  *pairDefMP         = new AliRsnPairDef((AliRsnDaughter::ESpecies)pType1, '-', (AliRsnDaughter::ESpecies)pType2, '+', pdgMother, massMother);
   AliRsnPairDef  *pairDefPP         = new AliRsnPairDef((AliRsnDaughter::ESpecies)pType1, '+', (AliRsnDaughter::ESpecies)pType2, '+', pdgMother, massMother);
   AliRsnPairDef  *pairDefMM         = new AliRsnPairDef((AliRsnDaughter::ESpecies)pType1, '-', (AliRsnDaughter::ESpecies)pType2, '-', pdgMother, massMother);

   // loop object creation
   AliRsnLoopPair *lp = 0;

   // sets +-
   lp = new AliRsnLoopPair(Form("%s_PM", name.Data()), pairDefPM, kFALSE);
   listLoops->Add(lp);

   if (!typeSame) {
      lp = new AliRsnLoopPair(Form("%s_MP", name.Data()), pairDefMP, kFALSE);
      listLoops->Add(lp);
   }

   // sets +- TRUE pairs
   if (isMC) {
      lp = new AliRsnLoopPair(Form("%s_PM_TRUE", name.Data()), pairDefPM, kFALSE);
      lp->SetOnlyTrue(kTRUE);
      lp->SetCheckDecay(kTRUE);
      listLoops->Add(lp);
      if (!typeSame) {
         lp = new AliRsnLoopPair(Form("%s_MP_TRUE", name.Data()), pairDefMP, kFALSE);
         lp->SetOnlyTrue(kTRUE);
         lp->SetCheckDecay(kTRUE);
         listLoops->Add(lp);
      }
      // sets +- TRUE paris (MC is used for momentum)
      lp = new AliRsnLoopPair(Form("%s_PM_TRUE_MC", name.Data()), pairDefPM, kFALSE);
      lp->SetTrueMC(kTRUE);
      listLoops->Add(lp);
      if (!typeSame) {
         // sets +- TRUE paris (MC is used for momentum)
         lp = new AliRsnLoopPair(Form("%s_MP_TRUE_MC", name.Data()), pairDefMP, kFALSE);
         lp->SetTrueMC(kTRUE);
         listLoops->Add(lp);
      }
   }

   // sets ++
   lp = new AliRsnLoopPair(Form("%s_PP", name.Data()), pairDefPP, kFALSE);
   listLoops->Add(lp);

   // sets --
   lp = new AliRsnLoopPair(Form("%s_MM", name.Data()), pairDefMM, kFALSE);
   listLoops->Add(lp);

   if (isMixing) {
      // sets +- Mixing (NOT mini)
      lp = new AliRsnLoopPair(Form("%s_PM", name.Data()), pairDefPM, kTRUE);
      listLoops->Add(lp);

      // sets -+ Mixing (NOT mini)
      lp = new AliRsnLoopPair(Form("%s_MP", name.Data()), pairDefMP, kTRUE);
      listLoops->Add(lp);

      // sets ++ Mixing (NOT mini)
      lp = new AliRsnLoopPair(Form("%s_PP", name.Data()), pairDefPP, kTRUE);
      listLoops->Add(lp);

      // sets -- Mixing (NOT mini)
      lp = new AliRsnLoopPair(Form("%s_MM", name.Data()), pairDefMM, kTRUE);
      listLoops->Add(lp);
   }


   // loops over all AliRsnLoops and sets everything (don't touch it if you don't know what you are doing)
   TIter next(listLoops);
   while ((lp = (AliRsnLoopPair *)next.Next())) {
      lp->SetListID(0, listID1);
      lp->SetListID(1, listID2);
      lp->SetMCRefInfo(gRsnUseMCMomentum);
      if (commonPairCuts) lp->SetPairCuts(commonPairCuts);
      if (commonEventCuts) lp->SetEventCuts(commonEventCuts);
      if (name.Contains("phi")) AddPairOutputPhi(lp);
      else if (name.Contains("kstar")) AddPairOutputKStar(lp);
      else continue;
      ((AliRsnAnalysisTask *)task)->AddLoop(lp);
   }
   return kTRUE;
}

// void AddMonitorOutput(AliRsnLoopDaughter *mon)

void AddMonitorOutput(TObjArray *mon)
{
// mcinfo is not supported yet
   if (gRsnUseMCMomentum) return ;

   // dEdx tpc
   AliRsnValueDaughter *axisMomTPC = new AliRsnValueDaughter("pTPC", AliRsnValueDaughter::kPtpc);
   AliRsnValueDaughter *axisSigTPC = new AliRsnValueDaughter("sTPC", AliRsnValueDaughter::kTPCsignal);
   axisMomTPC->SetBins(0.0,5.0,0.01);
   axisSigTPC->SetBins(0.0,500.0,2.0);

   // output: 2D histogram of TPC signal vs. TPC momentum
   AliRsnListOutput *outMonitordEdxTPC = new AliRsnListOutput("dEdx", AliRsnListOutput::kHistoDefault);
   outMonitordEdxTPC->AddValue(axisMomTPC);
   outMonitordEdxTPC->AddValue(axisSigTPC);

   // add outputs to loop
   mon->Add(outMonitordEdxTPC);

   // dEdx tpc
   AliRsnValueDaughter *axisMomTPCForTOF = new AliRsnValueDaughter("pTPC", AliRsnValueDaughter::kPtpc);
   AliRsnValueDaughter *axisSigTOF = new AliRsnValueDaughter("sTOF", AliRsnValueDaughter::kTOFsignal);
   axisMomTPCForTOF->SetBins(0.0,5.0,0.01);
   axisSigTOF->SetBins(0.0,500.0,2.0);

   // output: 2D histogram of TPC signal vs. TPC momentum
   AliRsnListOutput *outMonitordEdxTOF = new AliRsnListOutput("Edx", AliRsnListOutput::kHistoDefault);
   outMonitordEdxTOF->AddValue(axisMomTPCForTOF);
   outMonitordEdxTOF->AddValue(axisSigTOF);

   // add outputs to loop
//    mon->Add(outMonitordEdxTOF);


   // Momentum
   AliRsnValueDaughter *axisMomP = new AliRsnValueDaughter("p", AliRsnValueDaughter::kP);
   axisMomP->SetBins(0.0,5.0,0.01);

   // output: 2D histogram of TPC signal vs. TPC momentum
   AliRsnListOutput *outMonitorP = new AliRsnListOutput("P", AliRsnListOutput::kHistoDefault);
   outMonitorP->AddValue(axisMomP);

   // add outputs to loop
   mon->Add(outMonitorP);

   // Momentum Pt
   AliRsnValueDaughter *axisMomPt = new AliRsnValueDaughter("pt", AliRsnValueDaughter::kPt);
   axisMomPt->SetBins(0.0,5.0,0.01);

   // output: 2D histogram of TPC signal vs. TPC momentum
   AliRsnListOutput *outMonitorPt = new AliRsnListOutput("Pt", AliRsnListOutput::kHistoDefault);
   outMonitorPt->AddValue(axisMomPt);

   // add outputs to loop
   mon->Add(outMonitorPt);

   // Eta
   AliRsnValueDaughter *axisMomEta = new AliRsnValueDaughter("eta", AliRsnValueDaughter::kEta);
   axisMomEta->SetBins(-1.0,1.0,0.01);

   // output: 2D histogram of TPC signal vs. TPC momentum
   AliRsnListOutput *outMonitorEta = new AliRsnListOutput("Eta", AliRsnListOutput::kHistoDefault);
   outMonitorEta->AddValue(axisMomEta);

   // add outputs to loop
   mon->Add(outMonitorEta);

   // kTOFnsigmaK
   AliRsnValueDaughter *axisTPCnsigmaK = new AliRsnValueDaughter("K", AliRsnValueDaughter::kTPCnsigmaK);
   axisTPCnsigmaK->SetBins(1000,0,100);

   // output: 2D histogram of TPC signal vs. TPC momentum
   AliRsnListOutput *outMonitorTPCnsigmaK = new AliRsnListOutput("TPC_nsigma", AliRsnListOutput::kHistoDefault);
   outMonitorTPCnsigmaK->AddValue(axisTPCnsigmaK);

   // add outputs to loop
   mon->Add(outMonitorTPCnsigmaK);

   // kTPCnsigmaPi
   AliRsnValueDaughter *axisTPCnsigmaPi = new AliRsnValueDaughter("pi", AliRsnValueDaughter::kTPCnsigmaPi);
   axisTPCnsigmaPi->SetBins(1000,0,100);

   // output: 2D histogram of TPC signal vs. TPC momentum
   AliRsnListOutput *outMonitorTPCnsigmaPi = new AliRsnListOutput("TPC_nsigma", AliRsnListOutput::kHistoDefault);
   outMonitorTPCnsigmaPi->AddValue(axisTPCnsigmaPi);

   // add outputs to loop
   mon->Add(outMonitorTPCnsigmaPi);

   // kTPCnsigmaP
   AliRsnValueDaughter *axisTPCnsigmaP = new AliRsnValueDaughter("p", AliRsnValueDaughter::kTPCnsigmaP);
   axisTPCnsigmaP->SetBins(1000,0,100);

   // output: 2D histogram of TPC signal vs. TPC momentum
   AliRsnListOutput *outMonitorTPCnsigmaP = new AliRsnListOutput("TPC_nsigma", AliRsnListOutput::kHistoDefault);
   outMonitorTPCnsigmaP->AddValue(axisTPCnsigmaP);

   // add outputs to loop
   mon->Add(outMonitorTPCnsigmaP);


   // kTOFnsigmaK
   AliRsnValueDaughter *axisTOFnsigmaK = new AliRsnValueDaughter("K", AliRsnValueDaughter::kTOFnsigmaK);
   axisTOFnsigmaK->SetBins(1000,0,100);

   // output: 2D histogram of TPC signal vs. TPC momentum
   AliRsnListOutput *outMonitorTOFnsigmaK = new AliRsnListOutput("TOF_nsigma", AliRsnListOutput::kHistoDefault);
   outMonitorTOFnsigmaK->AddValue(axisTOFnsigmaK);

   // add outputs to loop
   mon->Add(outMonitorTOFnsigmaK);

   // kTOFnsigmaPi
   AliRsnValueDaughter *axisTOFnsigmaPi = new AliRsnValueDaughter("pi", AliRsnValueDaughter::kTOFnsigmaPi);
   axisTOFnsigmaPi->SetBins(1000,0,100);

   // output: 2D histogram of TPC signal vs. TPC momentum
   AliRsnListOutput *outMonitorTOFnsigmaPi = new AliRsnListOutput("TOF_nsigma", AliRsnListOutput::kHistoDefault);
   outMonitorTOFnsigmaPi->AddValue(axisTOFnsigmaPi);

   // add outputs to loop
   mon->Add(outMonitorTOFnsigmaPi);

   // kTOFnsigmaP
   AliRsnValueDaughter *axisTOFnsigmaP = new AliRsnValueDaughter("p", AliRsnValueDaughter::kTOFnsigmaP);
   axisTOFnsigmaP->SetBins(1000,0,100);

   // output: 2D histogram of TPC signal vs. TPC momentum
   AliRsnListOutput *outMonitorTOFnsigmaP = new AliRsnListOutput("TOF_nsigma", AliRsnListOutput::kHistoDefault);
   outMonitorTOFnsigmaP->AddValue(axisTOFnsigmaP);

   // add outputs to loop
   mon->Add(outMonitorTOFnsigmaP);


   AliRsnListOutput *outMonitorPTvsMult = new AliRsnListOutput("PTvsMult",AliRsnListOutput::kHistoDefault);
   AliRsnValueDaughter *vd1 = new AliRsnValueDaughter("pt",AliRsnValueDaughter::kPt);
   vd1->SetBins(0.0,5.0,0.01);
   outMonitorPTvsMult->AddValue(vd1);

   AliRsnValueEvent *ve1 = new AliRsnValueEvent("mult",AliRsnValueEvent::kMult);
   ve1->SetBins(0.0,100.0,1);
   outMonitorPTvsMult->AddValue(ve1);
   mon->Add(outMonitorPTvsMult);

//    mon->SetMCRefInfo(gRsnUseMCMomentum);

}

void AddMonitorOutputMini(AliRsnMiniAnalysisTask *task,Int_t listID1,TString name = "",Char_t charge='0')
{
   TString chargeName="all";
   if ( charge == '+' ) chargeName = "pos";
   if ( charge == '-' ) chargeName = "neg";

   AliRsnMiniMonitor *mondEdx = new AliRsnMiniMonitor(Form("%s_dEdxTPCvsP_%s", name.Data(),chargeName.Data()),AliRsnMiniMonitor::kdEdxTPCvsP, listID1);
   mondEdx->SetCharge(charge);

   AliRsnMiniMonitor *monPt = new AliRsnMiniMonitor(Form("%s_Pt_%s", name.Data(),chargeName.Data()),AliRsnMiniMonitor::kTrackPt, listID1);
   monPt->SetCharge(charge);

}

void AddParticleMonitor(AliAnalysisTaskSE *task, Bool_t isMC, Int_t listID1,AliRsnCutSet *commonEventCuts=0,AliRsnCutSet *cutPair=0,TString name = "")
{

   if (gRsnUseMiniPackage) {
      Printf("Monitoring by mini is not supported now. It will be soon !!!");
      return ;
//       AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
//
//       AddMonitorOutputMini((AliRsnMiniAnalysisTask *)task,listID1,name);
//       AddMonitorOutputMini((AliRsnMiniAnalysisTask *)task,listID1,name,'+');
//       AddMonitorOutputMini((AliRsnMiniAnalysisTask *)task,listID1,name,'-');
   } else {
      TList *listLoops = new TList;
      // monitor definition
      AliRsnDaughterDef *tracksAll = new AliRsnDaughterDef(AliRsnDaughter::kTrack /*'+' or '-'*/);
      AliRsnDaughterDef *tracksPos = new AliRsnDaughterDef(AliRsnDaughter::kTrack,'+');
      AliRsnDaughterDef *tracksNeg = new AliRsnDaughterDef(AliRsnDaughter::kTrack,'-');

      AliRsnLoopDaughter *lm =0;
      // loop object
      listLoops->Add(new AliRsnLoopDaughter(Form("%s_all", name.Data()), listID1, tracksAll));
      listLoops->Add(new AliRsnLoopDaughter(Form("%s_pos", name.Data()), listID1, tracksPos));
      listLoops->Add(new AliRsnLoopDaughter(Form("%s_neg", name.Data()), listID1, tracksNeg));

      TIter next(listLoops);
      while ((lm = (AliRsnLoopDaughter *)next.Next()))) {
         if (commonEventCuts) lm->SetEventCuts(commonEventCuts);
         AddMonitorOutput(lm);
         ((AliRsnAnalysisTask *)task)->AddLoop(lm);
      }

   }

}



