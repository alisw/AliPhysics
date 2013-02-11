#ifndef __CINT__
#include <AliAnalysisManager.h>
#include <PWGLF/RESONANCES/AliRsnCutPrimaryVertex.h>
#include <PWGLF/RESONANCES/AliRsnValuePair.h>
#include <PWGLF/RESONANCES/AliRsnListOutput.h>
#include <PWGLF/RESONANCES/AliRsnPairDef.h>
#include <PWGLF/RESONANCES/AliRsnLoopPair.h>
#include <PWGLF/RESONANCES/AliRsnAnalysisTask.h>
#include <PWGLF/RESONANCES/AliRsnLoopDaughter.h>
#include <PWGLF/RESONANCES/AliRsnValueDaughter.h>
#include <PWGLF/RESONANCES/AliRsnMiniAnalysisTask.h>
#include <PWGLF/RESONANCES/AliRsnCutMiniPair.h>
#include <PWGLF/RESONANCES/AliRsnInputHandler.h>
#include <PWGLF/RESONANCES/AliRsnMiniMonitor.h>
#include <PWGLF/RESONANCES/AliRsnValueEvent.h>
#include <AliRsnMiniMonitorTask.h>
#include <TROOT.h>
#endif

Bool_t RsnConfig(AliAnalysisTaskSE *task,AliRsnInputHandler *rsnIH=0,TList *listRsn=0) {
// ,Bool_t isMC,Bool_t isMixing,AliRsnInputHandler *rsnIH=0,TList *listRsn=0) {

   if (!task) return kFALSE;

   if (!listRsn) {
      return kFALSE;
   }
   Bool_t valid;
   Int_t isMC = AliRsnTrainManager::GetGlobalInt("IsMC",valid);
   Int_t isRsnMini = AliRsnTrainManager::GetGlobalInt("IsRsnMini",valid);
   Int_t isMixing = AliRsnTrainManager::GetGlobalInt("IsMixing",valid);
   // TODO
   TString rsnCutOptCommon = AliRsnTrainManager::GetGlobalStr("RsnLegoTrainCommonCutOption",valid);


   Int_t cutIndex = 0;
   Int_t numOfCuts = 0;

   // set commol eventCuts

   if (!RsnLoadMacroFromConfig("AddRsnCommonEventCuts.C")) return kFALSE;
   AliRsnCutSet *commonEventCuts = AddRsnCommonEventCuts(task);

   if (!RsnLoadMacroFromConfig("AddRsnCommonPairCuts.C")) return kFALSE;
   AliRsnCutSet *commonPairCuts = AddRsnCommonPairCuts();

   TIter next(listRsn);
   TNamed *rsnObj=0;
   TString rsnName,rsnNameOpt,rsnNameOptFull,rsnCutName,rsnCutOpt,rsnCutNameOptFull;
   while ((rsnObj = (TNamed *)next())) {
      GetOptionFromString(rsnObj->GetName(),rsnName,rsnNameOpt);
      rsnNameOptFull=rsnName;
      if (!rsnNameOpt.IsNull()) rsnNameOptFull += Form("_%s",rsnNameOpt.Data());

      GetOptionFromString(rsnObj->GetTitle(),rsnCutName,rsnCutOpt);
      rsnCutNameOptFull=rsnCutName;
      if (!rsnCutOpt.IsNull())rsnCutNameOptFull += Form("_%s",rsnCutOpt.Data());

      if (!rsnCutOptCommon.IsNull()) {
         if (!rsnCutOpt.IsNull()) rsnCutOpt += "_";
         rsnCutOpt += rsnCutOptCommon.Data();
      }

      if (!RsnLoadMacroFromConfig(Form("AddRsnDaughterCuts%s.C",rsnCutName.Data()))) return kFALSE;
      if (!RsnLoadMacroFromConfig(Form("AddRsnPairs%s.C",rsnName.Data()))) return kFALSE;

      rsnNameOptFull.ToLower();
      rsnName.ToLower();
      // add cuts
      if (!rsnName.CompareTo("phi")) {
         numOfCuts = gROOT->ProcessLine(Form("AddRsnDaughterCuts%s(AliPID::kKaon,AliPID::kKaon,\"%s\",(AliRsnInputHandler*)%p,(AliAnalysisTaskSE*)%p)",rsnCutName.Data(), rsnCutOpt.Data(),rsnIH, task));
         if (numOfCuts) {
            if (rsnNameOpt.Contains("mon")) AddParticleMonitor(task,isMC,cutIndex,commonEventCuts,commonPairCuts,Form("%s_%s_K",rsnNameOptFull.Data(),rsnCutNameOptFull.Data()));
            AddRsnPairsPhi(task,isMC,isMixing,AliPID::kKaon,cutIndex,AliPID::kKaon,cutIndex,commonEventCuts,commonPairCuts,Form("%s.%s",rsnNameOptFull.Data(),rsnCutNameOptFull.Data()));
            cutIndex+=numOfCuts;
         }
      } else if (!rsnName.CompareTo("kstar")) {
         numOfCuts = gROOT->ProcessLine(Form("AddRsnDaughterCuts%s(AliPID::kKaon,AliPID::kPion,\"%s\",(AliRsnInputHandler*)%p,(AliAnalysisTaskSE*)%p)",rsnCutName.Data(), rsnCutOpt.Data(),rsnIH,task));
         if (numOfCuts) {
            if (rsnNameOpt.Contains("mon")) AddParticleMonitor(task,isMC,cutIndex,commonEventCuts,commonPairCuts,Form("%s_%s_K",rsnNameOptFull.Data(),rsnCutNameOptFull.Data()));
            if (rsnNameOpt.Contains("mon")) AddParticleMonitor(task,isMC,cutIndex+1,commonEventCuts,commonPairCuts,Form("%s_%s_pi",rsnNameOptFull.Data(),rsnCutNameOptFull.Data()));
            AddRsnPairsKStar(task,isMC,isMixing,AliPID::kKaon,cutIndex,AliPID::kPion,cutIndex+1,commonEventCuts,commonPairCuts,Form("%s.%s",rsnNameOptFull.Data(),rsnCutNameOptFull.Data()));
            cutIndex+=numOfCuts;
         }
      } else if (!rsnName.CompareTo("rho")) {
         numOfCuts = gROOT->ProcessLine(Form("AddRsnDaughterCuts%s(AliPID::kPion,AliPID::kPion,\"%s\",(AliRsnInputHandler*)%p,(AliAnalysisTaskSE*)%p)",rsnCutName.Data(), rsnCutOpt.Data(),rsnIH,task));
         if (numOfCuts) {
            if (rsnNameOpt.Contains("mon")) AddParticleMonitor(task,isMC,cutIndex,commonEventCuts,commonPairCuts,Form("%s_%s_pi",rsnNameOptFull.Data(),rsnCutNameOptFull.Data()));
            AddRsnPairsRho(task,isMC,isMixing,AliPID::kPion,cutIndex,AliPID::kPion,cutIndex,commonEventCuts,commonPairCuts,Form("%s.%s",rsnNameOptFull.Data(),rsnCutNameOptFull.Data()));
            cutIndex+=numOfCuts;
         }
      } else if (!rsnName.CompareTo("lambda")) {
         numOfCuts = gROOT->ProcessLine(Form("AddRsnDaughterCuts%s(AliPID::kProton,AliPID::kKaon,\"%s\",(AliRsnInputHandler*)%p,(AliAnalysisTaskSE*)%p)",rsnCutName.Data(), rsnCutOpt.Data(),rsnIH,task));
         if (numOfCuts) {
            if (rsnNameOpt.Contains("mon")) AddParticleMonitor(task,isMC,cutIndex,commonEventCuts,commonPairCuts,Form("%s_%s_p",rsnNameOptFull.Data(),rsnCutNameOptFull.Data()));
            if (rsnNameOpt.Contains("mon")) AddParticleMonitor(task,isMC,cutIndex+1,commonEventCuts,commonPairCuts,Form("%s_%s_K",rsnNameOptFull.Data(),rsnCutNameOptFull.Data()));
            AddRsnPairsLambda(task,isMC,isMixing,AliPID::kProton,cutIndex,AliPID::kKaon,cutIndex+1,commonEventCuts,commonPairCuts,Form("%s.%s",rsnNameOptFull.Data(),rsnCutNameOptFull.Data()));
            cutIndex+=numOfCuts;
         }
      } else {
         Printf("Error : Particle %s is not supported !!!!",rsnName.Data());
         return kFALSE;
      }

   }

   return kTRUE;
}

Bool_t RsnLoadMacroFromConfig(TString macro,TString path="") {

   Bool_t valid;
   TString lego_path = AliAnalysisManager::GetGlobalStr("RsnLegoTrainPath",valid);
   if (!valid) lego_path = "$ALICE_ROOT/PWGLF/RESONANCES/macros/lego_train";

   if (!gSystem->AccessPathName(macro.Data())) {
      Int_t ret = gROOT->LoadMacro(macro.Data());
      Printf("Macro loaded from %s/%s [%d]...",gSystem->pwd(),macro.Data(),ret);
      return kTRUE;
   }

   if (!gSystem->AccessPathName(gSystem->ExpandPathName(TString::Format("%s/%s",lego_path.Data(),macro.Data()).Data()))) {
      Printf("Loading macro %s ...",gSystem->ExpandPathName(TString::Format("%s/%s",lego_path.Data(),macro.Data()).Data()));
      gROOT->LoadMacro(gSystem->ExpandPathName(TString::Format("%s/%s",lego_path.Data(),macro.Data()).Data()));
      Printf("Macro loaded from %s ...",gSystem->ExpandPathName(TString::Format("%s/%s",lego_path.Data(),macro.Data()).Data()));
      return kTRUE;
   }

   Printf("Error loading %s",macro.Data());

   return kFALSE;
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

   Bool_t valid;
   Int_t useMCMomentum = AliRsnTrainManager::GetGlobalInt("RsnUseMCMomentum",valid);
   Bool_t typeSame = (pType1 == pType2);

   Printf("------------- id1=%d id2=%d",pType1,pType2);

   TList *listLoops = new TList;

   // pair definition
   AliRsnPairDef  *pairDefPM         = new AliRsnPairDef((AliRsnDaughter::ESpecies)pType1, '+', (AliRsnDaughter::ESpecies)pType2, '-', pdgMother, massMother);
   AliRsnPairDef  *pairDefMP         = new AliRsnPairDef((AliRsnDaughter::ESpecies)pType1, '-', (AliRsnDaughter::ESpecies)pType2, '+', pdgMother, massMother);
   AliRsnPairDef  *pairDefPP         = new AliRsnPairDef((AliRsnDaughter::ESpecies)pType1, '+', (AliRsnDaughter::ESpecies)pType2, '+', pdgMother, massMother);
   AliRsnPairDef  *pairDefMM         = new AliRsnPairDef((AliRsnDaughter::ESpecies)pType1, '-', (AliRsnDaughter::ESpecies)pType2, '-', pdgMother, massMother);

   // loop object creation
   AliRsnLoopPair *lp = 0;

   // sets +-
   lp = new AliRsnLoopPair(Form("%s.RecPM", name.Data()), pairDefPM, kFALSE);
   listLoops->Add(lp);

   if (!typeSame) {
      lp = new AliRsnLoopPair(Form("%s.RecMP", name.Data()), pairDefMP, kFALSE);
      listLoops->Add(lp);
   }

   // sets +- TRUE pairs
   if (isMC) {
      lp = new AliRsnLoopPair(Form("%s.RecPM_RecMother", name.Data()), pairDefPM, kFALSE);
      lp->SetOnlyTrue(kTRUE);
      lp->SetCheckDecay(kTRUE);
      listLoops->Add(lp);
      if (!typeSame) {
         lp = new AliRsnLoopPair(Form("%s.RecMP_RecMother", name.Data()), pairDefMP, kFALSE);
         lp->SetOnlyTrue(kTRUE);
         lp->SetCheckDecay(kTRUE);
         listLoops->Add(lp);
      }
      // sets +- TRUE paris (MC is used for momentum)
      lp = new AliRsnLoopPair(Form("%s.GenPM_RecMother", name.Data()), pairDefPM, kFALSE);
      lp->SetTrueMC(kTRUE);
      listLoops->Add(lp);
      if (!typeSame) {
         // sets +- TRUE paris (MC is used for momentum)
         lp = new AliRsnLoopPair(Form("%s.GenMP_RecMother", name.Data()), pairDefMP, kFALSE);
         lp->SetTrueMC(kTRUE);
         listLoops->Add(lp);
      }
   }

   // sets ++
   lp = new AliRsnLoopPair(Form("%s.RecPP", name.Data()), pairDefPP, kFALSE);
   listLoops->Add(lp);

   // sets --
   lp = new AliRsnLoopPair(Form("%s.RecMM", name.Data()), pairDefMM, kFALSE);
   listLoops->Add(lp);

   if (isMixing) {
      // sets +- Mixing (NOT mini)
      lp = new AliRsnLoopPair(Form("%s.RecPM_mix", name.Data()), pairDefPM, kTRUE);
      listLoops->Add(lp);

      // sets -+ Mixing (NOT mini)
      lp = new AliRsnLoopPair(Form("%s.RecMP_mix", name.Data()), pairDefMP, kTRUE);
      listLoops->Add(lp);

      // sets ++ Mixing (NOT mini)
      lp = new AliRsnLoopPair(Form("%s.RecPP_mix", name.Data()), pairDefPP, kTRUE);
      listLoops->Add(lp);

      // sets -- Mixing (NOT mini)
      lp = new AliRsnLoopPair(Form("%s.RecMM_mix", name.Data()), pairDefMM, kTRUE);
      listLoops->Add(lp);
   }


   // loops over all AliRsnLoops and sets everything (don't touch it if you don't know what you are doing)
   TIter next(listLoops);
   while ((lp = (AliRsnLoopPair *)next.Next())) {
      lp->SetListID(0, listID1);
      lp->SetListID(1, listID2);
      lp->SetMCRefInfo(useMCMomentum);
      if (commonPairCuts) lp->SetPairCuts(commonPairCuts);
      if (commonEventCuts) lp->SetEventCuts(commonEventCuts);
      if (name.Contains("phi")) AddPairOutputPhi(lp);
      else if (name.Contains("kstar")) AddPairOutputKStar(lp);
      else if (name.Contains("rho")) AddPairOutputRho(lp);
      else if (name.Contains("lambda")) AddPairOutputLambda(lp);
      else if (name.Contains("sigma")) AddPairOutputSigma(lp);
      else continue;
      ((AliRsnAnalysisTask *)task)->AddLoop(lp);
   }
   return kTRUE;
}

void AddMonitorOutput(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *lm=0)
{

   Bool_t valid;
   Int_t useMCMon = AliRsnTrainManager::GetGlobalInt("RsnUseMCMonitoring",valid);
//     if (useMCMon) return;

   // dEdx tpc
   AliRsnValueDaughter *axisMomTPC = new AliRsnValueDaughter("pTPC", AliRsnValueDaughter::kPtpc);
   AliRsnValueDaughter *axisSigTPC = new AliRsnValueDaughter("sTPC", AliRsnValueDaughter::kTPCsignal);
   //axisMomTPC->SetBins(0.0,5.0,0.01);
   axisMomTPC->SetBins(0.0,12.0,0.05);
   axisSigTPC->SetBins(0.0,500.0,2.0);

   // output: 2D histogram of TPC signal vs. TPC momentum
   AliRsnListOutput *outMonitordEdxTPC = new AliRsnListOutput("dEdx", AliRsnListOutput::kHistoDefault);
   outMonitordEdxTPC->AddValue(axisMomTPC);
   outMonitordEdxTPC->AddValue(axisSigTPC);

   // add outputs to loop
   if (mon) mon->Add(outMonitordEdxTPC);
   if (lm) lm->AddOutput(outMonitordEdxTPC);

   // dEdx tpc
   AliRsnValueDaughter *axisMomTPCForTOF = new AliRsnValueDaughter("pTPC", AliRsnValueDaughter::kPtpc);
   AliRsnValueDaughter *axisSigTOF = new AliRsnValueDaughter("sTOF", AliRsnValueDaughter::kTOFsignal);
   //axisMomTPCForTOF->SetBins(0.0,5.0,0.01);
   //axisSigTOF->SetBins(0.0,500.0,2.0);
   axisMomTPCForTOF->SetBins(0.0,12.0,0.05);
   axisSigTOF->SetBins(0.0,5.e5,1.e3);

   // output: 2D histogram of TPC signal vs. TPC momentum
   AliRsnListOutput *outMonitordEdxTOF = new AliRsnListOutput("TOF", AliRsnListOutput::kHistoDefault);
   outMonitordEdxTOF->AddValue(axisMomTPCForTOF);
   outMonitordEdxTOF->AddValue(axisSigTOF);

   // add outputs to loop
   if (mon) mon->Add(outMonitordEdxTOF);
   if (lm) lm->AddOutput(outMonitordEdxTOF);


   // Momentum
   AliRsnValueDaughter *axisMomP = new AliRsnValueDaughter("p", AliRsnValueDaughter::kP);
   //axisMomP->SetBins(0.0,5.0,0.01);
   axisMomP->SetBins(0.0,12.0,0.05);

   // output: 2D histogram of TPC signal vs. TPC momentum
   AliRsnListOutput *outMonitorP = new AliRsnListOutput("P", AliRsnListOutput::kHistoDefault);
   outMonitorP->AddValue(axisMomP);

   // add outputs to loop
   if (mon) mon->Add(outMonitorP);
   if (lm) lm->AddOutput(outMonitorP);


   if (useMCMon) {
      AliRsnValueDaughter *axisMomPMC = new AliRsnValueDaughter("pMC", AliRsnValueDaughter::kP);
      axisMomPMC->SetUseMCInfo(kTRUE);
      //axisMomPMC->SetBins(0.0,5.0,0.01);
      axisMomPMC->SetBins(0.0,12.0,0.05);

      // output: 2D histogram of TPC signal vs. TPC momentum
      AliRsnListOutput *outMonitorPMC = new AliRsnListOutput("PMC", AliRsnListOutput::kHistoDefault);
      outMonitorPMC->AddValue(axisMomPMC);

      // add outputs to loop
      if (mon) mon->Add(outMonitorPMC);
      if (lm) lm->AddOutput(outMonitorPMC);
   }


   // Momentum Pt
   AliRsnValueDaughter *axisMomPt = new AliRsnValueDaughter("pt", AliRsnValueDaughter::kPt);
   //axisMomPt->SetBins(0.0,5.0,0.01);
   axisMomPt->SetBins(0.0,12.0,0.05);

   // output: 2D histogram of TPC signal vs. TPC momentum
   AliRsnListOutput *outMonitorPt = new AliRsnListOutput("Pt", AliRsnListOutput::kHistoDefault);
   outMonitorPt->AddValue(axisMomPt);

   // add outputs to loop
   if (mon) mon->Add(outMonitorPt);
   if (lm) lm->AddOutput(outMonitorPt);
   if (useMCMon) {
      // Momentum Pt
      AliRsnValueDaughter *axisMomPtMC = new AliRsnValueDaughter("ptMC", AliRsnValueDaughter::kPt);
      axisMomPtMC->SetUseMCInfo(kTRUE);
      //axisMomPtMC->SetBins(0.0,5.0,0.01);
      axisMomPtMC->SetBins(0.0,12.0,0.05);

      // output: 2D histogram of TPC signal vs. TPC momentum
      AliRsnListOutput *outMonitorPtMC = new AliRsnListOutput("PtMC", AliRsnListOutput::kHistoDefault);
      outMonitorPtMC->AddValue(axisMomPtMC);

      // add outputs to loop
      if (mon) mon->Add(outMonitorPtMC);
      if (lm) lm->AddOutput(outMonitorPtMC);
   }

   // Eta
   AliRsnValueDaughter *axisMomEta = new AliRsnValueDaughter("eta", AliRsnValueDaughter::kEta);
   axisMomEta->SetBins(-1.0,1.0,0.01);

   // output: 2D histogram of TPC signal vs. TPC momentum
   AliRsnListOutput *outMonitorEta = new AliRsnListOutput("Eta", AliRsnListOutput::kHistoDefault);
   outMonitorEta->AddValue(axisMomEta);

   // add outputs to loop
   if (mon) mon->Add(outMonitorEta);
   if (lm) lm->AddOutput(outMonitorEta);

   if (useMCMon) {
      // Eta
      AliRsnValueDaughter *axisMomEtaMC = new AliRsnValueDaughter("etaMC", AliRsnValueDaughter::kEta);
      axisMomEtaMC->SetUseMCInfo(kTRUE);
      axisMomEtaMC->SetBins(-1.0,1.0,0.01);

      // output: 2D histogram of TPC signal vs. TPC momentum
      AliRsnListOutput *outMonitorEtaMC = new AliRsnListOutput("EtaMC", AliRsnListOutput::kHistoDefault);
      outMonitorEtaMC->AddValue(axisMomEtaMC);

      // add outputs to loop
      if (mon) mon->Add(outMonitorEtaMC);
      if (lm) lm->AddOutput(outMonitorEtaMC);
   }

//   AliRsnValueDaughter *axisPtBig = new AliRsnValueDaughter("pt", AliRsnValueDaughter::kPt);
   AliRsnValueDaughter *axisPtBig = new AliRsnValueDaughter("p", AliRsnValueDaughter::kP);
   axisPtBig->SetBins(0.0,12.0,0.1);

   // kTOFnsigmaK
   AliRsnValueDaughter *axisTPCnsigmaK = new AliRsnValueDaughter("K", AliRsnValueDaughter::kTPCnsigmaK);
   axisTPCnsigmaK->SetBins(1000,-100,100);

   // output: 2D histogram of TPC signal vs. TPC momentum
   AliRsnListOutput *outMonitorTPCnsigmaK = new AliRsnListOutput("TPC_nsigma", AliRsnListOutput::kHistoDefault);
   outMonitorTPCnsigmaK->AddValue(axisPtBig);
   outMonitorTPCnsigmaK->AddValue(axisTPCnsigmaK);

   // add outputs to loop
   if (mon) mon->Add(outMonitorTPCnsigmaK);
   if (lm) lm->AddOutput(outMonitorTPCnsigmaK);

   // kTPCnsigmaPi
   AliRsnValueDaughter *axisTPCnsigmaPi = new AliRsnValueDaughter("pi", AliRsnValueDaughter::kTPCnsigmaPi);
   axisTPCnsigmaPi->SetBins(1000,-100,100);

   // output: 2D histogram of TPC signal vs. TPC momentum
   AliRsnListOutput *outMonitorTPCnsigmaPi = new AliRsnListOutput("TPC_nsigma", AliRsnListOutput::kHistoDefault);
   outMonitorTPCnsigmaPi->AddValue(axisPtBig);
   outMonitorTPCnsigmaPi->AddValue(axisTPCnsigmaPi);

   // add outputs to loop
   if (mon) mon->Add(outMonitorTPCnsigmaPi);
   if (lm) lm->AddOutput(outMonitorTPCnsigmaPi);

   // kTPCnsigmaP
   AliRsnValueDaughter *axisTPCnsigmaP = new AliRsnValueDaughter("p", AliRsnValueDaughter::kTPCnsigmaP);
   axisTPCnsigmaP->SetBins(1000,-100,100);

   // output: 2D histogram of TPC signal vs. TPC momentum
   AliRsnListOutput *outMonitorTPCnsigmaP = new AliRsnListOutput("TPC_nsigma", AliRsnListOutput::kHistoDefault);
   outMonitorTPCnsigmaP->AddValue(axisPtBig);
   outMonitorTPCnsigmaP->AddValue(axisTPCnsigmaP);

   // add outputs to loop
   if (mon) mon->Add(outMonitorTPCnsigmaP);
   if (lm) lm->AddOutput(outMonitorTPCnsigmaP);


   if (!opt.Contains("NoTOFSIGMA")) {

      // kTOFnsigmaK
      AliRsnValueDaughter *axisTOFnsigmaK = new AliRsnValueDaughter("K", AliRsnValueDaughter::kTOFnsigmaK);
      axisTOFnsigmaK->SetBins(1000,-100,100);

      // output: 2D histogram of TPC signal vs. TPC momentum
      AliRsnListOutput *outMonitorTOFnsigmaK = new AliRsnListOutput("TOF_nsigma", AliRsnListOutput::kHistoDefault);
      outMonitorTOFnsigmaK->AddValue(axisPtBig);
      outMonitorTOFnsigmaK->AddValue(axisTOFnsigmaK);

      // add outputs to loop
      if (mon) mon->Add(outMonitorTOFnsigmaK);
      if (lm) lm->AddOutput(outMonitorTOFnsigmaK);

      // kTOFnsigmaPi
      AliRsnValueDaughter *axisTOFnsigmaPi = new AliRsnValueDaughter("pi", AliRsnValueDaughter::kTOFnsigmaPi);
      axisTOFnsigmaPi->SetBins(1000,-100,100);

      // output: 2D histogram of TPC signal vs. TPC momentum
      AliRsnListOutput *outMonitorTOFnsigmaPi = new AliRsnListOutput("TOF_nsigma", AliRsnListOutput::kHistoDefault);
      outMonitorTOFnsigmaPi->AddValue(axisPtBig);
      outMonitorTOFnsigmaPi->AddValue(axisTOFnsigmaPi);

      // add outputs to loop
      if (mon) mon->Add(outMonitorTOFnsigmaPi);
      if (lm) lm->AddOutput(outMonitorTOFnsigmaPi);

      // kTOFnsigmaP
      AliRsnValueDaughter *axisTOFnsigmaP = new AliRsnValueDaughter("p", AliRsnValueDaughter::kTOFnsigmaP);
      axisTOFnsigmaP->SetBins(1000,-100,100);

      // output: 2D histogram of TPC signal vs. TPC momentum
      AliRsnListOutput *outMonitorTOFnsigmaP = new AliRsnListOutput("TOF_nsigma", AliRsnListOutput::kHistoDefault);
      outMonitorTOFnsigmaP->AddValue(axisPtBig);
      outMonitorTOFnsigmaP->AddValue(axisTOFnsigmaP);

      // add outputs to loop
      if (mon) mon->Add(outMonitorTOFnsigmaP);
      if (lm) lm->AddOutput(outMonitorTOFnsigmaP);
   }


   // nITCcluster
   AliRsnValueDaughter *axisITSnClusters = new AliRsnValueDaughter("ITS", AliRsnValueDaughter::kNITSclusters);
   axisITSnClusters->SetBins(8,-0.5,7.5);

   AliRsnListOutput *outMonitorITSnClusters = new AliRsnListOutput("nClusters", AliRsnListOutput::kHistoDefault);
   outMonitorITSnClusters->AddValue(axisPtBig);
   outMonitorITSnClusters->AddValue(axisITSnClusters);
   // add outputs to loop
   if (mon) mon->Add(outMonitorITSnClusters);
   if (lm) lm->AddOutput(outMonitorITSnClusters);

   // nTPCcluster
   AliRsnValueDaughter *axisTPCnClusters = new AliRsnValueDaughter("TPC", AliRsnValueDaughter::kNTPCclusters);
   axisTPCnClusters->SetBins(166,-0.5,165.5);

   AliRsnListOutput *outMonitorTPCnClusters = new AliRsnListOutput("nClusters", AliRsnListOutput::kHistoDefault);
   outMonitorTPCnClusters->AddValue(axisPtBig);
   outMonitorTPCnClusters->AddValue(axisTPCnClusters);
   // add outputs to loop
   if (mon) mon->Add(outMonitorTPCnClusters);
   if (lm) lm->AddOutput(outMonitorTPCnClusters);

   // ITSchi2
   AliRsnValueDaughter *axisITSchi2 = new AliRsnValueDaughter("ITS", AliRsnValueDaughter::kITSchi2);
   axisITSchi2->SetBins(100,0,10);

   AliRsnListOutput *outMonitorITSchi2 = new AliRsnListOutput("chi2", AliRsnListOutput::kHistoDefault);
   outMonitorITSchi2->AddValue(axisPtBig);
   outMonitorITSchi2->AddValue(axisITSchi2);
   // add outputs to loop
   if (mon) mon->Add(outMonitorITSchi2);
   if (lm) lm->AddOutput(outMonitorITSchi2);

   // TPCchi2
   AliRsnValueDaughter *axisTPCchi2 = new AliRsnValueDaughter("TPC", AliRsnValueDaughter::kTPCchi2);
   axisTPCchi2->SetBins(100,0,10);

   AliRsnListOutput *outMonitorTPCchi2 = new AliRsnListOutput("chi2", AliRsnListOutput::kHistoDefault);
   outMonitorTPCchi2->AddValue(axisPtBig);
   outMonitorTPCchi2->AddValue(axisTPCchi2);
   // add outputs to loop
   if (mon) mon->Add(outMonitorTPCchi2);
   if (lm) lm->AddOutput(outMonitorTPCchi2);

   // DCAXY
   AliRsnValueDaughter *axisDCAXY = new AliRsnValueDaughter("XY", AliRsnValueDaughter::kDCAXY);
   axisDCAXY->SetBins(200,-1,1);

   AliRsnListOutput *outMonitorDCAXY = new AliRsnListOutput("DCA", AliRsnListOutput::kHistoDefault);
   outMonitorDCAXY->AddValue(axisPtBig);
   outMonitorDCAXY->AddValue(axisDCAXY);
   // add outputs to loop
   if (mon) mon->Add(outMonitorDCAXY);
   if (lm) lm->AddOutput(outMonitorDCAXY);

   // DCAZ
   AliRsnValueDaughter *axisDCAZ = new AliRsnValueDaughter("Z", AliRsnValueDaughter::kDCAZ);
   axisDCAZ->SetBins(200,-1,1);

   AliRsnListOutput *outMonitorDCAZ = new AliRsnListOutput("DCA", AliRsnListOutput::kHistoDefault);
   outMonitorDCAZ->AddValue(axisPtBig);
   outMonitorDCAZ->AddValue(axisDCAZ);
   // add outputs to loop
   if (mon) mon->Add(outMonitorDCAZ);
   if (lm) lm->AddOutput(outMonitorDCAZ);

   AliRsnListOutput *outMonitorPTvsMult = new AliRsnListOutput("PTvsMult",AliRsnListOutput::kHistoDefault);
   AliRsnValueDaughter *vd1 = new AliRsnValueDaughter("pt",AliRsnValueDaughter::kPt);
   vd1->SetBins(0.0,5.0,0.01);
   outMonitorPTvsMult->AddValue(vd1);

   AliRsnValueEvent *ve1 = new AliRsnValueEvent("mult",AliRsnValueEvent::kMult);
   ve1->SetBins(0.0,100.0,1);
   outMonitorPTvsMult->AddValue(ve1);
   if (mon) mon->Add(outMonitorPTvsMult);
   if (lm) lm->AddOutput(outMonitorPTvsMult);

//    AliRsnListOutput *outMonitorMult = new AliRsnListOutput("EventMult",AliRsnListOutput::kHistoDefault);
//
//    AliRsnValueEvent *ve1Multi = new AliRsnValueEvent("centrality",AliRsnValueEvent::kCentralityV0);
//    ve1Multi->SetBins(0.0,100,10.0);
//    outMonitorMult->AddValue(ve1Multi);
//    if (mon) mon->Add(outMonitorMult);
//    if (lm) lm->AddOutput(outMonitorMult);
}

void AddMonitorOutputMini(AliRsnMiniMonitorTask *task,Int_t listID1,TString name = "",Char_t charge='0')
{
   TString chargeName="all";
   if ( charge == '+' ) chargeName = "pos";
   if ( charge == '-' ) chargeName = "neg";

   AliRsnMiniMonitor *mondEdx = task->CreateMonitor(Form("%s_dEdx_pTPC_%s", name.Data(),chargeName.Data()),AliRsnMiniMonitor::kdEdxTPCvsP, listID1);
   mondEdx->SetCharge(charge);
   AliRsnMiniMonitor *monPt = task->CreateMonitor(Form("%s_Pt_%s", name.Data(),chargeName.Data()),AliRsnMiniMonitor::kTrackPt, listID1);
   monPt->SetCharge(charge);
}

void AddParticleMonitor(AliAnalysisTaskSE *task, Bool_t isMC, Int_t listID1,AliRsnCutSet *commonEventCuts=0,AliRsnCutSet *cutPair=0,TString name = "")
{
   Bool_t valid;
   Int_t isRsnMini = AliRsnTrainManager::GetGlobalInt("IsRsnMini",valid);
   Int_t useMCMon = AliRsnTrainManager::GetGlobalInt("RsnUseMCMonitoring",valid);

   if (isRsnMini) {
//       Printf("Monitoring by mini is not supported now. It will be soon !!!");
//       return ;
      AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
      AliRsnMiniMonitorTask *monTask = new AliRsnMiniMonitorTask(name.Data(),useMCMon);
      AddMonitorOutputMini(monTask,listID1,name);
//       AddMonitorOutputMini(monTask,listID1,name,'+');
//       AddMonitorOutputMini(monTask,listID1,name,'-');
      mgr->AddTask(monTask);
      // connect input container according to source choice
      mgr->ConnectInput(monTask, 0, mgr->GetCommonInputContainer());

      // create paths for the output in the common file
      TString commonPath = AliAnalysisManager::GetCommonFileName();

      // create containers for output
      AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnMonMini%s", name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, commonPath.Data());

      mgr->ConnectOutput(monTask, 1, output);



   } else {

      TList *listLoops = new TList;
      // monitor definition
      AliRsnDaughterDef *tracksAll = new AliRsnDaughterDef(AliRsnDaughter::kTrack /*'+' or '-'*/);
// //       AliRsnDaughterDef *tracksPos = new AliRsnDaughterDef(AliRsnDaughter::kTrack,'+');
// //       AliRsnDaughterDef *tracksNeg = new AliRsnDaughterDef(AliRsnDaughter::kTrack,'-');
//
      AliRsnLoopDaughter *lm =0;
//       // loop object
      listLoops->Add(new AliRsnLoopDaughter(Form("ALL_%s", name.Data()), listID1, tracksAll));
//
// //       listLoops->Add(new AliRsnLoopDaughter(Form("%s_pos", name.Data()), listID1, tracksPos));
// //       listLoops->Add(new AliRsnLoopDaughter(Form("%s_neg", name.Data()), listID1, tracksNeg));
//
      TIter next(listLoops);
      while ((lm = (AliRsnLoopDaughter *)next.Next())) {
//          if (commonEventCuts) lm->SetEventCuts(commonEventCuts);
         AddMonitorOutput(0,"mc_loop",lm);
         ((AliRsnAnalysisTask *)task)->AddLoop(lm);
      }
   }
}

