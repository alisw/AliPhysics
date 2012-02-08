#ifndef __CINT__
#include <TSystem.h>
#include <TROOT.h>
#include <Rtypes.h>
#include <TString.h>
#include <TNamed.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TList.h>
#include <TStopwatch.h>
#endif

Bool_t RunALICE(TString anSrc = "grid",
                TString anMode = "terminate",
                TString input="aod" /*or "esd"*/,
                TString inputMC="" /*or "mc"*/,
                Long64_t nEvents = 1e10,
                Long64_t nSkip = 0,
                TString dsName="",
                TString alirsnliteManagers ="AddAMRsn",
                Bool_t useMultiHandler=kTRUE,
                TString alirsnlitesrc ="$ALICE_ROOT",
                TString alirsnlitetasks =""
               ) {

   // some init work
   anSrc.ToLower(); anMode.ToLower(); input.ToLower(); inputMC.ToLower();

   // loads libs and setup include paths
   if (LoadLibsBase(alirsnlitesrc)) return kFALSE;

   // reset manager if already exists
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (mgr) delete mgr;
   mgr = new AliAnalysisManager("AliRsnLiteAM","AliRsnLite Analysis Manager");

   Bool_t useAODOut = kFALSE;
   CreateInputHandlers(input,inputMC,useAODOut,useMultiHandler);

   // add default grid handler
   gROOT->LoadMacro("SetupAnalysisPlugin.C");
   AliAnalysisGrid *analysisPlugin = SetupAnalysisPlugin(anMode.Data());
   if (!analysisPlugin) { Printf("Error : analysisPlugin is null !!!"); return kFALSE; }
   mgr->SetGridHandler(analysisPlugin);
   if (!dsName.IsNull()) {
      if (!anSrc.CompareTo("proof") && !anMode.CompareTo("full")) {
         analysisPlugin->SetProofDataSet(dsName.Data());
         Printf(Form("Using DataSet %s ...",dsName.Data()));
      } else {
         analysisPlugin->SetFileForTestMode(dsName.Data());
         Printf(Form("Using Test file %s ...",dsName.Data()));
      }
   }

   TList *listManagers = CreateListOfManagersFromDir(alirsnliteManagers,alirsnlitetasks);
   if (!listManagers) { Printf("Error : CreateListOfManagersFromDir failed !!!"); return kFALSE;}

   // adds all tasks
   if (!AddAllManagers(listManagers, anSrc, anMode,input,inputMC)) { Printf("Error : AddAllManagers failed !!!"); return kFALSE;}

   TStopwatch timer;
   timer.Start();
   // runs analysis
   if (!RunAnalysisManager(anSrc, anMode.Data(), nEvents, nSkip)) { Printf("Error : RunAnalysisManager failed !!!"); return kFALSE;}

   timer.Stop();
   timer.Print();
   Printf("Working directory is %s ...", gSystem->WorkingDirectory());
   Printf("Done OK");
   return kTRUE;

}

Int_t LoadLibsBase(TString alirsnlitesrc) {
   Int_t num = 0;
   if (gSystem->Load("libTree.so") < 0) {num++; return num;}
   if (gSystem->Load("libGeom.so") < 0) {num++; return num;}
   if (gSystem->Load("libVMC.so") < 0) {num++; return num;}
   if (gSystem->Load("libMinuit.so") < 0) {num++; return num;}
   if (gSystem->Load("libPhysics.so") < 0) {num++; return num;}
   if (gSystem->Load("libSTEERBase.so") < 0) {num++; return num;}
   if (gSystem->Load("libESD.so") < 0) {num++; return num;}
   if (gSystem->Load("libAOD.so") < 0) {num++; return num;}
   if (gSystem->Load("libANALYSIS.so") < 0) {num++; return num;}
   if (gSystem->Load("libOADB.so") < 0) {num++; return num;}
   if (gSystem->Load("libANALYSISalice.so") < 0) {num++; return num;}

   gSystem->AddIncludePath(Form("-I\"%s/include\"", gSystem->ExpandPathName(alirsnlitesrc.Data())));
   gROOT->ProcessLine(Form(".include %s/include", gSystem->ExpandPathName(alirsnlitesrc.Data())));

   return 0;
}

Bool_t CreateInputHandlers(TString input,TString inputMC,Bool_t useAODOut=kFALSE,Bool_t useMultiHandler=kTRUE) {

   Bool_t useMC = !inputMC.CompareTo("mc");

   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) { Printf("Error [CreateInputHandlers] : mgr is null !!!"); return kFALSE; }

   if (useMultiHandler) {
      AliMultiInputEventHandler *inputHandler = new AliMultiInputEventHandler();
      if (!input.CompareTo("esd")) {
         inputHandler->AddInputEventHandler(new AliESDInputHandler());
         if (useMC) inputHandler->AddInputEventHandler(new AliMCEventHandler());
      } else if (!input.CompareTo("aod")) {
         inputHandler->AddInputEventHandler(new AliAODInputHandler());
      }

      mgr->SetInputEventHandler(inputHandler);
   } else {
      if (!input.CompareTo("esd")) {
         mgr->SetInputEventHandler(new AliESDInputHandler());
         if (useMC) mgr->SetMCtruthEventHandler(new AliMCEventHandler());
      } else if (!input.CompareTo("aod")) {
         mgr->SetInputEventHandler(new AliAODInputHandler());
      }
      mgr->SetInputEventHandler(inputHandler);
   }

   if (useAODOut) {
      AliAODHandler *aodHandler   = new AliAODHandler();
      aodHandler->SetOutputFileName("AliAOD.root");
      mgr->SetOutputEventHandler(aodHandler);
   }

   return kTRUE;

}

TList *CreateListOfManagersFromDir(TString listManagersNames="",TString dir="") {

   TList *listManagers = new TList;
   TString dirsStr;
   TObjArray *dirs=0;

   if (listManagersNames.IsNull()) {
      if (dir.IsNull() || gSystem->AccessPathName(gSystem->ExpandPathName(dir.Data()))) {
         Printf("Error [CreateListOfManagersFromDir] : Dir '%s' doesn't exists !!!",dir.Data());
         return 0;
      }
      dirsStr = gSystem->GetFromPipe(Form("ls %s",dir.Data()));
      dirs = dirsStr.Tokenize("\n");
   } else {
      dirsStr = listManagersNames;
      dirs = dirsStr.Tokenize(" ");
   }

   TIter next(dirs);
   Int_t counter=0;
   TObjString *str,*strtmp;
   TObjArray *mydirstrTok;
   TString mydirstr,main,prefix;
   while ((str = (TObjString *)next.Next())) {
      // TODO add direcotry
      mydirstr = str->GetString();
      if (mydirstr.IsNull()) continue;

      Printf("Adding %s .,,",mydirstr.Data());
      mydirstrTok = mydirstr.Tokenize("_");

      main = ((TObjString *)mydirstrTok->At(0))->GetString();

      strtmp = (TObjString *)mydirstrTok->At(1);
      if (strtmp) prefix = strtmp->GetString(); else prefix="";

      listManagers->Add(new TNamed(main,prefix));
   }

   return listManagers;
}

Bool_t AddAllManagers(TList *listManagers,TString anSrc, TString anMode,TString input,TString inputMC) {
   TIter next(listManagers);
   Int_t counter=0;
   TNamed *name;
   while ((name = (TNamed *)next.Next())) {
      if (!AddAnalysisManager(name->GetName(), anSrc, anMode,input,inputMC,name->GetTitle(),Form("%d",counter++))) {
         Printf("Error: Problem adding %s",name->GetName());
         return kFALSE;
      }
   }

   return kTRUE;
}

Bool_t AddAnalysisManager(TString managerMacro, TString anSrc, TString anMode,TString input,TString inputMC, TString postfix,TString idStr) {
   gROOT->LoadMacro(Form("%s.C", managerMacro.Data()));
   return gROOT->ProcessLine(Form("%s\(\"%s\",\"%s\",\"%s\"\,\"%s\",\"%s\",\"%s\"\);", managerMacro.Data(), anSrc.Data(), anMode.Data(),input.Data(),inputMC.Data(), postfix.Data(),idStr.Data()));
}

Bool_t RunAnalysisManager(TString anSrc = "proof", TString anMode = "test", Long64_t nEvents = 1e10, Long64_t nSkip = 0) {

   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

   if (!mgr) { Printf("Error [RunAnalysisManager] : mgr is null !!!"); return kFALSE; }

   // Run analysis
   mgr->InitAnalysis();
   mgr->PrintStatus();

   if ((!anSrc.CompareTo("proof")) || (!anSrc.CompareTo("local"))) {
      mgr->StartAnalysis(anSrc.Data(), nEvents, nSkip);
   } else {
      mgr->StartAnalysis(anSrc.Data());
   }

   return kTRUE;
}
