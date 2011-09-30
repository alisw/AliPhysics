// TODO:
// 1. run with Many centrality bins at once
#include <string.h>

enum { kMyRunModeLocal = 0, kMyRunModeCAF, kMyRunModeGRID};

TList * listToLoad = new TList();

TChain * GetAnalysisChain(const char * incollection);

void runProof(Char_t* data = "/alice/data/LHC10h_000138653_p2_AOD049", Long64_t nev = -1, Long64_t offset = 0, Bool_t debug = kFALSE, Int_t runMode = 1, Bool_t isMC = 0,
              const char* option = "SAVE", Int_t workers = -1)
{
   // runMode:
   //
   // 0 local
   // 1 proof
   // 2 grid

   if (nev < 0)
      nev = 1234567890;

   InitAndLoadLibs(runMode, workers, debug);

   // Create the analysis manager
   mgr = new AliAnalysisManager;

   // Add AOD handler
   AliAODInputHandler* AODH = new AliAODInputHandler;
   mgr->SetInputEventHandler(AODH);
   mgr->SetDebugLevel(10);

   // Is this needed for AOD?
   // if(isMC) {
   //   AliMCEventHandler* handler = new AliMCEventHandler;
   //   handler->SetPreReadMode(AliMCEventHandler::kLmPreRead);
   //   mgr->SetMCtruthEventHandler(handler);
   // }


   // If we are running on grid, we need the alien handler
   if (runMode == kMyRunModeGRID)
   {
      // Create and configure the alien handler plugin
      gROOT->LoadMacro("CreateAlienHandler.C");
      AliAnalysisGrid *alienHandler = CreateAlienHandler(data, listToLoad, "full", isMC);
      if (!alienHandler)
      {
         cout << "Cannot create alien handler" << endl;
         exit(1);
      }
      mgr->SetGridHandler(alienHandler);
   }

   // Parse option strings
   TString optionStr(option);

   // remove SAVE option if set
   // This  is copied from a macro by Jan. The reason I kept it is that I may want to pass textual options to the new task at some point
   Bool_t doSave = kFALSE;
   TString optionStr(option);
   if (optionStr.Contains("SAVE"))
   {
      optionStr = optionStr(0, optionStr.Index("SAVE")) + optionStr(optionStr.Index("SAVE") + 4, optionStr.Length());
      doSave = kTRUE;
   }

   //  AliLog::SetClassDebugLevel("AliESDtrackCuts", AliLog::kDebug);// FIXME

   // load my task
   gROOT->ProcessLine(".L AddTaskSpectraAOD.C");
   AliAnalysisTaskSpectraAOD * task = AddTaskSpectraAOD("SpectraAOD.root");
   task->SetIsMC(isMC);
   // Init and run the analy
   if (!mgr->InitAnalysis()) return;

   mgr->PrintStatus();

   if (runMode == kMyRunModeLocal)
   {
      // If running in local mode, create chain of ESD files
      cout << "RUNNING LOCAL, CHAIN" << endl;
      TChain * chain = GetAnalysisChain(data);
      //    chain->Print();
      mgr->StartAnalysis("local", chain, nev);
   }
   else if (runMode == kMyRunModeCAF)
   {
      mgr->StartAnalysis("proof", TString(data) + "#aodTree", nev);
   }
   else if (runMode == kMyRunModeGRID)
   {
      mgr->StartAnalysis("grid");
   }
   else
   {
      cout << "ERROR: unknown run mode" << endl;
   }

   TString pathsuffix = "";
   if (doSave) MoveOutput(data, pathsuffix.Data());

}


void MoveOutput(const char * data, const char * suffix = "")
{

   TString path("output/");
   path = path + TString(data).Tokenize("/")->Last()->GetName() + suffix;

   TString fileName = "SpectraAOD.root";
   gSystem->mkdir(path, kTRUE);
   gSystem->Rename(fileName, path + "/" + fileName);
   // for(Int_t ibin = 0; ibin < 20; ibin++){
   //   TString fileBin = fileName;
   //   fileBin.ReplaceAll(".root",Form("_%2.2d.root",ibin));
   //   gSystem->Rename(fileBin, path + "/" + fileBin);
   // }
   Printf(">>>>> Moved files to %s", path.Data());
}



TChain * GetAnalysisChain(const char * incollection)
{
   // Builds a chain of esd files
   // incollection can be
   // - a single root file
   // - an xml collection of files on alien
   // - a ASCII containing a list of local root files
   TChain* analysisChain = 0;
   // chain
   analysisChain = new TChain("aodTree");
   if (TString(incollection).Contains(".root"))
   {
      analysisChain->Add(incollection);
   }
   else if (TString(incollection).Contains("xml"))
   {
      TGrid::Connect("alien://");
      TAlienCollection * coll = TAlienCollection::Open(incollection);
      while (coll->Next())
      {
         analysisChain->Add(TString("alien://") + coll->GetLFN());
      }
   }
   else
   {
      ifstream file_collect(incollection);
      TString line;
      while (line.ReadLine(file_collect))
      {
         analysisChain->Add(line.Data());
      }
   }
   analysisChain->GetListOfFiles()->Print();

   return analysisChain;
}


void InitAndLoadLibs(Int_t runMode = kMyRunModeLocal, Int_t workers = 0, Bool_t debug = 0)
{
   // Loads libs and par files + custom task and classes

   // Custom stuff to be loaded

   listToLoad->Add(new TObjString("AliSpectraAODHistoManager.cxx+"));
   listToLoad->Add(new TObjString("AliSpectraAODEventCuts.cxx+"));
   listToLoad->Add(new TObjString("AliSpectraAODTrackCuts.cxx+"));
   listToLoad->Add(new TObjString("AliAnalysisTaskSpectraAOD.cxx+"));

   if (runMode == kMyRunModeCAF)
   {
      cout << "Init in CAF mode" << endl;

      gEnv->SetValue("XSec.GSI.DelegProxy", "2");
      TProof * p = TProof::Open("alice-caf.cern.ch", workers > 0 ? Form("workers=%d", workers) : "1x");
      //TProof * p = TProof::Open("skaf.saske.sk", workers>0 ? Form("workers=%d",workers) : "");
      p->Exec("TObject *o = gEnv->GetTable()->FindObject(\"Proof.UseMergers\"); gEnv->GetTable()->Remove(o);", kTRUE);

      gProof->EnablePackage("VO_ALICE@AliRoot::v4-21-29-AN");
      gSystem->Load("libCore.so");
      gSystem->Load("libTree.so");
      gSystem->Load("libGeom.so");
      gSystem->Load("libVMC.so");
      gSystem->Load("libPhysics.so");
      gSystem->Load("libSTEERBase");
      gSystem->Load("libESD");
      gSystem->Load("libAOD");
      gSystem->Load("libANALYSIS");
      gSystem->Load("libOADB");
      gSystem->Load("libANALYSISalice");

      // Enable the needed package
      // gProof->UploadPackage("$ALICE_ROOT/obj/STEERBase");
      // gProof->EnablePackage("$ALICE_ROOT/obj/STEERBase");
      // gProof->UploadPackage("$ALICE_ROOT/obj/ESD");
      // gProof->EnablePackage("$ALICE_ROOT/obj/ESD");
      // gProof->UploadPackage("$ALICE_ROOT/obj/AOD");
      // gProof->EnablePackage("$ALICE_ROOT/obj/AOD");
      // gProof->UploadPackage("$ALICE_ROOT/obj/ANALYSIS");
      // gProof->EnablePackage("$ALICE_ROOT/obj/ANALYSIS");
      // gProof->UploadPackage("$ALICE_ROOT/obj/OADB");
      // gProof->EnablePackage("$ALICE_ROOT/obj/OADB");
      // gProof->UploadPackage("$ALICE_ROOT/obj/ANALYSISalice");
      // gProof->EnablePackage("$ALICE_ROOT/obj/ANALYSISalice");
      // gProof->UploadPackage("$ALICE_ROOT/obj/PWG0base");
      // gProof->EnablePackage("$ALICE_ROOT/obj/PWG0base");
      // gROOT->ProcessLine(gSystem->ExpandPathName(".include $ALICE_ROOT/PWG0/multPb"));
      // gROOT->ProcessLine(gSystem->ExpandPathName(".include $ALICE_ROOT/PWG1/background"));
   }
   else
   {
      cout << "Init in Local or Grid mode" << endl;
      gSystem->Load("libCore.so");
      gSystem->Load("libTree.so");
      gSystem->Load("libGeom.so");
      gSystem->Load("libVMC.so");
      gSystem->Load("libPhysics.so");
      gSystem->Load("libSTEERBase");
      gSystem->Load("libESD");
      gSystem->Load("libAOD");
      gSystem->Load("libANALYSIS");
      gSystem->Load("libOADB");
      gSystem->Load("libANALYSISalice");
      // Use AliRoot includes to compile our task
      gROOT->ProcessLine(".include $ALICE_ROOT/include");

      // gSystem->Load("libVMC");
      // gSystem->Load("libTree");
      // gSystem->Load("libSTEERBase");
      // gSystem->Load("libESD");
      // gSystem->Load("libAOD");
      // gSystem->Load("libANALYSIS");
      // gSystem->Load("libANALYSISalice");
      // gSystem->Load("libPWG0base");

      // gROOT->ProcessLine(gSystem->ExpandPathName(".include $ALICE_ROOT/PWG0/multPb"));
      // gROOT->ProcessLine(gSystem->ExpandPathName(".include $ALICE_ROOT/PWG1/background"));
      //    gROOT->ProcessLine(gSystem->ExpandPathName(".include $ALICE_ROOT/PWG1/background/"));
   }
   // Load helper classes
   TIterator * iter = listToLoad->MakeIterator();
   TObjString * name = 0;
   while (name = (TObjString *)iter->Next())
   {
      gSystem->ExpandPathName(name->String());
      cout << name->String().Data();
      if (runMode == kMyRunModeCAF)
      {
         gProof->Load(name->String() + (debug ? "+g" : ""));
      }
      else
      {
         gROOT->LoadMacro(name->String() + (debug ? "+g" : ""));
      }
   }

}
