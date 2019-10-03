#include "FILTER_ESD2AODMUON.h"

#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "TChain.h"
#include "AliAODHandler.h"
#include "Riostream.h"
#include "AliLog.h"
#include "TROOT.h"
#include "TSystem.h"

namespace AAF {

int FILTER_ESD2AODMUON(const char* from, const char* to)
{
  AliLog::SetGlobalLogLevel(AliLog::kWarning);
  
  TString cdir = gSystem->WorkingDirectory();
  
  TString odir;
  
  odir.Form("/tmp/ALIPHYSICS_AAF_FILTER_ESD2AODMUON_%d",gSystem->GetPid());
  
  gSystem->MakeDirectory(odir.Data());

  gSystem->ChangeDirectory(odir.Data());

  AliAnalysisManager *mgr = new AliAnalysisManager("AOD2MUONAOD");
  
  AliInputEventHandler* input = new AliESDInputHandler;
  
  mgr->SetInputEventHandler(input);
  
  AliAODHandler* aodHandler = new AliAODHandler;
  
  aodHandler->SetOutputFileName("AliAOD.root");
  mgr->SetOutputEventHandler(aodHandler);

  gROOT->ProcessLine(".L $ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
  
  gROOT->ProcessLine("AddTaskMultSelection();");
  
  gROOT->ProcessLine(".L $ALICE_ROOT/ANALYSIS/ESDfilter/macros/AddTaskESDFilter.C");

  TString cmd;
  
  Bool_t useKineFilter=kFALSE;
  Bool_t writeMuonAOD=kTRUE;
  Bool_t writeDimuonAOD=kFALSE; /*obsolete*/
  Bool_t usePhysicsSelection=kTRUE;
  Bool_t useCentralityTask=kFALSE; /*obsolete*/
  Bool_t enableTPCOnlyAODTracks=kFALSE;
  Bool_t disableCascades=kTRUE;
  Bool_t disableKinks=kTRUE;
  Int_t runFlag = 1500; // The first 2 digits are the year, the second
  Int_t  muonMCMode = 3;
  Bool_t useV0Filter=kTRUE;
  Bool_t muonWithSPDTracklets=kFALSE;

  cmd.Form("AliAnalysisTaskESDfilter* task = AddTaskESDFilter(%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d);",
           useKineFilter,
           writeMuonAOD,
           writeDimuonAOD,
           usePhysicsSelection,
           useCentralityTask,
           enableTPCOnlyAODTracks,
           disableCascades,
           disableKinks,
           runFlag,
           muonMCMode,
           useV0Filter,
           muonWithSPDTracklets);
  
  gROOT->ProcessLine(cmd.Data());

  gROOT->ProcessLine("task->DisableV0s();");
  gROOT->ProcessLine("task->DisableTracks();");
  gROOT->ProcessLine("task->DisablePmdClusters();");
  gROOT->ProcessLine("task->DisableCaloClusters();");
  gROOT->ProcessLine("task->DisableCells();");
  gROOT->ProcessLine("task->DisableCaloTrigger(\"PHOS\");");
  gROOT->ProcessLine("task->DisableCaloTrigger(\"EMCAL\");");
  gROOT->ProcessLine("task->DisableHMPID();");
  gROOT->ProcessLine("task->SetPropagateTrackToEMCal(kFALSE);");
  
  if (!mgr->InitAnalysis())
  {
    std::cout << "Could not InitAnalysis" << std::endl;
    return -1;
  }

  // mgr->PrintStatus();
  
  TChain* chain = new TChain("esdTree");
  
  chain->Add(from);
  
  mgr->StartAnalysis("local",chain);
  
  AliInfoGeneral("FILTER_ESD2AODMUON","successfull");

  gSystem->ChangeDirectory(cdir.Data());
  
  gSystem->Rename(Form("%s/AliAOD.Muons.root",odir.Data()),to);

  // gSystem->Rename(Form("%s/AliAOD.root",odir.Data()),Form("%s/AliAOD.root",cdir.Data()));

  gSystem->Unlink(Form("%s/AliAOD.root",odir.Data()));
  gSystem->Unlink(Form("%s/AnalysisResults.root",odir.Data()));

  gSystem->Unlink(odir.Data());
  
  return 0;
}

}