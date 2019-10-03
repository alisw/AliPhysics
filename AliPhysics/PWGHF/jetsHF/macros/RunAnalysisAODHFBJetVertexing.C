/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: $ */

///////////////////////////////////////////////////////////////////
//                                                               //
// Macro to run Bjets Analysis Task on GRID                      //
//                                                               //
// Origin:                                                       //
// Sarah LaPointe                                                //
// Elena Bruna                                                   //
//                                                               //
// Modified:                                                     //
// Y. Corrales Morales, Torino, corrales@to.infn.it              //
//                                                               //
// Last change: Apr 20, 2015                                     //
//                                                               //
///////////////////////////////////////////////////////////////////

#if !defined(__CINT__) || defined(__MAKECINT__)
/* If you want to compile the macro put here all include                  */
#include <TSystem.h>
#include <TObject.h>
#include <TROOT.h>
#include <TGrid.h>

#include <TChain.h>

#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisGrid.h"
#include "AliAnalysisAlien.h"
/*                                                                        */
#endif

void LoadAllLibraries(TString &, Bool_t);
void AddWagonsTask(const Bool_t bIsMC);

AliAnalysisGrid *CreateAlienHandler(TString pluginmode = "test", Bool_t useParFiles = kFALSE);

//
//_____________________________________________________________________________
//
void RunAnalysisAODHFBJetVertexing()
{
  gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS -I$ALICE_PHYSICS/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS  -I$ALICE_PHYSICS/OADB -I$ALICE_PHYSICS/PWGHF -I$ALICE_PHYSICS/PWGHF/base -I$ALICE_PHYSICS/PWGHF/vertexingHF -I$ALICE_PHYSICS/PWG/FLOW/Base -I$ALICE_PHYSICS/PWG/FLOW/Tasks -I-I$ALICE_PHYSICS/JETAN -I$ALICE_PHYSICS/PWG/Tools -g");
  
  TString trainName     = "HFCJ";   //
  TString gridMode      = "test";   //stage
  TString loadMacroPath = "$ALICE_PHYSICS/PWGHF/vertexingHF/macros/";
  
  Bool_t useParFiles = kFALSE;
  
  LoadAllLibraries(loadMacroPath, useParFiles);

//  gROOT->ProcessLine(".L AliRDHFJetsCuts.cxx++g");
//  gROOT->ProcessLine(".L AliRDHFJetsCutsVertex.cxx++g");
//  gROOT->ProcessLine(".L AliHFJetsTagging.cxx++g");
//  gROOT->ProcessLine(".L AliHFJetsTaggingVertex.cxx++g");
//  gROOT->ProcessLine(".L AliHFJetsContainer.cxx++g");
//  gROOT->ProcessLine(".L AliHFJetsContainerVertex.cxx++g");
//  gROOT->ProcessLine(".L AliAnalysisTaskEmcalJetBtagSV.cxx++g");
  
  // Create the analysis manager
  AliAnalysisManager *mgr  = new AliAnalysisManager("My Manager", "My Manager");
  mgr->SetDebugLevel(10);

  // Create Alien plugin, if requested
  TGrid::Connect("alien://");
  AliAnalysisGrid *alienHandler = CreateAlienHandler(gridMode, useParFiles);
  if (!alienHandler) return;
  
  // Connect plug-in to the analysis manager
  mgr->SetGridHandler(alienHandler);
  
  AliAODInputHandler *inputHandler = new AliAODInputHandler("handler", "handler for HFCJ");
  //  inputHandler->AddFriend("./AliAOD.Jets.root");
  inputHandler->SetNeedField();
  mgr->SetInputEventHandler(inputHandler);
  
  Bool_t bIsMC = kTRUE;
  Bool_t doBkgRej = kTRUE;
  AddWagonsTask(bIsMC, doBkgRej);

  //-------------------------------------------------------------------
  //
  // Run the analysis
  //

  if ( !mgr->InitAnalysis() ) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("grid");

  return;
}
//
//_____________________________________________________________________________
//
void LoadAllLibraries(TString &loadMacroPath, Bool_t useParFiles) {
  
  TString loadLibraries="LoadLibraries.C";
  loadLibraries.Prepend(loadMacroPath.Data());
  gROOT->LoadMacro(loadLibraries.Data());
  LoadLibraries(useParFiles);
  
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
  gROOT->ProcessLine(".include $FASTJET/include");
  
  gSystem->Load("/usr/local/lib/libCGAL");
  gSystem->Load("libfastjet");
  gSystem->Load("libsiscone");
  gSystem->Load("libsiscone_spherical");
  gSystem->Load("libfastjetplugins");
  gSystem->Load("libfastjetcontribfragile");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libJETAN");
  gSystem->Load("libFASTJETAN");
  gSystem->Load("libCDB");
  gSystem->Load("libTENDER");
  gSystem->Load("libPWGTools");
  gSystem->Load("libPWGEMCAL");
  gSystem->Load("libPWGmuon");
  gSystem->Load("libPWGJE");
  gSystem->Load("libPWGJEEMCALJetTasks");
  
  return;
}
//
//_____________________________________________________________________________
//
AliAnalysisGrid *CreateAlienHandler(TString pluginmode, Bool_t useParFiles)
{
  // Check if user has a valid token, otherwise make one. This has limitations.
  // One can always follow the standard procedure of calling alien-token-init then
  // source /tmp/gclient_env_$UID in the current shell.
  
  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  
  // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
  plugin->SetRunMode(pluginmode.Data());
  plugin->SetUser("ycorrale");
  
  // Set versions of used packages
  plugin->SetAPIVersion("V1.1x");
  
  plugin->SetROOTVersion("v5-34-08-7");
  plugin->SetAliROOTVersion("v5-06-07");
  plugin->SetAliPhysicsVersion("vAN-20150305");
  plugin->AddExternalPackage("boost::v1_53_0");
  plugin->AddExternalPackage("cgal::v4.4");
  plugin->AddExternalPackage("fastjet::v3.0.6_1.012");
  
  plugin->SetNtestFiles(3);
  // gROOT->LoadMacro("$ALICE_ROOT/PWGHF/vertexingHF/AddGoodRuns.C");
  
  // Declare input data to be processed.
  //************************************************
  // Set data search pattern for DATA
  //************************************************
  
  //************************************************
  // Set data search pattern for MONTECARLO
  //************************************************
  plugin->SetGridDataDir("/alice/sim/LHC11a1b"); // specify MC/data sample /alice/sim/2014/LHC14g3a
  plugin->SetDataPattern("*/AliAOD.root");
  // plugin->SetRunPrefix("000");   // real data
  
  // [9] LHC14g3a large number of events - 195677, 195673, 195644, 195593, 195568, 195566, 195567, 195483, 195479
  
  // [14] LHC14g3a small number of events - 195675, 195635, 195633, 195596, 195592, 195531, 195529, 195482, 195480, 195481, 195389, 195390, 195391, 195478
  
  // small sample plus 4 runs with a larger no. of events
  //  const Int_t nRuns=14;
  //  Int_t runlist[nRuns] = {195675, 195635, 195633, 195596, 195592, 195531, 195529, 195482, 195480, 195481, 195389, 195390, 195391, 195478};
  
  const Int_t nRuns = 1;
  Int_t runlist[nRuns] = {126007};
  
  for (Int_t k = 0; k < nRuns; k++) {
    //if (runlist[k] < firstrun || runlist[k] > lastrun) continue;
    plugin->AddRunNumber(runlist[k]);
  }
  plugin->SetNrunsPerMaster(1);
  plugin->SetOutputToRunNo(1);
  
  plugin->SetGridWorkingDir("PWGLH-HFCJ_Bjets/TestTask");

  // Name of executable
  plugin->SetExecutable("BTag.sh");

  // Declare alien output directory. Relative to working directory.
  plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output
  
  // Declare the analysis source files names separated by blancs. To be compiled runtime
  // using ACLiC on the worker nodes.
  //plugin->SetAnalysisSource("AliAnalysisTaskEmcalJetBtagSV.cxx+");
  
  plugin->SetAdditionalLibs("/usr/local/lib/libCGAL.so libfastjet.so libsiscone.so libsiscone_spherical.so libfastjetplugins.so libfastjettools.so libfastjetcontribfragile.so libGeom.so libGui.so libJETAN.so libFASTJETAN.so libTENDER.so libCORRFW.so libPWGTools.so libCDB.so libProof.so libRAWDatabase.so  libSTEER.so libJETAN.so libFASTJETAN.so libPWGEMCAL.so libEMCALUtils.so libPWGJEEMCALJetTasks.so libPWGflowBase.so libPWGflowTasks.so libPWGHFbase.so libPWGHFvertexingHF.so libANALYSISalice.so");
  
  plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS -I$ALICE_PHYSICS/include -I$FASTJET/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS  -I$ALICE_PHYSICS/OADB -I$ALICE_ROOT/EMCAL  -I$ALICE_PHYSICS/PWGJE -I$ALICE_PHYSICS/PWGJE/EMCALJetTasks -I$ALICE_ROOT/PWGHF -I$ALICE_ROOT/PWGHF/base -I$ALICE_ROOT/PWGHF/vertexingHF -I$ALICE_ROOT/PWG/FLOW/Base -I$ALICE_ROOT/PWG/FLOW/Tasks -I$ALICE_ROOT/PWGJE -I$ALICE_ROOT/JETAN -g");

  plugin->SetDefaultOutputs(kTRUE);
  
  // merging via jdl
  plugin->SetMergeViaJDL(kFALSE);
  plugin->SetOneStageMerging(kFALSE);
  plugin->SetMaxMergeStages(2);
  plugin->SetSplitMaxInputFileNumber(25);
  
  plugin->SetAnalysisMacro("BtagAnalysis.C");
  plugin->SetJDLName("TaskHF.jdl");
  
  return plugin;
}
//
//_____________________________________________________________________________
//
void AddWagonsTask(const Bool_t bIsMC, const Bool_t doBkgRej) {
  // Analysis tasks (wagons of the train)
  //
  //  gROOT->ProcessLine(".L AliRDHFJetsCuts.cxx+g");
  //  gROOT->ProcessLine(".L AliRDHFJetsCutsVertex.cxx+g");
  //  gROOT->ProcessLine(".L AliHFJetsTagging.cxx+g");
  //  gROOT->ProcessLine(".L AliHFJetsTaggingVertex.cxx+g");
  //  gROOT->ProcessLine(".L AliHFJetsContainer.cxx+g");
  //  gROOT->ProcessLine(".L AliHFJetsContainerVertex.cxx+g");
  //  gROOT->ProcessLine(".L AliAnalysisTaskEmcalJetBtagSV.cxx+g");
  
  UInt_t iPhysSelectionFlag = AliVEvent::kAny;
  const char *dataType = "AOD";
  
  TString kTracksInName = "HybridTracks";

  TString kTracksName      = "PicoTracks";
  TString kEmcalTracksName = "EmcalTracks_";
  kEmcalTracksName += kTracksName;
  
  TString kClusName      = "EmcCaloClusters";
  TString kEmcalClusName = "EmcalClusters_";
  kEmcalClusName += kClusName;
  
  TString kCorrClusName = "CaloClustersCorr";
  TString kMCTracksName = "MCParticlesSelected";

  TString kRhoName   = "Rho";
  TString kRhoNameMC = "RhoMC";
  //const char *runPeriod  = "lhc13c";
  
  // EMCal framework
  // ################# Now: Add some basic tasks
  // Physics selection task
  if ( !bIsMC ) {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalPhysicsSelection.C");
    AliVEvent::EOfflineTriggerTypes pSel = AliVEvent::kAny;
    AliPhysicsSelectionTask *physSelTask = AddTaskEmcalPhysicsSelection(kTRUE, kTRUE, pSel, 0., 0., 10, kTRUE, -1, -1, -1, -1);
    if ( !physSelTask ) {
      ::Error("AddWagonsTask", "no physSelTask but running on data");
      return;
    }
  }
  
  // Setup task
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalSetup.C");
  AliEmcalSetupTask *setupTask = AddTaskEmcalSetup();
  
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/jetsHF/macros/AddTaskJetPreparationHF.C");
  AliAnalysisTaskSE *clusm = AddTaskJetPreparationHF("lhc10b" /*Period*/,
                                                     "PicoTracks" ,
                                                     bIsMC /*isMC*/,
                                                     "MCParticlesSelected",
                                                     "",
                                                     "",
                                                     2., 0., 0.03, 0.015, 0.15, 0 /*AliVEvent::kAny*/,
                                                     kFALSE,
                                                     kFALSE,
                                                     kTRUE,
                                                     kFALSE,
                                                     kFALSE,
                                                     1.,
                                                     kFALSE,
                                                     kFALSE
                                                     );
  // #### ADD NECESSARY JET FINDER TASKS
  enum AlgoType {kKT, kANTIKT};
  enum JetType  {kFULLJETS, kCHARGEDJETS, kNEUTRALJETS};
  
  // ################# Now: Add jet finders+analyzers
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");
  AliEmcalJetTask *jetFinderTask   = AddTaskEmcalJet("PicoTracks", "", kANTIKT, 0.4, kCHARGEDJETS, 0.150, 0.300);

  if ( bIsMC ) {
    AliEmcalJetTask *jetFinderTaskMC = AddTaskEmcalJet("MCParticlesSelected", "", kANTIKT, 0.4, kCHARGEDJETS, 0.150, 0.300);
    
    //    gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskJetResponseMaker.C");
    //    TString kJetAkTName   = jetFinderTask->GetJetsName();
    //    TString kJetAkTNameMC = jetFinderTaskMC->GetJetsName();
    
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJetTagger.C");
    AliAnalysisTaskEmcalJetTagger *tagger = AddTaskEmcalJetTagger(jetFinderTask->GetName(),
                                                                  jetFinderTaskMC->GetName(),
                                                                  0.4, "", "", "PicoTracks", "", 0,
                                                                  "V0M", AliVEvent::kAny, "", "");
    tagger->SetNCentBins(1);
    tagger->SetForceBeamType(AliAnalysisTaskEmcal::kpA);
    tagger->SetIsPythia(kTRUE);
    tagger->SetJetTaggingType(AliAnalysisTaskEmcalJetTagger::kClosest);
    tagger->SetJetTaggingMethod(AliAnalysisTaskEmcalJetTagger::kGeo);
  }
  
  //IF DO BACKGROUND REJECTION
  if (doBkgRej) {
    TString myRhoName("ExternalRhoTask");
    AlgoType aType = kKT; //or kANTIKT
    AliEmcalJetTask *jetFinderRhoKT = AddTaskEmcalJet("PicoTracks",
                                                      "",
                                                      aType,
                                                      0.4,
                                                      kCHARGEDJETS,
                                                      0.150,
                                                      0.300,
                                                      0.005,
                                                      1,
                                                      "Jet",
                                                      0., 0, 0);
    jetFinderRhoKT->SetMinJetPt(0);
    
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskRhoSparse.C");
    AliAnalysisTaskRhoSparse* rhotask = AddTaskRhoSparse(jetFinderRhoKT->GetName(),
                                                         jetFinderTask->GetName(),
                                                         "PicoTracks",   //pico trakcs
                                                         "",   //calo clusters
                                                         myRhoName.Data(),
                                                         0.4,  //jet radius
                                                         "TPC",  //cut type
                                                         0.,  //jet area cut
                                                         0.,  //jet pt cut ????????
                                                         0,  //enareacut
                                                         0,  //sfunc
                                                         2,  //excl Jets
                                                         kFALSE,   //no histo
                                                         myRhoName.Data(),  //task name
                                                         kTRUE); //calculate rho CMS
  
    if ( bIsMC ) {
      TString myRhoNameMC("ExternalRhoTaskMC");
      AliEmcalJetTask *jetFinderRhoKTMC = AddTaskEmcalJet("MCParticlesSelected",
                                                          "",
                                                          aType,
                                                          0.4,
                                                          kCHARGEDJETS,
                                                          0.150,
                                                          0.300,
                                                          0.005,
                                                          1,
                                                          "Jet",
                                                          0., 0, 0);
      jetFinderRhoKTMC->SetMinJetPt(0);
      
      AliAnalysisTaskRhoSparse *rhotaskMC = AddTaskRhoSparse(jetFinderRhoKTMC->GetName(),
                                                             jetFinderTaskMC->GetName(),
                                                             "MCParticlesSelected",   //pico trakcs
                                                             "",   //calo clusters
                                                             myRhoNameMC.Data(),
                                                             0.4,  //jet radius
                                                             "TPC",  //cut type
                                                             0.,  //jet area cut
                                                             0.,  //jet pt cut ????????
                                                             0,  //enareacut
                                                             0,  //sfunc
                                                             2,  //excl Jets
                                                             kFALSE,   //no histo
                                                             myRhoNameMC.Data(),  //task name
                                                             kTRUE); //claculate rho CMS
      
    } //isMC
  }//doBkgRej
  TString taskName;
  taskName="AddTaskEmcalJetBtagSV.C";
  gROOT->LoadMacro(taskName.Data());
  
  Int_t taskHardPtCut[4] = {10, 18, 30, 50};
  
  Int_t nBJetTask = 1;
  for (Int_t iTask = 0; iTask < nBJetTask; ++iTask) {
    TString ptHardName = Form("%d", taskHardPtCut[iTask]);
    
    AliAnalysisTaskEmcalJetBtagSV *taskJet = AddTaskEmcalJetBtagSV("PicoTracks",
                                                                   jetFinderTask->GetName(),
                                                                   "MCParticlesSelected",
                                                                   jetFinderTaskMC->GetName(),
                                                                   0.4 /*R*/,
                                                                   "TPC",
                                                                   "standard",
                                                                   bIsMC /*IsCorrMode*/,
                                                                   kTRUE /*PtHard*/,
                                                                   kTRUE /*DoBkgRej*/,
                                                                   ptHardName.Data(),
                                                                   "bbbar",
                                                                   0.4,
                                                                   "",
                                                                   0.,
                                                                   100.,
                                                                   "HFjetsContainer",
                                                                   kFALSE);

    if (doBkgRej) {
      taskJet->SetExternalRhoTaskName(myRhoName.Data());
      if (bIsMC) taskJet->SetMcExternalRhoTaskName(myRhoNameMC.Data());
    }
  }
  
  return;
}