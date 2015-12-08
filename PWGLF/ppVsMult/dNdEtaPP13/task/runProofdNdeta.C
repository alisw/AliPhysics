#if !defined (__CINT__) || (defined(__MAKECINT__))
#include <iostream>
#include "AliAnalysisGrid.h"
#include "TSystem.h"
#include "TROOT.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisGrid.h"
#include "AliVEventHandler.h"
#include "AliAODInputHandler.h"
#include "AliPhysicsSelection.h"
#include "AliAnalysisAlien.h"
#include "AliMCEventHandler.h"
#include "AliPhysicsSelectionTask.h"
#include "TRegexp.h"
#include "TProof.h"
#include "AliESDInputHandler.h"
#include "AliOADBPhysicsSelection.h"
#include "AliMultSelectionTask.h"

#endif

// Parameters of the add macro, listed here for reference
// Some of those should be moved as args for the run macro itself

// TODO: make the following settable?

const Char_t *listname = "clist"; 
//
Float_t etaMin     =-1;          // min eta range to fill in histos
Float_t etaMax     = 1;           // max eta range to fill in histos
Float_t zMin       = -7;         // process events with Z vertex min
Float_t zMax       =  7;          //                     max positions
const char* useCentVar = "V0M";         // centrality variable to use
//
Float_t cutSigNStd  = 1.5;         // cut on weighed distance used to extract signal
Float_t cutSigDPhiS = -1;         // cut on dPhi-phiBent used to extract signal (if negative -> dphi*sqrt(cutSigNStd)
//Bool_t  useMC  = kTRUE;         // fill MC info
//
Bool_t doRec  = kFALSE;            // fill data histos from new reco
Bool_t doInj  = kFALSE;            // create Inj. bg
Bool_t doRot  = kFALSE;            // create Rot. bg
//
                          // specific parameters for reconstruction
float  phiRot      = 3.14159e+00; // angle for bg. generation with rotation
float  injScale    = 1.;//0.7,    // inject injScale*Ncl(Lr1/Lr2) hits
Bool_t scaleDTheta = kTRUE;        // scale dTheta by 1/sin^2(theta) in trackleting
float  nStdDev     = 25.;         // number of st.dev. for tracklet cut to keep
float  dphi        = 0.06;        // dphi window (sigma of tracklet cut)
float  dtht        = 0.025;       // dtheta .... (if negative, abs will be used with additional cut on |dthetaX|, apart from w.distance
float  phishift    = 0.0045;      // bending shift
Bool_t remOvl      = kTRUE;     
float  ovlPhiCut   = 0.005;
float  ovlZetaCut  = 0.05;
Bool_t checkReconstructables = kFALSE;//kTRUE, // fill histos for reconstructable (needs useMC and doRec)

UInt_t trigSel = AliVEvent::kMB;


//______________________________________________________________________________

void runProofdNdeta(
                  TString dataset = "Find;"              
                  "BasePath=/alice/data/2015/LHC15f/000225768/pass1/%.%/;"
                  "FileName=root_archive.zip;"
                  "Anchor=AliESDs.root;"
                  "Tree=/esdTree;"
                  "Mode=remote;",

  // TString dataset = "Find;"
  //                   "BasePath=/alice/data/2015/LHC15f/000225106/cpass1_pass1/%.%/;"
  //                   "FileName=root_archive.zip;"
  //                   "Anchor=AliESDs_Barrel.root;"
  //                   "Tree=/esdTree;",
  // TString dataset = "Find;"
  //                   "BasePath=/alice/data/2013/LHC13e/000195949/ESDs/muon_pass2/AOD134/%/;"
  //                   "FileName=root_archive.zip;"
  //                   "Anchor=AliAOD.root;"
  //                   "Tree=/aodTree;",
                  Bool_t usePhysicsSelection = kFALSE,
                  Bool_t isMC = 0,
                  Int_t numEvents = 999999999,
                  Int_t firstEvent = 0,
                  const char * oadbMultSel = "LHC15f",                  
                  const char * addTaskString = "AddAnalysisTaskdNdEtaPP13(\"outfile.root\",\"%s\",%f,%f,%f,%f,\"%s\",%f,%f,%d,%d,%d,%d,%f,%f,%d,%f,%f,%f,%f,%d,%f,%f,%d,%u)"
                  
) {

  // VAF cheatsheet
  // vaf-enter
  // vafctl start
  // vafreq <NUM_OF_WORKERS>
  // vafcount
  // vafctl stop #release resources
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");


  TList *list = new TList(); 
  list->Add(new TNamed("ALIROOT_ENABLE_ALIEN", "1"));  // important: creates token on every PROOF worker

  TProof::Open("pod://");
  gProof->AddIncludePath("-I$ALICE_ROOT/include");
  gProof->AddIncludePath("-I$ALICE_PHYSICS/include");


  
  // Check the dataset before running the analysis!
  gProof->ShowDataSet( dataset.Data() );
  //return;  // <-- uncomment this to test search before running the analysis!


  // A single AliRoot package for *all* AliRoot versions: new on VAF
  TString aliceVafPar = "/afs/cern.ch/alice/offline/vaf/AliceVaf.par";
  gProof->UploadPackage(aliceVafPar.Data());
  gProof->EnablePackage(aliceVafPar.Data(), list);  // this "list" is the same as always

  AliAnalysisManager *mgr  = new AliAnalysisManager("Analysis Train");
  AliVEventHandler* iH = new AliESDInputHandler(); // WARNING: this needs to be changed for AOD!
  //    AliAODInputHandler* iH = new AliAODInputHandler();
  mgr->SetInputEventHandler(iH);
  //AliESDInputHandler *esdInputHandler = dynamic_cast<AliESDInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  //  if (!esdInputHandler) esdInputHandler = new AliESDInputHandlerRP();
  //mgr->SetInputEventHandler(esdInputHandler);

  if(isMC) {
    AliMCEventHandler* mchandler = new AliMCEventHandler();
    // Not reading track references
    mchandler->SetReadTR(kFALSE);
    mgr->SetMCtruthEventHandler(mchandler);
  }   
 // PHYSICS SELECTION
  if(usePhysicsSelection){
    std::cout << "WARNING! Custom Phuscs Selection" << std::endl;
    AliOADBPhysicsSelection *customPS = new AliOADBPhysicsSelection("customPS");
    customPS->AddCollisionTriggerClass(AliVEvent::EOfflineTriggerTypes(trigSel),"+CINT7-B-NOPF-ALLNOTRD","B",0);
    customPS->SetHardwareTrigger(0, "1");
    customPS->SetOfflineTrigger(0, "1");
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
    TString physSelMacro = ""; // comment for AODs"
    physSelMacro.Form("AddTaskPhysicsSelection(%d)",isMC);
    AliPhysicsSelectionTask *physSelTask = (AliPhysicsSelectionTask*)gROOT->ProcessLine( physSelMacro.Data() ); 
    physSelTask->GetPhysicsSelection()->SetCustomOADBObjects(customPS,0);
  }

  // MULT SELECTION
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
  //AddTask: Should take care of everything transparently...
    AliMultSelectionTask *taskMS = (AliMultSelectionTask*)gROOT->ProcessLine( "AddTaskMultSelection();");
  //  AliMultSelectionTask *taskMS = (AliMultSelectionTask*)gROOT->ProcessLine( "AddTaskMultSelection(kTRUE);"); // produces tree for centrality calibration with the true option
    
  //User Case FIXME??
  taskMS->SetAddInfo(kTRUE);

  taskMS->SetSelectedTriggerClass(AliVEvent::EOfflineTriggerTypes(trigSel));
  taskMS -> SetAlternateOADBforEstimators ( oadbMultSel );


 
  // create task

  // Load custom Tast; TODO: to be removed once this gets compiled with aliroot
  const char* loadTaskStr[] = {"AliITSMultRecBg.cxx+g", "AliAnalysisTaskdNdEtapp13.cxx+g","AddAnalysisTaskdNdEtaPP13.C+g",0};// FIXME: this assumes that the macro name will not change
  Int_t itask = -1;
  while (loadTaskStr[++itask]){
    if (gProof != NULL) {
      gProof->Load(loadTaskStr[itask]);
    }
    else {
      gROOT->LoadMacro(loadTaskStr[itask]); 
    }
  }
  AliAnalysisTaskSE *task;
  {
    TString buf1, buf2;
    buf1.Form("%s", addTaskString);
    buf2.Form(buf1.Data(), listname,etaMin, etaMax, zMin, zMax, useCentVar, cutSigNStd, cutSigDPhiS, isMC, doRec, doInj, doRot, phiRot, injScale, scaleDTheta, nStdDev, dphi, dtht, phishift, remOvl, ovlPhiCut, ovlZetaCut, checkReconstructables, trigSel);
    std::cout << "Add macro: " << buf2.Data() << std::endl;
    
    task = (AliAnalysisTaskSE *)gROOT->ProcessLine( buf2.Data() );
  }

  task->Print();
  if (!mgr->InitAnalysis()) return;
  mgr->Print();

  mgr->StartAnalysis("proof", dataset, numEvents, firstEvent);

}

