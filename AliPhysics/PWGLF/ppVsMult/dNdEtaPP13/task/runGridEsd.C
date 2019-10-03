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
#include "TGrid.h"
#include "AliMultSelectionTask.h"

#endif

// TString dataDir = "/alice/data/2015/LHC15j";
// //TString dataPattern = "/muon_calo_pass1/AOD/*/AliAOD.root";
// TString dataPattern = "/muon_calo_pass1/*/AliESDs.root";
// TString workingDir = "lhc15j_spdscan";

// // SPD scan
// Int_t nRuns=8;
// //Int_t runList[] = {238396,238399};
// Int_t runList[] = {238429,238431,238432,238395,238396,238398,238399,238400};

// // LHC15f
// TString dataDir = "/alice/data/2015/LHC15f";
// //TString dataPattern = "/muon_calo_pass1/AOD/*/AliAOD.root";
// TString dataPattern = "/pass2/*/AliESDs.root";
// TString workingDir = "dNdeta_lhc15f";
// // Int_t nRuns = 56;
// // Int_t runList[] = {226500, 226495, 226483, 226476, 226472, 226468, 226466, 226452, 226445, 226444, 226225, 226220, 226170, 226062, 225768, 225766, 225763, 225762, 225757, 225753, 225719, 225717, 225716, 225710, 225709, 225708, 225707, 225705, 225587, 225586, 225582, 225580, 225579, 225578, 225576, 225322, 225315, 225314, 225313, 225310, 225309, 225307, 225305, 225106, 225052, 225051, 225050, 225043, 225041, 225037, 225035, 225031, 225026, 225016, 225011, 225000};

// // LHC15f MC
// TString dataDir = "/alice/sim/2015/LHC15g3c2/";
// //TString dataPattern = "/muon_calo_pass1/AOD/*/AliAOD.root";
// TString dataPattern = "/*/AliESDs.root";
// TString workingDir = "dNdeta_LHC15g3c2_MC";



// Int_t nRuns = 79;
// Int_t runList[] = {
//     235683, 235684, 235685, 235687, 235694, 235710, 235714, 235721, 235759, 235811,
//     235813, 235839, 235841, 235892, 235893, 235895, 235896, 235897, 235898, 236062,
//     236137, 236138, 236150, 236151, 236153, 236158, 236159, 236161, 236163, 236164,
//     236203, 236204, 236222, 236227, 236234, 236238, 236240, 236242, 236244, 236246,
//     236248, 236281, 236284, 236285, 236331, 236334, 236337, 236348, 236349, 236352,
//     236353, 236354, 236356, 236357, 236359, 236360, 236386, 236389, 236393, 236395,
//     236397, 236441, 236443, 236444, 236446, 236453, 236459, 236462, 236541, 236554,
//     236556, 236558, 236562, 236563, 236564, 236565, 236569, 236575, 236866
// };
//TString dataDir = "/alice/data/2015/LHC15i";
//TString workingDir = "lhc15i";

// TString dataDir = "/alice/data/2015/LHC15j";
// TString workingDir = "lhc15j";

//Int_t nRuns = 79;
//
//Int_t runList[] = {
//    235683, 235684, 235685, 235687, 235694, 235710, 235714, 235721, 235759, 235811,
//    235813, 235839, 235841, 235892, 235893, 235895, 235896, 235897, 235898, 236062,
//    236137, 236138, 236150, 236151, 236153, 236158, 236159, 236161, 236163, 236164,
//    236203, 236204, 236222, 236227, 236234, 236238, 236240, 236242, 236244, 236246,
//    236248, 236281, 236284, 236285, 236331, 236334, 236337, 236348, 236349, 236352,
//    236353, 236354, 236356, 236357, 236359, 236360, 236386, 236389, 236393, 236395,
//    236397, 236441, 236443, 236444, 236446, 236453, 236459, 236462, 236541, 236554,
//    236556, 236558, 236562, 236563, 236564, 236565, 236569, 236575, 236866
//};

// Int_t nRuns = 24;

// Int_t runList[] = {
//     237029, 237645, 237670, 237671, 237675, 237676, 237678, 237707, 237708, 237710, 
//     237711, 237713, 237765, 237777, 237779, 237780, 237782, 237795, /*237844, 237845, 
//     237847,*/ 237945, 237948, 237969, 237978, 237982, 237983
// };
// TODO: clean up comments in this macro

TString aliPhysicsVersion="vAN-20151129-1";

// TODO: make the following settable?
Int_t * runList;
Int_t nRuns = 1;

const Char_t *listname = "clist"; 
//
Float_t etaMin     =-1;          // min eta range to fill in histos
Float_t etaMax     = 1;           // max eta range to fill in histos
Float_t zMin       = -7;         // process events with Z vertex min
Float_t zMax       =  7;          //                     max positions

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

AliVEvent::EOfflineTriggerTypes trigSel = AliVEvent::kMB;

void runGridEsd(TString dataDir = "/alice/sim/2015/LHC15g3c2/",
                TString strRunList = "226062",
                TString dataPattern = "/*/AliESDs.root",
                TString workingDir = "dNdeta_LHC15g3c2_MC",
                Bool_t usePhysicsSelection = kFALSE,
                Bool_t fIsMC = 0,
                Bool_t doMultSelTrees = kFALSE,
                const char * oadbMultSel = "LHC15f",
                const char * gridMode = "full",
                const char* useCentVar = "V0M",
                const char * addTaskString = "AddAnalysisTaskdNdEtaPP13(\"outfile.root\",\"%s\",%f,%f,%f,%f,\"%s\",%f,%f,%d,%d,%d,%d,%f,%f,%d,%f,%f,%f,%f,%d,%f,%f,%d,%u)"
                ){
  gSystem->AddIncludePath("-I. -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");

  // Build runList
  TObjArray * tokens = strRunList.Tokenize(",");
  nRuns = tokens->GetEntries();
  runList = new Int_t[nRuns];
  TIter iter (tokens);
  TObjString * runToken = 0;
  Int_t irun=0;
  while((runToken = (TObjString*)iter.Next())){
    runList[irun] = runToken->String().Atoi();
    std::cout << runList[irun] << "; ";
    
    irun++;
  }
  std::cout << "" << std::endl;
  delete tokens;
  

  
  if (!TGrid::Connect("alien://")) return;

  AliAnalysisManager *mgr = new AliAnalysisManager("Analysis");
//  AliAODInputHandler* handler = new AliAODInputHandler();
  AliESDInputHandler* handler = new AliESDInputHandler();
  handler->SetReadFriends(kTRUE);
  mgr->SetInputEventHandler(handler);

  if(fIsMC) {
    AliMCEventHandler* mchandler = new AliMCEventHandler();
    // Not reading track references
    mchandler->SetReadTR(kFALSE);
    mgr->SetMCtruthEventHandler(mchandler);
  }   

  
  // PHYSICS SELECTION
  if(usePhysicsSelection){
    std::cout << "WARNING! Custom Physics Selection" << std::endl;
    AliOADBPhysicsSelection *customPS = new AliOADBPhysicsSelection("customPS");
    customPS->AddCollisionTriggerClass(trigSel,"+CINT7-B-NOPF-ALLNOTRD","B",0);
    customPS->SetHardwareTrigger(0, "1");
    customPS->SetOfflineTrigger(0, "1");
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
    TString physSelMacro = ""; // comment for AODs"
    physSelMacro.Form("AddTaskPhysicsSelection(%d)",fIsMC);
    AliPhysicsSelectionTask *physSelTask = (AliPhysicsSelectionTask*)gROOT->ProcessLine( physSelMacro.Data() ); 
    physSelTask->GetPhysicsSelection()->SetCustomOADBObjects(customPS,0);
  }


  // MULT SELECTION
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
  //AddTask: Should take care of everything transparently...
  //  AliMultSelectionTask *taskMS = (AliMultSelectionTask*)gROOT->ProcessLine( "AddTaskMultSelection(kTRUE);");
  TString multSelMacro = ""; // comment for AODs"
  multSelMacro.Form("AddTaskMultSelection(%d)",doMultSelTrees);
  AliMultSelectionTask *taskMS = (AliMultSelectionTask*)gROOT->ProcessLine(multSelMacro.Data());// With the true option it saves the trees for calibration

  //  alien:///Users/mfloris/Work/ALICE/ANALYSIS/current/HMTF/dNdeta/task/
  if(fIsMC) {
    if(TString(oadbMultSel).EndsWith(".root")){
      // use micro-OADB calibrated on a specific MC production
      taskMS -> SetAlternateOADBFullManualBypassMC(oadbMultSel);
    }
    else {
      // Use boundaries as in data
      taskMS -> SetAlternateOADBforEstimators ( oadbMultSel );
    }
  }
  //  taskMS -> SetAlternateOADBFullManualBypassMC("alien:///alice/cern.ch/user/m/mfloris/dNdeta13TeV/OADB-LHC15g3c2_plus.root");
  taskMS->SetSelectedTriggerClass(trigSel);
    //taskMS -> SetAlternateOADBFullManualBypassMC("LHC15f");
  //User Case FIXME??
  taskMS->SetAddInfo(kTRUE);
  //  taskMS->SetSaveCalibInfo(kTRUE); //cross-check information for debugging

  // MY TASK
  const TString loadTaskStr[] = {"AliITSMultRecBg.cxx", "AliAnalysisTaskdNdEtapp13.cxx","AddAnalysisTaskdNdEtaPP13.C","END"};// FIXME: this assumes that the macro name will not change
  Int_t itask = -1;
  TString listOfGridSource = "";
  TString listOfGridAdditionalLibs ="";
  while (loadTaskStr[++itask]!="END"){
    if (gProof != NULL) {
      gProof->Load(loadTaskStr[itask]+"+g");
    }
    else {
      gROOT->LoadMacro(loadTaskStr[itask]+"+g"); 
    }
    // Prepare the lists of additional libs for the pugin
    if(loadTaskStr[itask].EndsWith(".cxx")){
      listOfGridSource = listOfGridSource + loadTaskStr[itask] + " ";
      listOfGridAdditionalLibs = listOfGridAdditionalLibs + TString(loadTaskStr[itask]).ReplaceAll(".cxx", ".h") + " " + loadTaskStr[itask] + " ";
    }
    
  }
  AliAnalysisTaskSE *task;
  {
    TString buf1, buf2;
    buf1.Form("%s", addTaskString);
    buf2.Form(buf1.Data(), listname,etaMin, etaMax, zMin, zMax, useCentVar, cutSigNStd, cutSigDPhiS, fIsMC, doRec, doInj, doRot, phiRot, injScale, scaleDTheta, nStdDev, dphi, dtht, phishift, remOvl, ovlPhiCut, ovlZetaCut, checkReconstructables, trigSel);
    std::cout << "Add macro: " << buf2.Data() << std::endl;
    
    task = (AliAnalysisTaskSE *)gROOT->ProcessLine( buf2.Data() );
  }

  
  
  // FIXME: Remove this stuff
  // gROOT->LoadMacro("AliITSMultRecBg.cxx+");
  // gROOT->LoadMacro("AliAnalysisTaskdNdEtaRuben.cxx+");
  // gROOT->LoadMacro("AddAnalysisTaskdNdEtaRuben.C");
  // AddAnalysisTaskdNdEtaRuben();
  //                                                               // "AnalysisResults.root",
  //                                                               //  "clist",
                                                                
  //                                                               //   /*  Float_t etaMin     = */ -1,          // min eta range to fill in histos
  //                                                               // /*Float_t etaMax     = */ 1,          // max eta range to fill in histos
  //                                                               // /*Float_t zMin       = */ -7,         // process events with Z vertex min
  //                                                               // /*Float_t zMax       = */  7,         //                     max positions
  //                                                               // /*const char* useCentVar = */ "V0M",          // centrality variable to use
  //                                                               // //
  //                                                               // /*Float_t cutSigNStd  = */ 1.5,        // cut on weighed distance used to extract signal
  //                                                               // /*Float_t cutSigDPhiS = */ -1,        // cut on dPhi-phiBent used to extract signal (if negative -> dphi*sqrt(cutSigNStd)
  //                                                               // /*Bool_t  useMC  = */ fIsMC,          // fill MC info
  //                                                               // //
  //                                                               // /*Bool_t doRec  = */ kFALSE,          // fill data histos from new reco
  //                                                               // /*Bool_t doInj  = */ kFALSE,          // create Inj. bg
  //                                                               // /*Bool_t doRot  = */ kFALSE,          // create Rot. bg
  //                                                               // //
  //                                                               // // specific parameters for reconstruction
  //                                                               // /*float  phiRot      = */ 3.14159e+00, // angle for bg. generation with rotation
  //                                                               // /*float  injScale    = */ 1.,//0.7,    // inject injScale*Ncl(Lr1/Lr2) hits
  //                                                               // /*Bool_t scaleDTheta = */ kTRUE,       // scale dTheta by 1/sin^2(theta) in trackleting
  //                                                               // /*float  nStdDev     = */ 25.,         // number of st.dev. for tracklet cut to keep
  //                                                               // /*float  dphi        = */ 0.06,        // dphi window (sigma of tracklet cut)
  //                                                               // /*float  dtht        = */ 0.025,       // dtheta .... (if negative, abs will be used with additional cut on |dthetaX|, apart from w.distance
  //                                                               // /*float  phishift    = */ 0.0045,      // bending shift
  //                                                               // /*Bool_t remOvl      = */ kTRUE,
  //                                                               // /*float  ovlPhiCut   = */ 0.005,
  //                                                               // /*float  ovlZetaCut  = */ 0.05,
  //                                                               // /*Bool_t checkReconstructables = */ kFALSE //kTRUE, // fill histos for reconstructable (needs useMC and doRec)
  //                                                               // //
  //                                                               // );

  if (!mgr->InitAnalysis()) return;
  
  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  plugin->SetRunMode(gridMode);
  plugin->AddIncludePath("-I. -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
  plugin->SetNtestFiles(2);
  plugin->SetAPIVersion("V1.1x");
  plugin->SetAliPhysicsVersion(aliPhysicsVersion);
  plugin->SetGridDataDir(dataDir.Data());
  plugin->SetDataPattern(dataPattern.Data());
  plugin->SetGridWorkingDir(workingDir.Data());
  if(!fIsMC) plugin->SetRunPrefix("000");
  for (Int_t i=0;i<nRuns;i++)  plugin->AddRunNumber(runList[i]);
  plugin->SetGridOutputDir("output");
  plugin->SetAnalysisSource(listOfGridSource.Data());
  plugin->SetAdditionalLibs(listOfGridAdditionalLibs.Data());
  plugin->SetNrunsPerMaster(1);
  plugin->SetSplitMaxInputFileNumber(300);
  plugin->SetOutputToRunNo(1);
  mgr->SetGridHandler(plugin);
  mgr->PrintStatus();
  mgr->StartAnalysis("grid");
 
//  TChain *chain = new TChain("esdTree");
//  for (Int_t i=1;i<=1;i++) chain->AddFile(Form("/data/esd/LHC12h/189616/AliESDs.%03i.root",i));
//  mgr->StartAnalysis("local",chain);
}
