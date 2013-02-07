/**
 * @file   AODTrain.C
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Wed Feb  6 23:10:29 2013
 * 
 * @brief  Setup of train like the `official' AODtrain.C script
 * 
 * 
 * @ingroup pwglf_forward_trains_examples
 */
#ifndef __CINT__
# include <AliAnalysisManager.h>
#else 
class AliAnalysisManager;
#endif
#include "TrainSetup.C"

/**
 * Setup of train like the `official' AODtrain.C script
 * 
 * @ingroup pwglf_forward_trains_examples
 */
class AODTrain : public TrainSetup
{
public:
  //------------------------------------------------------------------
  /** 
   * Constructor 
   * 
   * @param name Name of job
   */
  AODTrain(const char* name="myTest") : TrainSetup(name) 
  {
    fOptions.Set("type", "ESD");
    // over all options 
    fOptions.Add("collision", "TYPE", "Collision system", 0);
    fOptions.Add("year", "YEAR", "Year", 2011);
    fOptions.Add("cdb", "Connect to CDB", true);
    fOptions.Add("physics-selection", "Use PhysicsSelection", true);
    fOptions.Add("tender", "Use of Tenders", false);
    fOptions.Add("centrality", "Use centrality", false);
    fOptions.Add("v0-tender", "Use of V0 tender", false);
    fOptions.Add("track-refs", "Also read track references");
    fOptions.Add("kinematics-filter", "Use kinematics fitler");
    fOptions.Add("sys-info", "Use sys info");
    fOptions.Add("run", "RUNNUMBER", "Set initial run number", 0);
    
    // Task options 
    fOptions.Add("esd-filter", "Use ESD fitler", true);
    fOptions.Add("muon-copy", "Make additional MUON specific AOD", true);
    fOptions.Add("jetan", "Use Jet analysis", true);
    fOptions.Add("jetan-delta", "Make Jet analysis delta AOD", true);
    fOptions.Add("vertexing", "Use PWGHF vertexing", true);
    fOptions.Add("jpsi-filter", "Use PWGDQ J/Psi filter", false);
    fOptions.Add("d0-decay", "Use PWGHF D0->h+h", true);
    fOptions.Add("high-pt", "Use high pt", true);
    fOptions.Add("forward-nch", "Forward charged particle", true);
    // Other 
    fOptions.Add("pid-response", "Use PID response", true);
    fOptions.Add("pid-qa", "Do PID QA", true);
  }
  //------------------------------------------------------------------
  /** 
   * Create our tasks 
   * 
   * @param mgr Manager 
   */
  void CreateTasks(AliAnalysisManager* mgr)
  {

    fHelper->LoadLibrary("libCORRFW");

    Bool_t cdb              = fOptions.Has("cdb");
    Bool_t tender           = fOptions.Has("tender");
    Bool_t physicsSelection = fOptions.Has("physics-selection");

    if (cdb || tender) {
      fHelper->LoadLibrary("libCDB");
      fHelper->LoadLibrary("libProof");
      fHelper->LoadLibrary("libRAWDatabase");
      fHelper->LoadLibrary("libSTEER");
      fHelper->LoadLibrary("libSTAT");
      fHelper->LoadLibrary("libTRDbase");
      fHelper->LoadLibrary("libVZERObase");
      fHelper->LoadLibrary("libTPCbase");
      fHelper->LoadLibrary("libITSbase");
      fHelper->LoadLibrary("libHMPIDbase");
      fHelper->LoadLibrary("libTENDER");
      fHelper->LoadLibrary("libTENDERSupplies");
    }

    if (fOptions.Has("sys-info")) mgr->SetNSysInfo(100);
    
    // AliAnalysisTask*   task = 0;
    AliAnalysisTaskSE* seTask = 0;

    AliAnalysisManager::SetCommonFileName("AODQA.root");
    if (tender) 
      AddTask("$ALICE_ROOT/ANALYSIS/TenderSupplies/AddTaskTender.C",
	      Form("%d", fOptions.Has("v0-tender")));
    
    if (fOptions.Has("pid-response")) 
      AddTask("AddTaskPIDResponse.C");

    if (fOptions.Has("pid-qa")) {
      seTask = AddSETask("AddTaskPIDqa.C");
      seTask->SelectCollisionCandidates(AliVEvent::kAny);
    }
    
    if (cdb && !tender) {
      fHelper->LoadLibrary("libRAWDatarec");
      fHelper->LoadLibrary("libTRDrec");
      fHelper->LoadLibrary("libVZEROrec");
      fHelper->LoadLibrary("libTPCrec");
      fHelper->LoadLibrary("libITSrec");
      fHelper->LoadLibrary("libPWGPP");

      // AddTask("$ALICE_ROOT/PWGPP/PilotTrain/AddTaskCDBconnect.C");
      gROOT->LoadMacro("$ALICE_ROOT/PWGPP/PilotTrain/AddTaskCDBconnect.C");
      gROOT->ProcessLine(Form("AddTaskCDBconnect(%d)", fOptions.AsInt("run")));
      gROOT->ProcessLine("AliCDBManager* cdb = AliCDBManager::Instance();"
			 "cdb->SetDefaultStorage(\"raw://\");");
    }
    
    if (fOptions.Has("high-pt"))
      AddTask("AddTaskFilteredTree.C");
  
    if (fOptions.Has("esd-filter")) {
      fHelper->LoadLibrary("PWGmuon");
      if (fOptions.Has("muon-copy")) {
	mgr->RegisterExtraFile("AliAOD.Muons.root");
	mgr->RegisterExtraFile("AliAOD.Dimuons.root");
      }
      Int_t runFlag = (fOptions.AsInt("year") % 100) * 100;
      TString args = 
	TString::Format("%d,%d,false,false,false,false,true,false,false,%d",
			fOptions.Has("kinematics-filter"), 
			fOptions.Has("muon-copy"), runFlag);
      AddTask("AddTaskESDFilter.C", args);
    }
  
    Bool_t vertexing = fOptions.Has("vertexing");
    Bool_t d0Decay   = fOptions.Has("d0-decay");
    if (vertexing || d0Decay) {
      TString src = "$ALICE_ROOT/PWGHF/vertexingHF/ConfigVertexingHF.C";
      if (fOptions.AsInt("collision") == 1) 
	src.ReplaceAll(".C", "highmult.C");
      TFile::Cp(src, "ConfigVertexingHF.C");
      fHelper->LoadAux("ConfigVertexingHF.C");
    }    

    if (vertexing) {
      fHelper->LoadLibrary("PWGflowBase");
      fHelper->LoadLibrary("PWGflowTasks");
      fHelper->LoadLibrary("PWGHFvertexingHF");
      seTask = 
	AddSETask("$ALICE_ROOT/PWGHF/vertexingHF/macros/AddTaskVertexingHF.C");
      seTask->SelectCollisionCandidates(0);
      mgr->RegisterExtraFile("AliAOD.VertexingHF.root");
    }

    if (fOptions.Has("jpsi-filtering")) { 
      fHelper->LoadLibrary("PWGDQdielectron");
      seTask = 
	AddSETask("$ALICE_ROOT/PWGDQ/dielectron/macros/AddTaskJPSIFilter.C");
      seTask->SelectCollisionCandidates(0);
      mgr->RegisterExtraFile("AliAOD.Dielectron.root");
    }
  
    if (d0Decay) 
      AddTask("$ALICE_ROOT/PWGHF/vertexingHF/AddD2HTrain.C",
	      "kFALSE, 1,0,0,0,0,0,0,0,0,0,0");
  
    if (fOptions.Has("jetan")) {
      fHelper->LoadLibrary("JETAN");
      
      Bool_t jetanDelta = fOptions.Has("jetan-delta");      
      if (jetanDelta) {
	fHelper->LoadLibrary("CGAL");
	fHelper->LoadLibrary("fastjet");
	fHelper->LoadLibrary("siscone");
	fHelper->LoadLibrary("SISConePlugin");
	fHelper->LoadLibrary("FASTJETAN");
	// --- For newer fastjet ---
	// fHelper->LoadLibrary("CGAL");
	// fHelper->LoadLibrary("fastjet");
	// fHelper->LoadLibrary("fastjettools");
	// fHelper->LoadLibrary("siscone");
	// fHelper->LoadLibrary("siscone_spherical");
	// fHelper->LoadLibrary("fastjetplugins");
	// fHelper->LoadLibrary("FASTJETAN");
      }
      
      // Write script to do this - avoid useless links against loadable
      // library
      std::ofstream o("AddJets.C");
      o << "void AddJets(AliAnalysisManager*  mgr,\n"
	<< "             UInt_t               highPt,\n"
	<< "             const char*          deltaAOD,\n"
	<< "             Bool_t               useDelta,\n"
	<< "             Bool_t               usePS,\n"
	<< "             Bool_t               pbpb,\n"
	<< "             Float_t              etaCut,\n"
	<< "             Float_t              centLow,\n"
	<< "             Float_t              centHigh)\n"
	<< "{\n"
	<< "  gROOT->SetMacroPath(Form(\"%s:$ALICE_ROOT/PWGJE/macros\",\n"
	<< "                           gROOT->GetMacroPath()));\n"
	<< "  \n"                      
	<< "  gROOT->LoadMacro(\"AddTaskJets.C\");\n"
	<< "  AliAnalysisTaskJets *jTask = 0;\n"
	<< "  jTask = AddTaskJets(\"AOD\",\"UA1\",0.4,highPt,1.,0);\n"
	<< "  jTask->SetNonStdOutputFile(deltaAOD);\n"
	<< "  \n"
	<< "  if (useDelta) {\n"
	<< "    mgr->RegisterExtraFile(deltaAOD);\n"
	<< "    if (pbpb) {\n"
	<< "      jTask = AddTaskJets(\"AOD\",\"UA1\",0.4,highPt,1.,2);\n"
	<< "      jTask->SetNonStdOutputFile(deltaAOD);\n"
	<< "    }\n"
	<< "  }\n"
	<< "  \n";
      o << "  jTask =AddTaskJets(\"AOD\",\"SISCONE\",0.4,highPt,0.15,0);\n"
	<< "  jTask->SetNonStdOutputFile(deltaAOD);\n"
	<< "  TString subBranches = jTask->GetNonStdBranch();\n"
	<< "  \n";
      o << "  gROOT->LoadMacro(\"AddTaskJetCluster.C\");\n"
	<< "  AliAnalysisTaskJetCluster* cTask = 0;\n"
	<< "  cTask = AddTaskJetCluster(\"AOD\",\"\",highPt,usePS,\n"
	<< "                          \"KT\",0.4,0,1,deltaAOD,0.15,etaCut,0);\n"
	<< "  cTask->SetBackgroundCalc(kTRUE);\n"
	<< "  cTask->SetNRandomCones(10);\n"
	<< "  cTask->SetCentralityCut(centLow,centHigh);\n"
	<< "  cTask->SetGhostEtamax(etaCut);\n"
	<< "  TString bgBranch = Form(\"%s_%s\",\n"
	<< "    AliAODJetEventBackground::StdBranchName(),\n"
	<< "    cTask->GetJetOutputBranch());\n"
	<< "  \n";
      o << "  cTask = AddTaskJetCluster(\"AOD\",\"\",highPt,usePS,\n"
	<< "                          \"ANTIKT\",0.4,2,1,deltaAOD,0.15);\n"
	<< "  cTask->SetNRandomCones(10);\n"
	<< "  cTask->SetCentralityCut(centLow,centHigh);\n"
	<< "  if (pbpb) cTask->SetBackgroundBranch(bgBranch);\n"
	<< "  subBranches += Form(\" %s\",cTask->GetJetOutputBranch());\n"
	<< "  \n";
      o << "  cTask = AddTaskJetCluster(\"AOD\",\"\",highPt,usePS,\n"
	<<                           "\"ANTIKT\",0.2,0,1,deltaAOD,0.15);\n"
	<< "  cTask->SetCentralityCut(centLow,centHigh);\n"
	<< "  if (pbpb) cTask->SetBackgroundBranch(bgBranch);\n"
	<< "  subBranches += Form(\" %s\",cTask->GetJetOutputBranch());\n"
	<< "  \n";
      o	<< "  if (pbpb) {\n"
	<< "    gROOT->LoadMacro(\"AddTaskJetBackgroundSubtract.C\");\n"
	<< "    AliAnalysisTaskJetBackgroundSubtract* sTask = 0;\n"
	<< "    sTask = AddTaskJetBackgroundSubtract(subBranches,1,\"B0\",\n"
	<< "                                         \"B%d\");\n"
	<< "    sTask->SetBackgroundBranch(bgBranch);\n"
	<< "    sTask->SetNonStdOutputFile(deltaAOD);\n"
	<< "  }\n"
	<< "}\n"
	<< std::endl;
      o.close();

      UInt_t  highPtFilterMask = 272;
      TString deltaAODJetName  = "AliAOD.Jets.root";
      Float_t trackEtaWindow   = 0.9;
      Float_t centLow          = 0;
      Float_t centHigh         = 0;
      gROOT->Macro(Form("./AddJets.C((AliAnalysisManager*)%p,"
			"%d,\"%s\",%d,%d,%d,%f,%f,%f);",
			mgr, highPtFilterMask, deltaAODJetName.Data(),
			jetanDelta, physicsSelection,
			fOptions.AsInt("collision") == 1, // pbpb
			trackEtaWindow, centLow, centHigh));
    }
    if (fOptions.Has("forward-nch") && physicsSelection) {
      Bool_t fwdMC = (mgr->GetMCtruthEventHandler() != 0 && 
		      fOptions.Has("track-refs"));
    
      AddTask("$ALICE_ROOT/PWGLF/FORWARD/analysis2/AddTaskForwardMult.C",
	      Form("%d", fwdMC));
    }
  }
  //------------------------------------------------------------------
  /** 
   * Create physics selection, and add to manager
   * 
   * @param mc Whether this is for MC 
   * @param mgr Manager
   */
  virtual void CreatePhysicsSelection(Bool_t mc, AliAnalysisManager* mgr)
  {
    if (fOptions.Has("physics-selection")) 
      TrainSetup::CreatePhysicsSelection(mc, mgr);
  }
  //------------------------------------------------------------------
  /** 
   * Create centrality selection, and add to manager
   * 
   * @param mc Whether this is for MC 
   * @param mgr Manager
   */
  virtual void CreateCentralitySelection(Bool_t mc, AliAnalysisManager* mgr)
  {
    if (fOptions.Has("centrality")) 
      TrainSetup::CreateCentralitySelection(mc, mgr);
  }
  //------------------------------------------------------------------
  /** 
   * Create MC input handler 
   * 
   * @param mc    Assume monte-carlo input 
   * 
   * @return 
   */
  virtual AliVEventHandler* CreateMCHandler(UShort_t /*type*/, bool mc)
  {
    if (!mc) return 0;
    AliMCEventHandler* mcHandler = new AliMCEventHandler();
    mcHandler->SetReadTR(fOptions.Has("track-refs")); 
    return mcHandler;
  }
  //------------------------------------------------------------------
  /** 
   * Get the name of this set-up 
   * 
   * @return Name of this setup 
   */
  const char* ClassName() const { return "AODTrain"; }
};
//
// EOF
//
