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
    fOptions.Add("centrality", "Use centrality", false);
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
  //==================================================================
  void CouplePidResponse()
  {
    CoupleCar("AddTaskPIDResponse.C");
  }
  //------------------------------------------------------------------
  void CouplePidQA()
  {
    AliAnalysisTaskSE* seTask = CoupleSECar("AddTaskPIDqa.C");
    seTask->SelectCollisionCandidates(AliVEvent::kAny);
  }
  //------------------------------------------------------------------
  void CoupleHighPt()
  {
    CoupleCar("AddTaskFilteredTree.C");
  }
  //------------------------------------------------------------------
  void CoupleFilter()
  {
    fRailway->LoadLibrary("PWGmuon");
    if (fOptions.Has("muon-copy")) {
      mgr->RegisterExtraFile("AliAOD.Muons.root");
      mgr->RegisterExtraFile("AliAOD.Dimuons.root");
    }
    Int_t runFlag = (fOptions.AsInt("year") % 100) * 100;
    TString args = 
      TString::Format("%d,%d,false,false,false,false,true,false,false,%d",
		      fOptions.Has("kinematics-filter"), 
		      fOptions.Has("muon-copy"), runFlag);
    CoupleCar("AddTaskESDFilter.C", args);
  }
  //------------------------------------------------------------------
  void CoupleVertexing()
  {
    fRailway->LoadLibrary("PWGflowBase");
    fRailway->LoadLibrary("PWGflowTasks");
    fRailway->LoadLibrary("PWGHFvertexingHF");
    AliAnalysisTaskSE* seTask = 
      CoupleSECar("$ALICE_PHYSICS/PWGHF/vertexingHF/macros/AddTaskVertexingHF.C");
    seTask->SelectCollisionCandidates(0);
    mgr->RegisterExtraFile("AliAOD.VertexingHF.root");
  }
  //------------------------------------------------------------------
  void CoupleJPsiFilter()
  {
    fRailway->LoadLibrary("PWGDQdielectron");
    AliAnalysisTaskSE* seTask = 
      CoupleSECar("$ALICE_PHYSICS/PWGDQ/dielectron/macros/AddTaskJPSIFilter.C");
    seTask->SelectCollisionCandidates(0);
    mgr->RegisterExtraFile("AliAOD.Dielectron.root");
  }
  //------------------------------------------------------------------
  void CoupleD0Decay()
  {
    CoupleCar("$ALICE_PHYSICS/PWGHF/vertexingHF/AddD2HTrain.C",
	      "kFALSE, 1,0,0,0,0,0,0,0,0,0,0");
  }
  //------------------------------------------------------------------
  void CoupleJetAN()
  {
    fRailway->LoadLibrary("JETAN");
      
    Bool_t jetanDelta = fOptions.Has("jetan-delta");      
    if (jetanDelta) {
      fRailway->LoadLibrary("CGAL");
      fRailway->LoadLibrary("fastjet");
      fRailway->LoadLibrary("siscone");
      fRailway->LoadLibrary("SISConePlugin");
      fRailway->LoadLibrary("FASTJETAN");
      // --- For newer fastjet ---
      // fRailway->LoadLibrary("CGAL");
      // fRailway->LoadLibrary("fastjet");
      // fRailway->LoadLibrary("fastjettools");
      // fRailway->LoadLibrary("siscone");
      // fRailway->LoadLibrary("siscone_spherical");
      // fRailway->LoadLibrary("fastjetplugins");
      // fRailway->LoadLibrary("FASTJETAN");
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
  //------------------------------------------------------------------
  void CoupleForward()
  {
    Bool_t fwdMC = (mgr->GetMCtruthEventHandler() != 0 && 
		    fOptions.Has("track-refs"));
    
    CoupleCar("$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/AddTaskForwardMult.C",
	      Form("%d", fwdMC));
  }    
  //------------------------------------------------------------------
  /** 
   * Create our tasks 
   * 
   * @param mgr Manager 
   */
  void CreateTasks(AliAnalysisManager* mgr)
  {

    fRailway->LoadLibrary("libCORRFW");
    Bool_t physicsSelection = fOptions.Has("physics-selection");

    if (fOptions.Has("sys-info")) mgr->SetNSysInfo(100);

    const char* aliPhys = "$(ALICE_PHYSICS)";
    // AliAnalysisTask*   task = 0;
    AliAnalysisTaskSE* seTask = 0;

    AliAnalysisManager::SetCommonFileName("AODQA.root");
    if (fOptions.Has("pid-response")) CouplePidResponse();
    if (fOptions.Has("pid-qa"))       CouplePidQA();    
    if (fOptions.Has("high-pt"))      CoupleHighPt();
    if (fOptions.Has("esd-filter"))   CoupleFilter();

    Bool_t vertexing = fOptions.Has("vertexing");
    Bool_t d0Decay   = fOptions.Has("d0-decay");
    if (vertexing || d0Decay) {
      TString src = "$ALICE_PHYSICS/PWGHF/vertexingHF/ConfigVertexingHF.C";
      if (fOptions.AsInt("collision") == 1) src.ReplaceAll(".C", "highmult.C");
      TFile::Cp(src, "ConfigVertexingHF.C");
      fRailway->LoadAux("ConfigVertexingHF.C");
    }    

    if (vertexing)                      CoupleVertexing();
    if (fOptions.Has("jpsi-filtering")) CoupleJPsiFilter();
    if (d0Decay)                        CoupleD0Decay();
    if (fOptions.Has("jetan"))          CoupleJetAN();
    if (fOptions.Has("forward-nch") && physicsSelection) CoupleForward();
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
  virtual void CreateCentralitySelection(Bool_t mc)
  {
    if (fOptions.Has("centrality")) 
      TrainSetup::CreateCentralitySelection(mc);
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
