#include <string>

#include "TROOT.h"
#include "TString.h"
#include "AliAODHandler.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskValidation.h"
#include "AliOADBPhysicsSelection.h"
#include "AliPhysicsSelectionTask.h"
//#include "AddTaskForwardFlowRun2.h"
#include "AliVEvent.h"
#include "AliAnalysisDataSlot.h"
#include "AliForwardFlowRun2Task.h"
#include <sstream>
#include "AliMultSelectionTask.h"

#include "AliAnalysisTaskValidation.h"
enum mode {kRECON, kTRUTH};

void ConfigureTrain() {

  // Add mult selection Task
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
  gROOT->ProcessLine("AddTaskMultSelection()");

  // PhysicsSelection Configuration
  //gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  //AliPhysicsSelectionTask* ps = reinterpret_cast<AliPhysicsSelectionTask*>
  // Signature: Bool_t mCAnalysisFlag, Bool_t applyPileupCuts
  //gROOT->ProcessLine("AddTaskPhysicsSelection(false, true)");


  AliAnalysisTaskValidation* validation_task = AliAnalysisTaskValidation::ConnectTask("", false) ;

  gROOT->LoadMacro("$ALICE_PHYSICS/PWGCF/FLOW/macros/AddTaskForwardFlowRun2.C");

  bool doNUA, nua_mode, doetagap, mc, esd, prim_cen, prim_fwd;
  double gap;
  int tracktype;
  std::string centrality, name, nua_file;

  // ___________________________________________________________________________
  doNUA = true;
  nua_file = "/home/thoresen/Programs/alice/AliPhysics/PWGCF/FLOW/Forward/corrections/";
  nua_mode = false;
  doetagap = true;
  gap = 0.0;
  mc = false;
  esd = false;
  prim_cen = false;
  prim_fwd = false;
  tracktype = 768;
  centrality = "\"V0M\"";
  name = "\"Hybrid\"";

  std::stringstream inputstring;
  inputstring << "AddTaskForwardFlowRun2("
  << std::boolalpha << doNUA << ", "
  << nua_file << ", "
  << std::boolalpha << nua_mode << ", "
  << std::boolalpha << doetagap << ", "
  << gap << ", "
  << std::boolalpha << mc << ", "
  << std::boolalpha << esd << ", "
  << std::boolalpha << prim_cen << ", "
  << std::boolalpha << prim_fwd << ", "
  << tracktype << ", "
  << centrality << ", "
  << name << ")";

  AliForwardFlowRun2Task* mytask1 =
    reinterpret_cast<AliForwardFlowRun2Task*>
    (gROOT->ProcessLine(inputstring.str().c_str()));

  mytask1->ConnectInput(1,validation_task->GetOutputSlot(2)->GetContainer());

  // ___________________________________________________________________________
  tracktype = 128;
  nua_file = "/home/thoresen/Programs/alice/AliPhysics/PWGCF/FLOW/Forward/corrections/";
  name = "\"TPConly\"";

  std::stringstream inputstring2;
  inputstring2 << "AddTaskForwardFlowRun2("
  << std::boolalpha << doNUA << ", "
  << nua_file << ", "
  << std::boolalpha << nua_mode << ", "
  << std::boolalpha << doetagap << ", "
  << gap << ", "
  << std::boolalpha << mc << ", "
  << std::boolalpha << esd << ", "
  << std::boolalpha << prim_cen << ", "
  << std::boolalpha << prim_fwd << ", "
  << tracktype << ", "
  << centrality << ", "
  << name << ")";

  //mc prim esd
  AliForwardFlowRun2Task* mytask2 =
    reinterpret_cast<AliForwardFlowRun2Task*>
    (gROOT->ProcessLine(inputstring2.str().c_str()));
  mytask2->ConnectInput(1,validation_task->GetOutputSlot(2)->GetContainer());


  // ___________________________________________________________________________
  tracktype = 32;
  nua_file = "/home/thoresen/Programs/alice/AliPhysics/PWGCF/FLOW/Forward/corrections/";
  name = "\"Global\"";

  std::stringstream inputstring3;
  inputstring2 << "AddTaskForwardFlowRun2("
  << std::boolalpha << doNUA << ", "
  << nua_file << ", "
  << std::boolalpha << nua_mode << ", "
  << std::boolalpha << doetagap << ", "
  << gap << ", "
  << std::boolalpha << mc << ", "
  << std::boolalpha << esd << ", "
  << std::boolalpha << prim_cen << ", "
  << std::boolalpha << prim_fwd << ", "
  << tracktype << ", "
  << centrality << ", "
  << name << ")";

  //mc prim esd
  AliForwardFlowRun2Task* mytask3 =
    reinterpret_cast<AliForwardFlowRun2Task*>
    (gROOT->ProcessLine(inputstring3.str().c_str()));
  mytask3->ConnectInput(1,validation_task->GetOutputSlot(2)->GetContainer());

// ___________________________________________________________________________
  tracktype = 64;
  nua_file = "/home/thoresen/Programs/alice/AliPhysics/PWGCF/FLOW/Forward/corrections/";
  name = "\"GlobalLoose\"";

  std::stringstream inputstring4;
  inputstring4 << "AddTaskForwardFlowRun2("
  << std::boolalpha << doNUA << ", "
  << nua_file << ", "
  << std::boolalpha << nua_mode << ", "
  << std::boolalpha << doetagap << ", "
  << gap << ", "
  << std::boolalpha << mc << ", "
  << std::boolalpha << esd << ", "
  << std::boolalpha << prim_cen << ", "
  << std::boolalpha << prim_fwd << ", "
  << tracktype << ", "
  << centrality << ", "
  << name << ")";

  //mc prim esd
  AliForwardFlowRun2Task* mytask4 =
    reinterpret_cast<AliForwardFlowRun2Task*>
    (gROOT->ProcessLine(inputstring4.str().c_str()));
  mytask4->ConnectInput(1,validation_task->GetOutputSlot(2)->GetContainer());

  // ___________________________________________________________________________
    tracktype = 96;
    nua_file = "/home/thoresen/Programs/alice/AliPhysics/PWGCF/FLOW/Forward/corrections/";
    name = "\"GlobalCombined\"";

    std::stringstream inputstring5;
    inputstring5 << "AddTaskForwardFlowRun2("
    << std::boolalpha << doNUA << ", "
    << nua_file << ", "
    << std::boolalpha << nua_mode << ", "
    << std::boolalpha << doetagap << ", "
    << gap << ", "
    << std::boolalpha << mc << ", "
    << std::boolalpha << esd << ", "
    << std::boolalpha << prim_cen << ", "
    << std::boolalpha << prim_fwd << ", "
    << tracktype << ", "
    << centrality << ", "
    << name << ")";

    //mc prim esd
    AliForwardFlowRun2Task* mytask5 =
      reinterpret_cast<AliForwardFlowRun2Task*>
      (gROOT->ProcessLine(inputstring5.str().c_str()));
    mytask5->ConnectInput(1,validation_task->GetOutputSlot(2)->GetContainer());
}
