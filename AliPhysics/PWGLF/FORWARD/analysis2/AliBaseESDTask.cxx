#include "AliBaseESDTask.h"
#include "AliFMDEventInspector.h"
#include "AliForwardCorrectionManager.h"
#include "AliForwardUtil.h"
#include "AliFMDCorrELossFit.h"
#include <AliAnalysisManager.h>
#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include <AliAODHandler.h>
#include <AliLog.h>
#include <AliESDEvent.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TFile.h>
#include <iostream>
#include <iomanip>

//____________________________________________________________________
AliBaseESDTask::AliBaseESDTask()
  : AliAnalysisTaskSE(), 
    fFirstEvent(true),
    fList(0),
    fResults(0),
    fNeededCorrections(0),
    fExtraCorrections(0),
    fCloneList(false),
    fCorrManager(0)
{}

//____________________________________________________________________
AliBaseESDTask::AliBaseESDTask(const char* name, const char* title,
			       AliCorrectionManagerBase* manager)
  : AliAnalysisTaskSE(name), 
    fFirstEvent(true),
    fList(0),
    fResults(0),
    fNeededCorrections(0),
    fExtraCorrections(0),
    fCloneList(false),
    fCorrManager(0)
{
  // The line below doesn't actually do the job - when we're
  // constructing, the derived class ins't set yet and this explicitly
  // points to an object of _this_ class.
  // SetTitle(title && title[0] != '\0' ? title : this->ClassName());
  SetTitle(title && title[0] != '\0' ? title : "");
  fCorrManager = manager;
  // if (!manager) 
  //   AliFatal("Must pass in a valid correction manager object!");
  fBranchNames = 
    "ESD:AliESDRun.,AliESDHeader.,AliMultiplicity.,"
    "AliESDFMD.,SPDVertex.,TPCVertex.,PrimaryVertex.";

  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
}
//____________________________________________________________________
Bool_t
AliBaseESDTask::Connect(const char* sumFile, 
			const char* resFile,
			Bool_t      old)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskForwardMult", "No analysis manager to connect to.");
    return false;
  }   

  // Add to the manager 
  mgr->AddTask(this);
  
  // Create and connect output containers 
  TString sumOut;
  TString resOut;
  if      (sumFile && sumFile[0] != '\0') sumOut = sumFile;
  if      (resFile && resFile[0] != '\0') resOut = resFile;
  else if (sumFile && sumFile[0] != '\0') resOut = sumFile;
  // If the string is null or 'default' connect to standard output file 
  if (sumOut.IsNull() || sumOut.EqualTo("default", TString::kIgnoreCase)) 
    sumOut = AliAnalysisManager::GetCommonFileName();
  // If the string is null or 'default' connect to standard output file 
  if (resOut.IsNull() || resOut.EqualTo("default", TString::kIgnoreCase)) 
    resOut = AliAnalysisManager::GetCommonFileName();

  // Always connect input 
  mgr->ConnectInput(this, 0, mgr->GetCommonInputContainer());

  // Connect sum list unless the output 'none' is specified
  if (!sumOut.EqualTo("none", TString::kIgnoreCase)) {
    TString sumName(Form("%s%s", old ? "Forward" : GetName(),
			 old ? "" : "Sums"));
    AliAnalysisDataContainer* sumCon = 
      mgr->CreateContainer(sumName, TList::Class(), 
			   AliAnalysisManager::kOutputContainer, sumOut);
    mgr->ConnectOutput(this, 1, sumCon);
  }
  // Connect the result list unless the output 'none' is specified
  if (!resOut.EqualTo("none", TString::kIgnoreCase)) {
    TString resName(Form("%sResults", GetName()));
    AliAnalysisDataContainer* resCon = 
      mgr->CreateContainer(resName, TList::Class(), 
			   AliAnalysisManager::kParamContainer, resOut);
    mgr->ConnectOutput(this, 2, resCon);
  }
  
  return true;
}

//____________________________________________________________________
TAxis*
AliBaseESDTask::DefaultEtaAxis() const
{
  static TAxis* a = new TAxis(200, -4, 6);
  return a;
}
//____________________________________________________________________
TAxis*
AliBaseESDTask::DefaultVertexAxis() const
{
  static TAxis* a = AliForwardUtil::MakeFullIpZAxis(20);
  return a;
}
//____________________________________________________________________
void
AliBaseESDTask::SetDebug(Int_t dbg)
{
  // 
  // Set debug level 
  // 
  // Parameters:
  //    dbg debug level
  //
  GetEventInspector().SetDebug(dbg);
}
//____________________________________________________________________
void
AliBaseESDTask::SetIPzMethod(const char* str)
{
  AliFMDEventInspector::EVtxType meth = AliFMDEventInspector::kNormal;
  TString s(str);
  s.ToLower();
  if      (s.Contains("normal")) meth = AliFMDEventInspector::kNormal;
  else if (s.Contains("pa2012")) meth = AliFMDEventInspector::kpA2012;
  else if (s.Contains("pa2013")) meth = AliFMDEventInspector::kpA2013;
  else if (s.Contains("ud"))     meth = AliFMDEventInspector::kPWGUD;
  else if (s.Contains("disp") || s.Contains("sat"))
    meth = AliFMDEventInspector::kDisplaced;
  else
    AliWarningF("Unknown IPz method: %s, using normal", str);
  GetEventInspector().SetVertexMethod(meth);
}

//____________________________________________________________________
Bool_t 
AliBaseESDTask::Configure(const char* macro)
{
  // --- Configure the task ------------------------------------------
  TString macroPath(gROOT->GetMacroPath());
  if (!macroPath.Contains("$(ALICE_PHYSICS)/PWGLF/FORWARD/analysis2")) { 
    macroPath.Append(":$(ALICE_PHYSICS)/PWGLF/FORWARD/analysis2");
    gROOT->SetMacroPath(macroPath);
  }
  TString mac(macro);
  if (mac.EqualTo("-default-")) mac = DefaultConfig();
  const char* config = gSystem->Which(gROOT->GetMacroPath(), mac.Data());
  if (!config) {
    AliWarningF("%s not found in %s", mac.Data(), gROOT->GetMacroPath());
    return false;
  }
  if (fTitle.IsNull()) fTitle = this->ClassName();

  AliInfoF("Loading configuration of '%s' from %s",  ClassName(), config);
  gROOT->Macro(Form("%s((%s*)%p)", config, GetTitle(), this));

  AliInfoF("Unloading configuration script %s", config);
  gInterpreter->UnloadFile(config);

  delete config;
 
 return true;
}

//____________________________________________________________________
void 
AliBaseESDTask::LocalInit() 
{ 
  fFirstEvent = true; 
  DGUARD(fDebug,1,"Doing local initialization");
  Setup(); 
}

//____________________________________________________________________
void
AliBaseESDTask::UserCreateOutputObjects()
{
  // 
  // Create output objects 
  // 
  //
  DGUARD(fDebug,1,"Create user ouput");
  fList = new TList;
  fList->SetName(Form("%sSums", GetName()));
  fList->SetOwner();
  
  AliAnalysisManager* am = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler*      ah = 
    dynamic_cast<AliAODHandler*>(am->GetOutputEventHandler());
  //if (!ah) AliFatal("No AOD output handler set in analysis manager");
  if (ah)  CreateBranches(ah);
   
  GetEventInspector().CreateOutputObjects(fList);

  if (!Book()) AliFatalF("Failed to book output objects for %s", GetName());

  // gSystem->Exec("root-config --version --prefix");
  PostData(1, fList);
}

//____________________________________________________________________
void
AliBaseESDTask::UserExec(Option_t*)
{
  // Call pre-event setup 
  PreEvent();

  // Read in selected branches 
  LoadBranches();

  // Get the input data 
  AliESDEvent* esd = GetESDEvent();
  if (!esd) return;

  // Call the user code with our event passed in 
  Event(*esd);
  // if (!Event(*esd)) {
  //   AliWarningF("Failed to process the event for %s", GetName());
  //   return;
  // }

  // Post data 
  PostData(1, fList);

  // Call post-event processing 
  PostEvent();
}

//____________________________________________________________________
Bool_t
AliBaseESDTask::StoreTrainName(Int_t no)
{
  AliAnalysisDataSlot* slot = GetOutputSlot(no);
  if (!slot) return false;

  AliAnalysisDataContainer* cont = slot->GetContainer();
  if (!cont) return false;

  TFile* file = cont->GetFile();
  if (!file || !file->IsWritable()) return false;

  TDirectory* save  = gDirectory;
  TNamed* tag = new TNamed("trainName",
			   AliAnalysisManager::GetAnalysisManager()->GetName());
  file->cd();
  tag->Write();
  save->cd();

  return true;
}

//____________________________________________________________________
void
AliBaseESDTask::Terminate(Option_t*)
{
  TList* list = dynamic_cast<TList*>(GetOutputData(1));
  if (!list) {
    AliError(Form("No output list defined (%p)", GetOutputData(1)));
    if (GetOutputData(1)) GetOutputData(1)->Print();
    return;
  }

  // Assign to our internal variable for use by sub-classes 
  fList = list;

  // Create our output container 
  TString resName(Form("%sResults", GetName()));
  if (fCloneList)
    fResults = static_cast<TList*>(fList->Clone(resName));
  else {
    fResults = new TList;
    fResults->SetName(resName);
  }
  fResults->SetOwner();

  // Now call user defined routines 
  if (!Finalize()) {
    AliErrorF("Failed to finalize this task (%s)", GetName());
    return;
  }

  // Store name in output
  if (!StoreTrainName(2)) 
    fResults->Add(new TNamed("trainName",
			     AliAnalysisManager::GetAnalysisManager()
			     ->GetName()));
  PostData(2, fResults);
}

//____________________________________________________________________
Bool_t 
AliBaseESDTask::PreData(const TAxis&, const TAxis&) 
{ 
  return true; 
}


//____________________________________________________________________
Bool_t 
AliBaseESDTask::CheckCorrections(UInt_t what) const
{
  // 
  // Check if all needed corrections are there and accounted for.  If not,
  // do a Fatal exit 
  // 
  // Parameters:
  //    what Which corrections is needed
  // 
  // Return:
  //    true if all present, false otherwise
  //  
  DGUARD(fDebug,1,"Checking corrections 0x%x", what);
  if (what == 0) return true;

  AliCorrectionManagerBase* cm = GetManager();
  if (!cm) {
    AliErrorF("Check corrections=0x%x not null, "
	      "but no correction manager defined!", 
	      what);
    return false;
  }
  Bool_t ret = cm->CheckCorrections(what);
  return ret;
}
//____________________________________________________________________
Bool_t
AliBaseESDTask::ReadCorrections(const TAxis*& pe, 
				const TAxis*& pv, 
				Bool_t        mc,
				Bool_t        sat)
{
  //
  // Read corrections
  //
  //
  UInt_t what = fNeededCorrections|fExtraCorrections;
  
  DGUARD(fDebug,1,"Read corrections 0x%x", what);

  AliCorrectionManagerBase* cm = GetManager();
  if (!cm && fNeededCorrections) {
    AliErrorF("Needed/extra corrections=0x%x/0x%x not null, "
	      "but no correction manager defined!", 
	      fNeededCorrections, fExtraCorrections);
    return false;
  }
  if (!cm || !what) {
    // In case we have no needed corrections, we can return here 
    if (!pe) pe = DefaultEtaAxis();
    if (!pv) pv = DefaultVertexAxis();
    return true;
  }
  cm->EnableCorrections(what);
  if (!cm->InitCorrections(GetEventInspector().GetRunNumber(),
			   GetEventInspector().GetCollisionSystem(),
			   GetEventInspector().GetEnergy(),
			   GetEventInspector().GetField(),
			   mc,
			   sat,
			   false)) { 
    AliWarning("Failed to read in some corrections, making task zombie");
    return false;
  }
  if (!CheckCorrections(fNeededCorrections)) return false;

  // Sett our persistency pointer 
  // fCorrManager = &fcm;

  // Get the eta axis from the secondary maps - if read in
  if (!pe) {
    pe = cm->GetEtaAxis();
    if (!pe) pe = DefaultEtaAxis();
  }
  // Get the vertex axis from the secondary maps - if read in
  if (!pv) {
    pv = cm->GetVertexAxis();
    if (!pv) pv = DefaultVertexAxis();
  }

  return true;
}
//____________________________________________________________________
AliESDEvent*
AliBaseESDTask::GetESDEvent()
{
  //
  // Get the ESD event. IF this is the first event, initialise
  //
  DGUARD(fDebug,1,"Get the ESD event");

  // If we're marked as a zombie, do nothing and return a null
  if (IsZombie()) return 0;

  // Try to get the ESD event 
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!esd) {
    AliWarning("No ESD event found for input event");
    return 0;
  }

  // --- Load the data -----------------------------------------------
  LoadBranches();

  if (!fFirstEvent || !esd->GetESDRun()) return esd;

  // On the first event, initialize the parameters
  GetEventInspector().SetMC(MCEvent());
  GetEventInspector().ReadRunDetails(esd);
  
  AliInfoF("Initializing with parameters from the ESD:\n"
	   "         AliESDEvent::GetBeamEnergy()   ->%f\n"
	   "         AliESDEvent::GetBeamType()     ->%s\n"
	   "         AliESDEvent::GetCurrentL3()    ->%f\n"
	   "         AliESDEvent::GetMagneticField()->%f\n"
	   "         AliESDEvent::GetRunNumber()    ->%d",
	   esd->GetBeamEnergy(),
	   esd->GetBeamType(),
	   esd->GetCurrentL3(),
	   esd->GetMagneticField(),
	   esd->GetRunNumber());

  PreCorrections(esd);

  fFirstEvent = false;
  
  const   TAxis* pe = 0;
  const   TAxis* pv = 0;
  Bool_t  mc        = IsMC();
  Bool_t  sat       = false;
  Bool_t  ret       = ReadCorrections(pe, pv, mc, sat);
  if (!ret) {
    AliError("Failed to read corrections, making this a zombie");
    SetZombie(true);
    return 0;
  }
  Printf("Vertex axis: %p   Eta axis: %p", pv, pe);
  if (!pv) AliFatal("No vertex axis defined");
  if (!pe) AliFatal("No eta axis defined");

  // Initialize the event inspector 
  GetEventInspector().SetupForData(*pv);
  
  // Initialize the remaining stuff 
  if (!PreData(*pv, *pe)) {
    AliError("Failed to initialize sub-algorithms, making this a zombie");
    SetZombie(true);
    return 0;
  }
  
  this->Print("R");

  return esd;
}

//____________________________________________________________________
void
AliBaseESDTask::MarkEventForStore() const
{
  // Make sure the AOD tree is filled 
  DGUARD(fDebug,3,"Mark AOD event for storage");
  AliAnalysisManager* am = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler*      ah = 
    dynamic_cast<AliAODHandler*>(am->GetOutputEventHandler());
  if (ah) ah->SetFillAOD(kTRUE);
}
#define PF(N,V,...)					\
  AliForwardUtil::PrintField(N,V, ## __VA_ARGS__)
#define PFB(N,FLAG)				\
  do {									\
    AliForwardUtil::PrintName(N);					\
    std::cout << std::boolalpha << (FLAG) << std::noboolalpha << std::endl; \
  } while(false)
#define PFV(N,VALUE)					\
  do {							\
    AliForwardUtil::PrintName(N);			\
    std::cout << (VALUE) << std::endl; } while(false)

//____________________________________________________________________
void
AliBaseESDTask::Print(Option_t* option) const
{
  // 
  // Print information 
  // 
  // Parameters:
  //    option Not used
  //  
  std::cout << std::setfill('=') << std::setw(75) << "=" 
	    << std::setfill(' ') << std::endl;
  AliForwardUtil::PrintTask(*this);
  gROOT->IncreaseDirLevel();
  PF("Off-line trigger mask", "0x%0x", fOfflineTriggerMask);
  if (GetManager()) GetManager()->Print(option);
  else  PF("No correction manager","");

  GetEventInspector().Print(option);
  gROOT->DecreaseDirLevel();
}
//
// EOF
//
