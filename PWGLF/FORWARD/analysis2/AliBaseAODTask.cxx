#include "AliBaseAODTask.h"
#include "AliForwardUtil.h"
#include "AliAODForwardMult.h"
#include "AliAODCentralMult.h"
#include <AliAnalysisManager.h>
#include <AliLog.h>
#include <AliAODEvent.h>
#include <TROOT.h>

//____________________________________________________________________
AliBaseAODTask::AliBaseAODTask()
  : AliAnalysisTaskSE(),
    fTriggerMask(0xFFFFFFFF), 
    fMinIpZ(0), 
    fMaxIpZ(-1), 
    fCentAxis(0, 0, -1), 
    fTriggers(0), 
    fEventStatus(0), 
    fVertex(0),
    fCent(0),
    fAccVertex(0),
    fAccCent(0),
    fFirstEvent(true),
    fCloneList(false),
    fSums(0), 
    fResults(0)
{
}
//____________________________________________________________________
AliBaseAODTask::AliBaseAODTask(const char* name)
  : AliAnalysisTaskSE(name),
    fTriggerMask(0xFFFFFFFF), 
    fMinIpZ(0), 
    fMaxIpZ(-1), 
    fCentAxis(0, 0, -1), 
    fTriggers(0), 
    fEventStatus(0), 
    fVertex(0),
    fCent(0),
    fAccVertex(0),
    fAccCent(0),
    fFirstEvent(true), 
    fCloneList(false),
    fSums(0), 
    fResults(0)
{
  fCentAxis.SetName("centAxis");
  fCentAxis.SetTitle("Centrality [%]");
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
}
//________________________________________________________________________
void 
AliBaseAODTask::SetTriggerMask(const char* mask)
{
  // 
  // Set the trigger maskl 
  // 
  // Parameters:
  //    mask Trigger mask
  //
  DGUARD(fDebug,3,"Set the trigger mask: %s", mask);
  SetTriggerMask(AliAODForwardMult::MakeTriggerMask(mask));
}
//________________________________________________________________________
void 
AliBaseAODTask::SetTriggerMask(UShort_t mask) 
{ 
  DGUARD(fDebug,3,"Set the trigger mask: 0x%0x", mask);
  fTriggerMask = mask; 
  // if (fTriggerString) delete fTriggerString;
  // fTriggerString = AliForwardUtil::MakeParameter("trigger", fTriggerMask);
}
//________________________________________________________________________
void 
AliBaseAODTask::SetCentralityAxis(UShort_t n, Short_t* bins)
{
  DGUARD(fDebug,3,"Set centrality axis, %d bins", n);
  TArrayD dbins(n+1);
  for (UShort_t i = 0; i <= n; i++) 
    dbins[i] = (bins[i] == 100 ? 100.1 : bins[i]);
  fCentAxis.Set(n, dbins.GetArray());
}
//________________________________________________________________________
void 
AliBaseAODTask::SetCentralityAxis(Short_t low, Short_t high)
{
  Short_t a[] = { low, high };
  SetCentralityAxis(1, a);
}

//____________________________________________________________________
Bool_t
AliBaseAODTask::Connect(const char* sumFile, 
			const char* resFile)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskForwardMult", "No analysis manager to connect to.");
    return false;
  }   

  // --- Check that we have an AOD input handler ---------------------
  UShort_t aodInput = 0;
  if (!(aodInput = AliForwardUtil::CheckForAOD())) {
    AliError("Cannot proceed without and AOD handler");
    return false;
  }
  if (aodInput == 2 &&
      !AliForwardUtil::CheckForTask("AliForwardMultiplicityBase")) {
    AliError("The relevant task wasn't added to the train");
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
  if (sumOut.IsNull()) sumOut = AliAnalysisManager::GetCommonFileName();
  if (resOut.IsNull()) resOut = AliAnalysisManager::GetCommonFileName();

  AliAnalysisDataContainer* sumCon = 
    mgr->CreateContainer(Form("%sSums", GetName()), TList::Class(), 
			 AliAnalysisManager::kOutputContainer, sumOut);
  AliAnalysisDataContainer* resCon = 
    mgr->CreateContainer(Form("%sResults", GetName()), TList::Class(), 
			 AliAnalysisManager::kParamContainer, resOut);
  mgr->ConnectInput(this, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(this, 1, sumCon);
  mgr->ConnectOutput(this, 2, resCon);
  
  return true;
}
//____________________________________________________________________
void
AliBaseAODTask::UserCreateOutputObjects()
{
  // 
  // Create output objects 
  // 
  //
  DGUARD(fDebug,1,"Create user ouput");
  fSums = new TList;
  fSums->SetName(Form("%sSums", GetName()));
  fSums->SetOwner();

  fTriggers = AliAODForwardMult::MakeTriggerHistogram("triggers",fTriggerMask);
  fTriggers->SetDirectory(0);

  fEventStatus = AliAODForwardMult::MakeStatusHistogram("status");
  fEventStatus->SetDirectory(0);

  fSums->Add(fTriggers);
  fSums->Add(fEventStatus);

  TAxis* vA = AliForwardUtil::MakeFullIpZAxis(20);
  fVertex = new TH1D("vertex", "IP_{z} of all events", 
		     vA->GetNbins(), vA->GetXbins()->GetArray());
  fVertex->SetXTitle("IP_{z} [cm]");
  fVertex->SetYTitle("Events");
  fVertex->SetDirectory(0);
  fVertex->SetFillColor(kRed+2);
  fVertex->SetFillStyle(3001);
  fVertex->SetLineColor(kRed+2);
  fSums->Add(fVertex);
  fAccVertex = static_cast<TH1*>(fVertex->Clone("vertexAcc"));
  fAccVertex->SetTitle("IP_{z} of accepted events");
  fAccVertex->SetDirectory(0);
  fAccVertex->SetFillColor(kGreen+2);
  fAccVertex->SetLineColor(kGreen+2);
  fSums->Add(fAccVertex);

  fCent = new TH1D("cent","Centrality of all events",100, 0, 100);
  fCent->SetXTitle("Centrality [%]");
  fCent->SetYTitle("Events");
  fCent->SetFillColor(kRed+2);
  fCent->SetFillStyle(3001);
  fCent->SetLineColor(kRed+2);
  fCent->SetDirectory(0);
  fSums->Add(fCent);
  fAccCent = static_cast<TH1*>(fCent->Clone("centAcc"));
  fAccCent->SetTitle("Centrality of accepted events");
  fAccCent->SetDirectory(0);
  fAccCent->SetFillColor(kGreen+2);
  fAccCent->SetLineColor(kGreen+2);
  fSums->Add(fAccCent);


  if (!Book()) AliFatalF("Failed to book output objects for %s", GetName());

  // Store centrality axis as a histogram - which can be merged
  fSums->Add(fCentAxis.Clone("centAxis"));
  fSums->Add(AliForwardUtil::MakeParameter("trigger", ULong_t(fTriggerMask)));
  fSums->Add(AliForwardUtil::MakeParameter("count", 1));
  fSums->Add(AliForwardUtil::MakeParameter("alirootRev", 
					   AliForwardUtil::AliROOTRevision()));
  fSums->Add(AliForwardUtil::MakeParameter("alirootBranch", 
					   AliForwardUtil::AliROOTBranch()));

  Print();

  PostData(1, fSums);
}

//____________________________________________________________________
AliAODForwardMult*
AliBaseAODTask::GetForward(const AliAODEvent& aod, Bool_t mc, Bool_t verb)
{
  // Get the forward object that contains our event selection stuff 
  TObject* obj = 0;
  if (mc) obj = aod.FindListObject("ForwardMC");
  else    obj = aod.FindListObject("Forward");
  if (!obj) { 
    if (verb) AliWarning("No forward object found");
    return 0;
  }
  AliAODForwardMult* forward = static_cast<AliAODForwardMult*>(obj);
  return forward;
}
//____________________________________________________________________
AliAODCentralMult*
AliBaseAODTask::GetCentral(const AliAODEvent& aod, Bool_t mc, Bool_t verb)
{
  // Get the forward object that contains our event selection stuff 
  TObject* obj = 0;
  if (mc) obj = aod.FindListObject("CentralClustersMC");
  else    obj = aod.FindListObject("CentralClusters");
  if (!obj) { 
    if (verb) AliWarning("No central object found");
    return 0;
  }
  AliAODCentralMult* central = static_cast<AliAODCentralMult*>(obj);
  return central;
}
//____________________________________________________________________
TH2D*
AliBaseAODTask::GetPrimary(const AliAODEvent& aod)
{
  TObject* obj = aod.FindListObject("primary");
  // We should have a forward object at least 
  if (!obj) {
    return 0;
  }
  TH2D* ret = static_cast<TH2D*>(obj);
  return ret;
}
  
//____________________________________________________________________
void 
AliBaseAODTask::UserExec(Option_t *) 
{
  // 
  // Process a single event 
  // 
  // Parameters:
  //    option Not used
  //
  // Main loop
  DGUARD(fDebug,1,"Analyse the AOD event");
  if (!PreEvent()) return;

  AliAODEvent* aod = AliForwardUtil::GetAODEvent(this);
  if (!aod) return;

  // Get the forward object that contains our event selection stuff 
  AliAODForwardMult* forward = GetForward(*aod);
  if (!forward) return;

  if (fFirstEvent) { 
    if (!PreData()) return;
    StoreInformation(*forward);
    fFirstEvent = false;
  }

  // Get our ip_z and centrality 
  Double_t vtx   = forward->GetIpZ();
  Float_t  cent  = forward->GetCentrality();
  fVertex->Fill(vtx);
  fCent->Fill(cent);

  // Now check our event selectio up front 
  if (!CheckEvent(*forward)) return;

  // Let user defined code do the job 
  Bool_t taken = Event(*aod);

  // Fill our histograms 
  if (taken) {
    fAccVertex->Fill(vtx);
    fAccCent->Fill(cent);
  }

  PostData(1, fSums);
  
  PostEvent();
}

//____________________________________________________________________
Bool_t
AliBaseAODTask::CheckEvent(const AliAODForwardMult& forward)
{
  if (HasCentrality())
    return forward.CheckEvent(fTriggerMask, fMinIpZ, fMaxIpZ, 
			      UShort_t(fCentAxis.GetXmin()), 
			      UShort_t(fCentAxis.GetXmax()), 
			      fTriggers, fEventStatus);
 return forward.CheckEvent(fTriggerMask, fMinIpZ, fMaxIpZ, 
			   0, 0, fTriggers, fEventStatus);
}

//____________________________________________________________________
void
AliBaseAODTask::StoreInformation(AliAODForwardMult& forward)
{
  fSums->Add(AliForwardUtil::MakeParameter("sNN", forward.GetSNN()));
  fSums->Add(AliForwardUtil::MakeParameter("sys", forward.GetSystem()));
}

//____________________________________________________________________
void
AliBaseAODTask::Terminate(Option_t*)
{
  TList* list = dynamic_cast<TList*>(GetOutputData(1));
  if (!list) {
    AliError(Form("No output list defined (%p)", GetOutputData(1)));
    if (GetOutputData(1)) GetOutputData(1)->Print();
    return;
  }

  // Assign to our internal variable for use by sub-classes 
  fSums = list;

  // Create our output container 
  TString resName(Form("%sResults", GetName()));
  if (fCloneList)
    fResults = static_cast<TList*>(fSums->Clone(resName));
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

  PostData(2, fResults);
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
AliBaseAODTask::Print(Option_t* /*option=""*/) const 
{
  /** 
   * Print this task 
   * 
   * @param option Not used
   */
  AliForwardUtil::PrintTask(*this);
  gROOT->IncreaseDirLevel();
  PFV("Trigger mask",  AliAODForwardMult::GetTriggerString(fTriggerMask));
  PF("IP z range", "%++6.1f - %+6.1f", fMinIpZ, fMaxIpZ);
  PFV("Centrality bins", (HasCentrality() ? "" : "none"));
  if (HasCentrality()) {
    Int_t           nBins = fCentAxis.GetNbins();
    const Double_t* bins  = fCentAxis.GetXbins()->GetArray();
    for (Int_t i = 0; i <= nBins; i++) 
      std::cout << (i==0 ? "" : "-") << bins[i];
    std::cout << std::endl;
  }
  gROOT->DecreaseDirLevel();
}
//
// EOF
//
