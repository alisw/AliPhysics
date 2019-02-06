#include "AliBaseAODTask.h"
#include "AliForwardUtil.h"
#include "AliAODForwardMult.h"
#include "AliAODCentralMult.h"
#include "AliAODMultEventClass.h"
#include <AliAnalysisManager.h>
#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include <AliLog.h>
#include <AliAODEvent.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TFile.h>
#include <TH2.h>
#include <iostream>

//____________________________________________________________________
AliBaseAODTask::AliBaseAODTask()
  : AliAnalysisTaskSE(),
    fTriggerMask(0xFFFFFFFF),
    fFilterMask(AliAODForwardMult::kDefaultFilter),
    fCentAxis(0, 0, -1),
    fAbsMinCent(-1),
    fIPzAxis(1,-10,10),
    fTriggers(0),
    fEventStatus(0),
    fVertex(0),
    fCent(0),
    fAccVertex(0),
    fAccVertexXY(0),
    fAccCent(0),
    fFirstEvent(true),
    fCloneList(false),
    fSums(0),
    fResults(0)
{
}
//____________________________________________________________________
AliBaseAODTask::AliBaseAODTask(const char* name,
			       const char* title)
  : AliAnalysisTaskSE(name),
    fTriggerMask(0xFFFFFFFF),
    fFilterMask(AliAODForwardMult::kDefaultFilter),
    fCentAxis(0, 0, -1),
    fAbsMinCent(-1),
    fIPzAxis(1,-10,10),
    fTriggers(0),
    fEventStatus(0),
    fVertex(0),
    fCent(0),
    fAccVertex(0),
    fAccVertexXY(0),
    fAccCent(0),
    fFirstEvent(true),
    fCloneList(false),
    fSums(0),
    fResults(0)
{
  SetTitle(title && title[0] != '\0' ? title : this->ClassName());
  fCentAxis.SetName("centAxis");
  fCentAxis.SetTitle("Centrality [%]");
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
}

//____________________________________________________________________
Bool_t
AliBaseAODTask::Configure(const char* macro)
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
  // if (gInterpreter->IsLoaded(config))
  // gInterpreter->UnloadFile(config);

  AliInfoF("Loading configuration of '%s' from %s",  ClassName(), config);
  gROOT->Macro(Form("%s((%s*)%p)", config, GetTitle(), this));

  Info("Configure", "Unloading configuration script");
  gInterpreter->UnloadFile(config);
  delete config;

 return true;
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
AliBaseAODTask::SetTriggerMask(UInt_t mask)
{
  DGUARD(fDebug,3,"Set the trigger mask: 0x%0x", mask);
  fTriggerMask = mask;
}
//________________________________________________________________________
void
AliBaseAODTask::SetFilterMask(const char* mask)
{
  //
  // Set the trigger maskl
  //
  // Parameters:
  //    mask Trigger mask
  //
  DGUARD(fDebug,3,"Set the filter mask: %s", mask);
  SetFilterMask(AliAODForwardMult::MakeTriggerMask(mask, "|"));
}
//________________________________________________________________________
void
AliBaseAODTask::SetFilterMask(UInt_t mask)
{
  DGUARD(fDebug,3,"Set the filter mask: 0x%0x", mask);
  fFilterMask = mask;
}
//====================================================================
void AliBaseAODTask::FixAxis(TAxis& axis, const char* title)
{
  if (title && title[0] != '\0') axis.SetTitle(title);
  axis. SetNdivisions(210);
  axis. SetLabelFont(42);
  axis. SetLabelSize(0.03);
  axis. SetLabelOffset(0.005);
  axis. SetLabelColor(kBlack);
  axis. SetTitleOffset(1);
  axis. SetTitleFont(42);
  axis. SetTitleSize(0.04);
  axis. SetTitleColor(kBlack);
  axis. SetTickLength(0.03);
  axis. SetAxisColor(kBlack);
}
//____________________________________________________________________
void AliBaseAODTask::SetAxis(TAxis& axis, Int_t n, Double_t* borders)
{
  axis.Set(n, borders);
  FixAxis(axis);
}
//____________________________________________________________________
void AliBaseAODTask::SetAxis(TAxis&         axis,
			     const TString& spec,
			     const char*    sep)
{
  TString s(spec);
  Bool_t isRange = false, isUnit = false;
  if (s[0] == 'r' || s[0] == 'R') {
    isRange = true;
    s.Remove(0,1);
  }
  if (s[0] == 'u' || s[0] == 'U') {
    isUnit = true;
    s.Remove(0,1);
  }
  TObjArray*  tokens = s.Tokenize(sep);
  TArrayD     bins(tokens->GetEntries());
  TObjString* token = 0;
  TIter       next(tokens);
  Int_t       i = 0;
  while ((token = static_cast<TObjString*>(next()))) {
    Double_t v = token->String().Atof();
    bins[i] = v;
    i++;
  }
  delete tokens;
  if (isUnit) {
    if (bins.GetSize() > 1)
      SetAxis(axis, Int_t(bins[1]-bins[0]), bins[0], bins[1]);
    else
      SetAxis(axis, 2*Int_t(bins[0]), bins[0]);
  }      
  else if (isRange) {
    Int_t    nBins = Int_t(bins[0]);
    if (bins.GetSize() > 2) 
      SetAxis(axis, nBins, bins[1], bins[2]);
    else
      SetAxis(axis, nBins, bins[1]);
  }
  else 
    SetAxis(axis, bins.GetSize()-1,bins.GetArray());
}
//____________________________________________________________________
void AliBaseAODTask::SetAxis(TAxis&   axis,
			     Int_t    n,
			     Double_t l,
			     Double_t h)
{
  axis.Set(n, l, h);
  FixAxis(axis);
}
//____________________________________________________________________
void AliBaseAODTask::SetAxis(TAxis& axis, Int_t n, Double_t m)
{
  SetAxis(axis, n, -TMath::Abs(m), +TMath::Abs(m));
}
//____________________________________________________________________
void AliBaseAODTask::PrintAxis(const TAxis& axis,
			       Int_t nSig,
			       const char* alt)
{
  gROOT->IndentLevel();
  TString tit(alt ? alt : axis.GetTitle());
  tit.Append(" axis:");
  printf(" %-26s ", tit.Data());
  if (axis.GetXbins() && axis.GetXbins()->GetArray()) {
    printf("%.*f", nSig, axis.GetBinLowEdge(1));
    for (Int_t i = 1; i <= axis.GetNbins(); i++)
      printf(":%.*f", nSig, axis.GetBinUpEdge(i));
  }
  else
    printf("%d bins between %.*f and %.*f",
	   axis.GetNbins(), nSig, axis.GetXmin(),nSig,axis.GetXmax());
  printf("\n");
}

//________________________________________________________________________
void
AliBaseAODTask::SetCentralityAxis(UShort_t n, Short_t* bins)
{
  DGUARD(fDebug,3,"Set centrality axis, %d bins", n);
  TArrayD dbins(n+1);
  for (UShort_t i = 0; i <= n; i++)
    dbins[i] = (bins[i] == 100 ? 100.1 : bins[i]);
  SetAxis(fCentAxis, n, dbins.GetArray());
}

//________________________________________________________________________
void
AliBaseAODTask::SetCentralityAxis(const char* bins)
{
  DGUARD(fDebug,3,"Set centrality axis: %s", bins);
  if (!bins || bins[0] == '\0') return;

  TString     spec(bins);
  if (spec.EqualTo("none", TString::kIgnoreCase))
    return;

  TArrayD edges;
  if (spec.EqualTo("default", TString::kIgnoreCase) ||
      spec.EqualTo("pbpb", TString::kIgnoreCase)) {
    //                 1  2  3   4   5   6   7   8   9   10  11  12
    Double_t tmp[] = { 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };
    SetAxis(fCentAxis,11,tmp);
    return;
  }
  if (spec.EqualTo("ppb", TString::kIgnoreCase) ||
      spec.EqualTo("pbp", TString::kIgnoreCase)) {
    //                 1  2  3   4   5   6   7   8
    Double_t tmp[] = { 0, 5, 10, 20, 40, 60, 80, 100 };
    SetAxis(fCentAxis,7, tmp);
    return;
  }

  SetAxis(fCentAxis,bins);
}

//________________________________________________________________________
void
AliBaseAODTask::SetCentralityAxis(UShort_t n, Double_t* bins)
{
  DGUARD(fDebug,3,"Set centrality axis, %d bins", n);
  SetAxis(fCentAxis, n, bins);
}
//________________________________________________________________________
void
AliBaseAODTask::SetCentralityAxis(Short_t low, Short_t high)
{
  // Short_t a[] = { low, high };
  SetAxis(fCentAxis, (high-low), low, high);
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
  fSums->SetOwner(true);

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
  fVertex->SetFillStyle(3002);
  fVertex->SetLineColor(kRed+2);
  fSums->Add(fVertex);
  fAccVertex = static_cast<TH1D*>(fVertex->Clone("vertexAcc"));
  fAccVertex->SetTitle("IP_{z} of accepted events");
  fAccVertex->SetDirectory(0);
  fAccVertex->SetFillColor(kGreen+2);
  fAccVertex->SetLineColor(kGreen+2);
  fSums->Add(fAccVertex);

  fAccVertexXY = new TH2D("vertexAccXY", "IP_{x,y} of accepted events",
			  1000,-2,2,1000,-2,2);
  fAccVertexXY->SetXTitle("IP_{x} [cm]");
  fAccVertexXY->SetYTitle("IP_{y} [cm]");
  fAccVertexXY->SetDirectory(0);
  fSums->Add(fAccVertexXY);
  
  fCent = new TH1D("cent","Centrality of all events",102, -1, 101);
  fCent->SetXTitle("Centrality [%]");
  fCent->SetYTitle("Events");
  fCent->SetFillColor(kRed+2);
  fCent->SetFillStyle(3002);
  fCent->SetLineColor(kRed+2);
  fCent->SetDirectory(0);
  fSums->Add(fCent);
  fAccCent = static_cast<TH1D*>(fCent->Clone("centAcc"));
  fAccCent->SetTitle("Centrality of accepted events");
  fAccCent->SetDirectory(0);
  fAccCent->SetFillColor(kGreen+2);
  fAccCent->SetLineColor(kGreen+2);
  fSums->Add(fAccCent);

  // Store centrality axis as a histogram - which can be merged
  TH1* cH = 0;
  if (fCentAxis.GetXbins() && fCentAxis.GetXbins()->GetSize() > 0)
    cH = new TH1I(fCentAxis.GetName(), fCentAxis.GetTitle(),
      fCentAxis.GetNbins(), fCentAxis.GetXbins()->GetArray());
  else
    cH = new TH1I(fCentAxis.GetName(), fCentAxis.GetTitle(),
      fCentAxis.GetNbins(), fCentAxis.GetXmin(),
      fCentAxis.GetXmax());
  cH->SetBinContent(1,1);
  cH->GetXaxis()->SetTitle(fCentAxis.GetTitle());
  cH->GetXaxis()->SetName(fCentAxis.GetName());

  fSums->Add(cH);

  fSums->Add(AliForwardUtil::MakeParameter("trigger", ULong_t(fTriggerMask)));
  fSums->Add(AliForwardUtil::MakeParameter("filter",  ULong_t(fFilterMask)));
  fSums->Add(AliForwardUtil::MakeParameter("count", 1));
  fSums->Add(AliForwardUtil::MakeParameter("alirootRev",
             AliForwardUtil::AliROOTRevision()));
  fSums->Add(AliForwardUtil::MakeParameter("alirootBranch",
             AliForwardUtil::AliROOTBranch()));



  if (!Book()) AliFatalF("Failed to book output objects for %s", GetName());

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
AliAODMultEventClass*
AliBaseAODTask::GetMultClass(const AliAODEvent& aod, Bool_t verb)
{
  // Get the forward object that contains our event selection stuff
  TObject* obj = aod.FindListObject("MultClass");
  if (!obj) {
    if (verb) AliWarning("No multiplicity event class object found");
    return 0;
  }
  AliAODMultEventClass* multClass = static_cast<AliAODMultEventClass*>(obj);
  return multClass;
}
//____________________________________________________________________
Double_t
AliBaseAODTask::GetCentrality(AliAODEvent&,
			      AliAODForwardMult* forward,
			      Int_t&             qual)
{
  qual          = 0;
  Double_t cent = forward->GetCentrality();
  Double_t max  = (HasCentrality() ? fCentAxis.GetXmax() : 100);
  if (cent < 0)   { cent = -.5; qual = 0xFFFF; }
  if (cent > max) { cent = TMath::Max(max+.1,100.5); qual = 198; }
  return cent;
}
//____________________________________________________________________
Double_t
AliBaseAODTask::GetCentrality(AliAODEvent& event,
			      AliAODForwardMult* forward)
{
  Int_t    qual = 0;
  Double_t cent = GetCentrality(event, forward, qual);
  if (qual > 0)   forward->SetTriggerBits(AliAODForwardMult::kCentNoCalib);
  return cent;
}

//____________________________________________________________________
Double_t
AliBaseAODTask::GetIpZ(AliAODEvent&,
		       AliAODForwardMult* forward)
{
  return forward->GetIpZ();
}
//____________________________________________________________________
Bool_t
AliBaseAODTask::GetIpXY(AliAODEvent& aod, Double_t& x, Double_t& y)
			
{
  x = -10000;
  y = -10000;
  AliVVertex* gen   = aod.GetPrimaryVertex();
  AliVVertex* vtx[] = { aod.GetPrimaryVertexSPD(),
			aod.GetPrimaryVertexTPC(),
			gen };
  Bool_t ret = false;
  for (Int_t i = 0; i < 3; i++) { 
    if (!vtx[i] || (vtx[i] != gen && !vtx[i]->IsFromVertexer3D()))
      continue;
	
    x = vtx[i]->GetX();
    y = vtx[i]->GetY();
    ret = true;
    break;
  }
  return ret;
}

//____________________________________________________________________
AliAODCentralMult*
AliBaseAODTask::GetCentral(const AliAODEvent& aod, Bool_t mc, Bool_t verb)
{
  // Get the central object that contains our event selection stuff
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
  if (!obj) return 0;
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
  DGUARD(fDebug,1,"Analyse the AOD event in UserExec");
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
  Double_t vtx   = GetIpZ(*aod, forward);
  Float_t  cent  = GetCentrality(*aod, forward);
  fVertex->Fill(vtx);
  fCent->Fill(cent);
  if (HasCentrality() && fAbsMinCent >= 0 && cent < fAbsMinCent) return;

  // Now check our event selectio up front
  if (!CheckEvent(*forward)) return;

  // Let user defined code do the job
  Bool_t taken = Event(*aod);

  // Fill our histograms
  if (taken) {
    fAccVertex->Fill(vtx);
    fAccCent->Fill(cent);

    Double_t ipX, ipY;
    GetIpXY(*aod, ipX,  ipY);
    fAccVertexXY->Fill(ipX, ipY);
  }

  PostData(1, fSums);

  PostEvent();
}

//____________________________________________________________________
Bool_t
AliBaseAODTask::CheckEvent(const AliAODForwardMult& forward)
{
  if (HasCentrality()) {
    return forward.CheckEvent(fTriggerMask,
			      fIPzAxis.GetXmin(),
			      fIPzAxis.GetXmax(),
			      fCentAxis.GetXmin(),
			      fCentAxis.GetXmax(),
			      fTriggers,
			      fEventStatus,
			      fFilterMask);
}
  return forward.CheckEvent(fTriggerMask, 
			    fIPzAxis.GetXmin(),
			    fIPzAxis.GetXmax(),
			    0,
			    0,
			    fTriggers,
			    fEventStatus,
			    fFilterMask);
}

//____________________________________________________________________
void
AliBaseAODTask::StoreInformation(AliAODForwardMult& forward)
{
  fSums->Add(AliForwardUtil::MakeParameter("sNN", forward.GetSNN()));
  fSums->Add(AliForwardUtil::MakeParameter("sys", forward.GetSystem()));
}

//____________________________________________________________________
Bool_t
AliBaseAODTask::StoreTrainName(Int_t no)
{
  AliAnalysisDataSlot* slot = GetOutputSlot(no);
  if (!slot) {
    AliWarningF("No output slot # %d defined", no);
    return false;
  }

  AliAnalysisDataContainer* cont = slot->GetContainer();
  if (!cont) {
    AliWarningF("Output slot # %d has no container", no);
    return false;
  }

  TFile* file = cont->GetFile();
  if (!file) {
    AliWarningF("Container of output slot %d has no file", no);
    return false;
  }
  if (!file->IsWritable()) {
    AliWarningF("File (%s) associated with output slot # %d not writable",
		file->GetName(), no);
    return false;
  }

  TDirectory* save  = gDirectory;
  TNamed* tag = new TNamed("trainName",
			   AliAnalysisManager::GetAnalysisManager()->GetName());
  file->cd();
  tag->Write();
  AliInfoF("Wrote train name (%s) to output file %s",
	   tag->GetTitle(), file->GetName());
  save->cd();

  return true;
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
    fResults = static_cast<TList*>(fSums->Clone(resName.Data()));
  else {
    fResults = new TList;
    fResults->SetName(resName.Data());
  }
  fResults->SetOwner(true);

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
  PFV("Trigger mask",   AliAODForwardMult::GetTriggerString(fTriggerMask));
  PFV("Filter mask",    AliAODForwardMult::GetTriggerString(fFilterMask,"|"));
  PFV("Abs.least cent", fAbsMinCent);
  PrintAxis(fIPzAxis,  1, "IP z");
  PrintAxis(fCentAxis, 0, "Centraltiy");
  gROOT->DecreaseDirLevel();
}
//
// EOF
//
