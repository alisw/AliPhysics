//====================================================================
#include "AliForwarddNdetaTask.h"
#include <TMath.h>
#include <TH2D.h>
#include <TH1D.h>
#include <THStack.h>
#include <TList.h>
#include <TFile.h>
#include <AliAnalysisManager.h>
#include <AliAODEvent.h>
#include <AliAODHandler.h>
#include <AliAODInputHandler.h>
#include "AliForwardUtil.h"
#include "AliAODForwardMult.h"

//____________________________________________________________________
AliForwarddNdetaTask::AliForwarddNdetaTask()
  : AliBasedNdetaTask()
{
  //
  // Constructor 
  // 
}

//____________________________________________________________________
AliForwarddNdetaTask::AliForwarddNdetaTask(const char* /* name */)
  : AliBasedNdetaTask("Forward")
{
  // 
  // Constructor
  // 
  // Paramters
  //   name    Name of task 
}

//____________________________________________________________________
AliForwarddNdetaTask::AliForwarddNdetaTask(const AliForwarddNdetaTask& o)
  : AliBasedNdetaTask(o)
{
  // 
  // Copy constructor
  // 
}

//____________________________________________________________________
AliBasedNdetaTask::CentralityBin*
AliForwarddNdetaTask::MakeCentralityBin(const char* name, Short_t l, Short_t h) 
  const 
{
  // 
  // Make a new centrality bin
  // 
  // Parameters:
  //    name   Histogram names
  //    l      Lower cut
  //    h      Upper cut
  // 
  // Return:
  //    Newly allocated object (of our type)
  //
  return new AliForwarddNdetaTask::CentralityBin(name, l, h);
}

//____________________________________________________________________
void
AliForwarddNdetaTask::UserExec(Option_t* option)
{
  // 
  // Called at each event 
  //
  // This is called once in the master 
  // 
  // Parameters:
  //    option Not used 
  //
  AliBasedNdetaTask::UserExec(option);
  
  AliAODEvent* aod = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!aod) {
    AliError("Cannot get the AOD event");
    return;
  } 

  TObject* obj = aod->FindListObject("Forward");
  if (!obj) { 
    AliWarning("No forward object found");
    return;
  }
  AliAODForwardMult* forward = static_cast<AliAODForwardMult*>(obj);
 
  TObject* oPrimary   = aod->FindListObject("primary");
  if (!oPrimary) return;
  
  TH2D* primary   = static_cast<TH2D*>(oPrimary);

  // Loop over centrality bins 
  TIter next(fListOfCentralities);
  CentralityBin* bin = 0;
  while ((bin = static_cast<CentralityBin*>(next()))) 
    bin->ProcessPrimary(forward, fTriggerMask, fVtxMin, fVtxMax, primary);
}  
  
//________________________________________________________________________
void 
AliForwarddNdetaTask::Terminate(Option_t *option) 
{
  // 
  // Called at end of event processing.. 
  //
  // This is called once in the master 
  // 
  // Parameters:
  //    option Not used 
  AliBasedNdetaTask::Terminate(option);

  THStack* truth      = new THStack("dndetaTruth",      "dN/d#eta MC Truth");
  THStack* truthRebin = new THStack(Form("dndetaTruth_rebin%02d", fRebin), 
				    "dN/d#eta MC Truth");

  TIter next(fListOfCentralities);
  CentralityBin* bin = 0;
  while ((bin = static_cast<CentralityBin*>(next()))) {
    if (fCentAxis && bin->IsAllBin()) continue;

    TList* results = bin->GetResults();
    if (!results) continue; 

    TH1* dndeta      = static_cast<TH1*>(results->FindObject("dndetaTruth"));
    TH1* dndetaRebin = 
      static_cast<TH1*>(results->FindObject(Form("dndetaTruth_rebin%02d",
						fRebin)));
    if (dndeta)      truth->Add(dndeta);
    if (dndetaRebin) truthRebin->Add(dndetaRebin);
  }
  // If available output rebinned stack 
  if (!truth->GetHists() || 
      truth->GetHists()->GetEntries() <= 0) {
    AliWarning("No MC truth histograms found");
    delete truth;
    truth = 0;
  }
  if (truth) fOutput->Add(truth);

  // If available output rebinned stack 
  if (!truthRebin->GetHists() || 
      truthRebin->GetHists()->GetEntries() <= 0) {
    AliWarning("No rebinned MC truth histograms found");
    delete truthRebin;
    truthRebin = 0;
  }
  if (truthRebin) fOutput->Add(truthRebin);
  
}

//____________________________________________________________________
TH2D*
AliForwarddNdetaTask::GetHistogram(const AliAODEvent* aod, Bool_t mc)
{
  // 
  // Retrieve the histogram 
  // 
  // Parameters:
  //    aod AOD event 
  //    mc  Whether to get the MC histogram or not
  // 
  // Return:
  //    Retrieved histogram or null
  //
  TObject* obj = 0;
  if (mc) obj = aod->FindListObject("ForwardMC");
  else    obj = aod->FindListObject("Forward");

  // We should have a forward object at least 
  if (!obj) {
    if (!mc) AliWarning("No Forward object found AOD");
    return 0;
  }
  AliAODForwardMult* forward = static_cast<AliAODForwardMult*>(obj);
  return &(forward->GetHistogram());
}

//========================================================================
void
AliForwarddNdetaTask::CentralityBin::ProcessPrimary(const AliAODForwardMult* 
						    forward, 
						    Int_t triggerMask,
						    Double_t vzMin, 
						    Double_t vzMax, 
						    const TH2D* primary)
{ 
  // Check the centrality class unless this is the 'all' bin 
  if (!IsAllBin()) { 
    Double_t centrality = forward->GetCentrality();
    if (centrality < fLow || centrality >= fHigh) return;
  }

  // Create sum histogram 
  if (!fSumPrimary) { 
    fSumPrimary = static_cast<TH2D*>(primary->Clone("truth"));
    fSumPrimary->SetDirectory(0);
    fSumPrimary->Reset();
    fSums->Add(fSumPrimary);
  }
  
  // translate real trigger mask to MC trigger mask
  Int_t mask = AliAODForwardMult::kB;
  if (triggerMask == AliAODForwardMult::kNSD) {
    mask ^= AliAODForwardMult::kNSD;
    mask =  AliAODForwardMult::kMCNSD;
  }

  // Now use our normal check, but with the new mask, except 
  vzMin = vzMax = -10000; // ignore vertex 
  if (!forward->CheckEvent(mask, vzMin, vzMax, 0, 0, 0)) return;

  fSumPrimary->Add(primary);
  Int_t n = fSumPrimary->GetBinContent(0,0);
  fSumPrimary->SetBinContent(0,0, ++n);
}

//________________________________________________________________________
void
AliForwarddNdetaTask::CentralityBin::End(TList*      sums, 
					 TList*      results,
					 UShort_t    scheme,
					 const TH1*  shapeCorr, 
					 Double_t    trigEff,
					 Bool_t      symmetrice,
					 Int_t       rebin, 
					 Bool_t      rootProj,
					 Bool_t      corrEmpty, 
					 Bool_t      cutEdges,
					 Int_t       triggerMask,
					 Int_t       marker)
{
  AliInfo(Form("In End of %s with corrEmpty=%d, cutEdges=%d, rootProj=%d", 
	       GetName(), corrEmpty, cutEdges, rootProj));
  AliBasedNdetaTask::CentralityBin::End(sums, results, scheme, 
					shapeCorr, trigEff, 
					symmetrice, rebin, 
					rootProj, corrEmpty, cutEdges,
					triggerMask, marker);

  fSumPrimary     = static_cast<TH2D*>(fSums->FindObject("truth"));

  if (!fSumPrimary) { 
    AliWarning("No MC truth histogram found");
  }
  else {
#if 0
    Int_t n = fSumPrimary->GetBinContent(0,0);
#else
    Int_t n = (triggerMask == AliAODForwardMult::kNSD ? 
	       Int_t(fTriggers->GetBinContent(AliAODForwardMult::kBinMCNSD)) : 
	       Int_t(fTriggers->GetBinContent(AliAODForwardMult::kBinAll)));
#endif
    AliInfo(Form("Normalising MC truth to %d", n));
    
    TH1D* dndetaTruth = fSumPrimary->ProjectionX("dndetaTruth",1,
					       fSumPrimary->GetNbinsY(),"e");
    dndetaTruth->SetDirectory(0);
    dndetaTruth->Scale(1./n, "width");

    
    SetHistogramAttributes(dndetaTruth, GetColor(), 30, "Monte-Carlo truth");
    
    fOutput->Add(dndetaTruth);
    fOutput->Add(Rebin(dndetaTruth, rebin, cutEdges));
  }

  if (!IsAllBin()) return;
  TFile* file = TFile::Open("forward.root", "READ");
  if (!file) return;
  
  TList* forward = static_cast<TList*>(file->Get("Forward"));
  if (!forward) { 
    AliError("List Forward not found in forward.root");
    return;
  }
  TList* rings = static_cast<TList*>(forward->FindObject("ringResults"));
  if (!rings) { 
    AliError("List ringResults not found in forward.root");
    return;
  }
  THStack* res = static_cast<THStack*>(rings->FindObject("all"));
  if (!res) { 
    AliError(Form("Stack all not found in %s", rings->GetName()));
    return;
  }
  if (!fTriggers) { 
    AliError("Triggers histogram not set");
    return;
  }
  Double_t ntotal   = 0;
  Double_t epsilonT = trigEff;
  // TEMPORARY FIX
  if (triggerMask == AliAODForwardMult::kNSD) {
    // This is a local change 
    epsilonT = 0.92; 
    AliWarning(Form("Using hard-coded NSD trigger efficiency of %f",epsilonT));
  }
  AliInfo("Adding per-ring histograms to output");
  Double_t scaler = Normalization(*fTriggers, scheme, epsilonT, ntotal);
  TIter next(res->GetHists());
  TH1*  hist = 0;
  while ((hist = static_cast<TH1*>(next()))) hist->Scale(scaler);
  res->SetName("dndetaRings");
  fOutput->Add(res);
}

//________________________________________________________________________
//
// EOF
//
