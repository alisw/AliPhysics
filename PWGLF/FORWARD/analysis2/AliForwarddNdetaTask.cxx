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
  DGUARD(fDebug, 3, "Default CTOR of AliForwarddNdetaTask");
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
  SetTitle("FMD");
  DGUARD(fDebug, 3, "Named CTOR of AliForwarddNdetaTask");
}

//____________________________________________________________________
AliForwarddNdetaTask::AliForwarddNdetaTask(const AliForwarddNdetaTask& o)
  : AliBasedNdetaTask(o)
{
  // 
  // Copy constructor
  // 
  DGUARD(fDebug, 3, "Copy CTOR of AliForwarddNdetaTask");
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
  DGUARD(fDebug, 3,"Make a centrality bin for AliForwarddNdetaTask: %s [%d,%d]",
	 name, l, h);
  return new AliForwarddNdetaTask::CentralityBin(name, l, h);
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
AliForwarddNdetaTask::CentralityBin::End(TList*      sums, 
					 TList*      results,
					 UShort_t    scheme,
					 const TH2F* shapeCorr, 
					 Double_t    trigEff,
					 Double_t    trigEff0,
					 Bool_t      symmetrice,
					 Int_t       rebin, 
					 Bool_t      rootProj,
					 Bool_t      corrEmpty, 
					 Bool_t      cutEdges,
					 Int_t       triggerMask,
					 Int_t       marker,
					 Int_t       color,
					 TList*      mclist,
					 TList*      truthlist )
{
  DGUARD(fDebug, 1, "In End of %s with corrEmpty=%d, cutEdges=%d, rootProj=%d", 
	 GetName(), corrEmpty, cutEdges, rootProj);
  AliBasedNdetaTask::CentralityBin::End(sums, results, scheme, 
					shapeCorr, trigEff, trigEff0,
					symmetrice, rebin, 
					rootProj, corrEmpty, cutEdges,
					triggerMask, marker, color, mclist, 
					truthlist);

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
#if 0
  // TEMPORARY FIX
  if (triggerMask == AliAODForwardMult::kNSD) {
    // This is a local change 
    epsilonT = 0.92; 
    AliWarning(Form("Using hard-coded NSD trigger efficiency of %f",epsilonT));
  }
#endif
  AliInfo("Adding per-ring histograms to output");
  TString text;
  Double_t scaler = Normalization(*fTriggers, scheme, epsilonT, ntotal, &text);
  TIter next(res->GetHists());
  TH1*  hist = 0;
  while ((hist = static_cast<TH1*>(next()))) hist->Scale(scaler);
  res->SetName("dndetaRings");
  fOutput->Add(res);
  fOutput->Add(new TNamed("normCalc", text.Data()));
}

//________________________________________________________________________
//
// EOF
//
