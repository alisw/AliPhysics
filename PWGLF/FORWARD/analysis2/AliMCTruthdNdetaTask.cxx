//====================================================================
#include "AliMCTruthdNdetaTask.h"
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
AliMCTruthdNdetaTask::AliMCTruthdNdetaTask()
  : AliBasedNdetaTask(),
    fHasData(true)
{
  //
  // Constructor 
  // 
}

//____________________________________________________________________
AliMCTruthdNdetaTask::AliMCTruthdNdetaTask(const char* /* name */)
  : AliBasedNdetaTask("MCTruth"), 
    fHasData(true)
{
  // 
  // Constructor
  // 
  // Paramters
  //   name    Name of task 
}


//____________________________________________________________________
AliBasedNdetaTask::CentralityBin*
AliMCTruthdNdetaTask::MakeCentralityBin(const char* name, Float_t l, Float_t h) 
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
  return new AliMCTruthdNdetaTask::CentralityBin(name, l, h);
}

//____________________________________________________________________
TH2D*
AliMCTruthdNdetaTask::GetHistogram(const AliAODEvent& aod, Bool_t mc)
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
  if (!fHasData) return 0;
  if (mc) return 0;

  TH2D* ret = GetPrimary(aod);
  if (!ret) {
    fHasData = false;
    return 0;
  }
  Int_t nY  = ret->GetNbinsY();
  // Need to fill under-/overflow bin with 1's 
  for (Int_t i = 1; i <= ret->GetNbinsX(); i++)  {
    ret->SetBinContent(i, 0,    1);
    ret->SetBinContent(i, nY+1, 1);
  }
  return ret;
}

//________________________________________________________________________
Bool_t
AliMCTruthdNdetaTask::Finalize() 
{
  // 
  // Called at end of event processing.. 
  //
  // This is called once in the master 
  // 
  // Parameters:
  //    option Not used 
  if (!fHasData) {
    AliInfo("The MC truth dN/deta task didn't get any data");
    return false;
  }
  AliBasedNdetaTask::Finalize();

  THStack*       truth = new THStack("dndetaTruth", "dN/d#eta MC Truth");
  CentralityBin* bin   = 0;
  TIter          next(fListOfCentralities);
  while ((bin = static_cast<CentralityBin*>(next()))) {
    if (HasCentrality() && bin->IsAllBin()) continue;

    TList* results = bin->GetResults();
    if (!results) continue; 

    TH1* dndeta      = static_cast<TH1*>(results->FindObject("dndetaTruth"));
    if (dndeta)      truth->Add(dndeta);
  }
  // If available output rebinned stack 
  if (!truth->GetHists() || 
      truth->GetHists()->GetEntries() <= 0) {
    AliWarning("No MC truth histograms found");
    delete truth;
    truth = 0;
  }
  if (truth) fResults->Add(truth);

  return true;
}

//========================================================================
Bool_t
AliMCTruthdNdetaTask::CentralityBin::ProcessEvent(const AliAODForwardMult* 
						  forward, 
						  UInt_t triggerMask,
						  Bool_t isZero,
						  Double_t vzMin, 
						  Double_t vzMax, 
						  const TH2D* primary,
						  const TH2D*,
						  UInt_t filter,
						  Double_t weight)
{ 
  // Check the centrality class unless this is the 'all' bin 
  if (!primary) return false;

  if (!IsAllBin()) { 
    Double_t centrality = forward->GetCentrality();
    if (centrality < fLow || centrality >= fHigh) return false;
  }

  if (!fSum) CreateSums(primary, 0);
  if (!fSumTruth) { 
    fSumTruth = static_cast<TH2D*>(primary->Clone("truth"));
    fSumTruth->SetDirectory(0);
    fSumTruth->Reset();
    fSums->Add(fSumTruth);
  }

  // translate real trigger mask to MC trigger mask
  Int_t mask = AliAODForwardMult::kB;
  if (triggerMask == AliAODForwardMult::kNSD) {
    mask ^= AliAODForwardMult::kNSD;
    mask =  AliAODForwardMult::kMCNSD;
  }

  // Now use our normal check, but with the new mask, except ignore vertex
  if (forward->CheckEvent(mask, -10000, -10000, 0, 0, 0, 0, filter)) {
    fSumTruth->Add(primary, weight);

    // Store event count in left-most under- underflow bin 
    Int_t cnt = Int_t(fSumTruth->GetBinContent(0,0));
    fSumTruth->SetBinContent(0,0, ++cnt);
  }

  // Now use our normal check with the full trigger mask and vertex
  if (CheckEvent(forward, triggerMask, vzMin, vzMax, filter)) 
    fSum->Add(primary, isZero, weight);
  return true;
}

//________________________________________________________________________
bool
AliMCTruthdNdetaTask::CentralityBin::End(TList*      sums, 
					 TList*      results,
					 UShort_t    scheme,
					 Double_t    trigEff,
					 Double_t    trigEff0,
					 Bool_t      rootProj,
					 Bool_t      corrEmpty, 
					 Int_t       triggerMask,
					 Int_t       marker,
					 Int_t       color,
					 TList*      mclist,
					 TList*      truthlist)
{
#if 0
  AliInfo(Form("At end with sums=%p, results=%p, scheme=%d, "
	       "shapeCorr=%p, trigEff=%f, symmetrice=%d, rebin=%d, "
	       "rootProj=%d, corrEmpty=%d, cutEdges=%d, triggerMask=0x%08x, "
	       "marker=%d (%d)", 
	       sums, results, scheme, shapeCorr, trigEff, trigEff0, symmetrice, 
	       rebin, rootProj, corrEmpty, cutEdges, triggerMask, marker,
	       GetMarkerStyle(kStar)));
#endif

  if (!AliBasedNdetaTask::CentralityBin::End(sums, results, scheme, trigEff, 
					     trigEff0, rootProj, corrEmpty,
					     triggerMask,
					     marker, color, mclist, 
					     truthlist))
    return false;

  fSumTruth     = static_cast<TH2D*>(fSums->FindObject("truth"));
  

  if (fSumTruth) { 
    Int_t n0 = Int_t(fSumTruth->GetBinContent(0,0));
    Int_t n  = (triggerMask == AliAODForwardMult::kNSD ? 
		Int_t(fTriggers->GetBinContent(AliAODForwardMult::kBinMCNSD)) : 
		Int_t(fTriggers->GetBinContent(AliAODForwardMult::kBinAll)));
    DMSG(fDebug,0,"Normalising MC truth to %d (%d additions)", n, n0);
    
    TH1D* dndetaTruth = fSumTruth->ProjectionX("dndetaTruth",1,
					       fSumTruth->GetNbinsY(),"e");
    dndetaTruth->SetDirectory(0);
    dndetaTruth->Scale(1./n, "width");
    
    SetHistogramAttributes(dndetaTruth, GetColor(color)+1, 
			   GetMarkerStyle(kCross), 
			   "Monte-Carlo truth");

    fOutput->Add(dndetaTruth);
  }
  TH1* dndeta          =                       GetResult("");
  if (dndeta)    
    dndeta->SetTitle("Monte-Carlo truth (selected)");
  return true;
}

//________________________________________________________________________
//
// EOF
//
