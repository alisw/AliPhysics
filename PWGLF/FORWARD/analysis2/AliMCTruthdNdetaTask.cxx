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
  : AliBasedNdetaTask()
{
  //
  // Constructor 
  // 
}

//____________________________________________________________________
AliMCTruthdNdetaTask::AliMCTruthdNdetaTask(const char* /* name */)
  : AliBasedNdetaTask("MCTruth")
{
  // 
  // Constructor
  // 
  // Paramters
  //   name    Name of task 
}

//____________________________________________________________________
AliMCTruthdNdetaTask::AliMCTruthdNdetaTask(const AliMCTruthdNdetaTask& o)
  : AliBasedNdetaTask(o)
{
  // 
  // Copy constructor
  // 
}

//____________________________________________________________________
AliBasedNdetaTask::CentralityBin*
AliMCTruthdNdetaTask::MakeCentralityBin(const char* name, Short_t l, Short_t h) 
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
AliMCTruthdNdetaTask::GetHistogram(const AliAODEvent* aod, Bool_t mc)
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
  if (mc) return 0;
  TObject* obj = aod->FindListObject("primary");
  // We should have a forward object at least 
  if (!obj) return 0;
  TH2D* ret = static_cast<TH2D*>(obj);
  // Need to fill underflow bin with 1's 
  for (Int_t i = 1; i <= ret->GetNbinsX(); i++)  
    ret->SetBinContent(i, 0, 1);
  return ret;
}

//________________________________________________________________________
void 
AliMCTruthdNdetaTask::Terminate(Option_t *option) 
{
  // 
  // Called at end of event processing.. 
  //
  // This is called once in the master 
  // 
  // Parameters:
  //    option Not used 
  AliBasedNdetaTask::Terminate(option);

  THStack* truth      = new THStack("dndetaTruth", "dN/d#eta MC Truth");
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

//========================================================================
void
AliMCTruthdNdetaTask::CentralityBin::ProcessEvent(const AliAODForwardMult* 
						  forward, 
						  Int_t triggerMask,
						  Bool_t isZero,
						  Double_t vzMin, 
						  Double_t vzMax, 
						  const TH2D* primary,
						  const TH2D*)
{ 
  // Check the centrality class unless this is the 'all' bin 
  if (!primary) return;

  if (!IsAllBin()) { 
    Double_t centrality = forward->GetCentrality();
    if (centrality < fLow || centrality >= fHigh) return;
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
  if (forward->CheckEvent(mask, -10000, -10000, 0, 0, 0)) {
    fSumTruth->Add(primary);

    // Store event count in left-most under- underflow bin 
    Int_t cnt = Int_t(fSumTruth->GetBinContent(0,0));
    fSumTruth->SetBinContent(0,0, ++cnt);
  }

  // Now use our normal check with the full trigger mask and vertex
  if (CheckEvent(forward, triggerMask, vzMin, vzMax)) 
    fSum->Add(primary, isZero);
}

//________________________________________________________________________
void
AliMCTruthdNdetaTask::CentralityBin::End(TList*      sums, 
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
	       sums, results, scheme, shapeCorr, trigEff, symmetrice, 
	       rebin, rootProj, corrEmpty, cutEdges, triggerMask, marker,
	       GetMarkerStyle(kStar)));
#endif

  AliBasedNdetaTask::CentralityBin::End(sums, results, scheme, 
					shapeCorr, trigEff, 
					symmetrice, rebin, 
					rootProj, corrEmpty, cutEdges,
					triggerMask, marker, color, mclist, 
					truthlist);

  fSumTruth     = static_cast<TH2D*>(fSums->FindObject("truth"));
  

  if (fSumTruth) { 
    Int_t n0 = Int_t(fSumTruth->GetBinContent(0,0));
    Int_t n  = (triggerMask == AliAODForwardMult::kNSD ? 
		Int_t(fTriggers->GetBinContent(AliAODForwardMult::kBinMCNSD)) : 
		Int_t(fTriggers->GetBinContent(AliAODForwardMult::kBinAll)));
    AliInfo(Form("Normalising MC truth to %d (%d additions)", n, n0));
    
    TH1D* dndetaTruth = fSumTruth->ProjectionX("dndetaTruth",1,
					       fSumTruth->GetNbinsY(),"e");
    dndetaTruth->SetDirectory(0);
    dndetaTruth->Scale(1./n, "width");
    
    SetHistogramAttributes(dndetaTruth, GetColor(color)+1, 
			   GetMarkerStyle(kCross), 
			   "Monte-Carlo truth");

    fOutput->Add(dndetaTruth);
    fOutput->Add(Rebin(dndetaTruth, rebin, cutEdges));
  }
  TH1* dndeta          = GetResult(0,     false, "");
  TH1* dndetaSym       = GetResult(0,     true,  "");
  TH1* dndeta_rebin    = GetResult(rebin, false, "");
  TH1* dndetaSym_rebin = GetResult(rebin, true,  "");
  if (dndeta)    
    dndeta->SetTitle("Monte-Carlo truth (selected)");
  if (dndetaSym) 
    dndetaSym->SetTitle("Monte-Carlo truth (selected,mirrored)");
  if (dndeta_rebin)    
    dndeta_rebin->SetTitle("Monte-Carlo truth (selected)");
  if (dndetaSym_rebin) 
    dndetaSym_rebin->SetTitle("Monte-Carlo truth (selected,mirrored)");

}

//________________________________________________________________________
//
// EOF
//
