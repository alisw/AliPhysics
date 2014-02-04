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
  // SetTitle("FMD");
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
AliForwarddNdetaTask::MakeCentralityBin(const char* name, Short_t l,Short_t h) 
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
  DGUARD(fDebug, 3,
	 "Make a centrality bin for AliForwarddNdetaTask: %s [%d,%d]",
	 name, l, h);
  return new AliForwarddNdetaTask::CentralityBin(name, l, h);
}


//____________________________________________________________________
TH2D*
AliForwarddNdetaTask::GetHistogram(const AliAODEvent& aod, Bool_t mc)
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
  // We should have a forward object at least 
  AliAODForwardMult* forward = GetForward(aod, mc, !mc);
  if (!forward) return 0;
  return &(forward->GetHistogram());
}
//____________________________________________________________________
void
AliForwarddNdetaTask::CheckEventData(Double_t vtx, 
				     TH2*     data, 
				     TH2*     dataMC)
{
  // Check if this is satellite
  // if (!fSatelliteVertices) return;
  Double_t aVtx = TMath::Abs(vtx);
  if (aVtx < 37.5 || aVtx > 400) return;

  TH2* hists[] = { data, dataMC };

  // In satellite vertices FMD2i is cut away manually at this point
  // for certain vertices. It could be done in the ESDs, but as of
  // this writing not for specific vertices.
  // 
  // cholm comment: It would be difficult to setup the filter in the
  // reconstruction pass, but it could perhaps be done in the AOD
  // filtering.
  // 
  // This is what was done for
  // the Pb-Pb paper (arXiv:1304.0347).
  for (Int_t iX = 0; iX<=data->GetNbinsX(); iX++) {
    // Do all checks up front - as soon as we can - branching is
    // expensive!
    Double_t x    = data->GetXaxis()->GetBinCenter(iX);
    Bool_t   zero = false;
    if (((vtx >  60 && vtx <  90) && x < 3) ||
	((vtx > 330 && vtx < 350) && x > -2.5) ||
	((vtx < 100 || vtx > 305) && TMath::Abs(x) < 4.5) || 
	(vtx < 50                 && TMath::Abs(x) < 4.75))
      zero = true;
    if (!zero) continue;
    
    for (Int_t iH = 0; iH < 2; iH++) {
      if (!hists[iH]) continue;
      // if (iX > hists[iH]->GetNbinsX()+1) continue;
      // Also zero coverage and phi acceptance for this 
      for (Int_t iY = 0; iY<=hists[iH]->GetNbinsY()+1; iY++) {	  
	hists[iH]->SetBinContent(iX, iY, 0);
	hists[iH]->SetBinError(iX, iY, 0);
      }
    }
  }

  if (fCorrEmpty) {
    // Now, since we have some dead areas in FMD2i (sectors 16 and
    // 17), we need to remove the corresponding bins from the
    // histogram. However, it is not obvious which bins (in eta) to
    // remove, so remove everything starting from the most negative to
    // the middle of the histogram.
    // 
    // This hack was first introduced by HHD, but was done at the end of
    // the event processing (CentralityBin::MakeResults).  That is,
    // however, not very practical, as we'd like to normalize to the phi
    // acceptance rather than the eta coverage and then correct for
    // empty bins. Since the only way to really update the phi
    // acceptance stored in the overflow bin is on the event level, we
    // should really do it here.
    const Int_t phiBin1 = 17; // Sector 16
    const Int_t phiBin2 = 18; // Sector 17
    for (Int_t iH = 0; iH < 2; iH++) { 
      if (!hists[iH]) continue;
      
      Int_t midX = hists[iH]->GetNbinsX() / 2;
      // Int_t nY   = hists[iH]->GetNbinsY();
      for (Int_t i = 1; i <= midX; i++) { 
	hists[iH]->SetBinContent(i, phiBin1, 0);
	hists[iH]->SetBinContent(i, phiBin2, 0);
	hists[iH]->SetBinError(i, phiBin1, 0);
	hists[iH]->SetBinError(i, phiBin2, 0);
	
	// Here, we should also modify the overflow bin to reflect the
	// new phi acceptance.  First get the old phi acceptance -
	// then multiply this on the number of bins. This gives us -
	// roughly - the number of sectors we had.  Then take out two
	// from that number, and then calculate the new phi
	// Acceptance. Note, if the sectors where already taken out in
	// the AOD production, we _will_ end up with a wrong number,
	// so we should _not_ do that in the AOD production.  This is
	// tricky and may not work at all.  For now, we should rely on
	// the old way of correcting to the eta coverage and
	// correcting for empty bins.
      }
    }
  }
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
  DGUARD(fDebug, 1,"In End of %s with corrEmpty=%d, cutEdges=%d, rootProj=%d", 
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
  
  TList* forward = static_cast<TList*>(file->Get("ForwardSums"));
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
