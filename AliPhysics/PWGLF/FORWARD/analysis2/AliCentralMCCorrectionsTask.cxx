// 
// Calculate the corrections in the central regions
// 
// Inputs: 
//   - AliESDEvent 
//
// Outputs: 
//   - AliAODCentralMult 
// 
// Histograms 
//   
// Corrections used 
// 
#include "AliCentralMCCorrectionsTask.h"
#include "AliCentralCorrectionManager.h"
#include "AliTriggerAnalysis.h"
#include "AliPhysicsSelection.h"
#include "AliLog.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliESDEvent.h"
#include "AliAODHandler.h"
#include "AliMultiplicity.h"
#include "AliInputEventHandler.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliAODForwardMult.h"
#include "AliCentralCorrSecondaryMap.h"
#include "AliCentralCorrAcceptance.h"
#include "AliForwardUtil.h"
#include "AliMultiplicity.h"
#include <TVector3.h>
#include <TH1.h>
#include <TH2D.h>
#include <TDirectory.h>
#include <TList.h>
#include <TROOT.h>
#include <iostream>

//====================================================================
AliCentralMCCorrectionsTask::AliCentralMCCorrectionsTask()
  : AliBaseMCCorrectionsTask(),
    fTrackDensity(),
    fSecCorr(0), 
    fAccCorr(0),
    fNPhiBins(20),
    fEffectiveCorr(true),
    fEtaCut(1.9),
    fCorrCut(0.6)
{
  // 
  // Constructor 
  // 
  // Parameters:
  //    name Name of task 
  //
  DGUARD(fDebug, 3,"Default CTOR of AliCentralMCCorrectionsTask");
}

//____________________________________________________________________
AliCentralMCCorrectionsTask::AliCentralMCCorrectionsTask(const char* name)
  : AliBaseMCCorrectionsTask(name, &(AliCentralCorrectionManager::Instance())),
    fTrackDensity("trackDensity"),
    fSecCorr(0), 
    fAccCorr(0),
    fNPhiBins(20),
    fEffectiveCorr(true),
    fEtaCut(1.9),
    fCorrCut(0.6)
{
  // 
  // Constructor 
  // 
  // Parameters:
  //    name Name of task 
  //
  DGUARD(fDebug, 3,"Named CTOR of AliCentralMCCorrectionsTask: %s",name);
}

//____________________________________________________________________
AliBaseMCCorrectionsTask::VtxBin*
AliCentralMCCorrectionsTask::CreateVtxBin(Double_t low, Double_t high)
{
  return new AliCentralMCCorrectionsTask::VtxBin(low,high,fEtaAxis,fNPhiBins);
}

//____________________________________________________________________
Bool_t
AliCentralMCCorrectionsTask::ProcessESD(const AliESDEvent& esd, 
					const AliMCEvent& mc, 
					AliBaseMCCorrectionsTask::VtxBin& bin,
					const TVector3& ip)
{
  AliCentralMCCorrectionsTask::VtxBin& vb = 
    static_cast<AliCentralMCCorrectionsTask::VtxBin&>(bin);

  // Now process our input data and store in own ESD object 
  fTrackDensity.Calculate(mc, ip, *vb.fHits, bin.fPrimary);
  
  // Get the ESD object
  const AliMultiplicity* spdmult = esd.GetMultiplicity();
 
  // Count number of tracklets per bin 
  for(Int_t j = 0; j< spdmult->GetNumberOfTracklets();j++) 
    vb.fClusters->Fill(spdmult->GetEta(j),spdmult->GetPhi(j));
  //...and then the unused clusters in layer 1 
  for(Int_t j = 0; j< spdmult->GetNumberOfSingleClusters();j++) {
    Double_t eta = -TMath::Log(TMath::Tan(spdmult->GetThetaSingle(j)/2.));
    vb.fClusters->Fill(eta, spdmult->GetPhiSingle(j));
  }
  
  // Count events  
  bin.fCounts->Fill(0.5);
  
  return true;
}
//____________________________________________________________________
void
AliCentralMCCorrectionsTask::CreateCorrections(TList* results)
{
  fSecCorr = new AliCentralCorrSecondaryMap;
  fSecCorr->SetVertexAxis(fVtxAxis);

  fAccCorr = new AliCentralCorrAcceptance;
  fAccCorr->SetVertexAxis(fVtxAxis);

  results->Add(fSecCorr);
  results->Add(fAccCorr);
}

//____________________________________________________________________
Bool_t 
AliCentralMCCorrectionsTask::FinalizeVtxBin(AliBaseMCCorrectionsTask::VtxBin* 
					    bin,  UShort_t iVz) 
{
  
  AliCentralMCCorrectionsTask::VtxBin* vb = 
    static_cast<AliCentralMCCorrectionsTask::VtxBin*>(bin);
  vb->Terminate(fList, fResults, iVz, fEffectiveCorr, 
		fEtaCut, fCorrCut, fSecCorr,fAccCorr);
  return true;
}
					    

//____________________________________________________________________
void
AliCentralMCCorrectionsTask::Print(Option_t* option) const
{
  AliBaseMCCorrectionsTask::Print(option);
  std::cout << "  # of phi bins:    " << fNPhiBins << "\n"
	    << "  Effective corr.:  " << fEffectiveCorr << "\n"
	    << "  Eta cut-off:      " << fEtaCut << "\n"
	    << "  Acceptance cut:   " << fCorrCut 
	    << std::endl;
  gROOT->IncreaseDirLevel();
  fTrackDensity.Print(option);
  gROOT->DecreaseDirLevel();
}

//====================================================================
AliCentralMCCorrectionsTask::VtxBin::VtxBin()
  : AliBaseMCCorrectionsTask::VtxBin(),
    fHits(0), 
    fClusters(0)
{
}
//____________________________________________________________________
AliCentralMCCorrectionsTask::VtxBin::VtxBin(Double_t     low, 
					    Double_t     high, 
					    const TAxis& axis,
					    UShort_t     nPhi)
  : AliBaseMCCorrectionsTask::VtxBin(low, high, axis, nPhi),
    fHits(0), 
    fClusters(0)
{
  fHits = static_cast<TH2D*>(fPrimary->Clone("hits"));
  fHits->SetTitle("Hits");
  fHits->SetDirectory(0);

  fClusters = static_cast<TH2D*>(fPrimary->Clone("clusters"));
  fClusters->SetTitle("Clusters");
  fClusters->SetDirectory(0);
}


//____________________________________________________________________
TList*
AliCentralMCCorrectionsTask::VtxBin::CreateOutputObjects(TList* l)
{
  TList* d = AliBaseMCCorrectionsTask::VtxBin::CreateOutputObjects(l);

  d->Add(fHits);
  d->Add(fClusters);

  return d;
}
//____________________________________________________________________
void
AliCentralMCCorrectionsTask::VtxBin::Terminate(const TList* input, 
					       TList* output, 
					       UShort_t iVz, 
					       Bool_t effectiveCorr,
					       Double_t etaCut,
					       Double_t accCut,
					       AliCentralCorrSecondaryMap* map,
					       AliCentralCorrAcceptance* acorr)
{
  TList* out = new TList;
  out->SetName(GetName());
  out->SetOwner();
  output->Add(out);

  TList* l = static_cast<TList*>(input->FindObject(GetName()));
  if (!l) { 
    AliError(Form("List %s not found in %s", GetName(), input->GetName()));
    return;
  }

  // Get the sums 
  TH2D*   hits  = static_cast<TH2D*>(l->FindObject("hits"));
  TH2D*   clus  = static_cast<TH2D*>(l->FindObject("clusters"));
  TH2D*   prim  = static_cast<TH2D*>(l->FindObject("primary"));
  if (!hits || !prim) {
    AliError(Form("Missing histograms: %p, %p", hits, prim));
    return;
  }

  // Clone cluster and hit map
  TH2D* secMapEff = static_cast<TH2D*>(clus->Clone("secMapEff"));
  TH2D* secMapHit = static_cast<TH2D*>(hits->Clone("secMapHit"));
  secMapEff->SetTitle("2^{nd} map from clusters");
  secMapEff->SetDirectory(0);
  secMapHit->SetTitle("2^{nd} map from MC hits");
  secMapHit->SetDirectory(0);

  // Divide cluster and hit map with number of primaries 
  secMapEff->Divide(prim);
  secMapHit->Divide(prim);

  // Create acceptance histograms 
  TH1D* accEff = new TH1D("accEff",
			  "Acceptance correction for SPD (from 2^{nd} map)" ,
			  fPrimary->GetXaxis()->GetNbins(), 
			  fPrimary->GetXaxis()->GetXmin(), 
			  fPrimary->GetXaxis()->GetXmax());
  TH1D* accHit = static_cast<TH1D*>(accEff->Clone("accHit"));
  accHit->SetTitle("Acceptance correction for SPD (from clusters)");

  // Diagnostics histogra, 
  TH2*  dia    = static_cast<TH2D*>(clus->Clone("diagnostics"));
  dia->SetTitle("Scaled cluster density");

  // Total number of channels along phi and # of eta bins
  Int_t nTotal = secMapHit->GetNbinsY();
  Int_t nEta   = secMapHit->GetNbinsX();

  for(Int_t xx = 1; xx <= nEta; xx++) {
    Double_t eta = secMapHit->GetXaxis()->GetBinCenter(xx);
    Bool_t   ins = TMath::Abs(eta) <= etaCut;
    Double_t mm  = 0;
    if (ins) {
      // Find the maximum cluster signal in this phi range 
      for (Int_t yy = 1; yy <= nTotal; yy++) { 
	Double_t c = clus->GetBinContent(xx,yy);
	mm         = TMath::Max(mm, c);
      }
    }
    // Count number of phi bins with enough clusters or high enough
    // correction.
    Int_t nOKEff    = 0;
    Int_t nOKHit    = 0;
    for(Int_t yy = 1; yy <=nTotal; yy++) {
      if (!ins) { // Not inside Eta cut
	secMapEff->SetBinContent(xx,yy,0.); 
	secMapEff->SetBinError(xx,yy,0.); 
	secMapHit->SetBinContent(xx,yy,0.); 
	secMapHit->SetBinError(xx,yy,0.); 
	dia->SetBinContent(xx,yy,0);
	continue;
      }

      // Check if the background correction is large enough, or zero map
      if(secMapEff->GetBinContent(xx,yy) > 0.9) {
	// acc->Fill(h->GetXaxis()->GetBinCenter(xx));
	nOKEff++;
      }
      else {
	secMapEff->SetBinContent(xx,yy,0.); 
	secMapEff->SetBinError(xx,yy,0.); 
      }

      // Check if the number of cluster is large enough, or zero map
      Double_t c = clus->GetBinContent(xx,yy);
      Double_t s = (mm < 1e-6) ? 0 : c / mm;
      dia->SetBinContent(xx,yy,s);
      if (s >= accCut) {
	nOKHit++;
      }
      else {
	secMapHit->SetBinContent(xx,yy,0);
	secMapHit->SetBinError(xx,yy,0);
      }
    }

    // Calculate acceptance as ratio of bins with enough clusters and
    // total number of phi bins.
    Double_t accXX = float(nOKHit) / nTotal;
    if (accXX < 0.2) accXX = 0;
    accHit->SetBinContent(xx, accXX);

    // Calculate acceptance as ratio of bins with large enough
    // correction and total number of phi bins.
    accXX = float(nOKEff) / nTotal;
    if (accXX < 0.2) accXX = 0;
    accEff->SetBinContent(xx, accXX);
  }

  TH2D* secMap    = (effectiveCorr ? secMapEff : secMapHit);
  TH2D* secMapAlt = (effectiveCorr ? secMapHit : secMapEff);
  TH1D* acc       = (effectiveCorr ? accEff    : accHit);
  TH1D* accAlt    = (effectiveCorr ? accHit    : accEff);
  out->Add(secMap->Clone("secMap"));
  out->Add(secMapAlt->Clone());
  out->Add(acc->Clone("acc"));
  out->Add(accAlt->Clone());
  out->Add(dia->Clone());

  map->SetCorrection(iVz, secMap);
  acorr->SetCorrection(iVz, acc);
}

//
// EOF
//
