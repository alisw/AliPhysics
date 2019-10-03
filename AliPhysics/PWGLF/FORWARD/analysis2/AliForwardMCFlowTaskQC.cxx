//
// Calculate flow on MC data in the forward and central regions using the Q cumulants method.
//
// Inputs:
//  - AliAODEvent
//
// Outputs:
//  - AnalysisResults.root
//
#include "AliForwardMCFlowTaskQC.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "TGraph.h"
#include "TF1.h"
#include "AliAODEvent.h"
#include "AliAODForwardMult.h"
#include "AliAODCentralMult.h"
#include "AliGenEventHeader.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliGenEventHeaderTunedPbPb.h"
#include "AliCollisionGeometry.h"

ClassImp(AliForwardMCFlowTaskQC)
#if 0
;
#endif
//_____________________________________________________________________
AliForwardMCFlowTaskQC::AliForwardMCFlowTaskQC() 
  : AliForwardFlowTaskQC(), 
    fBinsForwardTR(),      // List of FMDTR analysis objects
    fBinsCentralTR(),      // List of SPDTR analysis objects
    fBinsMC(),             // List of MC particle analysis objects
    fAODMCHeader(),        // MC Header
    fHistdNdedpMC(),       // MC particle d^2N/detadphi histogram
    fHistFMDMCCorr(),      // FMD MC correlation
    fHistSPDMCCorr(),      // SPD MC correlation
    fWeights(0),            // Flow weights
    fImpactParToCent(),    // Impact parameter to centrality graph
    fUseImpactPar(0),      // Use impact par for centrality
    fUseMCVertex(0),       // Get vertex from MC header?
    fUseFlowWeights(0)     // Add flow to MC particles
{} 
  //
  // Default Constructor
  //
//_____________________________________________________________________
AliForwardMCFlowTaskQC::AliForwardMCFlowTaskQC(const char* name) 
  : AliForwardFlowTaskQC(name),
    fBinsForwardTR(),      // List of FMDTR analysis objects
    fBinsCentralTR(),      // List of SPDTR analysis objects
    fBinsMC(),             // List of MC particles analysis objects
    fAODMCHeader(0),       // MC Header
    fHistdNdedpMC(),       // MC particles d^2N/detadphi histogram
    fHistFMDMCCorr(),      // FMD MC correlation
    fHistSPDMCCorr(),      // SPD MC correlation
    fWeights(0),            // Flow weights
    fImpactParToCent(),    // Impact parameter to centrality graph
    fUseImpactPar(0),      // Use impact par for centrality
    fUseMCVertex(0),       // Get vertex from MC header?
    fUseFlowWeights(0)     // Add flow to MC particles
{ 
  // 
  // Constructor
  // Parameters:
  //  name: Name of task
  //
  
  //  Double_t impactParam[] = {0.,1.75,4.225,5.965,7.765,9.215,10.46,
  //                            11.565,12.575,13.515,16.679};
  //  Double_t centrality[] = {0.,2.5,7.5,15,25,35,45,55,65,75,90};
  Double_t impactParam[] = { 0.00,  3.72,  5.23,  7.31,  8.88, 10.20, 
			    11.38, 12.47, 13.50, 14.51, 16.679};
  Double_t centrality[]  = { 0.,    5.,   10.,   20.,   30.,   40., 
			    50.,   60.,   70.,   80.,  100.};

  Int_t nPoints = sizeof(impactParam)/sizeof(Double_t);
  fImpactParToCent = new TGraph(nPoints, impactParam, centrality);

}
//_____________________________________________________________________
AliForwardMCFlowTaskQC::AliForwardMCFlowTaskQC(const AliForwardMCFlowTaskQC& o) 
  : AliForwardFlowTaskQC(o), 
    fBinsForwardTR(),                      // List of FMDTR analysis objects
    fBinsCentralTR(),                      // List of SPDTR analysis objects
    fBinsMC(),                             // List of MC particles analysis objects
    fAODMCHeader(o.fAODMCHeader),          // MC Header
    fHistdNdedpMC(o.fHistdNdedpMC),        // MC particles d^2N/detadphi histogram
    fHistFMDMCCorr(o.fHistFMDMCCorr),      // FMD MC correlation
    fHistSPDMCCorr(o.fHistSPDMCCorr),      // SPD MC correlation
    fWeights(o.fWeights),                  // Flow weights
    fImpactParToCent(o.fImpactParToCent),  // Impact parameter to centrality
    fUseImpactPar(o.fUseImpactPar),        // Use impact par for centrality
    fUseMCVertex(o.fUseMCVertex),          // Get vertex from MC header?
    fUseFlowWeights(o.fUseFlowWeights)     // Add flow to MC particles
{
  //
  // Copy Constructor
  //
} 
//_____________________________________________________________________
AliForwardMCFlowTaskQC&
AliForwardMCFlowTaskQC::operator=(const AliForwardMCFlowTaskQC& o)
{
  // 
  // Assignment operator
  // Parameters:
  //    o Object to copy from 
  //
  if (&o == this) return *this;
  fAODMCHeader     = o.fAODMCHeader;
  fHistdNdedpMC    = o.fHistdNdedpMC;
  fHistFMDMCCorr   = o.fHistFMDMCCorr;
  fHistSPDMCCorr   = o.fHistSPDMCCorr;
  fWeights         = o.fWeights;
  fImpactParToCent = o.fImpactParToCent;
  fUseImpactPar    = o.fUseImpactPar;
  fUseMCVertex     = o.fUseMCVertex;
  fUseFlowWeights  = o.fUseFlowWeights;
  return *this;
}
//_____________________________________________________________________
void AliForwardMCFlowTaskQC::InitVertexBins()
{
  //
  // Initiate VertexBin lists
  //
  AliForwardFlowTaskQC::InitVertexBins();

  Bool_t isNUA = (fFlowFlags & kNUAcorr);
  for (Int_t v = 1; v <= fVtxAxis->GetNbins(); v++) {
    Int_t vL = Int_t(fVtxAxis->GetBinLowEdge(v));
    Int_t vH = Int_t(fVtxAxis->GetBinUpEdge(v));
    // FMD
    if ((fFlowFlags & kFMD)) {
      fBinsForwardTR.Add(new VertexBin(vL, vH, fMaxMoment, "FMDTR", fFlowFlags, fFMDCut, fEtaGap));
      if (!(fFlowFlags & k3Cor)) 
        fBinsCentralTR.Add(new VertexBin(vL, vH, fMaxMoment, "SPDTR", fFlowFlags|kSPD, fSPDCut, fEtaGap));
      if (isNUA) fFlowFlags ^= kNUAcorr;
      fBinsMC.Add(new VertexBin(vL, vH, fMaxMoment, "MC-FMD", fFlowFlags|kMC, -1, fEtaGap));
      if ((fFlowFlags & kStdQC)) 
	fBinsMC.Add(new VertexBin(vL, vH, fMaxMoment, "MC-SPD", fFlowFlags|kMC|kSPD, -1, fEtaGap));
      if (isNUA) fFlowFlags ^= kNUAcorr;
    }
    // VZERO
    else if ((fFlowFlags & kVZERO)) {
      fBinsMC.Add(new VertexBin(vL, vH, fMaxMoment, "MC-VZERO", fFlowFlags|kMC, -1, fEtaGap));
    }
  }
}
//_____________________________________________________________________
void AliForwardMCFlowTaskQC::InitHists()
{
  //
  // Initiate diagnostics hists and add to outputlist
  //
  AliForwardFlowTaskQC::InitHists();

  TString subDetName = ((fFlowFlags & kFMD) ? "FMD" : ((fFlowFlags & kVZERO) ? "VZERO" : "none"));
  fHistdNdedpMC = TH2D(Form("fdNdedpMC%s%s", subDetName.Data(), GetQCType(fFlowFlags)),
		   Form("fdNdedpMC%s%s", subDetName.Data(), GetQCType(fFlowFlags)),
		   240, -6., 6., 200, 0., TMath::TwoPi());
  
  fHistFMDMCCorr = new TH2D("hFMDMCCorr", "hFMDMCCorr", 200, 0., 15000., 200, 0, 20000);
  fHistSPDMCCorr = new TH2D("hSPDMCCorr", "hSPDMCCorr", 200, 0., 7500., 200, 0, 20000);
  TList* dList = (TList*)fSumList->FindObject("Diagnostics");
  if (!dList) {
    dList = new TList();
    dList->SetName("Diagnostics");
    fSumList->Add(dList);
  }
  dList->Add(fHistFMDMCCorr);
  dList->Add(fHistSPDMCCorr);

  TIter nextForwardTR(&fBinsForwardTR);
  VertexBin* bin = 0;
  while ((bin = static_cast<VertexBin*>(nextForwardTR()))) {
    bin->AddOutput(fSumList, fCentAxis);
  }
  TIter nextCentralTR(&fBinsCentralTR);
  while ((bin = static_cast<VertexBin*>(nextCentralTR()))) {
    bin->AddOutput(fSumList, fCentAxis);
  }
  TIter nextMC(&fBinsMC);
  while ((bin = static_cast<VertexBin*>(nextMC()))) {
    bin->AddOutput(fSumList, fCentAxis);
  }
  if (fWeights) {
    TList* wList = new TList();
    wList->SetName("FlowWeights");
    fWeights->Init(wList);
    fSumList->Add(wList);
  }

}
//_____________________________________________________________________
Bool_t AliForwardMCFlowTaskQC::Analyze() 
{
  // 
  // Load FMD and SPD MC objects from aod tree and call the cumulants 
  // calculation for the correct vertexbin
  //
  if (!AliForwardFlowTaskQC::Analyze()) return kFALSE;

  // Run analysis on trackrefs from FMD and SPD
  AliAODForwardMult* aodfmult = static_cast<AliAODForwardMult*>(fAOD->FindListObject("ForwardMC"));
  AliAODCentralMult* aodcmult = static_cast<AliAODCentralMult*>(fAOD->FindListObject("CentralClustersMC"));
  Int_t vtx = fVtxAxis->FindBin(fVtx)-1;
  
  // if objects are present, get histograms
  if (aodfmult) {
    TH2D& fmdTRdNdetadphi = aodfmult->GetHistogram();
    if ((fFlowFlags & kStdQC)) {
      FillVtxBinList(fBinsForwardTR, fmdTRdNdetadphi, vtx);
    } else if ((fFlowFlags & kEtaGap)) {
      FillVtxBinListEtaGap(fBinsForwardTR, fmdTRdNdetadphi, fmdTRdNdetadphi, vtx/*, kDoVtxCut*/);
    }
    if (aodcmult) {
      TH2D& spdTRdNdetadphi = aodcmult->GetHistogram();
      if ((fFlowFlags & kStdQC)) {
	FillVtxBinList(fBinsCentralTR, spdTRdNdetadphi, vtx);
      } else if ((fFlowFlags & kEtaGap)) {
	FillVtxBinListEtaGap(fBinsCentralTR, fmdTRdNdetadphi, spdTRdNdetadphi, vtx/*, kDoVtxCut*/);
      } else if ((fFlowFlags & k3Cor)) {
	FillVtxBinList3Cor(fBinsForwardTR, spdTRdNdetadphi, fmdTRdNdetadphi, vtx);
      }
    }
  }
  // Run analysis on MC branch
  if (!FillMCHist()) return kFALSE;

  if ((fFlowFlags & kStdQC)) {
    FillVtxBinList(fBinsMC, fHistdNdedpMC, vtx, kMC);
  } else if ((fFlowFlags & kEtaGap)) {
    FillVtxBinListEtaGap(fBinsMC, fHistdNdedpMC, fHistdNdedpMC, vtx, kMC);
  } else if ((fFlowFlags & k3Cor)) {
    FillVtxBinList3Cor(fBinsMC, fHistdNdedpMC, fHistdNdedpMC, vtx);
  }

  // Mult correlation diagnostics
  if (aodfmult && aodcmult) {
    AliAODForwardMult* fmult = static_cast<AliAODForwardMult*>(fAOD->FindListObject("Forward"));
    AliAODCentralMult* cmult = static_cast<AliAODCentralMult*>(fAOD->FindListObject("CentralClusters"));
    const TH2D& fhist = fmult->GetHistogram();
    const TH2D& chist = cmult->GetHistogram();

    Double_t totForward = fhist.Integral(1, fhist.GetNbinsX(), 1, fhist.GetNbinsY());
    Double_t totSPD = chist.Integral(1, chist.GetNbinsX(), 1, chist.GetNbinsY());
    Double_t totMC = fHistdNdedpMC.Integral(1, fHistdNdedpMC.GetNbinsX(), 1, fHistdNdedpMC.GetNbinsY());
    fHistFMDMCCorr->Fill(totForward, totMC);
    fHistSPDMCCorr->Fill(totSPD, totMC);
  }

  return kTRUE;
}
//_____________________________________________________________________
void AliForwardMCFlowTaskQC::Finalize() 
{
  //
  // Finalize command, called by Terminate()
  //
  AliForwardFlowTaskQC::Finalize();

  EndVtxBinList(fBinsForwardTR);
  EndVtxBinList(fBinsCentralTR);
  EndVtxBinList(fBinsMC);

  return;
}
//_____________________________________________________________________
Bool_t AliForwardMCFlowTaskQC::CheckEvent(const AliAODForwardMult* aodfm) 
{
  // 
  //  Function to check that an AOD event meets the cuts
  //
  //  Parameters: 
  //   AliAODForwardMult: forward mult object with trigger and vertex info
  //
  //  Return: false if there is no trigger or if the centrality or vertex
  //  is out of range. Otherwise true.
  //
  fAODMCHeader = static_cast<AliAODMCHeader*>(fAOD->FindListObject(AliAODMCHeader::StdBranchName()));
  return AliForwardFlowTaskQC::CheckEvent(aodfm);
}
//_____________________________________________________________________
Bool_t AliForwardMCFlowTaskQC::CheckTrigger(const AliAODForwardMult* aodfm) const 
{
  //
  // Function to look for a trigger string in the event.
  //
  // Parameters: 
  //  AliAODForwardMult: forward mult object with trigger and vertex info
  //
  // Returns true if B trigger is present - for some reason this is the one we use in MC
  //
  if (aodfm) return aodfm->IsTriggerBits(AliAODForwardMult::kB);
  else return (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))
                 ->IsEventSelected() & AliVEvent::kMB);

}
// _____________________________________________________________________
Bool_t AliForwardMCFlowTaskQC::GetCentrality(const AliAODForwardMult* aodfm)
{
  // 
  // Function to use centrality parametrization from impact parameter
  // if flag is not set call AliForwardFlowTaskQC::GetCentrality
  //
  // Parameters:
  //  AliAODForwardMult: forward mult object with trigger and vertex info
  //
  // Returns true when centrality is set.
  //
  AliGenEventHeaderTunedPbPb* header = 
    dynamic_cast<AliGenEventHeaderTunedPbPb*>(fAODMCHeader->GetCocktailHeader(0));
  if (header) fCent = header->GetCentrality();
  else if (fUseImpactPar) fCent = GetCentFromB();
  else return AliForwardFlowTaskQC::GetCentrality(aodfm);
  
  if (fCentAxis->GetXmin() > fCent || fCent >= fCentAxis->GetXmax()) {
    fHistEventSel->Fill(kInvCent);
    return kFALSE;
  }
  if (fCent != 0) return kTRUE;
  else {
    fHistEventSel->Fill(kNoCent);
    return kFALSE;
  }
}
//_____________________________________________________________________
Bool_t AliForwardMCFlowTaskQC::GetVertex(const AliAODForwardMult* aodfm)
{
  //
  // Function to look for vertex determination in the event using the MC header.
  //
  // Parameters: 
  //  AliAODForwardMult: Not used
  //
  // Returns true if vertex is determined
  //
  if (fUseMCVertex) {
    if (fAODMCHeader) {
      fVtx = fAODMCHeader->GetVtxZ();
      if (fVtx < fVtxAxis->GetXmin() || fVtx > fVtxAxis->GetXmax()) {
        fHistEventSel->Fill(kInvVtx);
	return kFALSE;
      }
      return kTRUE;
    } else {
      fHistEventSel->Fill(kNoVtx);
      return kFALSE;
    }
  } else 
    return AliForwardFlowTaskQC::GetVertex(aodfm);
}
//_____________________________________________________________________
Bool_t AliForwardMCFlowTaskQC::FillMCHist()  
{
  // 
  // Loop over AliAODParticle branch and fill d^2N/detadphi-histograms.
  // Add flow if set to do so in AddTask function
  //
  fHistdNdedpMC.Reset();
  Double_t minEta = -3.75;
  Double_t maxEta = 5.;

  //retreive MC particles from event
  TClonesArray* mcArray = 
    static_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
  if(!mcArray){
    AliWarning("No MC array found in AOD. Try making it again.");
    return kFALSE;
  }

  if (!fAODMCHeader) AliWarning("No MC header found.");

  Int_t ntracks = mcArray->GetEntriesFast();
  Double_t rp = -1, b = -1;
  if (fAODMCHeader) {
    rp = fAODMCHeader->GetReactionPlaneAngle();
    b = fAODMCHeader->GetImpactParameter();
    if (fAODMCHeader->GetNCocktailHeaders() > 1) {
      ntracks = fAODMCHeader->GetCocktailHeader(0)->NProduced();
    }
  }
  
  for (Int_t it = 0; it < ntracks; it++) { // Track loop
    Double_t weight = 1;
    AliAODMCParticle* particle = (AliAODMCParticle*)mcArray->At(it);
    if (!particle) {
      AliError(Form("Could not receive track %d", it));
      continue;
    }
    if (!particle->IsPrimary()) continue;
    if (particle->Charge() == 0) continue;
    Double_t pT = particle->Pt();
    Double_t eta = particle->Eta();
    Double_t phi = particle->Phi();
    if (eta >= minEta && eta < maxEta) {
      // Add flow if it is in the argument
      if (fUseFlowWeights && fWeights) { 
        weight = fWeights->CalcWeight(eta, pT, phi, particle->PdgCode(), rp, b); 
//      Printf("%f", weight);
      }
      fHistdNdedpMC.Fill(eta, phi, weight);
    }
  }
  // Set underflow bins for acceptance checks
  Int_t sBin = fHistdNdedpMC.GetXaxis()->FindBin(minEta);
  Int_t eBin = fHistdNdedpMC.GetXaxis()->FindBin(maxEta);
  for ( ; sBin <= eBin; sBin++) {
    fHistdNdedpMC.SetBinContent(sBin, 0, 1);
  } // End of eta bin loop

  return kTRUE;
}
//_____________________________________________________________________
Double_t AliForwardMCFlowTaskQC::GetCentFromB() const
{
  //
  // Get centrality from MC impact parameter.
  //
  Double_t cent = -1.;
  if (!fAODMCHeader) return cent;
  Double_t b = fAODMCHeader->GetImpactParameter();
  cent = fImpactParToCent->Eval(b);

  return cent;
}
//_____________________________________________________________________
//
//
// EOF
