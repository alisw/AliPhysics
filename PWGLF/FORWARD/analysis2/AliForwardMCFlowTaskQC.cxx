//
// Calculate flow on MC data in the forward and central regions using the Q cumulants method.
//
// Inputs:
//  - AliAODEvent
//
// Outputs:
//  - AnalysisResults.root
//
/**
 * @file   AliForwardMCFlowTaskQC.cxx
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Thu Feb  7 01:09:19 2013
 * 
 * @brief  
 * 
 * 
 * @ingroup pwglf_forward_flow
 */
#include "AliForwardMCFlowTaskQC.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "TGraph.h"
#include "TF1.h"
#include "TProfile2D.h"
#include "AliAODEvent.h"
#include "AliAODForwardMult.h"
#include "AliAODCentralMult.h"
#include "AliGenEventHeader.h"

ClassImp(AliForwardMCFlowTaskQC)
#if 0
;
#endif
//_____________________________________________________________________
AliForwardMCFlowTaskQC::AliForwardMCFlowTaskQC() 
  : AliForwardFlowTaskQC(), 
    fBinsFMDTR(),          // List of FMDTR analysis objects
    fBinsSPDTR(),          // List of SPDTR analysis objects
    fBinsMC(),             // List of MC truth analysis objects
    fdNdedpMC(),           // MC truth d^2N/detadphi histogram
    fWeights(),            // Flow weights
    fImpactParToCent(),    // Impact parameter to centrality graph
    fUseImpactPar(0),      // Use impact par for centrality
    fFMDMinEta(-6),        // FMD min eta coverage for this vtx
    fFMDMaxEta(6),         // FMD max eta coverage for this vtx
    fAddFlow(0),           // Add flow to MC truth
    fAddType(0),           // Add type of flow to MC truth
    fAddOrder(0)           // Add order of flow to MC truth        
{} 
  //
  // Default Constructor
  //
//_____________________________________________________________________
AliForwardMCFlowTaskQC::AliForwardMCFlowTaskQC(const char* name) 
  : AliForwardFlowTaskQC(name),
    fBinsFMDTR(),          // List of FMDTR analysis objects
    fBinsSPDTR(),          // List of SPDTR analysis objects
    fBinsMC(),             // List of MC truth analysis objects
    fdNdedpMC(),           // MC truth d^2N/detadphi histogram
    fWeights(),            // Flow weights
    fImpactParToCent(),    // Impact parameter to centrality graph
    fUseImpactPar(0),      // Use impact par for centrality
    fFMDMinEta(-6),        // FMD min eta coverage for this vtx
    fFMDMaxEta(6),         // FMD max eta coverage for this vtx
    fAddFlow(0),           // Add flow to MC truth
    fAddType(0),           // Add type of flow to MC truth
    fAddOrder(0)           // Add order of flow to MC truth        
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
    fBinsFMDTR(),                          // List of FMDTR analysis objects
    fBinsSPDTR(),                          // List of SPDTR analysis objects
    fBinsMC(),                             // List of MC truth analysis objects
    fdNdedpMC(o.fdNdedpMC),                // MC truth d^2N/detadphi histogram
    fWeights(o.fWeights),                  // Flow weights
    fImpactParToCent(o.fImpactParToCent),  // Impact parameter to centrality
    fUseImpactPar(o.fUseImpactPar),        // Use impact par for centrality
    fFMDMinEta(o.fFMDMinEta),              // FMD min eta coverage for this vtx
    fFMDMaxEta(o.fFMDMaxEta),              // FMD max eta coverage for this vtx
    fAddFlow(o.fAddFlow),                  // Add flow to MC truth
    fAddType(o.fAddType),                  // Add type of flow to MC truth
    fAddOrder(o.fAddOrder)                 // Add order of flow to MC truth
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
  fdNdedpMC        = o.fdNdedpMC;
  fWeights         = o.fWeights;
  fImpactParToCent = o.fImpactParToCent;
  fUseImpactPar    = o.fUseImpactPar;
  fFMDMinEta       = o.fFMDMinEta;
  fFMDMaxEta       = o.fFMDMaxEta;
  fAddFlow         = o.fAddFlow;
  fAddType         = o.fAddType;
  fAddOrder        = o.fAddOrder;
  return *this;
}
//_____________________________________________________________________
void AliForwardMCFlowTaskQC::InitVertexBins()
{
  //
  // Initiate VertexBin lists
  //
  AliForwardFlowTaskQC::InitVertexBins();

  Int_t moment = 0;
  for(UShort_t n = 0; n < fV.GetSize(); n++) {
    moment = fV.At(n);
    for (Int_t v = 1; v <= fVtxAxis->GetNbins(); v++) {
      Int_t vL = Int_t(fVtxAxis->GetBinLowEdge(v));
      Int_t vH = Int_t(fVtxAxis->GetBinUpEdge(v));
      fBinsFMDTR.Add(new VertexBin(vL, vH, moment, "FMDTR", fFlowFlags, fFMDCut, fEtaGap));
      fBinsSPDTR.Add(new VertexBin(vL, vH, moment, "SPDTR", fFlowFlags, fSPDCut, fEtaGap));
      fBinsMC.Add(new VertexBin(vL, vH, moment, "MC", fFlowFlags, -1, fEtaGap));
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

  fdNdedpMC = TH2D(Form("fdNdedpMC%s", ((fFlowFlags & kEtaGap) ? "_etaGap" : "")),
		   Form("fdNdedpMC%s", ((fFlowFlags & kEtaGap) ? "_etaGap" : "")),
		   48, -6., 6., 200, 0., 2.*TMath::Pi());
  fdNdedpMC.Sumw2();

  TIter nextFMDTR(&fBinsFMDTR);
  VertexBin* bin = 0;
  while ((bin = static_cast<VertexBin*>(nextFMDTR()))) {
    bin->AddOutput(fSumList);
  }
  TIter nextSPDTR(&fBinsSPDTR);
  while ((bin = static_cast<VertexBin*>(nextSPDTR()))) {
    bin->AddOutput(fSumList);
  }
  TIter nextMC(&fBinsMC);
  while ((bin = static_cast<VertexBin*>(nextMC()))) {
    bin->AddOutput(fSumList);
  }

  TList* wList = new TList();
  wList->SetName("FlowWeights");
  fWeights.Init(wList);
  fSumList->Add(wList);

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
  const AliAODForwardMult* aodfmult = 
    static_cast<AliAODForwardMult*>(fAOD->FindListObject("ForwardMC"));
  const AliAODCentralMult* aodcmult = 
    static_cast<AliAODCentralMult*>(fAOD->FindListObject("CentralClustersMC"));
  Int_t vtx = fVtxAxis->FindBin(fVtx)-1;
  
  // if objects are present, get histograms
  if (aodfmult) {
    const TH2D& fmdTRdNdetadphi = aodfmult->GetHistogram();
    if ((fFlowFlags & kEtaGap)) {
      FillVtxBinListEtaGap(fBinsFMDTR, fmdTRdNdetadphi, fmdTRdNdetadphi, vtx);
    } else {
      FillVtxBinList(fBinsFMDTR, fmdTRdNdetadphi, vtx);
    }
    
    if (aodcmult) {
      const TH2D& spdTRdNdetadphi = aodcmult->GetHistogram();
      if ((fFlowFlags & kEtaGap)) {
	FillVtxBinListEtaGap(fBinsSPDTR, fmdTRdNdetadphi, spdTRdNdetadphi, vtx);
      } else {
	FillVtxBinList(fBinsSPDTR, spdTRdNdetadphi, vtx);
      }
    }
  }
  // Run analysis on MC branch
  if (!LoopAODMC()) return kFALSE;
  if ((fFlowFlags & kEtaGap)) {
    FillVtxBinListEtaGap(fBinsMC, fdNdedpMC, fdNdedpMC, vtx);
  } else {
    FillVtxBinList(fBinsMC, fdNdedpMC, vtx);
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

  EndVtxBinList(fBinsFMDTR);
  EndVtxBinList(fBinsSPDTR);
  EndVtxBinList(fBinsMC);

  return;
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
  return aodfm->IsTriggerBits(AliAODForwardMult::kB);
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
  if (fUseImpactPar) {
    fCent = GetCentFromB();
    if (fCent != -1) return kTRUE;
  }
  return AliForwardFlowTaskQC::GetCentrality(aodfm);
}
//_____________________________________________________________________
void AliForwardMCFlowTaskQC::GetFMDLimits()
{

  const AliAODForwardMult* aodfmult = 
    static_cast<AliAODForwardMult*>(fAOD->FindListObject("Forward"));
  const TH2D& h = aodfmult->GetHistogram();

  for (Int_t e = 1; ; e++) {
    if (h.GetBinContent(e, 0) != 0) { 
      fFMDMinEta = h.GetXaxis()->GetBinLowEdge(e);
      break;
    }
  }
  for (Int_t e = h.GetNbinsX(); ; e--) {
    if (h.GetBinContent(e, 0) != 0) { 
      fFMDMaxEta = h.GetXaxis()->GetBinLowEdge(e);
      break;
    }
  }

  return;
}
//_____________________________________________________________________
Bool_t AliForwardMCFlowTaskQC::LoopAODMC()  
{
  // 
  // Loop over AliAODParticle branch and fill d^2N/detadphi-histograms.
  // Add flow if set to do so in AddTask function
  //
  fdNdedpMC.Reset();
  GetFMDLimits();

  //retreive MC particles from event
  TClonesArray* mcArray = 
    static_cast<TClonesArray*>(fAOD->FindListObject(
                               AliAODMCParticle::StdBranchName()));
  if(!mcArray){
    AliWarning("No MC array found in AOD. Try making it again.");
    return kFALSE;
  }

  AliAODMCHeader* header = 
    dynamic_cast<AliAODMCHeader*>(fAOD->FindListObject(
                                  AliAODMCHeader::StdBranchName()));
  if (!header) 
    AliWarning("No header file found.");
  
  // Double_t rp = header->GetReactionPlaneAngle();

  Int_t ntracks = mcArray->GetEntriesFast();
  // TODO: Make this bit smarter...
  if (header->GetNCocktailHeaders() > 1) {
    ntracks = header->GetCocktailHeader(0)->NProduced();
  }

  UShort_t flowFlags = 0;
  if (fAddFlow.Length() > 1) {
    if (fAddFlow.Contains("pt"))   flowFlags |= AliForwardFlowWeights::kPt;
    if (fAddFlow.Contains("pid"))  flowFlags |= AliForwardFlowWeights::kPt;
    if (fAddFlow.Contains("eta"))  flowFlags |= AliForwardFlowWeights::kEta;
    if (fAddFlow.Contains("cent")) flowFlags |= AliForwardFlowWeights::kCent;
  }
  // Double_t b = header->GetImpactParameter();

  // Track loop: chek how many particles will be accepted
  Double_t weight = 0;
  for (Int_t it = 0; it < ntracks; it++) {
    weight = 1;
    AliAODMCParticle* particle = (AliAODMCParticle*)mcArray->At(it);
    if (!particle) {
      AliError(Form("Could not receive track %d", it));
      continue;
    }
    if (!particle->IsPrimary()) continue;
    if (particle->Charge() == 0) continue;
    //Double_t pT = particle->Pt();
    Double_t eta = particle->Eta();
    Double_t phi = particle->Phi();
    if (eta >= fFMDMinEta && eta <= fFMDMaxEta) {
      // Add flow if it is in the argument
      /* FLOW WEIGHTS DISABLED IN THE VERSION - COMING BACK SOON
      if (flowFlags != 0) { 
//	weight = fWeights.CalcWeight(eta, pT, phi, particle->PdgCode(), 
//				     rp, fCent, fAddType, fAddOrder, 
//				     flowFlags) + 1;
        weight = fWeights.CalcWeight(eta, pT, phi, particle->PdgCode(),
                                    rp, b); 
//      Printf("%f", weight);
      }*/
        fdNdedpMC.Fill(eta, phi, weight);
    }
  }
  Int_t sBin = fdNdedpMC.GetXaxis()->FindBin(fFMDMinEta);
  Int_t eBin = fdNdedpMC.GetXaxis()->FindBin(fFMDMaxEta);
  for ( ; sBin < eBin; sBin++) fdNdedpMC.SetBinContent(sBin, 0, 1.);

  return kTRUE;
}
//_____________________________________________________________________
Double_t AliForwardMCFlowTaskQC::GetCentFromB() const
{
  //
  // Get centrality from MC impact parameter.
  //
  Double_t cent = -1.;
  AliAODMCHeader* header = 
    static_cast<AliAODMCHeader*>(fAOD->FindListObject(AliAODMCHeader::
						      StdBranchName()));
  if (!header) return cent;
  Double_t b = header->GetImpactParameter();

  cent = fImpactParToCent->Eval(b);

  return cent;
}
//_____________________________________________________________________
void AliForwardMCFlowTaskQC::PrintFlowSetup() const  
{
  //
  // Print the setup of the flow task
  //
  Printf("AliForwardMCFlowTaskQC::Print");
  Printf("Number of bins in vertex axis:\t%d", fVtxAxis->GetNbins());
  Printf("Range of vertex axis         :\t[%3.1f,%3.1f]", 
			  fVtxAxis->GetXmin(), fVtxAxis->GetXmax());
  printf("Doing flow analysis for      :\t");
  for (Int_t n  = 0; n < fV.GetSize(); n++) printf("v%d ", fV.At(n));
  printf("\n");
  Printf("Satellite vertex flag           :\t%s", ((fFlowFlags & kSatVtx) ? 
						   "true" : "false"));
  Printf("Symmetrize ref. flow wrt eta = 0:\t%s", ((fFlowFlags & kSymEta) ? 
						   "true" : "false"));
  Printf("Use an eta-gap for ref. flow    :\t%s", ((fFlowFlags & kEtaGap) ? 
						   "true" : "false"));
  Printf("FMD sigma cut:               :\t%f", fFMDCut);
  Printf("SPD sigma cut:               :\t%f", fSPDCut);
  if ((fFlowFlags & kEtaGap)) 
    Printf("Eta gap:                     :\t%f", fEtaGap);
}
//_____________________________________________________________________
//
//
// EOF

