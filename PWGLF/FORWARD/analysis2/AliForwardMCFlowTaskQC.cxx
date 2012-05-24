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
#include "TProfile2D.h"
#include "AliAODEvent.h"
#include "AliAODForwardMult.h"
#include "AliAODCentralMult.h"

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
    fAddFlow(0),           // Add flow to MC truth
    fAddType(0),           // Add type of flow to MC truth
    fAddOrder(0)           // Add order of flow to MC truth        
{ 
  // 
  // Constructor
  // Parameters:
  //  name: Name of task
  //
  fdNdedpMC = TH2D("fdNdedpMC", "fdNdedpMC", 
		   48, -6., 6., 200, 0., 2.*TMath::Pi());
  fdNdedpMC.Sumw2();
  
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

  for(UShort_t n = 1; n <= 6; n++) {
    if (!fv[n]) continue;
    for (Int_t v = 1; v <= fVtxAxis->GetNbins(); v++) {
      Int_t vL = Int_t(fVtxAxis->GetBinLowEdge(v));
      Int_t vH = Int_t(fVtxAxis->GetBinUpEdge(v));
      
      fBinsFMDTR.Add(new VertexBin(vL, vH, n, "FMDTR"));
      fBinsSPDTR.Add(new VertexBin(vL, vH, n, "SPDTR"));
      fBinsMC.Add(new VertexBin(vL, vH, n, "MC"));
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

  fWeights.Init(fSumList);
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
  if (!aodfmult || !aodcmult) return kFALSE;
  
  // if objects are present, get histograms
  const TH2D& fmdTRdNdetadphi = aodfmult->GetHistogram();
  const TH2D& spdTRdNdetadphi = aodcmult->GetHistogram();

  // Run analysis on tr refs
  Int_t vtx = fVtxAxis->FindBin(fVtx)-1;
  if (!FillVtxBinList(fBinsFMDTR, fmdTRdNdetadphi, vtx)) return kFALSE;
  if (!FillVtxBinList(fBinsSPDTR, spdTRdNdetadphi, vtx)) return kFALSE;

  // Run analysis on MC branch
  if (!LoopAODMC()) return kFALSE;
  if (!FillVtxBinList(fBinsMC, fdNdedpMC, vtx)) return kFALSE;

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
    fHistCent->Fill(fCent);
    return kTRUE;
  }
  else  return AliForwardFlowTaskQC::GetCentrality(aodfm);
}
//_____________________________________________________________________
Bool_t AliForwardMCFlowTaskQC::LoopAODMC()  
{
  // 
  // Loop over AliAODParticle branch and fill d^2N/detadphi-histograms.
  // Add flow if set to do so in AddTask function
  fdNdedpMC.Reset();

  //retreive MC particles from event
  TClonesArray* mcArray = 
    static_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::
						    StdBranchName()));
  if(!mcArray){
    //    AliWarning("No MC array found in AOD. Try making it again.");
    return kFALSE;
  }

  Double_t rp = 0;
  AliAODMCHeader* header = 
    dynamic_cast<AliAODMCHeader*>(fAOD->FindListObject(AliAODMCHeader::
						       StdBranchName()));
  if (!header) 
    AliWarning("No header file found.");
  else 
    rp = header->GetReactionPlaneAngle();

  Int_t ntracks = mcArray->GetEntriesFast();

  UShort_t flowFlags = 0;
  if (fAddFlow.Length() > 1) {
    if (fAddFlow.Contains("pt"))   flowFlags |= AliForwardFlowWeights::kPt;
    if (fAddFlow.Contains("pid"))  flowFlags |= AliForwardFlowWeights::kPt;
    if (fAddFlow.Contains("eta"))  flowFlags |= AliForwardFlowWeights::kEta;
    if (fAddFlow.Contains("cent")) flowFlags |= AliForwardFlowWeights::kCent;
  }


  // Track loop: chek how many particles will be accepted
  Double_t weight = 0;
  for (Int_t it = 0; it < ntracks; it++) {
    weight = 1;
    AliAODMCParticle* particle = (AliAODMCParticle*)mcArray->At(it);
    if (!particle) {
      AliError(Form("Could not receive track %d", it));
      continue;
    }
    if (!particle->IsPhysicalPrimary()) continue;
    if (particle->Charge() == 0) continue;
    Double_t pT = particle->Pt();
    Double_t eta = particle->Eta();
    Double_t phi = particle->Phi();
    if (eta > -4. && eta < 5.) {
      // Add flow if it is in the argument
      if (flowFlags != 0) 
	weight = fWeights.CalcWeight(eta, pT, phi, particle->PdgCode(), 
				     rp, fCent, fAddType, fAddOrder, 
				     flowFlags) + 1;
      fdNdedpMC.Fill(eta, phi, weight);
    }
  }

  return kTRUE;
}
//_____________________________________________________________________
Double_t AliForwardMCFlowTaskQC::GetCentFromB() const
{
  //
  // Get centrality from MC impact parameter.
  //
  Double_t cent = -1.;
  Double_t b = -1.;
  AliAODMCHeader* header = 
    static_cast<AliAODMCHeader*>(fAOD->FindListObject(AliAODMCHeader::
						      StdBranchName()));
  if (!header) return cent;
  b = header->GetImpactParameter();

  cent = fImpactParToCent->Eval(b);

  return cent;
}
//_____________________________________________________________________
//
//
// EOF

