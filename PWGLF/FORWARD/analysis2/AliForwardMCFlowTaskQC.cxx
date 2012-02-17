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

//_____________________________________________________________________
AliForwardMCFlowTaskQC::AliForwardMCFlowTaskQC() :
  AliForwardFlowTaskQC(), 
  fBinsFMDTR(),          // List of FMDTR analysis objects
  fBinsSPDTR(),          // List of SPDTR analysis objects
  fBinsMC(),             // List of MC truth analysis objects
  fdNdedpMC(),           // MC truth d^2N/detadphi histogram
  fAliceCent4th(),       // Alice QC4 vs. centrality data points
  fAlicePt2nd4050(),     // Alice QC2 vs. pT data points
  fAlicePt4th3040(),     // Alice QC4 vs. pT data points
  fAlicePt4th4050(),     // Alice QC4 vs. pT data points
  fImpactParToCent(),    // Impact parameter to centrality graph
  fAddFlow(0),           // Add flow to MC truth
  fAddType(0),           // Add type of flow to MC truth
  fAddOrder(0)           // Add order of flow to MC truth        
   {} 
  //
  // Default Constructor
  //
//_____________________________________________________________________
AliForwardMCFlowTaskQC::AliForwardMCFlowTaskQC(const char* name) :
  AliForwardFlowTaskQC(name),
  fBinsFMDTR(),          // List of FMDTR analysis objects
  fBinsSPDTR(),          // List of SPDTR analysis objects
  fBinsMC(),             // List of MC truth analysis objects
  fdNdedpMC(),           // MC truth d^2N/detadphi histogram
  fAliceCent4th(),       // Alice QC4 vs. centrality data points
  fAlicePt2nd4050(),     // Alice QC2 vs. pT data points
  fAlicePt4th3040(),     // Alice QC4 vs. pT data points
  fAlicePt4th4050(),     // Alice QC4 vs. pT data points
  fImpactParToCent(),    // Impact parameter to centrality graph
  fAddFlow(0),           // Add flow to MC truth
  fAddType(0),           // Add type of flow to MC truth
  fAddOrder(0)           // Add order of flow to MC truth        
{ 
  // 
  // Constructor
  // Parameters:
  //  name: Name of task
  //
  fdNdedpMC = TH2D("fdNdedpMC", "fdNdedpMC", 48, -6., 6., 200, 0., 2.*TMath::Pi());
  fdNdedpMC.Sumw2();
  
  // Add parametrizations of ALICE data
  Double_t xCumulant4thTPCrefMultTPConlyAll[] = {2.5,7.5,15,25,35,45,55,65};
  Double_t yCumulant4thTPCrefMultTPConlyAll[] = {0.017855,0.032440,0.055818,0.073137,0.083898,0.086690,0.082040,0.077777};
  Int_t nPointsCumulant4thTPCrefMultTPConlyAll = sizeof(xCumulant4thTPCrefMultTPConlyAll)/sizeof(Double_t);
  fAliceCent4th = new TGraph(nPointsCumulant4thTPCrefMultTPConlyAll,xCumulant4thTPCrefMultTPConlyAll,
                                                        yCumulant4thTPCrefMultTPConlyAll);
 
  Double_t xCumulant2nd4050ALICE[] = {0.000000,0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.850000,0.950000,
  1.100000,1.300000,1.500000,1.700000,1.900000,2.250000,2.750000,3.250000,3.750000,4.500000};
  Double_t yCumulant2nd4050ALICE[] = {0.000000,0.043400,0.059911,0.073516,0.089756,0.105486,0.117391,0.128199,0.138013,
  0.158271,0.177726,0.196383,0.208277,0.216648,0.242954,0.249961,0.240131,0.269006,0.207796};
  Int_t nPointsCumulant2nd4050ALICE = sizeof(xCumulant2nd4050ALICE)/sizeof(Double_t);                                      
  fAlicePt2nd4050 = new TGraph(nPointsCumulant2nd4050ALICE,xCumulant2nd4050ALICE,yCumulant2nd4050ALICE);

  Double_t xCumulant4th3040ALICE[] = {0.00000,0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.850000,0.950000,
  1.100000,1.300000,1.500000,1.700000,1.900000,2.250000,2.750000,3.250000,3.750000,4.500000,
  5.500000,7.000000,9.000000};
  Double_t yCumulant4th3040ALICE[] = {0.000000,0.037071,0.048566,0.061083,0.070910,0.078831,0.091396,0.102026,0.109691,
  0.124449,0.139819,0.155561,0.165701,0.173678,0.191149,0.202015,0.204540,0.212560,0.195885,
  0.000000,0.000000,0.000000};
  Int_t nPointsCumulant4th3040ALICE = sizeof(xCumulant4th3040ALICE)/sizeof(Double_t);                                      
  fAlicePt4th3040 = new TGraph(nPointsCumulant4th3040ALICE,xCumulant4th3040ALICE,yCumulant4th3040ALICE);

  Double_t xCumulant4th4050ALICE[] = {0.000000,0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.850000,0.950000,
  1.100000,1.300000,1.500000,1.700000,1.900000,2.250000,2.750000,3.250000,3.750000,4.500000};
  Double_t yCumulant4th4050ALICE[] = {0.000000,0.038646,0.049824,0.066662,0.075856,0.081583,0.099778,0.104674,0.118545,
  0.131874,0.152959,0.155348,0.169751,0.179052,0.178532,0.198851,0.185737,0.239901,0.186098};
  Int_t nPointsCumulant4th4050ALICE = sizeof(xCumulant4th4050ALICE)/sizeof(Double_t);   
  fAlicePt4th4050 = new TGraph(nPointsCumulant4th4050ALICE, xCumulant4th4050ALICE, yCumulant4th4050ALICE);

  Double_t impactParam[] = {0.,1.75,4.225,5.965,7.765,9.215,10.46,11.565,12.575,13.515,16.679};
  Double_t centrality[] = {0.,2.5,7.5,15,25,35,45,55,65,75,90};

  Int_t nPoints = sizeof(impactParam)/sizeof(Double_t);
  fImpactParToCent = new TGraph(nPoints, impactParam, centrality);

}
//_____________________________________________________________________
AliForwardMCFlowTaskQC::AliForwardMCFlowTaskQC(const AliForwardMCFlowTaskQC& o) :
  AliForwardFlowTaskQC(o), 
  fBinsFMDTR(),                          // List of FMDTR analysis objects
  fBinsSPDTR(),                          // List of SPDTR analysis objects
  fBinsMC(),                             // List of MC truth analysis objects
  fdNdedpMC(o.fdNdedpMC),                // MC truth d^2N/detadphi histogram
  fAliceCent4th(o.fAliceCent4th),        // Alice QC4 vs. centrality data points
  fAlicePt2nd4050(o.fAlicePt2nd4050),    // Alice QC2 vs. pT data points
  fAlicePt4th3040(o.fAlicePt4th3040),    // Alice QC4 vs. pT data points
  fAlicePt4th4050(o.fAlicePt4th4050),    // Alice QC4 vs. pT data points
  fImpactParToCent(o.fImpactParToCent),  // Impact parameter to centrality graph
  fAddFlow(o.fAddFlow),                  // Add flow to MC truth
  fAddType(o.fAddType),                  // Add type of flow to MC truth
  fAddOrder(o.fAddOrder)                 // Add order of flow to MC truth        
   {} 
  //
  // Copy Constructor
  //
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
  fAliceCent4th    = o.fAliceCent4th;
  fAlicePt2nd4050  = o.fAlicePt2nd4050;
  fAlicePt4th3040  = o.fAlicePt4th3040;
  fAlicePt4th4050  = o.fAlicePt4th4050;
  fImpactParToCent = o.fImpactParToCent;
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

  for(Int_t n = 1; n <= 6; n++) {
    if (!fv[n]) continue;
    for (Int_t v = -10; v < 10; v++) {
      fBinsFMDTR.Add(new VertexBin(v, v+1, n, "FMDTR"));
      fBinsSPDTR.Add(new VertexBin(v, v+1, n, "SPDTR"));
      fBinsMC.Add(new VertexBin(v, v+1, n, "MC"));
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
  if (!aodfmult) return kFALSE;
  TH2D fmdTRdNdetadphi = aodfmult->GetHistogram();

  AliAODCentralMult* aodcmult = static_cast<AliAODCentralMult*>(fAOD->FindListObject("CentralClustersMC"));
  if (!aodcmult) return kFALSE;
  TH2D spdTRdNdetadphi = aodcmult->GetHistogram();
  
  TIter nextFMDTR(&fBinsFMDTR);
  VertexBin* bin = 0;
  while ((bin = static_cast<VertexBin*>(nextFMDTR()))) {
    if (bin->CheckVertex(fZvertex)) {
      if (!bin->FillHists(&fmdTRdNdetadphi)) return kFALSE;
      bin->CumulantsAccumulate(fCent);
    }
  }

  TIter nextSPDTR(&fBinsSPDTR);
  while ((bin = static_cast<VertexBin*>(nextSPDTR()))) {
    if (bin->CheckVertex(fZvertex)) {
      if (!bin->FillHists(&spdTRdNdetadphi)) return kFALSE;
      bin->CumulantsAccumulate(fCent);
    }
  }

  // Run analysis on MC branch
  if (!LoopAODMC()) return kFALSE;

  TIter nextMC(&fBinsMC);
  while ((bin = static_cast<VertexBin*>(nextMC()))) {
    if (bin->CheckVertex(fZvertex)) {
      if (!bin->FillHists(&fdNdedpMC)) return kFALSE;
      bin->CumulantsAccumulate(fCent);
    }
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

  TIter nextFMDTR(&fBinsFMDTR);
  VertexBin* bin = 0;
  while ((bin = static_cast<VertexBin*>(nextFMDTR()))) {
    bin->CumulantsTerminate(fSumList, fOutputList);
  }
  TIter nextSPDTR(&fBinsSPDTR);
  while ((bin = static_cast<VertexBin*>(nextSPDTR()))) {
    bin->CumulantsTerminate(fSumList, fOutputList);
  }
  TIter nextMC(&fBinsMC);
  while ((bin = static_cast<VertexBin*>(nextMC()))) {
    bin->CumulantsTerminate(fSumList, fOutputList);
  }

  TProfile2D* fmdHist = 0;
  TProfile2D* spdHist = 0;
  TProfile2D* mcHist = 0;

  for (Int_t i = 2; i <= 4; i += 2) {
    for (Int_t n = 1; n <= 6; n++) {
      if (!fv[n]) continue;
      fmdHist = (TProfile2D*)fOutputList->FindObject(Form("FMDQC%d_v%d_unCorr", i, n))
		  ->Clone(Form("FMDQC%d_v%d_Correction", i, n));
      spdHist = (TProfile2D*)fOutputList->FindObject(Form("SPDQC%d_v%d_unCorr", i, n))
		  ->Clone(Form("SPDQC%d_v%d_Correction", i, n));
      mcHist = (TProfile2D*)fOutputList->FindObject(Form("MCQC%d_v%d_unCorr", i, n));
     
      if (!fmdHist || !spdHist || !mcHist) {
	AliError(Form("Histogram missing, correction object not created for v%d", n));
	continue;
      }

      fmdHist->Divide(mcHist);
      spdHist->Divide(mcHist);
      fmdHist->SetTitle(Form("FMD QC{%d} v_{%d} Correction Object", i, n));
      fmdHist->SetTitle(Form("SPD QC{%d} v_{%d} Correction Object", i, n));

      fOutputList->Add(fmdHist);
      fOutputList->Add(spdHist);
    }
  }

}
//_____________________________________________________________________
Bool_t AliForwardMCFlowTaskQC::LoopAODMC()  
{
  // 
  // Loop over AliAODParticle branch and fill d^2N/detadphi-histograms.
  // Add flow if set to do so in AddTask function
  fdNdedpMC.Reset();

  //retreive MC particles from event
  TClonesArray* mcArray = (TClonesArray*)fAOD->FindListObject(AliAODMCParticle::StdBranchName());
  if(!mcArray){
//    AliWarning("No MC array found in AOD. Try making it again.");
    return kFALSE;
  }

  Double_t rp = 0;
  AliAODMCHeader* header = dynamic_cast<AliAODMCHeader*>(fAOD->FindListObject(AliAODMCHeader::StdBranchName()));
  if (!header) {
    AliWarning("No header file found.");
  }
  else {
    rp = header->GetReactionPlaneAngle();
  }

  Int_t ntracks = mcArray->GetEntriesFast();

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
    if (TMath::Abs(eta) < 6.) {
      // Add flow if it is in the argument
      if (fAddFlow.Length() > 1) {
        if (fAddFlow.Contains("pt"))
          weight *= AddptFlow(pT);
        if (fAddFlow.Contains("pid"))
          weight *= AddpidFlow(particle->PdgCode());
        if (fAddFlow.Contains("eta"))
          weight *= AddetaFlow(eta);
        if (fAddFlow.Contains("cent")) 
          weight *= fAliceCent4th->Eval(fCent)/fAliceCent4th->Eval(45);
        
        weight *= 20*2.*TMath::Cos((Double_t)fAddOrder*(phi-rp)); 
        weight += 1;
      }
      fdNdedpMC.Fill(eta, phi, weight);
    }
  }

  return kTRUE;
}
//_____________________________________________________________________
Double_t AliForwardMCFlowTaskQC::AddptFlow(Double_t pt = 0) const 
{
  //
  // Add pt dependent flow factor
  // Parameters:
  //  pt: pT parametrization to use
  //
  Double_t weight = 0;

  switch(fAddType) 
  {
    case 1: weight = fAlicePt2nd4050->Eval(pt)*0.5+fAlicePt4th4050->Eval(pt)*0.5;
            break;
    case 2: weight = fAlicePt2nd4050->Eval(pt);
            break;
    case 3: weight = fAlicePt4th3040->Eval(pt);
            break;
    case 4: weight = fAlicePt4th4050->Eval(pt);
            break;
  } 
  
  return weight;
}
//_____________________________________________________________________
Double_t AliForwardMCFlowTaskQC::AddpidFlow(Int_t id = 0) const 
{
  //
  // Add pid dependent flow factor 
  // Parameters:
  //  id: choose PID dependent setup
  //
  Double_t weight = 0;

  switch(fAddType)
  {
    case 1: if (TMath::Abs(id) ==  211) // pion flow
              weight = 1.3;
            else if (TMath::Abs(id) ==  2212) // proton flow
              weight = 1.;
            else 
              weight = 0.7;
            break;
    case 2: weight = 1.207;
            break;
  }

  return weight;
}
//_____________________________________________________________________
Double_t AliForwardMCFlowTaskQC::AddetaFlow(Double_t eta = 0) const 
{
  //
  // Add eta dependent flow factor 
  // Parameters:
  //  eta: choose v_n(eta) shape
  //
  Double_t weight = 0;

  TF1 gaus = TF1("gaus", "gaus", -6, 6);

  switch(fAddType)
  {
     case 1: gaus.SetParameters(0.1, 0., 9);
             break;
     case 2: gaus.SetParameters(0.1, 0., 3);
             break;
     case 3: gaus.SetParameters(0.1, 0., 15);
             break;
  }

  weight = gaus.Eval(eta);

  return weight;
}
//_____________________________________________________________________
Double_t AliForwardMCFlowTaskQC::GetCentFromB() const
{
  //
  // Get centrality from MC impact parameter.
  // Values taken from: https://twiki.cern.ch/twiki/bin/viewauth/ALICE/CentStudies
  //
  Double_t cent = -1.;
  Double_t b = -1.;
  AliAODMCHeader* header = (AliAODMCHeader*)fAOD->FindListObject(AliAODMCHeader::StdBranchName());
  if (!header) return cent;
  b = header->GetImpactParameter();

  cent = fImpactParToCent->Eval(b);

  return cent;
}
//_____________________________________________________________________
//
//
// EOF

